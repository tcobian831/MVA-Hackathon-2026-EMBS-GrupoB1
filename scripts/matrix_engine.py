#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Corrected gene matrix adapted from Lucia's wmatrix.py.
Install: python -m pip install numpy pandas scipy
Run: python wmatrix_corrected.py
Input: embedded original table of 31 variants in 28 genes.
Complete MANUAL_SCORES to obtain the final ranking and sensitivity analysis.
Each run creates a new timestamped output directory to avoid stale results.
All listed variants contribute to gene averages; none is excluded.
"""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
from datetime import datetime
OUTPUT_DIR = Path(__import__('os').environ.get('MVA_OUTPUT_DIR', str(ROOT / 'results' / ('wmatrix_results_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f')))))
import io
import numpy as np
import pandas as pd
from scipy.stats import binomtest

# ===================== MANUAL EDITING SECTION =====================
# Key: exact gene symbol.
# Value: a LIST of THREE numeric scores, in this exact order:
#     [MVA link, phenotype fit, literature evidence]
# Each completed score must be a finite float between 0.0 and 2.0 inclusive.
# Integers in this range are accepted and converted to floats.
# None is an EMPTY placeholder: replace every None with your own score.
# Completed-format example only: "GENE_SYMBOL": [1.5, 2.0, 0.5]
# All 28 genes intentionally start empty. Missing values are NEVER set to zero.
MANUAL_SCORES = {
    "ACLY": [None, None, None],
    "ADAMTS1": [None, None, None],
    "ASAP3": [None, None, None],
    "BUB1B": [None, None, None],
    "CD82": [None, None, None],
    "CITED2": [None, None, None],
    "CNOT10": [None, None, None],
    "CTU2": [None, None, None],
    "DIDO1": [None, None, None],
    "FRZB": [None, None, None],
    "IGBP1": [None, None, None],
    "LRP4": [None, None, None],
    "LZTR1": [None, None, None],
    "MEP1B": [None, None, None],
    "PIK3C3": [None, None, None],
    "PRMT6": [None, None, None],
    "ROCK1": [None, None, None],
    "SIM2": [None, None, None],
    "SLC8A1": [None, None, None],
    "SMC4": [None, None, None],
    "SOCS5": [None, None, None],
    "SVIL": [None, None, None],
    "TACC2": [None, None, None],
    "THUMPD3": [None, None, None],
    "TRDMT1": [None, None, None],
    "TSSK6": [None, None, None],
    "ZNF197": [None, None, None],
    "ZNF418": [None, None, None],
}

# Joint phase configuration; edit only when supported by evidence.
# "trans": two heterozygous variants at different sites, confirmed in trans.
# "unknown": two heterozygous variants, phase unknown.
# "cis": two heterozygous variants on the same haplotype.
# "trans_split_site": two different ALT alleles at one site (original GT 1/2).
# These preserve the original assumptions. Split rows alone do not establish
# the original multiallelic 1/2 genotype.
PHASE = {"BUB1B": "unknown", "ASAP3": "unknown",
         "ADAMTS1": "cis"}
# Preserve the original assumption of a single X copy in an XY sample.
ASSUME_XY = True
W = {"mva": 5.0, "pheno": 4.5, "effect": 3.5, "zyg": 2.5,
     "freq": 2.0, "cadd": 1.5, "lit": 1.0, "call": 0.25}
# =================== END OF MANUAL EDITING SECTION ===================

# Data transcribed from document. grpmax_pct = max group frequency in % (blank = ND)
CSV = (ROOT / "results" / "curated_variants.csv").read_text()
v = pd.read_csv(io.StringIO(CSV))

# Sanity checks against the document
assert len(v) == 31 and v["gene"].nunique() == 28
assert v["grpmax_pct"].notna().sum() == 16 and v["cadd"].isna().sum() == 1

# Allele balance: fraction of reads supporting the ALT allele
v["ad_total"] = v["ad_ref"] + v["ad_alt"]
v["alt_frac"] = (v["ad_alt"] / v["ad_total"]).round(2)

# For heterozygous calls we expect ~50% ALT reads. Binomial test p-value against 0.5
def het_p(r):
    if r["gt"] in ("0/1", "1/0") and r["ad_total"] > 0:
        return binomtest(int(r["ad_alt"]), int(r["ad_total"]), 0.5).pvalue
    return np.nan
v["het_balance_p"] = v.apply(het_p, axis=1).round(3)

# Call-quality flags
v["flag_low_depth"] = v["dp"] <= 25
v["flag_low_gq"]    = v["gq"] <= 50
v["flag_ab_skew"]   = v["het_balance_p"] < 0.05         # ~1-2 of 30 flagged by chance alone
v["flag_hom_with_ref_reads"] = (v["gt"] == "1/1") & (v["ad_ref"] > 0)
# Missing CADD and repeat/long-allele context are NOT quality flags.
v["n_flags"] = v[[c for c in v.columns if c.startswith("flag_")]].sum(axis=1)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", None)
print(v[["id", "gene", "variant", "gt", "dp", "gq", "ad_ref", "ad_alt", "alt_frac",
         "het_balance_p", "n_flags"]].to_string(index=False))
#----------------------------------------------- Step 2 Weigthed matrix
# ---------- Automatic scores (0-2) ----------
LOF = ["stop_gained", "frameshift_del", "frameshift_ins", "splice_donor_del"]
v["s_effect"] = v["type"].apply(lambda t: 2.0 if t in LOF else (1.0 if t == "missense" else 0.0))
# Team-agreed adjustments; these are scoring decisions, not molecular validation.
EFFECT_OVERRIDE = {2: 1, 3: 1,    # ADAMTS1: cis pair, net deletion of 21 bp may restore the frame
                   11: 1,         # CTU2: possible splice rescue
                   12: 1.5,       # DIDO1: revised adjustment
                   26: 1.5,       # SOCS5: revised adjustment
                   32: 1.5}       # ZNF197: revised adjustment
for i, s in EFFECT_OVERRIDE.items():
    v.loc[v["id"] == i, "s_effect"] = s

# Anchors use the ENTIRE list, ignoring missing data.
CADD_MAX = float(v["cadd"].max())        # Currently 49.
FREQ_MIN = float(v["grpmax_pct"].min())  # Currently 0.00147 PERCENT.
if not np.isfinite(CADD_MAX) or CADD_MAX <= 20:
    raise ValueError("CADD normalization requires a maximum above 20.")
if not np.isfinite(FREQ_MIN) or not 0 <= FREQ_MIN < 0.1:
    raise ValueError("Frequency minimum must be in [0, 0.1) percent.")

def cadd_score(c):
    """Missing -> 1; CADD 20 -> 0; maximum -> 2; bounded to [0, 2]."""
    return 1.0 if pd.isna(c) else float(np.clip(2 * (c - 20) / (CADD_MAX - 20), 0, 2))

def frequency_score(f):
    """f in PERCENT: missing -> 1; minimum -> 2; f >= 0.1% -> 0."""
    return 1.0 if pd.isna(f) else float(np.clip(2 * (0.1 - f) / (0.1 - FREQ_MIN), 0, 2))

v["s_cadd"] = v["cadd"].map(cadd_score)
v["s_freq"] = v["grpmax_pct"].map(frequency_score)
v["s_call"] = v["n_flags"].apply(lambda n: 2 if n == 0 else (1 if n == 1 else 0))

# Validate the manual dictionary and preserve empty entries as NaN.
if set(MANUAL_SCORES) != set(v["gene"]):
    raise ValueError("MANUAL_SCORES must contain exactly the 28 gene symbols.")
for gene, values in MANUAL_SCORES.items():
    if not isinstance(values, list) or len(values) != 3:
        raise ValueError(f"{gene}: provide a list of exactly three scores.")
    for value in values:
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, (int, float))
            or not np.isfinite(value) or not 0 <= value <= 2
        ):
            raise ValueError(f"{gene}: use finite scores in [0, 2], or None.")
for index, column in enumerate(["s_mva", "s_pheno", "s_lit"]):
    v[column] = v["gene"].map(
        lambda gene: np.nan if MANUAL_SCORES[gene][index] is None
        else float(MANUAL_SCORES[gene][index])
    )

# Configuration is joint evidence about a gene, not a per-variant score.
CONFIG_TEXT = {   # configuration -> (description, are the variants on different alleles?)
    "single_het":       ("single heterozygous variant", "n/a (one allele only)"),
    "hom":              ("homozygous (same variant on both alleles)", "both alleles (same variant)"),
    "X_hemizygous":     ("X-linked, 1/1 likely hemizygous if XY", "n/a (single X copy)"),
    "trans":            ("two heterozygous variants at different sites, confirmed in trans", "yes (trans)"),
    "unknown":          ("two heterozygous variants, phase unknown", "unknown"),
    "trans_split_site": ("two different ALT alleles at one site (GT 1/2)", "yes by genotype, weak read support"),
    "cis":              ("two heterozygous variants on the same haplotype", "no (cis)"),
}
ZYG_SCORE = {"hom": 2, "trans_split_site": 2, "trans": 2, "unknown": 1, "X_hemizygous": 1, "single_het": 0.5, "cis": 0.5}

def vdesc(r):
    af = "AF ND" if pd.isna(r["grpmax_pct"]) else f"AF {r['grpmax_pct']}%"
    cadd = "CADD ND" if pd.isna(r["cadd"]) else f"CADD {r['cadd']}"
    return f"{r['variant']} ({r['chrom']}:{r['pos']}, {r['gt']}, {cadd}, {af})"

# ---------- One unit per gene: include ALL listed variants ----------
rows = []
for gene, g in v.groupby("gene"):
    assert len(g) <= 2, f"{gene} has more than 2 variants: define how to pair them"
    if len(g) == 1:
        r = g.iloc[0]
        set_cfg = "X_hemizygous" if r["chrom"] == "chrX" and ASSUME_XY else ("hom" if r["gt"] == "1/1" else "single_het")
    else:
        set_cfg = PHASE[gene]
    if set_cfg not in ZYG_SCORE:
        raise ValueError(f"{gene}: unsupported configuration {set_cfg}.")
    if set_cfg == "trans" and g["pos"].nunique() != 2:
        raise ValueError(f"{gene}: trans requires two different positions.")
    if set_cfg == "trans_split_site" and g["pos"].nunique() != 1:
        raise ValueError(f"{gene}: trans_split_site requires a single position.")
    # ALL listed variants contribute equally; there are no exclusion gates.
    members, ranked_cfg, note = g, set_cfg, ""
    rows.append({
        "gene": gene,
        "n_variants_in_set": len(g),
        "all_variants": " + ".join(vdesc(r) for _, r in g.iterrows()),
        "set_configuration": CONFIG_TEXT[set_cfg][0],
        "different_alleles": CONFIG_TEXT[set_cfg][1],
        "ranked_unit": " + ".join(members["variant"]),
        "ranked_as": CONFIG_TEXT[ranked_cfg][0],
        "note": note,
        "u_effect": members["s_effect"].mean(),
        "u_cadd":   members["s_cadd"].mean(),
        "u_freq":   members["s_freq"].mean(),
        "u_call":   members["s_call"].mean(),        # average ALL listed variants
        "u_zyg":    ZYG_SCORE[ranked_cfg],
        "u_pheno":  members["s_pheno"].iloc[0],      # gene-level, same for all its variants
        "u_mva":    members["s_mva"].iloc[0],
        "u_lit":    members["s_lit"].iloc[0],
    })
units = pd.DataFrame(rows)
assert len(units) == 28

# ---------- TABLE BEFORE THE RANKING: genes, specific variants and allele configuration ----------
pd.set_option("display.max_colwidth", 220); pd.set_option("display.width", 300)
table = units[["gene", "n_variants_in_set", "all_variants", "set_configuration",
               "different_alleles", "ranked_unit", "ranked_as", "note"]]
print(table.sort_values(["n_variants_in_set", "gene"], ascending=[False, True]).to_string(index=False))
OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
table.to_csv(OUTPUT_DIR / "gene_table_before_ranking.csv", index=False)

# Weighted total: fixed denominator; any missing component propagates.
def unit_totals(w):
    return sum(w[k] * units[f"u_{k}"] for k in w) / (2 * sum(w.values())) * 100

units["manual_complete"] = units[["u_mva", "u_pheno", "u_lit"]].notna().all(axis=1)
units["total"] = unit_totals(W)
units.to_csv(OUTPUT_DIR / "gene_scores.csv", index=False)
v.to_csv(OUTPUT_DIR / "variant_matrix.csv", index=False)
pending = units.loc[~units["manual_complete"], ["gene", "u_mva", "u_pheno", "u_lit"]]
pending.to_csv(OUTPUT_DIR / "manual_scores_pending.csv", index=False)
print(f"Output folder: {OUTPUT_DIR.resolve()}")
if not pending.empty:
    print(f"Complete MANUAL_SCORES for {len(pending)} genes and rerun.")
    print("Final ranking and sensitivity are pending; missing scores are never zero.")
    raise SystemExit(0)
# Rank unrounded totals, so display rounding does not create artificial ties.
ranking = units.sort_values("total", ascending=False).reset_index(drop=True)
ranking["rank"] = ranking["total"].rank(ascending=False, method="min").astype(int)   # ties share the rank
ranking["in_top15"] = ranking["rank"] <= 15

print("\nFULL RANKING (all 28 genes):")
print(ranking[["rank", "gene", "ranked_unit", "ranked_as", "total", "in_top15"]].to_string(index=False))

# ---------- Sensitivity: how stable is the ranking if the weights change? ----------
rng = np.random.default_rng(1)
N = 10000
fields = list(W)
base_weights = np.array([W[k] for k in fields], dtype=float)
if not np.isfinite(base_weights).all() or (base_weights < 0).any() or base_weights.sum() <= 0:
    raise ValueError("Weights must be finite, nonnegative, and have a positive sum.")
# Draw each field's weight independently, uniformly from 50% to 150% of baseline.
simulation_weights = rng.uniform(0.5, 1.5, size=(N, len(fields))) * base_weights
field_scores = units[[f"u_{k}" for k in fields]].to_numpy(dtype=float)
# Every simulation score is a weighted mean on the original 0-2 scale.
simulation_scores = (field_scores @ simulation_weights.T) / simulation_weights.sum(axis=1)
if not np.isfinite(simulation_scores).all() or (simulation_scores < 0).any() or (simulation_scores > 2).any():
    raise ValueError("Sensitivity scores must be finite and between 0 and 2.")
S = pd.DataFrame(simulation_scores, index=units["gene"])
R = S.rank(axis=0, ascending=False, method="min")
sens = pd.DataFrame({"median_score": S.median(axis=1),
                     "pct_in_top5": (R <= 5).mean(axis=1) * 100,
                     "pct_in_top10": (R <= 10).mean(axis=1) * 100,
                     "pct_in_top15": (R <= 15).mean(axis=1) * 100,
                     "pct_rank20_or_later": (R >= 20).mean(axis=1) * 100,
                     "median_rank": R.median(axis=1)})
print("\nSENSITIVITY:")
print(sens.sort_values("median_rank").to_string())

ranking["weighted_mean"] = ranking["total"] / 50.0
ranking = ranking.join(sens, on="gene", validate="one_to_one")
# Combine scores on the SAME 0-2 scale; total remains the original 0-100 output.
ranking["definitive_score"] = 0.8 * ranking["weighted_mean"] + 0.2 * ranking["median_score"]
summary = ranking[["gene"] + [f"u_{k}" for k in fields] + [
    "weighted_mean", "total", "median_score", "pct_in_top5", "pct_in_top10",
    "pct_in_top15", "pct_rank20_or_later", "definitive_score"]].copy()
summary = summary.rename(columns={f"u_{k}": k for k in fields})
summary = summary.sort_values(["definitive_score", "gene"], ascending=[False, True])
print("\nDEFINITIVE RANKING (0-2):")
print(summary[["gene", "weighted_mean", "median_score", "definitive_score"]].to_string(index=False))
ranking.to_csv(OUTPUT_DIR / "gene_ranking.csv", index=False)
sens.to_csv(OUTPUT_DIR / "sensitivity.csv", index_label="gene")
summary.to_csv(OUTPUT_DIR / "gene_ranking_summary.csv", index=False)
