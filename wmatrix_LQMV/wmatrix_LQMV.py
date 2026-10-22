# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:52:57 2026

@author: mavec
"""
#Here what I do is one BIG table with the 33 variants and their quality data
import io
import numpy as np
import pandas as pd
from scipy.stats import binomtest

# Data transcribed from document. grpmax_pct = max group frequency in % (blank = ND)
CSV = """id,gene,variant,chrom,pos,type,cadd,grpmax_pct,gt,dp,gq,ad_ref,ad_alt
1,ACLY,Arg244Trp,chr17,41907459,missense,28.3,,0/1,23,99,13,10
2,ADAMTS1,del7bp,chr21,26844275,frameshift_del,20.8,0.008821,0/1,58,99,24,34
3,ADAMTS1,del14bp,chr21,26844283,frameshift_del,24.0,0.00882,0/1,62,99,28,34
4,ASAP3,del119bp_downstream,chr1,23425132,downstream_del,,0.03293,0/1,55,99,30,25
5,ASAP3,Ala430Val,chr1,23437183,missense,24.8,0.4837,0/1,37,99,21,16
6,BUB1B,Leu737Ter,chr15,40209701,stop_gained,36.0,0.00735,0/1,46,99,21,25
7,BUB1B,Asn1002Lys,chr15,40220612,missense,24.5,,0/1,28,99,15,13
8,CD82,Tyr57Asn,chr11,44605090,missense,24.9,,0/1,35,99,17,18
9,CITED2,His160Pro,chr6,139373466,missense,25.8,0.004414,0/1,38,99,25,13
10,CNOT10,Ser110Arg,chr3,32708720,missense,23.1,,0/1,44,99,18,26
11,CTU2,donor_del6bp,chr16,88714226,splice_donor_del,33.0,,1/1,27,82,0,27
12,DIDO1,donor_del7bp,chr20,62926430,splice_donor_del,23.4,,0/1,24,99,10,14
13,FRZB,Cys180Phe,chr2,182842531,missense,28.7,0.01764,0/1,44,99,28,16
14,IGBP1,Met151Leu,chrX,70134785,missense,22.9,0.005617,1/1,23,69,0,23
15,LRP4,Arg1365Gln,chr11,46874935,missense,27.7,0.00147,0/1,46,99,19,27
16,LZTR1,Tyr748Ter,chr22,20996720,stop_gained,49.0,,0/1,48,99,24,24
17,MAD2L1BP,ins184_A,chr6,43738633,window_insertion,,,1/0,8,99,7,1
18,MAD2L1BP,ins184_B,chr6,43738633,window_insertion,,,0/1,8,99,8,0
19,MEP1B,Ser165Thr,chr18,32204306,missense,25.9,,0/1,49,99,33,16
20,PIK3C3,Met777Leu,chr18,42057948,missense,24.8,,0/1,38,99,15,23
21,PRMT6,Gly12Arg,chr1,107056749,missense,24.6,0.0147,0/1,37,99,13,24
22,ROCK1,Gln827Arg,chr18,20984360,missense,23.8,0.0147,0/1,39,99,19,20
23,SIM2,Glu224Lys,chr21,36726245,missense,33.0,0.00294,0/1,60,99,35,25
24,SLC8A1,Ile611Thr,chr2,40178470,missense,25.0,,0/1,41,99,23,18
25,SMC4,Asp256Asn,chr3,160411998,missense,33.0,0.00147,0/1,30,99,18,12
26,SOCS5,Gln410Ter,chr2,46759758,stop_gained,43.0,,0/1,30,99,13,17
27,SVIL,Ala2072Thr,chr10,29463555,missense,26.9,,0/1,31,99,11,20
28,TACC2,Ser1562Ter,chr10,122087185,stop_gained,34.0,0.00147,0/1,55,99,23,32
29,THUMPD3,Glu288Ter,chr3,9374570,stop_gained,44.0,,0/1,38,99,20,18
30,TRDMT1,Leu377Phe,chr10,17149087,missense,28.9,,0/1,51,99,24,27
31,TSSK6,Thr17Ile,chr19,19515378,missense,29.4,,0/1,58,99,30,28
32,ZNF197,Tyr790fs,chr3,44643497,frameshift_del,33.0,0.007351,0/1,43,99,16,27
33,ZNF418,Leu634fs,chr19,57926279,frameshift_ins,22.4,0.9449,1/1,15,45,0,15
"""
v = pd.read_csv(io.StringIO(CSV))

# Sanity checks against the document
assert len(v) == 33 and v["gene"].nunique() == 29
assert v["grpmax_pct"].notna().sum() == 16 and v["cadd"].isna().sum() == 3

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
v["flag_low_depth"] = v["dp"] <= 10                     # same cutoff as the agreed DP > 10
v["flag_low_gq"]    = v["gq"] <= 50
v["flag_ab_skew"]   = v["het_balance_p"] < 0.05         # ~1-2 of 30 flagged by chance alone
v["flag_hom_with_ref_reads"] = (v["gt"] == "1/1") & (v["ad_ref"] > 0)
v["flag_no_cadd"]   = v["cadd"].isna()
# Repeat / long-allele context, taken from the document (ASAP3 deletion, CTU2, MAD2L1BP insertions)
v["flag_repeat_context"] = v["id"].isin([4, 11, 17, 18])
v["n_flags"] = v[[c for c in v.columns if c.startswith("flag_")]].sum(axis=1)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", None)
print(v[["id", "gene", "variant", "gt", "dp", "gq", "ad_ref", "ad_alt", "alt_frac",
         "het_balance_p", "n_flags"]].to_string(index=False))
#----------------------------------------------- Step 2 Weigthed matrix
# ---------- Automatic scores (0-2) ----------
LOF = ["stop_gained", "frameshift_del", "frameshift_ins", "splice_donor_del"]
v["s_effect"] = v["type"].apply(lambda t: 2 if t in LOF else (1 if t == "missense" else 0))
#doubts in the loss-of-function effect, lower the score to 1 (edited) PREGUNTAR en reunion
EFFECT_OVERRIDE = {2: 1, 3: 1,    # ADAMTS1: cis pair, net deletion of 21 bp may restore the frame
                   11: 1,         # CTU2: possible splice rescue
                   12: 1,         # DIDO1: isoform-dependent
                   26: 1,         # SOCS5: stop in the last exon
                   32: 1}         # ZNF197: frameshift in the last exon
for i, s in EFFECT_OVERRIDE.items():
    v.loc[v["id"] == i, "s_effect"] = s

v["s_cadd"] = v["cadd"].apply(lambda c: 0 if pd.isna(c) else (2 if c >= 30 else (1 if c >= 24 else 0)))
v["s_freq"] = v["grpmax_pct"].apply(lambda f: 1 if pd.isna(f) else (2 if f < 0.01 else (1 if f < 0.1 else 0)))  # missing = neutral
v["s_call"] = v["n_flags"].apply(lambda n: 2 if n == 0 else (1 if n == 1 else 0))

# Zygosity / second hit (edit after meeting): 2 hom autosomal, 1 pair of unknown phase or trans, 0 cis, 0.5 single het
ZYG = {"CTU2": 2, "ZNF418": 2, "BUB1B": 1, "ASAP3": 1, "MAD2L1BP": 1, "IGBP1": 1, "ADAMTS1": 0}
v["s_zyg"] = v["gene"].map(ZYG).fillna(0.5)

# ---------- Manual scores (LQMV @ZELTZ1N) ----------
# Manual gene-level annotation: 0-2.
#
# pheno_fit:
#   0 = no convincing phenotype match identified
#   1 = partial match with one or more clinical features
#   2 = strong match with the overall clinical presentation
#
# mva_link:
#   0 = no relevant link identified
#   1 = indirect mechanistic link to mitosis, chromosome stability,
#       DNA repair or ploidy
#   2 = established human association with MVA
#
# lit:
#   0 = insufficient relevant evidence
#   1 = limited, indirect or preliminary evidence
#   2 = substantial relevant clinical or mechanistic literature
#
# These scores are gene-level. They do not classify individual variants.

MANUAL = {
    "ACLY":      dict(pheno_fit=1, mva_link=1, lit=2),
    "ADAMTS1":   dict(pheno_fit=1, mva_link=0, lit=1),
    "ASAP3":     dict(pheno_fit=1, mva_link=0, lit=1),
    "BUB1B":     dict(pheno_fit=2, mva_link=2, lit=2),
    "CD82":      dict(pheno_fit=1, mva_link=0, lit=1),
    "CITED2":    dict(pheno_fit=1, mva_link=0, lit=2),
    "CNOT10":    dict(pheno_fit=1, mva_link=1, lit=2),
    "CTU2":      dict(pheno_fit=1, mva_link=0, lit=2),
    "DIDO1":     dict(pheno_fit=1, mva_link=1, lit=2),
    "FRZB":      dict(pheno_fit=1, mva_link=0, lit=2),
    "IGBP1":     dict(pheno_fit=0, mva_link=1, lit=1),
    "LRP4":      dict(pheno_fit=1, mva_link=0, lit=2),
    "LZTR1":     dict(pheno_fit=1, mva_link=0, lit=2),
    "MAD2L1BP":  dict(pheno_fit=2, mva_link=2, lit=2),
    "MEP1B":     dict(pheno_fit=1, mva_link=0, lit=1),
    "PIK3C3":    dict(pheno_fit=1, mva_link=0, lit=2),
    "PRMT6":     dict(pheno_fit=0, mva_link=1, lit=2),
    "ROCK1":     dict(pheno_fit=0, mva_link=1, lit=1),
    "SIM2":      dict(pheno_fit=1, mva_link=0, lit=1),
    "SLC8A1":    dict(pheno_fit=1, mva_link=0, lit=2),
    "SMC4":      dict(pheno_fit=1, mva_link=1, lit=1),
    "SOCS5":     dict(pheno_fit=1, mva_link=0, lit=1),
    "SVIL":      dict(pheno_fit=1, mva_link=1, lit=2),
    "TACC2":     dict(pheno_fit=1, mva_link=1, lit=2),
    "THUMPD3":   dict(pheno_fit=1, mva_link=0, lit=2),
    "TRDMT1":    dict(pheno_fit=1, mva_link=1, lit=2),
    "TSSK6":     dict(pheno_fit=0, mva_link=0, lit=1),
    "ZNF197":    dict(pheno_fit=0, mva_link=0, lit=1),
    "ZNF418":    dict(pheno_fit=0, mva_link=0, lit=1),
}

# Validate that the library covers all 29 target genes.
assert set(MANUAL) == set(v["gene"].unique()), (
    f"Gene mismatch. Missing: {set(v['gene'].unique()) - set(MANUAL)}; "
    f"Extra: {set(MANUAL) - set(v['gene'].unique())}"
)

assert all(
    set(scores) == {"pheno_fit", "mva_link", "lit"}
    and all(score in (0, 1, 2) for score in scores.values())
    for scores in MANUAL.values()
)

for k, col in [("pheno_fit", "s_pheno"), ("mva_link", "s_mva"), ("lit", "s_lit")]:
    v[col] = v["gene"].map(lambda g: MANUAL.get(g, {}).get(k, np.nan))
missing = sorted(v.loc[v["s_pheno"].isna(), "gene"].unique())
print(len(missing), "genes still without manual scores")
for col in ("s_pheno", "s_mva", "s_lit"):
    v[col] = v[col].fillna(0)

# ---------- Gates: now only warning flags, every gene is still ranked ----------
v["gate_no_effect"] = v["type"].isin(["downstream_del", "window_insertion"])   # no demonstrated effect on the gene
v["gate_bad_call"]  = v["flag_low_depth"] & v["flag_repeat_context"]
v["gate_flag"]      = v["gate_no_effect"] | v["gate_bad_call"]

# ---------- Phase of the genes with 2 variants (edit if groups disagree (preguntar en meet)) ----------
PHASE = {"BUB1B": "unknown", "ASAP3": "unknown", "MAD2L1BP": "trans_split_site", "ADAMTS1": "cis"}

CONFIG_TEXT = {   # configuration -> (description, are the variants on different alleles?)
    "single_het":       ("single heterozygous variant", "n/a (one allele only)"),
    "hom":              ("homozygous (same variant on both alleles)", "both alleles (same variant)"),
    "X_hemizygous":     ("X-linked, 1/1 likely hemizygous if XY", "n/a (single X copy)"),
    "unknown":          ("two heterozygous variants, phase unknown", "unknown"),
    "trans_split_site": ("two different ALT alleles at one site (GT 1/2)", "yes by genotype, weak read support"),
    "cis":              ("two heterozygous variants on the same haplotype", "no (cis)"),
}
ZYG_SCORE = {"hom": 2, "trans_split_site": 2, "unknown": 1, "X_hemizygous": 1, "single_het": 0.5, "cis": 0.5}

def vdesc(r):
    af = "AF ND" if pd.isna(r["grpmax_pct"]) else f"AF {r['grpmax_pct']}%"
    cadd = "CADD ND" if pd.isna(r["cadd"]) else f"CADD {r['cadd']}"
    return f"{r['variant']} ({r['chrom']}:{r['pos']}, {r['gt']}, {cadd}, {af})"

# ---------- One "unit" per gene: the variant, or the pair of variants, that supports the gene ----------
rows = []
for gene, g in v.groupby("gene"):
    assert len(g) <= 2, f"{gene} has more than 2 variants: define how to pair them"
    if len(g) == 1:
        r = g.iloc[0]
        set_cfg = "X_hemizygous" if r["chrom"] == "chrX" else ("hom" if r["gt"] == "1/1" else "single_het")
    else:
        set_cfg = PHASE[gene]
    members, ranked_cfg, note = g, set_cfg, ""
    ok = g[~g["gate_flag"]]
    # If one partner of a pair has no demonstrated effect or an invalid call, rank the gene on the other one
    if len(g) == 2 and set_cfg != "cis" and len(ok) == 1:
        members, ranked_cfg = ok, "single_het"
        bad = g[g["gate_flag"]].iloc[0]
        note = f"partner excluded: {bad['variant']} (no demonstrated effect on the gene or invalid call)"
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
        "u_call":   members["s_call"].min(),         # the weakest call limits the unit
        "u_zyg":    ZYG_SCORE[ranked_cfg],
        "u_pheno":  members["s_pheno"].iloc[0],      # gene-level, same for all its variants
        "u_mva":    members["s_mva"].iloc[0],
        "u_lit":    members["s_lit"].iloc[0],
        "any_gate_flag": bool(members["gate_flag"].any()),
    })
units = pd.DataFrame(rows)
assert len(units) == 29

# ---------- TABLE BEFORE THE RANKING: genes, specific variants and allele configuration ----------
pd.set_option("display.max_colwidth", 220); pd.set_option("display.width", 300)
table = units[["gene", "n_variants_in_set", "all_variants", "set_configuration",
               "different_alleles", "ranked_unit", "ranked_as", "note"]]
print(table.sort_values(["n_variants_in_set", "gene"], ascending=[False, True]).to_string(index=False))
table.to_csv("gene_table_before_ranking.csv", index=False)

# ---------- Weights (pending agree them with the team BEFORE looking at the totals) ----------
W = {"effect": 3, "pheno": 3, "mva": 3, "zyg": 2, "call": 2, "lit": 2, "freq": 1, "cadd": 1}

def unit_totals(w):
    return sum(w[k] * units[f"u_{k}"] for k in w) / (2 * sum(w.values())) * 100

units["total"] = unit_totals(W).round(1)
ranking = units.sort_values("total", ascending=False).reset_index(drop=True)
ranking["rank"] = ranking["total"].rank(ascending=False, method="min").astype(int)   # ties share the rank
ranking["in_top15"] = ranking["rank"] <= 15

print("\nFULL RANKING (all 29 genes):")
print(ranking[["rank", "gene", "ranked_unit", "ranked_as", "total", "any_gate_flag", "in_top15"]].to_string(index=False))

# ---------- Sensitivity: how stable is the ranking if the weights change? ----------
rng = np.random.default_rng(1)
N = 1000
res = []
for _ in range(N):
    w = {k: x * rng.uniform(0.5, 1.5) for k, x in W.items()}
    t = pd.Series(unit_totals(w).values, index=units["gene"])
    res.append(t.rank(ascending=False, method="min"))
R = pd.concat(res, axis=1)
sens = pd.DataFrame({"pct_in_top10": (R <= 10).mean(axis=1) * 100,
                     "pct_in_top15": (R <= 15).mean(axis=1) * 100,
                     "median_rank": R.median(axis=1)}).round(0)
print("\nSENSITIVITY:")
print(sens.sort_values("median_rank").to_string())

ranking.to_csv("gene_ranking.csv", index=False); sens.to_csv("sensitivity.csv")
v.to_csv("variant_matrix.csv", index=False)
