# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 10:47:58 2026

@author: mavec
"""

import gzip, re, os
import pandas as pd
# NOTE: every value of flag_freq_missing is True (and freq_max is empty) because the
# VEP VCF used here does not include population frequency data (no gnomAD / 1000 Genomes
# fields in CSQ). Once Quetzally re-runs VEP with the frequency options enabled,
# update VCF_PATH to the new file and re-run the script; the frequency fields will be
# detected automatically.

# ================= CONFIGURATION =================
VCF_PATH    = r"C:\Users\mavec\OneDrive\Documentos\MVA_HCKATON\filter1_2026-10-08-00-51_VEP.vcf\filter1_2026-10-08-00-51_VEP.vcf"  # TODO: swap for the new VEP file with frequencies (flag_freq_missing is all True until then)
FILTER1_VCF = r"C:\Users\mavec\OneDrive\Documentos\MVA_HCKATON\filter1.norm.vcf\filter1.norm.vcf"                             # filter 1 output: DP, GQ, F1_REASON
GENES_PATH  = r"C:\Users\mavec\OneDrive\Documentos\MVA_HCKATON\mva_genes (1).tsv"
OUT_DIR     = r"C:\Users\mavec\OneDrive\Documentos\MVA_HCKATON\salida_filtros"              # new folder, so earlier outputs are not overwritten

# --- Thresholds (Tomas proposal) ---
CADD_OUT = 20          # outside gene windows: CADD PHRED >= 20
CADD_IN  = 10          # inside gene windows: CADD PHRED >= 10
DP_MIN   = 10          # DP > 10 (outside windows only)
GQ_MIN   = 20          # GQ > 20 (outside windows only)
GQ_INTERMEDIATE_MAX = 50   # flag calls with GQ between 21 and 50
WINDOW   = 100_000     # bases on each side of the gene
SAMPLE_INDEX = 0       # first sample in the filter 1 VCF

# --- Population frequency ---
FREQ_MAX    = None     # None = do NOT filter, only add the column and flag. E.g. 0.01 to filter
FREQ_FIELDS = None     # None = autodetect in CSQ, or give a list, e.g. ["gnomADg_AF"]
FREQ_TEST   = [0.05, 0.01, 0.001, 0.0001]   # thresholds compared in the summary table
# =================================================

os.makedirs(OUT_DIR, exist_ok=True)

def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def norm_chr(c):
    # Make chromosome names comparable: "chr1" -> "1", "M" -> "MT"
    c = str(c)
    if c.lower().startswith("chr"):
        c = c[3:]
    return "MT" if c == "M" else c

def to_num(x):
    # Convert to float, or return None if it is not a number
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

# ---------- 1) Gene windows ----------
with open_any(GENES_PATH) as f:
    lines = f.read().splitlines()

# Find the header row of the table (the file may have text lines before it)
KEYS = ("gene", "symbol", "chr", "chrom", "start", "end")
hdr_i = None
for i, l in enumerate(lines):
    parts = [p.strip().lstrip("#").lower() for p in l.split("\t")]
    if len(parts) >= 4 and sum(any(k in p for k in KEYS) for p in parts) >= 3:
        hdr_i = i
        break
if hdr_i is None:
    raise SystemExit("Could not find the header row in mva_genes.tsv.")

genes = pd.read_csv(GENES_PATH, sep="\t", skiprows=hdr_i, dtype=str)
genes.columns = [c.strip().lstrip("#") for c in genes.columns]
cols = {c.lower(): c for c in genes.columns}

def pick(*names):
    # Return the first matching column name (case-insensitive)
    for n in names:
        if n in cols:
            return cols[n]
    raise SystemExit(f"Column not found {names}. Available columns: {list(genes.columns)}")

g_gene  = pick("gene", "gene_name", "symbol", "gen")
g_chr   = pick("chr", "chrom", "chromosome", "seqname")
g_start = pick("start", "gene_start", "start_position", "start_bp")
g_end   = pick("end", "gene_end", "end_position", "end_bp")

# Build one window per gene: gene interval +/- WINDOW bases
genes = genes.dropna(subset=[g_chr, g_start, g_end])
windows = {}
for _, r in genes.iterrows():
    try:
        s, e = int(float(r[g_start])), int(float(r[g_end]))
    except ValueError:
        continue
    windows.setdefault(norm_chr(r[g_chr]), []).append(
        (s - WINDOW, e + WINDOW, str(r[g_gene])))

print("Genes loaded:", sum(len(v) for v in windows.values()))

# ---------- 2) Filter 1 VCF: DP, GQ, F1_REASON and original lines ----------
qc = {}                                  # (chrom, pos, ref, alt) -> (dp, gq, f1_reason, line index)
f1_header, f1_raw, f1_samples = [], [], []
with open_any(FILTER1_VCF) as f:
    for line in f:
        if line.startswith("##"):
            f1_header.append(line); continue
        if line.startswith("#"):
            f1_header.append(line)
            f1_samples = line.rstrip("\n").split("\t")[9:]
            continue
        c = line.rstrip("\n").split("\t")
        info = {}
        for item in c[7].split(";"):
            k, _, v = item.partition("=")
            info[k] = v
        dp = gq = None
        if len(c) > 9 + SAMPLE_INDEX:
            # Genotype fields (DP, GQ) of the selected sample
            fmt = dict(zip(c[8].split(":"), c[9 + SAMPLE_INDEX].split(":")))
            dp, gq = to_num(fmt.get("DP")), to_num(fmt.get("GQ"))
        qc[(norm_chr(c[0]), int(c[1]), c[3], c[4])] = (
            dp, gq, info.get("F1_REASON", ""), len(f1_raw))
        f1_raw.append(line)

print("Variants in filter 1 VCF:", len(qc))
print("Samples in filter 1 VCF:", f1_samples, "| using index", SAMPLE_INDEX)

# ---------- 3) VEP VCF: CADD PHRED, frequency and gene ----------
csq_fields, freq_fields = None, []
rows = []
with open_any(VCF_PATH) as f:
    for line in f:
        if line.startswith("##"):
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format: ([^\">]+)", line)
                csq_fields = m.group(1).strip().split("|")
            continue
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                if csq_fields is None:
                    raise SystemExit("CSQ header not found.")
                # Frequency fields: those whose name ends in AF (gnomADg_AF, AFR_AF, MAX_AF...)
                freq_fields = FREQ_FIELDS if FREQ_FIELDS is not None else \
                    [x for x in csq_fields if re.search(r"(^|_)AF$", x, re.I)]
                print("Frequency fields used:", freq_fields)
                if not freq_fields:
                    print("WARNING: no frequency fields found in CSQ. "
                          "The column will be empty and every variant will be flagged as missing frequency.")
            continue

        c = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = c[0], int(c[1]), c[3], c[4]

        info = {}
        for item in c[7].split(";"):
            k, _, v = item.partition("=")
            info[k] = v

        # Each variant has one CSQ entry per transcript: take the max CADD and max frequency
        cadd_vals, freq_vals, gene = [], [], ""
        if "CSQ" in info:
            for tx in info["CSQ"].split(","):
                d = dict(zip(csq_fields, tx.split("|")))
                v = to_num(d.get("CADD_PHRED"))              # PHRED, not RAW
                if v is not None:
                    cadd_vals.append(v)
                for ff in freq_fields:
                    for piece in str(d.get(ff, "")).split("&"):   # "&" separates co-located values
                        fv = to_num(piece)
                        if fv is not None:
                            freq_vals.append(fv)
                if not gene:
                    gene = d.get("SYMBOL", "")
        cadd = max(cadd_vals) if cadd_vals else None
        freq = max(freq_vals) if freq_vals else None         # highest frequency observed

        # Match with filter 1 to get DP/GQ
        dp, gq, f1r, f1_idx = qc.get((norm_chr(chrom), pos, ref, alt),
                                     (None, None, "", None))
        # Gene windows that contain this position
        w_genes = [g for (s, e, g) in windows.get(norm_chr(chrom), []) if s <= pos <= e]

        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                     "gene": gene, "f1_reason": f1r,
                     "dp": dp, "gq": gq, "cadd_phred": cadd, "freq_max": freq,
                     "in_window": bool(w_genes),
                     "window_genes": ",".join(sorted(set(w_genes))),
                     "f1_idx": f1_idx})

df = pd.DataFrame(rows)
for col in ("dp", "gq", "cadd_phred", "freq_max"):
    df[col] = pd.to_numeric(df[col], errors="coerce")

print("\nVariants in VEP VCF:", len(df))
print("With DP (matched to filter 1):", int(df["dp"].notna().sum()), "of", len(df))
print("With CADD PHRED:", int(df["cadd_phred"].notna().sum()), "of", len(df))
print("With population frequency:", int(df["freq_max"].notna().sum()), "of", len(df))
if df["dp"].notna().sum() < 0.8 * len(df):
    print("WARNING: few variants matched the filter 1 VCF. Check chrom/pos/ref/alt formatting.")

# ---------- 4) Flags ----------
df["pass_qc"] = (df["dp"] > DP_MIN) & (df["gq"] > GQ_MIN)
df["flag_low_quality"]    = ~df["pass_qc"]                                    # fails DP > 10 or GQ > 20
df["flag_gq_intermediate"] = (df["gq"] > GQ_MIN) & (df["gq"] <= GQ_INTERMEDIATE_MAX)
df["flag_freq_missing"]   = df["freq_max"].isna()                             # no frequency data

# ---------- 5) Selection (Tomás's proposal) ----------
has_cadd = df["cadd_phred"].notna()
cond_out = (~df["in_window"]) & df["pass_qc"] & (df["cadd_phred"] >= CADD_OUT)
cond_in  = df["in_window"] & (df["cadd_phred"] >= CADD_IN)   # windows are exempt from DP/GQ
sel_base = has_cadd & (cond_out | cond_in)                    # variants without CADD are excluded

if FREQ_MAX is None:
    freq_ok = pd.Series(True, index=df.index)                 # no frequency filter
else:
    # Variants without frequency data are kept (and flagged)
    freq_ok = df["flag_freq_missing"] | (df["freq_max"] <= FREQ_MAX)

sel = sel_base & freq_ok

# Window variants kept only because of their location (would not pass the criteria used outside windows)
df["kept_by_window"] = df["in_window"] & (df["flag_low_quality"] | (df["cadd_phred"] < CADD_OUT))

# ---------- 6) Summary ----------
print("\nWithout CADD score (separate list, excluded):", int((~has_cadd).sum()),
      "| of which inside windows:", int(((~has_cadd) & df["in_window"]).sum()))
print("Inside gene windows:", int(df["in_window"].sum()))
print("Pass DP/GQ:", int(df["pass_qc"].sum()))

print("\n=== Selection (before any frequency filter) ===")
print("Total:", int(sel_base.sum()),
      "| outside windows:", int((sel_base & ~df["in_window"]).sum()),
      "| inside windows:", int((sel_base & df["in_window"]).sum()),
      "| inside windows with low quality:", int((sel_base & df["in_window"] & df["flag_low_quality"]).sum()))

# How many variants would remain with different frequency thresholds
print("\n=== Effect of the frequency threshold (to help choose one) ===")
rows_t = []
for t in [None] + FREQ_TEST:
    ok = pd.Series(True, index=df.index) if t is None else \
         (df["flag_freq_missing"] | (df["freq_max"] <= t))
    s = sel_base & ok
    rows_t.append({"freq_threshold": "no filter" if t is None else t,
                   "retained": int(s.sum()),
                   "outside_windows": int((s & ~df["in_window"]).sum()),
                   "inside_windows": int((s & df["in_window"]).sum()),
                   "missing_freq": int((s & df["flag_freq_missing"]).sum())})
print(pd.DataFrame(rows_t).to_string(index=False))
print("\nFrequency distribution (variants with data):")
print(df["freq_max"].describe(percentiles=[.5, .9, .95]).round(5))

# ---------- 7) Output files ----------
out_cols = ["chrom", "pos", "ref", "alt", "gene", "f1_reason", "dp", "gq",
            "cadd_phred", "freq_max", "in_window", "window_genes",
            "flag_low_quality", "flag_gq_intermediate", "flag_freq_missing",
            "kept_by_window"]
suffix = "" if FREQ_MAX is None else f"_freq{FREQ_MAX}"

# Final candidates, sorted by CADD (highest first)
final = df[sel]
final.sort_values("cadd_phred", ascending=False)[out_cols].to_csv(
    os.path.join(OUT_DIR, f"final_candidates{suffix}.tsv"), sep="\t", index=False)

# Same variants as a VCF, using the original filter 1 lines and headers
idxs = sorted(int(i) for i in final["f1_idx"].dropna())
with open(os.path.join(OUT_DIR, f"final_candidates{suffix}.vcf"), "w", encoding="utf-8") as out:
    out.writelines(f1_header)
    out.writelines(f1_raw[i] for i in idxs)

# Variants without CADD score: separate list for possible manual review
df[~has_cadd][[c for c in out_cols if c != "cadd_phred"]].to_csv(
    os.path.join(OUT_DIR, "no_CADD_for_review.tsv"), sep="\t", index=False)

# If a frequency filter was applied, also save what it removed (for auditing)
if FREQ_MAX is not None:
    df[sel_base & ~freq_ok].sort_values("cadd_phred", ascending=False)[out_cols].to_csv(
        os.path.join(OUT_DIR, f"removed_by_frequency{suffix}.tsv"), sep="\t", index=False)

print("\nFiles saved in:", OUT_DIR)