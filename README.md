# MVA Hackathon 2026 — EMBS Grupo B1

A modular reconstruction of the team's variant-prioritization workflow for Track 1.

## What this repository reproduces

This repository separates **automated filtering**, **documented human curation**, **a subsequent genome-wide partner search**, and **decision-matrix ranking**. Human biological review is recorded as an explicit input, rather than disguised as an automatic filter. Four assessors' recorded manual scores are reused without modification to recompute comparable rankings on the baseline variant set.

Group ranking input: **31 variant entries in 28 genes**, after excluding MAD2L1BP and its two window insertions. After the group ranking was established, two additional variants were retained for review: LR2 and T18. These later additions do not change the ranking input or consensus reported here. The combined review inventory has 33 variants; the ranking baseline has 31. Stable IDs are not record counts.

The original caller VCF has **5,012,204 records**; the complete VEP-annotated input has **4,943,948 records**. They are different files. The restricted VEP file is supplied in the separate **private inputs archive**, not in this public-repository ZIP. Access the official dataset through https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data. Do not upload the private archive to GitHub.

## Directory guide

| Path | Purpose |
|---|---|
| `scripts/01_annotation_filter.py` | One canonical implementation of the initial annotation/window filter. |
| `scripts/02_normalize.sh` | Reference-checked chromosome renaming, multiallelic splitting and left normalization. |
| `scripts/03_add_cadd.py` | Exact-allele replay of the frozen CADD v1.7 annotation checkpoint. |
| `scripts/04_quality_cadd.py` | Daniel's quality/CADD rules and missing-CADD branch. |
| `scripts/05_rarity.py` | Allele-matched six-population groupmax rarity filter. |
| `scripts/06_shortlist_filter.py` | Historical stricter shortlist branch; does not automatically produce the manually selected genes. |
| `scripts/07_curate.py` | Replay the documented human selection, validating positions, gene annotations, DP and GQ against the full VEP. |
| `scripts/08_decision_matrix.py` | **Final analytical stage**: run the shared matrix for one or all four assessors and calculate their consensus. |
| `tools/review_additional_variants.py` | Separate post-ranking review; never feeds additional variants back into the group consensus. |
| `scripts/matrix_engine.py` | Shared scoring implementation used by stage 08; no second scoring algorithm. |
| `scripts/quality_helpers.py` | Shared access to stage 04's functions; no duplicated rules. |
| `config/` | Gene windows, positive controls, explicit curated variants, manual scores, input hashes. |
| `data/annotations/` | Frozen exact-allele annotation snapshots needed to reproduce scores without changing online databases. |
| `results/` | Current four-assessor results, consensus and review audit; clearly dated reconstruction, not original participant reruns. |
| `docs/` | Selection rationale, numerical validation, migration decisions and limitations. |

## Setup and inputs

Python 3.12 was used for validation. Install the pinned Python dependencies:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Normalization additionally needs bash, bcftools and tabix (historical version 1.19), plus the exact reference FASTA named in `data/reference/README.md`. Use Linux, macOS or WSL for that stage. The reference is an external dependency; it is not bundled.

Place the authorized VEP input at `data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz`. Verify SHA256 against `config/input_sha256.json` before a fresh run. The separate private archive also supplies the historical normalized checkpoint, enabling a verified restart after normalization.

## Stages, in execution order

### 01 — Annotation-based retention

```bash
python scripts/01_annotation_filter.py data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz results/filter1.vcf
```

Input: complete VEP VCF. Output: VCF with all existing CSQ annotations preserved and `F1_REASON` added. The conditions are joined by **OR**:

* HIGH or MODERATE impact on a canonical transcript.
* One of `splice_region_variant`, `splice_donor_5th_base_variant`, `splice_donor_region_variant`, `splice_polypyrimidine_tract_variant` on a canonical transcript.
* Position inside a configured MVA gene interval extended by **100,000 bp on both sides**.
* ClinVar pathogenic/likely pathogenic, or conflicting classifications of pathogenicity, in a CSQ annotation.

No DP, GQ, CADD or frequency filtering occurs here. The supplied 12-gene window table includes MAD2L1BP for historical reproduction; keeping that historical rescue rule does not reintroduce MAD2L1BP into the final 28-gene matrix. **Verified: 4,943,948 → 26,167 VCF records.**

### 02 — Normalize and split alleles

```bash
bash scripts/02_normalize.sh results/filter1.vcf data/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna results/filter1.norm.vcf.gz
```

Input: stage 01 plus exact reference. Output: normalized compressed VCF and tabix index. Chromosomes become `chr`-prefixed. The command uses `--check-ref e -m -any --keep-sum AD --old-rec-tag ORIGINAL`. The ORIGINAL tag is essential for matching VEP allele annotations after splitting. The normalized checkpoint contains **26,375 records**, not 26,167: splitting can increase row count.

Fresh normalization was not rerun in this environment because the exact FASTA and bcftools were unavailable. Downstream reproduction was tested against the historical normalized checkpoint; this limit is recorded in `docs/validation.json`.

### 03 — Incorporate CADD

```bash
python scripts/03_add_cadd.py results/filter1.norm.vcf.gz results/filter1.cadd.vcf
```

Input: normalized VCF and bundled frozen annotation table. Output: same 26,375 records with `INFO/CADD_PHRED` when available. Matching is exact CHROM/POS/REF/ALT. The script fails on keys outside the snapshot rather than silently annotating them as missing.

Provenance: CADD v1.7 GRCh38, precomputed SNVs and gnomAD v4.0 indels, recorded in Daniel's 13 October checkpoint. An earlier 8 October web-VEP output was the first historical incorporation of CADD into project files. This consolidated branch uses the later complete checkpoint; it does not claim CADD was present in the original full VEP.

### 04 — Quality and CADD filter

```bash
python scripts/04_quality_cadd.py results/filter1.cadd.vcf results/filter2.vcf results/filter2_no_cadd.vcf
```

Outside MVA windows: **DP >10, GQ >20, CADD PHRED ≥20**. Inside windows: **CADD ≥10**; DP/GQ become flags rather than exclusion rules. Missing CADD is routed to a separate file after the relevant quality decision. The 159 known positive controls are exact-allele overrides; one is present in the retained output. Non-PASS, intermediate GQ 21–50 and heterozygous ALT balance outside 0.2–0.8 are flagged, not independent hard exclusions. Verified: **2,430 scored/override records and 41 missing-CADD records**.

### 05 — Population rarity

```bash
python scripts/05_rarity.py results/filter2.vcf results/filter3.vcf
python scripts/05_rarity.py results/filter2_no_cadd.vcf results/filter3_no_cadd.vcf
```

Both branches use allele-matched CSQ frequencies. Groupmax is the maximum available gnomAD genomes AF among **AFR, AMR, EAS, MID, NFE, SAS**. AMI, ASJ, FIN and REMAINING do not enter this maximum. Keep groupmax **<0.01 (1%)** or missing; exact positive controls are protected. Global AF and excluded-group AF can produce flags, but do not replace this rule. Verified outputs: **235 + 11 records**.

### 06 — Historical stricter shortlist branch

```bash
python scripts/06_shortlist_filter.py results/filter3.vcf --pending-vcf results/filter3_no_cadd.vcf --out-dir results/strict_shortlist
```

Default observed historical rules: protect BUB1B chr15:40209701 T>G; otherwise groupmax **<0.0002 (0.02%)** or missing; available CADD **strictly >23 outside** and **strictly >18 inside** windows; missing CADD remains pending. Functional high/intermediate/low labels are reported but do not themselves override those default numeric rules. Verified: **60 selected + 8 pending** from 246 inputs.

These are the actual defaults in the recovered script, not the later proposed thresholds ≥25/≥20. A proposal is not substituted for an observed implementation.

### 07 — Explicit manual curation and current ranking input

```bash
python scripts/07_curate.py data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz
```

Input: complete VEP plus `config/curated_variants.csv`. Output: `results/curated_variants.csv` and `curation_validation.json`. Verified: all **31 selected positions**, their gene-symbol associations and DP/GQ are present in the full VEP; resulting set contains **28 genes**.

**This is a manual-selection replay, not an automatic filter that discovers 28 genes.** The historical review assembled the original 29 genes from a broader Lucía cohort (2,429 scored + 68 without CADD) and Daniel's pair proposals, not solely from the strict stage 06 output. Its all-target partner search used groupmax ≤5% or absent and CADD >20 or absent, without another DP/GQ cutoff. See `docs/curation_history.md` for the documented merge and gene-level exclusions. The broader cohort's first-principles re-execution is not claimed here; the curated manifest makes its resulting decisions reproducible and inspectable.

The correction before ranking removes both MAD2L1BP window insertions. LR2 and T18 are retained separately in a later review. `docs/curation_decisions.csv` records these actions. Do not apply the strict stage 06 cutoffs to stage 07 and expect the same manually reviewed list: that would silently change the analysis.

### 08 — Decision matrix (last pipeline script)

For a fresh output directory:

```bash
python scripts/08_decision_matrix.py --all-participants
```

This runs the same matrix for Quetzally, Daniel, Lucía and Tomás using their three recorded manual scores per gene. Their original Drive outputs show that each participant ran the earlier matrix; the current four-person outputs here are **centrally recomputed with their unchanged manual inputs**, not falsely attributed as four new participant executions. Each output directory must be absent before rerunning; move old outputs aside rather than mixing versions.

For a new assessor, copy `config/manual_scores_template.json`, replace every null with a score, then run:

```bash
python scripts/08_decision_matrix.py --manual my_scores.json --name myname
```

The dictionary key is the exact gene symbol. Each value is a list of three finite numeric values in [0,2], in the order **[MVA link, phenotype fit, literature evidence]**. Fractional scores are permitted. Missing values are never silently replaced by zero.

Weights: MVA 5, phenotype 4.5, effect 3.5, zygosity 2.5, frequency 2, CADD 1.5, literature 1, call quality 0.25. Gene effect/CADD/frequency/call scores average all included variants. Manual fields remain gene-level; zygosity is a joint configuration. Therefore adding a weak second variant can reduce total score despite increasing zygosity.

CADD matrix score: clip[2(CADD−20)/(49−20),0,2], missing →1. Frequency matrix score uses frequencies in PERCENT: clip[2(0.1−AF)/(0.1−0.00147),0,2], missing →1. These normalization anchors remain unchanged in the current set. Call quality uses the engine's depth/GQ/allele-balance flags. Consequence overrides and their stable IDs are retained in the single engine implementation.

The baseline weighted mean is Σ(weight × score)/Σweight on [0,2]; `total` is that mean ×50. Each of **10,000 simulations** independently varies each weight uniformly from 50% to 150% of baseline and renormalizes the weighted sum. Seed=1. The median score and top5/top10/top15/rank≥20 percentages are reported; ties share the minimum rank. Definitive score = **0.8 × baseline mean + 0.2 × simulated median**, on [0,2]. These are prioritization scores, not causal probabilities or submission EPCR.

Final stage outputs include each assessor's `gene_ranking_summary.csv`, detailed `gene_ranking.csv`, `variant_matrix.csv`, `sensitivity.csv` and a validated four-assessor `results/consensus.csv`. Consensus is the equal-weight mean of the four definitive scores, converted to percent by multiplying by 50. Rounding occurs only for display.

## Group consensus and final top ten

| Rank | Gene | Mean definitive score (%) |
|---:|---|---:|
| 1 | BUB1B | 81.2583 |
| 2 | TACC2 | 53.8819 |
| 3 | SMC4 | 52.4485 |
| 4 | SVIL | 50.1593 |
| 5 | LZTR1 | 49.6770 |
| 6 | CTU2 | 46.5415 |
| 7 | DIDO1 | 45.8408 |
| 8 | PRMT6 | 44.6462 |
| 9 | LRP4 | 38.8009 |
| 10 | TRDMT1 | 37.5330 |
| 11 | ROCK1 | 37.2942 |
| 12 | CITED2 | 36.5069 |
| 13 | IGBP1 | 36.0472 |
| 14 | PIK3C3 | 35.1607 |
| 15 | THUMPD3 | 34.5290 |
| 16 | ACLY | 33.3806 |
| 17 | ZNF418 | 31.8201 |
| 18 | ZNF197 | 31.4685 |
| 19 | FRZB | 31.2763 |
| 20 | SIM2 | 31.0414 |
| 21 | ADAMTS1 | 29.9651 |
| 22 | SOCS5 | 29.9563 |
| 23 | CNOT10 | 29.2460 |
| 24 | SLC8A1 | 27.7550 |
| 25 | TSSK6 | 24.1948 |
| 26 | MEP1B | 23.8849 |
| 27 | CD82 | 22.9778 |
| 28 | ASAP3 | 20.6400 |


The group final top ten, obtained by averaging the four assessors' definitive scores, are **BUB1B, TACC2, SMC4, SVIL, LZTR1, CTU2, DIDO1, PRMT6, LRP4, TRDMT1**. LRP4 remains ninth; ROCK1 is eleventh. Subsequent additional-variant review does not modify this ranking.

## Subsequent review — additional partner search after the group ranking

```bash
python tools/review_additional_variants.py data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz
```

Input: full VEP, all matching transcripts and the frozen supplemental CADD snapshot. Output: one audit row per called ALT/gene in the historical 15-gene search panel. It retains original unsplit GT, AD and local phase metadata.

The latest exploratory screen is **DP ≥5, GQ ≥15, groupmax ≤5% or absent, CADD ≥4 or absent**. It is distinct from all previous filtering thresholds. Passing this screen makes an allele eligible for review, not pathogenic. Missing CADD means no value was recovered in the frozen sources; no damage is inferred.

Only two **additional** variants were retained for later review: **LRP4 LR2 chr11:46878979 T>C** and **TACC2 T18 chr10:122194203 CA>C**. Each gene's proposed pair has unknown phase. A separate experimental re-ranking assigned zygosity 1, but that experiment is not the group ranking reported above. Neither is represented as confirmed compound heterozygosity. The audit distinguishes pre-existing matrix variants from genuinely new additions. Other candidates fail quality/frequency/CADD, or lack a sufficiently supported additional gene-specific mechanism. Detailed considerations are in `docs/additional_variants.md`; they are not all claimed benign.

The chronological sequence is: select the 28 genes and 31 baseline variants → each assessor runs the matrix with personal manual scores → average the four definitive scores → obtain the group top ten → reopen the full VEP for additional variants. The additional review does not overwrite the ranked set or automatically change the top ten. The two additions are stored in `config/post_ranking_additions.csv`; they are not read by the matrix engine.

## Validation and publication

Stage 01 was rerun against all 4,943,948 VEP records. CADD replay and stages 04–06 were rerun on the complete normalized checkpoint. Stage 07 checked the full input; the separate review tool rescans that full input. All four matrix runs completed 10,000 simulations; their weighted means and definitive-score arithmetic were checked independently. Fresh reference normalization remains an explicit external validation limit. See docs/validation.json.

This is a replacement repository tree. Unzip, inspect it, and replace the old working tree deliberately; uploading it alongside old files would preserve the existing duplicates. No remote repository or Drive file was changed. Do not upload the separate private-input ZIP, raw VCF, reference or read data. GitHub blocks regular Git files larger than 100 MiB; the genomic file is also access-restricted.

Official sources: [challenge](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), [rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py), [GitHub size limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github). Participant submissions are subject to the challenge CC-BY terms; raw data are not covered by an unrestricted redistribution permission. Follow the challenge retention/deletion rules.
