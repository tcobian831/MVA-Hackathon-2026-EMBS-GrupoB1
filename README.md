# MVA Hackathon 2026 — EMBS Grupo B1

Variant prioritization for a child with confirmed mosaic variegated aneuploidy (MVA), integrating genomic annotations, clinical interpretation and four independent assessments.

## Analysis overview

The analysis has four components:

1. Annotation, normalization, CADD annotation and automated filtering.
2. Biological review and a recorded selection of **31 variants in 28 genes**.
3. A decision matrix completed independently by Quetzally, Daniel, Lucía and Tomás, followed by an equal-weight average of their definitive scores.
4. A subsequent search for additional variants in the full VEP dataset. Two additional candidates were retained for TACC2 and LRP4; this review does not alter the group ranking.

The original VCF contains **5,012,204 records**. The full VEP-annotated VCF contains **4,943,948 records**. These record counts refer to different files, and neither should be interpreted as a count of distinct candidate genes.

## Reproducibility and manual review

The biological selection is an explicit part of the analysis. `config/curated_variants.csv` records the accepted baseline alleles, and `docs/curation_history.md` documents the scientific review. Stage 07 checks these exact called alleles against the full VEP, including allele-specific gene annotations, genotype, depth and allele counts. It produces the 31-variant, 28-gene input to the matrix.

The 28 genes are therefore reproduced **with the recorded manual decisions as inputs**. Automated annotation and numeric filters alone do not determine this gene set. The stricter historical filter in stage 06 is a parallel candidate-review branch; its output is not the sole source of the manually selected variants. This distinction is necessary to reproduce the analysis as performed.

The recorded four-person mean is reproduced by stage 08. The additional TACC2 and LRP4 alleles are stored separately and are not included in this mean. Reference-based normalization requires the exact reference FASTA; the available validation uses the recorded normalized checkpoint for stages 03–06. Validation scope is detailed in `docs/reproducibility.md` and `docs/validation.json`.

## Repository structure

| Path | Contents |
|---|---|
| `scripts/01_annotation_filter.py` | Initial annotation/window filter. |
| `scripts/02_normalize.sh` | Chromosome renaming, reference checking, normalization and multiallelic splitting. |
| `scripts/03_add_cadd.py` | Exact-allele annotation using the fixed CADD v1.7 table. |
| `scripts/04_quality_cadd.py` | Quality/CADD filtering, with a separate missing-CADD branch. |
| `scripts/05_rarity.py` | Allele-matched population groupmax filtering. |
| `scripts/06_shortlist_filter.py` | Historical stricter shortlist branch. |
| `scripts/07_curate.py` | Validation and export of the manually curated 28-gene baseline. |
| `scripts/08_decision_matrix.py` | Final pipeline stage: individual rankings, group consensus and final gene/variant tables. |
| `scripts/matrix_engine.py`, `scripts/quality_helpers.py` | Shared implementation modules. |
| `tools/review_additional_variants.py` | Subsequent additional-variant review. |
| `config/curated_variants.csv` | The 31 baseline alleles and their scoring annotations. |
| `config/participant_manual_scores.json` | The four assessors' recorded manual scores. |
| `config/post_ranking_additions.csv` | LR2 and T18, retained after ranking. |
| `data/annotations/` | Fixed CADD annotation tables. |
| `results/` | Recorded individual rankings, consensus and variant-review tables. |
| `docs/` | Scientific curation, references and validation records. |

## Inputs and software

The genomic data require authorized access through the [official dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data). Raw genomic files and reference sequences are excluded from version control. The full annotated input location is:

`data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz`

Input SHA256 values are recorded in `config/input_sha256.json`. The normalization reference is `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`; details appear in `data/reference/README.md`.

Python dependencies are fixed in `requirements.txt`. Validation used Python 3.12. Normalization additionally requires bash, bcftools and tabix; the recorded historical command used version 1.19.

```bash
python -m venv .venv
# Linux/macOS activation: source .venv/bin/activate
# Windows PowerShell activation: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Pipeline stages

### 01 — Annotation filter

```bash
python scripts/01_annotation_filter.py data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz results/filter1.vcf
```

A record is retained if any of the following criteria holds:

- HIGH or MODERATE impact in a canonical transcript.
- A canonical-transcript consequence in `splice_region_variant`, `splice_donor_5th_base_variant`, `splice_donor_region_variant` or `splice_polypyrimidine_tract_variant`.
- Position inside a configured MVA-gene interval extended by **100 kb on each side**.
- A pathogenic/likely pathogenic or conflicting-pathogenicity ClinVar annotation.

Existing CSQ annotations are preserved and `F1_REASON` records the retention criteria. No quality, frequency or CADD cutoff is applied at this stage. The configured 12-gene window panel is retained as a historical filtering input; it is distinct from the final 28-gene set. Observed output: **26,167 records**.

### 02 — Normalization

```bash
bash scripts/02_normalize.sh results/filter1.vcf data/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna results/filter1.norm.vcf.gz
```

Chromosome names are standardized, reference alleles checked, and multiallelic records split with `--check-ref e -m -any --keep-sum AD --old-rec-tag ORIGINAL`. The recorded normalized output contains **26,375 records**. Splitting explains the increase in row count. `ORIGINAL` preserves the relationship to the unsplit source site.

### 03 — CADD annotation

```bash
python scripts/03_add_cadd.py results/filter1.norm.vcf.gz results/filter1.cadd.vcf
```

CADD PHRED values are joined by exact CHROM/POS/REF/ALT using the fixed v1.7 GRCh38 table. Missing values remain missing; unmatched input keys cause an error. All **26,375 records** are retained. The fixed table records the 13 October annotation checkpoint. The earlier 8 October web-VEP output was the first incorporation of CADD into project files; CADD was not present in the original full VEP.

### 04 — Quality and CADD

```bash
python scripts/04_quality_cadd.py results/filter1.cadd.vcf results/filter2.vcf results/filter2_no_cadd.vcf
```

Outside MVA windows, retention requires **DP >10, GQ >20 and CADD ≥20**. Inside the windows, **CADD ≥10** is required, while DP/GQ failures are flagged. Missing-CADD records satisfying the applicable quality criteria enter a separate branch. Exact matches to the 159 configured positive controls are protected.

Non-PASS calls, GQ 21–50 and heterozygous ALT balance outside 0.2–0.8 are flagged. Recorded outputs are **2,430 main records and 41 without CADD**.

### 05 — Population frequency

```bash
python scripts/05_rarity.py results/filter2.vcf results/filter3.vcf
python scripts/05_rarity.py results/filter2_no_cadd.vcf results/filter3_no_cadd.vcf
```

Groupmax is the largest available allele-matched gnomAD genomes AF among **AFR, AMR, EAS, MID, NFE and SAS**. AMI, ASJ, FIN and REMAINING do not enter this maximum. The rule is **groupmax <1% or unavailable**, with exact positive-control protection. Outputs are **235 main records and 11 without CADD**.

### 06 — Stricter historical shortlist branch

```bash
python scripts/06_shortlist_filter.py results/filter3.vcf --pending-vcf results/filter3_no_cadd.vcf --out-dir results/strict_shortlist
```

The recorded default rules are groupmax **<0.02% or unavailable**, CADD **>23 outside** or **>18 inside** MVA windows, and protection of BUB1B chr15:40209701 T>G. Missing CADD enters a pending list. Functional-priority labels are descriptive; they do not override these numeric defaults. The output is **60 selected and 8 pending variants**.

This branch supported candidate review. The final gene set also incorporated the broader cohort described in the manual-review record. Thresholds proposed during the project are distinguished from those implemented in this version.

### 07 — Manual curation and the 28-gene baseline

```bash
python scripts/07_curate.py data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz
```

The historical review combined a broader cohort of **2,429 CADD-scored and 68 unscored records**, gene-specific pair proposals, and biological evidence. Its partner search used groupmax **≤5% or unavailable** and CADD **>20 or unavailable**, without an additional DP/GQ filter. The review yielded 33 variants assigned to 29 genes. Removal of the two window-only MAD2L1BP assignments left the **31 baseline variants in 28 genes** used for ranking.

The accepted decisions are recorded in `config/curated_variants.csv`; exclusion rationales appear in `docs/curation_history.md`. Stage 07 validates every exact selected ALT as called in the sample, verifies its own gene annotation, and checks DP, GQ, GT and AD. The outputs are `results/curated_variants.csv` and `results/curation_validation.json`.

### 08 — Individual matrix and group consensus

Each member completed the same matrix with independent manual scores. The group ranking is the arithmetic mean of the four definitive scores for each gene. The recorded scores are provided in `config/participant_manual_scores.json`.

```bash
python scripts/08_decision_matrix.py --all-participants
```

Each execution creates a new directory under `results/runs/`, preserving existing recorded results. A specific new directory can be selected with `--output-dir`. The baseline CSV can be selected with `--variants`.

A separate assessor's scores can be supplied as JSON:

```bash
python scripts/08_decision_matrix.py --manual my_scores.json --name assessor
```

Each dictionary key is a gene symbol; its value is a list of three finite scores in [0,2], ordered **[MVA link, phenotype fit, literature evidence]**. Fractional values are allowed. The empty template is `config/manual_scores_template.json`.

| Component | Weight | Scope |
|---|---:|---|
| MVA link | 5.0 | Manual, gene-level |
| Phenotype fit | 4.5 | Manual, gene-level |
| Variant effect | 3.5 | Mean across included variants |
| Zygosity | 2.5 | Joint variant configuration |
| Population frequency | 2.0 | Mean across included variants |
| CADD | 1.5 | Mean across included variants |
| Literature evidence | 1.0 | Manual, gene-level |
| Call quality | 0.25 | Mean across included variants |

The baseline score is the normalized weighted mean, on [0,2]. Its percentage representation is the mean multiplied by 50. Zygosity scores are 0.5 for a single heterozygote or cis pair, 1 for unknown-phase heterozygous pairs, and 2 for homozygosity or a confirmed trans configuration; the engine records X-chromosome handling separately.

CADD is normalized as `clip(2*(CADD-20)/(49-20),0,2)`, with missing CADD assigned 1. Frequency uses AF in percent: `clip(2*(0.1-AF)/(0.1-0.00147),0,2)`, with missing frequency assigned 1. Effect overrides and call-quality flags are defined in the shared engine.

Sensitivity analysis uses **10,000 simulations**, seed 1. Each component weight is sampled independently and uniformly between 50% and 150% of its baseline, and each simulated score is renormalized to [0,2]. The engine reports the median score and percentages of simulations in the top 5, top 10, top 15 and rank 20 or later.

`definitive_score = 0.8*weighted_mean + 0.2*median_score`

The consensus percentage is `50*mean(definitive_score across the four members)`. These are prioritization scores, not calibrated causal probabilities or submission EPCR values.

Outputs include individual detailed rankings and sensitivity tables, `consensus.csv` for all 28 genes, `top10_genes.csv` and `top10_variants.csv`. The last table identifies exact GRCh38 alleles, genotype, phase status and baseline versus post-ranking inclusion.

## Final group ranking

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


The final top ten are **BUB1B, TACC2, SMC4, SVIL, LZTR1, CTU2, DIDO1, PRMT6, LRP4 and TRDMT1**. Exact allele assignments are listed in [the final candidate table](results/top10_variants.csv).

## Subsequent additional-variant review

```bash
python tools/review_additional_variants.py data/raw/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf.gz
```

The full VEP was revisited for the historical 15-gene search panel using all matching transcripts and the original sample genotypes. The exploratory criteria were **DP ≥5, GQ ≥15, groupmax ≤5% or unavailable, and CADD ≥4 or unavailable**. The output records each called allele, its annotations and its review decision.

Two additional variants were retained: **TACC2 T18, chr10:122194203 CA>C**, and **LRP4 LR2, chr11:46878979 T>C**. They are recorded in `config/post_ranking_additions.csv`. Their phase relative to the respective baseline variant is unknown. The gene order and consensus scores remain those obtained before this review.

| Gene | GRCh38 allele | GT | Inclusion | Phase |
|---|---|---|---|---|
| BUB1B | chr15:40209701 T>G | `0/1` | Ranking baseline | unknown between selected variants |
| BUB1B | chr15:40220612 T>G | `0/1` | Ranking baseline | unknown between selected variants |
| TACC2 | chr10:122087185 C>G | `0/1` | Ranking baseline | unknown between selected variants |
| TACC2 | chr10:122194203 CA>C | `0/1` | Subsequent review | unknown between selected variants |
| SMC4 | chr3:160411998 G>A | `0/1` | Ranking baseline | single heterozygous candidate |
| SVIL | chr10:29463555 C>T | `0/1` | Ranking baseline | single heterozygous candidate |
| LZTR1 | chr22:20996720 C>G | `0/1` | Ranking baseline | single heterozygous candidate |
| CTU2 | chr16:88714226 AGGTGTG>A | `1/1` | Ranking baseline | same variant on both copies (GT 1/1) |
| DIDO1 | chr20:62926430 ACCGCAGG>A | `0/1` | Ranking baseline | single heterozygous candidate |
| PRMT6 | chr1:107056749 G>C | `0/1` | Ranking baseline | single heterozygous candidate |
| LRP4 | chr11:46874935 C>T | `0/1` | Ranking baseline | unknown between selected variants |
| LRP4 | chr11:46878979 T>C | `0/1` | Subsequent review | unknown between selected variants |
| TRDMT1 | chr10:17149087 G>A | `0/1` | Ranking baseline | single heterozygous candidate |

The final top-ten inventory contains **13 distinct alleles**, including the two later additions. A homozygous variant occupies one allele-record row while being called on both gene copies. Multiple heterozygous variants are not described as confirmed trans without phase evidence. Reasons for non-inclusion of the other additional variants are provided in [additional_variants.md](docs/additional_variants.md), with the per-allele audit in `results/expanded_review.csv`.

## Validation and data access

Stage 01 was checked against all 4,943,948 VEP records. CADD annotation and stages 04–06 match the recorded normalized checkpoints. The curated baseline was validated against the source alleles, and the four matrix outputs and consensus arithmetic were checked independently. Normalization with the original reference remains a separately documented validation dependency. Details are in [reproducibility.md](docs/reproducibility.md).

The genomic input is access-restricted and is not distributed through this repository. Data use and retention are governed by the [challenge rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py). Scientific references appear in [curation_bibliography.md](docs/curation_bibliography.md).
