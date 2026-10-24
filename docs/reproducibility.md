# Reproducibility and validation

## Computational stages

The complete VEP input contains 4,943,948 records. Stage 01 retains 26,167 records. The recorded reference-normalized output contains 26,375 split records. CADD annotation preserves that count. Quality/CADD filtering produces 2,430 main and 41 missing-CADD records; groupmax filtering produces 235 and 11. The stricter shortlist branch produces 60 selected and 8 pending records.

The stage 02 command is documented with bcftools/htslib 1.19 and the exact reference filename. Current computational verification resumes from the recorded normalized file; it does not constitute a fresh reference-normalization validation. The source input hashes are recorded separately.

## Human curation as a required analysis input

The final gene selection combined biological review, a broader cohort of 2,429 scored plus 68 unscored records, and gene-specific genotype/pair proposals. This review was not equivalent to selecting only the strict 60-variant branch. Its accepted decisions are the 31 baseline entries in config/curated_variants.csv after removal of the two MAD2L1BP window-only assignments.

Stage 07 verifies the selected alleles directly in the full VEP, including the called ALT, matching CSQ gene, DP, GQ, GT and AD. It exports 31 variant records in 28 genes. This reproduces the recorded human selection; it does not model biological judgement as a numeric filter or independently infer a new shortlist. A first-principles rerun of the broader historical cohort is not part of the available computational verification.

## Four-member ranking

The decision matrix uses the four members' recorded MVA-link, phenotype-fit and literature scores. The 28 gene sets must agree across members. Each run performs 10,000 weight simulations with seed 1. Independently checked identities are:

- Weighted mean = sum(component score × weight) / sum(weights).
- Percentage score = 50 × weighted mean.
- Definitive score = 0.8 × weighted mean + 0.2 × median simulated score.
- Consensus percentage = 50 × the equal-weight mean of the four definitive scores.

The resulting ordering is BUB1B, TACC2, SMC4, SVIL, LZTR1, CTU2, DIDO1, PRMT6, LRP4 and TRDMT1 for the first ten positions. LR2 and T18 do not enter the ranking baseline.

## Final gene-to-variant mapping

Stage 08 exports both the ranked gene table and an exact-allele table for those genes. The latter includes the baseline variants and two separately labelled additions retained during the subsequent review. The final top-ten gene list has 13 distinct variant records. Unknown phase remains unknown; no inferred trans relationship is introduced.

## Verification scope

The fixed-input ranking and its arithmetic are reproducible. The validity of the recorded biological selection, the clinical pathogenicity of each candidate, and the phase of unphased heterozygous pairs are separate questions. The computational validation does not establish these as experimentally confirmed.
