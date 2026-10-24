# Reference required for fresh normalization

Use GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna plus its samtools faidx index, matching the bcftools command recorded in the historical normalized VCF header. A generic different GRCh38 FASTA is not guaranteed equivalent. bcftools/htslib historical version: 1.19. This large reference is an external dependency, not a bundled input.
