#!/usr/bin/env bash
set -euo pipefail
# Match the historical normalization: split ALT alleles, preserve ORIGINAL and AD sum.
# Input 1: Filter 1 VCF; input 2: exact GRCh38 FASTA; input 3: normalized output.
[[ $# == 3 ]] || { echo 'Usage: bash scripts/02_normalize.sh INPUT.vcf REF.fa OUTPUT.vcf.gz' >&2; exit 1; }
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
bcftools query -f '%CHROM\n' "$1" | sort -u | awk '{c=$1; if(c !~ /^chr/) c="chr" c; print $1 "\t" c}' > "$work/chromosomes.tsv"
bcftools annotate --rename-chrs "$work/chromosomes.tsv" -Ou "$1" | bcftools norm -f "$2" --check-ref e -m -any --keep-sum AD --old-rec-tag ORIGINAL -Oz -o "$3"
tabix -f -p vcf "$3"
