#!/usr/bin/env python3
"""Validate the manually curated baseline against exact called alleles in the full VEP VCF."""
import argparse
import csv
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('full_vep')
parser.add_argument('--out-dir', type=Path, default=ROOT / 'results')
args = parser.parse_args()
with (ROOT / 'config/curated_variants.csv').open() as stream:
    reader = csv.DictReader(stream)
    columns = reader.fieldnames
    curated = list(reader)
manuals = json.loads((ROOT / 'config/participant_manual_scores.json').read_text())
keys = {(r['chrom'].removeprefix('chr'), r['pos'], r['ref'], r['alt']) for r in curated}
if len(keys) != 31 or len({r['gene'] for r in curated}) != 28:
    raise ValueError('Expected 31 distinct baseline alleles in 28 genes')
if any(set(scores) != {r['gene'] for r in curated} for scores in manuals.values()):
    raise ValueError('Gene symbols disagree between curated variants and assessor scores')
positions = {(k[0], k[1]) for k in keys}
found = {}
record_count = 0
opener = gzip.open if args.full_vep.endswith('.gz') else open
with opener(args.full_vep, 'rt') as stream:
    for line in stream:
        if line.startswith('##INFO=<ID=CSQ,'):
            names = line.split('Format: ', 1)[1].split('"', 1)[0].split('|')
        if line.startswith('#'):
            continue
        record_count += 1
        c = line.rstrip('\r\n').split('\t')
        chrom = c[0].removeprefix('chr')
        if (chrom, c[1]) not in positions:
            continue
        info = dict(x.split('=', 1) if '=' in x else (x, '') for x in c[7].split(';'))
        annotations = [dict(zip(names, a.split('|'))) for a in info.get('CSQ', '').split(',')]
        sample = dict(zip(c[8].split(':'), c[9].split(':')))
        alts = c[4].split(',')
        trim = len({a[0] for a in [c[3]] + alts}) == 1
        called = re.split(r'[/|]', sample.get('GT', '.'))
        for index, alt in enumerate(alts, 1):
            key = (chrom, c[1], c[3], alt)
            if key not in keys:
                continue
            if str(index) not in called:
                raise ValueError('Curated ALT is not called in the sample: ' + str(key))
            vep_alt = (alt[1:] or '-') if trim else alt
            genes = {a.get('SYMBOL', '') for a in annotations if a.get('Allele') == vep_alt}
            if key in found:
                raise ValueError('Duplicate source allele: ' + str(key))
            found[key] = (sample, genes)
for r in curated:
    key = (r['chrom'].removeprefix('chr'), r['pos'], r['ref'], r['alt'])
    if key not in found:
        raise ValueError('Curated allele missing from VEP: ' + str(key))
    sample, genes = found[key]
    if r['gene'] not in genes:
        raise ValueError('Gene is not annotated to the selected ALT: ' + str(key))
    for field, value in [('DP', r['dp']), ('GQ', r['gq'])]:
        if int(sample[field]) != int(value):
            raise ValueError('Sample quality mismatch: ' + str(key))
    if sample['GT'].replace('|', '/') != r['gt'].replace('|', '/'):
        raise ValueError('Sample genotype mismatch: ' + str(key))
    if sample.get('AD') != r['ad_ref'] + ',' + r['ad_alt']:
        raise ValueError('Allele-depth mismatch: ' + str(key))
if record_count != 4943948:
    raise ValueError('Input record count differs from the documented full VEP')
args.out_dir.mkdir(parents=True, exist_ok=True)
with (args.out_dir / 'curated_variants.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n')
    writer.writeheader()
    writer.writerows(curated)
summary = {'source_records_scanned': record_count, 'variants': 31, 'genes': 28,
           'exact_called_alleles_verified': True, 'allele_specific_gene_annotations_verified': True,
           'DP_GQ_GT_AD_verified': True, 'selection': 'recorded manual curation; post-ranking additions excluded'}
(args.out_dir / 'curation_validation.json').write_text(json.dumps(summary, indent=2) + '\n')
print('Validated full VEP: 31 called alleles in 28 curated genes; DP, GQ, GT, AD and gene annotations agree.')
