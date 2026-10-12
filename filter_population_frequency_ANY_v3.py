#!/usr/bin/env python3
"""Keep records if ANY CSQ/gnomADg_AF is below cutoff; keep missing AF.

Usage: python filter_population_frequency_ANY_v3.py final_candidates.vcf
Requires Python 3.8+; no third-party packages.
Only the global gnomADg_AF field is used, not population-specific fields.
"""
import argparse
import gzip
import math
import re
from pathlib import Path

VERSION = "ANY_v3_2026-10-12"


def frequencies(info, index):
    values = []
    csq = next((x[4:] for x in info.split(';') if x.startswith('CSQ=')), '')
    for annotation in csq.split(','):
        fields = annotation.split('|')
        if index >= len(fields):
            if annotation:
                raise ValueError('CSQ annotation shorter than header schema')
            continue
        for token in fields[index].split('&'):
            if token in ('', '.'):
                continue
            value = float(token)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError('Invalid gnomADg_AF: ' + token)
            values.append(value)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--threshold', type=float, default=0.01)
    parser.add_argument('-o', '--output', type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.threshold) or not 0 < args.threshold <= 1:
        parser.error('threshold must be in (0, 1]')
    stem = args.input.name.removesuffix('.gz').removesuffix('.vcf') if hasattr(str, 'removesuffix') else re.sub(r'\.vcf(?:\.gz)?$', '', args.input.name)
    tag = format(args.threshold, '.12g').replace('.', 'p')
    output = args.output or args.input.with_name(stem + '_AF_ANY_v3_lt_' + tag + '_keep_missing.vcf')
    if output.resolve() == args.input.resolve():
        parser.error('output must differ from input')
    if output.exists():
        parser.error('output already exists; choose another name with -o')
    total = passed = missing = different = 0
    headers = []
    kept = []
    index = None
    opener = gzip.open if args.input.suffix == '.gz' else open
    with opener(args.input, 'rt', encoding='utf-8', newline='') as source:
        for line_number, line in enumerate(source, 1):
            if line.startswith('#'):
                if line.startswith('##INFO=<ID=CSQ,'):
                    match = re.search(r'Format:\s*([^"\r\n]+)', line)
                    if not match:
                        raise ValueError('CSQ header has no Format schema')
                    names = match.group(1).strip().split('|')
                    if 'gnomADg_AF' not in names:
                        raise ValueError('CSQ header does not contain gnomADg_AF')
                    index = names.index('gnomADg_AF')
                headers.append(line)
                continue
            if not line.strip():
                continue
            if index is None:
                raise ValueError('Missing CSQ/gnomADg_AF header')
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) < 8:
                raise ValueError('Invalid VCF at line ' + str(line_number))
            if ',' in columns[4]:
                raise ValueError('Multiallelic record: split alleles before applying this filter')
            af = frequencies(columns[7], index)
            total += 1
            different += int(len(set(af)) > 1)
            # ANY is the actual selection rule. Missing values are retained.
            if not af:
                missing += 1
                kept.append(line)
            elif any(value < args.threshold for value in af):
                passed += 1
                kept.append(line)
    if not any(line.startswith('#CHROM\t') for line in headers):
        raise ValueError('Missing #CHROM header')
    with output.open('x', encoding='utf-8', newline='') as destination:
        for line in headers:
            if line.startswith('#CHROM\t'):
                destination.write('##population_frequency_filter="Version=' + VERSION + ';Field=CSQ/gnomADg_AF;Rule=ANY_AF_LT_THRESHOLD;Threshold=' + str(args.threshold) + ';KeepMissing=True"\n')
            destination.write(line)
        destination.writelines(kept)
    print('Version: ' + VERSION)
    print('Frequency rule: minimum used (ANY value below cutoff passes).')
    print('Field: CSQ/gnomADg_AF (global only)')
    print('Threshold (strict <): ' + str(args.threshold))
    print('Input records: ' + str(total))
    print('Passing frequency filter: ' + str(passed))
    print('Without frequency, retained: ' + str(missing))
    print('Records with differing frequency values: ' + str(different))
    print('Total output records: ' + str(len(kept)))
    print('Excluded records: ' + str(total - len(kept)))
    print('Output: ' + str(output.resolve()))


if __name__ == '__main__':
    main()
