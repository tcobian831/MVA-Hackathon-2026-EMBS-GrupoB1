#!/usr/bin/env python3
"""MVA final shortlist. Python 3.9+, standard library only; see README.md."""
import argparse
import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

VERSION = 'MVA_FINAL_v1_2026-10-14'
GROUPS = ('AFR', 'AMR', 'EAS', 'MID', 'NFE', 'SAS')
CONTROL = ('chr15', '40209701', 'T', 'G')  # GRCh38: BUB1B p.Leu737Ter
HIGH = set('stop_gained frameshift_variant splice_donor_variant splice_acceptor_variant'.split())
INTERMEDIATE = set('missense_variant inframe_insertion inframe_deletion start_lost stop_lost protein_altering_variant'.split())


def number(value):
    if value in ('', '.', None):
        return None
    n = float(value)
    if not math.isfinite(n):
        raise ValueError('Non-finite numeric annotation: ' + str(value))
    return n


def chr_name(s):
    return s if s.startswith('chr') else 'chr' + s


def own_allele(cols, info):
    if ',' in cols[4]:
        raise ValueError('Unsplit multiallelic record; normalize/split before this script')
    if info.get('ORIGINAL'):
        _, _, ref, alt, index = info['ORIGINAL'].split('|')
        alleles, index = [ref] + alt.split(','), int(index)
    else:
        alleles, index = [cols[3], cols[4]], 1
    if not 0 < index < len(alleles) or any(not a for a in alleles):
        raise ValueError('Invalid ORIGINAL allele index/alleles')
    # Match the allele representation used by Daniel's VEP annotation pipeline.
    if len({a[0] for a in alleles}) == 1:
        alleles = [a[1:] or '-' for a in alleles]
    return alleles[index]


def read_vcf(path, source):
    opener = gzip.open if str(path).endswith('.gz') else open
    header, names, records = [], None, []
    with opener(path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                header.append(line)
                if line.startswith('##INFO=<ID=CSQ,'):
                    names = line.split('Format: ', 1)[1].split('"', 1)[0].split('|')
                continue
            if not line.strip():
                continue
            cols = line.rstrip('\r\n').split('\t')
            if names is None or len(cols) < 10:
                raise ValueError('Expected CSQ header and a sample column in ' + str(path))
            info = dict(x.split('=', 1) if '=' in x else (x, '') for x in cols[7].split(';'))
            allele = own_allele(cols, info)
            required = ['Allele', 'Consequence', 'BIOTYPE', 'SYMBOL'] + ['gnomADg_' + g + '_AF' for g in GROUPS]
            if any(k not in names for k in required):
                raise ValueError('Missing required CSQ fields in ' + str(path))
            anns = []
            for raw in info.get('CSQ', '').split(','):
                fields = raw.split('|')
                if len(fields) != len(names):
                    raise ValueError('CSQ length disagrees with its header at ' + cols[0] + ':' + cols[1])
                a = dict(zip(names, fields))
                if a['Allele'] == allele:
                    anns.append(a)
            if not anns:
                raise ValueError('No allele-matched CSQ at ' + cols[0] + ':' + cols[1])
            vectors = {tuple(a['gnomADg_' + g + '_AF'] for g in GROUPS) for a in anns}
            if len(vectors) != 1:
                raise ValueError('Conflicting frequencies for the SAME allele; inspect input')
            raw_freq = next(iter(vectors))
            if any('&' in s for s in raw_freq):
                raise ValueError('Multiple values within an AF field; inspect allele mapping')
            freq = {g: number(s) for g, s in zip(GROUPS, raw_freq)}
            if any(n is not None and not 0 <= n <= 1 for n in freq.values()):
                raise ValueError('AF must be between 0 and 1')
            available = {g: n for g, n in freq.items() if n is not None}
            gm = max(available.values()) if available else None
            if 'F3_GRPMAX' in info and (gm is None or not math.isclose(gm, float(info['F3_GRPMAX']), rel_tol=1e-7)):
                raise ValueError('Recalculated groupmax disagrees with F3_GRPMAX')
            key = (chr_name(cols[0]), cols[1], cols[3], cols[4])
            control = key == CONTROL
            coding = [a for a in anns if a['BIOTYPE'] == 'protein_coding']
            h = [a for a in coding if set(a['Consequence'].split('&')) & HIGH]
            m = [a for a in coding if set(a['Consequence'].split('&')) & INTERMEDIATE]
            priority = 'high' if h or control else 'intermediate' if m else 'low'
            primary = h or m or [a for a in anns if a.get('CANONICAL') == 'YES'] or anns
            gtdata = dict(zip(cols[8].split(':'), cols[9].split(':')))
            gt = gtdata.get('GT', '.')
            alleles_gt = re.split(r'[/|]', gt)
            gtclass = ('Y_ploidy_review' if key[0] == 'chrY' else
                       'heterozygous' if len(alleles_gt) == 2 and '.' not in alleles_gt and len(set(alleles_gt)) == 2 else
                       'homozygous_alternate' if len(alleles_gt) == 2 and alleles_gt[0] == alleles_gt[1] and alleles_gt[0] not in ('0', '.') else
                       'other_or_missing')
            records.append(dict(key=key, cols=cols, info=info, source=source, cadd=number(info.get('CADD_PHRED')),
                                groupmax=gm, groups=','.join(g for g in GROUPS if freq[g] is not None),
                                window='MVA_WINDOW' in info.get('F1_REASON', '').split(','), control=control,
                                priority=priority, gt=gt, gtclass=gtclass,
                                phase=';'.join(k + '=' + gtdata[k] for k in ('PS', 'PID', 'PGT') if k in gtdata),
                                genes=sorted({a['SYMBOL'] for a in primary if a['SYMBOL']}),
                                allgenes=sorted({a['SYMBOL'] for a in anns if a['SYMBOL']}),
                                consequences=sorted({a['Consequence'] for a in primary}),
                                split_site=info.get('ORIGINAL', ''), status=''))
    if not any(l.startswith('#CHROM\t') for l in header):
        raise ValueError('Missing #CHROM header')
    return header, records


def write_vcf(path, header, records, rule):
    with path.open('x', encoding='utf-8', newline='') as f:
        for line in header:
            if line.startswith('#CHROM\t'):
                f.write('##INFO=<ID=MVA_FINAL_PRIORITY,Number=1,Type=String,Description="Functional review priority, not clinical classification">\n')
                f.write('##INFO=<ID=MVA_FINAL_STATUS,Number=1,Type=String,Description="selected or pending_no_cadd">\n')
                f.write('##MVA_FINAL=' + json.dumps(rule, separators=(',', ':')) + '\n')
            f.write(line)
        for r in records:
            c = list(r['cols'])
            c[7] += ';MVA_FINAL_PRIORITY=' + r['priority'] + ';MVA_FINAL_STATUS=' + r['status']
            f.write('\t'.join(c) + '\n')


def write_tsv(path, rows, fields):
    with path.open('x', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('main_vcf', type=Path, help="Daniel's filter3.vcf")
    p.add_argument('--pending-vcf', type=Path, help="Daniel's filter3_no_cadd.vcf")
    p.add_argument('--out-dir', type=Path, default=Path('final_shortlist_v1'))
    p.add_argument('--cadd-out', type=float, default=23)
    p.add_argument('--cadd-in', type=float, default=18)
    p.add_argument('--af', type=float, default=.0002, help='Proportion, not percent (default .0002 = .02%%)')
    args = p.parse_args()
    if not (all(math.isfinite(x) for x in (args.cadd_out, args.cadd_in, args.af)) and 0 < args.af <= 1):
        p.error('Invalid numeric thresholds')
    if args.out_dir.exists():
        p.error('Output directory exists; choose a new --out-dir (no overwriting)')
    sources = [('main', args.main_vcf)] + ([('pending', args.pending_vcf)] if args.pending_vcf else [])
    headers, records, seen = {}, [], set()
    for source, path in sources:
        headers[source], rr = read_vcf(path, source)
        for r in rr:
            if r['key'] in seen:
                raise ValueError('Duplicate CHROM/POS/REF/ALT across inputs: ' + str(r['key']))
            seen.add(r['key'])
            cutoff = args.cadd_in if r['window'] else args.cadd_out
            if r['control']:
                r['status'] = 'selected'
            elif r['groupmax'] is not None and r['groupmax'] >= args.af:
                r['status'] = 'rejected_frequency'
            elif r['cadd'] is None:
                r['status'] = 'pending_no_cadd'
            elif r['cadd'] > cutoff:
                r['status'] = 'selected'
            else:
                r['status'] = 'rejected_cadd'
            records.append(r)
    selected = [r for r in records if r['status'] == 'selected']
    pending = [r for r in records if r['status'] == 'pending_no_cadd']
    rule = dict(version=VERSION, cadd_out_strict_gt=args.cadd_out, cadd_in_strict_gt=args.cadd_in,
                groupmax_strict_lt=args.af, groups=GROUPS, protect_grch38=CONTROL)
    args.out_dir.mkdir(parents=True)
    # Source headers differ; keep provenance and annotations of each output intact.
    for source, _ in sources:
        for label, rr in [('selected', selected), ('pending', pending)]:
            subset = [r for r in rr if r['source'] == source]
            if subset or (source == 'main' and label == 'selected'):
                write_vcf(args.out_dir / (label + '_' + source + '.vcf'), headers[source], subset, rule)
    audit = []
    for r in records:
        audit.append(dict(chrom=r['key'][0], pos=r['key'][1], ref=r['key'][2], alt=r['key'][3],
                          status=r['status'], source=r['source'], cadd_phred=r['cadd'],
                          in_window=r['window'], groupmax_af=r['groupmax'], available_groups=r['groups'],
                          priority=r['priority'], genes=','.join(r['genes']), consequences=';'.join(r['consequences']),
                          GT=r['gt'], genotype_class=r['gtclass'], phase=r['phase'], original_site=r['split_site'],
                          genotype_note=('Split-site GT: 0 may represent another ALT; inspect original unsplit GT'
                                         if r['split_site'] else 'Reference at this site does not establish a reference gene copy')))
    fields = list(audit[0]) if audit else ['status']
    write_tsv(args.out_dir / 'audit_all.tsv', audit, fields)
    write_tsv(args.out_dir / 'selected.tsv', [r for r in audit if r['status'] == 'selected'], fields)
    write_tsv(args.out_dir / 'pending.tsv', [r for r in audit if r['status'] == 'pending_no_cadd'], fields)
    gene_records = defaultdict(list)
    for r in selected + pending:
        if r['gtclass'] == 'heterozygous':
            for g in r['allgenes']:
                gene_records[g].append(r)
    pairs = [dict(gene=g, variant=':'.join(r['key']), status=r['status'], priority=r['priority'],
                  GT=r['gt'], phase=r['phase'], interpretation='Multiple retained variants annotated to this gene; cis/trans not determined')
             for g, rr in sorted(gene_records.items()) if len(rr) > 1 for r in rr]
    write_tsv(args.out_dir / 'gene_pairs.tsv', pairs,
              ['gene', 'variant', 'status', 'priority', 'GT', 'phase', 'interpretation'])
    summary = dict(rule=rule, records_in=len(records), selected=len(selected), pending=len(pending),
                   selected_with_frequency=sum(r['groupmax'] is not None for r in selected),
                   selected_without_frequency=sum(r['groupmax'] is None for r in selected),
                   pending_with_frequency=sum(r['groupmax'] is not None for r in pending),
                   pending_without_frequency=sum(r['groupmax'] is None for r in pending),
                   status_counts=dict(Counter(r['status'] for r in records)),
                   selected_priorities=dict(Counter(r['priority'] for r in selected)),
                   pending_priorities=dict(Counter(r['priority'] for r in pending)),
                   selected_genotypes=dict(Counter(r['gtclass'] for r in selected)),
                   pending_genotypes=dict(Counter(r['gtclass'] for r in pending)),
                   control_found=any(r['control'] for r in records))
    (args.out_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(VERSION)
    print('Selected:', len(selected), '| with AF:', summary['selected_with_frequency'], '| without AF:', summary['selected_without_frequency'])
    print('Selected priorities:', summary['selected_priorities'])
    print('Pending without CADD:', len(pending), '(NOT part of selected candidates)')
    print('Known BUB1B control found:', summary['control_found'])
    print('Output:', args.out_dir.resolve())


if __name__ == '__main__':
    main()

