#!/usr/bin/env python3
"""Replay exact CADD annotations from the audited Drive checkpoint, without network drift."""
import argparse,csv,gzip
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('input');p.add_argument('output');a=p.parse_args()
root=Path(__file__).resolve().parents[1]
with open(root/'data/annotations/cadd_v1.7_snapshot.tsv') as f:
 cache={(r['chrom'],r['pos'],r['ref'],r['alt']):r['cadd'] for r in csv.DictReader(f,delimiter='\t')}
seen=set(); count=0
with (gzip.open(a.input,'rt') if a.input.endswith('.gz') else open(a.input)) as f,open(a.output,'w') as out:
 for line in f:
  if line.startswith('#'):
   if line.startswith('#CHROM'):out.write('##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="Frozen CADD v1.7 exact-allele snapshot">\n')
   if not line.startswith('##INFO=<ID=CADD_PHRED,'):out.write(line)
   continue
  c=line.rstrip('\n').split('\t');key=(c[0],c[1],c[3],c[4])
  if key not in cache:raise ValueError('Allele absent from CADD checkpoint: '+str(key))
  seen.add(key);count+=1
  c[7]=';'.join(x for x in c[7].split(';') if not x.startswith('CADD_PHRED='))
  if cache[key]:c[7]+=';CADD_PHRED='+cache[key]
  out.write('\t'.join(c)+'\n')
if len(seen)!=len(cache):raise ValueError('Input does not cover complete annotation checkpoint')
print('Annotated records:',count)
