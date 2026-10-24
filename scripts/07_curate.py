#!/usr/bin/env python3
"""Replay an explicit human-selected candidate set; never call this an automatic discovery filter."""
import argparse,csv,gzip,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('full_vep');p.add_argument('--out-dir',default='results');a=p.parse_args()
root=Path(__file__).resolve().parents[1];out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
with open(root/'config/curated_variants.csv') as f: rows=list(csv.DictReader(f))
expected={(r['chrom'].removeprefix('chr'),r['pos']) for r in rows};found={};names=[];count=0
with (gzip.open(a.full_vep,'rt') if a.full_vep.endswith('.gz') else open(a.full_vep)) as f:
 for line in f:
  if line.startswith('##INFO=<ID=CSQ,'):names=line.split('Format: ',1)[1].split('"',1)[0].split('|')
  if line.startswith('#'):continue
  count+=1;c=line.rstrip('\n').split('\t');k=(c[0].removeprefix('chr'),c[1])
  if k not in expected:continue
  info=dict(x.split('=',1) if '=' in x else (x,'') for x in c[7].split(';')); syms={x.split('|')[names.index('SYMBOL')] for x in info.get('CSQ','').split(',')}
  found[k]=(c,syms)
for r in rows:
 k=(r['chrom'].removeprefix('chr'),r['pos'])
 if k not in found:raise ValueError('Curated position missing: '+str(k))
 c,syms=found[k]
 if r['gene'] not in syms:raise ValueError('Gene association not present: '+str(k)+' '+r['gene'])
 sample=dict(zip(c[8].split(':'),c[9].split(':')))
 if int(sample['DP'])!=int(r['dp']) or int(sample['GQ'])!=int(r['gq']):raise ValueError('Quality mismatch: '+str(k))
assert len(rows)==31 and len({r['gene'] for r in rows})==28
with open(out/'curated_variants.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(out/'curation_validation.json').write_text(json.dumps({'source_records_scanned':count,'variants':31,'genes':28,'all_positions_and_gene_symbols_verified':True,'selection':'explicit manual decisions; excludes two post-ranking additions'},indent=2)+'\n')
print('Validated',count,'source records; 31 curated variants in 28 genes')
