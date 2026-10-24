#!/usr/bin/env python3
"""Revisit all annotated alleles for the 15 historical finalist genes; preserve original GT."""
import argparse,csv,gzip,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('vcf');p.add_argument('--output',default='results/expanded_review.csv');a=p.parse_args()
# This is the actual historical search panel, distinct from the final top ten.
genes=set('BUB1B MAD2L1BP SMC4 SVIL DIDO1 PRMT6 LZTR1 TACC2 PIK3C3 LRP4 CTU2 CITED2 TRDMT1 ROCK1 ACLY'.split())
cache={}
with open(ROOT/'data/annotations/cadd_v1.7_snapshot.tsv') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['cadd']:cache[(r['chrom'].removeprefix('chr'),r['pos'],r['ref'],r['alt'])]=float(r['cadd'])
for r in json.loads((ROOT/'data/annotations/expanded_cadd_snapshot.json').read_text()):cache[tuple(r['key'])]=r['cadd']
with open(ROOT/'config/curated_variants.csv') as f: baseline={(x['gene'],x['pos']) for x in csv.DictReader(f) if x['id'] not in ['34','35']}
fields=['gene','chrom','pos','ref','alt','GT','DP','GQ','AD','PGT','PID','PS','CADD','groupmax_fraction','consequence','passes_latest_search','review_decision']
out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
with open(out,'w') as o,(gzip.open(a.vcf,'rt') if a.vcf.endswith('.gz') else open(a.vcf)) as f:
 w=csv.DictWriter(o,fieldnames=fields);w.writeheader();n=0;count=0
 for l in f:
  if l.startswith('##INFO=<ID=CSQ,'):names=l.split('Format: ',1)[1].split('"',1)[0].split('|')
  if l.startswith('#'):continue
  count+=1
  # Cheap prefilter avoids parsing most records; exact symbols are checked below.
  if not any(g in l for g in genes):continue
  c=l.rstrip('\n').split('\t');info=dict(x.split('=',1) if '=' in x else (x,'') for x in c[7].split(';'))
  anns=[dict(zip(names,x.split('|'))) for x in info.get('CSQ','').split(',')]
  sample=dict(zip(c[8].split(':'),c[9].split(':')));gt=re.split(r'[/|]',sample.get('GT','.'));alts=c[4].split(',');alleles=[c[3]]+alts
  trim=len({x[0] for x in alleles})==1
  def numeric(s):return None if s in ('','. ','.',None) else float(s)
  dp=numeric(sample.get('DP'));gq=numeric(sample.get('GQ'))
  for index,alt in enumerate(alts,1):
   if str(index) not in gt:continue
   va=(alt[1:] or '-') if trim else alt
   aa=[x for x in anns if x.get('Allele')==va];syms={x.get('SYMBOL') for x in aa}&genes
   if not syms:continue
   freq={float(x['gnomADg_'+g+'_AF']) for x in aa for g in ['AFR','AMR','EAS','MID','NFE','SAS'] if x.get('gnomADg_'+g+'_AF') not in ('',None,'.')}
   gm=max(freq) if freq else None;k=(c[0].removeprefix('chr'),c[1],c[3],alt);score=cache.get(k)
   reasons=[]
   if dp is None or dp<5:reasons.append('DP<5 or missing')
   if gq is None or gq<15:reasons.append('GQ<15 or missing')
   if gm is not None and gm>.05:reasons.append('groupmax>5%')
   if score is not None and score<4:reasons.append('CADD<4')
   for gene in sorted(syms):
    accepted=(gene=='LRP4' and c[1]=='46878979') or (gene=='TACC2' and c[1]=='122194203')
    decision='accepted additional variant in post-ranking review; phase unknown' if accepted else ('not eligible: '+'; '.join(reasons) if reasons else 'not added: no additional sufficiently supported variant mechanism established in this review; passing exploratory cutoffs is not pathogenicity evidence')
    w.writerow(dict(gene=gene,chrom='chr'+k[0],pos=k[1],ref=k[2],alt=k[3],GT=sample.get('GT'),DP=dp,GQ=gq,AD=sample.get('AD',''),PGT=sample.get('PGT',''),PID=sample.get('PID',''),PS=sample.get('PS',''),CADD=score,groupmax_fraction=gm,consequence=';'.join(sorted({x.get('Consequence','') for x in aa if x.get('SYMBOL')==gene})),passes_latest_search=not reasons,review_decision=('already in baseline ranking; not a new addition' if (gene,c[1]) in baseline else decision)));n+=1
print('Scanned',count,'VCF records; wrote',n,'gene/allele review rows')
