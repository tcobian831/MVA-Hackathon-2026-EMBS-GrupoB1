#!/usr/bin/env python3
"""Final ranking stage: run the shared engine for one assessor or the four recorded assessors."""
import argparse,ast,contextlib,json,os,io
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--all-participants',action='store_true');p.add_argument('--manual',type=Path);p.add_argument('--name',default='custom');p.add_argument('--output-dir',type=Path);p.add_argument('--variants',type=Path,default=ROOT/'results/curated_variants.csv');a=p.parse_args()
if a.all_participants and a.manual:p.error('Choose either --all-participants or --manual')
if not a.all_participants and not a.manual:p.error('Supply --manual or --all-participants')
manuals=json.loads((ROOT/'config/participant_manual_scores.json').read_text()) if a.all_participants else {a.name:json.loads(a.manual.read_text())}
run_dir=a.output_dir or ROOT/'results'/'runs'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
run_dir.mkdir(parents=True,exist_ok=False)
os.environ['MVA_CURATED_CSV']=str(a.variants.resolve())
all_results=[]
for name,scores in manuals.items():
 if not name.replace('_','').isalnum():raise ValueError('Use an alphanumeric assessor name')
 out=run_dir/name
 if out.exists():raise ValueError('Output already exists: '+str(out))
 tree=ast.parse((ROOT/'scripts/matrix_engine.py').read_text())
 node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='MANUAL_SCORES' for t in n.targets))
 expected=set(ast.literal_eval(node.value))
 if set(scores)!=expected:raise ValueError('Manual scores must contain exactly the 28 current genes')
 for values in scores.values():
  if not isinstance(values,list) or len(values)!=3 or any(isinstance(v,bool) or not isinstance(v,(float,int)) or not np.isfinite(v) or not 0<=v<=2 for v in values):raise ValueError('Supply three finite scores between 0 and 2 for every gene')
 node.value=ast.parse(repr(scores),mode='eval').body;ast.fix_missing_locations(tree)
 os.environ['MVA_OUTPUT_DIR']=str(out)
 log=io.StringIO();ns={'__name__':'__main__','__file__':str(ROOT/'scripts/matrix_engine.py')}
 with contextlib.redirect_stdout(log):exec(compile(tree,str(ROOT/'scripts/matrix_engine.py'),'exec'),ns)
 (out/'execution.log').write_text(log.getvalue())
 frame=ns['summary'].set_index('gene');all_results.append(frame.definitive_score.rename(name))
 w=ns['W'];u=ns['units']; calculated=sum(u['u_'+k]*v for k,v in w.items())/sum(w.values())
 np.testing.assert_allclose(frame.loc[u.gene,'weighted_mean'],calculated,atol=1e-14,rtol=0)
 np.testing.assert_allclose(frame.definitive_score,.8*frame.weighted_mean+.2*frame.median_score,atol=1e-14,rtol=0)
 print(name,'validated: 28 genes, 10000 simulations')
if a.all_participants:
 merged=pd.concat(all_results,axis=1)
 if merged.isna().any().any():raise ValueError('Assessors have different gene sets')
 merged=merged*50;merged['mean_definitive_percent']=merged.mean(axis=1)
 merged=merged.sort_values('mean_definitive_percent',ascending=False);merged.insert(0,'rank',range(1,len(merged)+1))
 merged.to_csv(run_dir/'consensus.csv',index_label='gene')
 top=merged.head(10).copy()
 top.to_csv(run_dir/'top10_genes.csv',index_label='gene')
 baseline=pd.read_csv(a.variants,dtype={'pos':str,'ref':str,'alt':str})
 additions=pd.read_csv(ROOT/'config/post_ranking_additions.csv',dtype={'pos':str,'ref':str,'alt':str})
 baseline['in_ranking_baseline']=True
 additions['in_ranking_baseline']=False
 final=pd.concat([baseline,additions],ignore_index=True)
 final=final.loc[final.gene.isin(top.index)].copy()
 final['rank']=final.gene.map(top['rank'])
 final['mean_definitive_percent']=final.gene.map(top['mean_definitive_percent'])
 final['phase_status']=final.gene.map(lambda g: 'same variant on both copies (GT 1/1)' if g=='CTU2' else 'unknown between selected variants' if g in ('BUB1B','TACC2','LRP4') else 'single heterozygous candidate')
 final=final.sort_values(['rank','id'])
 cols=['rank','gene','mean_definitive_percent','id','variant','chrom','pos','ref','alt','gt','dp','gq','cadd','grpmax_pct','phase_status','in_ranking_baseline']
 final[cols].to_csv(run_dir/'top10_variants.csv',index=False)
 if len(final)!=13 or final[['chrom','pos','ref','alt']].duplicated().any():raise ValueError('Unexpected final candidate allele inventory')
 print(merged.to_string())
 print('Output directory:',run_dir.resolve())
