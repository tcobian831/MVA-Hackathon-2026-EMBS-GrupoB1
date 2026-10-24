#!/usr/bin/env python3
"""Final ranking stage: run the shared engine for one assessor or the four recorded assessors."""
import argparse,ast,contextlib,json,os,io
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--all-participants',action='store_true');p.add_argument('--manual',type=Path);p.add_argument('--name',default='custom');a=p.parse_args()
if a.all_participants and a.manual:p.error('Choose either --all-participants or --manual')
if not a.all_participants and not a.manual:p.error('Supply --manual or --all-participants')
manuals=json.loads((ROOT/'config/participant_manual_scores.json').read_text()) if a.all_participants else {a.name:json.loads(a.manual.read_text())}
all_results=[]
for name,scores in manuals.items():
 if not name.replace('_','').isalnum():raise ValueError('Use an alphanumeric assessor name')
 out=ROOT/'results'/name
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
 merged.to_csv(ROOT/'results/consensus.csv',index_label='gene')
 print(merged.to_string())
