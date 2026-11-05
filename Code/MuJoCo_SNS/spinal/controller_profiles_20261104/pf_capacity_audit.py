"""Read-only PF representability diagnostic using mirrored nonnegative decoders.

This is not an SNS simulation or a conductance lower bound: real MN dynamics,
reversal potentials, saturation, and reflexes are omitted deliberately. Compare
the existing PF projection mask against all same-side PF features on held-out
cycles before proposing extra cells or connections.
"""
import argparse,json,os,sys
from pathlib import Path
for key in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import numpy as np
from scipy.optimize import nnls
from controller_profiles_20261104 import profile as P
from muscle_map import classify
def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 meta=json.loads((a.fit/'MANIFEST.json').read_text());cfg=meta['configuration'];z=np.load(a.fit/'BEST_TRACE.npz');acts=z['acts'].tolist();pf=z['pf_names'].tolist();t=z['t'];T=cfg['target_period_s'];train=(t>=2*T)&(t<7*T);test=t>=7*T;train&=np.arange(len(t))%4==0;test&=np.arange(len(t))%4==0
 n=P.build(cfg['family'],acts,options=cfg.get('hunt_options'),config=cfg['controller_profile']);allowed={act:set()for act in acts}
 for edge,source,target in P.connections(n):
  if target.startswith('MN_')and target[3:]in allowed and source in pf and edge['params']['max_conductance']>0:allowed[target[3:]].add(source)
 selected=[i for i,act in enumerate(acts)if classify(act).groups[0]in('hip_ext','hip_flex','knee_ext','knee_flex','ankle_pf','ankle_df')];bases=sorted({classify(acts[i]).base for i in selected});modes={};weights={}
 for mode in('existing_projection_mask','all_ipsilateral_PF'):
  sq=[];total=0;rows=[];weights[mode]={}
  for base in bases:
   ids=[i for i in selected if classify(acts[i]).base==base];features=None
   for i in ids:
    side=acts[i][-1];pool=allowed[acts[i]]if mode=='existing_projection_mask'else{source for source in pf if source.endswith('_'+side)};keys={source[:-2]for source in pool};features=keys if features is None else features&keys
   features=sorted(features or[])
   if not features:
    residual=np.concatenate([-z['desired'][test,i]for i in ids]);sq.append(float(residual@residual));total+=len(residual);rows.append(dict(muscle=base,features=[],unrepresented=True,test_RMSE=float(np.sqrt(np.mean(residual**2)))));continue
   def design(i):return np.clip(z['pf'][:,[pf.index(key+'_'+acts[i][-1])for key in features]]/5.,0.,1.)
   X=np.concatenate([design(i)[train]for i in ids]);Y=np.concatenate([z['desired'][train,i]for i in ids]);w,_=nnls(X,Y);residual=np.concatenate([design(i)[test]@w-z['desired'][test,i]for i in ids]);sq.append(float(residual@residual));total+=len(residual);rows.append(dict(muscle=base,features=features,test_RMSE=float(np.sqrt(np.mean(residual**2))),actuators=[acts[i]for i in ids]));weights[mode][base]=dict(features=features,weights=w.tolist())
  modes[mode]=dict(test_RMSE=float(np.sqrt(sum(sq)/max(total,1))),represented_muscles=sum(not r.get('unrepresented',False)for r in rows),total_muscles=len(rows),rows=rows)
 result=dict(family=cfg['family'],source_fit=str(a.fit),training_cycles=[2,7],test_cycles=[7,12],feature_encoding='PF voltage clipped to0–5mV then normalized; diagnostic convention, not exact synaptic activation',modes=modes,weights=weights,scope=__doc__)
 (a.output/'RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(family=cfg['family'],modes={key:{k:v for k,v in value.items()if k!='rows'}for key,value in modes.items()}),indent=2))
if __name__=='__main__':main()
