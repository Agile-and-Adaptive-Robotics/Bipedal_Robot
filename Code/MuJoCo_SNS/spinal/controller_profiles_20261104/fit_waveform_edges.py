"""Current-balance-guided, mirrored fitting of existing PF→MN edge ratios.

No added PF cells/edges, target clamping, or neural/physical equation changes.
Use the discrete membrane-current balance to propose bounded conductances,
then integrate the full circuit and retain the best training-window result.
"""
import argparse,json,time,os,sys,hashlib
from pathlib import Path
for key in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parent))
import numpy as np
from scipy.optimize import lsq_linear
from pose_replay import ContractReplay
from readout_keys import edge_key,multiplier
from muscle_map import classify

def main():
 p=argparse.ArgumentParser();p.add_argument('--source-fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--iterations',type=int,default=10);a=p.parse_args();assert 1<=a.iterations<=20;a.output.mkdir(parents=True,exist_ok=False)
 old=json.loads((a.source_fit/'FINAL.json').read_text());meta=json.loads((a.source_fit/'MANIFEST.json').read_text());cfg=meta['configuration'];assert cfg['family']in('s3k','syn6','w2lvar')and cfg.get('feedforward_only');timing=Path(meta['source_timing_block']);fixed=json.loads((timing/'FINAL.json').read_text())['parameters'];r=ContractReplay(cfg);c=r.c;core=r.core;original=r.base_g.copy()
 from pf_design_gate import require_pf_design_ready
 require_pf_design_ready(cfg['family'])
 selected=np.flatnonzero(r.selected);groups={};by_act={}
 for i,(pre,post)in enumerate(zip(c.pre,c.post)):
  source,target=r.names[pre],r.names[post]
  if target.startswith('MN_')and source.startswith('PF_')and not source.startswith('PF_IN')and original[i]>0 and target[3:]in r.acts and r.selected[r.acts.index(target[3:])]:
   key=edge_key(source,target);groups.setdefault(key,[]).append(i);by_act.setdefault(r.acts.index(target[3:]),[]).append(i)
 assert not np.any(core.g_ion[:,r.mn[selected]])
 columns=np.flatnonzero(np.any(core.input_connectivity[r.mn[selected]]!=0.,axis=0));assert not np.any(r.u[:,columns]),'Account for nonzero direct MN inputs in current balance'
 initial={key:multiplier(r.names[c.pre[ids[0]]],r.names[c.post[ids[0]]],old)for key,ids in groups.items()}
 for key,ids in groups.items():assert np.max(abs(original[ids]-original[ids[0]]))<1e-12,'Homologous base conductances differ'
 train=(r.t>=2*r.T)&(r.t<7*r.T);test=r.t>=7*r.T;sample=train&(np.arange(len(r.t))%4==0);Vd=5*r.desired;previous=np.vstack((Vd[-1:],Vd[:-1]));dV=(Vd-previous)/r.dt;cap=r.dt/core.time_factor_membrane
 manifest=dict(source_timing_block=str(timing),configuration=cfg,source_peak_fit=str(a.source_fit),method='Discrete MN current balance, bounded per-existing-edge proposals, 0.35 damped update, actual full-network verification',gain_bounds=[.02,20.],training_cycles=[2,7],validation_cycles=[7,12],scope='Validation replays the same reference cycle at later times, not independent subjects; no body simulation',source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in(Path(__file__),HERE/'readout_keys.py',HERE/'replay_adapter.py',HERE/'pose_replay.py')})
 (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));before=np.empty((len(r.t),len(core.V)));forward=c.forward;counter=0
 def observed(u):
  nonlocal counter
  before[counter]=core.V;result=forward(u);counter+=1;return result
 c.forward=observed;gains=initial.copy();best=None;rows=[]
 for iteration in range(a.iterations):
  r.base_g[:]=original
  for key,ids in groups.items():r.base_g[ids]=original[ids]*gains[key]
  counter=0;result=r.evaluate(fixed,time.time()+180.);error=r.output[:,selected]-r.desired[:,selected];training=float(np.sqrt(np.mean(error[train]**2)));validation=float(np.sqrt(np.mean(error[test]**2)));row=dict(iteration=iteration,gains={},edge_gains=gains.copy(),training_RMSE=training,validation_RMSE=validation,**result);rows.append(row)
  if best is None or(result['diagnostics']['accepted']and not best['diagnostics']['accepted'])or(result['diagnostics']['accepted']==best['diagnostics']['accepted']and training<best['training_RMSE']):
   best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2));np.savez_compressed(a.output/'BEST_TRACE.npz',t=r.t,rg=r.rg_trace,pf=r.pf_trace,mn=r.output,desired=r.desired,acts=np.array(r.acts),rg_names=np.array([r.names[i]for i in r.rg]),pf_names=np.array([r.names[i]for i in r.pf]))
  (a.output/'ITERATIONS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(iteration=iteration,training_RMSE=training,validation_RMSE=validation,peak_RMSE=result['diagnostics']['mn_peak_rmse'])),flush=True)
  proposal={};bases=sorted({classify(r.acts[i]).base for i in by_act})
  for base in bases:
   acts=[i for i in by_act if classify(r.acts[i]).base==base];keys=sorted({edge_key(r.names[c.pre[e]],r.names[c.post[e]])for i in acts for e in by_act[i]});Xs=[];Ys=[]
   for i in acts:
    post=r.mn[i];incoming=np.flatnonzero(c.post==post);pfedges=by_act[i];others=np.array([e for e in incoming if e not in pfedges],dtype=int);known=np.full(len(r.t),core.i_b[post])
    if len(others):
     activation=np.clip((before[:,c.pre[others]]-c.elo[others])/c.erange[others],0.,1.);known+=np.sum(c.gmax[others]*activation*(c.reversal[others]-previous[:,i,None]),axis=1)
    rhs=cap[post]*dV[:,i]+core.g_m[post]*(previous[:,i]-core.V_rest[post])-known;X=np.zeros((len(r.t),len(keys)))
    for edge in pfedges:
     key=edge_key(r.names[c.pre[edge]],r.names[post]);X[:,keys.index(key)]=original[edge]*np.clip((before[:,c.pre[edge]]-c.elo[edge])/c.erange[edge],0.,1.)*(c.reversal[edge]-previous[:,i])
    Xs.append(X[sample]);Ys.append(rhs[sample])
   X=np.vstack(Xs);Y=np.concatenate(Ys);ridge=.002;X=np.vstack((X,np.sqrt(ridge)*np.eye(len(keys))));Y=np.r_[Y,np.sqrt(ridge)*np.array([gains[k]for k in keys])];fit=lsq_linear(X,Y,bounds=(.02,20.),tol=1e-8)
   proposal.update(zip(keys,map(float,fit.x)))
  gains={key:float(np.clip(.65*gains[key]+.35*proposal.get(key,gains[key]),.02,20.))for key in gains}
 # Reconstruct independently; the subsequent D worker supplies a true fresh-process check.
 fresh=ContractReplay(cfg)
 for i,(pre,post)in enumerate(zip(fresh.c.pre,fresh.c.post)):
  source,target=fresh.names[pre],fresh.names[post]
  if target.startswith('MN_')and source.startswith('PF_')and not source.startswith('PF_IN'):fresh.base_g[i]*=multiplier(source,target,best)
 repeat=fresh.evaluate(fixed,time.time()+180.);matched=abs(repeat['value']-best['value'])<1e-10;accepted=bool(best['diagnostics']['accepted']and repeat['diagnostics']['accepted']and matched)
 (a.output/'FINAL.json').write_text(json.dumps(dict(accepted=accepted,fresh_reconstruction_matched=matched,best=best,repeat=repeat,scope=manifest['scope']),indent=2))
if __name__=='__main__':main()
