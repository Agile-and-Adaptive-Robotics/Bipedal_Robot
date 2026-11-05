"""Bounded experimental PF-input calibration; frozen RG parameters and MN weights.

Fits mirrored existing upstream conductances by discrete LIF current balance.
Not installed into the body adapter or eligible for automatic curriculum advance.
"""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
for name in('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parent)]
import numpy as np
from scipy.optimize import lsq_linear
import motor_adapter as M
from pf_timing_audit import features,primitive_target

def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path);p.add_argument('--output',type=Path);p.add_argument('--iterations',type=int,default=8);p.add_argument('--verify',type=Path);a=p.parse_args()
 if a.verify:
  out=a.verify;manifest=json.loads((out/'MANIFEST.json').read_text());fit=Path(manifest['source_D_fit']);verify=True
 else:
  if not a.fit or not a.output or not 1<=a.iterations<=12:p.error('Bounded fit/output required')
  out=a.output;fit=a.fit;out.mkdir(parents=True,exist_ok=False);verify=False
 final=json.loads((fit/'FINAL.json').read_text());dm=json.loads((fit/'MANIFEST.json').read_text());cf=Path(dm['fit']);cm=json.loads((cf/'MANIFEST.json').read_text());cfg=cm['configuration']
 if cfg['family']!='syn6':raise ValueError('Explicit six-channel primitive targets required')
 ready=M.handle_request(dict(op='preflight',fit=str(cf),fit_sha256=hashlib.sha256((cf/'FINAL.json').read_bytes()).hexdigest()))
 if not ready['ready_for_optimization']:raise ValueError('Frozen readout dependency failed')
 r=M.INSTANCE;fixed={**M.FIXED,**final['parameters']};c=r.c;core=r.core;names=r.names;pf=r.pf
 assert not np.any(core.g_ion[:,pf]),'Account for PF ionic currents before inversion'
 fsa_path=HERE.parent/'fsa_results/fsa_backsolve.npz';fsa=np.load(fsa_path);target=np.column_stack([primitive_target(fsa,cfg,name[-1],r.t,int(name.split('_')[1][1:])-1)for name in [names[i]for i in pf]])
 desired=5*target;previous=np.vstack((desired[-1:],desired[:-1]));derivative=(desired-previous)/r.dt;rhs=(r.dt/core.time_factor_membrane[pf])*derivative+core.g_m[pf]*(previous-core.V_rest[pf])-core.i_b[pf]
 edges=np.flatnonzero(np.isin(c.post,pf));post_column=np.array([np.where(pf==c.post[i])[0][0]for i in edges]);keys=[names[c.pre[i]][:-1]+'|'+names[c.post[i]][:-1]for i in edges];groups={key:np.array([j for j,k in enumerate(keys)if k==key])for key in sorted(set(keys))};key_order=list(groups);key_id=np.array([key_order.index(k)for k in keys]);static=~np.isin(edges,r.net.command_runtime.slots);base_static=r.base_g[edges].copy();buffer=np.empty((len(r.t),len(edges)));counter=0;forward=c.forward;coeff=np.ones(len(groups));capture=False
 def advance(u):
  nonlocal counter
  base=c.gmax[edges].copy();base[static]=base_static[static]
  if capture:
   activation=np.clip((core.V[c.pre[edges]]-c.elo[edges])/c.erange[edges],0.,1.);buffer[counter]=base*activation*(c.reversal[edges]-previous[counter,post_column])
  c.gmax[edges]=base*coeff[key_id];counter+=1;return forward(u)
 c.forward=advance
 phase=np.arange(200)/200.;keep=r.t>=2*r.T;fit_rows=keep&(r.t<7*r.T)&(np.arange(len(r.t))%4==0)
 def evaluate(values,observe=False):
  nonlocal coeff,counter,capture
  coeff=values.copy();counter=0;capture=observe;baseline=r.evaluate(fixed,time.time()+180.);rows=[]
  for col,name in enumerate([names[i]for i in pf]):
   actual=np.array([np.interp((phase+n)*r.T,r.t,np.clip(r.pf_trace[:,col]/5.,0.,1.))for n in range(2,11)]).mean(axis=0);wanted=primitive_target(fsa,cfg,name[-1],phase*r.T,int(name.split('_')[1][1:])-1);af=features(actual);tf=features(wanted);lag=(tf['fundamental_phase_deg']-af['fundamental_phase_deg']+180)%360-180;duty=abs(af['duty_fraction']-tf['duty_fraction']);counts=len(af['onset_phase_percent'])==len(tf['onset_phase_percent'])and len(af['offset_phase_percent'])==len(tf['offset_phase_percent']);passed=bool(abs(lag)<=30 and duty<=.15 and counts and af['excursion']>=.2)
   rows.append(dict(cell=name,lag_deg=float(lag),duty_error=float(duty),burst_counts_match=counts,excursion=af['excursion'],maximum_activation_error=float(max(abs(actual-wanted))),mean_absolute_activation_error=float(np.mean(abs(actual-wanted))),passed=passed))
  value=float(np.mean([abs(row['lag_deg'])/180+2*row['duty_error']+float(not row['burst_counts_match'])+row['maximum_activation_error']+.5*row['mean_absolute_activation_error']for row in rows]));rg_ok=baseline['diagnostics']['relative_period_error']<=.1 and baseline['diagnostics']['bilateral_phase_error']<=.15
  return dict(value=value,all_PF_timing_passed=bool(all(row['passed']for row in rows)and rg_ok),PF=rows,RG=baseline['diagnostics'])
 if verify:
  best=json.loads((out/'BEST.json').read_text());result=evaluate(np.array([best['multipliers'][key]for key in key_order]));matched=abs(result['value']-best['metrics']['value'])<1e-10;(out/'FRESH_REPEAT.json').write_text(json.dumps(dict(exact_objective_match=matched,metrics=result),indent=2));return
 manifest=dict(fit=str(cf),source_D_fit=str(fit),configuration=cfg,source_timing_block=cm['source_timing_block'],scope='Experimental PF upstream timing fit only; original RG/MN parameters fixed; existing edges and mirrored multipliers; no body integration or installed controller revision',gain_bounds=[.001,20.],iterations=a.iterations,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in(Path(__file__),HERE/'pf_timing_audit.py',HERE/'motor_adapter.py',HERE/'replay_adapter.py',fsa_path)})
 (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));best=None;history=[];deadline=time.time()+600
 for iteration in range(a.iterations):
  if time.time()>deadline:break
  result=evaluate(coeff,True);row=dict(iteration=iteration,multipliers=dict(zip(key_order,map(float,coeff))),metrics=result);history.append(row)
  if best is None or result['value']<best['metrics']['value']:
   best=row;(out/'BEST.json').write_text(json.dumps(best,indent=2));np.savez_compressed(out/'BEST_TRACE.npz',t=r.t,rg=r.rg_trace,pf=r.pf_trace,mn=r.output,desired=r.desired,acts=np.array(r.acts),rg_names=np.array([names[i]for i in r.rg]),pf_names=np.array([names[i]for i in pf]))
  (out/'ITERATIONS.json').write_text(json.dumps(history,indent=2));print(json.dumps(dict(iteration=iteration,PF_loss=result['value'],PF_passers=sum(row['passed']for row in result['PF']))),flush=True)
  proposal=coeff.copy()
  for base in sorted({names[i][:-1]for i in pf}):
   targets=[j for j,i in enumerate(pf)if names[i][:-1]==base];columns=[j for j,key in enumerate(key_order)if key.endswith('|'+base)];X=[];Y=[]
   for post in targets:
    block=np.zeros((fit_rows.sum(),len(columns)))
    for k,column in enumerate(columns):
     slots=[j for j in groups[key_order[column]]if post_column[j]==post]
     if slots:block[:,k]=buffer[fit_rows][:,slots].sum(axis=1)
    X.append(block);Y.append(rhs[fit_rows,post])
   X=np.vstack(X);Y=np.concatenate(Y);ridge=.01;X=np.vstack((X,np.sqrt(ridge)*np.eye(len(columns))));Y=np.r_[Y,np.sqrt(ridge)*np.ones(len(columns))];solution=lsq_linear(X,Y,bounds=(.001,20.),tol=1e-7,max_iter=200);proposal[columns]=solution.x
  coeff=.7*coeff+.3*proposal
 subprocess.run([sys.executable,str(Path(__file__)),'--verify',str(out)],check=True,timeout=180,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));repeat=json.loads((out/'FRESH_REPEAT.json').read_text());(out/'FINAL.json').write_text(json.dumps(dict(PF_timing_passed=best['metrics']['all_PF_timing_passed'],fresh_process_matched=repeat['exact_objective_match'],accepted=False,experimental_only=True,body_adapter_integration=False,best=best,scope=manifest['scope']),indent=2))
if __name__=='__main__':main()
