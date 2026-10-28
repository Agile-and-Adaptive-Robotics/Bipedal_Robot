from pathlib import Path
import os,sys,json,hashlib,importlib.util
for n in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[n]='1'
R=Path('D:/GitHub/Bipedal_Robot');S=R/'Code/MuJoCo_SNS/spinal';O=S/'reports_20261028_publication_flow/stimulus_long_controls';O.mkdir(exist_ok=True);sys.path.insert(0,str(S));os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from numba import njit
import deng_cpg_ode as D
spec=importlib.util.spec_from_file_location('coupled',S/'reports_20261026_dengmujoco_20261026/deng_mujoco_runner.py');M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)
import mujoco
(O/'PREREG.md').write_text('Long-window follow-up controls: apply each intervention for two baseline periods to test whole-burst suppression as well as within-window activity loss. Preserve shorter-pulse results. Preserve all initial negative results. Fresh 24 s baseline and full parity check retained. Add RG/PF extensor +50 nA for .5/1 s at original .25-cycle phase, and direct MN flexor inhibition -100/-500 nA for .5 s beginning .05 s after the first baseline activation-.5 rising crossing after 8 s. These are mechanistic controls, not paper-exact physiological-current replications. Same a-priori deletion/phase gates and suspended-body scope. No tuning of circuit/model parameters.\n')
@njit
def integrate(VR,src,dst,g,esyn,elo,ehi,na,idx,amp,t0,t1,dur=24.):
 dt=.1;n=int(dur*1000/dt);V=VR.copy();h=np.ones(len(V))*(1/(1+D.KH*np.exp(D.SH*(D.EH-D.VREST_HC))));h[12:]=0
 out=np.zeros((n//20+1,len(V)));out[0]=V
 for k in range(n):
  Is=np.zeros(len(V))
  for e in range(len(src)):
   sat=min(1.,max(0.,(V[src[e]]-elo[e])/(ehi[e]-elo[e])));Is[dst[e]]+=g[e]*sat*(esyn[e]-V[dst[e]])
  vv=V.copy();tt=k*dt/1000
  for j in range(len(V)):
   mi=1/(1+D.KM*np.exp(D.SM*(D.EM-V[j])));hi=1/(1+D.KH*np.exp(D.SH*(D.EH-V[j])));ina=D.G_NA*mi*h[j]**D.KH*(D.E_NA-V[j]) if na[j] else 0.
   inj=10. if j==0 and 1000<=k<1020 else 0.
   if j==idx and t0<=tt<t1:inj+=amp
   vv[j]=V[j]+(D.GM*(VR[j]-V[j])+inj+Is[j]+ina)*dt/D.CM
   h[j]+=dt*(hi-h[j])/D.TAU_H
  V=vv
  if (k+1)%20==0:out[(k+1)//20]=V
 return out
na=np.array([i in D.NA for i in range(len(D.NAMES))]);args=(D.VREST,D.src,D.dst,D.g,D.esyn,D.elo,D.ehi,na)
base=integrate(*args,-1,0.,0.,0.);t=np.arange(len(base))*.002
# Full reference trajectory from the existing NumPy integrator, same cadence/kick.
v=D.VREST.copy();h=np.ones(len(v))*(1/(1+D.KH*np.exp(D.SH*(D.EH-D.VREST_HC))));h[12:]=0;ref=np.empty_like(base);ref[0]=v
for k in range(len(base)-1):
 inj=np.zeros(len(v))
 if 50<=k<51:inj[0]=10.
 v,h=M.step_neural(v,h,inj);ref[k+1]=v
err=float(np.max(np.abs(ref-base)));(O/'backend_parity.json').write_text(json.dumps({'full_24s_max_voltage_difference_mV':err,'tolerance_mV':1e-7},indent=2));assert err<=1e-7,err
rg_on=M.onsets(t,base[:,0]);T=float(np.median(np.diff(rg_on[rg_on>4])));on=float(rg_on[rg_on>8][0]+.25*T)
conditions=[('baseline',-1,0.,0.),('RG_plus50_2cycles',0,50.,2*T),('PF_plus50_2cycles',4,50.,2*T),('MNflex_minus500_2cycles',13,-500.,2*T)]
act=lambda v:np.clip(1/(1+np.exp(.1532*(-70-v)))-.01,0,1)
def events(x,thr):return t[np.flatnonzero((x[1:]>thr)&(x[:-1]<=thr))+1]
base_motor=act(base[:,13]);bm=events(base_motor,.5);summ=[]
on_original=on
for name,idx,amp,duration in conditions:
 crossings=events(base_motor,.5);on=float(crossings[crossings>8][0]+.05) if name.startswith('MNflex') else on_original
 V=base if idx<0 else integrate(*args,idx,amp,on,on+duration);motor=act(V[:,13]);win=(t>=on)&(t<on+duration);active=win&(base_motor>.5);active_ms=float(active.sum()*2);lost_ms=float(np.sum(active&(motor<=.5))*2)
 deletion=active_ms>=20 and lost_ms>=.5*active_ms
 def phase(sig,baseline_events,thr):
  ev=events(sig,thr);ev=ev[ev>on+duration+T]
  if len(ev)<3:return None
  return float(np.median([np.min(np.abs(baseline_events-e))/T for e in ev[-4:]]))
 residual=phase(motor,bm,.5);rgres=phase(V[:,0],rg_on,-45.)
 status='no observable deletion'
 if deletion:status='no recovery' if residual is None else 'non-resetting deletion' if residual<=.05 else 'resetting deletion' if residual>=.15 else 'deletion; phase indeterminate'
 model=mujoco.MjModel.from_xml_path(str(M.XML));data=mujoco.MjData(model);model.dof_damping[:]*=.5
 aids={a:mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,a) for a in M.MUS};qadr=[model.jnt_qposadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,j)] for j in ['hip_L','knee_L','ankle_L']];q=np.zeros((len(t),3));q[0]=data.qpos[qadr]
 for k in range(len(t)-1):
  for a,ni in M.MUS.items():data.ctrl[aids[a]]=min(M.ACT(V[k+1,ni]),M.CAP)
  for z in range(2):mujoco.mj_step(model,data)
  q[k+1]=data.qpos[qadr]
 np.savez_compressed(O/(name+'.npz'),t=t,V=V,motor=act(V[:,12:14]),q=q)
 row={'condition':name,'target':D.NAMES[idx] if idx>=0 else None,'current_nA':amp,'start_s':on,'duration_s':duration,'baseline_period_s':T,'baseline_active_ms':active_ms,'lost_active_ms':lost_ms,'motor_phase_residual_cycles':residual,'RG_phase_residual_cycles':rgres,'classification':status,'finite':bool(np.isfinite(q).all() and np.isfinite(V).all())};summ.append(row);print(json.dumps(row),flush=True)
(O/'stimulus_results.json').write_text(json.dumps(summ,indent=2));(O/'source_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'deng_cpg_ode.py',Path(M.__file__),Path(M.XML)]},indent=2))
