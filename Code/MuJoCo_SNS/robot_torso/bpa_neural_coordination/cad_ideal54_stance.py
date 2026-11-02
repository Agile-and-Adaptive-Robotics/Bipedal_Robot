"""Mechanical-only all54 preload; preserves52floors and freepair0floor."""
import json,hashlib,argparse,time
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear,minimize
import cad_compliant_stance as B
OUT=B.OUT/'ideal54';mujoco=B.mujoco

class All54(B.Mechanical):
 def inverse_fit(self,x):
  heights=self.pose(x);m,d=self.m,self.d;d.qacc[:]=0;mujoco.mj_inverse(m,d);Qrig=np.zeros(m.nv);mujoco.mj_applyFT(m,d,[0,.85*self.weight,0],[0,0,0],d.xipos[self.pelvis],self.pelvis,Qrig);A=-self.s.dynamic.length_jacobian[self.s.tids].T;lower=np.full(m.nu,.05);lower[self.s.ids]=0.;target=d.qfrc_inverse-Qrig;matrix=A*self.s.caps;fit=lsq_linear(self.scale[:,None]*matrix,self.scale*target,bounds=(lower,np.ones(m.nu)),method='bvls',tol=1e-12,max_iter=200);self.force[:]=fit.x*self.s.caps;error=self.scale*(A@self.force+Qrig-d.qfrc_inverse);gaps=np.array(self.s.constraints());return error,np.maximum(gaps-.0005,0)/.00005,np.maximum(heights+.00001,0)/.0001,np.maximum(-.004-heights,0)/.0001,heights,fit

def solve():
 OUT.mkdir(exist_ok=True);p=dict(machine='easteregg2',scope='Mechanical-onlyall54preload, NO SNS/pressure case',correction='All54mechanicalforces freewithin originalcapacities; prior202Npair-fixed negatives preserved',force_floors='Other52=.05*MIF; previouslyunclampedpair0; no floorresize',source_physics_mass_constraints_unchanged=True,one_stance_family_only=True,loop_gap_bound_m=.0005,actual_qacc_gate=.05,actual_COMacc_gate=.01,actual_ground_error_gate_N=.1,method='Powell after quantified nonsmoothderivatives',maxfev=1000,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),Path(B.__file__))});path=OUT/'PREREGISTRATION.json'
 if path.exists():assert json.loads(path.read_text())==p
 else:path.write_text(json.dumps(p,indent=2))
 b=All54();calls=[0]
 def objective(x):
  e,g,a,d,_,_=b.inverse_fit(x);calls[0]+=1;r=np.r_[e,g,a,d,1e-8*x];return float(r@r)
 start=np.array(json.loads((B.OUT/'STATIC_RESULT_default_fd.json').read_text())['variables']);lo=np.r_[np.full(6,-10.),.1,np.full(4,-4.)];hi=np.r_[np.full(6,10.),4.,np.full(4,4.)];t=time.perf_counter();sol=minimize(objective,start,method='Powell',bounds=list(zip(lo,hi)),options={'maxfev':1000,'xtol':1e-7,'ftol':1e-10});e,g,a,d,heights,fit=b.inverse_fit(sol.x);acc,loads=b.actual();gaps=b.s.constraints();knees=list(b.s.dynamic.coordinates(b.d).values());headings=[]
 for side in ('r','l'):
  axis=b.d.xmat[b.m.body('foot_'+side).id].reshape(3,3)[:,0];headings.append(float(np.rad2deg(np.arctan2(axis[2],axis[0]))))
 checks=dict(qacc=np.max(np.abs(b.d.qacc))<=.05,COMacc=np.linalg.norm(acc)<=.01,ground=abs(loads.sum()-.15*b.weight)<=.1,bothfeet=bool(np.all(loads.sum(axis=1)>1)),loop=max(gaps)<=.0005,knees=max(knees)<=-.02,heading=max(np.abs(headings))<=10,heights=bool(np.all(heights<=0)&np.all(heights>=-.004)));checks={key:bool(value) for key,value in checks.items()};result=dict(accepted=all(checks.values()),checks=checks,solver_success=bool(sol.success),nfev=int(sol.nfev),wall_s=time.perf_counter()-t,qacc_max=float(np.max(np.abs(b.d.qacc))),COM_acceleration=acc.tolist(),ground_N=float(loads.sum()),target_ground_N=.15*b.weight,contact_loads_N=loads.tolist(),closure_gaps_m=gaps,knees_rad=knees,heading_deg=headings,foot_min_y_m=heights.tolist(),pair_ideal_preload_N=b.force[b.s.ids].tolist(),inverse_residual_max_scaled=float(np.max(np.abs(e))),neural_or_pressure_simulation=False,physics_changed=False)
 np.savez_compressed(OUT/'PRELOAD_STATE.npz',qpos=b.d.qpos.copy(),qvel=b.d.qvel.copy(),qacc=b.d.qacc.copy(),force_N=b.force.copy(),variables=sol.x,root_anchor=b.d.qpos[:3].copy());(OUT/'STATIC_RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--solve',action='store_true');args=p.parse_args()
 if args.solve:solve()
