"""One numericalcontinuation of the same frozen fullCADstaticcase."""
import json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,minimize
import cad_rigid_stance as C
OUT=C.OUT

def main():
 b=C.Native();start=C.P.load(OUT/'PRELOAD_STATE.npz')['variables'];calls=[0]
 def residual(x):
  e,g,a,d,_,_=b.inverse_fit(x);calls[0]+=1;return np.r_[e,g,a,d,1e-8*x]
 def jac(x,h):
  columns=[]
  for j in range(len(x)):
   step=np.zeros(len(x));step[j]=h;columns.append((residual(x+step)-residual(x-step))/(2*h))
  return np.array(columns).T
 J1=jac(start,.001);J2=jac(start,.01);difference=float(np.linalg.norm(J1-J2)/max(np.linalg.norm(J2),1e-12));method='centralJacleast_squares' if difference<.5 else 'Powellcontinuation';p=dict(scope='SamefullCADcase numericalcontinuation; no newphysics/controller/rig/forcebound',starting_state=str(OUT/'PRELOAD_STATE.npz'),prior_result_preserved=True,relative_J_difference=difference,method=method,max_nfev=50,max_Powell_fevals=1500,central_step_scaled=.001,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());file=OUT/'REFINEMENT_PREREGISTRATION.json'
 if file.exists():assert json.loads(file.read_text())==p
 else:file.write_text(json.dumps(p,indent=2))
 lo=np.r_[np.full(6,-10.),.1,np.full(4,-4.)];hi=np.r_[np.full(6,10.),4.,np.full(4,4.)];t=time.perf_counter()
 if difference<.5:sol=least_squares(residual,start,jac=lambda x:jac(x,.001),bounds=(lo,hi),max_nfev=50,ftol=1e-13,xtol=1e-13,gtol=1e-13)
 else:sol=minimize(lambda x:float(residual(x)@residual(x)),start,method='Powell',bounds=list(zip(lo,hi)),options={'maxfev':1500,'xtol':1e-9,'ftol':1e-12})
 e,g,a,d,heights,fit=b.inverse_fit(sol.x);acc,loads=b.actual();qacc=b.d.qacc.copy();gaps=b.s.constraints();knees=list(b.s.dynamic.coordinates(b.d).values());headings=[]
 for side in ('r','l'):
  axis=b.d.xmat[b.m.body('foot_'+side).id].reshape(3,3)[:,0];headings.append(float(np.rad2deg(np.arctan2(axis[2],axis[0]))))
 ledger=C.accounting(b);checks=dict(qacc=np.max(np.abs(qacc))<=.05,COMacc=np.linalg.norm(acc)<=.01,ground=abs(loads.sum()-.15*b.weight)<=.1,bothfeet=bool(np.all(loads.sum(axis=1)>1)),heeltoe=bool(np.all(loads>.5)),loop=max(gaps)<=.0005,knees=max(knees)<=-.02,heading=max(np.abs(headings))<=10,heights=bool(np.all(heights<=0)&np.all(heights>=-.004)),EOM=ledger['native_EOM_residual_max']<1e-6,rank4=ledger['physical_rank']==4);checks={k:bool(v) for k,v in checks.items()};result=dict(accepted=all(checks.values()),checks=checks,method=method,relative_J_difference=difference,solver_success=bool(sol.success),nfev=int(sol.nfev),total_inverse_evaluations=calls[0],wall_s=time.perf_counter()-t,qacc_max=float(np.max(np.abs(qacc))),COM_acceleration=acc.tolist(),ground_N=float(loads.sum()),target_ground_N=.15*b.weight,contact_loads_N=loads.tolist(),closure_gaps_m=gaps,knees_rad=knees,heading_deg=headings,foot_min_y_m=heights.tolist(),pair_ideal_preload_N=b.force[b.s.ids].tolist(),inverse_residual_max_scaled=float(np.max(np.abs(e))),native_accounting=ledger,neural_or_pressure_simulation=False,numerical_parameters_frozen=True)
 np.savez_compressed(OUT/'REFINED_PRELOAD_STATE.npz',qpos=b.d.qpos.copy(),qvel=b.d.qvel.copy(),qacc=qacc,force_N=b.force.copy(),variables=sol.x,root_anchor=b.d.qpos[:3].copy());(OUT/'REFINEMENT_RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
