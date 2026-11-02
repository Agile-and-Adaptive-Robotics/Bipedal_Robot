"""Separate full-CAD ideal54 static audit using frozenfixture-derived numerics."""
import argparse,json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import cad_ideal54_stance as F
import rigid_fourbar_numerics as V
P=F.B.P;R=F.B.R;mujoco=F.mujoco;OUT=F.B.OUT.parent/'rigid_numerics_stance'

class Native(F.All54):
 def __init__(self,dt=.002):
  super().__init__();self.numerical_manifest=V.apply(self.m,True,iterations=500,timestep=dt)

def accounting(b):
 m,d=b.m,b.d;M=np.zeros((m.nv,m.nv));mujoco.mj_fullM(m,M,d.qM);Qrig=np.zeros(m.nv);mujoco.mj_applyFT(m,d,d.xfrc_applied[b.pelvis,:3],d.xfrc_applied[b.pelvis,3:],d.xipos[b.pelvis],b.pelvis,Qrig);eom=M@d.qacc+d.qfrc_bias-d.qfrc_passive-d.qfrc_constraint-d.qfrc_applied-Qrig;rows=np.flatnonzero(d.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY));J=d.efc_J.reshape(d.nefc,m.nv)[rows];ranks={}
 for si,side in enumerate(('r','l')):
  dofs=[m.jnt_dofadr[m.joint(part+'_'+side+'_joint').id] for part in ('back','shank','front')];block=J[3*si:3*si+3,dofs];ranks[side]=dict(rank=int(np.linalg.matrix_rank(block,tol=1e-10)),singular_values=np.linalg.svd(block,compute_uv=False).tolist())
 force=b.force.copy();acc=d.qacc.copy();applied=d.qfrc_applied+Qrig;mujoco.mj_inverse(m,d);inverse_error=float(np.max(np.abs(d.qfrc_inverse-applied)));return dict(native_EOM_residual_max=float(np.max(np.abs(eom))),forward_inverse_force_residual_max=inverse_error,per_loop_hinge_Jacobians=ranks,physical_rank=sum(row['rank'] for row in ranks.values()),full_equality_J_singular_values=np.linalg.svd(J,compute_uv=False).tolist(),raw_full_rank=int(np.linalg.matrix_rank(J,tol=1e-10)),warning='Nearclosednonzerogaps can addtinyglobalrankdirections; no constraintrows were dropped or syntheticforcesadded',full_qacc=acc.tolist())

def solve():
 OUT.mkdir(exist_ok=True);selection=json.loads((V.DEFAULT_RECORD/'DERIVED_VARIANT.json').read_text());p=dict(machine='easteregg2',scope='NEW fullCAD mechanical-onlyideal54 static audit, not SNS or BPA',frozen_numerics=selection,iterations=500,dt=.002,stance='Same groundedzeroheadingfamily, sourcegeometry/mass/damping/capacities/floors unchanged',force_floors='52reference routes .05*MIF; formerlyfreepair0floor, all54preloads mayvary',maxfev=1500,method='Powell derivative-free after existingnonsmoothcontactaudit',loop_gap_bound_m=.0005,qacc_max_gate=.05,COMacc_norm_gate=.01,ground_error_gate_N=.1,each_foot_load_gate_N=1.,each_heel_toe_load_gate_N=.5,EOM_residual_gate=1e-6,forward5s_only_after_pass=True,fresh_then_dt1ms_after_hold_pass=True,no_stiffness_controller_rig_geometry_retune=True,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),Path(F.__file__),Path(V.__file__))});path=OUT/'PREREGISTRATION.json'
 if path.exists():assert json.loads(path.read_text())==p
 else:path.write_text(json.dumps(p,indent=2))
 b=Native();calls=[0]
 def objective(x):
  e,g,a,d,_,_=b.inverse_fit(x);calls[0]+=1;r=np.r_[e,g,a,d,1e-8*x];return float(r@r)
 start=np.array(json.loads((F.B.OUT/'STATIC_RESULT_default_fd.json').read_text())['variables']);lo=np.r_[np.full(6,-10.),.1,np.full(4,-4.)];hi=np.r_[np.full(6,10.),4.,np.full(4,4.)];t=time.perf_counter();sol=minimize(objective,start,method='Powell',bounds=list(zip(lo,hi)),options={'maxfev':1500,'xtol':1e-8,'ftol':1e-12});e,g,a,d,heights,fit=b.inverse_fit(sol.x);acc,loads=b.actual();actual_qacc=b.d.qacc.copy();gaps=b.s.constraints();knees=list(b.s.dynamic.coordinates(b.d).values());headings=[]
 for side in ('r','l'):
  axis=b.d.xmat[b.m.body('foot_'+side).id].reshape(3,3)[:,0];headings.append(float(np.rad2deg(np.arctan2(axis[2],axis[0]))))
 ledger=accounting(b);checks=dict(qacc=np.max(np.abs(actual_qacc))<=.05,COMacc=np.linalg.norm(acc)<=.01,ground=abs(loads.sum()-.15*b.weight)<=.1,bothfeet=bool(np.all(loads.sum(axis=1)>1)),heeltoe=bool(np.all(loads>.5)),loop=max(gaps)<=.0005,knees=max(knees)<=-.02,heading=max(np.abs(headings))<=10,heights=bool(np.all(heights<=0)&np.all(heights>=-.004)),EOM=ledger['native_EOM_residual_max']<1e-6,rank4=ledger['physical_rank']==4);checks={key:bool(value) for key,value in checks.items()};result=dict(accepted=all(checks.values()),checks=checks,solver_success=bool(sol.success),nfev=int(sol.nfev),wall_s=time.perf_counter()-t,qacc_max=float(np.max(np.abs(actual_qacc))),COM_acceleration=acc.tolist(),ground_N=float(loads.sum()),target_ground_N=.15*b.weight,contact_loads_N=loads.tolist(),closure_gaps_m=gaps,knees_rad=knees,heading_deg=headings,foot_min_y_m=heights.tolist(),pair_ideal_preload_N=b.force[b.s.ids].tolist(),inverse_residual_max_scaled=float(np.max(np.abs(e))),numerical_manifest=b.numerical_manifest,native_accounting=ledger,neural_or_pressure_simulation=False,source_geometry_mass_damping_unchanged=True)
 np.savez_compressed(OUT/'PRELOAD_STATE.npz',qpos=b.d.qpos.copy(),qvel=b.d.qvel.copy(),qacc=actual_qacc,force_N=b.force.copy(),variables=sol.x,root_anchor=b.d.qpos[:3].copy());(OUT/'STATIC_RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--solve',action='store_true');a=parser.parse_args()
 if a.solve:solve()
