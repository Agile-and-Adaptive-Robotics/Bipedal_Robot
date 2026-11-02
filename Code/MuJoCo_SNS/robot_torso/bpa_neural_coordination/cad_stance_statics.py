"""Actual-engine static-preload audit for one repaired stance, no integration."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import cad_stance_repair as S
P=S.P;C=S.C;R=S.R;mujoco=S.mujoco;OUT=S.OUT

def diagnose(s,force):
 m,d=s.m,s.d;J=-s.dynamic.length_jacobian[s.tids].T;d.qfrc_applied[:]=J@force;d.xfrc_applied[:]=0;d.xfrc_applied[m.body('pelvis').id,1]=.85*m.body_mass.sum()*9.81;s.dynamic.update(d);names=[]
 for i in np.argsort(-np.abs(d.qacc))[:6]:
  j=int(m.dof_jntid[i]);names.append(dict(dof=int(i),joint=m.joint(j).name,qacc=float(d.qacc[i]),mass_diagonal=None))
 domains=[]
 for j,name in enumerate(C.PAIR):
  p=s.bank.parameters[j];length=float(d.ten_length[s.tids[s.ids[j]]]+s.cfg['actuators'][name]['constant_remote_path_offset_m']);minimum=max(0.,p[4]*(length-p[2]-p[0]));domains.append(dict(name=name,length_m=length,zero_tension_restplusfixed_m=float(p[2]+p[0]),minimum_series_tension_for_domain_N=float(minimum),frozen_neural_reference_N=float(s.pairT[j]),max_zero_contraction_capacity_N=float(p[3])))
 return dict(dominant_accelerations=names,domain=domains,closure_gaps_m=s.constraints())

def solve():
 registration=dict(machine='easteregg2',scope='Offlineactual-mj_forward mechanicalpreload/contactdepthfit, not neuraltraining or newbehavior',pose='First registered feasibleexactstance, rightknee−.55rad and solvedleftheight',pair_force_frozen=True,other52_bounds=[.05,1.],penetration_bounds_m=[.0001,.004],harness_feedforward=.85,each_initial_height_anchor_explicit=True,physics_and_model_unchanged=True,max_nfev=20,forward_gate_max_qacc=.05,forward_gate_COM_acc_norm=.01,forward_gate_closure_m=1e-8,no_timeintegration=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());path=OUT/'ACTUAL_STATICS_PREREGISTRATION.json'
 if path.exists():assert json.loads(path.read_text())==registration
 else:path.write_text(json.dumps(registration,indent=2))
 s=S.Stance();report=json.loads((OUT/'OFFLINE_RESULTS.json').read_text());row=next(r for r in report['cases'] if r.get('geometric_height_root') and r['statics']['ideal_constraint_capacity_feasible']);angles=[row['right_knee_rad'],row['left_knee_rad']];s.pose(angles);initialforce=np.array(row['statics']['LP_force_N']);before=diagnose(s,initialforce);m,d=s.m,s.d;other=np.array([i for i in range(m.nu) if i not in s.ids]);qbase=d.qpos.copy();M=np.zeros((m.nv,m.nv));mujoco.mj_fullM(m,M,d.qM);scale=np.r_[np.ones(3),np.full(m.nv-3,.1)];calls=[0]
 def evaluate(x):
  d.qpos[:]=qbase;d.qpos[1]-=(x[-1]*.001-.0002);d.qvel[:]=0;s.dynamic.update(d);force=np.zeros(m.nu);force[other]=s.caps[other]*x[:-1];force[s.ids]=s.pairT;d.qfrc_applied[:]=-s.dynamic.length_jacobian[s.tids].T@force;d.xfrc_applied[:]=0;d.xfrc_applied[m.body('pelvis').id,1]=.85*m.body_mass.sum()*9.81;s.dynamic.update(d);calls[0]+=1;return force
 def residual(x):
  evaluate(x);return np.r_[(M@d.qacc)/scale,1e-6*x[:-1]]
 x0=np.r_[initialforce[other]/s.caps[other],.2];lo=np.r_[np.full(len(other),.05),.1];hi=np.r_[np.ones(len(other)),4.];x0=np.clip(x0,lo+1e-12,hi-1e-12);sol=least_squares(residual,x0,bounds=(lo,hi),max_nfev=20,ftol=1e-9,xtol=1e-9,gtol=1e-9);force=evaluate(sol.x);after=diagnose(s,force);jac=np.zeros((3,m.nv));mujoco.mj_jacSubtreeCom(m,d,jac,m.body('pelvis').id);acc=jac@d.qacc;passed=bool(np.max(np.abs(d.qacc))<=.05 and np.linalg.norm(acc)<=.01 and max(s.constraints())<=1e-8);result=dict(scope=registration['scope'],actual_statics_valid=passed,solver_success=bool(sol.success),nfev=int(sol.nfev),evaluations=calls[0],chosen_angles_rad=angles,penetration_m=float(sol.x[-1]*.001),qacc_max=float(np.max(np.abs(d.qacc))),COM_acceleration=acc.tolist(),ground_force_N=float(R.contact_loads(m,d).sum()),required_ground_force_N=float(.15*m.body_mass.sum()*9.81),before=before,after=after,pressure_live_allowed=False,no_closed_loop_run=True)
 (OUT/'ACTUAL_STATICS_RESULT.json').write_text(json.dumps(result,indent=2));np.savez_compressed(OUT/'ACTUAL_STATICS_STATE.npz',qpos=d.qpos.copy(),qvel=d.qvel.copy(),qacc=d.qacc.copy(),force_N=force,initial_force_N=initialforce,variables=sol.x);print(json.dumps(result,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--solve',action='store_true');a=p.parse_args()
 if a.solve:solve()
