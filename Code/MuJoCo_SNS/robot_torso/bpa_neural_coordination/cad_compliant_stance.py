"""Bounded native-compliance mechanical initialization; no pressure/SNS claim."""
import argparse,json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,lsq_linear,minimize
from scipy.spatial.transform import Rotation
import cad_stance_repair as S
P=S.P;R=S.R;mujoco=S.mujoco;OUT=S.OUT/'compliance05mm'

class Mechanical:
 def __init__(self):
  self.s=S.Stance();self.m=self.s.m;self.d=self.s.d;rows=json.loads((S.OUT/'OFFLINE_RESULTS.json').read_text())['cases'];self.source=next(row for row in rows if row.get('geometric_height_root') and row['statics']['ideal_constraint_capacity_feasible']);self.base=self.s.pose([self.source['right_knee_rad'],self.source['left_knee_rad']]);self.pelvis=self.m.body('pelvis').id;self.weight=self.m.body_mass.sum()*9.81;self.hinge=[int(self.m.jnt_qposadr[self.m.joint(p+'_'+s+'_joint').id]) for s in ('r','l') for p in ('back','shank','front')];self.other=np.array([i for i in range(self.m.nu) if i not in self.s.ids]);self.scale=np.r_[np.ones(3),np.full(self.m.nv-3,20.)];self.force=np.zeros(self.m.nu);self.x=np.r_[np.zeros(6),.2,np.zeros(4)]
 def pose(self,x):
  m,d=self.m,self.d;d.qpos[:]=self.base;d.qpos[1]=.9;d.qvel[:]=0;d.qfrc_applied[:]=0;d.xfrc_applied[:]=0;d.qpos[self.hinge]+=x[:6]*.001
  for side in ('r','l'):
   a=m.jnt_qposadr[m.joint('thigh_'+side+'_joint').id];d.qpos[a:a+4]=[1,0,0,0]
  mujoco.mj_kinematics(m,d);minima=[]
  for j,side in enumerate(('r','l')):
   goal=Rotation.from_rotvec([x[7+2*j]*.001,0,x[8+2*j]*.001]).as_matrix()@self.s.soles[side]['matrix'];forward=goal[:,0];goal=Rotation.from_rotvec([0,np.arctan2(forward[2],forward[0]),0]).as_matrix()@goal;bid=m.body('foot_'+side).id;cur=d.xmat[bid].reshape(3,3);quat=Rotation.from_matrix(goal@cur.T).as_quat();a=m.jnt_qposadr[m.joint('thigh_'+side+'_joint').id];d.qpos[a:a+4]=np.r_[quat[3],quat[:3]];minima.append(float((self.s.soles[side]['local_vertices']@goal.T)[:,1].min()))
  mujoco.mj_kinematics(m,d);ys=np.array([d.xpos[m.body('foot_'+side).id,1]+minima[j] for j,side in enumerate(('r','l'))]);d.qpos[1]-=ys.mean()+x[6]*.001;self.s.dynamic.update(d);return ys-ys.mean()-x[6]*.001
 def inverse_fit(self,x):
  heights=self.pose(x);m,d=self.m,self.d;d.qacc[:]=0;mujoco.mj_inverse(m,d);Qrig=np.zeros(m.nv);mujoco.mj_applyFT(m,d,[0,.85*self.weight,0],[0,0,0],d.xipos[self.pelvis],self.pelvis,Qrig);A=-self.s.dynamic.length_jacobian[self.s.tids].T;target=d.qfrc_inverse-Qrig-A[:,self.s.ids]@self.s.pairT;matrix=A[:,self.other]*self.s.caps[self.other];fit=lsq_linear(self.scale[:,None]*matrix,self.scale*target,bounds=(.05,1.),method='bvls',tol=1e-12,max_iter=200);self.force[:]=0;self.force[self.other]=fit.x*self.s.caps[self.other];self.force[self.s.ids]=self.s.pairT;error=self.scale*(A@self.force+Qrig-d.qfrc_inverse);gaps=np.array(self.s.constraints());barrier=np.maximum(gaps-.0005,0)/.00005;above=np.maximum(heights+.00001,0)/.0001;deep=np.maximum(-.004-heights,0)/.0001;return error,barrier,above,deep,heights,fit
 def actual(self):
  d=self.d;d.qfrc_applied[:]=-self.s.dynamic.length_jacobian[self.s.tids].T@self.force;d.xfrc_applied[:]=0;d.xfrc_applied[self.pelvis,1]=.85*self.weight;self.s.dynamic.update(d);jac=np.zeros((3,self.m.nv));mujoco.mj_jacSubtreeCom(self.m,d,jac,self.pelvis);return jac@d.qacc,R.contact_loads(self.m,d)

def protocol():
 OUT.mkdir(exist_ok=True);p=dict(machine='easteregg2',scope='NEW mechanical-only native-compliance initialization with IDEAL54 reference forces; no neural/pressureclaim',loop_gap_bound_m=.0005,source_constraint_law_unchanged=True,mass_kg=75.,original_force_floors_preserved='Other52>=.05*MIF; selectedpair observed202N idealreference, no imposedpairfloor',pose='First exactgroundedzeroheading family, protectedknees; only small hingecompliance and footcontactpreload/tilt change',hinge_component_window_rad=[-.01,.01],contact_penetration_window_m=[.0001,.004],foot_tilt_component_window_rad=[-.004,.004],harness_support=.85,newinitialheightanchor_explicit=True,max_outer_nfev=40,actual_gate_qacc_max_rad_s2=.05,actual_gate_COM_acc_norm_m_s2=.01,actual_gate_ground_load_error_N=.1,actual_gate_each_foot_load_N=1.,actual_gate_heading_abs_deg=10.,actual_gate_knee_max_rad=-.02,hold_duration_s=5.,fresh_only_after_pass=True,no_geometry_or_physics_edit=True,no_hidden_forces=True,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),Path(S.__file__),P.C.ROBOT/'ah001_neural.xml',P.C.ROBOT/'dynamic_routes.py')});file=OUT/'PREREGISTRATION.json'
 if file.exists():assert json.loads(file.read_text())==p
 else:file.write_text(json.dumps(p,indent=2))
 return p

def solve():
 p=protocol();b=Mechanical();calls=[0]
 def residual(x):
  e,g,a,d,heights,fit=b.inverse_fit(x);calls[0]+=1
  return np.r_[e,g,a,d,1e-8*x]
 lo=np.r_[np.full(6,-10.),.1,np.full(4,-4.)];hi=np.r_[np.full(6,10.),4.,np.full(4,4.)];start=time.perf_counter();startx=np.array(json.loads((OUT/'STATIC_RESULT_default_fd.json').read_text())['variables']);revision=dict(method='Powell derivative-free',maxfev=1000,same_mechanical_case_and_physical_gates=True,reason='Finite-difference inconsistency quantified in FINITE_DIFFERENCE_AUDIT.json',initial= 'Preservedfailed defaultgradientcandidate');(OUT/'NUMERICAL_REVISION.json').write_text(json.dumps(revision,indent=2))
 def objective(x):r=residual(x);return float(r@r)
 sol=minimize(objective,startx,method='Powell',bounds=list(zip(lo,hi)),options={'maxfev':1000,'xtol':1e-7,'ftol':1e-10});e,g,a,d,heights,fit=b.inverse_fit(sol.x);b.x=sol.x;COMacc,loads=b.actual();gaps=b.s.constraints();headings=[]
 for side in ('r','l'):
  axis=b.d.xmat[b.m.body('foot_'+side).id].reshape(3,3)[:,0];headings.append(float(np.rad2deg(np.arctan2(axis[2],axis[0]))))
 knees=list(b.s.dynamic.coordinates(b.d).values());checks=dict(qacc=np.max(np.abs(b.d.qacc))<=.05,COMacc=np.linalg.norm(COMacc)<=.01,ground=abs(loads.sum()-.15*b.weight)<=.1,bothfeet=bool(np.all(loads.sum(axis=1)>1)),loop=max(gaps)<=.0005,heading=max(np.abs(headings))<=10,knees=max(knees)<=-.02,heights=bool(np.all(heights<=0)&np.all(heights>=-.004)));checks={k:bool(v) for k,v in checks.items()};result=dict(accepted=all(checks.values()),checks=checks,solver_success=bool(sol.success),method='Powell derivative-free numericalrevision',nfev=int(sol.nfev),forward_or_inverse_evaluations=calls[0],wall_s=time.perf_counter()-start,variables=sol.x.tolist(),qacc_max=float(np.max(np.abs(b.d.qacc))),COM_acceleration=COMacc.tolist(),contact_loads_N=loads.tolist(),ground_total_N=float(loads.sum()),ground_target_N=.15*b.weight,closure_gaps_m=gaps,heading_deg=headings,knees_rad=knees,foot_min_y_m=heights.tolist(),inverse_force_residual_max_scaled=float(np.max(np.abs(e))),pressure_geometry_unchanged_not_tested=True,neural_controller_not_run=True,physics_unchanged=True)
 np.savez_compressed(OUT/'PRELOAD_STATE.npz',qpos=b.d.qpos.copy(),qvel=b.d.qvel.copy(),qacc=b.d.qacc.copy(),force_N=b.force.copy(),variables=sol.x,root_anchor=b.d.qpos[:3].copy());(OUT/'STATIC_RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--solve',action='store_true');args=parser.parse_args()
 if args.solve:solve()
