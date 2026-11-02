"""Source-inertia pinned-hip fixture and derived numerical rigid-link variant."""
import json,argparse,hashlib,xml.etree.ElementTree as ET,time
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
import cad_posture_audit as P
C=P.C;R=P.R;mujoco=P.mujoco;OUT=P.OUT.parent/'fourbar_numerical_fidelity'
def fmt(x):return ' '.join(format(float(v),'.17g') for v in np.ravel(x))

def build_fixture():
 source=mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));d=mujoco.MjData(source);d.qpos[:]=P.load(P.OUT.parent/'mechanical_stance_repair/STANCES.npz')['qpos'][0];mujoco.mj_forward(source,d);tree=ET.parse(C.ROBOT/'ah001_neural.xml');root=tree.getroot();root.set('model','Pinned source-inertia rightfourbar fixture')
 for tag in ('actuator','tendon','sensor','contact','asset','visual','equality'):
  node=root.find(tag)
  if node is not None:root.remove(node)
 compiler=root.find('compiler');compiler.set('inertiafromgeom','false');keep={'back_r_joint','shank_r_joint','front_r_joint'}
 for parent in list(root.iter()):
  for child in list(parent):
   if child.tag in ('geom','site','camera','light') or child.tag in ('joint','freejoint') and child.get('name') not in keep:parent.remove(child)
 for node in root.iter('body'):
  bid=source.body(node.get('name')).id;pid=source.body_parentid[bid];Rp=d.xmat[pid].reshape(3,3);relative=Rp.T@d.xmat[bid].reshape(3,3);quat=Rotation.from_matrix(relative).as_quat();node.set('pos',fmt(Rp.T@(d.xpos[bid]-d.xpos[pid])));node.set('quat',fmt(np.r_[quat[3],quat[:3]]));inertial=node.find('inertial')
  if inertial is not None:node.remove(inertial)
  ET.SubElement(node,'inertial',mass=str(float(source.body_mass[bid])),pos=fmt(source.body_ipos[bid]),quat=fmt(source.body_iquat[bid]),diaginertia=fmt(source.body_inertia[bid]))
 eq=ET.SubElement(root,'equality');original=ET.parse(C.ROBOT/'ah001_neural.xml').getroot().find('equality');eq.append(next(e for e in original if e.get('name')=='knee_loop_r'));text=ET.tostring(root,encoding='unicode');fixture=mujoco.MjModel.from_xml_string(text);(OUT/'PINNED_SOURCE_FIXTURE.xml').write_text(text,encoding='utf-8')
 proof=dict(body_mass_exact=bool(np.array_equal(source.body_mass,fixture.body_mass)),body_inertia_exact=bool(np.array_equal(source.body_inertia,fixture.body_inertia)),body_ipos_exact=bool(np.array_equal(source.body_ipos,fixture.body_ipos)),total_mass_kg=float(fixture.body_mass.sum()),fixture_NV=fixture.nv,source_mass_sha256=hashlib.sha256(source.body_mass.tobytes()).hexdigest(),source_inertia_sha256=hashlib.sha256(source.body_inertia.tobytes()).hexdigest(),external_fixture='All but three rightknee hinges pinned. No gravity/contact/muscle/NNrobotstanding claim; bodymass/inertia arrays retained exactly.',source_xml_sha256=hashlib.sha256((C.ROBOT/'ah001_neural.xml').read_bytes()).hexdigest())
 (OUT/'SOURCE_INERTIA_PROOF.json').write_text(json.dumps(proof,indent=2));return text,proof

class Fixture:
 def __init__(self,text,tc=None,imp=None,dt=.002,iterations=100):
  self.m=mujoco.MjModel.from_xml_string(text);self.m.opt.timestep=dt;self.m.opt.iterations=iterations
  if tc is not None:self.m.eq_solref[0]=[tc,1.]
  if imp is not None:self.m.eq_solimp[0,:2]=imp
  self.d=mujoco.MjData(self.m);self.c=np.array([self.m.jnt_axis[self.m.joint('back_r_joint').id,2],self.m.jnt_axis[self.m.joint('shank_r_joint').id,2],0.]);self.W=null_space(self.c[None,:]);self.q0=self.d.qpos.copy();self.base();self.J=self.d.efc_J.reshape(self.d.nefc,self.m.nv).copy();self.rows=np.flatnonzero(self.d.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY));self.sourceK=self.localK()
 def base(self):self.d.qvel[:]=0;self.d.qfrc_applied[:]=0;mujoco.mj_forward(self.m,self.d)
 def inverse(self,x,load):
  self.d.qpos[:]=self.q0+self.W@x;self.d.qvel[:]=0;mujoco.mj_forward(self.m,self.d);self.d.qacc[:]=0;mujoco.mj_inverse(self.m,self.d);tau=np.array([load,0.,0.]);return self.W.T@(self.d.qfrc_inverse-tau)
 def gap(self):return float(np.linalg.norm(self.d.efc_pos[self.rows]))
 def localK(self):
  h=1e-6;return np.column_stack([(self.inverse(np.eye(2)[i]*h,0)-self.inverse(-np.eye(2)[i]*h,0))/(2*h) for i in range(2)])
 def equilibrium(self,load):
  residual0=self.inverse(np.zeros(2),load);guess=np.linalg.lstsq(self.sourceK,-residual0,rcond=None)[0];prediction_q=self.q0+self.W@guess;self.d.qpos[:]=prediction_q;self.base();prediction_gap=self.gap();seed=np.clip(guess,-.5,.5);sol=least_squares(lambda x:self.inverse(x,load),seed,bounds=(-np.ones(2)*.6,np.ones(2)*.6),max_nfev=80,ftol=1e-13,gtol=1e-13,xtol=1e-13,diff_step=1e-4);res=self.inverse(sol.x,load);required=self.d.qfrc_inverse.copy();tau=np.array([load,0.,0.]);rho=float(self.c@(required-tau)/(self.c@self.c));applied=tau+self.c*rho;inverse_force=self.d.efc_force.copy();q=self.d.qpos.copy();gap=self.gap();self.d.qfrc_applied[:]=applied;mujoco.mj_forward(self.m,self.d);acc=self.d.qacc.copy();equilibrium_force=self.d.efc_force.copy();J=self.d.efc_J.reshape(self.d.nefc,3).copy();M=np.zeros((3,3));mujoco.mj_fullM(self.m,M,self.d.qM);balance=M@acc+self.d.qfrc_bias-self.d.qfrc_passive-self.d.qfrc_constraint-applied;self.d.qacc[:]=acc;mujoco.mj_inverse(self.m,self.d);parity=float(np.max(np.abs(self.d.qfrc_inverse-applied)))
  row=dict(load_back_native_Nm=load,static_residual_Nm=res.tolist(),solution_success=bool(sol.success),linear_predicted_gap_m=prediction_gap,actual_gap_m=gap,hinge_delta_rad=(q-self.q0).tolist(),fixture_hold_generalized_torque_Nm=(self.c*rho).tolist(),qacc_max=float(np.max(np.abs(acc))),forward_inverse_error_Nm=parity,EOM_balance_max_Nm=float(np.max(np.abs(balance))),constraint_force_N=equilibrium_force.tolist(),constraint_joint_torque_Nm=(J.T@equilibrium_force).tolist(),constraint_J_rank=int(np.linalg.matrix_rank(J,tol=1e-10)),solver_iterations=int(getattr(self.d,'solver_iter',-1)),physical_hardware_calibration=False)
  return row,q,applied
 def details(self):
  self.d.qpos[:]=self.q0;self.base();M=np.zeros((3,3));mujoco.mj_fullM(self.m,M,self.d.qM);J=self.d.efc_J.reshape(self.d.nefc,3);A=J@np.linalg.solve(M,J.T)
  return dict(eq_solref=self.m.eq_solref.tolist(),eq_solimp=self.m.eq_solimp.tolist(),hinge_mass_matrix=M.tolist(),constraint_J=J.tolist(),J_singular_values=np.linalg.svd(J,compute_uv=False).tolist(),constraint_inverse_inertia=A.tolist(),native_regularizer_R=self.d.efc_R.tolist(),native_diagApprox=self.d.efc_diagApprox.tolist(),local_inverse_stiffness_Nm_per_rad=self.sourceK.tolist(),loop_anchor_leverarms_m=[float(np.linalg.norm(self.m.eq_data[0,:3])),float(np.linalg.norm(self.m.eq_data[0,3:6]))],damping=self.m.dof_damping.tolist(),armature=self.m.dof_armature.tolist())

def run():
 OUT.mkdir(exist_ok=True);text,proof=build_fixture();source=Fixture(text);details=source.details();rows=[]
 for load in (0.,1.,8.,9.):row,_,_=source.equilibrium(load);rows.append(row)
 target=.0001;worst=max(row['linear_predicted_gap_m'] for row in rows);tc=.004;source_d=.95;factor=max(1.,worst/target);imp=1.-(1.-source_d)*(.01/tc)**2/factor
 if imp>.9999:tc=.002;imp=1.-(1.-source_d)*(.01/tc)**2/factor
 imp=float(min(.9999,max(.95,imp)));dt=min(.002,tc/2.);selection=dict(reason='Scale staticstiffness~1/[(1−impedance)*timeconst²] from measured source linearcompliance at9Nm to100um target; constantimpedance removesdistance-dependentsoftlayer.',source_linear_worst_gap_m=worst,target_gap_m=target,required_stiffening_factor=factor,timeconst_s=tc,constant_impedance=imp,timestep_s=dt,refsafe_preserved=True,not_hardware_calibration=True)
 (OUT/'DERIVED_VARIANT.json').write_text(json.dumps(selection,indent=2));variants={}
 for tag,h,iterations in [('selected',dt,100),('iterations500',dt,500),('halfstep',dt/2.,500)]:
  f=Fixture(text,tc,imp,h,iterations);out=[]
  for load in (0.,1.,8.,9.):row,_,_=f.equilibrium(load);out.append(row)
  variants[tag]=dict(details=f.details(),cases=out)
 result=dict(scope='Pinnedsource-inertiafixture only; no robottrial, neuralmodel orBPAchanges',proof=proof,source_details=details,source_cases=rows,derived=selection,variants=variants,primary_sources=['https://mujoco.readthedocs.io/en/2.3.7/modeling.html#solver-parameters','https://mujoco.readthedocs.io/en/2.3.7/computation.html#parameters']);(OUT/'FIXTURE_STATIC_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(source=rows,derived=selection,variants={key:v['cases'] for key,v in variants.items()}),indent=2));return result

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--static-audit',action='store_true');a=parser.parse_args()
 if a.static_audit:run()
