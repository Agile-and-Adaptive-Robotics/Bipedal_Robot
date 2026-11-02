"""Separate geometry-only CAD standing initialization; original models unchanged."""
import json,argparse,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from scipy.optimize import brentq,linprog
import cad_posture_audit as P
C=P.C;R=P.R;mujoco=P.mujoco
OUT=P.OUT.parent/'mechanical_stance_repair'

def sole_frame(m,side):
 d=mujoco.MjData(m);mujoco.mj_forward(m,d);bid=m.body('foot_'+side).id;cloud=[]
 for gid in range(m.ngeom):
  if not m.geom(gid).name.startswith('foot_'+side+'_'):continue
  mid=m.geom_dataid[gid]
  if mid>=0:
   v=m.mesh_vert[m.mesh_vertadr[mid]:m.mesh_vertadr[mid]+m.mesh_vertnum[mid]];world=v@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid];cloud.append((world-d.xpos[bid])@d.xmat[bid].reshape(3,3))
 v=np.vstack(cloud);h=ConvexHull(v);candidates=[];seen=set()
 for equation in h.equations:
  n=equation[:3]
  if n[1]>-.866:continue
  key=tuple(np.round(n,4))
  if key in seen:continue
  seen.add(key);rot,_=Rotation.align_vectors([[0,-1,0]],[n]);M=rot.as_matrix();forward=M[:,0];heading=np.arctan2(forward[2],forward[0]);M=Rotation.from_rotvec([0,heading,0]).as_matrix()@M;points=v@M.T;minimum=points[:,1].min();touch=points[points[:,1]<minimum+.0003][:,[0,2]]
  if len(touch)<3:continue
  try:area=ConvexHull(touch).volume
  except Exception:continue
  candidates.append((area,M,float(minimum),np.ptp(touch,axis=0)))
 area,M,minimum,span=max(candidates,key=lambda row:row[0]);return dict(matrix=M,local_vertices=v,minimum=minimum,area=float(area),span=span)

class Stance:
 def __init__(self):
  self.native=mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));self.caps=-self.native.actuator_gainprm[:,0].copy();self.m,self.dynamic=P.model();self.d=mujoco.MjData(self.m);self.soles={s:sole_frame(self.m,s) for s in ('r','l')};self.cfg=json.loads((P.S.N.BASE/'ASSUMED_PAIR_CONFIGURATION.json').read_text());self.bank=C.RoutedBPABank(C.PAIR,self.cfg,self.m.opt.timestep);self.ids=[self.m.actuator(n).id for n in C.PAIR];self.tids=self.m.actuator_trnid[:,0];a=P.load(P.S.ROOT/'mixed/mixed/ACTUATOR_BRIDGE.npz');self.pairT=self.caps[self.ids]*a['neural_output'][-1,self.ids]
 def relative_pose(self,angles):
  m,d=self.m,self.d;d.qpos[:]=m.qpos0;d.qpos[:3]=[0,.9,0];d.qpos[3:7]=[1,0,0,0];d.qvel[:]=0;d.qfrc_applied[:]=0;d.xfrc_applied[:]=0
  for side,angle in zip(('r','l'),angles):R.knee_pose(m,d,side,float(angle))
  mujoco.mj_kinematics(m,d)
  for side in ('r','l'):
   hip=m.jnt_qposadr[m.joint('thigh_'+side+'_joint').id];foot=m.body('foot_'+side).id;cur=d.xmat[foot].reshape(3,3);quat=Rotation.from_matrix(self.soles[side]['matrix']@cur.T).as_quat();d.qpos[hip:hip+4]=np.r_[quat[3],quat[:3]]
  mujoco.mj_kinematics(m,d);return np.array([d.xpos[m.body('foot_'+s).id,1]+self.soles[s]['minimum'] for s in ('r','l')])
 def pose(self,angles,penetration=.0002):
  y=self.relative_pose(angles);self.d.qpos[1]-=float(y.mean()+penetration);self.dynamic.update(self.d);return self.d.qpos.copy()
 def constraints(self):
  m,d=self.m,self.d;gaps=[]
  for side in ('r','l'):
   e=m.equality('knee_loop_'+side).id;one,two=m.eq_obj1id[e],m.eq_obj2id[e];gaps.append(float(np.linalg.norm(d.xpos[one]+d.xmat[one].reshape(3,3)@m.eq_data[e,:3]-d.xpos[two]-d.xmat[two].reshape(3,3)@m.eq_data[e,3:6])))
  return gaps
 def statics(self):
  m,d=self.m,self.d;dyn=self.dynamic;pelvis=m.body('pelvis').id;weight=m.body_mass.sum()*9.81;F=np.array([0,.85*weight,0]);tau=np.zeros(3);Qrig=np.zeros(m.nv);mujoco.mj_applyFT(m,d,F,tau,d.xipos[pelvis],pelvis,Qrig);N=R.anatomical_knee_tangent(m,d);T=np.zeros((m.nv,6+N.shape[1]));T[:6,:6]=np.eye(6);T[6:,6:]=N;A=T.T@dyn.moment_arms[self.tids].T;positions=[];bodies=[];mus=[]
  for c in d.contact:
   names=[m.geom(int(g)).name for g in (c.geom1,c.geom2)]
   if 'ground' in names and any(n.startswith('foot_') for n in names):positions.append(c.pos.copy());bodies.append(int(m.geom_bodyid[c.geom2 if names[0]=='ground' else c.geom1]));mus.append(float(c.friction[0])/np.sqrt(2))
  columns=[A]
  for point,bid in zip(positions,bodies):
   jac=np.zeros((3,m.nv));ang=np.zeros_like(jac);mujoco.mj_jac(m,d,jac,ang,point,bid);columns.append(T.T@jac.T)
  E=np.column_stack(columns);bounds=[(.05*cap,cap) for cap in self.caps]
  for j,idx in enumerate(self.ids):bounds[idx]=(float(self.pairT[j]),float(self.pairT[j]))
  bounds.extend([(None,None),(0,weight),(None,None)]*len(positions));ineq=[]
  for i,mu in enumerate(mus):
   for axis in (0,2):
    for sign in (-1,1):row=np.zeros(E.shape[1]);row[m.nu+3*i+axis]=sign;row[m.nu+3*i+1]=-mu;ineq.append(row)
  lp=linprog(np.r_[1/self.caps,np.zeros(3*len(positions))],A_ub=np.array(ineq),b_ub=np.zeros(len(ineq)),A_eq=E,b_eq=T.T@(d.qfrc_bias-Qrig),bounds=bounds,method='highs') if positions else None
  com=d.subtree_com[pelvis];point=d.xipos[pelvis];cop=(weight*com[[0,2]]-F[1]*point[[0,2]])/(weight-F[1]);margin=P.inside_margin(np.asarray(positions)[:,[0,2]],cop) if len(positions)>=3 else np.nan;press=[]
  for j,name in enumerate(C.PAIR):p,flag=C.inverse_pressure(float(self.pairT[j]),float(d.ten_length[self.tids[self.ids[j]]]+self.cfg['actuators'][name]['constant_remote_path_offset_m']),self.bank.parameters[j]);press.append(dict(name=name,kpa=float(p),flag=int(flag)))
  result=dict(ideal_constraint_capacity_feasible=bool(lp is not None and lp.success),capacity_message=lp.message if lp is not None else 'No geometriccontact',required_static_CoP_xz=cop.tolist(),CoP_geometric_margin_m=float(margin),contacts_per_foot={s:sum(b==m.body('foot_'+s).id for b in bodies) for s in ('r','l')},static_pair_pressure=press,COM=com.tolist(),harness_force_N=F.tolist(),harness_torque_Nm=tau.tolist())
  if lp is not None and lp.success:
   force=lp.x[:m.nu];d.qfrc_applied[:]=-dyn.length_jacobian[self.tids].T@force;d.xfrc_applied[pelvis,:3]=F;dyn.update(d);jac=np.zeros((3,m.nv));mujoco.mj_jacSubtreeCom(m,d,jac,pelvis);result.update(ideal_reduced_equilibrium_error=float(np.max(np.abs(E@lp.x-T.T@(d.qfrc_bias-Qrig)))),actual_mj_forward_qacc_max=float(np.max(np.abs(d.qacc))),actual_mj_forward_COM_acceleration=(jac@d.qacc).tolist(),actual_mj_forward_ground_N=float(R.contact_loads(m,d).sum()),actual_ground_target_N=float(weight-F[1]),LP_force_N=force.tolist())
  return result

def offline():
 OUT.mkdir(exist_ok=True);source=Path(__file__);registration=dict(machine='easteregg2',scope='Separateinitialization/reference only; offlineexactloop/soleplane/forceaudit, no behavior',model_unchanged=True,mass_unchanged=True,heading_goal_body_axis_deg=0.,maximum_acceptable_heading_deg=10.,right_knee_samples_rad=np.linspace(-.55,-.08,9).tolist(),left_knee_solve_bounds_rad=[-.55,-.08],sole_penetration_m=.0002,root_and_lumbar_world_orientation='Identity',harness_feedforward=.85,new_harness_height_anchor='Newexplicitgroundedinitialpose, not weightweaning comparison',old_results_unchanged=True,pair_weights_unchanged=P.S.N.U.weights(),pair_reference_force_N='Frozenobserved202N output for staticcapacity only, not neuralprediction',other52_force_bounds='Original .05*MIF..MIF referencecontext',no_physics_changes=True,no_hardware_ROM_certification=True,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest());path=OUT/'PREREGISTRATION.json'
 if path.exists():assert json.loads(path.read_text())==registration
 else:path.write_text(json.dumps(registration,indent=2))
 system=Stance();rows=[];states=[];soles={side:dict(matrix=value['matrix'].tolist(),support_area_m2=value['area'],support_AP_ML_span_m=value['span'].tolist()) for side,value in system.soles.items()}
 for right in registration['right_knee_samples_rad']:
  def mismatch(left):y=system.relative_pose([right,left]);return float(y[0]-y[1])
  grid=np.linspace(-.55,-.08,17);values=[mismatch(left) for left in grid];roots=[]
  for i in range(len(grid)-1):
   if values[i]*values[i+1]<=0:roots.append(brentq(mismatch,float(grid[i]),float(grid[i+1]),xtol=1e-12))
  if not roots:rows.append(dict(right_knee_rad=right,geometric_height_root=False,height_residual_range_m=[min(values),max(values)]));continue
  left=min(roots,key=lambda value:abs(value-right));q=system.pose([right,left]);m,d=system.m,system.d;details=system.statics();row=dict(right_knee_rad=right,left_knee_rad=left,geometric_height_root=True,closure_gaps_m=system.constraints(),root_height_m=float(q[1]),mass_kg=float(m.body_mass.sum()),statics=details);rows.append(row);states.append(q);print(json.dumps({k:v for k,v in row.items() if k!='statics'}),flush=True);print(json.dumps({k:v for k,v in details.items() if k!='LP_force_N'}),flush=True)
  (OUT/'OFFLINE_RESULTS.json').write_text(json.dumps(dict(registration=registration,sole_frames=soles,cases=rows),indent=2))
 result=dict(registration=registration,sole_frames=soles,cases=rows);(OUT/'OFFLINE_RESULTS.json').write_text(json.dumps(result,indent=2));np.savez_compressed(OUT/'STANCES.npz',qpos=np.array(states));return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--offline',action='store_true');a=p.parse_args()
 if a.offline:offline()
