"""Opt-in offline supported, in-place stepping reference preparation.

Stance feet are fixed in world coordinates. C2 swing joint bumps remain inside
the frozen source-IK envelopes. Native guide springs, original body geometry,
and source force laws are retained. Nothing is injected into a live SNS run.
"""
import os,sys,json,argparse,hashlib
from time import perf_counter
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'functional_gait2392')]
import numpy as np,mujoco
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from scipy.interpolate import PchipInterpolator
from mechanics import CoupledArms,JOINTS
from contact_standing_solver import _support_regions,_external

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def smooth(x):
 x=np.clip(x,0.,1.);return x*x*x*(10.+x*(-15.+6.*x))
def bump(x):return 64.*x**3*(1.-x)**3

class NativeFootPitchTeacher:
 """Bounded ankle targets for calibrated world-foot pitch; hip/knee unchanged."""
 def __init__(self,model,keyframes,spec):
  if not np.allclose(model.opt.gravity/np.linalg.norm(model.opt.gravity),[0.,0.,-1.]):raise ValueError('Foot-pitch adapter requires the explicit captured z-up human frame')
  self.model=model;self.feet=[model.body('calcn_'+s).id for s in ('r','l')];self.joints=[model.joint('ankle_angle_'+s).id for s in ('r','l')];self.rows=model.jnt_dofadr[self.joints];self.addresses=model.jnt_qposadr[self.joints];self.bounds=np.array([spec['joints']['ankle_angle_'+s]['native_bounds_rad'] for s in ('r','l')]);self.pitch=[];self.jp=np.zeros((3,model.nv));self.jr=np.zeros_like(self.jp);d=mujoco.MjData(model)
  for frame in keyframes:
   d.qpos[:]=frame['qpos'];mujoco.mj_fwdPosition(model,d);self.pitch.append(np.array([self.angle(d.xmat[b].reshape(3,3)[:,0]) for b in self.feet]))
  self.updates=0;self.max_error=0.;self.max_truncation=0.;self.bound_hits=0;self.min_sensitivity=np.inf

 @staticmethod
 def angle(axis):return float(np.arctan2(axis[2],np.hypot(axis[0],axis[1])))

 def solve(self,live,index,blend,rate,qref,vref):
  q=qref.copy();v=vref.copy()
  for side,(body,row,address) in enumerate(zip(self.feet,self.rows,self.addresses)):
   axis=live.xmat[body].reshape(3,3)[:,0];h=np.hypot(axis[0],axis[1]);gradient=np.array([-axis[2]*axis[0]/h,-axis[2]*axis[1]/h,h]);angular_gradient=np.cross(axis,gradient);mujoco.mj_jac(self.model,live,self.jp,self.jr,live.xpos[body],body);sensitivity=float(angular_gradient@self.jr[:,row])
   if abs(sensitivity)<.2:raise ValueError('Ankle cannot reliably control world-foot pitch at this pose')
   target=self.pitch[index][side]+blend*(self.pitch[index+1][side]-self.pitch[index][side]);target_rate=rate*(self.pitch[index+1][side]-self.pitch[index][side]);error=target-self.angle(axis);requested=live.qpos[address]+error/sensitivity;q[address]=np.clip(requested,*self.bounds[side]);v[row]=live.qvel[row]+(target_rate-angular_gradient@self.jr@live.qvel)/sensitivity;self.max_error=max(self.max_error,abs(error));self.max_truncation=max(self.max_truncation,abs(requested-q[address]));self.bound_hits+=int(requested!=q[address]);self.min_sensitivity=min(self.min_sensitivity,abs(sensitivity))
  self.updates+=1;return q,v

 def diagnostics(self):
  return dict(updates=self.updates,keyframe_pitch_deg=np.rad2deg(self.pitch).tolist(),maximum_actual_pitch_error_deg=float(np.rad2deg(self.max_error)),maximum_desired_ankle_bound_truncation_deg=float(np.rad2deg(self.max_truncation)),bound_hit_updates=self.bound_hits,minimum_ankle_pitch_sensitivity=self.min_sensitivity,scope='Original calibrated world-foot pitch using actual foot orientation/angular velocity. Only desired ankle targets change; hip/knee posture targets remain. Desired ankle saturation is explicit; no live qpos/qvel clamp. Frozen physical contact/clearance/ROM gates decide acceptance.')

class NativeCartesianTeacher:
 """Bounded desired-leg IK using actual root pose; never changes live state."""
 def __init__(self,model,keyframes,spec,stance_posture=False,preserve_loaded_followers=False,tangential_clearance=None):
  self.model=model;self.data=mujoco.MjData(model);self.arms=CoupledArms(model)
  self.feet=[model.body('calcn_'+side).id for side in ('r','l')]
  self.columns=[np.array([JOINTS.index(name+'_'+side) for name in ('hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle')]) for side in ('r','l')]
  self.bounds=[];self.position=[];self.rotation=[];self.gaps=[];self.last=None;self.stance_posture=stance_posture
  self.preserve_loaded_followers=preserve_loaded_followers;self.depaddr=np.array([model.jnt_qposadr[model.dof_jntid[k]] for k in self.arms.dependent]);self.max_loaded_offset=0.;self.max_loaded_offset_velocity=0.
  if preserve_loaded_followers and not np.allclose(model.opt.gravity/np.linalg.norm(model.opt.gravity),[0.,0.,-1.]):raise ValueError('Loaded adapter requires the explicit captured z-up human frame')
  if tangential_clearance is not None and (not preserve_loaded_followers or not np.isfinite(tangential_clearance) or tangential_clearance<=.001):raise ValueError('Tangential clearance requires loaded FK and an explicit margin above1mm')
  self.tangential_clearance=tangential_clearance;self.max_planned_loaded_tangent=0.;self.max_planned_loaded_rotation=0.
  self.jp=np.zeros((3,model.nv));self.jr=np.zeros_like(self.jp)
  self.max_position_error=0.;self.max_rotation_error=0.;self.bound_hits=0;self.updates=0;self.max_desired_speed=0.;self.first_failure=None;self.max_unweighted_vertical_error=0.;self.max_knee_posture_error=0.;self.max_normal_posture_task_error=0.
  for columns in self.columns:
   lower=[];upper=[]
   for column in columns:
    joint=int(self.arms.jids[column]);name=JOINTS[column]
    lo,hi=spec['joints'][name]['native_bounds_rad'] if name in spec['joints'] else (max(model.jnt_range[joint,0],np.deg2rad(-10.)),min(model.jnt_range[joint,1],np.deg2rad(10.)))
    if hi<=lo:raise ValueError('Missing native desired-IK bounds for '+name)
    if name.startswith('knee_angle_'):hi=min(hi,np.deg2rad(-2.))
    lower.append(lo);upper.append(hi)
   self.bounds.append((np.array(lower),np.array(upper)))
  for frame in keyframes:
   self.data.qpos[:]=frame['qpos']
   if not self.preserve_loaded_followers:self.arms.project(self.data)
   mujoco.mj_fwdPosition(model,self.data)
   self.position.append(self.data.xpos[self.feet].copy());self.rotation.append(self.data.xmat[self.feet].reshape(2,3,3).copy())
   regions,_=_support_regions(model,self.data);self.gaps.append(np.array([min(r['gap_m'] for r in regions[2*s:2*s+2]) for s in range(2)]))

 def target_kinematics(self,index,side,blend,rate):
  delta=self.position[index+1][side]-self.position[index][side];g0=self.gaps[index][side];g1=self.gaps[index+1][side];horizontal_blend=blend;horizontal_rate=rate
  if self.tangential_clearance is not None:
   margin=self.tangential_clearance
   if g0<=.001 and g1>margin:
    start=(margin-g0)/(g1-g0);u=np.clip((blend-start)/(1.-start),0.,1.);horizontal_blend=float(smooth(u));horizontal_rate=30.*u*u*(1.-u)**2*rate/(1.-start) if 0.<u<1. else 0.
   elif g1<=.001 and g0>margin:
    finish=(g0-margin)/(g0-g1);u=np.clip(blend/finish,0.,1.);horizontal_blend=float(smooth(u));horizontal_rate=30.*u*u*(1.-u)**2*rate/finish if 0.<u<1. else 0.
   elif min(g0,g1)<=.001 and max(g0,g1)>.005:raise ValueError('Tangential margin exceeds declared swing-keyframe clearance')
  target=self.position[index][side]+blend*delta;angular=Rotation.from_matrix(self.rotation[index+1][side]@self.rotation[index][side].T).as_rotvec();orientation=Rotation.from_rotvec(horizontal_blend*angular).as_matrix()@self.rotation[index][side];velocity=np.r_[rate*delta,horizontal_rate*angular]
  if self.tangential_clearance is not None:target[:2]=self.position[index][side][:2]+horizontal_blend*delta[:2];velocity[:2]=horizontal_rate*delta[:2]
  return target,orientation,velocity,float(g0+blend*(g1-g0))

 def solve(self,live,index,blend,rate,desired_q,desired_v):
  m,d=self.model,self.data;q=desired_q.copy();q[:6]=live.qpos[:6];v=desired_v.copy()
  if self.preserve_loaded_followers:
   # Retain compliant path-point displacement and velocity in scratch FK.
   # The ideal projection defines the offset; it never replaces live state.
   d.qpos[:]=live.qpos;self.arms.project(d);loaded_offset=live.qpos[self.depaddr]-d.qpos[self.depaddr];self.arms.update(live);loaded_velocity=live.qvel[self.arms.dependent]-self.arms.E[self.arms.dependent]@live.qvel[self.arms.rows];self.max_loaded_offset=max(self.max_loaded_offset,float(np.max(abs(loaded_offset))));self.max_loaded_offset_velocity=max(self.max_loaded_offset_velocity,float(np.max(abs(loaded_velocity))))
  for side,(body,columns,bounds) in enumerate(zip(self.feet,self.columns,self.bounds)):
   address=self.arms.addresses[columns];rows=self.arms.rows[columns];lo,hi=bounds
   target,orientation,target_velocity,reference_gap=self.target_kinematics(index,side,blend,rate)
   if reference_gap<=.001:self.max_planned_loaded_tangent=max(self.max_planned_loaded_tangent,float(np.linalg.norm(target_velocity[:2])));self.max_planned_loaded_rotation=max(self.max_planned_loaded_rotation,float(np.linalg.norm(target_velocity[3:])))
   stance_weight=float(smooth((.001-reference_gap)/.001)) if self.stance_posture else 0.;vertical_weight=1.-stance_weight;reference_knee=float(desired_q[address[3]])
   def residual(x):
    q[address]=x;d.qpos[:]=q;self.arms.project(d)
    if self.preserve_loaded_followers:d.qpos[self.depaddr]+=loaded_offset
    mujoco.mj_fwdPosition(m,d)
    position=d.xpos[body]-target;rotation=Rotation.from_matrix(d.xmat[body].reshape(3,3)@orientation.T).as_rotvec()
    if self.stance_posture:return np.r_[position[:2]/.0002,rotation/.001,vertical_weight*position[2]/.0002-stance_weight*(x[3]-reference_knee)/.01]
    return np.r_[position/.0002,rotation/.001]
   guess=q[address] if self.last is None else self.last[address]+desired_q[address]-self.previous_reference[address]
   fit=least_squares(residual,np.clip(guess,lo+1e-10,hi-1e-10),bounds=(lo,hi),max_nfev=12,xtol=1e-9,ftol=1e-9,gtol=1e-9)
   residual(fit.x);position_error=d.xpos[body]-target;normal_task=vertical_weight*position_error[2]/.0002-stance_weight*(fit.x[3]-reference_knee)/.01;self.max_unweighted_vertical_error=max(self.max_unweighted_vertical_error,float(abs(position_error[2])));self.max_normal_posture_task_error=max(self.max_normal_posture_task_error,float(abs(normal_task))) if self.stance_posture else 0.;position_error[2]*=0. if self.stance_posture else 1.;self.max_position_error=max(self.max_position_error,float(np.max(abs(position_error))));self.max_knee_posture_error=max(self.max_knee_posture_error,float(stance_weight*abs(fit.x[3]-reference_knee)))
   self.max_rotation_error=max(self.max_rotation_error,float(np.linalg.norm(Rotation.from_matrix(d.xmat[body].reshape(3,3)@orientation.T).as_rotvec())))
   if self.first_failure is None and (np.max(abs(position_error))>.0002 or (self.stance_posture and abs(normal_task)>.001) or np.linalg.norm(Rotation.from_matrix(d.xmat[body].reshape(3,3)@orientation.T).as_rotvec())>.001):self.first_failure=dict(time_s=float(live.time),side=('r','l')[side],position_error_world_m=(d.xpos[body]-target).tolist(),normal_posture_task_error=float(normal_task),stance_weight=stance_weight,root_qpos_native=live.qpos[:6].tolist(),desired_leg_angles_deg=np.rad2deg(fit.x).tolist(),nfev=int(fit.nfev),optimizer_status=int(fit.status),scope='Bounded local target solve did not meet pose criterion; not a global kinematic infeasibility proof.')
   self.bound_hits+=int(np.any(np.minimum(fit.x-lo,hi-fit.x)<1e-6))
   self.arms.update(d);mujoco.mj_jac(m,d,self.jp,self.jr,d.xpos[body],body)
   J=np.vstack([self.jp@self.arms.E[:,columns],self.jr@self.arms.E[:,columns]])
   root_velocity=np.vstack([self.jp[:,:6],self.jr[:,:6]])@live.qvel[:6]
   rhs=target_velocity-root_velocity
   if self.preserve_loaded_followers:rhs-=np.vstack([self.jp[:,self.arms.dependent],self.jr[:,self.arms.dependent]])@loaded_velocity
   if self.stance_posture:
    u=np.clip((.001-reference_gap)/.001,0.,1.);weight_rate=30.*u*u*(1.-u)**2*(-rate*(self.gaps[index+1][side]-self.gaps[index][side])/.001) if 0.<u<1. else 0.
    normal_J=vertical_weight*J[2]/.0002-stance_weight*np.eye(6)[3]/.01;normal_rhs=vertical_weight*rhs[2]/.0002-stance_weight*desired_v[rows[3]]/.01+weight_rate*((d.xpos[body,2]-target[2])/.0002+(fit.x[3]-reference_knee)/.01)
    J=np.vstack([J[:2]/.0002,J[3:]/.001,normal_J]);rhs=np.r_[rhs[:2]/.0002,rhs[3:]/.001,normal_rhs];v[rows]=np.linalg.lstsq(J,rhs,rcond=1e-5)[0]
   else:v[rows]=np.linalg.lstsq(J,rhs,rcond=1e-5)[0]
   self.max_desired_speed=max(self.max_desired_speed,float(np.max(abs(v[rows]))))
  if self.preserve_loaded_followers:
   d.qpos[:]=q;self.arms.project(d);d.qpos[self.depaddr]+=loaded_offset;q=d.qpos.copy();self.arms.update(d);v[self.arms.dependent]=self.arms.E[self.arms.dependent]@v[self.arms.rows]+loaded_velocity
  self.last=q.copy();self.previous_reference=desired_q.copy();self.updates+=1
  return q,v

 def diagnostics(self):
  return dict(updates=self.updates,maximum_target_position_error_m=self.max_position_error,maximum_normal_posture_task_error=self.max_normal_posture_task_error,maximum_unweighted_vertical_error_m=self.max_unweighted_vertical_error,maximum_weighted_knee_posture_error_rad=self.max_knee_posture_error,bound_hit_updates=self.bound_hits,maximum_target_rotation_error_rad=self.max_rotation_error,maximum_desired_leg_speed_rad_s=self.max_desired_speed,accepted_target_IK=bool(self.max_position_error<=.0002 and self.max_rotation_error<=.001 and self.max_normal_posture_task_error<=.001),first_failure=self.first_failure,stance_posture=self.stance_posture,preserve_loaded_followers=self.preserve_loaded_followers,maximum_preserved_loaded_offset_m=self.max_loaded_offset,maximum_preserved_loaded_offset_velocity_mps=self.max_loaded_offset_velocity,tangential_clearance_margin_m=self.tangential_clearance,maximum_planned_low_clearance_tangential_speed_mps=self.max_planned_loaded_tangent,maximum_planned_low_clearance_angular_speed_radps=self.max_planned_loaded_rotation,scope='Desired IK only; optional clearance margin finishes tangential/orientation motion before lowering and delays it until after lift. Full loaded FK/velocity retained; no live state changes. Physical contact/ROM gates unchanged.')

def prepare(model_folder,spec_path,out,period=2.4,frames=61,support=.5,lift_clearance=0.,cartesian=False,cartesian_height=None,toe_fraction=.5):
 out.mkdir(parents=True,exist_ok=True)
 manifest=dict(format='aarl_supported_step_teacher_v1',source_model_sha256=sha(model_folder/'runtime_model.mjb'),source_envelope_sha256=sha(spec_path),source_code_sha256=sha(__file__),period_s=period,frames=frames,support_fraction=support,duty=.65,stance_gap_m=-.001,swing_native_peak_deg=[-18.,-55.,10.],other_leg_axis_preparation_limit_deg=10.,acceptance=dict(stance_speed_max_mps=.005,stance_body_position_error_max_m=.0002,stance_rotation_error_max_rad=.001,mesh_penetration_max_m=.002,central_swing_clearance_min_m=.01,knee_extension_max_deg=-2.),scope='Declared engineering stepping-in-place teacher, not measured walking. Kinematic and inverse feasibility only; original IK and accepted SNS remain unchanged. Source6-axis envelopes are acceptance limits;10deg other-axis caps are explicit preparation bounds, not physiological claims.')
 manifest['swing_minimum_mesh_clearance_plateau_m']=lift_clearance;manifest['swing_clearance_profile']='C2 ramp over first/last20percent of swing,plateau in between; actual inverse kinematics,not a live trajectory clamp';manifest['swing_pose_tracking']='Soft position/orientation goal with hard actual mesh-clearance penalty; stance pose stays fixed';manifest['root_height_reference_peak_m']=.005 if lift_clearance and not cartesian else 0.;manifest['dense_validation_dt_s']=.002;manifest['interpolation']='Periodic tiled componentwise PCHIP,exact polynomial equality followers';manifest['acceptance']['periodic_endpoint_max_rad_or_m']=1e-5
 manifest['quasistatic_distal_toe_load_fraction']=toe_fraction
 if cartesian:manifest['cartesian_target']=dict(forward_excursion_m=.04,vertical_lift_m=.025 if cartesian_height is None else cartesian_height,flat_foot=True,stance='Exact fixed world body transform',swing='C2x-bump and C2z-plateau,returns to same stance transform; original joint envelopes constrain IK')
 mp=out/'MANIFEST.json'
 if mp.exists() and json.loads(mp.read_text())!=manifest:raise RuntimeError('Existing teacher manifest differs; select a new checkpoint folder')
 mp.write_text(json.dumps(manifest,indent=2));m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);ar=CoupledArms(m);spec=json.loads(spec_path.read_text());bw=float(m.body_mass.sum()*9.81);base=m.qpos0.copy();base[:6]=0.;base[m.joint('pelvis_ty').qposadr[0]]=.94;base[ar.addresses]=0.;vertical=int(m.joint('pelvis_ty').qposadr[0]);rootrot=np.array([int(m.joint(n).qposadr[0]) for n in ('pelvis_tilt','pelvis_list','pelvis_rotation')]);footids=[m.body('calcn_'+s).id for s in ('r','l')];legnames=[n+'_'+s for s in ('r','l') for n in ('hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle')];legaddr=np.array([int(m.joint(n).qposadr[0]) for n in legnames]);jp=np.zeros((3,m.nv));jr=np.zeros_like(jp)
 common=np.array([np.max([spec['joints'][n+'_'+s]['native_bounds_rad'][0] for s in ('r','l')]) for n in ('hip_flexion','knee_angle','ankle_angle')]);commonhi=np.array([np.min([spec['joints'][n+'_'+s]['native_bounds_rad'][1] for s in ('r','l')]) for n in ('hip_flexion','knee_angle','ankle_angle')]);commonhi[1]=min(commonhi[1],np.deg2rad(-2.))
 def put(q):d.qpos[:]=q;d.qvel[:]=0.;d.xfrc_applied[:]=0.;d.xfrc_applied[m.body('pelvis').id,2]=support*bw;ar.project(d);mujoco.mj_fwdPosition(m,d);mujoco.mj_fwdVelocity(m,d)
 put(base);fields=('xpos','xmat','geom_xpos','geom_xmat','actuator_length','actuator_moment','qM','qfrc_bias','qfrc_passive');fast={name:getattr(d,name).copy() for name in fields};mujoco.mj_forward(m,d);fk_errors={name:float(np.max(abs(getattr(d,name)-fast[name]))) for name in fields}
 if max(fk_errors.values())>1e-10:raise RuntimeError('Fast FK pipeline does not reproduce required native quantities')
 timing={}
 for name,call in [('full_forward',lambda:mujoco.mj_forward(m,d)),('required_position_velocity',lambda:(mujoco.mj_fwdPosition(m,d),mujoco.mj_fwdVelocity(m,d)))]:
  start=perf_counter()
  for _ in range(200):call()
  timing[name]=(perf_counter()-start)/200
 (out/'FK_PIPELINE.json').write_text(json.dumps(dict(scope='Identical held pose,200calls per pipeline; not complete simulation or SNS speedup',max_errors=fk_errors,seconds_per_call=timing,speedup=timing['full_forward']/timing['required_position_velocity']),indent=2))
 def initial_pose(x):
  q=base.copy();q[vertical]=x[3]
  for s in ('r','l'):
   for j,v in zip(('hip_flexion','knee_angle','ankle_angle'),x[:3]):q[m.joint(j+'_'+s).qposadr[0]]=v
  put(q);regions,_=_support_regions(m,d);return np.r_[np.array([r['gap_m']+.001 for r in regions])/.0002,.02*(x[:3]-np.deg2rad([-7.,-15.,8.]))/.1]
 sol=least_squares(initial_pose,np.r_[np.deg2rad([-7.,-15.,8.]),.94],bounds=(np.r_[common,.80],np.r_[commonhi,1.05]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=100)
 initial_pose(sol.x);base=d.qpos.copy();stance=sol.x[:3].copy();low=[];high=[]
 for name in legnames:
  if name in spec['joints']:lo,hi=spec['joints'][name]['native_bounds_rad'];hi=min(hi,np.deg2rad(-2.)) if name.startswith('knee') else hi
  else:lo,hi=np.deg2rad([-10.,10.])
  low.append(lo);high.append(hi)
 low=np.r_[low,np.deg2rad([-10.]*3)];high=np.r_[high,np.deg2rad([10.]*3)];address=np.r_[legaddr,rootrot];qout=[];targets=[];phases=[];weights=[];errors=[];t=np.linspace(0.,period,frames);last=None
 if cartesian and toe_fraction!=.5:
  put(base);fixedpos=d.xpos[footids].copy();fixedrot=d.xmat[footids].reshape(2,3,3).copy();regions,_=_support_regions(m,d);local=[d.xmat[r['body']].reshape(3,3).T@(r['point']-d.xpos[r['body']]) for r in regions];calibration_address=np.r_[address,vertical];oldheight=float(base[vertical]);kneeaddr=[int(m.joint('knee_angle_'+s).qposadr[0]) for s in ('r','l')]
  def calibration(x):
   q=base.copy();q[calibration_address]=x;put(q);contact=np.zeros(m.nv)
   for k,r in enumerate(regions):
    point=d.xpos[r['body']]+d.xmat[r['body']].reshape(3,3)@local[k];mujoco.mj_jac(m,d,jp,jr,point,r['body']);contact+=jp[2]*(1.-support)*bw*.5*(toe_fraction if k%2 else 1.-toe_fraction)
   return np.r_[(d.xpos[footids]-fixedpos).ravel()/.0002,np.concatenate([Rotation.from_matrix(fixedrot[k].T@d.xmat[b].reshape(3,3)).as_rotvec() for k,b in enumerate(footids)])/.001,(d.qfrc_bias-d.qfrc_passive-_external(m,d)-contact)[rootrot]/.2,(d.qpos[kneeaddr]+np.deg2rad(15.))/.01]
  fit=least_squares(calibration,base[calibration_address],bounds=(np.r_[low,oldheight-.03],np.r_[high,oldheight+.005]),max_nfev=120,xtol=1e-11,ftol=1e-11,gtol=1e-11);r=calibration(fit.x);base=d.qpos.copy();stance=np.mean([[base[m.joint(j+'_'+s).qposadr[0]] for j in ('hip_flexion','knee_angle','ankle_angle')] for s in ('r','l')],axis=0);(out/'STANDING_GEOMETRY_CALIBRATION.json').write_text(json.dumps(dict(original_height_m=oldheight,calibrated_height_m=float(base[vertical]),root_angles_native_deg=np.rad2deg(base[rootrot]).tolist(),knees_native_deg=np.rad2deg(base[kneeaddr]).tolist(),position_error_m=float(np.max(abs(r[:6]))*.0002),rotation_error_rad=float(np.max(abs(r[6:12]))*.001),guide_moment_error_Nm=float(np.max(abs(r[12:15]))*.2),toe_fraction=toe_fraction,scope='One constant-height standing calibration,not a per-frame free-height branch or altered spring reference'),indent=2))
 checkpoint=out/'knot_checkpoint.npz'
 if checkpoint.exists():
  saved=np.load(checkpoint);qout=list(saved['qpos']);targets=list(saved['targets']);phases=list(saved['phase']);weights=list(saved['weights']);errors=json.loads(str(saved['errors']));last=qout[-1][address].copy() if qout else None
 for time in t[len(qout):]:
  phase=(time/period+np.array([0.,.5]))%1.;q=base.copy();weight=np.array([smooth(p/.1)*smooth((.65-p)/.1) if p<=.65 else 0. for p in phase]);weight/=weight.sum()
  for s,p in zip(('r','l'),phase):
   b=bump((p-.65)/.35) if p>.65 else 0.
   for j,v in zip(('hip_flexion','knee_angle','ankle_angle'),stance+(np.deg2rad([-18.,-55.,10.])-stance)*b):q[m.joint(j+'_'+s).qposadr[0]]=v
  if lift_clearance>0. and not cartesian:
   for side,p in zip(('r','l'),phase):
    if p<=.65:continue
    u=(p-.65)/.35;goal=-.001+lift_clearance*smooth(min(u/.2,(1.-u)/.2));addresses=np.array([int(m.joint(j+'_'+side).qposadr[0]) for j in ('hip_flexion','knee_angle','ankle_angle')]);original=q[addresses].copy();put(q);original_x=float(d.xpos[m.body('calcn_'+side).id,0]);seed=q.copy()
    def swing_residual(x):
     qq=seed.copy();qq[addresses]=x;put(qq);rr,_=_support_regions(m,d);ss=0 if side=='r' else 2;gap=min(r['gap_m'] for r in rr[ss:ss+2]);return np.r_[(gap-goal)/.0002,(d.xpos[m.body('calcn_'+side).id,0]-original_x)/.05,.005*(x-original)/.1]
    fit_swing=least_squares(swing_residual,np.clip(original,common+1e-8,commonhi-1e-8),bounds=(common,commonhi),max_nfev=100,xtol=1e-10,ftol=1e-10,gtol=1e-10);q[addresses]=fit_swing.x
  if cartesian:
   put(base);targetpos=d.xpos[footids].copy();targetrot=d.xmat[footids].reshape(2,3,3).copy();regions,_=_support_regions(m,d);desired=base.copy()
   for k,p in enumerate(phase):
    if p>.65:
     u=(p-.65)/.35;targetpos[k,0]+=.04*bump(u);targetpos[k,2]+=manifest['cartesian_target']['vertical_lift_m']*smooth(min(u/.2,(1.-u)/.2))
  else:put(q);targetpos=d.xpos[footids].copy();targetrot=d.xmat[footids].reshape(2,3,3).copy();regions,_=_support_regions(m,d);desired=q.copy()
  localpoints=[d.xmat[r['body']].reshape(3,3).T@(r['point']-d.xpos[r['body']]) for r in regions]
  desired[vertical]+=manifest['root_height_reference_peak_m']*sum(smooth(min(((p-.65)/.35)/.2,(1.-((p-.65)/.35))/.2)) for p in phase if p>.65)
  def residual(x):
   q=desired.copy();q[address]=x;put(q);position=np.concatenate([(d.xpos[b]-targetpos[k])/(.0002 if cartesian or phase[k]<=.65 else .002) for k,b in enumerate(footids)]);rotation=np.concatenate([Rotation.from_matrix(targetrot[k].T@d.xmat[b].reshape(3,3)).as_rotvec()/(.001 if cartesian or phase[k]<=.65 else .02) for k,b in enumerate(footids)]);contact=np.zeros(m.nv)
   for k,r in enumerate(regions):
    point=d.xpos[r['body']]+d.xmat[r['body']].reshape(3,3)@localpoints[k];mujoco.mj_jac(m,d,jp,jr,point,r['body']);contact+=jp[2]*((1.-support)*bw*weight[k//2]*(toe_fraction if k%2 else 1.-toe_fraction))
   moment=(d.qfrc_bias-d.qfrc_passive-_external(m,d)-contact)[rootrot]/.2
   clearance=[]
   if lift_clearance>0. and not cartesian:
    actual_regions,_=_support_regions(m,d)
    for k,p in enumerate(phase):
     if p<=.65:clearance.append(0.);continue
     minimum=.0102 if .72<p<.93 else -.0015;gap=min(r['gap_m'] for r in actual_regions[2*k:2*k+2]);clearance.append(max(0.,minimum-gap)/.0002)
   return np.r_[position,rotation,moment,clearance,.005*(x-desired[address])/.1]
  guess=desired[address].copy() if last is None else last.copy();guess[:12]=desired[legaddr] if last is None else guess[:12];guess=np.clip(guess,low+1e-9,high-1e-9);fit=least_squares(residual,guess,bounds=(low,high),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=70);residual(fit.x);last=fit.x.copy();qout.append(d.qpos.copy());targets.append(np.concatenate([targetpos.ravel(),targetrot.ravel()]));phases.append(phase);weights.append(weight);errors.append(dict(position_m=float(np.max(abs(d.xpos[footids]-targetpos))),stance_position_m=float(max(np.max(abs(d.xpos[b]-targetpos[k])) for k,b in enumerate(footids) if phase[k]<=.65)),rotation_rad=float(max(np.max(abs(Rotation.from_matrix(targetrot[k].T@d.xmat[b].reshape(3,3)).as_rotvec())) for k,b in enumerate(footids) if phase[k]<=.65)),root_static_moment_residual_Nm=float(np.max(abs(residual(fit.x)[12:15]))*.2),nfev=int(fit.nfev)))
  pending=out/'knot_checkpoint.pending'
  with pending.open('wb') as handle:np.savez_compressed(handle,qpos=qout,targets=targets,phase=phases,weights=weights,errors=json.dumps(errors))
  pending.replace(checkpoint)
 qout=np.array(qout);endpoint=float(np.max(abs(qout[-1]-qout[0])));qout[-1]=qout[0];spline=PchipInterpolator(np.r_[t[:-1]-period,t,t[1:]+period],np.concatenate([qout[:-1],qout,qout[1:]]),axis=0);np.savez_compressed(out/'knots.npz',t=t,qpos=qout);t=np.linspace(0.,period,round(period/.002)+1);qout=spline(t);qvel=spline(t,1);qacc=spline(t,2);frames=len(t);phases=np.array([(value/period+np.array([0.,.5]))%1. for value in t]);weights=np.array([[smooth(p/.1)*smooth((.65-p)/.1) if p<=.65 else 0. for p in phase] for phase in phases]);weights/=weights.sum(axis=1)[:,None];mesh=[];speeds=[];strict=[];central=[]
 for k,time in enumerate(t):
  put(qout[k]);qout[k]=d.qpos;d.qvel[:]=qvel[k]
  for dep,ind,c,offset,qa in zip(ar.dependent,ar.independent,ar.coefficients,ar.offsets,ar.qadr):
   x=d.qpos[qa]-offset;first=c[1]+x*(2*c[2]+x*(3*c[3]+4*x*c[4]));second=2*c[2]+x*(6*c[3]+12*x*c[4]);qvel[k,dep]=first*qvel[k,ind];qacc[k,dep]=first*qacc[k,ind]+second*qvel[k,ind]**2
  d.qvel[:]=qvel[k];mujoco.mj_forward(m,d);regions,_=_support_regions(m,d);mesh.append([r['gap_m'] for r in regions]);ss=[]
  for r in regions:mujoco.mj_jac(m,d,jp,jr,r['point'],r['body']);ss.append(float(np.linalg.norm(jp@d.qvel)))
  speeds.append(ss);strict.append(max(float(np.rad2deg(d.qpos[m.joint('knee_angle_'+s).qposadr[0]])) for s in ('r','l')))
 mesh=np.array(mesh);speeds=np.array(speeds);phases=np.array(phases);stance_mask=phases<=.65;central_mask=(phases>.72)&(phases<.93);footgap=mesh.reshape(frames,2,2).min(axis=2);foot_speed=speeds.reshape(frames,2,2).max(axis=2);metrics=dict(max_foot_pose_error_m=max(e['position_m'] for e in errors),max_foot_rotation_error_rad=max(e['rotation_rad'] for e in errors),max_root_static_moment_residual_Nm=max(e['root_static_moment_residual_Nm'] for e in errors),max_stance_point_speed_mps=float(foot_speed[stance_mask].max()),minimum_all_mesh_gap_m=float(mesh.min()),minimum_central_swing_gap_m=float(footgap[central_mask].min()),maximum_knee_extension_deg=max(strict),periodic_endpoint_correction_max_rad_or_m=endpoint,stance_native_hip_knee_ankle_deg=np.rad2deg(stance).tolist())
 metrics['max_stance_body_position_error_m']=max(e['stance_position_m'] for e in errors);metrics['dense_source_envelope_max_excess_rad']=max(max(float(np.max(row['native_bounds_rad'][0]-qout[:,row['qpos_address']])),float(np.max(qout[:,row['qpos_address']]-row['native_bounds_rad'][1])),0.) for row in spec['joints'].values());metrics['root_height_range_m']=[float(qout[:,vertical].min()),float(qout[:,vertical].max())]
 accepted=metrics['max_stance_body_position_error_m']<=.0002 and metrics['max_foot_rotation_error_rad']<=.001 and metrics['max_stance_point_speed_mps']<=.005 and metrics['minimum_all_mesh_gap_m']>=-.002 and metrics['minimum_central_swing_gap_m']>=.01 and metrics['maximum_knee_extension_deg']<=-2. and metrics['dense_source_envelope_max_excess_rad']<=1e-10 and endpoint<=1e-5
 np.savez_compressed(out/'teacher.npz',t=t,qpos=qout,qvel=qvel,qacc=qacc,phase=phases,nominal_foot_load_share=weights,mesh_gaps=mesh,point_speeds=speeds,foot_targets=targets);report=dict(accepted_kinematics=bool(accepted),metrics=metrics,per_frame_solve=errors,scope='Offline kinematic reference only. Full root/contact dynamic and source-muscle recruitment feasibility remain required before neural fitting.',manifest=manifest);(out/'KINEMATIC_RESULT.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(accepted_kinematics=bool(accepted),metrics=metrics),indent=2))

def calibrate_native_stance(model_folder,spec_path,source,out,phase_origin=.075,knee_target=15.,knee_side=None,support_override=None,unload_from=None,unload_side='l',unload_gap=.003):
 """Actual inverse-contact static calibration, with original native guide model."""
 out.mkdir(parents=True,exist_ok=True);cfg=json.loads((source/'MANIFEST.json').read_text());spec=json.loads(spec_path.read_text());z=np.load(source/'teacher.npz');seed=z['qpos'][int(np.argmin(abs(z['t']-phase_origin*cfg['period_s'])))].copy();m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);ar=CoupledArms(m);bw=float(m.body_mass.sum()*9.81);footids=[m.body('calcn_'+s).id for s in ('r','l')];legnames=[n+'_'+s for s in ('r','l') for n in ('hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle')];legaddr=np.array([int(m.joint(n).qposadr[0]) for n in legnames]);address=np.r_[legaddr,np.arange(6)];kneeaddr=[int(m.joint('knee_angle_'+s).qposadr[0]) for s in ('r','l')];force6=np.zeros(6);lo=[];hi=[]
 if unload_from is not None:seed=np.load(unload_from/'source_recruitment.npz')['qpos'].copy()
 for name in legnames:
  if name in spec['joints']:lower,upper=spec['joints'][name]['native_bounds_rad'];upper=min(upper,np.deg2rad(-2.)) if name.startswith('knee') else upper
  else:lower,upper=np.deg2rad([-10.,10.])
  lo.append(lower);hi.append(upper)
 lo=np.r_[lo,seed[0]-.03,seed[1]-.02,seed[2]-.003,np.deg2rad([-10.]*3),-.002,np.deg2rad(-.5)];hi=np.r_[hi,seed[0]+.03,seed[1]+.02,seed[2]+.003,np.deg2rad([10.]*3),.002,np.deg2rad(.5)]
 if support_override is not None:cfg['support_fraction']=support_override
 def update(q):d.qpos[:]=q;ar.project(d);d.qvel[:]=0.;d.qacc[:]=0.;d.xfrc_applied[:]=0.;d.xfrc_applied[m.body('pelvis').id,2]=cfg['support_fraction']*bw;mujoco.mj_inverse(m,d)
 update(seed);positions=d.xpos[footids].copy();orientations=d.xmat[footids].reshape(2,3,3).copy();last={};manifest=dict(format='aarl_native_stance_calibration_v1',source=str(source),source_code_sha256=sha(__file__),support_fraction=cfg['support_fraction'],distal_toe_fraction_target=.08,knee_flexion_target_deg=knee_target,knee_target_side=knee_side,source_phase=phase_origin,scope='Original native model/soft contact,static qvel/qacc zero. Foot vertical/pitch microdeformation is solved,not a geometry change. This is an engineering calibration,not human pressure-distribution data.',acceptance=dict(root_force_max_BW=.01,root_moment_max_BWm=.01,foot_pose_position_max_m=.0002,foot_pose_rotation_max_rad=.001,penetration_max_m=.002,toe_fraction_error_max=.02,knee_extension_max_deg=-2.))
 if unload_from is not None:
  side=0 if unload_side=='r' else 1;regions,_=_support_regions(m,d);gap=min(r['gap_m'] for r in regions[2*side:2*side+2]);positions[side,2]+=unload_gap-gap;manifest['unload_before_advancing']=dict(source_recruitment=str(unload_from),unloaded_side=unload_side,nominal_gap_m=unload_gap,constant_fore_aft_position=True,constant_foot_orientation=True,scope='Unloaded near-flat intermediate keyframe; common native contact microdeformation may change the final gap slightly.')
 (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
 def residual(x):
  q=seed.copy();q[address]=x[:18];update(q);target=positions.copy();target[:,2]+=x[18];R=Rotation.from_rotvec([0.,x[19],0.]).as_matrix();rotation=np.concatenate([Rotation.from_matrix((R@orientations[k])@d.xmat[b].reshape(3,3).T).as_rotvec() for k,b in enumerate(footids)]);root=d.qfrc_inverse[:6]-_external(m,d)[:6];normal=toe=0.;contact=[]
  for i in range(d.ncon):
   c=d.contact[i];a,b=int(c.geom1),int(c.geom2)
   if m.geom_type[a]==mujoco.mjtGeom.mjGEOM_PLANE:g=b
   elif m.geom_type[b]==mujoco.mjtGeom.mjGEOM_PLANE:g=a
   else:continue
   name=m.body(int(m.geom_bodyid[g])).name
   if not name.startswith(('calcn_','toes_')):continue
   mujoco.mj_contactForce(m,d,i,force6);normal+=force6[0];toe+=force6[0] if name.startswith('toes_') else 0.;contact.append(dict(body=name,normal_N=float(force6[0]),gap_m=float(c.dist)))
  ratio=toe/max(normal,1e-12);regions,_=_support_regions(m,d);gaps=np.array([r['gap_m'] for r in regions]);last.update(root=root.copy(),contacts=contact,gaps=gaps,ratio=ratio,normal=normal,target=target,rotation=rotation);knee=d.qpos[kneeaddr].mean() if knee_side is None else d.qpos[kneeaddr[0 if knee_side=='r' else 1]];return np.r_[(d.xpos[footids]-target).ravel()/.0002,rotation/.001,root[:3]/(bw*.005),root[3:]/(bw*.005),(knee+np.deg2rad(knee_target))/.01,(ratio-.08)/.01,np.maximum(0.,-.002-gaps)/.00005]
 x=np.r_[seed[address],0.,0.];fit=least_squares(residual,x,bounds=(lo,hi),max_nfev=150,xtol=1e-10,ftol=1e-10,gtol=1e-10);r=residual(fit.x);metrics=dict(root_force_max_BW=float(np.max(abs(last['root'][:3]))/bw),root_moment_max_BWm=float(np.max(abs(last['root'][3:]))/bw),foot_position_error_m=float(np.max(abs(d.xpos[footids]-last['target']))),foot_rotation_error_rad=float(np.max(abs(last['rotation']))),minimum_mesh_gap_m=float(last['gaps'].min()),actual_distal_toe_load_fraction=float(last['ratio']),native_normal_total_BW=float(last['normal']/bw),knees_native_deg=np.rad2deg(d.qpos[kneeaddr]).tolist(),root_qpos_native=d.qpos[:6].tolist(),foot_common_vertical_shift_m=float(fit.x[18]),foot_pitch_adjustment_deg=float(np.rad2deg(fit.x[19])),nfev=int(fit.nfev));accepted=metrics['root_force_max_BW']<=.01 and metrics['root_moment_max_BWm']<=.01 and metrics['foot_position_error_m']<=.0002 and metrics['foot_rotation_error_rad']<=.001 and metrics['minimum_mesh_gap_m']>=-.002 and abs(metrics['actual_distal_toe_load_fraction']-.08)<=.02 and max(metrics['knees_native_deg'])<=-2.;np.savez_compressed(out/'native_stance.npz',qpos=d.qpos.copy(),required_generalized_muscle_force=d.qfrc_inverse-_external(m,d));(out/'NATIVE_STANCE_RESULT.json').write_text(json.dumps(dict(accepted_native_static_pose=bool(accepted),metrics=metrics,contacts=last['contacts'],manifest=manifest),indent=2));print(json.dumps(dict(accepted_native_static_pose=bool(accepted),metrics=metrics),indent=2))

def refine(model_folder,spec_path,source,out,cartesian_height=None):
 """Dense constrained IK projection; never used inside a forward SNS simulation."""
 out.mkdir(parents=True,exist_ok=True);parent=json.loads((source/'MANIFEST.json').read_text());z=np.load(source/'teacher.npz');q=z['qpos'].copy();t=z['t'];phase=z['phase'];dt=float(np.median(np.diff(t)));m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);ar=CoupledArms(m);spec=json.loads(spec_path.read_text());jp=np.zeros((3,m.nv));jr=np.zeros_like(jp);position=[];rotation=[];bounds_failures=0;projection_failures=[]
 if cartesian_height is not None:
  if 'cartesian_target' not in parent:raise ValueError('Height override requires explicit Cartesian target')
  parent['cartesian_target']['vertical_lift_m']=cartesian_height
 manifest=dict(parent,refinement_source=str(source),refinement_source_sha256=sha(source/'teacher.npz'),source_code_sha256=sha(__file__),refinement='Exact dense stance-body inverse kinematics; C2 world-body pose anchoring over first/last15percent swing; separate hindfoot/forefoot inequalities avoid nonsmooth minimum ridge; source-anchored C2 floor reaches11mm by15percent swing; unchanged acceptance thresholds; periodic centered finite-difference independent rates plus exact equality-follower derivatives')
 (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
 def put(x):d.qpos[:]=x;d.qvel[:]=0.;ar.project(d);mujoco.mj_fwdPosition(m,d);mujoco.mj_fwdVelocity(m,d)
 put(q[0]);footids=[m.body('calcn_'+s).id for s in ('r','l')];targetpos=d.xpos[footids].copy();targetrot=d.xmat[footids].reshape(2,3,3).copy();initial_regions,_=_support_regions(m,d);basegap=[min(r['gap_m'] for r in initial_regions[2*i:2*i+2]) for i in range(2)]
 for k in range(len(t)):
  for side_index,s in enumerate(('r','l')):
   leg=[n+'_'+s for n in ('hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle')];addresses=np.array([int(m.joint(n).qposadr[0]) for n in leg]);columns=[JOINTS.index(n) for n in leg];lo=np.array([spec['joints'][n]['native_bounds_rad'][0] if n in spec['joints'] else np.deg2rad(-10.) for n in leg]);hi=np.array([min(spec['joints'][n]['native_bounds_rad'][1],np.deg2rad(-2.)) if n.startswith('knee') else spec['joints'][n]['native_bounds_rad'][1] if n in spec['joints'] else np.deg2rad(10.) for n in leg]);body=footids[side_index]
   if phase[k,side_index]<=.65 or 'cartesian_target' in parent:
    desired_position=targetpos[side_index].copy();desired_rotation=targetrot[side_index]
    if phase[k,side_index]>.65:
     u=(phase[k,side_index]-.65)/.35;desired_position[0]+=parent['cartesian_target']['forward_excursion_m']*bump(u);desired_position[2]+=parent['cartesian_target']['vertical_lift_m']*smooth(min(u/.2,(1.-u)/.2))
    for _ in range(5):
     put(q[k]);error=np.r_[desired_position-d.xpos[body],Rotation.from_matrix(desired_rotation@d.xmat[body].reshape(3,3).T).as_rotvec()]
     if np.max(abs(error))<1e-11:break
     ar.update(d);mujoco.mj_jacBody(m,d,jp,jr,body);J=np.vstack([jp,jr])@ar.E[:,columns];step=np.linalg.solve(J,error);candidate=q[k,addresses]+step
     if np.any(candidate<lo-1e-10) or np.any(candidate>hi+1e-10):
      seed=q[k].copy()
      def pose_objective(x):
       qq=seed.copy();qq[addresses]=x;put(qq);return np.r_[(d.xpos[body]-desired_position)/.0002,Rotation.from_matrix(desired_rotation@d.xmat[body].reshape(3,3).T).as_rotvec()/.001]
      fit=least_squares(pose_objective,q[k,addresses],bounds=(lo,hi),max_nfev=100,xtol=1e-11,ftol=1e-11,gtol=1e-11);q[k,addresses]=fit.x;rr=pose_objective(fit.x)
      if np.max(abs(rr[:3]))*.0002>parent['acceptance']['stance_body_position_error_max_m'] or np.max(abs(rr[3:]))*.001>parent['acceptance']['stance_rotation_error_max_rad']:bounds_failures+=1;projection_failures.append(dict(index=k,side=s,phase=float(phase[k,side_index]),position_error_m=float(np.max(abs(rr[:3]))*.0002),rotation_error_rad=float(np.max(abs(rr[3:]))*.001),angles_deg=np.rad2deg(fit.x).tolist()))
      break
     q[k,addresses]=candidate
   else:
    u=(phase[k,side_index]-.65)/.35;put(q[k]);blend=smooth(min(u/.15,(1.-u)/.15));desired_position=targetpos[side_index]+blend*(d.xpos[body]-targetpos[side_index]);relative=Rotation.from_matrix(targetrot[side_index].T@d.xmat[body].reshape(3,3)).as_rotvec();desired_rotation=targetrot[side_index]@Rotation.from_rotvec(blend*relative).as_matrix()
    if blend<1.:
     for _ in range(5):
      put(q[k]);error=np.r_[desired_position-d.xpos[body],Rotation.from_matrix(desired_rotation@d.xmat[body].reshape(3,3).T).as_rotvec()]
      if np.max(abs(error))<1e-11:break
      ar.update(d);mujoco.mj_jacBody(m,d,jp,jr,body);J=np.vstack([jp,jr])@ar.E[:,columns];candidate=q[k,addresses]+np.linalg.solve(J,error)
      if np.any(candidate<lo-1e-10) or np.any(candidate>hi+1e-10):break
      q[k,addresses]=candidate
    put(q[k]);regions,_=_support_regions(m,d);ss=side_index*2;gap=min(r['gap_m'] for r in regions[ss:ss+2]);goal=basegap[side_index]+(.011-basegap[side_index])*smooth(min(u/.15,(1.-u)/.15))
    if gap<goal:
     indices=np.array([0,3,4]);aa=addresses[indices];seed=q[k].copy();original=seed[aa].copy()
     def objective(x):
      qq=seed.copy();qq[aa]=x;put(qq);rr,_=_support_regions(m,d);g=np.array([r['gap_m'] for r in rr[ss:ss+2]]);return np.r_[np.maximum(0.,goal-g)/.0001,.02*(x-original)/.1]
     fit=least_squares(objective,original,bounds=(lo[indices],hi[indices]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=60);q[k,aa]=fit.x
  put(q[k]);q[k]=d.qpos.copy();position.append(max(float(np.max(abs(d.xpos[b]-targetpos[i]))) for i,b in enumerate(footids) if phase[k,i]<=.65));rotation.append(max(float(np.max(abs(Rotation.from_matrix(targetrot[i]@d.xmat[b].reshape(3,3).T).as_rotvec()))) for i,b in enumerate(footids) if phase[k,i]<=.65))
 endpoint=float(np.max(abs(q[-1]-q[0])));q[-1]=q[0];cycle=q[:-1];vcycle=(np.roll(cycle,-1,axis=0)-np.roll(cycle,1,axis=0))/(2*dt);acycle=(np.roll(cycle,-1,axis=0)-2*cycle+np.roll(cycle,1,axis=0))/dt**2;v=np.vstack([vcycle,vcycle[0]]);a=np.vstack([acycle,acycle[0]]);mesh=[];speed=[]
 for k in range(len(t)):
  put(q[k])
  for dep,ind,c,offset,qa in zip(ar.dependent,ar.independent,ar.coefficients,ar.offsets,ar.qadr):
   x=d.qpos[qa]-offset;first=c[1]+x*(2*c[2]+x*(3*c[3]+4*x*c[4]));second=2*c[2]+x*(6*c[3]+12*x*c[4]);v[k,dep]=first*v[k,ind];a[k,dep]=first*a[k,ind]+second*v[k,ind]**2
  d.qvel[:]=v[k];mujoco.mj_fwdVelocity(m,d);regions,_=_support_regions(m,d);mesh.append([r['gap_m'] for r in regions]);ss=[]
  for r in regions:mujoco.mj_jac(m,d,jp,jr,r['point'],r['body']);ss.append(float(np.linalg.norm(jp@d.qvel)))
  speed.append(ss)
 mesh=np.array(mesh);speed=np.array(speed);stance=phase<=.65;central=(phase>.72)&(phase<.93);footgap=mesh.reshape(-1,2,2).min(axis=2);footspeed=speed.reshape(-1,2,2).max(axis=2);excess=max(max(float(np.max(row['native_bounds_rad'][0]-q[:,row['qpos_address']])),float(np.max(q[:,row['qpos_address']]-row['native_bounds_rad'][1])),0.) for row in spec['joints'].values());knees=[int(m.joint('knee_angle_'+s).qposadr[0]) for s in ('r','l')];metrics=dict(max_stance_body_position_error_m=max(position),max_stance_rotation_error_rad=max(rotation),max_stance_point_speed_mps=float(footspeed[stance].max()),minimum_all_mesh_gap_m=float(mesh.min()),minimum_central_swing_gap_m=float(footgap[central].min()),maximum_knee_extension_deg=float(np.rad2deg(q[:,knees].max())),dense_source_envelope_max_excess_rad=excess,periodic_endpoint_correction_max_rad_or_m=endpoint,stance_projection_bound_failures=bounds_failures)
 accepted=max(position)<=.0002 and max(rotation)<=.001 and metrics['max_stance_point_speed_mps']<=.005 and mesh.min()>=-.002 and footgap[central].min()>=.01 and metrics['maximum_knee_extension_deg']<=-2.+1e-8 and excess<=1e-10 and endpoint<=1e-5 and bounds_failures==0
 np.savez_compressed(out/'teacher.npz',t=t,qpos=q,qvel=v,qacc=a,phase=phase,mesh_gaps=mesh,point_speeds=speed);(out/'KINEMATIC_RESULT.json').write_text(json.dumps(dict(accepted_kinematics=bool(accepted),metrics=metrics,projection_failures=projection_failures,manifest=manifest,scope='Dense constrained offline reference,not forward gait validation. Full root/contact source-dynamic inverse remains required before neural fitting.'),indent=2));print(json.dumps(dict(accepted_kinematics=bool(accepted),metrics=metrics),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',required=True,type=Path);p.add_argument('--spec',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--period',type=float,default=2.4);p.add_argument('--frames',type=int,default=61);p.add_argument('--support',type=float,default=.5);p.add_argument('--native-support',type=float);p.add_argument('--lift-clearance',type=float,default=0.);p.add_argument('--refine',type=Path);p.add_argument('--native-stance',type=Path);p.add_argument('--native-phase',type=float,default=.075);p.add_argument('--knee-target',type=float,default=15.);p.add_argument('--knee-target-side',choices=('r','l'));p.add_argument('--unload-from',type=Path);p.add_argument('--unload-side',choices=('r','l'),default='l');p.add_argument('--unload-gap',type=float,default=.003);p.add_argument('--cartesian',action='store_true');p.add_argument('--cartesian-height',type=float);p.add_argument('--toe-fraction',type=float,default=.5);a=p.parse_args()
 if a.native_stance:calibrate_native_stance(a.model,a.spec,a.native_stance,a.out,a.native_phase,a.knee_target,a.knee_target_side,a.native_support,a.unload_from,a.unload_side,a.unload_gap)
 elif a.refine:refine(a.model,a.spec,a.refine,a.out,a.cartesian_height)
 else:prepare(a.model,a.spec,a.out,a.period,a.frames,a.support,a.lift_clearance,a.cartesian,a.cartesian_height,a.toe_fraction)
