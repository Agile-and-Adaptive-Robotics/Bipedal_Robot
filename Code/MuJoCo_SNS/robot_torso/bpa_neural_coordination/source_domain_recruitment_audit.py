"""Offline all-route feasibility and signed equality-aware recruitment loads."""
import json,hashlib,csv
from pathlib import Path
import numpy as np
import gait_standing_family as F
from mechanics import CoupledArms,JOINTS
from valve_authority_diagnostics import load
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');DEST=ROOT/'source_domain_recruitment';SOURCE=ROOT/'valve_authority'

def model(reference):
 cp=json.loads(Path(reference['checkpoint']).read_text());q0=np.load(Path(cp['configuration']).with_suffix('.npz'))['qfull'][0].copy();m0=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));d0=F.mujoco.MjData(m0);F.mujoco.mj_resetDataKeyframe(m0,d0,0);d0.qpos[:]=q0;F.mujoco.mj_forward(m0,d0);m=F.R.apply_harness(m0,d0,kxy=1500.)
 force=json.loads(Path(reference['force_law']).read_text())
 for row in force['muscles'].values():
  for name in row['actuator_names']:
   i=m.actuator(name).id
   for key in ('gainprm','biasprm','dynprm','lengthrange'):getattr(m,'actuator_'+key)[i]=row[key]
 for side in ('r','l'):
  j=m.joint('knee_angle_'+side);m.jnt_limited[j.id]=True;m.jnt_range[j.id]=np.radians([-90.,0.]);m.jnt_margin[j.id]=np.radians(2.);m.jnt_solref[j.id]=[.004,1.];m.jnt_solimp[j.id]=[.999,.9999,.001,.5,2.]
 d=F.mujoco.MjData(m);F.mujoco.mj_forward(m,d)
 for axis in ('tx','ty','tz'):
  j=m.joint('pelvis_'+axis)
  if abs(d.xaxis[j.id,2])>.8:m.jnt_stiffness[j.id]*=.01;m.dof_damping[j.dofadr[0]]*=.42
 # Pure offline force attribution: no native Hill force mixed into BPA demand.
 m.actuator_gaintype[:]=F.mujoco.mjtGain.mjGAIN_FIXED;m.actuator_gainprm[:]=0.;m.actuator_biastype[:]=F.mujoco.mjtBias.mjBIAS_NONE
 return m,d,q0

def main():
 DEST.mkdir(exist_ok=True);cfg=json.loads((SOURCE/'PREREGISTRATION.json').read_text());trials=[json.loads(s) for s in (SOURCE/'pressure/rustuna/TRIALS.jsonl').read_text().splitlines()];ref=json.loads(Path(trials[0]['diagnostics']['conditions'][0]['result']).read_text());m,d,q0=model(ref);arms=CoupledArms(m);names=[m.actuator(i).name for i in range(m.nu)];rows=[];snapshots=[]
 for t in trials:
  for c in t['diagnostics']['conditions']:
   folder=Path(c['result']).parent;final=load(folder,'w2lvar_case.final_actuator_mechanics_state.npz');spec=next(s for s in cfg['training'] if s['name']==c['condition']);actual=json.loads(Path(spec['path']).read_text());d.qpos[:]=final['qpos'];d.qvel[:]=0.;d.qfrc_applied[:]=0.;F.mujoco.mj_forward(m,d);A=arms.update(d).copy();E=arms.E.copy();saved=d.qpos.copy()
   if m.nq!=m.nv:raise ValueError('Tangent finite difference requires scalar-coordinate model')
   fd=np.zeros_like(A)
   for k in range(len(JOINTS)):
    for sign in (-1.,1.):d.qpos[:]=saved+sign*1e-6*E[:,k];F.mujoco.mj_forward(m,d);fd[k]-=sign*d.actuator_length/(2e-6)
   d.qpos[:]=saved;F.mujoco.mj_forward(m,d);error=float(np.max(np.abs(A-fd)))
   if error>1e-6:raise ValueError('Equality-aware pulling torque sign/FD mismatch:'+str(error))
   length_error=None
   if c['failure']:
    observed=load(folder,'w2lvar_case.failure_state.npz');length_error=float(np.max(np.abs(d.actuator_length-observed['actuator_length'])))
    if length_error>1e-10:raise ValueError('Reconstructed observed lengths mismatch')
   length=d.actuator_length.copy();qerror=(saved-q0)[arms.addresses];gravity=E.T@d.qfrc_bias;passive=E.T@d.qfrc_passive;frozen_contact=E.T@d.qfrc_constraint;demand=gravity-passive-frozen_contact;target=np.zeros(m.nu);minimum=np.zeros(m.nu)
   for name,p in actual['actuators'].items():
    i=names.index(name);target[i]=float(np.clip(final['ctrl'][i],0,1))*p['force_command_scale_N'];minimum[i]=max(0.,p['series_stiffness_N_m']*(length[i]-p['resting_length_m']-p['nonmuscle_length_m']))
   deficit=np.maximum(0,minimum-target);tau=A@target;increment=A@deficit;residual=demand-tau;improvement=float(np.linalg.norm(residual)-np.linalg.norm(residual-increment));posture_restoring=-float(increment@qerror)
   for i in np.flatnonzero(deficit>1e-7):
    unit=A[:,i];delta=unit*deficit[i];rows.append(dict(trial=t['number'],cap_kpa=t['params']['pressure_cap_kpa'],condition=c['condition'],route=names[i],length_m=float(length[i]),neural_target_N=float(target[i]),geometric_minimum_N=float(minimum[i]),deficit_N=float(deficit[i]),posture_restoring_work_J=-float(delta@qerror),frozen_contact_torque_residual_improvement_Nm=float(np.linalg.norm(residual)-np.linalg.norm(residual-delta)),signed_moment_arms_m=dict(zip(JOINTS,unit.tolist())),incremental_torques_Nm=dict(zip(JOINTS,delta.tolist()))))
   snapshots.append(dict(trial=t['number'],condition=c['condition'],finite_difference_max_arm_error_m=error,observed_length_max_error_m=length_error,deficient_routes=[names[i] for i in np.flatnonzero(deficit>1e-7)],total_missing_requested_tension_N=float(deficit.sum()),incremental_torque_Nm=dict(zip(JOINTS,increment.tolist())),posture_restoring_work_J=posture_restoring,frozen_contact_torque_residual_improvement_Nm=improvement,reference_demand_Nm=dict(zip(JOINTS,demand.tolist())),original_requested_torque_Nm=dict(zip(JOINTS,tau.tolist())),joint_error_rad=dict(zip(JOINTS,qerror.tolist()))))
 report=dict(machine='easteregg2',joints=JOINTS,snapshots=snapshots,all_deficient_routes=rows,controller_signal_derivation='Tmin=max(0,K*(L-rest-fixed)); e_target=max(0,Tmin+declaredmargin-uMN*Fscale). Optional measured tension/pressure reserve only, no true gain. This is proposed neural sensory error, never enforced force/recruitment.',scope='Offline observed poses, equality-aware pulling moment arms verified against tangent finite difference. Static zero-native-force contact recomputation is only a frozen-load diagnostic, not actual-contact reconstruction or proof of dynamic balance. Posture restoring work uses actual recorded q minus same initial corrected pose. Incremental requested torque is NOT equal to actual force change: nearest-feasible nominal inversion already exceeds original NN requests when they are infeasible. No geometry/model modification, no time integration or fit yet.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (DEST/'OFFLINE_RECRUITMENT_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(snapshots=len(snapshots),deficient_names=sorted({r['route'] for r in rows}),snapshots_summary=[dict(trial=s['trial'],condition=s['condition'],deficit_N=s['total_missing_requested_tension_N'],posture_work_J=s['posture_restoring_work_J'],frozen_residual_improvement_Nm=s['frozen_contact_torque_residual_improvement_Nm']) for s in snapshots]),indent=2))

def measured_load_probe():
 cfg=json.loads((SOURCE/'PREREGISTRATION.json').read_text());trials=[json.loads(s) for s in (SOURCE/'pressure/rustuna/TRIALS.jsonl').read_text().splitlines()];reference=json.loads(Path(trials[0]['diagnostics']['conditions'][0]['result']).read_text());m,d,q0=model(reference);ar=CoupledArms(m);names=[m.actuator(i).name for i in range(m.nu)];d.qpos[:]=q0;F.mujoco.mj_forward(m,d);origin_com=d.subtree_com[0].copy();rows=[]
 for t in trials:
  for c in t['diagnostics']['conditions']:
   folder=Path(c['result']).parent;end=load(folder,'w2lvar_case.final_actuator_mechanics_state.npz');pressure=load(folder,'w2lvar_case.pressure.npz');spec=next(s for s in cfg['training'] if s['name']==c['condition']);actual=json.loads(Path(spec['path']).read_text());d.qpos[:]=end['qpos'];d.qvel[:]=end['qvel'];d.qfrc_applied[:]=0.;F.mujoco.mj_forward(m,d);A=ar.update(d).copy();T=np.zeros(m.nu);missing=[];delta=np.zeros(m.nu);error=(d.qpos-q0)[ar.addresses];comerror=d.subtree_com[0]-origin_com
   for i,name in enumerate(actual['actuators']):
    aid=names.index(name);value=end['force'][i]
    if not np.isfinite(value):value=pressure['samples'][-1,2,i];missing.append(name)
    T[aid]=value;p=actual['actuators'][name];minimum=max(0.,p['series_stiffness_N_m']*(d.actuator_length[aid]-p['resting_length_m']-p['nonmuscle_length_m']));target=float(np.clip(end['ctrl'][aid],0,1))*p['force_command_scale_N']
    if minimum>target+1e-7:delta[aid]=max(0.,minimum+5.-value)
   moment=d.actuator_moment.copy();d.qfrc_applied[:]=-moment.T@T;F.mujoco.mj_forward(m,d);base=d.qacc.copy();base_constraint=d.qfrc_constraint.copy();jac=np.zeros((3,m.nv));F.mujoco.mj_jacSubtreeCom(m,d,jac,0);points=[];normalforces=[]
   plane=m.geom('ground-plane').id
   for k in range(d.ncon):
    contact=d.contact[k]
    if plane not in (contact.geom1,contact.geom2):continue
    other=contact.geom2 if contact.geom1==plane else contact.geom1;body=m.body(int(m.geom_bodyid[other])).name
    if body not in ('calcn_r','calcn_l','toes_r','toes_l','talus_r','talus_l'):continue
    force=np.zeros(6);F.mujoco.mj_contactForce(m,d,k,force)
    if force[0]>1.:points.append(contact.pos[:2].copy());normalforces.append(float(force[0]))
   supportmargin=None
   if len(points)>=3:
    from scipy.spatial import ConvexHull
    hull=ConvexHull(np.asarray(points));supportmargin=float(np.min(-(hull.equations[:,:2]@d.subtree_com[0,:2]+hull.equations[:,2])/np.linalg.norm(hull.equations[:,:2],axis=1)))
   d.qfrc_applied[:]=-moment.T@(T+delta);F.mujoco.mj_forward(m,d);da=d.qacc-base;torque=A@delta;residual=ar.E.T@(d.qfrc_bias-d.qfrc_passive-base_constraint)-A@T
   rows.append(dict(trial=t['number'],condition=c['condition'],temporal_force_fallback_routes=missing,force_source_scope='Final actual bank tension where finite. Failed route uses previous valid measured tension at most one2ms physics step older; this uncertainty is explicit, not filled from minimum.',proposed_physical_headroom_N=5.,added_measured_headroom_force_N={names[i]:float(delta[i]) for i in np.flatnonzero(delta>0)},posture_restoring_work_J=-float(torque@error),incremental_mechanical_power_W=float(torque@end['qvel'][ar.rows]),loaded_contact_COM_margin_m=supportmargin,loaded_contact_points=len(points),incremental_selected_coordinate_acceleration_rad_s2=dict(zip(JOINTS,da[ar.rows].tolist())),incremental_COM_acceleration_world_m_s2=(jac@da).tolist(),COM_error_world_m=comerror.tolist(),COM_horizontal_restoring_metric_m2_s2=-float(comerror[:2]@(jac@da)[:2]),frozen_actual_force_residual_improvement_Nm=float(np.linalg.norm(residual)-np.linalg.norm(residual-torque)),scope='Explicit offline generalized-force probes only, no integration/neuraloverride/forced physicalminimum. Original measured forces plus hypothetical positive5N sleeve headroom deficit. COM metric relative initial standing reference, legacyxanterior/ylateral/zsuperior. Loaded contact hull uses recomputedforces>1N at observedq, not groundtruth archivedcontactwrenches; instantaneous results do not prove closedloop stability or anatomical jointreaction safety.'))
 report=dict(machine='easteregg2',rows=rows,architecture_chosen=False,scope='Torque/load compatibility audit only, new fitting/validation unstarted. Hypothetical headroom is not an enforced force floor.')
 (DEST/'MEASURED_TENSION_LOAD_PROBE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':
 import sys
 measured_load_probe() if '--measured-load' in sys.argv else main()
