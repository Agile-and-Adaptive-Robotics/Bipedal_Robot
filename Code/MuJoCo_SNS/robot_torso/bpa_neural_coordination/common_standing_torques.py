"""Equality-aware shared fixed-pose torque demand and observed recruitment projections."""
import json
from pathlib import Path
import numpy as np
import gait_standing_family as F
from mechanics import CoupledArms,JOINTS
from muscle_map import classify
from physics_fidelity.standing_solver import solve_standing_activations
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis'

def main():
 baseline=json.loads((OUT/'s3k_s3k_v1.result.json').read_text());cp=json.loads(Path(baseline['checkpoint']).read_text());q0=np.load(Path(cp['configuration']).with_suffix('.npz'))['qfull'][0].copy();m0=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));d0=F.mujoco.MjData(m0);F.mujoco.mj_resetDataKeyframe(m0,d0,0);d0.qpos[:]=q0;F.mujoco.mj_forward(m0,d0);m=F.R.apply_harness(m0,d0,kxy=1500.)
 force=json.loads(Path(baseline['force_law']).read_text())
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
 d.qpos[:]=q0;F.mujoco.mj_forward(m,d);activation=solve_standing_activations(m,d,preload=False,bilateral=True);match=np.array_equal(activation,baseline['standing_activation_targets'])
 if not np.allclose(activation,baseline['standing_activation_targets'],atol=1e-12,rtol=0):raise ValueError('Common mechanical reconstruction does not match actual standing targets')
 arms=CoupledArms(m);arms.update(d);E=arms.E.copy();gravity=E.T@d.qfrc_bias;support=E.T@(d.qfrc_constraint+d.qfrc_passive);passive=E.T@(d.actuator_moment.T@d.actuator_force);target=gravity-support-passive;names=[m.actuator(i).name for i in range(m.nu)];rows=[]
 for family,label in (('s3k','s3k'),('syn6','syn6'),('w2lvar','w2lvar'),('w2lvar','fanin_mean')):
  # The legacy native log retains only20selected activations. Reconstruct all92
  # through the actual native activation kernel and gate it on those measured20.
  import xml.etree.ElementTree as ET
  os_family=__import__('os');os_family.environ['AARL_NET']=family;n=F.B.build(names,dt=.002,interleg=True);nnames=F.canonical_neuron_names(n);mn=np.array([n.idx['MN_'+a] for a in names])
  with np.load(OUT/(family+'_'+label+'_v1.apparatus.npz')) as z:commands=np.clip(z['neuro'][:,mn]/F.R.E_HI,0,1)
  tree=ET.fromstring('<mujoco><option timestep=".002" gravity="0 0 0"/><worldbody><body><joint name="q"/><geom type="sphere" size=".01" mass="1"/></body></worldbody><actuator/></mujoco>')
  for i in range(m.nu):ET.SubElement(tree.find('actuator'),'general',dict(name=names[i],joint='q',dyntype='muscle',dynprm=' '.join(format(x,'.17g') for x in m.actuator_dynprm[i]),gainprm='0',biastype='none',ctrllimited='true',ctrlrange='0 1'))
  am=F.mujoco.MjModel.from_xml_string(ET.tostring(tree,encoding='unicode'));ad=F.mujoco.MjData(am);activation_rows=[]
  for u in commands:
   ad.ctrl[:]=u;F.mujoco.mj_forward(am,ad);ad.act[:]+=am.opt.timestep*ad.act_dot;np.clip(ad.act,0,1,out=ad.act);activation_rows.append(ad.act.copy())
  activation_rows=np.asarray(activation_rows)
  state_name=family+'_'+label+('_v1.motor_state.npz' if label=='fanin_mean' else '_measured_v2.motor_state.npz')
  with np.load(OUT/state_name) as z:
   activation_error=float(np.max(abs(activation_rows-z['activation'])));control_error=float(np.max(abs(commands-z['excitation'])));observed_activation=z['activation'][-250:].mean(0)
  if activation_error>1e-12 or control_error>1e-12:raise ValueError('Native activation reconstruction mismatch: '+str((activation_error,control_error)))
  act=observed_activation
  d.act[:]=act;F.mujoco.mj_forward(m,d);contribution=E.T@(d.actuator_moment.T*d.actuator_force);grouped={}
  for group in sorted({classify(n).groups[0] for n in names}):grouped[group]=contribution[:,[i for i,n in enumerate(names) if classify(n).groups[0]==group]].sum(1).tolist()
  rows.append(dict(family=family if label==family else 'w2lvar_meanIIIN',mean_observed_activation=act.tolist(),native_all92_activation_reconstruction_error=activation_error,native_all92_control_error=control_error,common_pose_total_muscle_torque_Nm=contribution.sum(1).tolist(),primary_group_torque_Nm=grouped))
 role_dof={'hip_ext':('hip_flexion','glut_max2'),'hip_flex':('hip_flexion','psoas'),'knee_ext':('knee_angle','vas_int'),'knee_flex':('knee_angle','bifemsh'),'ankle_pf':('ankle_angle','soleus'),'ankle_df':('ankle_angle','tib_ant'),'hip_abd':('hip_adduction','glut_med1'),'hip_add':('hip_adduction','add_long')};roles=[]
 for i,name in enumerate(names):
  group=classify(name).groups[0]
  if group not in role_dof:continue
  joint,reference=role_dof[group];side=name[-1];row=JOINTS.index(joint+'_'+side);a=float(arms.arms[row,i]);ref=float(arms.arms[row,names.index(reference+'_'+side)])
  roles.append(dict(actuator=name,source_primary_group=group,joint=JOINTS[row],common_pose_moment_arm_m=a,reference_muscle=reference+'_'+side,reference_arm_m=ref,opposite_reference_action=bool(a*ref<0 and abs(a)>.0005)))
 report=dict(machine='easteregg2',joints=JOINTS,actuators=names,static_target_reconstruction_bit_exact=match,static_target_reconstruction_max_abs_error=float(np.max(abs(activation-baseline['standing_activation_targets']))),gravity_torque_Nm=gravity.tolist(),frozen_zero_activation_contact_and_passive_joint_torque_Nm=support.tolist(),passive_muscle_torque_Nm=passive.tolist(),desired_active_torque_Nm=target.tolist(),families=rows,primary_role_sign_audit=roles,scope='Identical corrected pose, zero velocity, equal biological force law and apparatus. Equality-aware projection includes moving-point followers. Actual late mean activation is replayed only for this OFFLINE force attribution at the common pose; not used as a simulation control override or teacher. Native activation replay is checked against all92 directly logged states. Legacy NPZ act is selected excitation, not activation. Primary-group allocation counts each muscle once while retaining all its joint-torque components, including RF hip action. The sign audit measures this pose only, relative to an explicit reference muscle; it does not silently reclassify muscles across the full RoM. Frozen contacts are not a proof of dynamic equilibrium.')
 (OUT/'COMMON_JOINT_TORQUE_REFERENCE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(targets_bit_exact=match,joints=len(JOINTS),families=[r['family'] for r in rows]),indent=2))

if __name__=='__main__':main()
