"""Worst validated left9Nm loopnormalcouple loadstep convergence, fixtureonly."""
import json
import numpy as np
from scipy.linalg import null_space
import cad_fourbar_fixture as F
import cad_fourbar_signed_check as L
mujoco=F.mujoco;OUT=F.OUT

def run(tag,dt,iterations):
 sel=json.loads((OUT/'DERIVED_VARIANT.json').read_text());normal=json.loads((OUT/'NORMAL_LOAD_WRENCH_AUDIT.json').read_text());external=np.array(normal['cases']['left'][-1]['selected']['known_joint_torque_Nm']);f=L.Left((OUT/'PINNED_SOURCE_FIXTURE_LEFT.xml').read_text(),sel['timeconst_s'],sel['constant_impedance'],dt,iterations);_,q,Q=f.equilibrium(0);f.d.qpos[:]=q;f.d.qvel[:]=0;f.d.qfrc_applied[:]=Q;mujoco.mj_forward(f.m,f.d);J=f.d.efc_J.reshape(f.d.nefc,3);n=null_space(J,rcond=1e-8).ravel();v=n/(f.c@n);M=np.zeros((3,3));mujoco.mj_fullM(f.m,M,f.d.qM);I=float(v@M@v);kp=I/.05**2;kd=2*I/.05;rows=[];failure=None
 for k in range(round(1/dt)):
  t=k*dt;load=external if t>=.1 else np.zeros(3);mujoco.mj_forward(f.m,f.d);n=null_space(f.d.efc_J.reshape(f.d.nefc,3),rcond=1e-8).ravel();rho=float(n@(f.d.qfrc_bias-load)/(n@f.c));error=float(f.c@(f.d.qpos-f.q0));hold=rho-kp*error-kd*float(f.c@f.d.qvel);f.d.qfrc_applied[:]=load+f.c*hold;mujoco.mj_step(f.m,f.d);mujoco.mj_forward(f.m,f.d);rows.append(dict(t=float(f.d.time),gap=f.gap(),qpos=f.d.qpos.copy(),qvel=f.d.qvel.copy(),force=f.d.efc_force.copy(),applied=f.d.qfrc_applied.copy(),angle_error=float(f.c@(f.d.qpos-f.q0))))
  if not np.isfinite(f.d.qpos).all() or f.gap()>.005:failure='Nonfiniteor5mmfixturegap';break
 a={key:np.asarray([r[key] for r in rows]) for key in rows[0]};tail=a['t']>=.8;result=dict(tag=tag,scope='Left pinnedfixture9Nmback/13.05Nmshankconstantjointcouple, external1DOFhold; notrobot',dt=dt,iterations=iterations,failure=failure,maxgap_m=float(a['gap'].max()),tailmean_gap_m=float(a['gap'][tail].mean()) if tail.any() else None,tailmax_qvel=float(np.abs(a['qvel'][tail]).max()) if tail.any() else None,known_joint_load_Nm=external.tolist(),hold_kp=kp,hold_kd=kd);np.savez_compressed(OUT/(tag+'_TRACE.npz'),**a);return result

def main():
 sel=json.loads((OUT/'DERIVED_VARIANT.json').read_text());h=sel['timestep_s'];cases=[run('normal_left_dt',h,100),run('normal_left_500iter',h,500),run('normal_left_halfstep',h/2,500),run('normal_left_quarterstep',h/4,500)];diff=abs(cases[-1]['tailmean_gap_m']-cases[-2]['tailmean_gap_m'])/abs(cases[-1]['tailmean_gap_m']);result=dict(cases=cases,source_variant_unchanged=True,physical_model_unchanged=True,half_quarter_relative_difference=diff,half_quarter_converged=diff<.05,parent_05mm_gate=all(row['failure'] is None and row['maxgap_m']<.0005 for row in cases),original_100um_target_gate=all(row['failure'] is None and row['maxgap_m']<.0001 for row in cases),no_robot_trial=True);(OUT/'NORMAL_LOAD_FORWARD_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
