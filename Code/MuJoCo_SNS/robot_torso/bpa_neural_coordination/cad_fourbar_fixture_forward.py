"""Bounded fixture load-step convergence, not a robot standing simulation."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
import cad_fourbar_fixture as F
OUT=F.OUT;mujoco=F.mujoco

def run_case(tag,dt,iterations,source=False):
 text=(OUT/'PINNED_SOURCE_FIXTURE.xml').read_text();selection=json.loads((OUT/'DERIVED_VARIANT.json').read_text());f=F.Fixture(text,None if source else selection['timeconst_s'],None if source else selection['constant_impedance'],dt,iterations);zero,q,applied=f.equilibrium(0.);f.d.qpos[:]=q;f.d.qvel[:]=0;f.d.qfrc_applied[:]=applied;mujoco.mj_forward(f.m,f.d);J=f.d.efc_J.reshape(f.d.nefc,3);n=null_space(J).ravel();vtheta=n/(f.c@n);M=np.zeros((3,3));mujoco.mj_fullM(f.m,M,f.d.qM);effective=float(vtheta@M@vtheta);kp=effective/.05**2;kd=2*effective/.05;rows=[];failure=None
 for k in range(round(1./dt)):
  t=k*dt;load=0. if t<.1 else 8.;mujoco.mj_forward(f.m,f.d);J=f.d.efc_J.reshape(f.d.nefc,3);n=null_space(J).ravel();external=np.array([load,0,0]);rho=float(n@(f.d.qfrc_bias-external)/(n@f.c));error=float(f.c@(f.d.qpos-f.q0));rate=float(f.c@f.d.qvel);hold=rho-kp*error-kd*rate;f.d.qfrc_applied[:]=external+f.c*hold;mujoco.mj_step(f.m,f.d);mujoco.mj_forward(f.m,f.d);rows.append(dict(t=float(f.d.time),gap=f.gap(),qpos=f.d.qpos.copy(),qvel=f.d.qvel.copy(),constraint_force=f.d.efc_force.copy(),applied=f.d.qfrc_applied.copy(),angle_error=float(f.c@(f.d.qpos-f.q0)),load=load))
  if not np.isfinite(f.d.qpos).all() or np.max(np.abs(f.d.qvel))>100 or f.gap()>.05:failure='Fixturefinite/velocity/50mmgap stop';break
 arrays={key:np.asarray([r[key] for r in rows]) for key in rows[0]};tail=arrays['t']>=.8;row8,_,_=f.equilibrium(8.);result=dict(tag=tag,scope='Pinnedsource-inertiafixture loadstep0→8Nm at0.1s, not robot/SNS/BPA',duration_s=float(arrays['t'][-1]),dt=dt,iterations=iterations,source_constraints=source,failure=failure,max_gap_m=float(arrays['gap'].max()),tail_mean_gap_m=float(arrays['gap'][tail].mean()) if tail.any() else None,tail_max_qvel=float(np.abs(arrays['qvel'][tail]).max()) if tail.any() else None,angle_error_max_rad=float(np.max(np.abs(arrays['angle_error']))),static8_gap_m=row8['actual_gap_m'],fixture_hold_kp_Nm_rad=kp,fixture_hold_kd_Nms_rad=kd,fixture_effective_inertia_kg_m2=effective,hold_policy='Explicitexternal1DOFfixturehold: tangentialgravity/load feedforward + criticalPD50mstimeconstant; sourcejointdamping/mass unchanged.')
 np.savez_compressed(OUT/(tag+'_TRACE.npz'),**arrays);(OUT/(tag+'_RESULT.json')).write_text(json.dumps(result,indent=2));return result

def main():
 selection=json.loads((OUT/'DERIVED_VARIANT.json').read_text());dt=selection['timestep_s'];p=dict(scope='Fixtureonly loadstepconvergence',known_load_Nm=8.,step_time_s=.1,duration_s=1.,sourcecase=True,variants=['selected_dt100iter','selected_dt500iter','halfdT500iter','quarterdT500iter'],source_mass_inertia_unchanged=True,derivedvariant_frozen=True,gate_max_gap_m=.0005,gate_tail_dt_relative_difference=.05,no_robot_trial=True,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),Path(F.__file__),OUT/'DERIVED_VARIANT.json')});file=OUT/'FORWARD_PREREGISTRATION.json'
 if file.exists():assert json.loads(file.read_text())==p
 else:file.write_text(json.dumps(p,indent=2))
 cases=[run_case('source_dt',.002,100,True),run_case('selected_dt',dt,100),run_case('selected_500iter',dt,500),run_case('halfstep',dt/2,500),run_case('quarterstep',dt/4,500)];selected=cases[1:];tail=[r['tail_mean_gap_m'] for r in selected];relative=abs(tail[-1]-tail[-2])/max(abs(tail[-1]),1e-12);result=dict(cases=cases,selected_fixture_gate=all(r['failure'] is None and r['max_gap_m']<=.0005 for r in selected),half_quarter_tail_gap_relative_difference=relative,time_convergence_pass=relative<=.05,no_robot_trial_performed=True);(OUT/'FIXTURE_FORWARD_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
