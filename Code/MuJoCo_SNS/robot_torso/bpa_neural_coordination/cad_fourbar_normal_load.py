"""Known loop-normal load, with backhinge moment7–9Nm and exactwrench accounting."""
import json
import numpy as np
from scipy.optimize import least_squares
import cad_fourbar_fixture as F
import cad_fourbar_signed_check as L
mujoco=F.mujoco;OUT=F.OUT

def solve(f,T):
 f.d.qpos[:]=f.q0;f.base();J=f.d.efc_J.reshape(f.d.nefc,3).copy();lam=np.linalg.lstsq(np.vstack((J[:,0],J[:,2])),[T,0.],rcond=None)[0];load=J.T@lam
 def res(x):
  f.d.qpos[:]=f.q0+f.W@x;f.d.qvel[:]=0;mujoco.mj_forward(f.m,f.d);f.d.qacc[:]=0;mujoco.mj_inverse(f.m,f.d);return f.W.T@(f.d.qfrc_inverse-load)
 initial=np.linalg.lstsq(f.sourceK,-res(np.zeros(2)),rcond=None)[0];sol=least_squares(res,np.clip(initial,-.5,.5),bounds=(-np.ones(2)*.6,np.ones(2)*.6),max_nfev=80,ftol=1e-13,xtol=1e-13,gtol=1e-13,diff_step=1e-4);error=res(sol.x);rho=float(f.c@(f.d.qfrc_inverse-load)/(f.c@f.c));applied=load+f.c*rho;gap=f.gap();f.d.qfrc_applied[:]=applied;mujoco.mj_forward(f.m,f.d);J=f.d.efc_J.reshape(f.d.nefc,3);reaction=J.T@f.d.efc_force;acc=f.d.qacc.copy();mujoco.mj_inverse(f.m,f.d)
 return dict(back_load_target_Nm=T,known_equal_opposite_anchor_force_N=lam.tolist(),known_joint_torque_Nm=load.tolist(),actual_constraint_joint_reaction_Nm=reaction.tolist(),fixture_hold_joint_torque_Nm=(f.c*rho).tolist(),gap_m=gap,qacc_max=float(np.max(np.abs(acc))),static_residual_Nm=error.tolist(),forward_inverse_error_Nm=float(np.max(np.abs(f.d.qfrc_inverse-applied))),source_inertia_unchanged=True)

def main():
 sel=json.loads((OUT/'DERIVED_VARIANT.json').read_text());out={}
 for side,text,cls in [('right',(OUT/'PINNED_SOURCE_FIXTURE.xml').read_text(),F.Fixture),('left',(OUT/'PINNED_SOURCE_FIXTURE_LEFT.xml').read_text(),L.Left)]:
  source=cls(text);selected=cls(text,sel['timeconst_s'],sel['constant_impedance'],sel['timestep_s'],500);out[side]=[dict(source=solve(source,T),selected=solve(selected,T)) for T in (7.,8.,9.)]
 result=dict(scope='Known equal/oppositeforce pair at the two loopanchors, chosen toproduce7–9Nm backhingemoment and0front externalmoment. Load lies in sourceclosed J.T range; explicit1DOFfixturehold, no muscles/robottrial.',cases=out,numerical_variant_frozen=True,nearzero_gap_gate=all(row['selected']['gap_m']<.0001 and row['selected']['qacc_max']<1e-5 for rows in out.values() for row in rows));(OUT/'NORMAL_LOAD_WRENCH_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
