"""Axis-sign audit: does each muscle pull its joint in the anatomically
correct direction in the converted MJCF?

Method: at the keyframe pose, MuJoCo's actuator_moment[m, dof] gives the
generalized torque per unit muscle force. A positive torque on a hinge DoF
drives the joint in its positive direction. We compare that sign against the
anatomical expectation (per this file's table) for muscles whose action is
unambiguous.

!! CAVEAT (2026-09-10, after _muscle_direction_test.py): this method
injects +moment*Fmax as qfrc_applied, which measures the MOMENT-ARM sign
(dL/dq), NOT the direction the joint turns under real muscle activation -
the two come out SIGN-INVERTED here (e.g. it reports vas_lat pulling
knee_angle negative, while real activation of vas_lat drives knee_angle
POSITIVE). The table's expected signs are calibrated to this method, so
the audit stays internally consistent for catching routing flips, but the
signs are NOT OpenSim coordinate conventions. Ground truth for the knee
(real ctrl=1 activation tests, unjammed followers): knee_angle NEGATIVE =
flexion (OpenSim flexion-negative convention, preserved by the converter;
gravity buckles the standing knee negative, semimem/bifemsh/med_gas drive
negative, vas_lat/rect_fem drive positive = extension, +10 deg limit).
The converted knee_angle joint also ships limited="false" (range inert).

Verdict per joint: which fraction of anchor muscles agree; a flipped joint
shows up as near-total disagreement.
"""
from pathlib import Path

import numpy as np
import mujoco

MODEL = Path(r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
             r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml")

# (muscle, joint, expected torque sign in the joint's positive direction)
# anatomy from the gait2392 documentation (Delp 1990 / OpenSim docs)
ANCHORS = [
    ("soleus_r",   "ankle_angle_r",    +1),  # plantarflexion (+ in gait2392)
    ("tib_ant_r",  "ankle_angle_r",    -1),  # dorsiflexion
    ("med_gas_r",  "ankle_angle_r",    +1),
    ("vas_lat_r",  "knee_angle_r",     -1),  # knee extension (=- flexion)
    ("rect_fem_r", "knee_angle_r",     -1),
    ("bifemsh_r",  "knee_angle_r",     +1),  # knee flexion
    ("semimem_r",  "knee_angle_r",     +1),
    ("glut_max2_r", "hip_flexion_r",   -1),  # hip extension
    ("iliacus_r",  "hip_flexion_r",    +1),  # hip flexion
    ("psoas_r",    "hip_flexion_r",    +1),
    ("sar_r",      "hip_flexion_r",    +1),
    ("add_long_r", "hip_adduction_r",  +1),  # adduction (+)
    ("glut_med1_r", "hip_adduction_r", -1),  # abduction
    ("tib_ant_r",  "subtalar_angle_r", -1),  # inversion
    ("ercspn_r",   "lumbar_extension", -1),  # lumbar extension (=- flexion)
    # left side (mirrored frames)
    ("soleus_l",   "ankle_angle_l",    +1),
    ("vas_lat_l",  "knee_angle_l",     -1),
    ("iliacus_l",  "hip_flexion_l",    +1),
    ("glut_max2_l", "hip_flexion_l",   -1),
    ("add_long_l", "hip_adduction_l",  +1),
    ("glut_med1_l", "hip_adduction_l", -1),
]


def main():
    model = mujoco.MjModel.from_xml_path(str(MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    # audit the REPAIRED model (hip axis flips, prunes; pathpoint couplings
    # kept — welding them breaks moment arms)
    import runner as R
    key_pose = R.capture_pose(model)
    model = R.apply_harness(model, data)
    data = mujoco.MjData(model)
    R.seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)

    aid = {mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i): i
           for i in range(model.nu)}

    def joint_acc_with(m_idx: int) -> np.ndarray:
        """Delta joint acceleration from fully activating one muscle."""
        base = data.qacc.copy()
        f = float(model.actuator_gainprm[m_idx, 2])
        data.qfrc_applied[:] = data.actuator_moment[m_idx, :] * f
        mujoco.mj_forward(model, data)
        dq = data.qacc - base
        data.qfrc_applied[:] = 0.0
        mujoco.mj_forward(model, data)
        return dq

    print(f"{'muscle':14s} {'joint':18s} {'dq/dt2':>10s}  {'expected':>8s}  verdict")
    by_joint = {}
    for mus, jnt, want in ANCHORS:
        j = model.joint(jnt)
        dof = j.dofadr[0]
        m = aid[mus]
        if float(model.actuator_gainprm[m, 2]) == 0.0:
            print(f"{mus:14s} {jnt:18s} {'pruned':>10s}")
            continue
        dq = joint_acc_with(m)[dof]
        sign = np.sign(dq) if abs(dq) > 1e-6 else 0.0
        ok = (sign == want)
        by_joint.setdefault(jnt, []).append(ok)
        print(f"{mus:14s} {jnt:18s} {dq:10.4f}  {'+' if want > 0 else '-':>8s}  "
              f"{'OK' if ok else 'FLIPPED?' if sign != 0 else 'zero'}")
    print("\n-- per-joint agreement --")
    for jnt, oks in by_joint.items():
        print(f"{jnt:18s} {sum(oks)}/{len(oks)} anchors agree "
              f"{'  <-- AXIS LIKELY FLIPPED' if sum(oks) < len(oks) / 2 else ''}")


if __name__ == "__main__":
    main()
