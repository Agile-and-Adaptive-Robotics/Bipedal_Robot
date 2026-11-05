"""Verify Gait2392 knee FK: closed-form CustomJoint vs converted MuJoCo model.

Closed form (from gait2392_simbody.osim knee_r, validated structure):
    child origin in femur frame: p(q) = (f1(q), f2(q), 0)   [SimmSpline tables]
    child orientation: Rz(q)
f1 = translation1 table (x), f2 = translation2 table (y).
Comparison uses the converted model's quartic equality fits (approximation of
the same splines), so expect sub-mm agreement, exact rotation.
"""
import numpy as np
import mujoco
from scipy.interpolate import CubicSpline

XML = (r"C:/Users/Ben/Documents/GitHub/Bipedal_Robot/"
       r"Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/"
       r"gait2392_simbody_cvt3.xml")

# SimmSpline tables from the osim (knee_r translation1/translation2)
x1 = np.array([-2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533,
               0.197344, 0.337395, 0.490178, 1.52146, 2.0944])
y1 = np.array([-0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031,
               -0.005227, -0.005435, -0.005574, -0.005435, -0.00525])
x2 = np.array([-2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944])
y2 = np.array([-0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396])
# OpenSim SimmSpline = natural cubic spline
f1 = CubicSpline(x1, y1, bc_type="natural")
f2 = CubicSpline(x2, y2, bc_type="natural")


def closed_form(q):
    p = np.array([float(f1(q)), float(f2(q)), 0.0])
    c, s = np.cos(q), np.sin(q)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return p, R


model = mujoco.MjModel.from_xml_path(XML)
data = mujoco.MjData(model)

jid_femur = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "femur_r")
jid_tibia = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "tibia_r")
j_tr1 = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "knee_r_translation1")
j_tr2 = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "knee_r_translation2")
j_ang = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "knee_angle_r")

# quartic fits from the equality constraints
pc1 = np.array([-0.003639, -0.006267, 0.002807, 0.00134, -0.0006776])
pc2 = np.array([-0.3957, 0.003413, -0.005587, 0.0006697, 0.0005674])
qadr1 = model.jnt_qposadr[j_tr1]
qadr2 = model.jnt_qposadr[j_tr2]
qadra = model.jnt_qposadr[j_ang]

mujoco.mj_resetData(model, data)
mujoco.mj_forward(model, data)
# pelvis may be a free joint: pin everything at defaults, then zero pelvis vel
qpos0 = data.qpos.copy()

print(" q(rad)   mujoco_tibia_in_femur                    closed_form                         |dp| mm  dR deg")
for q in (0.0, -0.3, -0.7, -1.2, -1.7, -2.0944):
    mujoco.mj_resetData(model, data)
    data.qpos[:] = qpos0
    data.qpos[qadr1] = np.polyval(pc1[::-1], q)
    data.qpos[qadr2] = np.polyval(pc2[::-1], q)
    data.qpos[qadra] = q
    mujoco.mj_forward(model, data)
    p_mjc_world = data.xpos[jid_tibia].copy()
    R_mjc_world = data.xmat[jid_tibia].reshape(3, 3).copy()
    p_f_world = data.xpos[jid_femur].copy()
    R_f_world = data.xmat[jid_femur].reshape(3, 3).copy()
    R_fT = R_f_world.T
    p_mjc = R_fT @ (p_mjc_world - p_f_world)
    R_mjc = R_fT @ R_mjc_world

    p_cf, R_cf = closed_form(q)
    dp = np.linalg.norm(p_mjc - p_cf) * 1000
    dR = np.degrees(np.arccos(np.clip((np.trace(R_mjc.T @ R_cf) - 1) / 2, -1, 1)))
    print(f"{q:+7.3f}  {np.array2string(p_mjc, precision=5)}  {np.array2string(p_cf, precision=5)}  "
          f"{dp:7.3f}  {dR:7.3f}")

print("\nDONE_FK_CHECK")
