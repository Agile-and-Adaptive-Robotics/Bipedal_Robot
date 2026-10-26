r"""Sagittal-locked gait2392 plant for the Di Russo replication (D2).

"Constrain to sagittal plane like she did": her model's world is 2D —
fore-aft, vertical, and pitch free; everything else locked. The converted
model exposes named joints, so the lock is EXACT (equality joint
constraints to zero, no force leashes):
  locked: pelvis_list, pelvis_rotation, pelvis_ty (frontal/yaw/lateral),
          hip_adduction_r/l, hip_rotation_r/l, subtalar_angle_r/l,
          lumbar_extension/bending/rotation (trunk = single segment).
  free:   pelvis_tx (fore-aft), pelvis_tz (vertical), pelvis_tilt (pitch),
          hip_flexion, knee_angle, ankle_angle, mtp_angle per side.
Muscle restriction: her 9 groups/leg keep force; every other actuator is
zeroed (active gainprm[2] AND passive biasprm[2] — the runner prune-repair
pattern) so the model carries exactly her muscle set.
"""
from __future__ import annotations

import re

LOCK_JOINTS = [
    # NOTE axis bake (DESIGN trap #1): pelvis_ty = world VERTICAL (free,
    # harness channel); pelvis_tz = world LATERAL (locked).
    "pelvis_list", "pelvis_rotation", "pelvis_tz",
    "hip_adduction_r", "hip_adduction_l",
    "hip_rotation_r", "hip_rotation_l",
    "subtalar_angle_r", "subtalar_angle_l",
    "lumbar_extension", "lumbar_bending", "lumbar_rotation",
]

HER_ACTUATORS = set()
for _grp, _names in {
    "ILPSO": ["iliacus", "psoas"],
    "GMAX": ["glut_max1", "glut_max2", "glut_max3"],
    "HAMS": ["semimem", "semiten", "bifemlh"],
    "BFSH": ["bifemsh"],
    "RF": ["rect_fem"],
    "VAS": ["vas_med", "vas_int", "vas_lat"],
    "GAS": ["med_gas", "lat_gas"],
    "SOL": ["soleus"],
    "TA": ["tib_ant"],
}.items():
    for _n in _names:
        HER_ACTUATORS |= {_n + "_r", _n + "_l"}


def build_sagittal_xml(src_xml: str) -> str:
    txt = open(src_xml, encoding="utf-8").read()
    # 1) exact sagittal lock via joint equalities (append into the existing
    #    <equality> section)
    eq_add = "".join(
        f'  <joint name="lock_{j}" joint1="{j}"/>' for j in LOCK_JOINTS)
    m = re.search(r"<equality>", txt)
    if not m:
        txt = txt.replace("</mujoco>", "<equality>\n" + eq_add +
                          "\n</equality>\n</mujoco>")
    else:
        txt = txt[:m.end()] + "\n" + eq_add + txt[m.end():]
    # 2) muscle restriction: zero force+passive on every non-her actuator
    n_zero = 0

    def _zero(match):
        nonlocal n_zero
        tag = match.group(0)
        name = re.search(r'name="([^"]+)"', tag).group(1)
        if name in HER_ACTUATORS:
            return tag
        n_zero += 1
        for prm in ("gainprm", "biasprm"):
            mm = re.search(prm + r'="[^"]+"', tag)
            if mm:
                vals = mm.group(0).split('"')[1].split()
                vals[2] = "0"
                tag = tag[:mm.start()] + prm + '="' + " ".join(vals) + \
                    '"' + tag[mm.end():]
        return tag

    txt = re.sub(r'<general [^>]*name="[^"]+"[^>]*/>', _zero, txt)
    return txt, n_zero


def index_maps(model, mujoco):
    """Group->actuator ids, foot geoms, pelvis body, tilt joint qpos adr."""
    from controller import MAP, MUSCLES
    ids = {}
    for m in MUSCLES:
        for s in ("r", "l"):
            names = [t.format(s=s) for t in MAP[m]]
            ids[(m, s.upper())] = [mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_ACTUATOR, n) for n in names]
    geoms = {}
    for s in ("r", "l"):
        geoms[("heel", s.upper())] = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_GEOM, f"calcn_{s}_geom_1")
        geoms[("toe", s.upper())] = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_GEOM, f"toes_{s}_geom_1")
    pelvis = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "pelvis")
    jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "pelvis_tilt")
    tilt_adr = model.jnt_qposadr[jid]
    tx = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "pelvis_tx")
    tx_adr = model.jnt_qposadr[tx]
    return ids, geoms, pelvis, tilt_adr, tx_adr
