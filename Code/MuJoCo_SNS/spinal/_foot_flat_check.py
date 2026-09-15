"""Foot-flat diagnostic v2: mesh-aware sole height at the normal.mot pose
as a function of pelvis height."""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import mujoco

import runner as R

model = mujoco.MjModel.from_xml_path(str(R.MODEL))
data = mujoco.MjData(model)


def geom_bottom_z(g):
    t = model.geom_type[g]
    if t == 7:  # mesh: min over transformed vertices
        m = model.geom_dataid[g]
        va, vn = model.mesh_vertadr[m], model.mesh_vertnum[m]
        verts = model.mesh_vert[va:va + vn].astype(float)
        Rm = data.geom_xmat[g].reshape(3, 3)
        wz = (verts @ Rm.T)[:, 2] + data.geom_xpos[g][2]
        return float(wz.min())
    r = model.geom_size[g][2] if model.geom_type[g] in (2, 3, 6) else 0.0
    return float(data.geom_xpos[g][2] - r)


def body_bottom_z(body_name):
    bid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, body_name)
    if bid < 0:
        return float("nan")
    z = 1e9
    for g in range(model.ngeom):
        if model.geom_bodyid[g] == bid:
            z = min(z, geom_bottom_z(g))
    # include child bodies' geoms (toes may parent nothing, but calcn may)
    for b in range(model.nbody):
        p = model.body_parentid[b]
        if (p == bid or b == bid) and bid >= 0:
            for g in range(model.ngeom):
                if model.geom_bodyid[g] == b and \
                        not (model.body_weldid[b] == 0 and b != bid):
                    pass
    return z


def side_bottom(side):
    z = 1e9
    for bn in (f"toes_{side}", f"calcn_{side}", f"talus_{side}",
               f"foot_{side}", f"bofoot_{side}"):
        bid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, bn)
        if bid < 0:
            continue
        for g in range(model.ngeom):
            if model.geom_bodyid[g] == bid:
                z = min(z, geom_bottom_z(g))
    return z


def set_pose(height):
    mujoco.mj_resetDataKeyframe(model, data, 0)
    jadr = {jn: model.joint(jn).qposadr[0] for jn in R.START_POSE_DEG
            if model.joint(jn).id >= 0}
    for jn, deg in R.START_POSE_DEG.items():
        if jn in jadr:
            data.qpos[jadr[jn]] = np.radians(deg)
    data.qpos[model.joint("pelvis_ty").qposadr[0]] = height
    R._project_followers(model, data)
    mujoco.mj_forward(model, data)


print("height  sole_r   sole_l   (target: soles at 0, ankles -1.7/+9.8)")
rows = []
for h in np.arange(0.96, 0.79, -0.02):
    set_pose(h)
    sr, sl = side_bottom("r"), side_bottom("l")
    rows.append((h, sr, sl))
    print(f"{h:5.2f}  {sr:+7.3f}  {sl:+7.3f}")
# interpolate the height where the lower sole crosses 0
lo, hi = rows[-1], rows[0]
for (h1, s1, _), (h2, s2, _) in zip(rows, rows[1:]):
    if s1 < 0 <= s2:
        frac = (0 - s1) / (s2 - s1)
        print(f"\nright sole touches 0 at pelvis_ty ~ "
              f"{h1 + frac * (h2 - h1):.3f}")
        break
