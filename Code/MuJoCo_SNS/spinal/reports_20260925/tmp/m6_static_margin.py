"""M6 static-balance margin: where is the COM relative to the ground
support polygon (contact plates + foot/toe body geoms) at the grounded
rest pose? A rigid statue stands only if the COM x-y is inside the
polygon."""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, HERE)

import numpy as np
import mujoco

m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "w2l_ground.xml"))
d = mujoco.MjData(m)
jadr = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "root")]
d.qpos[jadr + 2] -= 0.026
mujoco.mj_forward(m, d)

com = d.subtree_com[0][:2]
print(f"COM world (x, y) = ({com[0]:+.4f}, {com[1]:+.4f})")

# every geom whose bottom face plausibly supports: the 4 plates + 2 foot
# boxes + 2 toe boxes; project each box's 8 corners to ground plane and
# keep those within 3 mm of the lowest surface.
geoms = ("foot_L_contact", "foot_R_contact", "toe_L_contact", "toe_R_contact",
         "foot_L", "foot_R", "toe_L", "toe_R")
corners = np.zeros((0, 2))
for g in geoms:
    gi = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g)
    R = d.geom_xmat[gi].reshape(3, 3)
    sz = m.geom_size[gi]
    loc = np.stack([np.array([-1, -1, -1, -1, 1, 1, 1, 1], float) * sz[0],
                    np.array([-1, -1, 1, 1, -1, -1, 1, 1], float) * sz[1],
                    np.array([-1, 1, -1, 1, -1, 1, -1, 1], float) * sz[2]])
    world = d.geom_xpos[gi][:, None] + R @ loc
    print(f"  {g:16s} lowest corner z = {world[2].min():+.4f}")
    keep = world[2] < world[2].min() + 0.003
    corners = np.vstack([corners, world[:2, keep].T])

x0, x1 = corners[:, 0].min(), corners[:, 0].max()
y0, y1 = corners[:, 1].min(), corners[:, 1].max()
print(f"\nsupport bbox (plates+feet+toes, both sides): "
      f"x [{x0:+.4f}, {x1:+.4f}] ({(x1-x0)*100:.1f} cm)  "
      f"y [{y0:+.4f}, {y1:+.4f}] ({(y1-y0)*100:.1f} cm)")
print(f"COM x inside x-range? {x0 <= com[0] <= x1}   "
      f"margin to front edge {100*(x1-com[0]):+.2f} cm, "
      f"to rear edge {100*(com[0]-x0):+.2f} cm")
print(f"COM y inside y-range? {y0 <= com[1] <= y1}   "
      f"margin {min(100*(com[1]-y0), 100*(y1-com[1])):+.2f} cm")

# per-side support (single-stance check)
for side, sfx in (("L", "L"), ("R", "R")):
    sel = [g for g in geoms if g.endswith(sfx)]
    cc = np.zeros((0, 2))
    for g in sel:
        gi = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g)
        R = d.geom_xmat[gi].reshape(3, 3)
        sz = m.geom_size[gi]
        loc = np.stack([np.array([-1, -1, -1, -1, 1, 1, 1, 1], float) * sz[0],
                        np.array([-1, -1, 1, 1, -1, -1, 1, 1], float) * sz[1],
                        np.array([-1, 1, -1, 1, -1, 1, -1, 1], float) * sz[2]])
        w = d.geom_xpos[gi][:, None] + R @ loc
        keep = w[2] < w[2].min() + 0.003
        cc = np.vstack([cc, w[:2, keep].T])
    print(f"  {side} support: x [{cc[:,0].min():+.4f}, {cc[:,0].max():+.4f}] "
          f"({(cc[:,0].max()-cc[:,0].min())*100:.1f} cm)  "
          f"y [{cc[:,1].min():+.4f}, {cc[:,1].max():+.4f}]")

# also: pelvis pitch at rest (rest pose is pitched vs world)
root = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, "Root")
from mujoco import mju_mat2Quat
q = np.zeros(4)
mju_mat2Quat(q, d.xmat[root].reshape(9))
print(f"root quat at rest: {np.round(q, 5)} (yaw about MJ y = pitch)")
