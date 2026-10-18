"""Probe the top band of the render + check skybox texture state."""
from pathlib import Path
import numpy as np
from PIL import Image
import mujoco

p = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\30-torque\opensim_robot_isometric.png")
a = np.asarray(Image.open(p).convert("RGB")).astype(int)

print("== row profile: first non-near-white row per column sample ==")
nonwhite = ~(a > 246).all(axis=2)
for r in [0, 2, 5, 10, 20, 40, 80, 120, 160, 200, 240, 300, 400]:
    row = nonwhite[r]
    if row.any():
        cmin = int(np.argmax(row)); cmax = int(len(row) - 1 - np.argmax(row[::-1]))
        seg = a[r, cmin:cmax + 1]
        print(f"row {r:4d}: nonwhite cols {cmin}..{cmax} ({row.sum()} px), mean color there {seg.reshape(-1,3).mean(axis=0).round(0).tolist()}")
    else:
        print(f"row {r:4d}: all near-white")

print("== color at exact probe points ==")
for (r, c) in [(5, 5), (5, 720), (100, 5), (150, 720), (60, 200)]:
    print(f"({r},{c}): {a[r, c].tolist()}")

print("== model texture table ==")
m = mujoco.MjModel.from_xml_path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_robot\gait2392_robot_cvt3.xml")
print("ntex:", m.ntex)
off = 0
for i in range(m.ntex):
    n = int(m.tex_height[i]) * int(m.tex_width[i]) * 3
    seg = m.tex_rgb[off:off + n]
    print(f"tex {i}: type={m.tex_type[i]} w={m.tex_width[i]} h={m.tex_height[i]} "
          f"rgb range [{seg.min():.3f},{seg.max():.3f}] mean={seg.mean():.3f}")
    off += n
