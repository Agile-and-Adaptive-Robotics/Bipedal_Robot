"""Full-image scan: any bluish or strongly non-white pixels anywhere?"""
from pathlib import Path
import numpy as np
from PIL import Image

p = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\30-torque\opensim_robot_isometric.png")
a = np.asarray(Image.open(p).convert("RGB")).astype(int)
h, w, _ = a.shape
R, G, B = a[:, :, 0], a[:, :, 1], a[:, :, 2]

bluish = (B > R + 15)  # any pixel where blue dominates red = colored band
print("bluish pixel count:", int(bluish.sum()))
if bluish.any():
    rows = np.where(bluish.any(axis=1))[0]
    cols = np.where(bluish.any(axis=0))[0]
    print("  rows", rows.min(), "..", rows.max(), "cols", cols.min(), "..", cols.max())
    print("  sample colors:", a[bluish][:10].tolist())

nonwhite = ~(a > 246).all(axis=2)
rows = np.where(nonwhite.any(axis=1))[0]
print("rows containing ANY non-near-white pixel:", rows.min(), "..", rows.max(),
      f"({len(rows)} of {h} rows)")
# column extent of nonwhite in top half (band check)
top = nonwhite[: h // 2]
tr = np.where(top.any(axis=1))[0]
if len(tr):
    band = top[tr[0]]
    print(f"first nonwhite row {tr[0]}: {int(band.sum())} px, cols "
          f"{int(np.argmax(band))}..{len(band)-1-int(np.argmax(band[::-1]))}")
