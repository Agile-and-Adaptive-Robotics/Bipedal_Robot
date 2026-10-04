# Probe Vas_Pam_20mm_Result_20260925.mat for the -p8 prose-fix round (2026-10-03).
# Verifies: top-level removedText (+ plotting leftovers), routeInfo.active rows
# for p3/p8 (active from full flexion?), and eliminatedAngleD.
import numpy as np
import scipy.io as sio

MAT = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_Result_20260925.mat"
m = sio.loadmat(MAT, struct_as_record=False, squeeze_me=True)

keys = sorted(k for k in m if not k.startswith("__"))
print("top-level keys:", keys)

for k in ("removedText", "qPlot", "plotIdx", "tileTitle"):
    if k in m:
        print(f"{k} = {m[k]!r}")
    else:
        print(f"{k}: ABSENT")

ri = m.get("routeInfo")
if ri is None:
    raise SystemExit("no routeInfo")
print("routeInfo fields:", ri._fieldnames)

act = np.atleast_2d(np.asarray(ri.active))
print("active shape (rows=points, cols=frames):", act.shape)
elim = np.asarray(ri.eliminatedAngleD, dtype=float)
print("eliminatedAngleD (p1..p9):", np.round(elim, 2).tolist())

for f in ri._fieldnames:
    if "FullFlexion" in f or "activeAt" in f:
        v = getattr(ri, f)
        arr = np.asarray(v)
        print(f"{f}: shape={arr.shape} all-ones={bool(np.all(arr.astype(bool)))}")

# p8 = row index 7, p3 = row index 2 (1-based pN -> index N-1)
for name, row in (("p3", 2), ("p8", 7)):
    a = act[row, :].astype(bool)
    off = np.flatnonzero(~a)
    first_off = int(off[0]) if off.size else None
    n_on = int(a.sum())
    print(f"{name}: active frames {n_on}/{a.size}, first-off frame {first_off}, "
          f"active at frame 0 (full flexion) = {bool(a[0])}, "
          f"active at final frame = {bool(a[-1])}")

# fraction of grid with p8 active (frames 0..first_off-1)
a8 = act[7, :].astype(bool)
off8 = np.flatnonzero(~a8)
if off8.size:
    frac = off8[0] / a8.size
    print(f"p8 inactive starting at frame {off8[0]} of {a8.size} "
          f"-> active over {100.0 * (1.0 - frac):.1f}% of frames")
