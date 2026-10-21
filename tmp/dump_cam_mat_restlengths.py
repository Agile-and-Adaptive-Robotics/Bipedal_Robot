"""Dump BPA resting lengths + tendon lengths from the cam result mat of record."""
import numpy as np
from scipy.io import loadmat

MAT = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_Result_cam_20261014_0416.mat"

m = loadmat(MAT, squeeze_me=True, struct_as_record=False)
print("== top-level keys ==")
for k, v in m.items():
    if k.startswith("__"):
        continue
    t = type(v).__name__
    try:
        shape = np.asarray(v).shape
    except Exception:
        shape = "?"
    print(f"  {k}: {t} shape={shape}")

def show(obj, name, depth=0, maxdepth=2):
    ind = "  " * (depth + 1)
    for attr in dir(obj):
        if attr.startswith("_") or callable(getattr(obj, attr, None)) and not isinstance(getattr(obj, attr, None), np.ndarray):
            continue
        try:
            val = getattr(obj, attr)
        except Exception:
            continue
        if isinstance(val, np.ndarray):
            s = np.array2string(np.asarray(val, dtype=float), precision=6, separator=", ")
            if len(s) > 200:
                s = s[:200] + " ..."
            print(f"{ind}{attr} = {s}")
        elif hasattr(val, "_fieldnames"):
            print(f"{ind}{attr}: struct")
            if depth < maxdepth:
                show(val, attr, depth + 1, maxdepth)
        else:
            print(f"{ind}{attr} = {val}")

for var in ("xFullBest", "xBestCam", "xFixed", "tendonRef"):
    if var in m:
        v = m[var]
        print(f"\n== {var} ==")
        print(np.array2string(np.asarray(v, dtype=float), precision=6, separator=", "))

for var in ("predBest", "ctxCamBest", "ctx"):
    if var in m:
        print(f"\n== {var} fields ==")
        show(m[var], var)
