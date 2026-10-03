"""Focused: dump XiUsed / lock fields from the redesign mats, recursively but skipping MCOS."""
import scipy.io as sio
import numpy as np
from pathlib import Path

MO = Path(r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results")
KEYS = ("xiused", "locksrc", "flxrow", "xi1lock", "xi2lock", "xi3", "removedtext", "fbest")

def walk(name, obj, depth=0):
    if depth > 4:
        return
    if isinstance(obj, np.ndarray) and obj.dtype.names:
        for f in obj.dtype.names:
            if obj.size:
                walk(f, obj.flat[0][f], depth + 1)
    elif isinstance(obj, np.ndarray) and obj.dtype == object and obj.size:
        for i, e in enumerate(obj.ravel()[:3]):
            walk(f"{name}[{i}]", e, depth + 1)
    else:
        ln = name.lower()
        if any(k in ln for k in KEYS):
            print(f"   {name} = {obj}")

for f in ["Vas_Pam_20mm_Result_20260920_1519.mat", "Bifemsh_20mm_Result.mat", "Vas_Pam_20mm_Result.mat"]:
    p = MO / f
    print(f"== {f}" + ("" if p.exists() else "  MISSING"))
    if not p.exists():
        continue
    r = sio.loadmat(p, squeeze_me=True, struct_as_record=True)
    for k in r:
        if not k.startswith("__"):
            walk(k, r[k])
