"""Verify the picks of record from the result mats:
- flexor front row 77 of minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat
- XiUsed in Bifemsh_20mm_Result.mat (flexor redesign of record)
- XiUsed in Vas_Pam_20mm_Result_20260920_1519.mat (extensor, flxr77 row 32)
"""
import scipy.io as sio
import numpy as np
from pathlib import Path

TD = Path(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo")
MO = Path(r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results")

def peek(name, obj, depth=0, maxdepth=2):
    pad = "  " * depth
    if depth > maxdepth:
        return
    if isinstance(obj, np.ndarray) and obj.dtype.names:
        print(f"{pad}{name}: struct{obj.shape} fields={obj.dtype.names}")
        for f in obj.dtype.names:
            if obj.size:
                peek(f, obj.flat[0][f] if isinstance(obj.flat[0][f], np.ndarray) else obj.flat[0][f], depth + 1, maxdepth)
    elif isinstance(obj, np.ndarray):
        s = str(obj.shape)
        v = ""
        if obj.size <= 12:
            v = f" = {obj.ravel()[:12]}"
        elif obj.size > 12:
            v = f" first={obj.ravel()[:6]}"
        print(f"{pad}{name}: {obj.dtype} {s}{v}")
    else:
        print(f"{pad}{name}: {type(obj).__name__} = {obj}")

m = sio.loadmat(TD / "minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat", squeeze_me=False, struct_as_record=True)
print("== FLX noT3 mat top keys:", [k for k in m if not k.startswith("__")])
for k in m:
    if not k.startswith("__"):
        peek(k, m[k], 0, 1)

for f in ["Bifemsh_20mm_Result.mat", "Vas_Pam_20mm_Result_20260920_1519.mat", "Vas_Pam_20mm_Result.mat"]:
    p = MO / f
    if not p.exists():
        print(f"== {f}: MISSING at {p}")
        continue
    r = sio.loadmat(p, squeeze_me=True, struct_as_record=True)
    print(f"\n== {f} keys:", [k for k in r if not k.startswith("__")])
    for k in r:
        if k.startswith("__"):
            continue
        v = r[k]
        if hasattr(v, "dtype") and v.dtype.names:
            for fld in v.dtype.names:
                try:
                    x = v[fld].item() if np.ndim(v[fld]) == 0 else v[fld]
                except Exception:
                    x = v[fld]
                print(f"   {k}.{fld} = {x}")
        else:
            print(f"   {k} = {v}")
