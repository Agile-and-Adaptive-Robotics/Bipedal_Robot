"""Round 2: list all field names; print rest/tendon-related fields."""
import numpy as np
from scipy.io import loadmat

MAT = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_Result_cam_20261014_0416.mat"
m = loadmat(MAT, squeeze_me=True, struct_as_record=False)

KEYWORDS = ("rest", "tendon", "slack", "l0", "lt", "seed")

def brief(val, maxlen=300):
    a = np.asarray(val)
    if a.dtype.kind in "OSU" or a.dtype.names:
        return f"<nested {a.dtype}>"
    if a.size > 12:
        flat = a.ravel()
        return f"array {a.shape}: first={np.array2string(flat[:6], precision=6)}, last={np.array2string(flat[-3:], precision=6)}"
    return np.array2string(a, precision=6, separator=", ")

def walk(obj, path, depth=0):
    if depth > 3 or not hasattr(obj, "_fieldnames"):
        return
    for f in obj._fieldnames:
        try:
            val = getattr(obj, f)
        except Exception as e:
            print(f"{path}.{f}: <error {e}>")
            continue
        if hasattr(val, "_fieldnames"):
            print(f"{path}.{f}: STRUCT fields={val._fieldnames}")
            walk(val, f"{path}.{f}", depth + 1)
        elif isinstance(val, np.ndarray) and val.dtype.names:
            print(f"{path}.{f}: STRUCT-ARRAY fields={val.dtype.names}")
        else:
            print(f"{path}.{f} = {brief(val)}")

for var in ("predBest", "ctxCamBest", "ctx"):
    print(f"\n================ {var} ================")
    walk(m[var], var)

print("\n================ keyword hits (values) ================")
def kw_print(obj, path, depth=0):
    if depth > 3 or not hasattr(obj, "_fieldnames"):
        return
    for f in obj._fieldnames:
        try:
            val = getattr(obj, f)
        except Exception:
            continue
        low = f.lower()
        if hasattr(val, "_fieldnames"):
            kw_print(val, f"{path}.{f}", depth + 1)
        elif any(k in low for k in KEYWORDS) and not (isinstance(val, np.ndarray) and val.dtype.names):
            try:
                print(f"{path}.{f} = {brief(val)}")
            except Exception:
                pass

for var in ("predBest", "ctxCamBest", "ctx"):
    kw_print(m[var], var)
