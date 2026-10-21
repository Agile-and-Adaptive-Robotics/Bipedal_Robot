"""Round 3: dump predBest.bpa cells + ctxCamBest/ctx BPA-related scalars."""
import numpy as np
from scipy.io import loadmat

MAT = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_Result_cam_20261014_0416.mat"
m = loadmat(MAT, squeeze_me=True, struct_as_record=False)
pred = m["predBest"]

def unwrap(v):
    """scipy cell-array representation -> python list."""
    if isinstance(v, np.ndarray) and v.dtype.names:
        return [unwrap(v[n]) for n in v.dtype.names]
    return v

bpas = unwrap(pred.bpa)
print(f"bpa container type: {type(bpas)}, len={len(bpas) if hasattr(bpas,'__len__') else 'n/a'}")
if isinstance(bpas, list):
    for i, b in enumerate(bpas):
        print(f"\n-- bpa[{i}] --")
        if hasattr(b, "_fieldnames"):
            for f in b._fieldnames:
                val = getattr(b, f)
                if hasattr(val, "_fieldnames"):
                    print(f"   {f}: STRUCT {val._fieldnames}")
                else:
                    a = np.asarray(val) if isinstance(val, np.ndarray) else val
                    if isinstance(a, np.ndarray) and a.size > 8:
                        print(f"   {f} = array {a.shape}")
                    else:
                        print(f"   {f} = {a}")
        else:
            print(f"   raw: {b}")

for var in ("ctxCamBest", "ctx"):
    obj = m[var]
    print(f"\n== {var} scalar/keyword fields ==")
    for f in obj._fieldnames:
        low = f.lower()
        if any(k in low for k in ("bpa", "count", "pressure", "kpa")):
            val = getattr(obj, f)
            if hasattr(val, "_fieldnames"):
                print(f"   {f}: STRUCT {val._fieldnames}")
            elif isinstance(val, np.ndarray) and val.dtype.names:
                print(f"   {f}: CELL/STRUCT-ARRAY {val.dtype.names}")
            else:
                print(f"   {f} = {val}")

print("\nrest      =", pred.rest, "m")
print("tendon    =", pred.tendon, "m")
print("restLmt   =", pred.restLmt, "m")
print("cRestLength =", pred.cRestLength, "m")
print("BPAcount  =", pred.BPAcount)
print("maxPathLength0 =", pred.maxPathLength0, "m (at", pred.maxPathAngleD0, "deg)")
print("Lmt0 range:", pred.Lmt0.min(), "to", pred.Lmt0.max())
print("Lmt_p range:", pred.Lmt_p.min(), "to", pred.Lmt_p.max())
