"""Extract torque spans from flexor redesign + pulley mats (READ-ONLY)."""
import scipy.io as sio
import numpy as np

MESH = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

for f, lab in [
    (r"\Bifemsh_20mm_Result_20260921_2037.mat", "FLX DIRECT 0921"),
    (r"\Bifemsh_20mm_Result_pulley_20260926_1652.mat", "FLX PULLEY 0926"),
]:
    print(f"\n===== {lab} =====")
    m = sio.loadmat(MESH + f, squeeze_me=True, struct_as_record=False)
    for k in sorted([k for k in m if not k.startswith("__")]):
        v = m[k]
        if isinstance(v, np.ndarray) and v.dtype.kind != "O":
            if v.ndim >= 1 and v.size > 10:
                # candidate torque/angle arrays
                lk = k.lower()
                if any(t in lk for t in ("torque", "human", "phi", "angle", "tau")):
                    a = np.atleast_1d(v.astype(float)).ravel()
                    a = a[np.isfinite(a)]
                    if a.size:
                        print(f"  {k}: n={a.size} min={a.min():.4f} max={a.max():.4f} mean={a.mean():.4f}")
    # struct fields for ctx-like objects
    for k in m:
        if k.startswith("__"):
            continue
        v = m[k]
        if hasattr(v, "_fieldnames"):
            print(f"  struct {k}: fields {v._fieldnames}")
