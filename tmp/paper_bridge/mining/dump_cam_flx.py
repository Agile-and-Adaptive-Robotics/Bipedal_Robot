"""Dump cam mat + flexor direct-route mats + routeInfo fields (READ-ONLY)."""
import scipy.io as sio
import numpy as np

MESH = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"


def scalars(path, label, maxel=12):
    print(f"\n===== {label} =====")
    try:
        m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    except Exception as e:
        print("  LOAD FAILED:", e)
        return
    for k in sorted([k for k in m if not k.startswith("__")]):
        v = m[k]
        if np.isscalar(v) or isinstance(v, str):
            print(f"  {k} = {v!r}")
        elif hasattr(v, "shape") and v.size <= maxel and v.dtype.kind != "O":
            print(f"  {k} = {np.array2string(np.atleast_1d(v), precision=5)}")


scalars(MESH + r"\Vas_Pam_20mm_Result_cam_20261014_0416.mat", "CAM EXT")

# routeInfo of the extensor-of-record + cam (nested struct)
for f, lab in [
    (r"\Vas_Pam_20mm_Result_20260925.mat", "VAS0925 routeInfo"),
    (r"\Vas_Pam_20mm_Result_cam_20261014_0416.mat", "CAM routeInfo"),
]:
    print(f"\n===== {lab} =====")
    m = sio.loadmat(MESH + f, squeeze_me=True, struct_as_record=False)
    ri = m.get("routeInfo")
    if ri is None:
        print("  no routeInfo")
        continue
    names = ri._fieldnames
    print("  fields:", names)
    for n in names:
        v = getattr(ri, n)
        if isinstance(v, np.ndarray) and v.dtype.kind != "O" and v.size <= 16:
            print(f"  {n} = {np.array2string(np.atleast_1d(v), precision=5)}")
        elif np.isscalar(v) or isinstance(v, str):
            print(f"  {n} = {v!r}")

# flexor direct-route candidates (newest first)
for f in [
    r"\Bifemsh_20mm_Result_20260921_2037.mat",
    r"\Bifemsh_20mm_Result_20260920_1443.mat",
]:
    scalars(MESH + f, f, maxel=10)
