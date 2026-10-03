import scipy.io as sio
import numpy as np

RES = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
FILES = [
    "Bifemsh_20mm_Result.mat",
    "Vas_Pam_20mm_Result.mat",
    "Vas_Pam_20mm_Result_20260920_1519.mat",
    "Bifemsh_20mm_Result_pulley_20260926_1652.mat",
    "BiPulley_opensim_bifemsh_r_20260926_1702.mat",
]

def describe(v, depth=0, name=""):
    pad = "  " * depth
    if isinstance(v, np.ndarray) and v.dtype.names:
        print(f"{pad}{name}: struct{v.shape} fields={v.dtype.names}")
        if depth < 2:
            for f in v.dtype.names[:25]:
                sub = v.flat[0][f] if v.size else None
                describe(sub, depth + 1, f)
        return
    if isinstance(v, np.ndarray):
        flat = v.ravel()
        if v.size <= 12:
            print(f"{pad}{name}: array{v.shape} dtype={v.dtype} = {flat.tolist()}")
        else:
            print(f"{pad}{name}: array{v.shape} dtype={v.dtype} min={np.nanmin(flat):.6g} max={np.nanmax(flat):.6g}")
        return
    if isinstance(v, np.generic) or isinstance(v, (int, float, str, bool)):
        print(f"{pad}{name}: {type(v).__name__} = {v}")
        return
    print(f"{pad}{name}: {type(v)} = {v}")

for fn in FILES:
    print("=" * 80)
    print("FILE:", fn)
    try:
        m = sio.loadmat(f"{RES}\\{fn}", squeeze_me=False, struct_as_record=True)
    except Exception as e:
        print("  LOAD ERROR:", e)
        continue
    for k, v in m.items():
        if k.startswith("__"):
            continue
        describe(v, 0, k)
