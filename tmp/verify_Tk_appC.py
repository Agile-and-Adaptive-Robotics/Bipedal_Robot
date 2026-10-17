"""Verify 94-AppendixC.tex prose claims against FlxPinBPASet.mat (kf(1)):
- per-frame transform T_k has rotation block R_z(theta_k) sweeping -120..+10 deg
- translation column constant [4.5, -395.8, 0] mm
Also dump ExtPinBPASet.mat ke(1) transform summary for contrast."""
import scipy.io as sio
import numpy as np

BASE = r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo"

def walk(obj, path, depth=0, out=None):
    if out is None:
        out = []
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        # struct array
        if obj.size:
            names = getattr(obj.flat[0], 'dtype', None)
        flat = obj.ravel()
        if flat.size and isinstance(flat[0], np.ndarray) and flat[0].dtype.names:
            pass
    return out

def describe(name, val, maxitems=6):
    if isinstance(val, np.ndarray) and val.dtype.names:
        print(f"{name}: struct fields={val.dtype.names}")
    elif isinstance(val, np.ndarray):
        print(f"{name}: ndarray shape={val.shape} dtype={val.dtype}")
        if val.size and val.size <= 24:
            print(f"   value={np.asarray(val).ravel()}")
    else:
        r = repr(val)
        print(f"{name}: {r[:120]}")

for fname, key in (("FlxPinBPASet.mat", "kf"), ("ExtPinBPASet.mat", "ke")):
    print("=" * 70)
    print(fname)
    d = sio.loadmat(rf"{BASE}\{fname}", squeeze_me=False, struct_as_record=True)
    top = [k for k in d if not k.startswith("__")]
    print("top-level keys:", top)
    if key in d:
        arr = d[key]
        print(f"{key}: shape={arr.shape}")
        s0 = arr.flat[0]
        print(f"{key}(1) fields:", s0.dtype.names)
        for nm in s0.dtype.names:
            v = s0[nm]
            describe(f"  {key}(1).{nm}", v)
