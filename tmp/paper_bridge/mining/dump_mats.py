"""Dump key result mats for the numbers ledger (READ-ONLY)."""
import scipy.io as sio
import numpy as np

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"
MESH = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"


def show(name, path, fields=None):
    print(f"\n===== {name}: {path} =====")
    try:
        m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    except Exception as e:
        print("  LOAD FAILED:", e)
        return None
    keys = [k for k in m if not k.startswith("__")]
    print("  top keys:", keys)
    return m


def row_dump(m, matvar, row):
    """Dump one front row from a results mat."""
    for k in m:
        if k.startswith("__"):
            continue
        v = m[k]
        if isinstance(v, np.ndarray) and v.dtype == object:
            continue
        try:
            if hasattr(v, "shape") and v.ndim >= 2 and v.shape[0] > row:
                print(f"  {k}[{row},:] = {np.array2string(np.atleast_1d(v[row]), precision=6)}")
            elif hasattr(v, "shape") and v.ndim == 1 and v.shape[0] > row:
                print(f"  {k}[{row}] = {v[row]!r}")
        except Exception as e:
            print(f"  {k}: ({e})")


# 1. Flexor two-bracket front of record, row 77
m = show("FLX 2brkt 2trans noT3 (row 77)", BASE + r"\minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat")
if m is not None:
    for cand in ("Xi", "X", "front", "pareto", "sol", "results"):
        hits = [k for k in m if not k.startswith("__") and cand.lower() in k.lower()]
        if hits:
            print("  candidate-key hits:", hits)
    row_dump(m, None, 77)

# 2. Extensor locked refit, row 32
m = show("EXT flxr77 (row 32)", BASE + r"\minimizeExt10mmX3_results_20260920_flxr77.mat")
if m is not None:
    row_dump(m, None, 32)
    for k in m:
        if not k.startswith("__"):
            v = m[k]
            if np.isscalar(v) or (hasattr(v, "size") and v.size <= 4):
                print(f"  {k} = {v!r}")

# 3. Extensor redesign of record
m = show("VAS 20260925", MESH + r"\Vas_Pam_20mm_Result_20260925.mat")
if m is not None:
    for k in m:
        if k.startswith("__"):
            continue
        v = m[k]
        if np.isscalar(v):
            print(f"  {k} = {v!r}")
        elif hasattr(v, "shape") and v.size <= 8:
            print(f"  {k} = {np.array2string(np.atleast_1d(v), precision=5)}")

# 4. Cam redesign
m = show("VAS cam 20261014", MESH + r"\Vas_Pam_20mm_Result_cam_20261014_0416.mat")
if m is not None:
    for k in m:
        if k.startswith("__"):
            continue
        v = m[k]
        if np.isscalar(v):
            print(f"  {k} = {v!r}")
        elif hasattr(v, "shape") and v.size <= 10:
            print(f"  {k} = {np.array2string(np.atleast_1d(v), precision=5)}")

# 5. Pulley flexor
m = show("BIFEMSH pulley 20260926", MESH + r"\Bifemsh_20mm_Result_pulley_20260926_1652.mat")
if m is not None:
    for k in m:
        if k.startswith("__"):
            continue
        v = m[k]
        if np.isscalar(v):
            print(f"  {k} = {v!r}")
        elif hasattr(v, "shape") and v.size <= 10:
            print(f"  {k} = {np.array2string(np.atleast_1d(v), precision=5)}")

print("\nDONE")
