# inspect_mats.py - Xi/GoF audit (easteregg2 myo env, scipy.io.loadmat)
# Read-only inspection of the front mats + optimized mats + result mats.
import scipy.io as sio
import numpy as np
import os, sys

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"
RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def show_keys(path):
    try:
        m = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
        print(f"  keys: {[k for k in m.keys() if not k.startswith('__')]}")
        return m
    except NotImplementedError as e:
        print(f"  V7.3/HDF5 MAT (loadmat cannot read): {e}")
        return None
    except Exception as e:
        print(f"  ERROR: {type(e).__name__}: {e}")
        return None

def dump_front(path, picks):
    print(f"\n=== {os.path.basename(path)} ===")
    m = show_keys(path)
    if m is None:
        return None
    fr = m.get('filtered_results')
    xc = m.get('xCols')
    if fr is None:
        print("  no filtered_results")
        return None
    fr = np.asarray(fr); xc = np.asarray(xc).flatten().astype(int)
    print(f"  filtered_results shape = {fr.shape}, xCols = {xc.tolist()}")
    out = {}
    for label, row in picks:
        g = fr[row-1, :].flatten()
        xs = ", ".join(f"{v:.6g}" for v in g[xc-1])
        print(f"  row {row} ({label}): xCols values [{xs}]")
        out[label] = g
    # look for GoF-ish sibling variables
    for k in m.keys():
        if k.startswith('__'):
            continue
        if k in ('filtered_results', 'xCols'):
            continue
        v = m[k]
        try:
            print(f"  var {k}: type={type(v).__name__} shape={getattr(v,'shape',None)}")
        except Exception:
            pass
    return m, out

print("### 1. flexor front (pick-77 source)")
fm = dump_front(os.path.join(BASE, "minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat"),
                [("pick77", 77), ("row107", 107), ("row1", 1)])

print("\n### 2. extensor flxr77 front (pick-32 source)")
em = dump_front(os.path.join(BASE, "minimizeExt10mmX3_results_20260920_flxr77.mat"),
                [("pick32", 32)])
