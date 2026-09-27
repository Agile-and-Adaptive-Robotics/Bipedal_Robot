# dump_more.py - baselineScores from flxr77 mat + optimized-mat structure
import scipy.io as sio
import numpy as np
import os

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

e = sio.loadmat(os.path.join(BASE, "minimizeExt10mmX3_results_20260920_flxr77.mat"),
                squeeze_me=True, struct_as_record=False)
for k in ['baselineScores', 'baselineScores1', 'baselineScores2', 'grpNames', 'grpTests', 'Xi0', 'Xi1', 'Xi2', 'Xi3']:
    if k in e:
        v = np.asarray(e[k])
        print(f"{k} = {np.array2string(v, precision=6, max_line_width=200)}")
rcvE = np.atleast_1d(e.get('results_cv'))
print(f"results_cv n = {rcvE.size}")
for i, item in enumerate(rcvE.flatten()):
    if hasattr(item, '_fieldnames'):
        fns = item._fieldnames
        out = {}
        for fn in fns:
            try:
                out[fn] = np.asarray(getattr(item, fn)).flatten()[:6]
            except Exception:
                out[fn] = '?'
        print(f" results_cv[{i}]: {out}")

print("\n=== KneeFlx_20mm_optimized.mat ===")
try:
    m = sio.loadmat(os.path.join(BASE, "KneeFlx_20mm_optimized.mat"), squeeze_me=True, struct_as_record=False)
    for k in m.keys():
        if k.startswith('__'):
            continue
        v = m[k]
        if hasattr(v, '_fieldnames'):
            print(f" {k}: struct fields={v._fieldnames[:20]}{'...' if len(v._fieldnames)>20 else ''}")
        else:
            a = np.asarray(v)
            print(f" {k}: {a.dtype} {a.shape}" + (f" val={a.flatten()[:8]}" if a.size <= 12 else ""))
except Exception as ex:
    print(f" loadmat failed: {type(ex).__name__}: {ex}")

print("\n=== KneeExt_20mm_optimized.mat ===")
try:
    m2 = sio.loadmat(os.path.join(BASE, "KneeExt_20mm_optimized.mat"), squeeze_me=True, struct_as_record=False)
    for k in m2.keys():
        if k.startswith('__'):
            continue
        v = m2[k]
        if hasattr(v, '_fieldnames'):
            print(f" {k}: struct fields={v._fieldnames[:20]}{'...' if len(v._fieldnames)>20 else ''}")
        else:
            a = np.asarray(v)
            print(f" {k}: {a.dtype} {a.shape}" + (f" val={a.flatten()[:8]}" if a.size <= 12 else ""))
except Exception as ex:
    print(f" loadmat failed: {type(ex).__name__}: {ex}")

print("\n=== Mesh_Optimization\\Results result mats ===")
RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
for f in ["Bifemsh_20mm_Result.mat", "Vas_Pam_20mm_Result_20260920_1519.mat"]:
    p = os.path.join(RES, f)
    if not os.path.isfile(p):
        print(f" {f}: MISSING")
        continue
    print(f" {f}:")
    try:
        mm = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
        for k in mm.keys():
            if k.startswith('__'):
                continue
            v = mm[k]
            if hasattr(v, '_fieldnames'):
                print(f"   {k}: struct fields={v._fieldnames[:25]}")
            else:
                a = np.asarray(v)
                if a.size <= 8:
                    print(f"   {k}: {a.flatten()}")
                else:
                    print(f"   {k}: {a.dtype} {a.shape}")
    except Exception as ex:
        print(f"   loadmat failed: {type(ex).__name__}: {ex}")
