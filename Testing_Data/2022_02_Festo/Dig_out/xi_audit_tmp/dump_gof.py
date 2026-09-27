# dump_gof.py - extract provenance + GoF fields from the two front mats
import scipy.io as sio
import numpy as np
import os

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

def s(x):
    if hasattr(x, 'dtype') and x.dtype.kind in 'US':
        return str(x)
    return None

def scalarish(v):
    try:
        a = np.asarray(v).flatten()
        return a
    except Exception:
        return v

print("=== FLEXOR FRONT minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat ===")
m = sio.loadmat(os.path.join(BASE, "minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat"),
                squeeze_me=True, struct_as_record=False)
for k in ['TAG', 'TRANSMODE', 'SOLVER', 'PICK', 'MAXGEN', 'POP', 'NUMHOLD', 'USE_BRACKET2', 'W', 'a0', 'f']:
    if k in m:
        v = m[k]
        print(f" {k} = {v!r}" if not hasattr(v, '_fieldnames') else f" {k}: fields={v._fieldnames}")
rcv = m.get('results_cv')
print(f" results_cv: type={type(rcv).__name__}", end="")
if hasattr(rcv, '_fieldnames'):
    print(f" fields={rcv._fieldnames}")
    for fn in rcv._fieldnames:
        val = getattr(rcv, fn)
        print(f"   .{fn} -> {np.asarray(val).shape if not hasattr(val,'_fieldnames') else 'struct'}")
else:
    print(f" value={rcv!r}")
fr = np.asarray(m['filtered_results'])
print(f" filtered_results full rows (all 13 cols):")
for row in (1, 77, 107):
    print(f"  row {row}: {np.array2string(fr[row-1,:], precision=6, max_line_width=200)}")
# results_cv struct array contents (per-fold GoF of save-time pick)
if hasattr(rcv, '_fieldnames') or (hasattr(rcv, 'dtype') and rcv.dtype == np.object_):
    arr = np.atleast_1d(rcv)
    for i, item in enumerate(arr):
        if hasattr(item, '_fieldnames'):
            print(f" results_cv[{i}]: " + ", ".join(f"{fn}={np.asarray(getattr(item,fn)).flatten()[:6]}" for fn in item._fieldnames))

print()
print("=== EXTENSOR FRONT minimizeExt10mmX3_results_20260920_flxr77.mat ===")
e = sio.loadmat(os.path.join(BASE, "minimizeExt10mmX3_results_20260920_flxr77.mat"),
                squeeze_me=True, struct_as_record=False)
for k in ['RESULTFILE', 'LOCKSRC', 'KMODE', 'BOUNDS', 'tag', 'flxRow', 'xi1lock', 'xi2lock',
          'pick', 'fold', 'numHold', 'exitflag']:
    if k in e:
        print(f" {k} = {e[k]!r}")
frE = np.asarray(e['filtered_results'])
print(f" filtered_results rows (all 15 cols):")
for row in (1, 32):
    print(f"  row {row}: {np.array2string(frE[row-1,:], precision=6, max_line_width=250)}")
for k in ['f_all', 'scores_cv', 'baselineScores', 'fvals', 'f']:
    if k in e:
        v = np.asarray(e[k])
        print(f" {k}: shape={v.shape}")
        if v.size <= 30:
            print(f"   values: {np.array2string(v, precision=6, max_line_width=200)}")
rcvE = e.get('results_cv')
arrE = np.atleast_1d(rcvE)
print(f" results_cv: n={arrE.size}")
for i, item in enumerate(arrE.flatten()):
    if hasattr(item, '_fieldnames'):
        vals = {}
        for fn in item._fieldnames:
            try:
                a = np.asarray(getattr(item, fn)).flatten()
                vals[fn] = a[:8]
            except Exception:
                vals[fn] = '?'
        print(f" results_cv[{i}]: {vals}")
