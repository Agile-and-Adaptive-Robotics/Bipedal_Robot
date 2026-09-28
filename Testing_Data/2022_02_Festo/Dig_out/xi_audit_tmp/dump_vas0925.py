# dump_vas0925.py - full inspection of Vas_Pam_20mm_Result_20260925 + optimized plot mats
import scipy.io as sio
import numpy as np
import os, datetime

RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

def mt(p):
    return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M")

p = os.path.join(RES, "Vas_Pam_20mm_Result_20260925.mat")
print(f"--- Vas_Pam_20mm_Result_20260925.mat (mtime {mt(p)}) ---")
m = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
ks = sorted(k for k in m.keys() if not k.startswith('__'))
print("keys:", ks)
for k in ks:
    v = m[k]
    if hasattr(v, '_fieldnames'):
        print(f" {k}: struct fields={v._fieldnames[:30]}")
    else:
        a = np.asarray(v)
        if a.size <= 10:
            print(f" {k} = {a.flatten()}")
        else:
            print(f" {k}: {a.dtype} {a.shape}")

# margins if predBest-like structures exist under any name
for cand in ['predBest', 'pred', 'predX0']:
    if cand in m and hasattr(m[cand], '_fieldnames'):
        pb = m[cand]
        print(f" {cand}.Torque shape:", np.asarray(pb.Torque).shape if hasattr(pb, 'Torque') else 'n/a')

print("\n--- optimized plot mats: mtimes + curve ratio checks ---")
for f, lbl in [("KneeFlx_20mm_optimized.mat", "FLEXOR"), ("KneeExt_20mm_optimized.mat", "EXTENSOR")]:
    p = os.path.join(BASE, f)
    m2 = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
    print(f"\n {f} (mtime {mt(p)}):")
    ang = np.asarray(m2['Angle']).flatten()
    print(f"  Angle: [{ang.min():.2f}, {ang.max():.2f}] deg, N={ang.size}")
    # flexor plot: Bifemsh_T (robot?) vs H (human); extensor: Torque0/2/3 vs HumanLocation?
    for k in m2.keys():
        if k.startswith('__') or k == 'Angle':
            continue
        v = np.asarray(m2[k])
        if v.dtype.kind in 'fiu' and v.size == ang.size:
            print(f"  {k}: [{np.nanmin(v):.4f}, {np.nanmax(v):.4f}]")
