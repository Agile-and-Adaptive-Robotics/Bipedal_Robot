# verify_appendix_data.py - check appendix tables vs BPA set mats + extensor release schedule
import scipy.io as sio
import numpy as np
import os

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

print("=== FlxPinBPASet.mat kf(1..5) vs appendix tab:app_flxpin_tests ===")
m = sio.loadmat(os.path.join(BASE, "FlxPinBPASet.mat"), squeeze_me=True, struct_as_record=False)
kf = np.atleast_1d(m['kf'])
print(f" n = {kf.size}")
for i, s in enumerate(kf, 1):
    ak = np.asarray(s.Ak).flatten()
    print(f" kf({i}): rest={float(s.rest):.4f} ten={float(s.ten):.4f} P={float(s.P):.2f} "
          f"Fm={float(s.Fm):.2f} dBPA={float(s.dBPA)*1000:.1f}mm Ak[{ak.min():.1f},{ak.max():.1f}] M={np.asarray(s.Loc).shape[0]}")

print("\n=== ExtPinBPASet.mat ke(1..9) vs appendix tab:app_extpin_tests ===")
e = sio.loadmat(os.path.join(BASE, "ExtPinBPASet.mat"), squeeze_me=True, struct_as_record=False)
ke = np.atleast_1d(e['ke'])
print(f" n = {ke.size}")
for i, s in enumerate(ke, 1):
    ak = np.asarray(s.Ak).flatten()
    print(f" ke({i}): rest={float(s.rest):.4f} ten={float(s.ten):.4f} P={float(s.P):.2f} "
          f"Fm={float(s.Fm):.2f} Ak[{ak.min():.1f},{ak.max():.1f}] M={np.asarray(s.Loc).shape[0]}")

print("\n=== extensor route release schedule from Vas_Pam_20mm_Result_20260925.mat ===")
RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
v = sio.loadmat(os.path.join(RES, "Vas_Pam_20mm_Result_20260925.mat"), squeeze_me=True, struct_as_record=False)
ri = v['routeInfo']
for fn in ['eliminatedPoint', 'eliminatedAngleD', 'eliminatedSweepIndex', 'addedPoint', 'addedAngleD']:
    if hasattr(ri, fn):
        print(f" routeInfo.{fn} = {np.asarray(getattr(ri, fn)).flatten()}")
if hasattr(v, 'transitionIdx'):
    print(f" transitionIdx = {np.asarray(v['transitionIdx']).flatten()}")
