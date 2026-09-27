"""SUPERVISOR GATE: independent scipy re-derivation of the Xi numbers of
record (read-only). Dumps:
  - flexor front rows 77 + 107 (xCols [4 5 6] 1-based) from
    minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat
  - extensor front rows 32 + 1 (xCols [5 6 7 8] 1-based) from
    minimizeExt10mmX3_results_20260920_flxr77.mat + LOCKSRC/flxRow/xi1lock/
    xi2lock/BOUNDS fields
  - Bifemsh_20mm_Result.mat XiUsed + min torque margin
Expected (claims): g77 = [3.93736985e-03, 3.99787285e+04, 1.47337523e+04];
g107 = [8.85741e-3, 56237.8, 18542.3]; g32 = [-6.378022e-3, 39978.7,
14733.8, 0.158289]; ext row1 = [-0.0101197, 43535.7, 17014.3, 0.620932].
"""
import io
import os
import sys

import numpy as np
from scipy.io import loadmat

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"
RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"


def row(mat, name, r, cols0):
    arr = mat[name]
    print(f"{name} shape={arr.shape}")
    return np.asarray(arr[r - 1, cols0], dtype=float).ravel()


m1 = loadmat(os.path.join(BASE, "minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat"),
             squeeze_me=True, struct_as_record=False)
print("flexor mat keys:", [k for k in m1 if not k.startswith("__")])
g77 = row(m1, "filtered_results", 77, [3, 4, 5])
g107 = row(m1, "filtered_results", 107, [3, 4, 5])
print("g77  =", [f"{v:.10g}" for v in g77])
print("g107 =", [f"{v:.10g}" for v in g107])

m2 = loadmat(os.path.join(BASE, "minimizeExt10mmX3_results_20260920_flxr77.mat"),
             squeeze_me=True, struct_as_record=False)
print("ext mat keys:", [k for k in m2 if not k.startswith("__")])
g32 = row(m2, "filtered_results", 32, [4, 5, 6, 7])
g1 = row(m2, "filtered_results", 1, [4, 5, 6, 7])
print("g32  =", [f"{v:.10g}" for v in g32])
print("g1   =", [f"{v:.10g}" for v in g1])
for fld in ["LOCKSRC", "flxRow", "xi1lock", "xi2lock", "BOUNDS", "KMODE"]:
    if fld in m2:
        print(f"{fld} = {m2[fld]!r}")

m3 = loadmat(os.path.join(RES, "Bifemsh_20mm_Result.mat"),
             squeeze_me=True, struct_as_record=False)
print("bifemsh keys:", [k for k in m3 if not k.startswith("__")])
if "XiUsed" in m3:
    print("XiUsed =", np.asarray(m3["XiUsed"], dtype=float).ravel())
for fld in ["minMargin", "minTorqueMarginFraction", "stamp", "liveRun",
            "fBest"]:
    if fld in m3:
        print(f"{fld} = {m3[fld]!r}")
if "torqueMarginFraction" in m3:
    tm = np.asarray(m3["torqueMarginFraction"], dtype=float).ravel()
    print(f"min torqueMarginFraction = {tm.min():.6g} at idx {int(tm.argmin())}")
if "phiD" in m3:
    phi = np.asarray(m3["phiD"], dtype=float).ravel()
    tm = np.asarray(m3["torqueMarginFraction"], dtype=float).ravel()
    print(f"at phiD = {phi[int(tm.argmin())]:.3f}")
