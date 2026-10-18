"""Supervisor spot-check: Appendix C adopted-values table vs the two front mats.
Reads flexor pooled front row 77 and extensor flxr77 row 32 + lock metadata."""
import scipy.io as sio
import numpy as np

flx = sio.loadmat(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat", squeeze_me=True, struct_as_record=False)
print("=== FLEXOR MAT keys ===")
print([k for k in flx.keys() if not k.startswith("__")])
for k in flx.keys():
    if k.startswith("__"):
        continue
    v = flx[k]
    if isinstance(v, np.ndarray):
        print(k, "shape", v.shape, "dtype", v.dtype)
    else:
        print(k, type(v), v if np.isscalar(v) or isinstance(v, str) else "")

ext = sio.loadmat(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\minimizeExt10mmX3_results_20260920_flxr77.mat", squeeze_me=True, struct_as_record=False)
print("=== EXTENSOR MAT keys ===")
print([k for k in ext.keys() if not k.startswith("__")])
for k in ext.keys():
    if k.startswith("__"):
        continue
    v = ext[k]
    if isinstance(v, np.ndarray):
        print(k, "shape", v.shape, "dtype", v.dtype)
    else:
        print(k, type(v), v if np.isscalar(v) or isinstance(v, str) else "")
