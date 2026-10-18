# Read-only dump of the archived result mats that Table 4.2/4.3 numbers must
# trace to: flexor front mat a0/baseline + save-time pick GoF, extensor flxr77
# baselineScores + per-fold records. scipy only; nothing written.
import scipy.io as sio
import numpy as np
import os

BASE = r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo"

def load(name):
    return sio.loadmat(os.path.join(BASE, name), squeeze_me=True, struct_as_record=False)

def show(label, v):
    a = np.asarray(v)
    print(f"{label}: shape={a.shape}")
    print(np.array2string(a, precision=6, max_line_width=240))

print("=== flexor front: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat ===")
F = load("minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat")
print("keys:", sorted(k for k in F if not k.startswith('__')))
for k in ("labels", "a0", "f", "k1", "k2", "k3", "ALLBPA"):
    if k in F:
        show(k, F[k])
fr = np.asarray(F["filtered_results"])
print("filtered_results shape:", fr.shape, "xCols:", np.asarray(F["xCols"]).flatten() if "xCols" in F else "n/a")
print("row77:", np.array2string(fr[76, :], precision=6, max_line_width=240))

print()
print("=== extensor front: minimizeExt10mmX3_results_20260920_flxr77.mat ===")
E = load("minimizeExt10mmX3_results_20260920_flxr77.mat")
print("keys:", sorted(k for k in E if not k.startswith('__')))
for k in ("baselineScores", "labels", "ALLBPA", "xi1lock", "xi2lock", "flxRow", "LOCKSRC"):
    if k in E:
        show(k, E[k])
frE = np.asarray(E["filtered_results"])
print("filtered_results shape:", frE.shape)
print("row32:", np.array2string(frE[31, :], precision=6, max_line_width=240))
