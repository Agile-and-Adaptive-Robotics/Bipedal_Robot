import numpy as np
import scipy.io as sio
import re

R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

# 1. OpenSim target txt column spans
for f, col in [(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\OpenSim_Vasti_Results.txt", 5),
               (r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\OpenSim_Bifem_Results.txt", 4)]:
    H = np.loadtxt(f, skiprows=7)
    print(f.split("\\")[-1], "shape", H.shape,
          "col%d span %.4f..%.4f" % (col, H[:, col-1].min(), H[:, col-1].max()),
          "col2 (angle) %.2f..%.2f n=%d" % (H[:, 1].min(), H[:, 1].max(), H.shape[0]))

# 2. Flexor + pulley ctx constants
for f in ["Bifemsh_20mm_Result.mat", "Bifemsh_20mm_Result_pulley_20260926_1652.mat"]:
    d = sio.loadmat(R + "\\" + f, squeeze_me=True, struct_as_record=False)
    c = d["ctx"]
    print(f, "Dia=%s BPAcount=%s P=%s KMAX=%s reqMargin=%s" %
          (c.Dia, c.BPAcount, c.targetPressure, c.KMAX, c.requiredTorqueMargin))

# 3. muscle_map pool counts (primary-group partition)
src = open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\muscle_map.py").read()
body = src.split("_GROUPS_BY_NAME")[1].split("}\n")[0]
entries = re.findall(r'"([a-z0-9_]+)":\s*\(([^)]*)\)', body)
pools = {}
for name, groups in entries:
    gs = re.findall(r'"([a-z_]+)"', groups)
    if gs:
        pools.setdefault(gs[0], []).append(name)
total = sum(len(v) for v in pools.values())
print("muscle_map primary-group pools:")
for g in ["hip_ext", "hip_flex", "knee_ext", "knee_flex", "ankle_pf", "ankle_df",
          "hip_abd", "hip_add", "trunk_ext", "trunk_flex"]:
    print("  %-10s %2d  %s" % (g, len(pools.get(g, [])), sorted(pools.get(g, []))))
print("TOTAL per side =", total)
