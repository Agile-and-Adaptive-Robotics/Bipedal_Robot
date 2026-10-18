# Read-only key dump of the per-test archive mats to look for stored GoF
# fields (RMSE/FVU) of any model vintage.
import scipy.io as sio
import os
BASE = r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo"
for f in ["KneeFlxPin_10mm_48cm.mat", "Plot_KneeFlxPin10mm_48cm.mat",
          "KneeExtPin_10mm_all.mat", "Plot_KneeExtPin10mm_40cm.mat"]:
    p = os.path.join(BASE, f)
    try:
        m = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
        keys = sorted(k for k in m if not k.startswith("__"))
        print(f, "->", keys)
    except Exception as e:
        print(f, "ERR", e)
