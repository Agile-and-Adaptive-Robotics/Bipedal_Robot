# Verify ZCODE_HANDOFF.md syn6-section numbers by re-deriving them (read-only).
import sys
sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
import numpy as np
import build_network_syn6 as b

t = b._load_synergy_tables()
for side in ("r", "l"):
    tab = t[side]
    print(side, "| W shape:", np.shape(tab["W"]),
          "| family:", tab.get("family"),
          "| tau mult:", tab.get("tau"))
sf_r = np.asarray(t["r"]["stance_frac"], dtype=float)
sf_l = np.asarray(t["l"]["stance_frac"], dtype=float)
mean = 0.5 * (sf_r + sf_l)
print("bilateral-mean stance_frac S1..S6:", [round(float(x), 6) for x in mean])
print("handoff claims:                    [0.195713, 0.952928, 0.922159, 0.24633, 0.410065, 0.889709]")
