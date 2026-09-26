import io, sys
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
HERE = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
d = np.load(HERE / "reports_20260925" / "logs" / "syn6_air_smoke.npz", allow_pickle=True)
t = np.asarray(d["t"], float); neuro = np.asarray(d["neuro"], float)
names = [str(x) for x in d["neuro_names"]]
q = np.asarray(d["q"], float)
msk = t >= 5.0
for nm in ("RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l"):
    v = neuro[:, names.index(nm)]
    pk, _ = find_peaks(v[msk], prominence=0.3)
    print(nm, "peaks:", np.round(t[msk][pk], 2).tolist())
ev = []
for nm, tag in (("RG_E_r", "E"), ("RG_F_r", "F")):
    v = neuro[:, names.index(nm)]
    pk, _ = find_peaks(v[msk], prominence=0.3)
    ev += [(x, tag) for x in t[msk][pk]]
ev.sort()
print("interleave:", "".join(s for _, s in ev))
jn = [str(x) for x in d["key_joints"]].index("knee_angle_r")
walk = (t >= 5.0) & (t <= 15.0)
print("knee_angle_r 5-15s:", round(float(q[walk, jn].min()), 1), "..",
      round(float(q[walk, jn].max()), 1), "deg; full-run:",
      round(float(q[:, jn].min()), 1), "..", round(float(q[:, jn].max()), 1))
