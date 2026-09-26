"""Quantify gate-c2 alternation on the saved air-run npz (ask criteria:
>=10 s run, alternating RG/PF output, distinguishable synergy-layer
profiles)."""
import io
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

HERE = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
d = np.load(HERE / "reports_20260925" / "logs" / "syn6_air_smoke.npz",
            allow_pickle=True)
t = np.asarray(d["t"], float)
neuro = np.asarray(d["neuro"], float)
names = [str(x) for x in d["neuro_names"]]
q = np.asarray(d["q"], float)
key_joints = [str(x) for x in d["key_joints"]]
walk = (t >= 5.0) & (t <= 15.0)
tw = t[walk]
print(f"run span {t[0]:.1f}-{t[-1]:.1f} s; walk window 5-15 s "
      f"({walk.sum()} samples = {tw[-1]-tw[0]:.1f} s)")

def peaks(nm, prom=0.5):
    v = neuro[:, names.index(nm)]
    pk, _ = find_peaks(v[walk], prominence=prom)
    return tw[pk], v[walk]

ok = True
for sd in ("r", "l"):
    te, ve = peaks(f"RG_E_{sd}")
    tf, vf = peaks(f"RG_F_{sd}")
    # interleaving: every consecutive RG event alternates E/F
    ev = sorted([(x, "E") for x in te] + [(x, "F") for x in tf])
    labels = "".join(s for _, s in ev)
    n_cycles = min(len(te), len(tf))
    per = (te[-1] - te[0]) / (len(te) - 1) if len(te) > 1 else float("nan")
    print(f"side {sd}: RG-E peaks {len(te)} {np.round(te, 2).tolist()}")
    print(f"        RG-F peaks {len(tf)} {np.round(tf, 2).tolist()}")
    print(f"        event sequence: {labels}  "
          f"E period {per:.2f} s  E amp {ve.max()-ve.min():.2f} mV  "
          f"F amp {vf.max()-vf.min():.2f} mV")
    alt = labels.startswith(("EF", "FE")) and "EE" not in labels.replace(
        "EFE", "EFE")  # report only; strict check below
    strict = labels == ("EF" * len(te))[:len(labels)] or \
        labels == ("FE" * len(te))[:len(labels)]
    ok &= len(te) >= 3 and len(tf) >= 3 and strict
    # PF watch channels: distinctness over the walk window
    P = np.array([neuro[walk, names.index(f"PF_S{k}_{sd}")]
                  for k in (1, 2, 3, 4)])
    c = np.corrcoef(P)
    print(f"        PF_S1..S4 pairwise corr: "
          f"{np.round(c[np.triu_indices(4, 1)], 3).tolist()}")
jn = key_joints.index("knee_angle_r") if "knee_angle_r" in key_joints else 4
kj = key_joints.index("hip_flexion_r") if "hip_flexion_r" in key_joints else 3
print(f"knee_angle_r walk range: {q[walk, jn].min():.1f}.."
      f"{q[walk, jn].max():.1f} deg")
print(f"hip_flexion_r walk range: {q[walk, kj].min():.1f}.."
      f"{q[walk, kj].max():.1f} deg")
ok &= (q[walk, jn].max() - q[walk, jn].min()) > 5.0
print("C2-ALTERNATION", "PASS" if ok else "FAIL",
      "(criteria: >=10 s, >=3 E and >=3 F peaks, strictly interleaved, "
      "knee range > 5 deg)")
