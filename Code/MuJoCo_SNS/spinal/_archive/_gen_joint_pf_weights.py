"""Generate joint_pf_weights.json from the T1 joint-layer fit.

One masked fit per side (seed-42 warm start, same protocol as
fsa_jointlayers.py T1); H rows normalized to max 1.0 per HC; merged
across sides by max.  Muscles that never entered the fit (inactive in
the SO, e.g. sartorius) fall back to the anatomical group mapping at
1.0 so the wiring stays complete for all 46 muscles.
"""
import io
import json
import sys

import numpy as np
from scipy.optimize import nnls
from scipy.signal import butter, filtfilt
from sklearn.decomposition import NMF

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import gait_phase as gp  # noqa: F401  (phase not needed here, keep parity)
import muscle_map

HC6 = ("HIP-E", "HIP-F", "KNEE-E", "KNEE-F", "ANK-E", "ANK-F")
GROUP2HC = {"hip_ext": ("HIP-E",), "hip_abd": ("HIP-E",),
            "hip_add": ("HIP-E", "HIP-F"), "hip_flex": ("HIP-F",),
            "knee_ext": ("KNEE-E",), "knee_flex": ("KNEE-F",),
            "ankle_pf": ("ANK-E",), "ankle_df": ("ANK-F",)}
CROSS = {"rect_fem": ("HIP-F",), "semimem": ("HIP-E",),
         "semiten": ("HIP-E",), "bifemsh": ("HIP-E",),
         "gas_med": ("KNEE-F",), "gas_lat": ("KNEE-F",),
         "grac": ("HIP-F",), "sart": ("KNEE-F",)}

d = np.load("bsolve_out.npz", allow_pickle=True)
acts = np.asarray(d["acts"], float)
names = [str(s) for s in d["act_names"]]
fmax = np.asarray(d["Fmax"], float)
t = np.asarray(d["t"], float)
dt = float(np.median(np.diff(t)))
bb, aa = butter(4, 6.0 / (0.5 / dt))
SM = np.clip(filtfilt(bb, aa, acts, axis=0), 0.0, 1.0)


def solve_w(A, H):
    W = np.zeros((A.shape[0], H.shape[0]))
    for i in range(A.shape[0]):
        W[i], _ = nnls(H.T, A[i])
    return W


W = {hc: {} for hc in HC6}
for side in ("r", "l"):
    keep, lay = [], []
    for i, nm in enumerate(names):
        if not nm.endswith(f"_{side}") or fmax[i] <= 5.0 \
                or acts[:, i].max() <= 0.05:
            continue
        base = nm[:-2]
        g = muscle_map._GROUPS_BY_NAME.get(base)
        hcs = list(GROUP2HC.get(g[0], ())) if g else []
        if not hcs:
            continue
        for sub, extra in CROSS.items():
            if base.startswith(sub):
                hcs += [x for x in extra if x not in hcs]
        keep.append(i)
        lay.append(hcs)
    A = SM[:, keep]
    M = np.zeros((6, len(keep)))
    for j, allowed in enumerate(lay):
        for L in allowed:
            M[HC6.index(L), j] = 1.0
    tr = np.arange(0, A.shape[0], 2)
    nmf = NMF(n_components=6, init="nndsvda", max_iter=5000, tol=1e-7,
              random_state=42)
    nmf.fit(np.maximum(A, 0.0))
    H = np.maximum(nmf.components_, 1e-6) * M
    Wm = solve_w(A[tr], H) + 1e-6
    eps = 1e-10
    for _ in range(4000):
        H *= (Wm.T @ A[tr]) / (Wm.T @ Wm @ H + eps)
        H *= M
        Wm *= (A[tr] @ H.T) / (Wm @ (H @ H.T) + eps)
    Hn = H / np.maximum(H.max(axis=1, keepdims=True), 1e-12)
    for j, gi in enumerate(keep):
        base = names[gi][:-2]
        for ci, hc in enumerate(HC6):
            if M[ci, j] > 0 and Hn[ci, j] > 0.05:
                W[hc][base] = max(W[hc].get(base, 0.0), float(Hn[ci, j]))

# fallback: every model muscle gets its anatomical mapping at 1.0 unless
# the fit supplied a weight
for base, g in muscle_map._GROUPS_BY_NAME.items():
    for hc in GROUP2HC.get(g[0], ()):
        for hc_key in (hc,):
            pass
    for hc in GROUP2HC.get(g[0], ()):
        if base not in W[hc]:
            W[hc][base] = 1.0
    for sub, extra in CROSS.items():
        if base.startswith(sub):
            for hc in extra:
                if base not in W[hc]:
                    W[hc][base] = 1.0

out = {"source": "T1 joint-layer fit (fsa_jointlayers protocol, seed 42, "
                  "max over sides; rows max-normalized)",
       "hcs": list(HC6), "w": W}
open("joint_pf_weights.json", "w", encoding="utf-8").write(
    json.dumps(out, indent=1))
n = sum(len(v) for v in W.values())
print(f"joint_pf_weights.json: {n} (hc,muscle) entries")
for hc in HC6:
    top = sorted(W[hc].items(), key=lambda kv: -kv[1])[:4]
    print(f"  {hc}: {len(W[hc])} muscles; top: "
          + ", ".join(f"{k}={v:.2f}" for k, v in top))
