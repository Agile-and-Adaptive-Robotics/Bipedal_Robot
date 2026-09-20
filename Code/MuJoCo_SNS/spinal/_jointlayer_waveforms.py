"""Fit the T1 joint-layer model (one seed) and plot the six fitted
PF half-center waveforms per side on the stance-rescaled gait-phase axis.
Also reports KNEE-F bimodality (the double-knee discriminator: a second
burst in 45-75% beyond the loading burst).
Output: fsa_results/joint_layer_hcs.{png,pdf} + stdout metrics.
"""
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import nnls
from scipy.signal import butter, filtfilt
from sklearn.decomposition import NMF

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import gait_phase as gp
import muscle_map
from draw_circuit import OI_BLUE, OI_VERM, OI_GREEN, OI_ORANGE, OI_PURPLE

HC6 = ("HIP-E", "HIP-F", "KNEE-E", "KNEE-F", "ANKLE-E", "ANKLE-F")
COLORS = {"HIP-E": OI_BLUE, "HIP-F": OI_VERM, "KNEE-E": OI_GREEN,
          "KNEE-F": OI_ORANGE, "ANKLE-E": OI_PURPLE, "ANKLE-F": "#8c6d31"}
GROUP2HC = {"hip_ext": ("HIP-E",), "hip_abd": ("HIP-E",),
            "hip_add": ("HIP-E", "HIP-F"), "hip_flex": ("HIP-F",),
            "knee_ext": ("KNEE-E",), "knee_flex": ("KNEE-F",),
            "ankle_pf": ("ANKLE-E",), "ankle_df": ("ANKLE-F",)}
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

def phase_arr(side):
    hs_s, to_s = gp.contact_events(side)
    ph = np.full(len(t), np.nan)
    for h0, h1 in zip(hs_s[:-1], hs_s[1:]):
        m = (t >= h0) & (t <= h1)
        ph[m] = (t[m] - h0) / (h1 - h0) * 100.0
        for kk in to_s:
            if h0 < kk < h1:
                mm = (t >= h0) & (t <= kk)
                ph[mm] = (t[mm] - h0) / (kk - h0) * 50.0
                mm = (t > kk) & (t <= h1)
                ph[mm] = 50.0 + (t[mm] - kk) / (h1 - kk) * 50.0
    return ph

def solve_w(A, H):
    W = np.zeros((A.shape[0], H.shape[0]))
    for i in range(A.shape[0]):
        W[i], _ = nnls(H.T, A[i])
    return W

fig, axes = plt.subplots(2, 1, figsize=(7.0, 7.2), sharex=True)
summary = {}
for row, side in enumerate(("r", "l")):
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
    ph = phase_arr(side)
    tr = np.arange(0, A.shape[0], 2)
    nmf = NMF(n_components=6, init="nndsvda", max_iter=5000, tol=1e-7,
              random_state=42)
    _Wf = nmf.fit_transform(np.maximum(A, 0.0))
    H = np.maximum(nmf.components_, 1e-6) * M
    W = solve_w(A[tr], H) + 1e-6
    eps = 1e-10
    for _ in range(4000):
        H *= (W.T @ A[tr]) / (W.T @ W @ H + eps)
        H *= M
        W *= (A[tr] @ H.T) / (W @ (H @ H.T) + eps)
    scale = np.maximum(W.max(axis=0), 1e-12)
    Wn = W / scale
    ax = axes[row]
    ordc = np.argsort(ph[tr])
    for ci, hc in enumerate(HC6):
        ax.plot(ph[tr][ordc], Wn[ordc, ci], color=COLORS[hc], lw=1.6,
                label=hc)
    ax.axvline(50.0, color="0.6", ls=(0, (3, 3)), lw=0.9)
    ax.set_xlim(0, 100)
    ax.set_ylabel(f"{side.upper()} leg\nHC amplitude (norm)")
    ax.text(25, 1.13 * Wn.max(), "stance", ha="center", fontsize=7,
            color="0.35")
    ax.text(75, 1.13 * Wn.max(), "swing", ha="center", fontsize=7,
            color="0.35")
    if row == 0:
        ax.legend(fontsize=6.5, ncol=6, loc="upper center",
                  bbox_to_anchor=(0.5, 1.28), frameon=False)
    # knee-F bimodality: primary burst vs 45-75% secondary
    ki = HC6.index("KNEE-F")
    pw = Wn[:, ki]
    m_swing = (ph[tr] >= 45) & (ph[tr] <= 75)
    m_all = np.isfinite(ph[tr])
    prim = float(pw[m_all].max())
    swing_peak = float(pw[m_swing].max()) if m_swing.any() else 0.0
    summary[side] = {"kneef_primary": prim,
                     "kneef_swing_second": swing_peak,
                     "bimodal_ratio": swing_peak / max(prim, 1e-9)}
    print(f"{side}: KNEE-F primary {prim:.2f}, swing-window second "
          f"{swing_peak:.2f}, ratio {summary[side]['bimodal_ratio']:.2f} "
          f"(>0.5 = genuine double-knee drive)")
axes[1].set_xlabel("stance-rescaled gait phase (%)")
fig.suptitle("Fitted joint-layer PF half-center waveforms (T1, SO "
             "activations; ankle timings inherit the known SO artifact)",
             fontsize=8.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = Path("fsa_results/joint_layer_hcs.png")
fig.savefig(out, dpi=200)
fig.savefig(out.with_suffix(".pdf"))
print(f"written: {out}")
