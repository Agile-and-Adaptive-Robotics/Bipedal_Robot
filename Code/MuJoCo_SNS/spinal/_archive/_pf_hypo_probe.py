"""PF-layer design probe (Ben's 2026-09-18 hypotheses, read-only).

(1) Where does each key muscle peak within the measured gait cycle?
(2) Does rect_fem/sartorius timing support driving swing-phase knee
    re-extension from a HIP-layer cell (biarticular convergence)?
Phase = stance-rescaled percent: heel strike 0, toe-off 50, next HS 100
(same normalization as plot_gait_joint_angles.py / gait_phase.py).
"""
import io
import json
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

d = np.load("bsolve_out.npz", allow_pickle=True)
acts = d["acts"]                       # (frames, 92) converted-model SO
names = [str(s) for s in d["act_names"]]
t = np.asarray(d["t"], float)
print(f"acts {acts.shape}, t {t[0]:.3f}..{t[-1]:.3f}s")

# --- measured events (same source as the joint-angle figure) ----------
import gait_phase as gp

side = "r"
hs, to = gp.contact_events(side)
print(f"events: HS {np.round(hs, 3)}  TO {np.round(to, 3)}")

# stance-rescaled phase per frame (right leg)
ph = np.full(len(t), np.nan)
cycles = list(zip(hs[:-1], hs[1:])) if len(hs) > 1 else []
for h0, h1 in cycles:
    m = (t >= h0) & (t <= h1)
    ph[m] = (t[m] - h0) / (h1 - h0) * 100.0
    for k in to:
        if h0 < k < h1:
            mm = (t >= h0) & (t <= k)
            ph[mm] = (t[mm] - h0) / (k - h0) * 50.0
            mm = (t > k) & (t <= h1)
            ph[mm] = 50.0 + (t[mm] - k) / (h1 - k) * 50.0
ok = np.isfinite(ph)
print(f"frames with phase: {ok.sum()}/{len(t)} "
      f"({ph[ok].min():.1f}..{ph[ok].max():.1f}%)")

# --- S-component muscle weights from the FSA backsolve ----------------
try:
    fb = json.load(open("fsa_results/fsa_backsolve.json"))
    comp = fb["sides"]["r"]["components"]
except Exception:
    comp = None

def in_comp(ci, muscle):
    if comp is None:
        return "?"
    c = comp[ci] if isinstance(comp, list) else comp[str(ci)]
    wdict = c.get("muscle_weights", c.get("weights", {}))
    for nm, w in wdict.items():
        if muscle.replace("_r", "") in nm and float(w) > 0.02:
            return f"{float(w):.2f}"
    return "."

WIND = [("loading flexion 5-20%", 5, 20),
        ("late-stance ext 30-48%", 30, 48),
        ("swing flexion 50-65%", 50, 65),
        ("swing re-ext 65-100%", 65, 100)]
MUS = ["rect_fem", "sart", "vas_lat", "vas_med", "semimem", "bifemsh",
       "grac", "glut_max", "glut_med", "iliacus", "psoas",
       "tib_ant", "per_tert", "sol", "gas_med", "gas_lat", "tib_post"]
print("\nmuscle          peak%  |  " + " | ".join(f"{w[0]:<21s}" for w in WIND)
      + " | S1..S6 weights (RF/sart row = Ben's hip-layer test)")
hdr = None
for mus in MUS:
    idx = [i for i, n in enumerate(names) if n.startswith(mus) and
           n.endswith("_r")]
    for i in idx:
        a = acts[:, i]
        if a.max() < 0.03:
            continue
        a_ph = a[ok]
        p_ph = ph[ok]
        pk = p_ph[int(np.argmax(a_ph))]
        cells = []
        for _, lo, hi in WIND:
            m = (p_ph >= lo) & (p_ph < hi)
            cells.append(f"{a_ph[m].mean():.3f}")
        wrow = ""
        if comp is not None:
            try:
                wrow = " ".join(in_comp(ci, mus) for ci in range(6))
            except Exception:
                wrow = "?"
        print(f"{names[i]:<15s} {pk:5.1f}  |  " +
              " | ".join(f"{c:<21s}" for c in cells) + f" | {wrow}")

# --- RF/sartorius vs swing knee re-extension --------------------------
knee_i = [i for i, n in enumerate(names) if n == "rect_fem_r"][0]
print("\nRF (rect_fem_r) activation, per 5% phase bin (right leg):")
bins = np.arange(0, 101, 5)
a = acts[:, knee_i][ok]
p = ph[ok]
for lo in bins[:-1]:
    m = (p >= lo) & (p < lo + 5)
    if m.any():
        bar = "#" * int(round(a[m].mean() * 60))
        print(f"  {lo:3.0f}-{lo+5:3.0f}%  {a[m].mean():.3f} {bar}")
