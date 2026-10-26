r"""Evidence figure for Order B1: original Deng 2019 response to the Deng 2023
perturbation protocol. Four panels from the saved chart files (mV, 1 ms).

A: LH RG flexor - baseline vs paper-exact P1b (rapid short step at 2.0 s)
B: LH RG extensor - same runs (phase reset, rhythm continues)
C: LH RG extensor - baseline vs paper-exact P2b (knee-PF kick: no RG effect)
D: LH knee extensor MN - baseline vs long-window P2 (deep local suppression
   only appears when the PF neuron is held for 0.5 s)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
COLS = ["LH_RG_E", "LH_RG_F", "RH_RG_E", "RH_RG_F",
        "LH_KneeMN_E", "LH_KneeMN_F", "RH_KneeMN_E", "RH_KneeMN_F",
        "LH_KAPF_E", "LH_KAPF_F", "LH_HipPF_E"]

def load(run, col):
    raw = np.genfromtxt(os.path.join(HERE, run, "deng_rg_mn.txt"),
                        delimiter="\t", skip_header=1)
    nonzero = np.where(np.abs(np.nan_to_num(raw)).sum(axis=1) > 0)[0]
    raw = raw[:nonzero[-1] + 1]
    raw = raw[~np.isnan(raw).any(axis=1)]
    t = raw[:, 1]
    v = raw[:, 2 + COLS.index(col)] * 1000.0
    return t, v

def style(ax):
    ax.set_ylabel("membrane voltage (mV)")
    for run, colr, lab in []:
        pass

fig, axs = plt.subplots(4, 1, figsize=(7.5, 9.0), sharex=True)
C0, C1 = "#4477AA", "#EE6677"   # Tol bright: blue = baseline, red = perturbed
T = (1.6, 3.6)

panels = [
    ("A  LH RG flexor - 2 nA, 2.0-2.1 s (paper-exact)", "run_B0_baseline", "LH_RG_F",
     "run_P1b_rgflexor_2p1", "LH_RG_F", "P1b: rapid short step at 2.0 s"),
    ("B  LH RG extensor - same runs", "run_B0_baseline", "LH_RG_E",
     "run_P1b_rgflexor_2p1", "LH_RG_E", "P1b: phase reset, rhythm continues"),
    ("C  LH RG extensor - knee-PF 2 nA, 2.0-2.1 s", "run_B0_baseline", "LH_RG_E",
     "run_P2b_kneepfflexor_2p1", "LH_RG_E", "P2b: no RG effect (2019 family)"),
    ("D  LH knee extensor MN - knee-PF 2 nA, 2.0-2.5 s", "run_B0_baseline", "LH_KneeMN_E",
     "run_P2_kneepfflexor", "LH_KneeMN_E", "P2 long window: local suppression"),
]
for ax, (title, r0, c0, r1, c1, note) in zip(axs, panels):
    t0, v0 = load(r0, c0)
    t1, v1 = load(r1, c1)
    ax.plot(t0, v0, color=C0, lw=0.8, label="baseline")
    ax.plot(t1, v1, color=C1, lw=0.8, alpha=0.85, label=note)
    ax.axvspan(2.0, 2.5, color="0.85", zorder=0)
    ax.axvline(2.1, color="0.55", lw=0.7, ls=":")
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_ylim(-105, 5)
    ax.legend(loc="lower right", fontsize=7, frameon=False)
axs[-1].set_xlabel("time (s)")
axs[-1].set_xlim(*T)
for ax in axs:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.tight_layout()
out = os.path.join(HERE, "deng_original_perturbation_20261026.png")
fig.savefig(out, dpi=180)
print("wrote", out)
