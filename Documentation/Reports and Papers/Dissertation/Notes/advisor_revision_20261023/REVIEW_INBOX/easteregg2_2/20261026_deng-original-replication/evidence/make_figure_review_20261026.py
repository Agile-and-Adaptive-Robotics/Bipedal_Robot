r"""Review-batch figure for the Deng-original replication proposal.

Dissertation conventions: Arial; (A)-(D) as MARGIN letters left of each tile
straddling the plot-box top (FlxPin_group convention); Tol-bright blue
(baseline) / red (perturbed); vector PDF + PNG preview. Data: the saved
deng_rg_mn.txt chart files (unchanged from the campaign analysis).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..",
                   "Documentation", "Reports and Papers", "Dissertation",
                   "Notes", "advisor_revision_20261023", "REVIEW_INBOX",
                   "easteregg2_2", "20261026_deng-original-replication",
                   "figures"))
os.makedirs(OUT, exist_ok=True)

COLS = ["LH_RG_E", "LH_RG_F", "RH_RG_E", "RH_RG_F",
        "LH_KneeMN_E", "LH_KneeMN_F", "RH_KneeMN_E", "RH_KneeMN_F",
        "LH_KAPF_E", "LH_KAPF_F", "LH_HipPF_E"]

def load(run, col):
    raw = np.genfromtxt(os.path.join(HERE, run, "deng_rg_mn.txt"),
                        delimiter="\t", skip_header=1)
    nonzero = np.where(np.abs(np.nan_to_num(raw)).sum(axis=1) > 0)[0]
    raw = raw[:nonzero[-1] + 1]
    raw = raw[~np.isnan(raw).any(axis=1)]
    return raw[:, 1], raw[:, 2 + COLS.index(col)] * 1000.0

fig, axs = plt.subplots(4, 1, figsize=(5.2, 6.8), sharex=True)
C0, C1 = "#4477AA", "#EE6677"
LETTERS = ["(A)", "(B)", "(C)", "(D)"]

panels = [
    ("run_B0_baseline", "LH_RG_F", "run_P1b_rgflexor_2p1", "LH_RG_F",
     "left rhythm-generator flexor, \\textbf{2 nA at the RG flexor} (2.0--2.1 s)"),
    ("run_B0_baseline", "LH_RG_E", "run_P1b_rgflexor_2p1", "LH_RG_E",
     "left rhythm-generator extensor, same stimulus"),
    ("run_B0_baseline", "LH_RG_E", "run_P2b_kneepfflexor_2p1", "LH_RG_E",
     "left rhythm-generator extensor, \\textbf{2 nA at the knee PF flexor} (2.0--2.1 s)"),
    ("run_B0_baseline", "LH_KneeMN_E", "run_P2_kneepfflexor", "LH_KneeMN_E",
     "left knee extensor motoneuron, knee PF flexor held 2.0--2.5 s"),
]
for ax, (r0, c0, r1, c1, desc), letter in zip(axs, panels, LETTERS):
    t0, v0 = load(r0, c0)
    t1, v1 = load(r1, c1)
    ax.plot(t0, v0, color=C0, lw=0.7, label="baseline")
    ax.plot(t1, v1, color=C1, lw=0.7, alpha=0.85, label="stimulated")
    ax.axvspan(2.0, 2.5, color="0.90", zorder=0)
    ax.axvline(2.1, color="0.45", lw=0.6, ls=(0, (2, 2)))
    ax.set_ylim(-105, 5)
    ax.set_yticks([-100, -75, -50, -25, 0])
    ax.text(-0.155, 0.97, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top", ha="left")
    ax.text(0.995, 0.05, desc, transform=ax.transAxes, fontsize=7,
            va="bottom", ha="right", color="0.25")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
axs[0].legend(loc="lower right", fontsize=6.5, frameon=False, ncol=2)
for ax in axs:
    ax.set_ylabel("voltage (mV)", fontsize=8)
axs[-1].set_xlabel("time (s)", fontsize=8)
axs[-1].set_xlim(1.6, 3.6)
fig.subplots_adjust(left=0.17, right=0.98, top=0.97, bottom=0.07, hspace=0.32)
base = os.path.join(OUT, "deng_original_perturbation_20261026")
fig.savefig(base + ".pdf")
fig.savefig(base + ".png", dpi=300)
print("wrote", base + ".pdf/.png")
