"""Make the three missing SCONE figures for the 2026-09-25 walker campaign.

Read-only on logs/ and the goal1 report (numbers transcribed from its tables).
Outputs PNGs into ../  (reports_20260925/figs/). Regenerate with:
  "C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe" make_scone_figs.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from scone_sto import col, read_sto  # noqa: E402

LOGS = os.path.normpath(os.path.join(HERE, "..", "..", "logs"))
FIGS = os.path.normpath(os.path.join(HERE, ".."))

# Okabe-Ito CVD-safe palette (house style)
C_INIT = "#0072B2"   # blue  - pretrained / default
C_D4 = "#D55E00"     # vermillion - short-horizon (overfit) best
C_D20 = "#009E73"    # bluish green - full-horizon best
C_40G = "#E69F00"    # orange - 40-gen balance
C_MAX = "#D55E00"
C_MIN = "#0072B2"


def load(fname):
    names, data = read_sto(os.path.join(LOGS, fname))
    t = col(names, data, "time")
    return t, col(names, data, "pelvis.pos.x"), col(names, data, "pelvis.pos.y")


# ----------------------------------------------------------------------------
# Figure 1 - gait overfit: forward progression + pelvis height, 3 pars
# ----------------------------------------------------------------------------
def fig_gait_overfit():
    t_i, x_i, y_i = load("motion_reeval_gait4a_init_d20.sto")
    t_s, x_s, y_s = load("motion_reeval_gait4a_best_d20.sto")
    t_f, x_f, y_f = load("motion_reeval_gait4a_bestD20_d20.sto")

    t_fall = t_s[-1]  # 4.68 s (report: 4.685)
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True, constrained_layout=True
    )

    ax1.plot(t_i, x_i, color=C_INIT, lw=1.6,
             label="pretrained par (20-s fitness 0.883, 36 steps)")
    ax1.plot(t_s, x_s, color=C_D4, lw=1.6,
             label="4-s-window CMA-ES best (fitness 0.858 @ 4 s)")
    ax1.plot(t_f, x_f, color=C_D20, lw=1.6,
             label="20-s full-horizon CMA-ES best (fitness 0.873, 36 steps)")
    for ax in (ax1, ax2):
        ax.axvline(t_fall, color=C_D4, ls="--", lw=1.0, alpha=0.6)

    ax1.annotate(
        "short-horizon best FALLS at t = %.2f s\n(5.05 m; 20-s fitness 83.3 = fell)"
        % t_fall,
        xy=(t_fall, x_s[-1]), xytext=(6.8, 2.6),
        arrowprops=dict(arrowstyle="->", color=C_D4), color=C_D4, fontsize=10,
    )
    ax1.annotate("21.48 m", xy=(t_i[-1], x_i[-1]), xytext=(-2, 10),
                 textcoords="offset points", ha="right", color=C_INIT)
    ax1.annotate("21.08 m", xy=(t_f[-1], x_f[-1]), xytext=(-2, -16),
                 textcoords="offset points", ha="right", color=C_D20)
    ax1.set_ylabel("pelvis forward distance (m)")
    ax1.set_title("SCONE Tutorial 4a Gait (Hyfydy): short-window optimization overfits")
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, 20.6)
    ax1.set_ylim(0, 23)

    ax2.plot(t_i, y_i, color=C_INIT, lw=1.6)
    ax2.plot(t_s, y_s, color=C_D4, lw=1.6)
    ax2.plot(t_f, y_f, color=C_D20, lw=1.6)
    ax2.annotate("fall: pelvis height collapses", xy=(t_fall, y_s[-1]),
                 xytext=(6.2, 0.45), arrowprops=dict(arrowstyle="->", color=C_D4),
                 color=C_D4, fontsize=10)
    ax2.set_xlabel("time (s)")
    ax2.set_ylabel("pelvis height (m)")
    ax2.set_ylim(0.3, 1.05)
    ax2.grid(alpha=0.3)

    out = os.path.join(FIGS, "scone_gait_overfit.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


# ----------------------------------------------------------------------------
# Figure 2 - balance: pelvis height, init vs 40-gen vs 300-gen
# ----------------------------------------------------------------------------
def fig_balance_gains():
    t0, _, y0 = load("motion_reeval_bal_init_d30.sto")
    t4, _, y4 = load("motion_reeval_bal40_best_d30.sto")
    t3, _, y3 = load("motion_reeval_bal300_best_d30.sto")

    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    ax.plot(t0, y0, color=C_INIT, lw=1.8,
            label="default/init par (fitness 101.8 @ 30 s)")
    ax.plot(t4, y4, color=C_40G, lw=1.8,
            label="CMA-ES 40-gen best (fitness 94.5 @ 30 s)")
    ax.plot(t3, y3, color=C_D20, lw=1.8,
            label="CMA-ES 300-gen best (fitness 2.19 @ 30 s)")

    ax.annotate("falls at t = %.2f s" % t0[-1], xy=(t0[-1], y0[-1]),
                xytext=(2.2, 0.62), arrowprops=dict(arrowstyle="->", color=C_INIT),
                color=C_INIT, fontsize=10)
    ax.annotate("falls at t = %.2f s" % t4[-1], xy=(t4[-1], y4[-1]),
                xytext=(4.6, 0.50), arrowprops=dict(arrowstyle="->", color=C_40G),
                color=C_40G, fontsize=10)
    ax.annotate("stands the full 30 s\n(pelvis %.2f m)" % y3[-1],
                xy=(t3[-1], y3[-1]), xytext=(21.5, 0.80),
                arrowprops=dict(arrowstyle="->", color=C_D20),
                color=C_D20, fontsize=10)

    ax.set_xlabel("time (s)")
    ax.set_ylabel("pelvis height (m)")
    ax.set_xlim(0, 31)
    ax.set_ylim(0.3, 1.05)
    ax.set_title("SCONE Tutorial 3a Balance (Hyfydy): generations buy standing balance")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(alpha=0.3)

    out = os.path.join(FIGS, "scone_balance_gains.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


# ----------------------------------------------------------------------------
# Figure 3 - tutorial score bars (15 never-run tutorials, default params)
# Numbers transcribed verbatim from goal1_scone_tutorials_and_optimization.md §1.
# ----------------------------------------------------------------------------
TUT = [
    # (label, result, style)  style: 'max' = maximize objective, 'min' = minimize
    ("1   Introduction (JumpMeasure)", 42.3804, "max"),
    ("2a  High Jump", 42.3804, "max"),
    ("2b  High Jump Polynomial", 42.3804, "max"),
    ("2c  Straight Pose Jump", 13.7741, "max"),
    ("3b  Motor Noise Balance", 102.629, "min"),
    ("4b  Fast Gait", 89.3969, "min"),
    ("4c  Perturbed Gait", 91.3812, "min"),
    ("4d  Slippery Slope", 97.7414, "min"),
    ("5a  Plantarflexor Weakness", 79.4282, "min"),
    ("5b  Short Hamstrings", 96.933, "min"),
    ("5c  Hyper-reflexia", 91.6098, "min"),
    ("6a  Script Body Height", 0.024403, "max"),
    ("6b  Script Gyro Balance", 95.8342, "min"),
    ("6c  Script Reflex Modulation", 87.3623, "min"),
    ("6d  Script Neural Delays", 0.00946241, "max"),
]


def fig_tutorial_scores():
    labels = [r[0] for r in TUT]
    vals = np.array([r[1] for r in TUT])
    styles = [r[2] for r in TUT]
    colors = [C_MAX if s == "max" else C_MIN for s in styles]

    ypos = np.arange(len(TUT))
    fig, ax = plt.subplots(figsize=(10, 7.8))
    fig.subplots_adjust(top=0.89, bottom=0.17, left=0.26, right=0.985)
    ax.barh(ypos, vals, color=colors, edgecolor="black", lw=0.4, height=0.68)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontfamily="monospace")
    ax.invert_yaxis()  # tutorial 1 on top
    ax.set_xlabel("objective score at default parameters (as printed by sconecmd -l 2)")
    ax.set_xlim(0, 118)
    for y, v in zip(ypos, vals):
        ax.text(v + 1.2, y, "%.4g" % v, va="center", fontsize=8.5, color="0.25")

    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(fc=C_MIN, label="minimize-style (gait/balance; lower = better)"),
            Patch(fc=C_MAX, label="maximize-style (Jump/ScriptMeasure; higher = better)"),
        ],
        loc="lower right", fontsize=9,
    )
    ax.set_title(
        "SCONE tutorials 1, 2a-2c, 3b, 4b-4d, 5a-5c, 6a-6d (Hyfydy, default params,\n"
        "first evaluation of each on EB475WS4, 2026-09-25)"
    )
    fig.text(
        0.01, 0.045,
        "Jump-style tutorials (1, 2a-2c) and ScriptMeasure scenarios (6a, 6d) use MAXIMIZE objectives; "
        "all gait/balance rows use MINIMIZE.\nScores are not comparable across objective styles. "
        "Source: goal1_scone_tutorials_and_optimization.md table, §1.",
        fontsize=8, color="0.35", va="bottom",
    )
    ax.grid(axis="x", alpha=0.3)

    out = os.path.join(FIGS, "scone_tutorial_scores.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    fig_gait_overfit()
    fig_balance_gains()
    fig_tutorial_scores()
