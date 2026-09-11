"""Figure suite for spinal runs - Ben's requested plots (2026-09-11):

  fig1  hip/knee/ankle OVER the neural stack (RG, PF, MN layers)
  fig2  joint movement over stimulus (DRIVE)
  fig3  RG/PF/MN activity over stimulus (DRIVE)
  fig4  joint-motion cycle overlays (stand->swing), cycle-normalized
  fig5  our mean gait cycle vs the OpenSim IK benchmark (if available)

Usage (myo env, from Code/MuJoCo_SNS/spinal):
    python plot_run.py [run.npz] [--side r]

Reads spinal_run.npz by default (the newest run) and writes fig1..fig5 PNGs
next to it. The benchmark panel needs `subject01_walk1_ik.mot` (or any
*_ik.mot) under Solid_Models\\OpenSim\\ - the converted model preserved
OpenSim coordinate conventions, so the curves overlay directly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
REPO = HERE.parents[2]

C_JOINTS = {"hip": "#0072B2", "knee": "#D55E00", "ankle": "#009E73"}
C_STANCE = "0.85"

JOINT_KEYS = {"hip": "hip_flexion_{s}", "knee": "knee_angle_{s}",
              "ankle": "ankle_angle_{s}"}
MN_GROUPS = {
    "extensors": ("vas_lat_{s}", "soleus_{s}", "med_gas_{s}", "glut_max2_{s}"),
    "flexors": ("semimem_{s}", "tib_ant_{s}", "psoas_{s}", "sar_{s}"),
}


def load(path: str, side: str):
    d = np.load(path, allow_pickle=True)
    t = d["t"]
    q = d["q"]                      # degrees, saved by runner
    act = d["act"]
    neuro = d["neuro"]
    names = list(d["key_joints"])
    neuro_names = (list(d["neuro_names"]) if "neuro_names" in d
                   else ["DRIVE", "POSTURE", "RG_E_r", "RG_F_r",
                         "PF_E2_r", "PF_F1_r", "RG_E_l", "RG_F_l"])
    acts = list(d["key_acts"])
    cols = {jn: i for i, jn in enumerate(names)}
    acols = {a: i for i, a in enumerate(acts)}
    ncol = {nm: i for i, nm in enumerate(neuro_names)}
    j = {k: q[:, cols[v.format(s=side)]] for k, v in JOINT_KEYS.items()}
    drive = neuro[:, ncol["DRIVE"]]
    return d, t, j, act, acols, neuro, ncol, drive, side


def stance_spans(neuro, ncol, t, side="r"):
    """(t0, t1) spans where RG_E is high -> shade as stance."""
    rge = neuro[:, ncol[f"RG_E_{side}"]]
    thr = 0.5 * np.max(rge)
    on = rge > thr
    edges = np.diff(on.astype(int))
    starts = list(np.flatnonzero(edges == 1) + 1)
    ends = list(np.flatnonzero(edges == -1) + 1)
    if on[0]:
        starts = [0] + starts
    if on[-1]:
        ends = ends + [len(on)]
    return [(t[a], t[b - 1]) for a, b in zip(starts, ends) if b - 1 > a]


def shade_stance(ax, spans):
    for a, b in spans:
        ax.axvspan(a, b, color=C_STANCE, zorder=0, lw=0)


# ---------------------------------------------------------------- fig 1
def fig1(t, j, neuro, ncol, act, acols, side, spans):
    fig, ax = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    fig.subplots_adjust(hspace=0.14, left=0.08, right=0.97, top=0.95)
    for a in ax:
        shade_stance(a, spans)
    # joints
    for k, v in j.items():
        ax[0].plot(t, v, lw=1.2, color=C_JOINTS[k], label=k)
    ax[0].set_ylabel(f"{side} joint angle [deg]")
    ax[0].legend(fontsize=8, ncol=3, loc="upper right")
    ax[0].set_title("joints over the neural stack (shaded = RG-E stance)",
                    fontsize=11)
    # RG
    for nm, c in ((f"RG_E_{side}", "#D55E00"), (f"RG_F_{side}", "#0072B2"),
                  (f"RG_E_{'l' if side == 'r' else 'r'}", "#E69F00"),
                  (f"RG_F_{'l' if side == 'r' else 'r'}", "#56B4E9")):
        if nm in ncol:
            ax[1].plot(t, neuro[:, ncol[nm]], lw=1.0, color=c,
                       label=nm.replace("_", " "))
    ax[1].set_ylabel("RG [mV]")
    ax[1].legend(fontsize=8, ncol=4, loc="upper right")
    # PF
    for nm, c in (("PF_E1", "#D55E00"), ("PF_E2", "#E69F00"),
                  ("PF_F1", "#0072B2"), ("PF_F2", "#56B4E9")):
        k = f"{nm}_{side}"
        if k in ncol:
            ax[2].plot(t, neuro[:, ncol[k]], lw=1.0, color=c, label=nm)
    ax[2].set_ylabel("PF [mV]")
    ax[2].legend(fontsize=8, ncol=4, loc="upper right")
    # MN (activation)
    for a, c in (("vas_lat", "#D55E00"), ("soleus", "#E69F00"),
                 ("med_gas", "#8B4513"), ("semimem", "#0072B2"),
                 ("tib_ant", "#56B4E9"), ("psoas", "#009E73"),
                 ("glut_max2", "#CC79A7"), ("rect_fem", "#000000")):
        k = f"{a}_{side}"
        if k in acols:
            ax[3].plot(t, act[:, acols[k]], lw=1.0, color=c, label=a)
    ax[3].set_ylabel("MN activation")
    ax[3].set_xlabel("t [s]")
    ax[3].legend(fontsize=7, ncol=4, loc="upper right")
    fig.savefig(HERE / "fig1_joints_over_neural.png", dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------- fig 2/3
def fig23(t, j, neuro, ncol, act, acols, drive, side, spans):
    # fig2: joints over stimulus
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    fig.subplots_adjust(hspace=0.16, left=0.08, right=0.93, top=0.95)
    for a in ax:
        shade_stance(a, spans)
        at = a.twinx()
        at.plot(t, drive, lw=1.4, color="0.4", ls="--")
        at.set_ylabel("DRIVE [nA]", color="0.4")
        at.tick_params(axis="y", colors="0.4")
        at.set_ylim(0, max(1e-9, drive.max()) * 1.6)
    for i, k in enumerate(("hip", "knee", "ankle")):
        ax[i].plot(t, j[k], lw=1.2, color=C_JOINTS[k])
        ax[i].set_ylabel(f"{k} [deg]")
    ax[0].set_title(f"{side} joint movement over stimulus (grey dashed "
                    "= DRIVE; shaded = stance)", fontsize=11)
    ax[-1].set_xlabel("t [s]")
    fig.savefig(HERE / "fig2_joints_over_stimulus.png", dpi=140)
    plt.close(fig)

    # fig3: neural activity over stimulus
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    fig.subplots_adjust(hspace=0.16, left=0.08, right=0.93, top=0.95)
    for a in ax:
        shade_stance(a, spans)
        at = a.twinx()
        at.plot(t, drive, lw=1.4, color="0.4", ls="--")
        at.set_ylabel("DRIVE [nA]", color="0.4")
        at.tick_params(axis="y", colors="0.4")
        at.set_ylim(0, max(1e-9, drive.max()) * 1.6)
    for nm, c in ((f"RG_E_{side}", "#D55E00"), (f"RG_F_{side}", "#0072B2")):
        ax[0].plot(t, neuro[:, ncol[nm]], lw=1.1, color=c, label=nm)
    ax[0].set_ylabel("RG [mV]")
    ax[0].legend(fontsize=8, loc="upper right")
    for nm, c in (("PF_E1", "#D55E00"), ("PF_E2", "#E69F00"),
                  ("PF_F1", "#0072B2"), ("PF_F2", "#56B4E9")):
        k = f"{nm}_{side}"
        if k in ncol:
            ax[1].plot(t, neuro[:, ncol[k]], lw=1.0, color=c, label=nm)
    ax[1].set_ylabel("PF [mV]")
    ax[1].legend(fontsize=8, ncol=4, loc="upper right")
    ext = [f"MN {g}" for g in MN_GROUPS]
    for g, c in (("extensors", "#D55E00"), ("flexors", "#0072B2")):
        members = [a.format(s=side) for a in MN_GROUPS[g]]
        idx = [acols[m] for m in members if m in acols]
        ax[2].plot(t, act[:, idx].mean(axis=1), lw=1.2, color=c,
                   label=f"{g} (mean MN)")
    ax[2].set_ylabel("MN activation")
    ax[2].set_xlabel("t [s]")
    ax[2].legend(fontsize=8, loc="upper right")
    ax[0].set_title(f"neural activity over stimulus ({side} leg)", fontsize=11)
    fig.savefig(HERE / "fig3_neural_over_stimulus.png", dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------- fig 4/5
def gait_cycles(t, neuro, ncol, side="r"):
    """Index ranges of full gait cycles: RG_E rise -> next RG_E rise,
    within the walk window (DRIVE > 0.8*max)."""
    drive = neuro[:, ncol["DRIVE"]]
    walk = drive > 0.8 * drive.max()
    rge = neuro[:, ncol[f"RG_E_{side}"]]
    thr = 0.5 * np.max(rge)
    on = (rge > thr) & walk
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1) + 1
    if len(rises) < 3:
        return []
    return [(rises[i], rises[i + 1]) for i in range(len(rises) - 1)]


def fig45(t, j, neuro, ncol, side, bench):
    cycles = gait_cycles(t, neuro, ncol, side)
    fig, ax = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    fig.subplots_adjust(hspace=0.13, left=0.1, right=0.97, top=0.95)
    mean_cyc = {}
    for i, k in enumerate(("hip", "knee", "ankle")):
        traces = []
        for a, b in cycles:
            x = j[k][a:b]
            if len(x) < 10:
                continue
            ph = np.linspace(0, 100, len(x))
            traces.append(np.interp(np.linspace(0, 100, 200), ph, x))
            ax[i].plot(ph, x, lw=0.5, color=C_JOINTS[k], alpha=0.35)
        if traces:
            mean_cyc[k] = np.mean(traces, axis=0)
            ax[i].plot(np.linspace(0, 100, 200), mean_cyc[k], lw=2.2,
                       color=C_JOINTS[k], label=f"sim mean ({len(traces)} cyc)")
            if k in bench:
                ax[i].plot(bench[k][0], bench[k][1], lw=1.8, color="k",
                           ls="--", label="OpenSim IK")
                ax[i].set_ylim(min(ax[i].get_ylim()[0], bench[k][1].min() - 5),
                               max(ax[i].get_ylim()[1], bench[k][1].max() + 5))
        ax[i].set_ylabel(f"{k} [deg]")
        ax[i].legend(fontsize=8, loc="upper right")
        ax[i].axvspan(0, 60, color=C_STANCE, zorder=0, lw=0)
    ax[0].set_title(f"gait cycles, cycle-normalized ({side} leg; 0 = RG-E "
                    "rise ~ stance start)", fontsize=11)
    ax[-1].set_xlabel("gait cycle [%]")
    fig.savefig(HERE / "fig4_gait_cycles.png", dpi=140)
    plt.close(fig)
    if not bench:
        print("(fig4: no OpenSim benchmark file found - sim-only overlay. "
              "Drop subject01_walk1_ik.mot under Solid_Models\\OpenSim\\ "
              "and rerun for fig5 comparison.)")
    return mean_cyc


def load_benchmark(side):
    """Find an OpenSim IK .mot and extract hip/knee/ankle cycle curves.
    Returns {} if absent."""
    cands = list((REPO / "Solid_Models" / "OpenSim").rglob("*_ik.mot"))
    if not cands:
        return {}
    f = cands[0]
    print(f"benchmark: {f.name}")
    rows, headers, names = [], None, []
    with open(f, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if line == "endheader":
                headers = True
                continue
            if headers and not names:
                names = line.split()
                continue
            if names:
                try:
                    rows.append([float(x) for x in line.split()])
                except ValueError:
                    break
    data = np.array(rows)
    col = {n: i for i, n in enumerate(names)}
    out = {}
    for k, v in JOINT_KEYS.items():
        cn = v.format(s=side)
        if cn in col:
            out[k] = (None, np.degrees(data[:, col[cn]]))
    # cycle-normalize by time column (assume steady walking, full file =
    # integer cycles; simple fixed split at equal chunks)
    if "time" in col:
        tt = data[:, col["time"]]
        for k in out:
            y = out[k][1]
            n = len(y)
            ph = np.linspace(0, 100, n)
            out[k] = (np.linspace(0, 100, 200),
                      np.interp(np.linspace(0, 100, 200), ph, y))
    return out


def main(argv):
    path = "spinal_run.npz"
    side = "r"
    args = list(argv)
    while args:
        a = args.pop(0)
        if a == "--side":
            side = args.pop(0)
        elif not a.startswith("--"):
            path = a
    d, t, j, act, acols, neuro, ncol, drive, side = load(path, side)
    spans = stance_spans(neuro, ncol, t, side)
    fig1(t, j, neuro, ncol, act, acols, side, spans)
    fig23(t, j, neuro, ncol, act, acols, drive, side, spans)
    bench = load_benchmark(side)
    fig45(t, j, neuro, ncol, side, bench)
    print("saved fig1_joints_over_neural.png, fig2_joints_over_stimulus.png,")
    print("      fig3_neural_over_stimulus.png, fig4_gait_cycles.png")


if __name__ == "__main__":
    main(sys.argv[1:])
