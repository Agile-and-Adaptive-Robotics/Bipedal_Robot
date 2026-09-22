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
    contact = d["contact"] if "contact" in d.files else None
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
    d.close()   # NpzFile keeps the handle open (WinError-5 trap)
    return d, t, q, j, act, acols, neuro, ncol, drive, side, contact


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
    """(kept for reference/back-compat) RG_E-rise cycles. fig45 no
    longer uses this - the v1 neural phasing is what let the one-legged
    s3b 'winner' look like a walk (Ben 2026-09-21)."""
    drive = neuro[:, ncol["DRIVE"]]
    walk = drive > 0.8 * drive.max()
    rge = neuro[:, ncol[f"RG_E_{side}"]]
    thr = 0.5 * np.max(rge)
    on = (rge > thr) & walk
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1) + 1
    if len(rises) < 3:
        return []
    return [(rises[i], rises[i + 1]) for i in range(len(rises) - 1)]


def fig45(t, q, neuro, ncol, contact, side="r"):
    """fig4 v2 (2026-09-21, after Ben's critique of the stale render):
    BOTH legs, cycles cut at each foot's own CONTACT loading onsets
    (kine_ref v2), against that leg's OpenSim reference cycle. A leg
    with no cycles is annotated FROZEN - never silently dropped."""
    import kine_ref as KR

    drive = neuro[:, ncol["DRIVE"]]
    wmask = drive > 0.8 * drive.max()
    t0 = float(t[np.flatnonzero(wmask)[0]]) if wmask.any() else 0.0
    ref = KR.ref_cached()
    fig, ax = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    fig.subplots_adjust(hspace=0.14, left=0.09, right=0.97, top=0.93)
    labels = {"hip": "Hip flexion [deg]", "knee": "Knee angle [deg]",
              "ankle": "Ankle angle [deg]"}
    mean_cyc = {}
    for col_i, s in enumerate(("r", "l")):
        mean, per, duty_c, cfrac, used_c, n_cyc, on = KR.sim_side(
            t, q, neuro, t0, s, contact)
        note = (f"{n_cyc} cycles, contact duty {duty_c:.2f} "
                f"(ref {ref[f'duty_{s}']:.2f})" if mean is not None
                else f"FROZEN - no cycles (contact frac {cfrac:.2f})")
        for row, k in enumerate(("hip", "knee", "ankle")):
            a = ax[row, col_i]
            a.plot(KR.GRID, ref[s][k], "k--", lw=1.8,
                   label="OpenSim ref")
            if mean is not None:
                a.plot(KR.GRID, mean[k], lw=2.2, color=C_JOINTS[k],
                       label=f"SNS sim mean ({note})")
            else:
                a.text(0.5, 0.5, note, transform=a.transAxes,
                       ha="center", va="center", fontsize=10,
                       color="#B2182B", fontweight="bold")
            a.set_ylabel(labels[k], fontsize=9)
            if row == 0:
                a.set_title(f"{'Right' if s == 'r' else 'Left'} leg",
                            fontsize=11, fontweight="bold")
                a.legend(fontsize=8, loc="upper right")
            if row == 2:
                a.set_xlabel("gait cycle [%]  (0 = loading onset)")
        if mean is not None:
            mean_cyc[s] = mean
    fig.suptitle("Mean gait cycles vs OpenSim (contact-phased, both "
                 "legs)", fontsize=12)
    fig.savefig(HERE / "fig4_gait_cycles.png", dpi=140)
    plt.close(fig)
    return mean_cyc


def load_benchmark(side):
    """Find an OpenSim IK .mot and extract hip/knee/ankle cycle curves.
    Returns {} if absent. Units are honored via the header's inDegrees
    flag - subject01_walk1_ik.mot says inDegrees=yes, so values are used
    AS-IS (an earlier np.degrees() here double-converted, Ben 2026-09-11)."""
    cands = list((REPO / "Solid_Models" / "OpenSim").rglob("*_ik.mot"))
    if not cands:
        return {}
    f = cands[0]
    print(f"benchmark: {f.name}")
    lines = f.read_text(encoding="utf-8", errors="ignore").splitlines()
    in_deg, end = True, None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("inDegrees"):
            in_deg = s.split("=")[1].strip().lower() == "yes"
        if s == "endheader":
            end = i
            break
    names, rows = None, []
    for ln in lines[end + 1:]:
        s = ln.strip()
        if not s:
            continue
        if names is None:
            names = s.split()
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            break
    data = np.array(rows)
    col = {n: i for i, n in enumerate(names)}
    out = {}
    for k, v in JOINT_KEYS.items():
        cn = v.format(s=side)
        if cn in col:
            y = data[:, col[cn]]
            out[k] = (None, y if in_deg else np.degrees(y))
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


# ---------------------------------------------------------------- fig 6
def fig6(t, j, side, spans):
    """Ben's Figure-8 style: hip/knee/ankle time series, left (black) vs
    right (magenta), one panel per joint."""
    other = "l" if side == "r" else "r"
    d2 = load("spinal_run.npz", other)
    t2, j2 = d2[1], d2[3]   # new 11-tuple: [1]=t, [3]=j dict
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    fig.subplots_adjust(hspace=0.16, left=0.08, right=0.97, top=0.93)
    labels = {"hip": "Hip", "knee": "Knee", "ankle": "Ankle"}
    for i, k in enumerate(("hip", "knee", "ankle")):
        ax[i].plot(t, j[k], lw=1.3, color="m", label="right")
        ax[i].plot(t2, j2[k], lw=1.3, color="k", label="left")
        ax[i].set_ylabel(f"{labels[k]} [deg]")
        for a, b in spans:
            ax[i].axvspan(a, b, color=C_STANCE, zorder=0, lw=0)
    ax[0].legend(fontsize=9, ncol=2, loc="upper right")
    ax[0].set_title("limb joint motion, both legs (shaded = right RG-E "
                    "stance)", fontsize=11)
    ax[-1].set_xlabel("t [s]")
    fig.savefig(HERE / "fig6_limbs_timecourse.png", dpi=140)
    plt.close(fig)


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
    d, t, q, j, act, acols, neuro, ncol, drive, side, contact = \
        load(path, side)
    spans = stance_spans(neuro, ncol, t, side)
    fig1(t, j, neuro, ncol, act, acols, side, spans)
    fig23(t, j, neuro, ncol, act, acols, drive, side, spans)
    fig45(t, q, neuro, ncol, contact, side)
    fig6(t, j, side, spans)
    print("saved fig1_joints_over_neural.png, fig2_joints_over_stimulus.png,")
    print("      fig3_neural_over_stimulus.png, fig4_gait_cycles.png "
          "(v2: both legs, contact-phased),")
    print("      fig6_limbs_timecourse.png")


if __name__ == "__main__":
    main(sys.argv[1:])
