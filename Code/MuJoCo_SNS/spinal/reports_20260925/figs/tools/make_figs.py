"""FIGURES task (2026-09-25): build the goal-4 variant-campaign figures
from the stage-5 winner replay npz's (see replay_s5.py).

Outputs (into reports_20260925/figs/):
  w2lvar_s5_walk_overlay.png   sim mean cycles vs the OpenSim reference
  syn6_s5_walk_overlay.png
  w2lvar_s5_traces.png         joints + tilt/COM + contact + neural raster
  syn6_s5_traces.png           (syn6 lanes are the S1..S4 synergy channels)
  variant_stage_scores.png     the 10 curriculum_{variant}_stage{1..5}.json

Conventions (kine_ref.py): OpenSim/ISB throughout - hip_flexion
+flexion, knee_angle NEGATIVE = flexion, ankle_angle +dorsiflexion;
X anterior, Y up, Z right.

Usage: python make_figs.py [w2lvar|syn6|scores]   (default: everything)
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

import kine_ref  # noqa: E402

FIGS = os.path.join("reports_20260925", "figs")
# Okabe-Ito CVD-safe palette (lab convention)
C_SIM_R = "#0072B2"   # blue
C_SIM_L = "#D55E00"   # vermillion
C_REF = "#000000"
WALK_T0, WALK_T1 = 5.0, 15.0   # params.SCHEDULE["walk"]
SCORES = {"w2lvar": -165.0405055213092, "syn6": -197.21456977796157}
STAGE_NAMES = ("1 air\ndeaff", "2 air\nafferented", "3 standing\nbalance",
               "4 walk\nno contact", "5 walk\ncontact")


def load_npz(variant):
    z = np.load(os.path.join(FIGS, "tmp", f"figs_npz_{variant}.npz"),
                allow_pickle=True)
    return (z["t"], z["q"], z["neuro"], z["contact"], z["com"],
            [str(x) for x in z["neuro_names"]],
            [str(x) for x in z["key_joints"]])


def kine_metrics(variant):
    path = os.path.join(FIGS, "tmp", f"replay_{variant}_s5.log")
    out = {}
    for ln in open(path, encoding="utf-8", errors="replace"):
        if "kine_score=" in ln:
            for tok in ln.split():
                for key in ("kine_score", "kz", "tilt_max", "duty",
                            "T_r", "lag_rl", "n_r=", "n_l="):
                    if tok.startswith(key):
                        out[key.rstrip("=")] = float(
                            tok.split("=")[-1].rstrip(","))
                if tok.startswith("bilateral"):
                    out["bilateral"] = tok.split("=")[-1]
        if "RECOMPUTED" in ln:
            out["recomputed"] = float(ln.split("RECOMPUTED")[1]
                                      .split("vs json")[0])
    return out


# ---------------------------------------------------------------- overlay
def sim_cycles(t, q, neuro, contact, side, ref):
    """Mean cycle + per-cycle list for one side (kine_ref semantics).
    Also returns the side's real contact duty (frac of walk window
    loaded > 50 N)."""
    m = t >= WALK_T0
    tt = t[m]
    si = 0 if side == "r" else 1
    cols = (3, 4, 5) if side == "r" else (8, 9, 10)
    force = contact[m, si]
    duty = float(np.mean(force > kine_ref.LOAD_N))
    on = kine_ref._loading_onsets(tt, force)
    if len(on) < 3:
        rge_col = 2 if side == "r" else 4
        rge = neuro[m, rge_col]
        bon = (rge > 0.5 * max(np.max(rge), 1e-9)).astype(int)
        on = tt[np.flatnonzero(np.diff(bon) == 1)]
    mean, periods, cycles = kine_ref._cycles_from_onsets(
        t[m], tt, on, cols, q[m])
    return mean, periods, cycles, duty


def fig_overlay(variant):
    t, q, neuro, contact, com, names, joints = load_npz(variant)
    ref = kine_ref.ref_cached()
    km = kine_metrics(variant)
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 8.2), sharex=True)
    ttl = (f"{variant} stage-5 winner - ground walk mean cycles vs "
           f"OpenSim reference (subject01)\n"
           f"score {SCORES[variant]:.1f}   kine {km.get('kine_score', float('nan')):.1f}"
           f"   kz {km.get('kz', float('nan')):.2f}"
           f"   tilt_max {km.get('tilt_max', float('nan')):.1f} deg"
           f"   bilateral {km.get('bilateral', '?')}")
    fig.suptitle(ttl, fontsize=11)
    for col, side in enumerate(("r", "l")):
        mean, periods, cycles, duty_c = sim_cycles(t, q, neuro, contact,
                                                   side, ref)
        cc = C_SIM_R if side == "r" else C_SIM_L
        for row, j in enumerate(("hip", "knee", "ankle")):
            ax = axes[row, col]
            if cycles:
                arr = np.vstack([c[j] for c in cycles])
                ax.fill_between(kine_ref.GRID, arr.mean(0) - arr.std(0),
                                arr.mean(0) + arr.std(0),
                                color=cc, alpha=0.18, lw=0)
            if mean is not None:
                ax.plot(kine_ref.GRID, mean[j], color=cc, lw=2.0,
                        label=f"sim {side} (n={len(periods)})")
            ax.plot(kine_ref.GRID, ref[side][j], color=C_REF, lw=1.4,
                    ls="--", label=f"OpenSim ref {side}")
            ax.set_ylabel({"hip": "hip flexion (deg, + = flexion)",
                           "knee": "knee angle (deg, - = flexion)",
                           "ankle": "ankle angle (deg, + = dorsiflexion)"}[j],
                          fontsize=9)
            if row == 0:
                T = float(np.mean(periods)) if periods else float("nan")
                ax.set_title(
                    f"{'RIGHT' if side == 'r' else 'LEFT'} leg  "
                    f"T {T:.2f} s (ref {ref[f'T_{side}']:.2f})  "
                    f"duty {duty_c:.2f}"
                    f" (ref {ref[f'duty_{side}']:.2f})  "
                    f"cycles {len(periods)}", fontsize=10)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8, loc="best")
    for ax in axes[-1]:
        ax.set_xlabel("gait cycle (%)")
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    fig.text(0.5, 0.015,
             "Cycles cut at each foot's own contact-loading onsets "
             "(t >= 5 s walk window); mean +/- 1 sd across cycles. "
             "OpenSim conventions: hip +flexion, knee -flexion, "
             "ankle +dorsiflexion. Score = stage-5 objective "
             "(kine, clipped -315; kz<0.62 -20; tilt>40 deg -10).",
             ha="center", fontsize=8, color="0.3")
    out = os.path.join(FIGS, f"{variant}_s5_walk_overlay.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


# ----------------------------------------------------------------- traces
def fig_traces(variant):
    t, q, neuro, contact, com, names, joints = load_npz(variant)
    km = kine_metrics(variant)
    syn6 = variant == "syn6"
    fig = plt.figure(figsize=(12.5, 10.0))
    gs = fig.add_gridspec(4, 1, height_ratios=(3, 1.4, 1.2, 3.2),
                          hspace=0.32)
    ax_j = fig.add_subplot(gs[0])
    ax_p = fig.add_subplot(gs[1], sharex=ax_j)
    ax_c = fig.add_subplot(gs[2], sharex=ax_j)
    ax_n = fig.add_subplot(gs[3], sharex=ax_j)
    axs = (ax_j, ax_p, ax_c, ax_n)
    ttl = (f"{variant} stage-5 winner - ground walk traces (22 s eval; "
           f"grey = DRIVE walk window {WALK_T0:.0f}-{WALK_T1:.0f} s)\n"
           f"score {SCORES[variant]:.1f}   kine {km.get('kine_score', float('nan')):.1f}"
           f"   kz {km.get('kz', float('nan')):.2f}"
           f"   tilt_max {km.get('tilt_max', float('nan')):.1f} deg")
    fig.suptitle(ttl, fontsize=11)

    # 1) joint angles (deg) r + l
    for j, lab in ((3, "hip_flexion"), (4, "knee_angle"), (5, "ankle_angle")):
        ax_j.plot(t, q[:, j], color=C_SIM_R,
                  ls=("-", "--", ":")[("hip", "knee", "ankle").index(
                      lab.split("_")[0])],
                  lw=1.4, label=f"{lab}_r")
        ax_j.plot(t, q[:, j + 5], color=C_SIM_L,
                  ls=("-", "--", ":")[("hip", "knee", "ankle").index(
                      lab.split("_")[0])],
                  lw=1.4, alpha=0.85, label=f"{lab}_l")
    ax_j.set_ylabel("joint angle (deg)\nhip +flex / knee -flex / ank +dorsi")
    ax_j.legend(fontsize=7, ncol=3, loc="upper right")
    ax_j.set_title("OpenSim key joints (r solid-family blue, l vermillion)",
                   fontsize=9)

    # 2) pelvis tilt + COM height
    ax_p.plot(t, q[:, 0], color="#009E73", lw=1.2, label="pelvis_tilt")
    ax_p.set_ylabel("pelvis tilt (deg)")
    ax2 = ax_p.twinx()
    ax2.plot(t, com[:, 2], color="#CC79A7", lw=1.2, label="COM z (m)")
    ax2.set_ylabel("COM z (m)", color="#CC79A7")
    ax_p.axhline(0.0, color="0.6", lw=0.6, ls=":")
    ax_p.set_title("pelvis tilt (green; + = pitched back) and COM height "
                   "(pink, m)", fontsize=9)

    # 3) per-foot contact normal force
    ax_c.plot(t, contact[:, 0], color=C_SIM_R, lw=1.1, label="contact_r (N)")
    ax_c.plot(t, contact[:, 1], color=C_SIM_L, lw=1.1, label="contact_l (N)")
    ax_c.axhline(kine_ref.LOAD_N, color="0.5", lw=0.7, ls=":",
                 label=f"load thr {kine_ref.LOAD_N:.0f} N")
    ax_c.set_ylabel("heel+toe\nnormal force (N)")
    ax_c.legend(fontsize=7, ncol=3, loc="upper right")
    ax_c.set_title("per-foot ground contact (real MuJoCo normal force)",
                   fontsize=9)

    # 4) neural raster: DRIVE/POSTURE + RG + PF watch lanes
    lanes = []
    for nm in names:
        lanes.append(nm)
    off = 0.0
    yticks, ylabels = [], []
    for nm in lanes:
        sig = neuro[:, names.index(nm)]
        rng = float(np.ptp(sig)) or 1.0
        ax_n.plot(t, (sig - np.mean(sig)) / rng + off, lw=0.8,
                  color="#56B4E9" if nm.startswith(("PF_S", "PF_")) else
                  ("#E69F00" if nm.startswith("BAL") else "#333333"))
        yticks.append(off)
        ylabels.append(nm)
        off -= 1.25
    ax_n.set_yticks(yticks)
    ax_n.set_yticklabels(ylabels, fontsize=7)
    ax_n.set_ylim(off + 0.6, 1.4)
    lane_title = ("PF_S1..S4_r = the syn6 synergy-PF channels (watch "
                  "selection in runner.py:1148)" if syn6 else
                  "PF_HIP-E/KNEE-E/KNEE-F/ANK-F_r = the w2lvar merged "
                  "joint-layer PF channels (runner.py:1150)")
    ax_n.set_title(f"neural traces per lane: ptp-normalised deviation "
                   f"from mean, stacked {lane_title}", fontsize=9)
    ax_n.set_xlabel("t (s)")

    for ax in axs:
        ax.axvspan(WALK_T0, WALK_T1, color="0.85", zorder=0)
        ax.grid(alpha=0.25)
        ax.set_xlim(0.0, float(t[-1]))
    for ax in axs[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)
    fig.tight_layout(rect=(0, 0.03, 1, 0.93))
    note = ("S1..S4 lanes: burst count visible per synergy channel; the "
            "chain report documents S2/S3 near-degenerate within-family "
            "channels." if syn6 else
            "KNEE-E drives knee_ext + ankle_pf and KNEE-F drives "
            "knee_flex + ankle_df (W2L biarticular synergy; deviation D4).")
    fig.text(0.5, 0.005, note, ha="center", fontsize=8, color="0.3")
    out = os.path.join(FIGS, f"{variant}_s5_traces.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


# ----------------------------------------------------------------- scores
def fig_scores():
    variants = ("w2lvar", "syn6")
    cols = {"w2lvar": C_SIM_R, "syn6": C_SIM_L}
    data = {v: [] for v in variants}
    for v in variants:
        for st in range(1, 6):
            J = json.loads(open(
                os.path.join(SPINAL, f"curriculum_{v}_stage{st}.json"),
                encoding="utf-8").read())
            data[v].append(float(J["score"]))
    x = np.arange(5)
    w = 0.38
    fig, ax = plt.subplots(figsize=(10.0, 5.6))
    for i, v in enumerate(variants):
        bars = ax.bar(x + (i - 0.5) * w, data[v], w, color=cols[v],
                      label=(f"w2lvar (W2L-on-gait2392, 888n)"
                             if v == "w2lvar" else
                             f"syn6 (6-synergy walker, 794n)"))
        for b, val in zip(bars, data[v]):
            ax.annotate(f"{val:.1f}", (b.get_x() + b.get_width() / 2, val),
                        textcoords="offset points",
                        xytext=(0, 3 if val >= 0 else -11),
                        ha="center", fontsize=8)
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(STAGE_NAMES)
    ax.set_ylabel("stage-winner objective score")
    ax.set_title("goal-4 variant curricula - five stage-winner scores per "
                 "variant (2026-09-25 campaign)")
    ax.grid(alpha=0.3, axis="y")
    ax.legend(loc="lower left", fontsize=9)
    ax.text(0.5, -0.30,
            "Stages 1-3: higher is better (air rhythm 3*rises + 0.5*(-knee_min), "
            "rhythm-gated; standing 100 - 400*sway - tilt - 40*|0.5-sym|).\n"
            "Stages 4-5: kine objective (<= 0, 0 = perfect; -315 clip, "
            "-320 frozen sentinel; stage 5 adds kz<0.62 (-20) and "
            "tilt>40 deg (-10) penalties). Dbs: optuna_{w2lvar,syn6}.db.",
            transform=ax.transAxes, ha="center", fontsize=8, color="0.3")
    fig.tight_layout()
    out = os.path.join(FIGS, "variant_stage_scores.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "w2lvar"):
        fig_overlay("w2lvar")
        fig_traces("w2lvar")
    if what in ("all", "syn6"):
        fig_overlay("syn6")
        fig_traces("syn6")
    if what in ("all", "scores"):
        fig_scores()
