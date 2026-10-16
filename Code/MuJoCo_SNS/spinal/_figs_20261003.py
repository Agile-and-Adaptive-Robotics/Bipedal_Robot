"""Defense-slideshow figures + gifs for the 2026-10-03 curriculum campaign
(Ben's ask: "results pictures and gifs").

Per model npz artifact:
  1. joint-angle overlay figure (hip/knee/ankle, both legs, sim vs the
     OpenSim subject01 reference cycle; Tol palette; journal-style)
  2. animated walking gif: 2-D sagittal stick figure from the joint
     angles, ~2.5 s per gait cycle, looping (matplotlib/Pillow)

Two artifact formats:
  * gait2392 runner npz (t/q/key_joints/neuro/com/contact) - s3k, syn6,
    w2lvar, pruned, spiking-mirror
  * BPA walker npz (t/q/h/f/p/rgE from _curriculum_bpa.run_bpa) - the
    6-joint w2l body

Usage:
  python _figs_20261003.py <npz> <model_label> [--format runner|bpa]
Writes  <label>_joints.png / <label>_walk.gif into
  D:\\Github\\Bipedal_Robot\\Documentation\\Reports and Papers\\Dissertation\\Defense\\slideshow
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = (r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers"
       r"\Dissertation\Defense\slideshow")
# Paul Tol bright palette (CVD-safe)
TOL = dict(blue="#4477AA", red="#EE6677", green="#228833",
           yellow="#CCBB44", cyan="#66CCEE", purple="#AA3377",
           grey="#BBBBBB")
plt.rcParams.update({
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
    "figure.dpi": 150, "savefig.dpi": 200,
    "axes.spines.top": False, "axes.spines.right": False,
})

SEG = dict(thigh=0.42, shank=0.43, foot=0.25)


def _cycles(t, sig, min_gap=0.3):
    """Onset indices of cycles (sig rising above mid-range)."""
    thr = 0.5 * (sig.max() + sig.min())
    on = sig > thr
    st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    keep = [s for i, s in enumerate(st)
            if i == 0 or (t[s] - t[st[i - 1]]) > min_gap]
    return np.array(keep)


def _runner_overlay(z, label):
    t = z["t"]
    q = z["q"]
    joints = [str(x) for x in z["key_joints"]]
    j = {name: i for i, name in enumerate(joints)}
    # walking window: everything after the first 4 s (stand->walk ramp)
    w = t >= 4.0
    fig, axes = plt.subplots(3, 1, figsize=(6.4, 6.4), sharex=True)
    rows = (("hip_flexion", "Hip flexion", "hip"),
            ("knee_angle", "Knee flexion", "knee"),
            ("ankle_angle", "Ankle angle", "ankle"))
    try:
        import kine_ref as KR
        ref = KR.ref_cached()
    except Exception:
        ref = None
    for ax, (key, title, rkey) in zip(axes, rows):
        for side, col in (("r", TOL["blue"]), ("l", TOL["red"])):
            ax.plot(t[w], q[w, j[f"{key}_{side}"]], color=col,
                    lw=1.1, label=f"SNS {side}")
        if ref is not None:
            T = ref[f"T_{side}"]
            for k, cyc in enumerate((0, 1)):
                tc = t[w][0] + k * T + np.linspace(0, T, 100)
                ax.plot(tc, ref[side][rkey], color=TOL["grey"],
                        lw=1.0, ls="--", alpha=0.8,
                        label="OpenSim ref" if (k == 0 and side == "r")
                        else None)
        ax.set_ylabel(f"{title} (deg)")
        ax.legend(loc="upper right", fontsize=7, ncol=3,
                  framealpha=0.8)
    axes[-1].set_xlabel("time (s)")
    axes[0].set_title(f"{label}: joint angles vs OpenSim reference",
                      fontsize=10)
    fig.tight_layout()
    out = os.path.join(OUT, f"{label}_joints.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def _runner_gif(z, label, fps=20, gif_cycle_s=2.5, n_cycles=3):
    t = z["t"]
    q = z["q"]
    joints = [str(x) for x in z["key_joints"]]
    j = {name: i for i, name in enumerate(joints)}
    com = z["com"] if "com" in z.files else None
    w = t >= 4.0
    tw, qw = t[w], q[w]
    knee = -qw[:, j["knee_angle_r"]]     # flexion-positive
    on = _cycles(tw, knee)
    if len(on) < 2:
        on = _cycles(tw, qw[:, j["hip_flexion_r"]])
    if len(on) >= 2:
        T = float(np.median(np.diff(tw[on[:6]]))) if len(on) > 2 else \
            float(tw[on[1]] - tw[on[0]])
    else:
        T = 1.2
    i0 = int(on[-2]) if len(on) >= 2 else 0
    span = min(int(T * n_cycles / (tw[1] - tw[0])), len(tw) - i0 - 1)
    stride = max(1, int(round(T / (gif_cycle_s * fps)
                                / (tw[1] - tw[0]))))
    frames = np.arange(i0, i0 + span, stride)

    def fk(qrow, side, hip_y):
        th = np.radians(qrow[j[f"hip_flexion_{side}"]])
        tk = np.radians(qrow[j[f"knee_angle_{side}"]])
        ta = np.radians(qrow[j[f"ankle_angle_{side}"]])
        hip = np.array([0.0, hip_y])
        knee = hip + SEG["thigh"] * np.array([np.sin(th), -np.cos(th)])
        sa = th + tk                       # shank abs angle (knee<0 flex)
        ank = knee + SEG["shank"] * np.array([np.sin(sa), -np.cos(sa)])
        # foot: neutral = perpendicular to shank pointing forward;
        # dorsiflexion (+) rotates the toe UP (CCW)
        fa = sa + ta
        toe = ank + SEG["foot"] * np.array([np.cos(fa), np.sin(fa)])
        return np.array([hip, knee, ank, toe])

    # self-calibrate the pelvis height over the animated window: shift
    # the figure so the LOWEST heel/toe of the window sits on the ground
    lows = []
    for k in frames:
        qrow = qw[k]
        base = float(com[k, 2]) if com is not None else 0.90
        for side in ("r", "l"):
            lows.append(fk(qrow, side, base)[-1][1])
    dy = float(np.min(lows)) - 0.005

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    ax.set_xlim(-0.55, 0.75)
    ax.set_ylim(-0.12, 1.15)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.plot([-0.6, 0.8], [0, 0], color="#555555", lw=2)
    lnR, = ax.plot([], [], color=TOL["blue"], lw=2.6)
    lnL, = ax.plot([], [], color=TOL["red"], lw=2.6, alpha=0.75)
    ttl = ax.set_title("", fontsize=9)

    def draw(fi):
        k = frames[fi]
        qrow = qw[k]
        hip_y = (float(com[k, 2]) if com is not None else 0.90) - dy
        pR, pL = fk(qrow, "r", hip_y), fk(qrow, "l", hip_y)
        lnR.set_data(pR[:, 0], pR[:, 1])
        lnL.set_data(pL[:, 0], pL[:, 1])
        ttl.set_text(f"{label}   t = {tw[k]:.2f} s")
        return lnR, lnL

    from matplotlib.animation import FuncAnimation, PillowWriter
    anim = FuncAnimation(fig, draw, frames=len(frames), interval=50,
                         blit=False)
    out = os.path.join(OUT, f"{label}_walk.gif")
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out


def _bpa_overlay(z, label):
    t, q = z["t"], z["q"]
    fig, axes = plt.subplots(3, 1, figsize=(6.4, 6.4), sharex=True)
    names = (("hip", "Hip flexion"), ("knee", "Knee flexion"),
             ("ankle", "Ankle"))
    cols = {"L": TOL["blue"], "R": TOL["red"]}
    for r, (key, title) in enumerate(names):
        for side, ci in (("L", 0), ("R", 3)):
            sig = -q[:, ci + r]           # flexion-positive
            axes[r].plot(t, sig, color=cols[side], lw=1.1,
                         label=f"{side}")
        axes[r].set_ylabel(f"{title} (deg)")
        axes[r].legend(loc="upper right", fontsize=7, ncol=2)
    axes[-1].set_xlabel("time (s)")
    axes[0].set_title(f"{label} (BPA walker): air-stepping -> "
                      "retained-support march", fontsize=10)
    fig.tight_layout()
    out = os.path.join(OUT, f"{label}_joints.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def _bpa_gif(z, label, fps=20, gif_cycle_s=2.5, n_cycles=3):
    t, q, h = z["t"], z["q"], z["h"]
    knee = -q[:, 4]
    on = _cycles(t, knee)
    T = float(np.median(np.diff(t[on[:6]]))) if len(on) > 2 else 1.5
    i0 = int(on[-2]) if len(on) >= 2 else 0
    span = min(int(T * n_cycles / (t[1] - t[0])), len(t) - i0 - 1)
    stride = max(1, int(round(T / (gif_cycle_s * fps) / (t[1] - t[0]))))
    frames = np.arange(i0, i0 + span, stride)
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    ax.set_xlim(-0.65, 0.85)
    ax.set_ylim(-0.12, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.plot([-0.7, 0.9], [0, 0], color="#555555", lw=2)
    lnL, = ax.plot([], [], color=TOL["blue"], lw=2.6)
    lnR, = ax.plot([], [], color=TOL["red"], lw=2.6, alpha=0.75)
    ttl = ax.set_title("", fontsize=9)

    def fk(qrow, side, hip_y):
        ci = 0 if side == "L" else 3
        th = np.radians(-qrow[ci])          # flexion-positive -> forward
        tk = np.radians(-qrow[ci + 1])
        ta = np.radians(-qrow[ci + 2])
        hip = np.array([0.0, hip_y])
        knee = hip + SEG["thigh"] * np.array([np.sin(th), -np.cos(th)])
        sa = th + tk
        ank = knee + SEG["shank"] * np.array([np.sin(sa), -np.cos(sa)])
        fa = sa + ta
        toe = ank + SEG["foot"] * np.array([np.cos(fa), np.sin(fa)])
        return np.array([hip, knee, ank, toe])

    # self-calibrate pelvis height: lowest toe of the window on ground
    lows = [fk(q[k], s, float(h[k, 0]))[-1][1]
            for k in frames for s in ("L", "R")]
    dy = float(np.min(lows)) - 0.005

    def draw(fi):
        k = frames[fi]
        hip_y = float(h[k, 0]) - dy
        lnL.set_data(*fk(q[k], "L", hip_y).T)
        lnR.set_data(*fk(q[k], "R", hip_y).T)
        ttl.set_text(f"{label}   t = {t[k]:.2f} s")
        return lnL, lnR

    from matplotlib.animation import FuncAnimation, PillowWriter
    anim = FuncAnimation(fig, draw, frames=len(frames), interval=50)
    out = os.path.join(OUT, f"{label}_walk.gif")
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out


def main():
    npz = sys.argv[1]
    label = sys.argv[2]
    fmt = "runner"
    if "--format" in sys.argv:
        fmt = sys.argv[sys.argv.index("--format") + 1]
    z = np.load(npz, allow_pickle=True)
    os.makedirs(OUT, exist_ok=True)
    if fmt == "bpa":
        p1 = _bpa_overlay(z, label)
        p2 = _bpa_gif(z, label)
    else:
        p1 = _runner_overlay(z, label)
        p2 = _runner_gif(z, label)
    print("wrote", p1)
    print("wrote", p2)


if __name__ == "__main__":
    main()
