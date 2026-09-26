r"""FIG BUILDER — assembles the 2026-09-25 walker-campaign deliverable
figures/GIFs from the w2l_runs.py outputs (data_*.npz + frames_*/\*.png).
STANDALONE, no simulation. Usage:
    C:\Users\Ben Bolen\.conda\envs\myo\python.exe w2l_figs.py m3|m4|m5|m6|all
"""
from __future__ import annotations

import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
from PIL import Image, ImageDraw, ImageEnhance

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TOOLS = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\figs\tools"
FIGS = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\figs"

RW, RH = 360, 470            # raw render size (w2l_runs.py)
IW = 430                     # inset panel width px
C_L, C_R = "#1f77b4", "#ff7f0e"
C_GL, C_GR = "#2ca02c", "#d62728"
C_HL, C_HR = "#9467bd", "#8c564b"
WIN = 4.0                    # sliding window, s


# --------------------------------------------------------------- inset draw
def inset_png(rows, overview, tcur, title, footer, vlines=(),
              w=IW, h=RH) -> Image.Image:
    """Render the right-hand inset panel to a PIL image.
    rows: list of (ax_title, [(label, t, y, color, lw), ...])
    overview: (label, t, y, color) or None — full-run strip with cursor.
    vlines: list of (t, label, color) vertical markers on every axes."""
    dpi = 100
    fw, fh = w / dpi, h / dpi
    n = len(rows)
    row_h, gap = 0.19, 0.052
    fig = plt.figure(figsize=(fw, fh), dpi=dpi)
    axes = []
    for k in range(n):
        bottom = 0.955 - 0.048 - row_h - k * (row_h + gap)
        a = fig.add_axes([0.12, bottom, 0.84, row_h])
        axes.append(a)
    a_ov = fig.add_axes([0.12, 0.085, 0.84, 0.068])
    for k, (atitle, traces) in enumerate(rows):
        a = axes[k]
        w0 = max(0.0, tcur - WIN)
        for label, t, y, color, lw in traces:
            m = (t >= w0) & (t <= max(tcur, WIN))
            if m.any():
                a.plot(t[m] - w0, y[m], color=color, lw=lw, label=label)
        a.set_xlim(0.0, WIN)
        a.set_title(atitle, fontsize=7, pad=1.5)
        a.tick_params(labelsize=6)
        a.grid(alpha=0.25, lw=0.4)
        if k == 0:
            a.legend(fontsize=6, loc="upper right", framealpha=0.7,
                     handlelength=1.2, borderpad=0.2)
        for tv, lab, color in vlines:
            if w0 <= tv <= max(tcur, WIN):
                a.axvline(tv - w0, color=color, ls="--", lw=1.0, alpha=0.9)
    if overview is not None:
        label, t, y, color = overview
        ds = max(1, len(t) // 8000)
        a_ov.plot(t[::ds], y[::ds], color=color, lw=0.6)
        a_ov.axvline(tcur, color="k", lw=1.0)
        for tv, _, color in vlines:
            a_ov.axvline(tv, color=color, ls="--", lw=0.9)
        a_ov.set_xlim(t[0], t[-1])
        a_ov.tick_params(labelsize=5)
        a_ov.set_title(label + " (full run, cursor = frame time)",
                       fontsize=6, pad=1.0)
    fig.suptitle(title, fontsize=8.5, y=0.995)
    fig.text(0.5, 0.004, footer, fontsize=5.6, ha="center", va="bottom",
             linespacing=1.25)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB").resize((w, h))


def build_gif(run: str, outfile: str, title: str, footer: str,
              rows_fn, overview_fn, vlines=(), data_name: str = "") -> None:
    fdir = os.path.join(TOOLS, "frames_" + run)
    frames = sorted(f for f in os.listdir(fdir) if f.endswith(".png"))
    data = np.load(os.path.join(TOOLS,
                                f"data_{data_name or run}.npz"))
    ftimes = data["ftimes"]
    pil_frames = []
    for k, fname in enumerate(frames):
        if k >= len(ftimes):
            break
        tcur = float(ftimes[k])
        base = Image.open(os.path.join(fdir, fname)).convert("RGB")
        base = ImageEnhance.Brightness(base).enhance(1.22)
        base = ImageEnhance.Contrast(base).enhance(1.05)
        ins = inset_png(rows_fn(tcur, data), overview_fn(data), tcur,
                        title, footer, vlines=vlines)
        comp = Image.new("RGB", (base.width + ins.width,
                                 max(base.height, ins.height)), "white")
        comp.paste(base, (0, 0))
        comp.paste(ins, (base.width, 0))
        ImageDraw.Draw(comp).text((4, 4), "t = %.1f s" % tcur,
                                  fill="yellow")
        pil_frames.append(comp)
    out = os.path.join(FIGS, outfile)
    pil_frames[0].save(out, save_all=True, append_images=pil_frames[1:],
                       duration=100, loop=0)
    mb = os.path.getsize(out) / 1e6
    print(f"[gif] {outfile}: {len(pil_frames)} frames, {mb:.1f} MB")


# ------------------------------------------------------------------ M3 gif
def m3_rows(tcur, data):
    t, q, rg = data["t"], data["q"], data["rgL"]
    hipL, hipR = -q[:, 0], -q[:, 3]
    r1 = ("hip flexion L/R (deg, window)",
          [("L", t, hipL, C_L, 1.0), ("R", t, hipR, C_R, 1.0)])
    r2 = ("L RG ext (mV, window)", [("RG-E L", t, rg, C_GL, 1.0)])
    return [r1, r2]


def m3_overview(data):
    t, q = data["t"], data["q"]
    return ("hip flexion L", t, -q[:, 0], C_L)


# ------------------------------------------------------------------ M4 gif
def m4_rows(tcur, data):
    t, q = data["t"], data["q"]
    rgL, rgR = data["rgL"], data["rgR"]
    r1 = ("hip flexion L/R (deg, window)",
          [("L", t, -q[:, 0], C_L, 1.0), ("R", t, -q[:, 3], C_R, 1.0)])
    r2 = ("RG ext L/R (mV, window)",
          [("RG-E L", t, rgL, C_GL, 1.0), ("RG-E R", t, rgR, C_GR, 1.0)])
    return [r1, r2]


def m4_overview(data):
    t, q = data["t"], data["q"]
    return ("hip flexion L/R", t, -q[:, 0], C_L)


# ------------------------------------------------------------------ M5 gif
def m5_rows(tcur, data):
    t, heel, rg = data["t"], data["heelSN"], data["rgL"]
    r1 = ("heel SN L/R (mV; scripted 1.5, extra 4.0 pulse)",
          [("heel L", t, heel[:, 0], C_HL, 1.0),
           ("heel R", t, heel[:, 1], C_HR, 1.0)])
    r2 = ("RG ext L/R (mV, window)",
          [("RG-E L", t, rg, C_GL, 1.0), ("RG-E R", t, data["rgR"], C_GR, 1.0)])
    return [r1, r2]


def m5_overview(data):
    t, heel = data["t"], data["heelSN"]
    return ("heel-L SN", t, heel[:, 0], C_HL)


# ------------------------------------------------------------ static figs
def fig_m3_joints() -> None:
    data = np.load(os.path.join(TOOLS, "data_m3.npz"))
    t, q = data["t"], data["q"]
    an = t >= 2.0
    flex = -q[an]
    tt = t[an]
    refs = [("hip", (0, 3), 38.0), ("knee", (1, 4), 61.0),
            ("ankle", (2, 5), 16.0)]
    fig, axs = plt.subplots(3, 1, figsize=(9.5, 8.0), sharex=True)
    for ax, (name, (cl, cr), ref) in zip(axs, refs):
        ax.plot(tt, flex[:, cl], color=C_L, lw=0.8, label="L")
        ax.plot(tt, flex[:, cr], color=C_R, lw=0.8, label="R")
        rl = float(flex[:, cl].ptp())
        rr = float(flex[:, cr].ptp())
        mid = 0.5 * (flex[:, cl].max() + flex[:, cl].min())
        ax.axhspan(mid - ref / 2, mid + ref / 2, color="gray", alpha=0.15)
        ax.text(0.995, 0.04,
                f"sim range L {rl:.1f} / R {rr:.1f} deg   "
                f"(AnimatLab ref ~{ref:.0f})",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                bbox=dict(fc="white", alpha=0.8, lw=0.5))
        ax.set_ylabel(f"{name} flexion (deg)")
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8, ncol=2, loc="upper right")
    axs[-1].set_xlabel("time (s)")
    fig.suptitle("M3: W2L 2023-original (single LH RG) air stepping, MuJoCo "
                 "body — joint angles vs AnimatLab references\n"
                 "measured cadence ~0.97 Hz (RG-E period 1.027 s); "
                 "2023-original ref 0.77 Hz (within 2x), modern W2L air ref "
                 "2.22 Hz\n"
                 "knobs: te=3, tf=4, tau_h=0.25 s, ctrl cap 0.5, joint "
                 "damping 3.0 N·m·s/rad, pelvis lifted +0.30 m (air rig)",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(FIGS, "w2l_air_joints_vs_reference.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("[png]", out)


def fig_m6() -> None:
    data = np.load(os.path.join(TOOLS, "data_m6.npz"))
    t, q = data["t"], data["q"]
    com, pel, tilt = data["com"], data["pel"], data["tilt"]
    heelF, toeF = data["heelF"], data["toeF"]
    rg = data["rgE"]
    heelVmax = float(data["heelVmax"])
    w = t >= 1.0
    # RG-E burst onsets (rising through 0.5 max, refractory 0.4 s) — same
    # definition as test_w2l_air.bursts_refractory, inlined (no path needed)
    DT_PHY = 0.001
    rgw = rg[w, 0]
    on = rgw > 0.5 * rgw.max()
    starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    keepL = []
    for x in starts:
        if not keepL or (x - keepL[-1]) * DT_PHY > 0.4:
            keepL.append(x)
    stL = np.array(keepL)
    per = (float(np.diff(stL).mean() * DT_PHY) if len(stL) >= 3
           else float("nan"))
    duty = {s: float((heelF[:, k] > 5.0).mean()) for k, s in enumerate("LR")}
    harness = 100.0 * data["rigFz"][w].mean() / 411.0

    fig, axs = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.2, 1.0, 1.0]))
    ax = axs[0]
    ax.plot(t, -q[:, 0], color=C_L, lw=0.8, label="hip L")
    ax.plot(t, -q[:, 3], color=C_R, lw=0.8, label="hip R")
    ax.plot(t, -q[:, 1], color=C_L, lw=0.8, ls="--", label="knee L")
    ax.plot(t, -q[:, 4], color=C_R, lw=0.8, ls="--", label="knee R")
    ax.plot(t, -q[:, 2], color=C_L, lw=0.8, ls=":", label="ankle L")
    ax.plot(t, -q[:, 5], color=C_R, lw=0.8, ls=":", label="ankle R")
    ax.set_ylabel("flexion (deg)")
    ax.set_title("joint angles: marches in place at the autonomous RG period "
                 f"({per:.3f} s, {1.0 / per:.2f} Hz)", fontsize=9)
    ax.legend(fontsize=7, ncol=6, loc="upper right")
    ax.grid(alpha=0.3)

    ax = axs[1]
    ax.plot(t, pel[:, 2], color="k", lw=0.9, label="pelvis z (m)")
    ax.plot(t, com[:, 2], color=C_GL, lw=0.9, label="COM z (m)")
    ax2 = ax.twinx()
    ax2.plot(t, tilt, color="#d62728", lw=0.6, alpha=0.7)
    ax2.set_ylabel("tilt (deg)", color="#d62728", fontsize=8)
    ax2.tick_params(labelsize=7, colors="#d62728")
    ax.set_ylabel("height (m)")
    ax.set_title(f"harness-supported (rig S=1.0): harness carries "
                 f"{harness:+.0f}% of body weight; pelvis z min "
                 f"{pel[w, 2].min():.3f} m, tilt max {tilt[w].max():.1f} deg",
                 fontsize=9)
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(alpha=0.3)

    ax = axs[2]
    ax.plot(t, heelF[:, 0], color=C_L, lw=0.8, label="heel L")
    ax.plot(t, heelF[:, 1], color=C_R, lw=0.8, label="heel R")
    ax.plot(t, toeF[:, 0], color=C_L, lw=0.7, ls="--", label="toe L")
    ax.plot(t, toeF[:, 1], color=C_R, lw=0.7, ls="--", label="toe R")
    ax.axhline(20.0, color="gray", lw=0.7, ls=":")
    ax.set_ylabel("plate normal force (N)")
    ax.set_xlabel("time (s)")
    ax.set_title(f"contact loading: ZERO heel loading (heel duty L/R "
                 f"{duty['L']:.2f}/{duty['R']:.2f}, heel SN max "
                 f"{heelVmax:.2f} mV) — heel stance reset never engages",
                 fontsize=9)
    ax.legend(fontsize=7, ncol=4, loc="upper right")
    ax.grid(alpha=0.3)

    fig.suptitle("M6: ground, pelvis FREE, harness-supported in-place march "
                 "(rig S=1.0, kxy=kz=2000 N/m, krot 400 N·m/rad)\n"
                 "the honest result: DOES NOT WALK — the rest pose cannot "
                 "free-stand (COM starts 3.09 cm outside the support polygon; "
                 "free stand topples at 0.89 s)\n"
                 "heels never load (duty 0.00, heel SN max 0.00 mV), so the "
                 "contact-driven layer stays open — what remains is the "
                 "autonomous M4 rhythm", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(FIGS, "w2l_ground_supported_walk.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("[png]", out)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("m3", "all"):
        build_gif("m3", "w2l_air_stepping.gif",
                  "M3: W2L 2023 original (single LH RG), AIR stepping",
                  "knobs te=3 tf=4 tau_h=0.25 s cap=0.5 damp=3.0 lift=+0.30 m\n"
                  "RG-E period 1.027 s (0.97 Hz), hips antiphase | backed by "
                  "goal2_m3_w2l_cpg_air.md", m3_rows, m3_overview)
        fig_m3_joints()
    if which in ("m4", "all"):
        build_gif("m4", "w2l_split_rg.gif",
                  "M4: SPLIT RGs, commissural coupling (c1 3.0 / V3 0.08)",
                  "each leg on its OWN RG; coupling LOCKS the pair at "
                  "1.357 s (0.74 Hz)\nL/R phase 0.507 cycle, band r -0.941 | "
                  "backed by goal2_m4_rg_split.md",
                  m4_rows, m4_overview)
    if which in ("m5", "all"):
        d5 = np.load(os.path.join(TOOLS, "data_m5p.npz"))
        tstar = float(d5["tstar"])
        build_gif("m5", "w2l_afferented_air.gif",
                  "M5: AFFERENTED split-RG air walk (Ben's rules)",
                  "scripted contacts amp 1.5 nA: heel [onset+0.02,+0.40T], "
                  "toe [0.62T,0.92T]; Ib 1 nA\nextra causal heel-L 4 nA x "
                  f"0.5 s at t*={tstar:.2f} s -> RG-E onset shift +54 ms "
                  "(gate b PASS)\nbacked by goal2_m5_afferents.md",
                  m5_rows, m5_overview,
                  vlines=[(tstar, "heel-L extra pulse", "#e377c2")],
                  data_name="m5p")
    if which in ("m6", "all"):
        fig_m6()
    print("done")
