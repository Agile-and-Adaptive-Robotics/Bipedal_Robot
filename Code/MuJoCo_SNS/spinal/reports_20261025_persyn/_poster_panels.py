"""Deng-poster-style perturbation panels on the W2L neuromechanical model
(frozen per-synapse w25 winner; wsup 0.25; drive pinned te=3/tf=4).

Panels: (A) control knee angles + heel contact; (B,C) inhibitory tonic
stimulus to the rhythm layer (mild/strong, 1 s window); (D,E) external
forward/backward push (40 N, 0.15 s at the pelvis); (F) pathway-deletion
summary (crossed heel->contra-RG, Ia chain, Renshaw, heel cascade,
commissural, Ib) on the same walker.

Style: Tol muted palette, Arial (DejaVu fallback), 10 pt floor, L solid /
R dashed (never hue alone). Outputs PNG+PDF+alt txt into Revision_D figs.
"""
from __future__ import annotations

import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np  # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")
import test_w2l_ground as T  # noqa: E402

W25 = json.load(open(os.path.normpath(os.path.join(
    HERE, "..", "reports_20261025_persyn", "round_w25.json")),
    encoding="utf-8"))["params"]
BASE = dict(T.KV)
BASE.update(dict(net="bilateralrg2", dur_walk=16.0, rig=0.25, wsup=0.25,
                 sagittal=50.0, te=3.0, tf=4.0, pushN=0.0, stimN=0.0,
                 legcoll=0.0, sink=0.06, walker_yield=1.5, **W25))

TOL = dict(blue="#4477AA", red="#EE6677", teal="#228833",
           yellow="#CCBB44", purple="#AA3377", grey="#BBBBBB")
plt.rcParams.update({"font.size": 10, "font.family": "sans-serif",
                     "axes.linewidth": 0.8})


def run(tag, **over):
    KV = dict(BASE)
    KV.update({k: v for k, v in over.items() if v is not None})
    T.KV.clear(); T.KV.update(KV)
    log = T.run_phase("walk", KV["dur_walk"], KV["rig"], None)
    t = np.array(log["t"]); q = np.array(log["q"])   # (n, nv) joint angles
    heelF = np.array(log["heelF"]); rgE = np.array(log["rgE"])
    com = np.array(log["com"]); tilt = np.array(log["tilt"])
    finite = bool(log["finite"]); fall_t = log["fall_t"]
    return dict(tag=tag, t=t, q=q, heelF=heelF, rgE=rgE, com=com,
                tilt=tilt, finite=finite, fall_t=fall_t)


def cycles(heelF, t, thr=20.0, merge=0.45):
    f = heelF[:, 0]
    on = f > thr
    st = t[1:][on[1:] & ~on[:-1]]
    m = []
    for x in st:
        if not m or (x - m[-1]) > merge:
            m.append(x)
    return len(m)


# ---- conditions -----------------------------------------------------------
IA_KEYS = [k for k in W25 if k.startswith("ia_")]
RC_KEYS = [k for k in W25 if k.startswith(("rc_own", "rc_inh", "rc_mut"))]
HEEL_KEYS = [k for k in W25 if k.startswith("heel_")]
IB_KEYS = [k for k in W25 if k.startswith("ib_")] + ["ibnA"]

runs = {}
runs["control"] = run("control")
runs["stim_mild"] = run("stim_mild", stimN=1.5, stimT=8.0, stimDur=1.0)
runs["stim_strong"] = run("stim_strong", stimN=3.0, stimT=8.0, stimDur=1.0)
runs["push_fwd"] = run("push_fwd", pushN=40.0, pushT=8.0, pushAxis="+x")
runs["push_back"] = run("push_back", pushN=40.0, pushT=8.0, pushAxis="-x")
dels = {
    "crossed heel edge": {"ch_contra": 0.0},
    "Ia chain (all)": {k: 0.0 for k in IA_KEYS},
    "Renshaw (all)": {k: 0.0 for k in RC_KEYS},
    "heel cascade": {k: 0.0 for k in HEEL_KEYS},
    "commissural": {"comm": 0.0},
    "Ib afferent": {k: 0.0 for k in IB_KEYS},
}
del_results = {}
for name, over in dels.items():
    r = run("del", **over)
    del_results[name] = dict(
        cycL=cycles(r["heelF"], r["t"]),
        cycR=cycles(np.nan_to_num(r["heelF"]), r["t"]),
        finite=r["finite"], fell=r["fall_t"] is not None,
        kneeL_deg=float(np.degrees(r["q"][:, 0]).max() -
                        np.degrees(r["q"][:, 0]).min()) if r["finite"] else 0,
        last_t=float(r["t"][-1]))
    print(f"  DEL {name:26s} cycL {del_results[name]['cycL']:3d} "
          f"cycR {del_results[name]['cycR']:3d} "
          f"fell={del_results[name]['fell']} "
          f"end={del_results[name]['last_t']:.1f}s", flush=True)

# joint angle column indices: find hip/knee/ankle in q via joint names
# (q follows qpos order; w2l joints: root(6) + hip_L,knee_L,ankle_L,...)
def col(name):
    import mujoco
    return None  # placeholder; use fixed offsets below

# w2l_mjcf joint order: root free(7 qpos) then per-side hip,knee,ankle
QOFF = {"hip_L": 0, "knee_L": 1, "ankle_L": 2,
        "hip_R": 3, "knee_R": 4, "ankle_R": 5}   # log["q"] = JNTS order, deg

# ---- figure ---------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.6), constrained_layout=True)

def knee_panel(ax, r, title, shade=None, note=""):
    t = r["t"]
    kL = -r["q"][:, QOFF["knee_L"]]   # negate: flexion-positive
    kR = -r["q"][:, QOFF["knee_R"]]
    ax.plot(t, kL, color=TOL["blue"], lw=1.4, label="left", ls="-")
    ax.plot(t, kR, color=TOL["red"], lw=1.4, label="right", ls="--")
    ax.fill_between(t, 0, 1, where=r["heelF"][:, 0] > 20,
                    color=TOL["blue"], alpha=0.10, transform=ax.get_xaxis_transform())
    ax.fill_between(t, 0, 1, where=r["heelF"][:, 1] > 20,
                    color=TOL["red"], alpha=0.10, transform=ax.get_xaxis_transform())
    if shade:
        ax.axvspan(shade[0], shade[1], color=TOL["yellow"], alpha=0.30, lw=0)
    ax.set_title(title + note, fontsize=10)
    ax.set_xlabel("time (s)"); ax.set_ylabel("knee flexion (deg, +)")
    ax.tick_params(labelsize=9)

knee_panel(axes[0, 0], runs["control"], "(A) control",
           note=f"  {cycles(runs['control']['heelF'], runs['control']['t'])} L cycles")
knee_panel(axes[0, 1], runs["stim_mild"], "(B) RG stimulus, mild",
           shade=(8.0, 9.0))
knee_panel(axes[0, 2], runs["stim_strong"], "(C) RG stimulus, strong",
           shade=(8.0, 9.0))
for ax, key, lab in ((axes[1, 0], "push_fwd", "(D) push forward"),
                     (axes[1, 1], "push_back", "(E) push backward")):
    r = runs[key]
    ax.plot(r["t"], r["com"][:, 0], color=TOL["teal"], lw=1.4, label="COM x (m)")
    ax2 = ax.twinx()
    ax2.plot(r["t"], r["tilt"], color=TOL["purple"], lw=1.2, ls="--",
             label="tilt (deg)")
    ax2.tick_params(labelsize=9)
    ax.axvline(8.0, color=TOL["grey"], lw=1.0, ls=":")
    ax.set_title(lab, fontsize=10)
    ax.set_xlabel("time (s)"); ax.set_ylabel("COM fore-aft (m)")
    ax2.set_ylabel("pelvis tilt (deg)", fontsize=9)
    ax.tick_params(labelsize=9)

axf = axes[1, 2]
axf.axis("off")
rows = [("pathway removed", "L cyc", "R cyc", "fell")]
ctrl_c = cycles(runs["control"]["heelF"], runs["control"]["t"])
rows.append(("control", str(ctrl_c), str(cycles(
    np.nan_to_num(runs["control"]["heelF"]), runs["control"]["t"])), "no"))
for name, d in del_results.items():
    rows.append((name, str(d["cycL"]), str(d["cycR"]),
                 "yes" if d["fell"] else "no"))
tbl = axf.table(cellText=rows[1:], colLabels=rows[0], loc="center",
                cellLoc="center", colWidths=[0.52, 0.14, 0.14, 0.12])
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.35)
axf.set_title("(F) pathway deletions (same walker)", fontsize=10)

OUT = r"D:\GitHub\_zcode_backups\diss_context_20261024\Revision_D\figs"
os.makedirs(OUT, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"w2l_perturbation_panels.{ext}"), dpi=300)
alt = ("Figure: perturbation and pathway-deletion panels for the corrected "
       "W2L walker (per-synapse conductances frozen, 25 percent harness "
       "weight support). A: control left/right knee angles with heel-contact "
       "shading. B and C: one-second inhibitory tonic stimulus to the rhythm "
       "layer, mild and strong (yellow band), with recovery. D and E: "
       "40 newton forward and backward pelvis pushes at 8 seconds (dotted "
       "line) shown as fore-aft center-of-mass position (solid) and pelvis "
       "tilt (dashed). F: strides per leg and fall status when individual "
       "pathway classes are removed: " +
       "; ".join(f"{n}: {d['cycL']}/{d['cycR']} cycles, "
                 f"{'fell' if d['fell'] else 'no fall'}"
                 for n, d in del_results.items()) + ".")
open(os.path.join(OUT, "w2l_perturbation_panels_alt.txt"), "w",
     encoding="utf-8").write(alt)
json.dump(del_results, open(os.path.join(
    HERE, "deletion_results.json"), "w"), indent=1)
print("figure written ->", OUT)
