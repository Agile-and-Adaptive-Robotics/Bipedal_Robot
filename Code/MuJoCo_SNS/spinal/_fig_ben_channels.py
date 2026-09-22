"""Full-channel figure (Ben 2026-09-21): per-side RG / PF / MN activity
and joint angles for pelvis (wrt ground), hip, knee, ankle, MTP -
including MTP + lumbar from the raw qpos log (qfull), which the named
KEY_JOINTS channels omit.

Columns: right leg | left leg. Rows: RG half-centers, PF cells, MN
activations (swing vs stance groups), joint angles (hip, knee, ankle,
MTP; pelvis tilt wrt ground on both columns).
Standards: AGENTS.md figure rules (Arial 10pt, Tol palette, no italics).
"""
import io
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import mujoco

import runner as R

HERE = Path(__file__).parent
FIGS = HERE / "figures"
DISS = (HERE.parents[2] / "Documentation" / "Reports and Papers"
        / "Dissertation" / "CPG_airstepping_figs")

z = np.load("spinal_run.npz", allow_pickle=True)
t, qf, neuro, act = z["t"], z["qfull"], z["neuro"], z["act"]
names = [str(s) for s in z["neuro_names"]]
acts = [str(s) for s in z["key_acts"]]
z.close()

model = mujoco.MjModel.from_xml_path(str(R.MODEL))


def jadr(name):
    j = model.joint(name)
    return j.qposadr[0] if j.id >= 0 else None


J = {nm: jadr(nm) for nm in
     ("pelvis_tilt", "lumbar_extension",
      "hip_flexion_r", "hip_flexion_l",
      "knee_angle_r", "knee_angle_l",
      "ankle_angle_r", "ankle_angle_l",
      "mtp_angle_r", "mtp_angle_l")}


def jr(name, m=t >= 2.0):
    i = J[name]
    return np.degrees(qf[m, i])


m = t >= 2.0
plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                     "font.style": "normal"})
fig, ax = plt.subplots(4, 2, figsize=(9.5, 10.5), sharex=True)
fig.suptitle("Per-side neural activity and joint angles "
             "(walk window; RG/PF right and left)", fontsize=12,
             fontweight="bold")
COL = {"r": "#d62728", "l": "#1f77b4"}

for col, side in enumerate(("r", "l")):
    a = ax[0, col]
    for ch, c, ls in ((f"RG_E_{side}", "#9D02D7", "-"),
                      (f"RG_F_{side}", "#0072B2", "-")):
        a.plot(t[m], neuro[m, names.index(ch)], color=c, ls=ls, lw=1.2,
               label=ch)
    a.set_ylabel("RG [mV]", fontsize=10)
    a.legend(fontsize=8, ncol=2, loc="upper right")
    a.set_title(f"{'Right' if side == 'r' else 'Left'} side",
                fontsize=11, fontweight="bold")

    a = ax[1, col]
    for ch, c in ((f"PF_E1_{side}", "#FFB14E"),
                  (f"PF_F1_{side}", "#EA5F94"),
                  (f"PF_E2_{side}", "#CD34B5"),
                  (f"PF_F2_{side}", "#0000FF")):
        if ch in names:
            a.plot(t[m], neuro[m, names.index(ch)], color=c, lw=1.0,
                   label=ch)
    a.set_ylabel("PF [mV]", fontsize=10)
    a.legend(fontsize=7, ncol=2, loc="upper right")

    a = ax[2, col]
    for base, c in (("vas_lat", "#FFB14E"), ("soleus", "#CD34B5"),
                    ("semimem", "#0072B2"), ("tib_ant", "#56B4E9"),
                    ("psoas", "#009E73"), ("glut_max2", "#CC79A7")):
        key = f"{base}_{side}"
        if key in acts:
            a.plot(t[m], act[m, acts.index(key)], color=c, lw=1.0,
                   label=base)
    a.set_ylabel("MN activation", fontsize=10)
    a.legend(fontsize=7, ncol=3, loc="upper right")

    a = ax[3, col]
    a.plot(t[m], jr("pelvis_tilt"), color="#555555", lw=1.2,
           label="pelvis tilt (wrt gnd)")
    for jn, c in ((f"hip_flexion_{side}", "#d62728"),
                  (f"knee_angle_{side}", "#0072B2"),
                  (f"ankle_angle_{side}", "#009E73"),
                  (f"mtp_angle_{side}", "#9D02D7")):
        a.plot(t[m], jr(jn), color=c, lw=1.1, label=jn.split("_")[0])
    a.set_ylabel("angle [deg]", fontsize=10)
    a.legend(fontsize=7, ncol=3, loc="upper right")
    a.set_xlabel("t [s]", fontsize=10)

fig.tight_layout(rect=(0, 0.005, 1, 0.97))
fig.savefig(FIGS / "curr3i_channels.png", dpi=200)
fig.savefig(DISS / "curr3i_channels.png", dpi=200)
plt.close(fig)

ALT = """Alt text - curr3i_channels figure. Four rows, two columns
(right side, left side) over the walk window. Row 1: rhythm-generator
half-center potentials RG-E and RG-F per side. Row 2: the four
pattern-formation phase-window cells per side. Row 3: motoneuron
activations of six representative muscles per side (knee extensor,
ankle plantarflexor, knee/hip flexor, ankle dorsiflexor, hip flexor,
hip extensor). Row 4: joint angles - pelvis tilt relative to ground
(common to both columns), hip flexion, knee angle (negative =
flexion), ankle angle (positive = dorsiflexion), and MTP angle, in
degrees. The right side shows rhythmic activity and joint motion; the
left side shows tonic activity with its joints near-static (the
frozen-leg limitation). Axes follow the OpenSim/ISB conventions."""
(FIGS / "curr3i_channels.alt.txt").write_text(ALT, encoding="utf-8")
(DISS / "curr3i_channels.alt.txt").write_text(ALT, encoding="utf-8")
print("saved curr3i_channels.png + alt (figures + Dissertation)")
