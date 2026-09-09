"""Plot the saved phase-1 recordings without rerunning or changing the model."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").exists())
DISS = REPO / "Documentation/Reports and Papers/Dissertation"
SOURCE = REPO / "Neuromechanical_Models/Biped_2xCPG_wSubs"
OUT = DISS / "ProofFinal/figs/Preliminary"
OUT.mkdir(parents=True, exist_ok=True)

def read_chart(name):
    path = SOURCE / name
    with path.open() as stream:
        names = stream.readline().strip().split("\t")
    raw = np.loadtxt(path, skiprows=1)
    assert raw.shape == (50010, len(names))
    assert np.isfinite(raw).all()
    assert np.allclose(np.diff(raw[:, 1]), 0.0002, atol=1e-9, rtol=0)
    keep = (raw[:, 1] >= 0) & (raw[:, 1] <= 10)
    assert keep.sum() == 50001
    assert np.all(raw[~keep, 2:] == 0), "Unexpected nonzero data beyond chart end"
    return names, raw[keep], {
        "path": str(path.relative_to(REPO)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "total_rows": len(raw), "retained_rows": int(keep.sum()),
        "excluded_padding_rows": int((~keep).sum()),
    }

n7, d7, info7 = read_chart("DataTool_7.txt")
n8, d8, info8 = read_chart("DataTool_8.txt")
assert np.array_equal(d7[:, :2], d8[:, :2])
model_path = SOURCE / "walk new new tester added 2 axis_phase1.asim"
model = ET.parse(model_path).getroot()
columns = {e.findtext("ColumnName"): e.findtext("DataType")
           for e in model.iter("DataColumn")}
for name in ["knee_L_deg", "knee_R_deg"]:
    assert columns[name] == "JointRotationDeg"
for name in ["RG_L_stance", "RG_R_stance", "L_foot ground contact", "R_foot ground contact"]:
    assert columns[name] == "MembraneVoltage"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "axes.linewidth": 0.7})
fig, axes = plt.subplots(3, 1, figsize=(6.5, 6.1), sharex=True,
                         layout="constrained")
blue, orange = "#17659b", "#b85314"
panels = [
    (n8, d8, ["knee_L_deg", "knee_R_deg"], 1,
     "A  Knee joint recordings", "Joint angle (deg)"),
    (n8, d8, ["RG_L_stance", "RG_R_stance"], 1000,
     "B  Stance rhythm-generator recordings", "Membrane voltage (mV)"),
    (n7, d7, ["L_foot ground contact", "R_foot ground contact"], 1000,
     "C  Foot-contact sensory-neuron recordings", "Membrane voltage (mV)"),
]
for ax, (names, data, channels, scale, title, ylabel) in zip(axes, panels):
    for channel, label, color, style in zip(channels, ["Left", "Right"],
                                           [blue, orange], ["-", "--"]):
        ax.plot(data[:, 1], scale * data[:, names.index(channel)],
                color=color, ls=style, lw=0.7, label=label)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.set_ylabel(ylabel)
    low, high = ax.get_ylim()
    ax.set_ylim(low, high + 0.22 * (high - low))
    ax.grid(axis="y", color="#dddddd", lw=0.5)
    ax.legend(loc="upper right", ncol=2, frameon=True, fontsize=8,
              framealpha=0.95, edgecolor="none")
    ax.set_xlim(0, 10)
axes[-1].set_xlabel("Time (s)")
axes[-1].set_xticks(np.arange(0, 11, 2))
fig.suptitle("AnimatLab phase 1: preliminary signal recording", fontsize=11)
pdf = OUT / "animatlab_phase1_preliminary.pdf"
fig.savefig(pdf)
plt.close(fig)

# Main dissertation result: bilateral hip, knee, and ankle recordings.
# Retain the original neural/contact diagnostic figure as a supporting asset.
fig, axes = plt.subplots(3, 1, figsize=(6.5, 5.6), sharex=True, layout="constrained")
for ax, joint, letter in zip(axes, ["hip", "knee", "ankle"], "ABC"):
    for side, label, color, style in zip(["L", "R"], ["Left", "Right"],
                                       [blue, orange], ["-", "--"]):
        channel = f"{joint}_{side}_deg"
        assert columns[channel] == "JointRotationDeg"
        ax.plot(d8[:, 1], d8[:, n8.index(channel)], color=color, ls=style,
                lw=0.9, label=label)
    ax.set_title(f"{letter}  {joint.capitalize()}", loc="left", fontweight="bold", fontsize=10)
    ax.set_ylabel("Joint angle (deg)")
    low, high = ax.get_ylim()
    ax.set_ylim(low, high + 0.22 * (high - low))
    ax.legend(loc="upper right", ncol=2, frameon=False, fontsize=8)
    ax.grid(axis="y", color="#dddddd", lw=0.5)
    ax.set_xlim(0, 10)
axes[-1].set_xlabel("Time (s)")
axes[-1].set_xticks(np.arange(0, 11, 2))
fig.suptitle("AnimatLab phase 1: preliminary joint motion", fontsize=11)
fig.savefig(OUT / "animatlab_phase1_joint_angles.pdf")
plt.close(fig)
report = {
    "source_model": str(model_path.relative_to(REPO)),
    "source_model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
    "charts": [info7, info8], "retained_time_seconds": [0, 10],
    "sample_interval_seconds": 0.0002,
    "checks": ["finite values", "uniform increasing times", "aligned chart times",
               "channel types checked against source model", "zero padding excluded only after chart end"],
    "plotted_channels": [c for panel in panels for c in panel[2]],
    "joint_figure_channels": [f"{joint}_{side}_deg" for joint in ["hip", "knee", "ankle"] for side in ["L", "R"]],
    "processing": "No filtering or resampling; volts multiplied by 1000 for mV. Within-window zero-voltage samples retained.",
    "limits": ["Contact channels are sensory-neuron voltages, not forces or Boolean contact states.",
               "Model includes an enabled external ForceX=40 input to Root during 0-1 s.",
               "No claim of stable walking or effective ground-contact gating.",
               "Body-position channels were recorded but are not plotted; their unit convention was not resolved."],
}
(DISS / "Notes/preliminary_animatlab_verification.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"pdf": str(pdf), "retained_samples_per_chart": 50001,
                  "padding_rows_excluded_per_chart": 9}, indent=2))
