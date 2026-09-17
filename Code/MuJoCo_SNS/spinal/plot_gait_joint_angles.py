"""Phase-normalized OpenSim IK joint angles for interpreting PF synergies.

Rows follow the requested trunk/hip/knee/ankle/MTP organization. The left
column contains sagittal coordinates. The right contains the available
OpenSim YZ-plane coordinates; structurally absent knee and MTP panels remain
blank rather than implying a zero-valued coordinate.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gait_phase import IK_MOT, PHASE, phase_normalize, read_storage


HERE = Path(__file__).parent
OUT = HERE / "fsa_results"
COLORS = {"r": "#0072B2", "l": "#D55E00"}

SAGITTAL = (
    ("Trunk", "lumbar_extension"),
    ("Hip", "hip_flexion_{side}"),
    ("Knee", "knee_angle_{side}"),
    ("Ankle", "ankle_angle_{side}"),
    ("MTP", "mtp_angle_{side}"),
)
YZ_PLANE = (
    ("Trunk", "lumbar_bending"),
    ("Hip", "hip_adduction_{side}"),
    ("Knee", None),
    ("Subtalar", "subtalar_angle_{side}"),
    ("MTP", None),
)


def decorate_phase_axis(ax):
    ax.axvspan(0, 50, color="0.92", zorder=-10)
    ax.axvline(50, color="0.35", ls="--", lw=0.8)
    ax.axhline(0, color="0.75", lw=0.6)
    ax.set_xlim(0, 100)
    ax.grid(alpha=0.20)


def main():
    OUT.mkdir(exist_ok=True)
    time, names, values = read_storage(IK_MOT)
    col = {name: index for index, name in enumerate(names)}
    phased = {
        side: phase_normalize(time, values, side)
        for side in ("r", "l")
    }

    fig, axes = plt.subplots(
        5, 2, figsize=(10.5, 12.5), sharex="col", constrained_layout=True)
    axes[0, 0].set_title("Sagittal-plane OpenSim coordinates", fontsize=11)
    axes[0, 1].set_title("OpenSim YZ-plane coordinates", fontsize=11)

    for column, specification in enumerate((SAGITTAL, YZ_PLANE)):
        for row, (region, template) in enumerate(specification):
            ax = axes[row, column]
            if template is None:
                ax.axis("off")
                ax.text(0.5, 0.5, f"{region}: no modeled coordinate",
                        transform=ax.transAxes, ha="center", va="center",
                        color="0.45", fontsize=9, style="italic")
                continue
            decorate_phase_axis(ax)
            for side in ("r", "l"):
                coordinate = template.format(side=side)
                index = col[coordinate]
                result = phased[side]
                mean = result["mean"][:, index]
                std = result["std"][:, index]
                label = (f"{'right' if side == 'r' else 'left'} "
                         f"(n={len(result['cycles'])})")
                ax.plot(PHASE, mean, color=COLORS[side], lw=1.5, label=label)
                if len(result["cycles"]) > 1:
                    ax.fill_between(PHASE, mean - std, mean + std,
                                    color=COLORS[side], alpha=0.15,
                                    linewidth=0)
            ax.set_ylabel(f"{region}\nangle (deg)", fontsize=8)
            ax.text(0.99, 0.94, template.replace("_{side}", ""),
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=6.5, color="0.38")
    for ax in axes[-1, :]:
        if ax.axison:
            ax.set_xlabel("stance-rescaled gait phase (%)\n"
                          "0–50 stance; 50–100 swing")
    axes[0, 0].legend(loc="best", fontsize=7)
    fig.suptitle(
        "OpenSim IK joint kinematics aligned to each ipsilateral gait cycle\n"
        "heel strike = 0%, toe-off = 50%, next heel strike = 100%",
        fontsize=12)
    fig.savefig(OUT / "gait_joint_angles_phase.png", dpi=220)
    fig.savefig(OUT / "gait_joint_angles_phase.pdf")
    plt.close(fig)

    arrays = {"phase": PHASE}
    metadata = {
        "source": str(IK_MOT),
        "phase_convention":
            "heel strike=0, toe off=50, next heel strike=100",
        "coordinates": {
            "sagittal": [value for _, value in SAGITTAL],
            "openSim_YZ_plane": [value for _, value in YZ_PLANE],
        },
        "sides": {},
    }
    for side, result in phased.items():
        arrays[f"{side}_cycles"] = result["cycles"]
        arrays[f"{side}_mean"] = result["mean"]
        arrays[f"{side}_std"] = result["std"]
        metadata["sides"][side] = {
            "n_complete_cycles": int(len(result["cycles"])),
            "events_s": result["events"],
            "measured_duty": result["duty"].tolist(),
        }
    arrays["coordinate_names"] = np.asarray(names)
    np.savez_compressed(OUT / "gait_joint_angles_phase.npz", **arrays)
    (OUT / "gait_joint_angles_phase.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata["sides"], indent=2))
    print(f"wrote {OUT / 'gait_joint_angles_phase.png'}")


if __name__ == "__main__":
    main()

