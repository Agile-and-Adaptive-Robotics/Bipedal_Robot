"""Compare the saved gracilis route with native signed human torque."""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import opensim as osim


def human_torque(model, state, muscle, curve, coord, poses, with_lengths=False):
    coordinates = model.getCoordinateSet()
    axis = coordinates.get(coord)
    out = []
    lengths = []
    for q in poses:
        for j in range(coordinates.getSize()):
            item = coordinates.get(j)
            if item.getName() in q:
                item.setValue(state, q[item.getName()], False)
        model.realizePosition(state)
        path_length = muscle.getLength(state)
        length = (path_length - muscle.getTendonSlackLength()) / muscle.getOptimalFiberLength()
        active = curve.calcValue(float(length))
        out.append(muscle.getMaxIsometricForce()*active*muscle.computeMomentArm(state, axis))
        if with_lengths:
            lengths.append(path_length)
    return (np.asarray(out), np.asarray(lengths)) if with_lengths else np.asarray(out)


def plot(folder):
    data = json.loads((folder/"gracilis_two_dof_sweeps.json").read_text())
    model = osim.Model(data["humanModel"])
    state = model.initSystem()
    muscle = model.getMuscles().get(data["muscle"])
    curve = osim.ActiveForceLengthCurve()
    labels = ("Knee angle (°) · flexion negative", "Hip angle (°) · adduction positive")
    titles = ("Knee flexion / extension", "Hip abduction / adduction")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.9), layout="constrained")
    for ax, sweep, label, title in zip(axes, data["sweeps"], labels, titles):
        x = np.asarray(sweep["angles_deg"])
        robot = np.asarray(sweep["robot_Nm"], dtype=float)
        human = human_torque(model, state, muscle, curve, sweep["coordinate"], sweep["poses"])
        slack = np.asarray(sweep["slack"], dtype=bool)
        invalid = np.asarray(sweep["invalid"], dtype=bool)
        meets = ~invalid & (robot*np.sign(human) >= np.abs(human)-1e-6)
        ax.axhline(0, color="#adb5bd", lw=.8)
        ax.plot(x, human, "--", color="#505763", lw=2, label="Human gracilis")
        ax.plot(x, robot, color="#1262b2", lw=2, label="Robot candidate")
        if slack.any():
            ax.scatter(x[slack], robot[slack], color="#ba3e22", s=14, label="Slack")
        if invalid.any():
            ax.scatter(x[invalid], np.zeros(invalid.sum()), marker="x", color="#ba3e22",
                       s=18, label="Stroke infeasible")
        ax.set(title=title, xlabel=label, ylabel="Signed torque (N·m)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
        with (folder/f"gracilis_{sweep['coordinate']}.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["angle_deg", "human_Nm", "robot_Nm", "slack", "invalid", "meets"])
            writer.writerows(zip(x, human, robot, slack, invalid, meets))
        print(f"{title}: {meets.sum()}/{len(x)} met; {slack.sum()} slack; {invalid.sum()} invalid")
    fig.suptitle("Gracilis · saved compliant route · independent joint sweeps")
    fig.savefig(folder/"gracilis_two_dof_torque.png", dpi=200)
    plot_coupled(folder, data["coupled"], model, state, muscle, curve)


def plot_coupled(folder, grid, model, state, muscle, curve):
    hip, knee = np.asarray(grid["hip_deg"]), np.asarray(grid["knee_deg"])
    shape = (knee.size, hip.size)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), layout="constrained", sharey=True)
    margins, curves = [], []
    for name, coordinate in (("knee", "knee_angle_r"), ("hip", "hip_adduction_r")):
        robot = np.asarray(grid[f"{name}_Nm"], dtype=float)
        human = human_torque(model, state, muscle, curve, coordinate, grid["poses"])
        margins.append(robot*np.sign(human)-np.abs(human))
        curves.extend([human, robot])
    scale = max(np.nanmax(np.abs(x)) for x in margins)
    for ax, margin, title in zip(axes, margins, ("Knee flexion", "Hip adduction")):
        mesh = ax.pcolormesh(hip, knee, margin.reshape(shape), cmap="RdBu",
                             shading="nearest", vmin=-scale, vmax=scale)
        ax.set(title=title, xlabel="Hip angle (°) · adduction positive")
        print(f"Coupled grid {title}: {np.count_nonzero(margin >= -1e-6)}/{margin.size} met")
    both = (margins[0] >= -1e-6) & (margins[1] >= -1e-6)
    print(f"Coupled grid both: {both.sum()}/{both.size}; "
          f"{sum(grid['slack'])} slack; {sum(grid['invalid'])} invalid")
    with (folder/"gracilis_coupled_torque.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["hip_deg", "knee_deg", "human_knee_Nm", "robot_knee_Nm",
                         "human_hip_Nm", "robot_hip_Nm", "slack", "invalid", "both_meet"])
        writer.writerows(zip(np.tile(hip,knee.size), np.repeat(knee,hip.size), *curves,
                             grid["slack"], grid["invalid"], both))
    axes[0].set_ylabel("Knee angle (°) · flexion negative")
    fig.colorbar(mesh, ax=axes, label="Signed torque margin (N·m)")
    fig.suptitle("Gracilis · one BiPulley equilibrium per combined joint pose")
    fig.savefig(folder/"gracilis_coupled_torque_margins.png", dpi=200)


if __name__ == "__main__":
    plot(Path(sys.argv[1]))
