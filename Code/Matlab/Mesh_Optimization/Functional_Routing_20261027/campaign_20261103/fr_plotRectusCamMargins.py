"""Plot signed torque margins for the two outputs of the rectus cam pilot."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import opensim as osim
from scipy.io import loadmat

from fr_plotPilotAnatomy import mesh_data


def plot(folder, branch):
    d = np.genfromtxt(folder/f"{branch}_pose_grid.csv", delimiter=",", names=True)
    hip, knee = np.unique(d["hip_deg"]), np.unique(d["knee_deg"])
    margins = []
    for joint in ("Knee", "Hip"):
        human = d[f"human{joint}_Nm"]
        margins.append(d[f"robot{joint}_Nm"]*np.sign(human)-abs(human))
    scale = max(np.max(abs(m)) for m in margins)
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4), layout="constrained", sharey=True)
    for ax, margin, label in zip(axes, margins, ("Knee extension", "Hip flexion")):
        mesh = ax.pcolormesh(hip, knee, margin.reshape(knee.size, hip.size), shading="nearest",
                             cmap="RdBu", vmin=-scale, vmax=scale)
        ax.set(title=label, xlabel="Hip flexion angle (°)")
    axes[0].set_ylabel("Knee angle (°) · flexion negative")
    fig.colorbar(mesh, ax=axes, label="Signed torque margin (N·m)")
    fig.suptitle(f"Rectus cam pilot · {branch} · each cell is a coupled pose")
    fig.savefig(folder/f"{branch}_coupled_margins.png", dpi=180)
    plt.close(fig)


def plot_route(folder, branch, data, geometry):
    b = data["bestByBranch"][("rigid", "compliant").index(branch)]
    hip, knee = data["H"].ravel(order="F"), data["K"].ravel(order="F")
    ix = np.argmin(hip**2+knee**2)
    route = data["path"][ix]
    fit = float(data.get("fit", 0.021))  # saved 20 mm extensor context
    direction = route[1]-route[0]
    direction /= np.linalg.norm(direction)
    start = route[0]+fit*direction
    sleeve = np.array([start,start+(b["L0_m"]-b["s_m"][ix])*direction])
    ref = json_model(folder)
    osim.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry))
    model = osim.Model(ref)
    state = model.initSystem()
    coordinates = model.getCoordinateSet()
    coordinates.get("hip_flexion_r").setValue(state, float(np.deg2rad(hip[ix])), False)
    coordinates.get("knee_angle_r").setValue(state, float(np.deg2rad(knee[ix])), False)
    model.realizePosition(state)
    femur = model.getBodySet().get("femur_r")
    tibia = model.getBodySet().get("tibia_r")
    fv, ff = mesh_data(geometry/"femur_r.vtp")
    tv, tf = mesh_data(geometry/"tibia_r.vtp")
    def to_femur(xyz):
        point = tibia.findStationLocationInAnotherFrame(state, osim.Vec3(*xyz), femur)
        return [point.get(i) for i in range(3)]
    tv = np.array([to_femur(v) for v in tv])
    fig, axes = plt.subplots(1,2,figsize=(7.5,5.4),layout="constrained",sharey=True)
    for ax, i, label in zip(axes,(0,2),("Anterior +x (mm)","Right +z (mm)")):
        for v, faces, color in ((fv,ff,"#b9c2ce"),(tv,tf,"#7e8da2")):
            ax.add_collection(PolyCollection(1000*v[faces][:,:,[i,1]],
                                             facecolors=color,edgecolors="none",alpha=.23))
        ax.plot(1000*route[:,i],1000*route[:,1],"-o",color="#1262b2",
                lw=1.5,ms=3,label="Cable and cam contacts")
        ax.plot(1000*sleeve[:,i],1000*sleeve[:,1],color="#dd7517",
                lw=4,label="Center sleeve axis")
        ax.set_ylim(1000*(route[:,1].min()-.055),1000*(route[:,1].max()+.04))
        ax.set_aspect("equal",adjustable="box")
        ax.set_xlabel(label)
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Superior +y (mm)")
    axes[0].set_title("Sagittal view"); axes[1].set_title("Frontal view")
    axes[1].legend(loc="upper center",fontsize=7)
    fig.suptitle(f"Rectus femoris · {branch} cam candidate · packaging failed")
    fig.savefig(folder/f"{branch}_route_anatomical.png",dpi=180)
    plt.close(fig)


def json_model(folder):
    import json
    record = folder.parents[1]/"run_20261105_rectus_sizing"/"rectus_route.json"
    return json.loads(record.read_text())["robot"]


if __name__ == "__main__":
    output = Path(sys.argv[1])
    data = loadmat(output/"rectus_cam_pilot.mat", simplify_cells=True)
    for name in ("rigid", "compliant"):
        plot(output, name)
        plot_route(output,name,data,Path(sys.argv[2]))
