"""Draw the solved sleeve and hardware envelopes over the bone surfaces."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image

from fr_plotAnkleClash import capsule
from fr_plotAssembledAnkle import load_meshes, triangle_tube
from fr_plotDistalStages import GOLD, INDIGO
from fr_wholeLegPack import read, unpack


def draw(folder, options, selected, poses, verification, before_options=None,
         headings=("Starting layout", "Updated layout")):
    scene = dict(anatomy=poses["anatomy"],
                 transforms={k: v[0] for k, v in poses["transforms"].items()})
    bones = load_meshes(scene, Path("C:/OpenSim 4.6/Geometry"))
    clear = {r["group"]: r["clearance_pass"] for r in verification}
    plt.rcParams.update({"font.family": "Arial", "font.size": 10})
    fig, axes = plt.subplots(1, 4, figsize=(12, 7.8), layout="constrained", sharey=True)
    if before_options is None:
        before_options = options
    for index, title in enumerate(headings):
        for ax, dim, name in zip(axes[2*index:2*index+2], (0, 2), ("Sagittal", "Posterior")):
            for tri in bones:
                ax.add_collection(PolyCollection(tri[:, :, [dim, 1]], facecolor="#bdc5d0",
                                                 edgecolor="none", alpha=.22))
            for g, data in (options if index else before_options).items():
                option = data[selected[g] if index else 1]
                color = INDIGO if not index or clear[g] else GOLD
                for row in option["parts"][0]:
                    p, q, radius = 1000*np.array(row[:3]), 1000*np.array(row[3:6]), 1000*row[6]
                    a, b = p[[dim, 1]], q[[dim, 1]]
                    if row[7] == 0:
                        capsule(ax, a, b, radius, color, .34)
                        ax.plot(*np.array([a, b]).T, color=color, alpha=.58, lw=.65)
                    elif row[7] == 1:
                        ax.plot(*np.array([a, b]).T, color="#5a626c", marker=".", alpha=.5, lw=1)
                    else:
                        ax.plot(*np.array([a, b]).T, color="#444a52", alpha=.33, lw=.5)
            ax.set(title=f"{title}\n{name}", aspect="equal",
                   xlabel=("Anterior +x (mm)" if dim == 0 else "Right +z (mm)"),
                   ylim=(-390, 930))
            ax.autoscale_view(scaley=False)
            ax.grid(alpha=.12)
    # Ground coordinates include the model's default pelvis height.
    y = np.concatenate([t[:, :, 1].ravel() for t in bones])
    vessel_y = [1000*r[1] for g, opts in options.items() for r in opts[selected[g]]["parts"][0]]
    axes[0].set_ylim(min(y.min(), min(vessel_y))-25, max(vessel_y)+55)
    axes[0].set_ylabel("Superior +y (mm)")
    fig.legend(handles=[Line2D([], [], lw=4, color=c, label=t) for c, t in (
        (INDIGO, "Starting layout or clears sampled checks"), (GOLD, "Clearance failure in updated layout"),
        ("#5a626c", "Hardware and cable centerlines"))],
        loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Whole-leg packaging at the reference pose · compliant branch", fontweight="bold")
    fig.savefig(folder / "whole_leg_packaging.png", dpi=190)
    plt.close(fig)
    (folder / "whole_leg_packaging.alt.txt").write_text(
        "Sagittal and posterior overlays compare the original and coordinated right-leg "
        "sleeve envelopes at the reference pose. Gold marks groups with a clearance "
        "failure at some sampled pose. The clearance tables give the governing margins.\n",
        encoding="utf-8")
    orbit(folder, options, selected, bones, clear)


def orbit(folder, options, selected, bones, clear):
    fig = plt.figure(figsize=(6.3, 7.1), dpi=92, facecolor="white")
    ax = fig.add_subplot(projection="3d", proj_type="ortho")
    allpoints = []
    for tris in bones:
        ax.add_collection3d(Poly3DCollection(tris[:, :, [0, 2, 1]], facecolor="#bac4ce",
                                            edgecolor="none", alpha=.27))
    for g, opts in options.items():
        rows = opts[selected[g]]["parts"][0]
        color = INDIGO if clear[g] else GOLD
        quads = []
        for row in rows:
            p, q = 1000*np.asarray(row[:3])[[0, 2, 1]], 1000*np.asarray(row[3:6])[[0, 2, 1]]
            allpoints.extend((p, q))
            if row[7] == 0:
                quads.extend(triangle_tube(p, q, 1000*row[6]))
            elif row[7] == 1:
                ax.plot(*np.array([p, q]).T, color="#454b55", marker=".", lw=1)
            else:
                ax.plot(*np.array([p, q]).T, color="#555a64", lw=.65, alpha=.45)
        if quads:
            ax.add_collection3d(Poly3DCollection(quads, facecolor=color, edgecolor="none", alpha=.62))
    points = np.asarray(allpoints)
    lower, upper = points.min(axis=0)-44, points.max(axis=0)+44
    ax.set(xlim=(lower[0], upper[0]), ylim=(lower[1], upper[1]), zlim=(lower[2], upper[2]),
           box_aspect=upper-lower)
    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=.96)
    fig.suptitle("Inflated sleeves · compliant packaging candidate", fontsize=12)
    frames = []
    for azimuth in np.linspace(-60, 120, 19):
        ax.view_init(elev=8, azim=azimuth)
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert("RGB"))
    frames[0].save(folder / "whole_leg_180deg_orbit.gif", save_all=True,
                   append_images=frames[1:], duration=115, loop=0, optimize=True)
    plt.close(fig)


if __name__ == "__main__":
    import csv
    folder = Path(sys.argv[1])
    options = unpack(folder / "dense", "compliant")
    selection = json.loads((folder / "selection.json").read_text())
    with (folder / "compliant_packed_groups.csv").open() as stream:
        verification = list(csv.DictReader(stream))
    for row in verification:
        row["clearance_pass"] = row["clearance_pass"] == "True"
    draw(folder, options, {g: selection[g+"_compliant"] for g in options},
         read(folder / "dense" / "poses.json"), verification)
