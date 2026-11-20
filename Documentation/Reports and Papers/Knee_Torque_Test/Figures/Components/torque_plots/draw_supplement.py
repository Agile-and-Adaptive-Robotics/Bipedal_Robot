"""Rebuild the two supplementary figures from the torque paper's saved designs."""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.io import loadmat


FIGURES = Path(__file__).resolve().parents[2]
ROOT = next(folder for folder in FIGURES.parents if (folder / "Code/Matlab").is_dir())
RESULTS = ROOT / "Code/Matlab/Mesh_Optimization/Results"
SURVEY = RESULTS / "KompareXXVII"
sys.path.insert(0, str(SURVEY / "Slide57"))
import route_visual as slides


def knee_routes():
    designs = [
        "Vas_Pam_20mm_Result_20260925.mat",
        "Vas_Pam_20mm_Result_cam_20261014_0416.mat",
        "Bifemsh_20mm_Result_20260921_2037.mat",
    ]
    extensor, cam, flexor = [
        loadmat(RESULTS / name, variable_names=["ctx", "predBest"], simplify_cells=True)
        for name in designs
    ]
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.95))
    fig.subplots_adjust(left=.13, right=.985, bottom=.40, top=.89, wspace=.37)
    for ax, direct, name, limits in zip(axes, [extensor, flexor], ["Extensor", "Flexor"],
                                        [(0, 120), (-90, 0)]):
        ctx = direct["ctx"]
        theta = ctx["phiD"]
        human = ctx["humanTorque"]
        ax.plot(theta, direct["predBest"]["TorqueZ"], color="#e27b69", lw=2.1, label="Direct route")
        if direct is extensor:
            ax.plot(cam["ctx"]["phiD"], cam["predBest"]["TorqueZ"], color="#b92ba1",
                    lw=2.2, ls="--", label="Cam extensor")
        ax.plot(theta, 1.05 * human, color="#666666", ls=":", lw=2.1, label="Required (+5%)")
        ax.plot(theta, human, color="black", ls="-.", lw=1.7, label="Human reference")
        ax.set(xlim=(-120, 10), ylim=limits, xticks=[-120, -90, -60, -30, 0, 10])
        ax.set_title(name, loc="center", fontsize=11, weight="bold", pad=8)
        ax.set_xlabel(r"$\mathbf{\theta_k}$ (°)", weight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=1.6, length=4)
    axes[0].set_ylabel("Torque (N m)", weight="bold")
    for letter, left in [("(A)", .004), ("(B)", .526)]:
        fig.text(left, .927, letter, weight="bold", fontsize=11)
    legend = fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=2,
                        bbox_to_anchor=(.5, .007), title="Legend", frameon=True, fancybox=False,
                        facecolor="white", edgecolor="black", framealpha=1, handlelength=3, columnspacing=2)
    legend.get_title().set_weight("bold")
    fig.savefig(FIGURES / "knee_route_torque.pdf")
    plt.close(fig)


def survey_routes(geometry):
    with (SURVEY / "route_configurations.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if row["side"] == "selected"]
    model = loadmat(slides.HERE / "route_projection.mat", squeeze_me=True, struct_as_record=False)
    assert [r["site"] for r in rows] == [s.site for s in model["sel"]]
    assert len(rows) == 27 and sum(int(row["nBPA"]) for row in rows) == 69
    assert all(row["points_identical_to_baseline"] == "true" for row in rows)
    neutral = json.loads(slides.GEOMETRY.read_text(encoding="utf-8"))["bodyTransformsInGround"]["neutral"]
    points = {row["site"]: rec.P for row, rec in zip(rows, model["Pel"])}
    points["psoas_r"] = slides.connected_psoas(points["psoas_r"], points["iliacus_r"],
                                               slides.relative_frame(neutral, "torso", "pelvis"))
    fig = plt.figure(figsize=(6.5, 5.1))
    left, right = [fig.add_axes(bounds) for bounds in [[.035, .18, .42, .74], [.49, .18, .47, .74]]]
    left.imshow(plt.imread(slides.HERE / "robot_routes_2019.png"))
    slides.draw_skeleton(right, model, neutral, geometry)
    for mesh in right.collections:
        mesh.set_rasterized(True)
    styles = {"single": (slides.COLORS["single"], "-."), "bundle": (slides.COLORS["bundle"], "-"),
              "tackle": (slides.COLORS["tackle"], "--")}
    for row in rows:
        color, style = styles[slides.group(row)]
        xy = slides.project(points[row["site"]])
        right.plot(*xy.T, color=color, linestyle=style, linewidth=1.3)
    right.set_aspect("equal")
    for ax in [left, right]:
        ax.axis("off")
    for x, label in [(.045, "(A)  Inherited routes"), (.495, "(B)  Selected transmissions")]:
        fig.text(x, .96, label, fontsize=11, weight="bold", va="top")
    labels = ["One BPA, direct", "Parallel BPAs, direct", "Parallel BPAs, tackle"]
    handles = [Line2D([], [], color=c, ls=s, lw=2) for c, s in styles.values()]
    legend = fig.legend(handles, labels, loc="lower center", ncol=3, title="Legend", frameon=True,
                        fancybox=False, facecolor="white", edgecolor="black", framealpha=1,
                        handlelength=2.3, columnspacing=1., handletextpad=.4, bbox_to_anchor=(.5, .008))
    legend.get_title().set_weight("bold")
    fig.savefig(FIGURES / "survey_routes.pdf", dpi=320)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=Path("C:/OpenSim 4.6/Geometry"))
    args = parser.parse_args()
    plt.rcParams.update({"font.family":"Arial", "font.size":10.5, "axes.linewidth":1.6,
                         "pdf.fonttype":42, "mathtext.fontset":"custom", "mathtext.rm":"Arial",
                         "mathtext.it":"Arial:italic", "mathtext.bf":"Arial:bold"})
    knee_routes()
    survey_routes(args.geometry)
