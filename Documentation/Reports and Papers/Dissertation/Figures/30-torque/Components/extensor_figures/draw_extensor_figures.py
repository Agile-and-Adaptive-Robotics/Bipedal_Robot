"""Draw the six pinned-extensor tests and the design-of-record cam endpoints.

Run export_extensor_data in MATLAB to reconstruct the input arrays. Use the
project's scientific Python with SciPy, Matplotlib, Pillow, and Arial fonts.
PDFs go to the chapter's figure root. SVGs keep their labels editable here.
"""
import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import numpy as np
from PIL import Image, ImageOps
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
CHAPTER = HERE.parents[1]
W = 468
INK = "#161b1e"
GRAY = "#555b5f"
ROUTE = "#213eb3"
FEMUR_CYL = "#a27a00"
FEMUR_WALL = "#d17a16"
CONDYLE = "#ae482a"
TIBIA_LOWER = "#d14b83"
TIBIA_UPPER = "#9a1872"
BEND = "#773ba5"
GOLD = "#efbf16"
MEASURED = "#253fae"
BASE = "#bb4d31"
FIT = "#b12a95"


def setup():
    fonts = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
    for filename in ("arial.ttf", "arialbd.ttf"):
        font_manager.fontManager.addfont(str(fonts / filename))
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 10.5, "axes.linewidth": 0.85,
        "axes.labelcolor": INK, "axes.edgecolor": INK, "xtick.color": INK,
        "ytick.color": INK, "text.color": INK, "pdf.fonttype": 42,
        "svg.fonttype": "none", "svg.hashsalt": "extensor-figures",
        "path.simplify": False, "axes.unicode_minus": True,
    })


def canvas(height):
    return plt.figure(figsize=(W / 72, height / 72), facecolor="white")


def label(fig, x, y, text, **kwargs):
    height = fig.get_figheight() * 72
    fig.text(x / W, 1 - y / height, text, va="center", **kwargs)


def axis(fig, x, top, width, height):
    total = fig.get_figheight() * 72
    return fig.add_axes([x / W, (total - top - height) / total, width / W, height / total])


def ticks(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", width=0.8, length=3, pad=2.5, labelsize=10)


def save(fig, stem, output):
    fig.savefig(output / (stem + ".pdf"), metadata={"Creator": "Ben Bolen", "CreationDate": None,
                                                   "ModDate": None, "Title": stem.replace("_", " ")})
    fig.savefig(HERE / (stem + ".svg"), metadata={"Creator": "Ben Bolen", "Date": None})
    plt.close(fig)


def torque(tests, output):
    fig = canvas(351)
    titles = ["43.6 cm", "43.2 cm + tendon", "45.7 cm", "40.0 cm", "40.0 cm + tendon", "46.5 cm"]
    for k, trial in enumerate(tests):
        row, col = divmod(k, 3)
        left = 40 + 148 * col
        top = 39 + 150 * row
        label(fig, left - 13, top - 30, f"({chr(65 + k)})", weight="bold", size=11)
        label(fig, left + 59, top - 30, titles[k], ha="center", weight="bold", size=11)
        label(fig, left + 59, top - 15, "Training" if row == 0 else "Validation",
              ha="center", size=10.5)
        ax = axis(fig, left, top, 123, 92)
        ax.axhline(0, color="#adb2b4", lw=0.6, zorder=1)
        ax.plot(trial["angle"], trial["unoptimized"], color=BASE, lw=1.4, ls=(0, (4, 2.5)),
                label="Unoptimized", zorder=2)
        ax.plot(trial["angle"], trial["optimized"], color=FIT, lw=1.6,
                label="Optimized", zorder=3)
        ax.plot(trial["hybrid_angle"], trial["hybrid"], ls="none", marker="^", ms=4.5,
                mfc=GOLD, mec="#765b00", mew=0.45, label="Hybrid", zorder=4)
        ax.plot(trial["measured_angle"], trial["measured"], ls="none", marker="o", ms=4.0,
                color=MEASURED, mec="white", mew=0.30, label="Measured", zorder=5)
        ax.set(xlim=(-135, 50), ylim=(-2, 25), xticks=[-120, -80, -40, 0, 40],
               yticks=[0, 5, 10, 15, 20, 25])
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(which="minor", length=1.8, width=0.55)
        ticks(ax)
    label(fig, 9, 159, "Extensor torque (N·m)", rotation=90, ha="center", size=11)
    label(fig, 237, 309, "Knee angle (°)", ha="center", size=11)
    handles, names = ax.get_legend_handles_labels()
    by_name = dict(zip(names, handles))
    names = ["Hybrid", "Measured", "Unoptimized", "Optimized"]
    fig.legend([by_name[n] for n in names], names, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=4, frameon=False, handlelength=2.0,
               columnspacing=2.4, fontsize=10.5)
    save(fig, "pinned_extensor_torque", output)


def circle(center, radius):
    angle = np.linspace(0, 2 * np.pi, 200)
    return center[:2] + radius * np.column_stack([np.cos(angle), np.sin(angle)])


def tibia_to_femur(xy, transform):
    xyz = np.column_stack([xy, np.zeros(len(xy))])
    return (xyz @ transform[:3, :3].T + transform[:3, 3])[:, :2]


def line(ax, xy, color, style="-", width=1.15, zorder=2):
    ax.plot(*np.asarray(xy).T, color=color, ls=style, lw=width, zorder=zorder)


def pose_plot(ax, pose, geo, left):
    points = pose["point"][:, :2]
    active = pose["active"].astype(bool)
    for bone in ("femur", "tibia"):
        ax.fill(*pose[bone].T, color="#e3e5e5", ec="#a5aaac", lw=0.55, zorder=0)
    line(ax, circle(geo["femurCylCenter"], geo["femurCylClearRadius"]), FEMUR_CYL)
    line(ax, np.column_stack([np.repeat(geo["femurLineX"], 2), geo["femurLineY"]]), FEMUR_WALL)
    curve = geo["femurOffsetBoundary"][:, :2]
    line(ax, np.vstack([curve, curve[0]]), CONDYLE)
    line(ax, np.column_stack([np.repeat(geo["femurCondyleClipX"], 2), geo["femurCondyleClipY"]]), CONDYLE)
    transform = pose["tibia_transform"]
    for center, radius, color in (("tibiaLowerCenter", "tibiaLowerClearRadius", TIBIA_LOWER),
                                   ("tibiaUpperCenter", "tibiaUpperClearRadius", TIBIA_UPPER)):
        line(ax, tibia_to_femur(circle(geo[center], geo[radius]), transform), color)
    if active[1]:
        line(ax, [geo["femurCylCenter"], points[1]], BEND, (0, (3, 2)), zorder=3)
    if active[7]:
        center = tibia_to_femur([geo["tibiaLowerCenter"]], transform)[0]
        line(ax, [center, points[7]], BEND, (0, (3, 2)), zorder=3)
    ax.plot(*points.T, ls="-", color=ROUTE, marker="o", mfc=ROUTE, mec=ROUTE,
            lw=1.6, ms=2.7, zorder=4)
    ax.plot(*points[5:8].T, ls="none", marker="o", color=INK, ms=4.2, zorder=5)
    ax.plot(*points[0], marker="s", mfc=GOLD, mec="#765b00", mew=0.55, ms=5.6, zorder=6)
    ax.plot(*points[8], marker="D", mfc="#f0a051", mec="#944b12", mew=0.55, ms=5, zorder=6)

    offsets = {1:(8,1), 2:(9,2), 3:(17,12), 4:(22,1), 5:(19,-10),
               6:(10,-18), 7:(-4,-28), 8:(-24,-15), 9:(-7,14)} if left else {
                   1:(7,2), 7:(10,9), 8:(10,-8), 9:(6,-1)}
    for index, offset in offsets.items():
        assert active[index-1]
        needs_leader = index in (3, 4, 5, 6, 7, 8, 9) and (left or index in (7, 8))
        leaders = dict(arrowstyle="-", lw=0.6, color=GRAY, shrinkA=0.5, shrinkB=2) if needs_leader else None
        ax.annotate(f"p{index}", xy=points[index-1], xytext=offset, textcoords="offset points",
                    ha="left" if offset[0] >= 0 else "right", va="center", size=10,
                    arrowprops=leaders, clip_on=False, zorder=8)
    ax.set(xlim=(-.29, .25), ylim=(-.67, .17), aspect="equal", xticks=[-.2, 0, .2],
           yticks=[-.6, -.4, -.2, 0, .1])
    ticks(ax)
    ax.set_xlabel("Femur-frame x (m)", fontsize=10.2, labelpad=2)
    if left:
        ax.set_ylabel("Femur-frame y (m)", fontsize=10.2, labelpad=1)


def cam(poses, geo, output):
    fig = canvas(372)
    for i, (x, title) in enumerate(((56, "Full flexion (−120°)"), (187, "Extension (+10°)"))):
        label(fig, 1 if i == 0 else x-15, 10, f"({chr(65+i)})", size=11, weight="bold")
        label(fig, x+5, 10, title, size=11, weight="bold")
        ax = axis(fig, x, 27, 128, 205)
        pose_plot(ax, poses[i], geo, i == 0)
    label(fig, 330, 10, "(C)", size=11, weight="bold")
    label(fig, 350, 10, "Test assembly", size=11, weight="bold")
    photo = ImageOps.exif_transpose(Image.open(HERE / "cam_hardware_photo.jpg")).convert("RGB")
    photo.thumbnail((900, 1600), Image.Resampling.LANCZOS)
    pic = axis(fig, 335, 24, 130, 231)
    pic.imshow(photo, interpolation="none")
    pic.set_axis_off()
    entries = [
        ("Femur cylinder", Line2D([], [], color=FEMUR_CYL)),
        ("Femur wall", Line2D([], [], color=FEMUR_WALL)),
        ("Condyle envelope", Line2D([], [], color=CONDYLE)),
        ("Lower tibia", Line2D([], [], color=TIBIA_LOWER)),
        ("Upper tibia", Line2D([], [], color=TIBIA_UPPER)),
        ("Bend radius", Line2D([], [], color=BEND, ls=(0,(3,2)))),
        ("Optimized route", Line2D([], [], color=ROUTE, marker="o", ms=3, lw=1.6)),
        ("Cam contacts", Line2D([], [], ls="none", color=INK, marker="o", ms=4.2)),
        ("Origin p1", Line2D([], [], ls="none", marker="s", mfc=GOLD, mec="#765b00", ms=5.6)),
        ("Insertion p9", Line2D([], [], ls="none", marker="D", mfc="#f0a051", mec="#944b12", ms=5)),
        ("Reference bones", Patch(facecolor="#e3e5e5", edgecolor="#a5aaac", lw=0.55)),
    ]
    legend = fig.legend([item[1] for item in entries], [item[0] for item in entries],
                        loc="lower center", bbox_to_anchor=(.5,.006), ncol=3, title="Legend",
                        frameon=True, fancybox=False, edgecolor="black", facecolor="white",
                        framealpha=1, fontsize=10.5, handlelength=2.4, columnspacing=3.0,
                        labelspacing=.55)
    legend.get_title().set_weight("bold")
    save(fig, "cam_extensor_route", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CHAPTER)
    parser.add_argument("--figure", choices=("both", "torque", "cam"), default="both")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    setup()
    source = loadmat(HERE / "extensor_data.mat", simplify_cells=True)
    if args.figure in ("both", "torque"):
        torque(source["tests"], args.output)
    if args.figure in ("both", "cam"):
        cam(source["poses"], source["geometry"], args.output)

