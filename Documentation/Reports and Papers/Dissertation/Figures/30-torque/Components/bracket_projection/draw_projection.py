"""Draw Figure 4.5 in its original point coordinates. Displacements are schematic."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

BLUE, ROSE, GREEN = "#4477AA", "#EE6677", "#228833"
DISPLACEMENT, PROJECTION, TENDON = "#882255", "#332288", "#CC6677"
GUIDE = "#8C8C8C"
SIZE = (416.592, 500.4)


def arrow(ax, start, end, color, width, style="-|>", scale=11):
    patch = FancyArrowPatch(start, end, arrowstyle=style, lw=width, color=color,
                            mutation_scale=scale, shrinkA=0, shrinkB=0,
                            capstyle="round", joinstyle="round", zorder=3)
    ax.add_patch(patch)


def text(ax, point, label, color="black", size=11.5):
    ax.text(*point, label, color=color, fontsize=size, ha="left", va="baseline")


def frame(ax, index, origin, ends, labels):
    for axis, (end, label, color) in enumerate(zip(ends, labels, (BLUE, ROSE, GREEN)), 1):
        arrow(ax, origin, end, color, 1.9)
        text(ax, label, rf"$\hat{{\mathbf{{e}}}}_{{{axis}}}^{{({index})}}$", color, 10.5)


def bracket(ax, index, origin, direction, displacement, offset, labels):
    # Both arrows use the same axial length. Only the dimension line is offset.
    axial = np.dot(displacement, direction)
    foot = origin + axial * direction
    tip = origin + displacement
    start = origin + offset
    end = foot + offset
    ax.plot([tip[0], end[0]], [tip[1], end[1]], color=GUIDE, lw=1.1,
            dashes=(3, 2.2), zorder=1)
    arrow(ax, origin, tip, DISPLACEMENT, 2.7, scale=12)
    arrow(ax, start, end, PROJECTION, 2.6, style="<->", scale=11)
    text(ax, labels[0], rf"$\mathbf{{\eta}}_{{br,{index}}}$", DISPLACEMENT)
    text(ax, labels[1], rf"$\delta_{{br,{index}}}$", PROJECTION)
    return end


def draw(destination):
    with plt.rc_context({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                         "pdf.fonttype": 42, "figure.facecolor": "white"}):
        fig = plt.figure(figsize=np.array(SIZE) / 72)
        ax = fig.add_axes([0, 0, 1, 1], xlim=(0, SIZE[0]), ylim=(0, SIZE[1]), aspect="equal")
        ax.axis("off")
        lower, upper = np.array([107.028, 72.108]), np.array([292.104, 393.372])
        u = (upper - lower) / np.linalg.norm(upper - lower)
        normal = np.array([u[1], -u[0]])
        ax.plot([lower[0], upper[0]], [lower[1], upper[1]], lw=1, dashes=(4, 3), color=GUIDE)

        frame(ax, 1, lower, [(149.84, 63.15), (82.48, 43.64), (107.03, 121.37)],
              [(153, 57), (65.3, 30.7), (98, 130.4)])
        frame(ax, 2, upper, [(324.15, 416.05), (261.18, 403.04), (292.10, 442.64)],
              [(325.4, 419.8), (240.2, 403.3), (283.1, 451.6)])
        text(ax, (37.7, 84.7), r"$\{br,1\}$")
        text(ax, (321.4, 348), r"$\{br,2\}$")
        for origin, direction in ((lower, u), (upper, -u)):
            arrow(ax, origin + 8 * direction, origin + 108 * direction, "black", 2.5, scale=12)
        text(ax, (170, 314), r"$\mathbf{F}'_2=-\mathbf{F}'_1$")
        text(ax, (133.4, 204), r"$\hat{\mathbf{u}}_1=-\hat{\mathbf{u}}_2$", "#595959", 10)

        tendon_start = bracket(ax, 1, lower, u, 77*u - 29*normal, 14*normal,
                              ((84, 155), (140, 87)))
        tendon_end = tendon_start + 62*u
        arrow(ax, tendon_start + 2.5*u, tendon_end, TENDON, 2.6, style="<->", scale=11)
        text(ax, (178, 143), r"$\delta_{\mathrm{tendon}}$", TENDON)
        bracket(ax, 2, upper, -u, -62*u - 23*normal, 10*normal,
                ((210, 358), (294, 365)))

        fig.savefig(destination, metadata={
            "Title": "Two-bracket deflection and force-path projection",
            "Subject": "At each bracket, displacement is oblique to the black force vector. "
                       "The parallel dimension arrow shows its axial projection. Dotted guides "
                       "connect each displacement tip to its projection. The two bracket "
                       "projections and tendon elongation contribute to the change in actuator length.",
            "CreationDate": None, "ModDate": None})
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, nargs="?",
                        default=Path(__file__).resolve().parents[2] / "xiProjection.pdf")
    draw(parser.parse_args().output)

