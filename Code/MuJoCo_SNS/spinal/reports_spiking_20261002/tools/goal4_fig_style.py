"""Shared figure style for the goal-4 dissertation figures.

Project figure standards (from the campaign ask):
  - letter-size page, 7.5 x 10 in usable area  -> every figure <= 7.5x10 in
  - Arial at 10 pt or larger, NO italics      -> font.family Arial, size 10,
                                                 no mathtext anywhere
  - Colors.m Paul Tol palette                 -> TOL list below (repo
                                                 Code/Matlab/Colors.m)
Line conventions follow the existing dissertation figures (e.g.
animatlab_phase1_joint_angles): left side = solid indigo/blue, right side =
dashed orange. Every figure gets a matching *_alt.txt alt-text file written
next to it.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Paul Tol rainbow palette, verbatim from Code/Matlab/Colors.m
TOL = ["#FFD700", "#FFB14E", "#FA8775", "#EA5F94", "#CD34B5", "#9D02D7",
       "#0000FF"]
INDIGO = TOL[6]   # left / non-spiking / baseline
ORANGE = TOL[1]   # right / spiking / twin
PINK = TOL[3]
MAGENTA2 = TOL[5]


def apply_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "mathtext.default": "regular",   # no italic math anywhere
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def save_fig(fig, stem):
    """Save pdf + png + alt-text file. stem = path WITHOUT extension."""
    fig.savefig(stem + ".pdf")
    fig.savefig(stem + ".png", dpi=300)
    plt.close(fig)
    print("wrote", stem + ".pdf", "and", stem + ".png")
