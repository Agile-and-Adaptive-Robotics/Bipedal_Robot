from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


OUT = Path(__file__).resolve().parents[1] / "ProofFinal" / "figs" / "Background"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "blue": "#3B4CC0",
    "orange": "#E69F00",
    "teal": "#009E73",
    "magenta": "#CC79A7",
    "dark": "#2C2C2C",
    "mid": "#707070",
    "light": "#E9E9EC",
    "white": "#FFFFFF",
}


def setup(ax, letter, title, subtitle):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(0.15, 6.73, letter, fontsize=12, fontweight="bold", va="top")
    ax.text(0.85, 6.73, title, fontsize=11, fontweight="bold", va="top")
    ax.text(0.85, 6.16, subtitle, fontsize=8.5, color=COLORS["mid"], va="top")


def torque_motor(ax):
    setup(ax, "A", "Electric torque motor", "Electromagnetic rotation")
    center = (4.8, 3.25)
    ax.add_patch(Circle(center, 2.0, facecolor=COLORS["light"], edgecolor=COLORS["dark"], lw=1.8))
    for angle in range(0, 360, 45):
        x = center[0] + 1.48 * __import__("math").cos(__import__("math").radians(angle))
        y = center[1] + 1.48 * __import__("math").sin(__import__("math").radians(angle))
        ax.add_patch(Rectangle((x - 0.27, y - 0.18), 0.54, 0.36, angle=angle,
                               facecolor=COLORS["blue"], edgecolor="none"))
    ax.add_patch(Circle(center, 0.92, facecolor=COLORS["orange"], edgecolor=COLORS["dark"], lw=1.5))
    ax.add_patch(Circle(center, 0.22, facecolor=COLORS["dark"], edgecolor="none"))
    ax.plot([5.02, 8.5], [3.25, 3.25], color=COLORS["dark"], lw=6, solid_capstyle="round")
    ax.add_patch(FancyArrowPatch((3.05, 5.0), (6.15, 5.0), connectionstyle="arc3,rad=-0.42",
                                 arrowstyle="-|>", mutation_scale=16, lw=2.2, color=COLORS["blue"]))
    ax.text(7.05, 3.72, "output shaft", fontsize=8.5, ha="center")
    ax.text(4.8, 0.48, "precise, mature; often needs gearing", fontsize=8.2, ha="center", color=COLORS["mid"])


def bpa(ax):
    setup(ax, "B", "Braided pneumatic actuator", "Pressure-driven axial contraction")
    x0, x1, y0, y1 = 2.8, 7.2, 1.65, 4.9
    braid_body = Rectangle((x0, y0), x1-x0, y1-y0, facecolor="#DDE8FF",
                           edgecolor=COLORS["dark"], lw=1.6)
    ax.add_patch(braid_body)
    for k in range(-3, 8):
        line_a, = ax.plot([x0, x1], [y0 + 0.52*k, y0 + 0.52*k + 2.7],
                          color=COLORS["blue"], lw=1.0)
        line_b, = ax.plot([x0, x1], [y1 - 0.52*k, y1 - 0.52*k - 2.7],
                          color=COLORS["orange"], lw=1.0)
        line_a.set_clip_path(braid_body)
        line_b.set_clip_path(braid_body)
    ax.add_patch(Rectangle((2.25, 2.25), 0.55, 2.05, facecolor=COLORS["dark"], edgecolor="none"))
    ax.add_patch(Rectangle((7.2, 2.25), 0.55, 2.05, facecolor=COLORS["dark"], edgecolor="none"))
    ax.plot([2.5, 1.25], [3.25, 3.25], color=COLORS["dark"], lw=4)
    ax.plot([7.5, 8.75], [3.25, 3.25], color=COLORS["dark"], lw=4)
    ax.add_patch(FancyArrowPatch((1.25, 3.25), (2.05, 3.25), arrowstyle="-|>", mutation_scale=15,
                                 lw=2.0, color=COLORS["teal"]))
    ax.add_patch(FancyArrowPatch((8.75, 3.25), (7.95, 3.25), arrowstyle="-|>", mutation_scale=15,
                                 lw=2.0, color=COLORS["teal"]))
    ax.add_patch(FancyArrowPatch((5.0, 5.55), (5.0, 4.95), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.8, color=COLORS["teal"]))
    ax.text(5.0, 5.67, "air pressure", fontsize=8.5, ha="center")
    ax.text(5.0, 0.48, "light, compliant, and muscle-like", fontsize=8.2, ha="center", color=COLORS["mid"])


def hydraulic(ax):
    setup(ax, "C", "Hydraulic cylinder", "Fluid pressure on a piston")
    ax.add_patch(Rectangle((1.65, 2.0), 5.5, 2.65, facecolor="#E4F5EF", edgecolor=COLORS["dark"], lw=1.7))
    ax.add_patch(Rectangle((2.1, 2.18), 2.55, 2.29, facecolor=COLORS["teal"], alpha=0.38, edgecolor="none"))
    ax.add_patch(Rectangle((4.55, 2.08), 0.32, 2.49, facecolor=COLORS["dark"], edgecolor="none"))
    ax.plot([4.86, 8.7], [3.32, 3.32], color=COLORS["dark"], lw=7, solid_capstyle="round")
    ax.add_patch(Rectangle((2.15, 4.65), 0.45, 0.65, facecolor=COLORS["dark"], edgecolor="none"))
    ax.add_patch(Rectangle((6.15, 4.65), 0.45, 0.65, facecolor=COLORS["dark"], edgecolor="none"))
    ax.add_patch(FancyArrowPatch((2.38, 5.55), (2.38, 5.05), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.8, color=COLORS["teal"]))
    ax.add_patch(FancyArrowPatch((7.25, 3.95), (8.85, 3.95), arrowstyle="-|>", mutation_scale=16,
                                 lw=2.2, color=COLORS["orange"]))
    ax.text(2.38, 5.67, "pressurized fluid", fontsize=8.5, ha="center")
    ax.text(7.95, 4.18, "linear force", fontsize=8.5, ha="center")
    ax.text(5.0, 0.48, "high force density; heavy infrastructure", fontsize=8.2, ha="center", color=COLORS["mid"])


def dea(ax):
    setup(ax, "D", "Dielectric elastomer actuator", "Electric-field-induced thinning")
    body = Polygon([[2.3, 2.35], [7.7, 2.35], [7.2, 4.35], [2.8, 4.35]], closed=True,
                   facecolor="#F6DFF0", edgecolor=COLORS["dark"], lw=1.6)
    ax.add_patch(body)
    ax.plot([2.8, 7.2], [4.35, 4.35], color=COLORS["magenta"], lw=6, solid_capstyle="round")
    ax.plot([2.3, 7.7], [2.35, 2.35], color=COLORS["magenta"], lw=6, solid_capstyle="round")
    ax.text(5.0, 4.75, "+ electrode", fontsize=8.5, ha="center")
    ax.text(5.0, 1.74, "− electrode", fontsize=8.5, ha="center")
    ax.add_patch(FancyArrowPatch((5.0, 5.65), (5.0, 4.75), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.8, color=COLORS["magenta"]))
    ax.add_patch(FancyArrowPatch((5.0, 1.05), (5.0, 1.9), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.8, color=COLORS["magenta"]))
    ax.add_patch(FancyArrowPatch((2.7, 3.3), (1.45, 3.3), arrowstyle="-|>", mutation_scale=15,
                                 lw=2.0, color=COLORS["blue"]))
    ax.add_patch(FancyArrowPatch((7.3, 3.3), (8.55, 3.3), arrowstyle="-|>", mutation_scale=15,
                                 lw=2.0, color=COLORS["blue"]))
    ax.text(5.0, 0.48, "large strain; requires high voltage", fontsize=8.2, ha="center", color=COLORS["mid"])


fig, axs = plt.subplots(2, 2, figsize=(10.6, 7.6))
fig.patch.set_facecolor("white")
torque_motor(axs[0, 0])
bpa(axs[0, 1])
hydraulic(axs[1, 0])
dea(axs[1, 1])

for ax in axs.flat:
    for spine in ax.spines.values():
        spine.set_visible(False)

fig.subplots_adjust(left=0.035, right=0.985, top=0.985, bottom=0.03, wspace=0.08, hspace=0.14)
for name in ("actuator_technology_comparison.pdf", "actuator_technology_comparison.png"):
    kwargs = {"dpi": 300} if name.endswith(".png") else {}
    fig.savefig(OUT / name, bbox_inches="tight", facecolor="white", **kwargs)
plt.close(fig)

print(OUT / "actuator_technology_comparison.pdf")
