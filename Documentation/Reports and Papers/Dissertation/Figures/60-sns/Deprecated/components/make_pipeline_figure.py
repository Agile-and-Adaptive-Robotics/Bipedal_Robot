"""Simulation-to-robot pipeline figure (Chapter 6, sec:preliminary_sim_methods).

Three phase bands left to right: model conversion, closed-loop simulation,
hardware deployment. The closed loop is a vertical exchange (MuJoCo above,
SNS-Toolbox below) with the two signal labels set beside their arrows at full
band width. Content follows the dissertation text.
Style: Arial, 10 pt floor, no italics, black outlines, minimal fills.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

OUT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns"
plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

W = 12.6
fig, ax = plt.subplots(figsize=(6.5, 2.85))

def box(x, y, w, h, text):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.026",
                                facecolor="#f2f2f2", edgecolor="#1f1f1f",
                                linewidth=0.9))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10,
            linespacing=1.3)

def har(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=12, color="#1f1f1f",
                                 linewidth=1.0))

def seg(x1, y1, x2, y2):
    ax.plot([x1, x2], [y1, y2], color="#1f1f1f", linewidth=1.0,
            solid_capstyle="butt")

def phase(x, w, title):
    ax.add_patch(Rectangle((x, 0.30), w, 4.35, facecolor="none",
                           edgecolor="#8a8a8a", linewidth=0.8,
                           linestyle=(0, (4, 3))))
    ax.text(x + w / 2, 4.80, title, ha="center", va="bottom", fontsize=10,
            fontweight="bold")

# Phase 1: model conversion
phase(0.15, 3.90, "Model conversion")
box(0.45, 3.05, 3.30, 1.10, "OpenSim robot-body model\n(.osim, Gait2392 routing)")
box(0.45, 1.75, 3.30, 0.85, "MyoConverter")
box(0.45, 0.55, 3.30, 0.85, "MuJoCo MJCF model\n(muscle actuators with\nactivation state)")
har(2.10, 3.00, 2.10, 2.65)
har(2.10, 1.70, 2.10, 1.45)

# Phase 2: closed-loop simulation (vertical exchange)
phase(4.30, 4.60, "Closed-loop simulation")
box(4.60, 2.95, 2.15, 1.20, "MuJoCo physics\nwith BPA force\nmodel (Ch. 3)")
box(4.60, 0.75, 2.15, 1.05, "SNS-Toolbox\nnetwork (Python)")
har(5.05, 1.83, 5.05, 2.92)                       # up: activation
ax.text(5.22, 2.38, "per-muscle activation", ha="left", va="center",
        fontsize=10)
har(6.30, 2.92, 6.30, 1.83)                       # down: afferents
ax.text(6.47, 2.38, "Ia, Ib, II afferents;\nheel and toe contact",
        ha="left", va="center", fontsize=10, linespacing=1.3)

# Phase 3: hardware deployment
phase(9.15, 3.30, "Hardware deployment")
box(9.45, 3.05, 2.70, 1.10, "Teensy-class MCUs\nper-joint valve, pressure,\nand contact loops")
box(9.45, 1.70, 2.70, 1.05, "Jetson Orin\nSNS network evaluation")
box(9.45, 0.55, 2.70, 0.85, "Treadmill biped")
har(10.80, 3.00, 10.80, 2.78)
har(10.80, 1.65, 10.80, 1.43)

# MJCF -> MuJoCo (elbow under the loop, entering the physics box)
seg(3.80, 0.98, 4.12, 0.98)
seg(4.12, 0.98, 4.12, 3.55)
har(4.12, 3.55, 4.55, 3.55)

# loop -> hardware (elbow up into the Teensy box)
seg(6.78, 1.28, 8.98, 1.28)
seg(8.98, 1.28, 8.98, 3.60)
har(8.98, 3.60, 9.42, 3.60)
ax.text(8.90, 2.45, "computer-in-the-loop,\nthen embedded execution",
        ha="right", va="center", fontsize=10, linespacing=1.3)

ax.set_xlim(0.0, W)
ax.set_ylim(0.0, 5.30)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\pipeline.pdf")
fig.savefig(OUT + r"\pipeline.png", dpi=300)
plt.close(fig)
print("wrote pipeline to", OUT)
