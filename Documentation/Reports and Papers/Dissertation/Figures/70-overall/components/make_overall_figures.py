"""Chapter 7 figures: first-pass design workflow and sensorium hierarchy.

design_workflow.pdf (sec:foundation): the designer's loop that the four
enablers make possible - targets plus force model in, actuator selection,
route design under the corrected torque model, torque prediction against
targets, commit to hardware, with the redesign loop and the simulation
verification step.
sensorium_hierarchy.pdf (sec:fw_senses): three-tier control stack - dashed
supraspinal additions (basal ganglia, cerebellar gateway, vestibular and
ocular sensing) over the implemented spinal network and BPA-actuated body.
Style: Arial, 10 pt floor, no italics, black outlines, minimal fills.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

OUT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall"
plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

def box(ax, x, y, w, h, text, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.026",
                                facecolor="#f2f2f2", edgecolor="#1f1f1f",
                                linewidth=0.9,
                                linestyle="-" if not dashed else (0, (5, 3))))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10,
            linespacing=1.3)

def har(ax, x1, y1, x2, y2, label=None, lx=0.0, ly=0.12, conn=None,
        ha="center", va="bottom"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=12, color="#1f1f1f",
                                 linewidth=1.0, connectionstyle=conn))
    if label:
        ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, ha=ha, va=va,
                fontsize=10, linespacing=1.3)

# ------------------------------------------------------------- workflow
fig, ax = plt.subplots(figsize=(6.5, 3.15))
ax.set_xlim(0, 12.8)
ax.set_ylim(0, 6.3)
box(ax, 0.30, 4.35, 2.70, 1.20, "Human torque and\nrange-of-motion targets\n(OpenSim benchmark)")
box(ax, 0.30, 2.55, 2.70, 1.20, "BPA force model\n(Ch. 3)")
box(ax, 3.90, 3.45, 2.30, 1.05, "Select actuator set,\ndiameter, and\nresting length")
box(ax, 3.90, 1.35, 2.30, 1.05, "Route design under\nthe corrected torque\nmodel (Ch. 4)")
D = Polygon([(7.75, 3.55), (9.05, 4.20), (10.35, 3.55), (9.05, 2.90)],
            closed=True, facecolor="#f2f2f2", edgecolor="#1f1f1f",
            linewidth=0.9)
ax.add_patch(D)
ax.text(9.05, 3.55, "Meets\ntargets?", ha="center", va="center", fontsize=10,
        linespacing=1.2)
box(ax, 10.90, 3.05, 1.80, 1.00, "Commit to\nhardware")
box(ax, 7.55, 0.30, 2.55, 0.95,
    "Closed-loop verification\nin simulation (Ch. 6)", dashed=True)
har(ax, 3.05, 4.95, 3.85, 4.10, "force envelope", lx=0.0, ly=0.05, conn="arc3,rad=0.12")
har(ax, 3.05, 3.15, 3.85, 2.00, "force surface", lx=0.0, ly=0.05, conn="arc3,rad=-0.12")
har(ax, 5.05, 3.40, 5.05, 2.45, "actuator set", lx=0.62, ly=0.0, ha="left",
    va="center")
har(ax, 6.25, 1.87, 7.85, 3.10, "predicted torque\nand range of motion",
    lx=-0.60, ly=-0.05, conn="arc3,rad=-0.15", va="top")
har(ax, 10.40, 3.55, 10.85, 3.55, "yes", ly=0.10)
har(ax, 8.35, 2.95, 6.30, 1.90, "no: redesign route\nor reselect actuator",
    lx=-1.05, ly=0.0, conn="arc3,rad=0.35", ha="right", va="center")
har(ax, 11.80, 3.00, 11.80, 1.30, conn="angle,angleA=90,angleB=0,rad=4")
ax.plot([11.80, 10.15], [1.30, 1.30], color="#1f1f1f", linewidth=1.0,
        solid_capstyle="butt")
har(ax, 10.15, 1.30, 10.12, 1.30)
ax.set_ylim(0.0, 6.0)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\design_workflow.pdf")
fig.savefig(OUT + r"\design_workflow.png", dpi=300)
plt.close(fig)
print("wrote design_workflow")

# ------------------------------------------------------- sensorium tiers
fig, ax = plt.subplots(figsize=(6.5, 3.4))
ax.set_xlim(0, 10.0)
ax.set_ylim(0, 7.4)
ax.text(0.25, 6.95, "Supraspinal additions (future work)", fontsize=10,
        fontweight="bold", ha="left", va="center")
box(ax, 0.25, 5.05, 2.90, 1.45, "Basal ganglia\naction and mode\nselection", dashed=True)
box(ax, 3.55, 5.05, 2.90, 1.45, "Cerebellar adaptive\ngateway: gait gating\nand phase tuning",
    dashed=True)
box(ax, 6.85, 5.05, 2.90, 1.45, "Vestibular and ocular\nsensing: IMU orientation,\nvisual drift correction",
    dashed=True)
ax.text(0.25, 4.55, "Spinal network (implemented, Ch. 6)", fontsize=10,
        fontweight="bold", ha="left", va="center")
box(ax, 0.25, 2.85, 4.35, 1.35, "Two-level CPG: rhythm\ngenerator and pattern\nformation")
box(ax, 5.15, 2.85, 4.60, 1.35, "Motor-neuron pools with\nIa, Ib, and II afferents\nand heel and toe contact")
ax.text(0.25, 2.35, "Body (demonstrated)", fontsize=10, fontweight="bold",
        ha="left", va="center")
box(ax, 2.20, 0.75, 5.60, 1.30, "BPA actuators and robot body\n(balance-platform precursor of the\nvestibular channel demonstrated)")
har(ax, 5.00, 5.00, 5.00, 4.25, "descending\ncommands", lx=0.14, ha="left")
har(ax, 2.40, 2.80, 2.40, 2.10)
har(ax, 5.00, 2.80, 5.00, 2.10, "activation", lx=0.14, ha="left")
har(ax, 7.60, 2.10, 7.60, 2.80, "proprioceptive and\ncontact feedback", lx=0.14,
    ha="left")
har(ax, 8.30, 2.80, 8.30, 5.00, "body orientation and\nvision (future)", lx=0.14,
    ha="left", conn="arc3,rad=0.0")
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\sensorium_hierarchy.pdf")
fig.savefig(OUT + r"\sensorium_hierarchy.png", dpi=300)
plt.close(fig)
print("wrote sensorium_hierarchy")

# ---------------------------------------------------------- BPA cutaway
fig, ax = plt.subplots(figsize=(6.5, 2.6))
ax.set_xlim(0, 12.6)
ax.set_ylim(0, 5.0)
X0, X1 = 1.60, 8.60          # tube span
YMID, RT = 2.55, 0.80        # centerline and outer half-height
# bladder
ax.add_patch(FancyBboxPatch((X0 + 0.55, YMID - RT * 0.55), X1 - X0 - 1.10,
                            2 * RT * 0.55, boxstyle="round,pad=0.02,rounding_size=0.30",
                            facecolor="#d9d9d9", edgecolor="#555555",
                            linewidth=0.8))
# sleeve outline
ax.add_patch(FancyBboxPatch((X0 + 0.35, YMID - RT), X1 - X0 - 0.70, 2 * RT,
                            boxstyle="round,pad=0.02,rounding_size=0.42",
                            facecolor="none", edgecolor="#1f1f1f",
                            linewidth=1.1))
# two-handed braid: diagonal segments clipped by hand
import numpy as np
xs = np.linspace(X0 + 0.45, X1 - 0.45, 13)
for k in range(len(xs) - 1):
    for sign in (1, -1):
        ax.plot([xs[k], xs[k + 1]],
                [YMID - sign * RT * 0.8, YMID + sign * RT * 0.8],
                color="#4a4a4a", linewidth=0.8, solid_capstyle="butt")
# end caps
for xc in (X0, X1):
    ax.add_patch(Rectangle((xc - 0.42, YMID - 1.05), 0.84, 2.10,
                           facecolor="#bbbbbb", edgecolor="#1f1f1f",
                           linewidth=0.9))
ax.plot([X0 - 0.9, X1 + 0.9], [YMID, YMID], color="#8a8a8a", linewidth=0.8,
        linestyle=(0, (5, 3)))
ax.text(X1 + 0.98, YMID, "symmetry axis", fontsize=10, ha="left", va="center",
        color="#333333")
def leader(ax, xy, xytext, text):
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=10, ha="left",
                va="center", linespacing=1.3,
                arrowprops=dict(arrowstyle="-", color="#1f1f1f",
                                linewidth=0.7, shrinkA=0, shrinkB=2))
leader(ax, (4.1, YMID + RT), (4.6, 4.35),
       "braided aramid sleeve,\ntwo handed weaves")
leader(ax, (6.9, YMID - RT * 0.4), (8.30, 4.35), "elastomer bladder")
leader(ax, (X0 + 0.1, YMID - 1.05), (0.35, 0.75),
       "end cap\n(Onyx or aluminum)")
ax.text(0.30, 4.55, "Longitudinal section (schematic, not to scale)",
        fontsize=10, fontweight="bold", ha="left", va="center")
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\bpa_construction.pdf")
fig.savefig(OUT + r"\bpa_construction.png", dpi=300)
plt.close(fig)
print("wrote bpa_construction")
print("all Chapter 7 figures written to", OUT)
