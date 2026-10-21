"""Five-stage training curriculum figure for Chapter 6 (sec:sns_curriculum).

One full-width row per stage: number chip, stage name, capability added,
stage score. Rows are linked by a seeding arrow chain; the chip fill encodes
the preparation's support level (suspended air -> stance -> ground).
Stages and scores follow the dissertation text exactly.
Style: Arial, 10 pt floor, no italics, CVD-safe lightness ramp.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns"
plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

stages = [
    ("1", "Deafferented\nair stepping", "half-center oscillation, legs\nsuspended, afferents silent",
     "rhythm and swing depth,\nrhythm-gated objective", "#dbe9f6"),
    ("2", "Afferented\nair stepping", "heel, toe, and load\nfeedback active",
     "rhythm and swing depth,\nrhythm-gated objective", "#b8d1ea"),
    ("3", "Standing\nbalance", "posture held under\npush perturbations",
     "posture maintenance\nunder push", "#8fb6dc"),
    ("4", "Ground\nwalking", "ground contact\nand forward gait",
     "cycle-shape kinematic\nscore vs human reference", "#6a9cc9"),
    ("5", "Contact-driven\nground walking", "contact-driven pathway\nvariant of the walker",
     "cycle-shape kinematic\nscore vs human reference", "#4a82b4"),
]

W, H = 10.0, 5.0            # axis units close to cm at \textwidth inclusion
ROW, GAP = 0.78, 0.30
fig, ax = plt.subplots(figsize=(6.5, 3.55))
x_name, x_adds, x_score = 1.35, 4.05, 7.15
top = 0.20 + 5 * ROW + 4 * GAP
for i, (num, name, adds, score, fc) in enumerate(stages):
    y = top - ROW - i * (ROW + GAP)
    ax.add_patch(FancyBboxPatch((0.20, y), W - 0.40, ROW,
                                boxstyle="round,pad=0.022",
                                facecolor="#f7f7f7", edgecolor="#1f1f1f",
                                linewidth=0.9))
    ax.add_patch(FancyBboxPatch((0.34, y + 0.13), 0.52, ROW - 0.26,
                                boxstyle="round,pad=0.018", facecolor=fc,
                                edgecolor="#1f1f1f", linewidth=0.7))
    ax.text(0.60, y + ROW / 2, num, ha="center", va="center", fontsize=11,
            fontweight="bold")
    ax.text(x_name, y + ROW / 2, name, ha="left", va="center", fontsize=10,
            fontweight="bold", linespacing=1.25)
    ax.text(x_adds, y + ROW / 2, adds, ha="left", va="center", fontsize=10,
            linespacing=1.25)
    ax.text(x_score, y + ROW / 2, score, ha="left", va="center", fontsize=10,
            linespacing=1.25)
    if i < 4:
        ax.add_patch(FancyArrowPatch((0.60, y - 0.035), (0.60, y - GAP + 0.035),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color="#1f1f1f", linewidth=1.0))
ax.text(0.20, 0.06,
        "Each stage's winner seeds the next; only the newly activated pathway gains are searched.",
        ha="left", va="top", fontsize=10)
ax.text(0.20, -0.22,
        "Finalists at every stage must pass the perturbation and basin gates. "
        "Darker chip = less harness support.",
        ha="left", va="top", fontsize=10, color="#333333")
ax.set_xlim(0.0, W)
ax.set_ylim(-0.62, top + 0.10)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\sns_curriculum.pdf")
fig.savefig(OUT + r"\sns_curriculum.png", dpi=300)
plt.close(fig)
print("wrote sns_curriculum to", OUT)
