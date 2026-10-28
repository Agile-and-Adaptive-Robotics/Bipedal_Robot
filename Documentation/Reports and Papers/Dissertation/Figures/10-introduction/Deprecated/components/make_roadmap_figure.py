"""Dissertation roadmap figure (Chapter 1, sec:organization end).

Two-row block diagram of the program flow, blocks tagged with the chapters
that establish them, feedback arc from the physical robot back to the
benchmark. Solid blocks are built and validated in the dissertation; the
dashed block is the future-work endpoint. Black outlines, minimal fills,
one arrow style.
Style: Arial, 10 pt floor, no italics.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\10-introduction"
plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

top_row = [
    ("Human biomechanics\nbenchmark\n(Ch. 2)", True),
    ("BPA force models\nvalidated on test data\n(Ch. 3)", True),
    ("Torque correction\nfactors and routes\n(Ch. 4)", True),
]
bot_row = [
    ("Sensory Afferent\nDatabase\n(Ch. 5)", True),
    ("Neuromechanical\nsimulation tools\n(Ch. 6)", True),
    ("Physical BPA robot\nwith SNS control", False),
]
arrow_labels = {
    ("t0", "t1"): "targets",
    ("t1", "t2"): "actuator model",
    ("t2", "b1"): "routes and torque",
    ("b0", "b1"): "sensory channels",
    ("b1", "b2"): "controller and body",
}

W, H = 10.0, 4.6
BW, BH = 2.55, 1.05
Y1, Y0 = 3.05, 1.15
X0 = 0.55
GAP = (W - 2 * X0 - 3 * BW) / 2

fig, ax = plt.subplots(figsize=(6.5, 3.0))
centers = {}
for row, y, pre in ((top_row, Y1, "t"), (bot_row, Y0, "b")):
    for i, (txt, solid) in enumerate(row):
        x = X0 + i * (BW + GAP)
        ax.add_patch(FancyBboxPatch((x, y), BW, BH, boxstyle="round,pad=0.028",
                                    facecolor="#f2f2f2" if solid else "white",
                                    edgecolor="#1f1f1f", linewidth=0.9,
                                    linestyle="-" if solid else (0, (5, 3))))
        ax.text(x + BW / 2, y + BH / 2, txt, ha="center", va="center",
                fontsize=10, linespacing=1.3)
        centers[pre + str(i)] = (x, y, x + BW, y + BH)

def arrow(a, b, label, below=True):
    xa, ya, xb, yb = centers[a]
    ua, va, ub, vb = centers[b]
    if abs(ya - va) < 0.01:                      # horizontal
        ax.add_patch(FancyArrowPatch((xb + 0.05, ya + BH / 2), (ua - 0.05, va + BH / 2),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color="#1f1f1f", linewidth=1.0))
        ax.text((xb + ua) / 2, ya + BH / 2 - 0.12, label, ha="center",
                va="top", fontsize=10)
    else:                                         # elbow down (t2 -> b1)
        sx, sy = xa + BW / 2, ya
        ex, ey = ua + BW / 2, yb + BH
        ax.add_patch(FancyArrowPatch((sx, sy - 0.05), (ex, ey + 0.05),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color="#1f1f1f", linewidth=1.0,
                                     connectionstyle="angle,angleA=-90,angleB=0,rad=6"))
        ax.text(sx + 0.16, (sy + ey) / 2, label, ha="left", va="center",
                fontsize=10)

for (a, b), lab in arrow_labels.items():
    arrow(a, b, lab)

# feedback arc: physical robot back to the benchmark, "biological insight"
fx0, fy0, fx1, fy1 = centers["b2"]
ax0_, ay0, ax1_, ay1 = centers["t0"]
ax.add_patch(FancyArrowPatch((fx0 + BW / 2, fy0 - 0.05), (ax0_ + BW / 2, ay0 + 0.05),
                             arrowstyle="-|>", mutation_scale=12,
                             color="#1f1f1f", linewidth=1.0,
                             connectionstyle="arc3,rad=-0.30"))
ax.text(W / 2, -0.02, "biological insight", ha="center", va="top", fontsize=10)

# framing key
ax.plot([X0, X0 + 0.42], [Y1 + BH + 0.30, Y1 + BH + 0.30], color="#1f1f1f",
        linewidth=0.9)
ax.text(X0 + 0.52, Y1 + BH + 0.30, "built and validated in this dissertation",
        ha="left", va="center", fontsize=10)
ax.plot([X0 + 4.30, X0 + 4.72], [Y1 + BH + 0.30, Y1 + BH + 0.30],
        color="#1f1f1f", linewidth=0.9, linestyle=(0, (5, 3)))
ax.text(X0 + 4.82, Y1 + BH + 0.30, "future work", ha="left", va="center",
        fontsize=10)

ax.set_xlim(0.0, W)
ax.set_ylim(-0.55, Y1 + BH + 0.60)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\roadmap.pdf")
fig.savefig(OUT + r"\roadmap.png", dpi=300)
plt.close(fig)
print("wrote roadmap to", OUT)
