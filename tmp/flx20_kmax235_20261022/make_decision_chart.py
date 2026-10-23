"""Margin-vs-tendon decision chart for the flexor 2x20 mm pair (KMAX 0.235 vs 0.255)."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
rows = []
with open(here / "flx20_kmax235_tendon_sweep.csv", newline="") as f:
    for r in csv.DictReader(f):
        if int(r["conv_idx"]) == 2:  # script convention (identical to physical within 0.1 Nm)
            rows.append((int(r["KMAX_idx"]), float(r["dTendon_mm"]),
                         float(r["minMargin_Nm"]), int(r["nFailAt620kPa"])))

fig, ax = plt.subplots(figsize=(7.0, 4.4), dpi=200)
for kidx, color, label in ((1, "#0072B2", "KMAX 0.255 (old)"), (2, "#D55E00", "KMAX 0.235 (corrected)")):
    pts = sorted((r for r in rows if r[0] == kidx), key=lambda r: r[1])
    xs = [p[1] for p in pts]
    ys = [p[2] if p[2] != float("-inf") else float("nan") for p in pts]
    ax.plot(xs, ys, "-o", ms=3.5, lw=2, color=color, label=label)

ax.axhline(0, color="k", lw=1.5)
ax.axvline(0, color="0.5", lw=1, ls=":")
ax.fill_between([-105, 50], -40, 0, color="0.92", zorder=0)
ax.text(-98, -6.5, "fails to meet human torque at 620 kPa", fontsize=8.5, color="0.35")
ax.text(22, 14.5, "meets / exceeds\n(620 kPa)", fontsize=8.5, color="0.35", ha="center")
ax.annotate("current tendon\n(5 mm in model)", xy=(0, 3.96), xytext=(12, 9.5),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=1))
ax.annotate("shorter tendon (left) = more margin", xy=(-60, 22.8), xytext=(-8, 26),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=1))
ax.set_xlabel("Change in tendon length from current (mm)  [negative = shorter]")
ax.set_ylabel("Worst-angle torque margin vs human at 620 kPa (N·m)")
ax.set_title("Flexor 2x20 mm pair: meet-or-exceed margin vs tendon length")
ax.set_xlim(-105, 50)
ax.set_ylim(-16, 36)
ax.legend(loc="upper right", frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(here / "flx20_tendon_decision_chart.png")
print("wrote", here / "flx20_tendon_decision_chart.png")
