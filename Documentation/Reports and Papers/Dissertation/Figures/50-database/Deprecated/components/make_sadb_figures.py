"""Generate the two Sensory Afferent Database example figures for Chapter 5.

Figure 1 (sadb_pivot): afferent class x animal preparation pivot (heatmap with
  counts and margins) over the retained papers, counted per paper-tag pair.
Figure 2 (sadb_clusters): Louvain cluster bubble chart - median publication
  year vs retained-paper count, bubble area = audited cluster size, fill depth
  = retained fraction.

Data: SADb_audit/export/sadb_export.csv, sadb_layout.json, cluster_labels.json
Style: Arial, 10 pt floor, no italics, lightness-ramp single-hue colors
(CVD-safe); matches the dissertation figure contract.
"""
import csv, json, math, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

EXPORT = r"D:\Github\Bipedal_Robot\SADb_audit\export"
OUT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\50-database"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 10, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

rows = list(csv.DictReader(open(EXPORT + r"\sadb_export.csv", encoding="utf-8-sig")))
retained = [r for r in rows if r["prune"] != "True" and r["archived"] != "True"]

def tags(s):
    return [t.strip() for t in s.replace(";", ",").split(",") if t.strip()]

# ---------------------------------------------------------------- Figure 1
AFF_ROWS = [
    ("Ia", ["Ia"]),
    ("Ib", ["Ib"]),
    ("II", ["II"]),
    ("Cutaneous", ["Cutaneous"]),
    ("Other classes", ["Mechanosensory", "Flexor reflex afferents", "III/IV", "Nociceptive"]),
]
ANIM_COLS = [
    ("Cat", ["Cat"]),
    ("Human", ["Human"]),
    ("Rodent", ["Mice", "Rat"]),
    ("Arthropod", ["Insects", "Stick Insect", "Cockroach", "Arthropods"]),
    ("Lower vertebrate", ["Lamprey", "Salamander", "Frog", "Zebrafish"]),
    ("Other or generic", ["Mammals", "Vertebrates", "Dog"]),
    ("No preparation tag", []),
]
classified = [r for r in retained if tags(r["afferents"])]
n_untagged_animal = sum(1 for r in classified if not tags(r["animals"]))

grid = [[0] * len(ANIM_COLS) for _ in AFF_ROWS]
col_papers = [set() for _ in ANIM_COLS]   # distinct papers per column
row_papers = [set() for _ in AFF_ROWS]    # distinct papers per row
col_all = set()
for r in classified:
    at, an = tags(r["afferents"]), tags(r["animals"])
    for i, (_, akeys) in enumerate(AFF_ROWS):
        if set(at) & set(akeys):
            row_papers[i].add(r["id"])
            for j, (_, nkeys) in enumerate(ANIM_COLS):
                hit = (set(nkeys) & set(an)) if nkeys else (not an)
                if hit:
                    grid[i][j] += 1
                    col_papers[j].add(r["id"])
    if an:
        col_all.add(r["id"])

nr, nc = len(AFF_ROWS), len(ANIM_COLS)
fig, ax = plt.subplots(figsize=(6.5, 2.9))
cmap = plt.get_cmap("Blues")
vmax = max(max(row) for row in grid)
for i in range(nr):
    for j in range(nc):
        v = grid[i][j]
        lum = 0.0 if vmax == 0 else v / vmax
        face = cmap(0.12 + 0.78 * (lum ** 0.75))
        ax.add_patch(Rectangle((j, nr - 1 - i), 0.92, 0.92, facecolor=face,
                               edgecolor="#777777", linewidth=0.5))
        if v > 0:
            tc = "white" if lum > 0.55 else "#1a1a1a"
            ax.text(j + 0.46, nr - 1 - i + 0.46, str(v), ha="center",
                    va="center", fontsize=10, color=tc)
        else:
            ax.text(j + 0.46, nr - 1 - i + 0.46, "0", ha="center",
                    va="center", fontsize=10, color="#a8a8a8")
# margins
for i in range(nr):
    ax.text(nc + 0.10, nr - 1 - i + 0.46, str(len(row_papers[i])), ha="left",
            va="center", fontsize=10, fontweight="bold")
for j in range(nc):
    ax.text(j + 0.46, -0.28, str(len(col_papers[j])), ha="center", va="center",
            fontsize=10, fontweight="bold")
ax.text(nc + 0.10, -0.28, str(len(classified)), ha="left", va="center",
        fontsize=10, fontweight="bold")
# headers
for i, (name, _) in enumerate(AFF_ROWS):
    ax.text(-0.10, nr - 1 - i + 0.46, name, ha="right", va="center", fontsize=10)
for j, (name, _) in enumerate(ANIM_COLS):
    ax.text(j + 0.46, nr + 0.10, name, ha="center", va="bottom", fontsize=10,
            rotation=32)
ax.text(nc + 0.10, nr + 0.30, "Papers", ha="left", va="bottom", fontsize=10,
        fontweight="bold")
ax.text(-0.10, nr + 0.10, "Afferent class", ha="right", va="bottom", fontsize=10,
        fontweight="bold")
ax.set_xlim(-1.9, nc + 1.15)
ax.set_ylim(-0.95, nr + 1.15)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT + r"\sadb_pivot.pdf")
fig.savefig(OUT + r"\sadb_pivot.png", dpi=300)
plt.close(fig)
print("pivot: base", len(classified), "of", len(retained),
      "retained; afferent-tagged papers lacking an animal tag:", n_untagged_animal)
for i, (name, _) in enumerate(AFF_ROWS):
    print("  row", name, "papers:", len(row_papers[i]))
for j, (name, _) in enumerate(ANIM_COLS):
    print("  col", name, "papers:", len(col_papers[j]))

# ---------------------------------------------------------------- Figure 2
layout = json.load(open(EXPORT + r"\sadb_layout.json", encoding="utf-8"))
labels = json.load(open(EXPORT + r"\cluster_labels.json", encoding="utf-8"))
by_id = {r["id"]: r for r in rows}
members = {}
for rid, info in layout.items():
    members.setdefault(info["cl"], []).append(by_id[rid])

def wrap(s, w=22):
    words, lines, cur = s.split(), [], ""
    for w_ in words:
        if len(cur) + len(w_) + 1 > w and cur:
            lines.append(cur); cur = w_
        else:
            cur = (cur + " " + w_).strip()
    if cur:
        lines.append(cur)
    return "\n".join(lines)

unclustered = len(members.pop(-1, []))
print("unclustered records (excluded from the chart):", unclustered)
stats = []
for cl, recs in members.items():
    ret = [r for r in recs if r["prune"] != "True" and r["archived"] != "True"]
    yrs = [int(r["year"]) for r in ret if r["year"].strip().isdigit()]
    if not yrs:
        yrs = [int(r["year"]) for r in recs if r["year"].strip().isdigit()]
    stats.append({"cl": cl, "label": labels[str(cl)], "audited": len(recs),
                  "retained": len(ret), "median_year": statistics.median(yrs)})
stats.sort(key=lambda s: s["retained"], reverse=True)
live = [s for s in stats if s["retained"] > 0]
dead = [s for s in stats if s["retained"] == 0]
dead_papers = sum(s["audited"] for s in dead)
dead_years = []
for s in dead:
    for r in members[s["cl"]]:
        if r["year"].strip().isdigit():
            dead_years.append(int(r["year"]))
dead_median = statistics.median(dead_years)
print("live clusters:", len(live), " fully pruned clusters:", len(dead),
      "(", dead_papers, "papers, median year", dead_median, ")")

fig, ax = plt.subplots(figsize=(6.5, 5.0))
smax = max(s["audited"] for s in stats)
for s in live:
    r = 0.62 * math.sqrt(s["audited"] / smax)       # area-proportional
    frac = s["retained"] / s["audited"]
    ax.scatter(s["median_year"], s["retained"], s=(r * 2 * 72) ** 2 / 4,
               color=cmap(0.15 + 0.75 * frac), edgecolor="#1f1f1f",
               linewidth=0.6, zorder=3)
r = 0.62 * math.sqrt(dead_papers / smax)
ax.scatter(dead_median, 0, s=(r * 2 * 72) ** 2 / 4, color="#c8c8c8",
           edgecolor="#1f1f1f", linewidth=0.6, zorder=3)
s_labels = {
    0: (0, 1), 1: (0, 1), 2: (0, 1), 3: (0, 1), 4: (0, 1), 5: (0, 1),
    6: (0, 1),
}
for s in live:
    r = 0.62 * math.sqrt(s["audited"] / smax)
    dy = s_labels.get(s["cl"], (0, 1))[1]
    ax.text(s["median_year"], s["retained"] + dy * (r + 2.2),
            wrap(s["label"]) + "\n" + str(s["retained"]) + " of " + str(s["audited"]),
            ha="center", va="bottom" if dy > 0 else "top", fontsize=10, zorder=4)
ax.text(dead_median, r + 2.2,
        "Eleven satellite clusters,\nfully pruned (" + str(dead_papers) + " papers)",
        ha="center", va="bottom", fontsize=10, zorder=4)
ax.set_xlabel("Median publication year of retained papers", fontsize=11)
ax.set_ylabel("Retained papers in cluster", fontsize=11)
ax.set_xlim(1988, 2032)
ax.set_ylim(-6, max(s["retained"] for s in live) * 1.32)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=10)
handles = [
    Line2D([], [], marker="o", ls="", markersize=9, markerfacecolor=cmap(0.9),
           markeredgecolor="#1f1f1f", label="fully retained"),
    Line2D([], [], marker="o", ls="", markersize=9, markerfacecolor=cmap(0.15),
           markeredgecolor="#1f1f1f", label="largely pruned"),
    Line2D([], [], marker="o", ls="", markersize=14, markerfacecolor="#9ecae1",
           markeredgecolor="#1f1f1f", label="200-paper cluster"),
    Line2D([], [], marker="o", ls="", markersize=7, markerfacecolor="#9ecae1",
           markeredgecolor="#1f1f1f", label="25-paper cluster"),
]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=10,
          labelspacing=1.0, borderpad=0.6, handletextpad=0.5)
fig.tight_layout()
fig.savefig(OUT + r"\sadb_clusters.pdf")
fig.savefig(OUT + r"\sadb_clusters.png", dpi=300)
plt.close(fig)
print("wrote sadb_pivot and sadb_clusters to", OUT)
