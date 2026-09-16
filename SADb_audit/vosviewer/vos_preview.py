"""Quick PNG preview of the SADb citation network (not a VOSviewer product).
Bubble size = OpenAlex cited_by_count, color = corpus source, labels on top-15 cited.
Run: myo python vos_preview.py  (needs networkx + matplotlib - both in the myo env).
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))

nodes = {}
with open(os.path.join(HERE, "sadb_map.txt"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        nodes[r["id"]] = r

G = nx.Graph()
for n in nodes.values():
    G.add_node(n["id"], **n)
with open(os.path.join(HERE, "sadb_network.txt"), encoding="utf-8") as f:
    hdr = f.readline()
    for line in f:
        s, t, w = line.rstrip("\n").split("\t")
        if s in nodes and t in nodes:
            G.add_edge(s, t, weight=float(w))

print("graph:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "edges")
# largest connected component for layout sanity
comps = sorted(nx.connected_components(G), key=len, reverse=True)
print("components:", [len(c) for c in comps][:8])

pos = nx.spring_layout(G, k=0.22, iterations=60, seed=11, weight=None)

colors = {"originals99": "#0072B2", "demo50": "#009E73", "rest383": "#CC79A7", "digest20": "#D55E00"}
src_color = [colors.get(nodes[n]["source"], "#999999") for n in G.nodes()]
sizes = []
for n in G.nodes():
    c = float(nodes[n].get("weight") or 0)
    sizes.append(18 + 4.5 * min(c, 400) ** 0.62)

fig, ax = plt.subplots(figsize=(17, 12), dpi=130)
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.055, width=0.7, edge_color="#333333")
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=src_color,
                       linewidths=0.25, edgecolors="white", alpha=0.9)

top = sorted(G.nodes(), key=lambda n: -float(nodes[n].get("weight") or 0))[:15]
for n in top:
    x, y = pos[n]
    ax.text(x, y + 0.016, nodes[n]["label"][:46] + ("..." if len(nodes[n]["label"]) > 46 else ""),
            fontsize=7.2, ha="center", color="#111111")

handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=8,
                      markerfacecolor=colors[k], label="%s (%d)" % (k, sum(1 for n in G.nodes() if nodes[n]["source"] == k)))
           for k in colors]
ax.legend(handles=handles, loc="lower left", fontsize=9, frameon=False, title="corpus source")
ax.set_title("SADb corpus citation network - 552 papers, %d citation links (OpenAlex) - bubble size = cited-by count" % G.number_of_edges(),
             fontsize=12)
ax.axis("off")
fig.tight_layout()
out = os.path.join(HERE, "sadb_preview.png")
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
