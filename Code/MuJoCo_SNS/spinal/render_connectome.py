"""Render the BilateralRG connectome diagram from Ben's canonical JSON.
Uses the connectome editor's own data format (nodes + edges) and renders
a layered schematic with proper left/right separation."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BEN_JSON = Path(r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Walker_2_Layer_CPG_BilateralRG_connectome.json")
OUT_PNG = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\slideshow\bilateralrg_connectome.png")
OUT_PDF = OUT_PNG.with_suffix(".pdf")


def classify_node(label):
    """Assign a layer based on the node label."""
    l = label.lower()
    if "rg" in l and ("ext" in l or "flx" in l):
        return "rg"
    if "pf" in l and ("hip" in l or "knee" in l or "ank" in l):
        return "pf"
    if "mn" in l:
        return "mn"
    if "in" in l or "renshaw" in l or "re " in l:
        return "in"
    if "heel" in l or "toe" in l or "ia" in l or "ib" in l or "ii" in l:
        return "sensory"
    return "other"


def classify_edge(edge):
    """Determine if an edge is excitatory or inhibitory."""
    etype = str(edge.get("type", edge.get("sign", ""))).lower()
    if "inh" in etype or "gaba" in etype or "gly" in etype:
        return "inh"
    return "exc"


def main():
    data = json.loads(BEN_JSON.read_text(encoding="utf8"))
    nodes = data.get("nodes", [])
    edges = data.get("edges", data.get("synapses", []))
    print(f"Loaded {len(nodes)} nodes, {len(edges)} edges")

    # Build node info
    node_info = {}
    for n in nodes:
        label = n.get("label", n.get("name", n.get("id", "?")))
        side = "L" if "L " in str(label) or str(label).startswith("L") else \
               "R" if "R " in str(label) or str(label).startswith("R") else "C"
        layer = classify_node(str(label))
        node_info[label] = {"side": side, "layer": layer, "label": label}

    # Layout: layers as horizontal bands, sides as columns
    LAYER_ORDER = ["sensory", "rg", "pf", "in", "mn", "other"]
    LAYER_Y = {layer: i * 2.0 for i, layer in enumerate(LAYER_ORDER)}
    SIDE_X = {"L": -3.0, "R": 3.0, "C": 0.0}

    # Position nodes
    positions = {}
    layer_counts = {}
    for label, info in node_info.items():
        layer = info["layer"]
        side = info["side"]
        key = (layer, side)
        layer_counts[key] = layer_counts.get(key, 0) + 1
        idx = layer_counts[key]
        # Spread nodes within their layer/side band
        spread = 0.8
        offset = (idx - 1) * spread - (layer_counts[key] - 1) * spread / 2
        positions[label] = (SIDE_X[side] + offset, LAYER_Y[layer])

    # Draw
    fig, ax = plt.subplots(1, 1, figsize=(16, 12))

    # Draw edges
    for e in edges:
        src = str(e.get("from", e.get("source", "")))
        dst = str(e.get("to", e.get("target", "")))
        if src not in positions or dst not in positions:
            continue
        x1, y1 = positions[src]
        x2, y2 = positions[dst]
        kind = classify_edge(e)
        if kind == "exc":
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                       arrowprops=dict(arrowstyle="->", color="#009E73", lw=0.8, alpha=0.6))
        else:
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                       arrowprops=dict(arrowstyle="-|>", color="#D55E00", lw=0.8, alpha=0.6,
                                    linestyle="--"))

    # Draw nodes
    LAYER_COLORS = {"sensory": "#F0E442", "rg": "#56B4E9", "pf": "#CC79A7",
                    "in": "#E69F00", "mn": "#0072B2", "other": "#999999"}
    for label, (x, y) in positions.items():
        info = node_info[label]
        color = LAYER_COLORS.get(info["layer"], "#999999")
        ax.scatter(x, y, s=300, c=color, edgecolors="black", linewidths=0.8, zorder=5)
        # Label (shortened)
        short = str(label).replace("contact", "").strip()
        if len(short) > 18:
            short = short[:15] + "..."
        ax.annotate(short, (x, y), textcoords="offset points", xytext=(0, 12),
                   ha="center", fontsize=5.5, zorder=6)

    # Layer labels
    for layer, y in LAYER_Y.items():
        if any(info["layer"] == layer for info in node_info.values()):
            ax.text(-9.5, y, layer.upper(), fontsize=9, fontweight="bold",
                   ha="left", va="center", color="#555555")

    # Legend
    exc_patch = mpatches.Patch(color="#009E00", label="Excitatory")
    inh_patch = mpatches.Patch(color="#D55E00", label="Inhibitory")
    ax.legend(handles=[exc_patch, inh_patch], loc="lower right", fontsize=9)

    ax.set_xlim(-10, 10)
    ax.set_ylim(-1, len(LAYER_ORDER) * 2 + 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("BilateralRG Connectome (Ben's canonical wiring, 2026-10-03)\n"
                "Toe → ipsilateral ankle ext only; Heel → all ext + ipsi RG reset + contra RG-flx",
                fontsize=11)

    plt.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT_PNG}")
    print(f"Saved {OUT_PDF}")


if __name__ == "__main__":
    main()
