"""Plot the two hardware envelopes at the saved governing pose."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Polygon
import numpy as np

from fr_plotDistalStages import COLORS, GOLD, INDIGO
from fr_plotPilotAnatomy import mesh_data


def capsule(ax, a, b, r, color, alpha, edge=None):
    a, b = np.asarray(a), np.asarray(b)
    v = b - a
    if np.linalg.norm(v) > 1e-8:
        side = r * np.array([-v[1], v[0]]) / np.linalg.norm(v)
        ax.add_patch(Polygon([a+side, b+side, b-side, a-side], color=color, alpha=alpha))
    for center in (a, b):
        ax.add_patch(Circle(center, r, facecolor=color, alpha=alpha, edgecolor=edge or "none"))
        if np.allclose(a, b):
            break


def draw(file, mesh):
    d = json.loads(file.read_text())
    vertices, faces = mesh_data(mesh)
    tris = 1000*vertices[faces]
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.4), layout="constrained", sharey=True)
    hit = (d["tib_ant"], d["ehl"])
    colors = (INDIGO, COLORS["magenta"])
    for ax, dim, label in zip(axes, (0, 2), ("Anterior +x", "Right +z")):
        axes2 = [dim, 1]
        ax.add_collection(PolyCollection(tris[:, :, axes2], facecolor="#a8b4c1", alpha=.28,
                                         edgecolor="none"))
        for key, color in zip(("tib_antParts", "ehlParts"), colors):
            for part in d[key]:
                capsule(ax, np.array(part["a"])[axes2], np.array(part["b"])[axes2], part["r"], color, .13)
        for part, color, group in zip(hit, colors, ("Tibialis anterior", "EHL")):
            a, b = np.array(part["a"])[axes2], np.array(part["b"])[axes2]
            capsule(ax, a, b, part["r"], color, .5)
            ax.plot([a[0], b[0]], [a[1], b[1]], "o-", color=color, lw=2,
                    label=f"{group} {part['name']}")
            if part["rigid"]:
                ax.add_patch(Circle(a, part["r"], facecolor="none", edgecolor=GOLD, lw=2))
        at = np.mean([p["a"] for p in hit]+[p["b"] for p in hit], axis=0)
        ax.set(xlim=(at[dim]-87, at[dim]+87), ylim=(at[1]-91, at[1]+91),
               xlabel=f"{label} (mm)", aspect="equal")
        ax.grid(alpha=.18)
    axes[0].set_ylabel("Superior +y in tibia frame (mm)")
    axes[1].legend(fontsize=8, loc="lower right")
    ankle, subtalar, mtp = d["angles_deg"]
    fig.suptitle(f"Tibialis anterior against EHL · {abs(d['margin_mm']):.2f} mm clearance shortfall\n"
                 f"Ankle {ankle:+g}° · subtalar {subtalar:+g}° · MTP {mtp:+g}°")
    fig.savefig(file.with_suffix(".png"), dpi=185)
    plt.close(fig)


if __name__ == "__main__":
    draw(Path(sys.argv[1]), Path(sys.argv[2]))
