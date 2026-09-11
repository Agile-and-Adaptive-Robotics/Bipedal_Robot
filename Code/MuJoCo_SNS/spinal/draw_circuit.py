"""Rybak-style schematic of one side's spinal circuit (RG + PF + MN +
proprioception), drawn to Ben's SNS diagram conventions:
open circle = neuron, solid black dot at the target = inhibitory synapse,
open triangle = excitatory synapse. Output: spinal_circuit.png.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon

FIG_W, FIG_H = 15.0, 9.0
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 15)
ax.set_ylim(0, 9)
ax.axis("off")

R = 0.28          # neuron radius


def neuron(x, y, label, sub=None, color="white"):
    ax.add_patch(Circle((x, y), R, fc=color, ec="black", lw=1.6, zorder=3))
    ax.text(x, y + (0.06 if sub else 0), label, ha="center", va="center",
            fontsize=8.5, zorder=4)
    if sub:
        ax.text(x, y - 0.12, sub, ha="center", va="center", fontsize=7,
                style="italic", zorder=4)
    return (x, y)


def syn(p, q, exc: bool, label=None, loff=(0.12, 0.12), lw=1.5):
    """Synapse arrow p->q. Excitatory: open triangle head; inhibitory:
    solid black dot at the target end (drawn as a short arrow + marker)."""
    arrow = FancyArrowPatch(p, q, arrowstyle="-", lw=lw, color="black",
                            shrinkA=R, shrinkB=R + (0.09 if exc else 0.12),
                            zorder=2)
    ax.add_patch(arrow)
    import numpy as np
    p, q = np.array(p, float), np.array(q, float)
    d = q - p
    d = d / np.linalg.norm(d)
    tip = q - d * (R + 0.12)
    perp = np.array([-d[1], d[0]])
    if exc:
        # open triangle at the target
        tri = Polygon([tip, tip - d * 0.18 + perp * 0.09,
                       tip - d * 0.18 - perp * 0.09],
                      closed=True, fc="white", ec="black", lw=1.4, zorder=3)
        ax.add_patch(tri)
    else:
        ax.add_patch(Circle(tip, 0.085, fc="black", ec="black", zorder=3))
    if label:
        mid = (p + q) / 2 + np.array(loff)
        ax.text(*mid, label, fontsize=7.5, ha="center",
                style="italic", color="0.25")


def box(x, y, w, h, label):
    ax.add_patch(plt.Rectangle((x - w / 2, y - h / 2), w, h, fc="0.93",
                               ec="black", lw=1.4, zorder=3))
    ax.text(x, y, label, ha="center", va="center", fontsize=8.5, zorder=4)
    return (x, y)


# ---------------- titles ----------------
ax.text(7.5, 8.6, "gait2392 spinal cord - one side (RG + PF + MN)",
        ha="center", fontsize=13, weight="bold")
ax.text(7.5, 8.25, "the other side is identical; RG-F cells mutually "
        "inhibit across sides (strong), RG-E weakly",
        ha="center", fontsize=8.5, style="italic", color="0.3")

# ---------------- descending inputs ----------------
drv = neuron(0.7, 6.6, "DRIVE", "MLR surr.")
pos = neuron(0.7, 5.6, "POSTURE", "tone")
balp = neuron(0.7, 4.6, "BAL_PF", "balance")
bald = neuron(0.7, 3.8, "BAL_DF", "balance")

# ---------------- rhythm generator (half-center + adaptation) ----------------
rge = neuron(3.0, 6.8, "RG-E", "extensor CPG")
rgf = neuron(3.0, 5.2, "RG-F", "flexor CPG")
ade = neuron(1.9, 7.4, "ADAP-E", "slow")
adf = neuron(1.9, 4.6, "ADAP-F", "slow")

syn(drv, rge, True, "0.9", (0.0, 0.15))
syn(drv, rgf, True, "0.7", (0.05, -0.18))
syn(pos, rge, True, "0.8")
syn(rge, rgf, False, "4.0", (0.16, 0.05))     # half-center mutual inhibition
syn(rgf, rge, False, "4.0", (0.16, 0.05))
syn(rge, ade, True, "2.5", (-0.14, 0.0))
syn(ade, rge, False, "2.5", (-0.24, 0.0))
syn(rgf, adf, True, "2.5", (-0.14, 0.0))
syn(adf, rgf, False, "2.5", (-0.24, 0.0))
ax.text(3.0, 4.35, "frequency rises with DRIVE;\nstance-biased drive sets "
        "duty", ha="center", fontsize=7.5, style="italic", color="0.3")

# ---------------- pattern formation (4 phase groups + adaptation) ----------------
pfe1 = neuron(6.0, 7.6, "PF-E1", "early stance")
pfe2 = neuron(6.0, 6.4, "PF-E2", "push-off")
pff1 = neuron(6.0, 5.0, "PF-F1", "early swing")
pff2 = neuron(6.0, 3.8, "PF-F2", "late swing")

syn(rge, pfe1, True, "1.2")
syn(rge, pfe2, True, "1.2")
syn(rgf, pff1, True, "1.2")
syn(rgf, pff2, True, "1.2")
syn(drv, pfe1, True, "0.2", (-0.05, 0.2))
syn(drv, pfe2, True, "0.2")
for px, py, nm in ((6.9, 7.6, "E1"), (6.9, 6.4, "E2"), (6.9, 5.0, "F1"),
                   (6.9, 3.8, "F2")):
    pfa = neuron(px, py, f"PFA-{nm}", "adapt")
    syn((6.0, py), pfa, True, "1.5", (0.0, 0.12))
    syn(pfa, (6.0, py), False, "1.5", (0.0, -0.14))
for a, b in ((pfe2, pff1), (pff1, pfe1), (pfe2, pff2)):
    syn(a, b, False, "3.0", (0.0, 0.1))
    syn(b, a, False)
ax.text(6.45, 2.9, "PF groups = phase windows within the step cycle;\n"
        "(tau, adaptation) shapes stagger their bursts",
        ha="center", fontsize=7.5, style="italic", color="0.3")

# ---------------- motoneuron pools ----------------
mns = {
    "hip_ext":  box(10.0, 7.4, 1.35, 0.55, "MN hip ext"),
    "knee_ext": box(10.0, 6.5, 1.35, 0.55, "MN knee ext"),
    "ankle_pf": box(10.0, 5.6, 1.35, 0.55, "MN ankle PF"),
    "hip_abd":  box(10.0, 4.7, 1.35, 0.55, "MN hip abd"),
    "hip_flex": box(12.0, 7.4, 1.35, 0.55, "MN hip flex"),
    "knee_flex": box(12.0, 6.5, 1.35, 0.55, "MN knee flex"),
    "ankle_df": box(12.0, 5.6, 1.35, 0.55, "MN ankle DF"),
}
# PF -> MN weights (primary group; biarticular muscles add half the
# secondary group's weight) - table from params.W_PF_MN
W = {
    ("PF-E1", "hip_ext"): 0.35, ("PF-E1", "knee_ext"): 0.25,
    ("PF-E1", "hip_abd"): 0.35,
    ("PF-E2", "ankle_pf"): 0.70, ("PF-E2", "hip_ext"): 0.45,
    ("PF-E2", "knee_ext"): 0.30, ("PF-E2", "hip_abd"): 0.30,
    ("PF-F1", "hip_flex"): 0.60, ("PF-F1", "knee_flex"): 1.20,
    ("PF-F1", "ankle_df"): 0.30,
    ("PF-F2", "ankle_df"): 0.40,
}
srcs = {"PF-E1": pfe1, "PF-E2": pfe2, "PF-F1": pff1, "PF-F2": pff2}
for (s, g), w in W.items():
    syn(srcs[s], mns[g], True, f"{w:.2f}", (0.06, 0.07), lw=1.1)
# POSTURE + BAL converge on MNs (weights per params.W_POSTURE / BAL)
syn(pos, mns["ankle_pf"], True, "posture", (-0.28, 0.14))
syn(balp, mns["ankle_pf"], True, "BAL")
syn(bald, mns["ankle_df"], True, "BAL")

# ---------------- proprioceptive reflex paths ----------------
ia = neuron(10.0, 3.3, "Ia", "velocity")
ii = neuron(11.2, 3.3, "II", "length")
ib = neuron(12.4, 3.3, "Ib", "force")
ax.text(11.2, 2.55, "afferent encoders: tendon velocity / length / force\n"
        "(gains speed-modulated + stance-gated in the runner)",
        ha="center", fontsize=7.5, style="italic", color="0.3")
syn(ia, mns["knee_ext"], True, "homonymous")
syn(ia, mns["knee_flex"], False, "reciprocal")
syn(ii, mns["knee_ext"], True)
syn(ib, mns["ankle_pf"], False, "autogenic")
ibexc = neuron(9.0, 2.9, "IB-EXC", "load sharing")
syn(ib, ibexc, True, "0.5")
syn((3.0, 6.8), ibexc, True, "RG-E gate", (0.1, -0.28))
syn(ibexc, mns["ankle_pf"], True, "stance reflex reversal")

# ---------------- legend ----------------
ly = 1.2
ax.add_patch(Circle((1.0, ly), R, fc="white", ec="black", lw=1.6))
ax.text(1.45, ly, "neuron", fontsize=8.5, va="center")
ax.add_patch(FancyArrowPatch((2.8, ly), (3.9, ly), arrowstyle="-", lw=1.5))
ax.add_patch(Polygon([(3.9, ly), (3.74, ly + 0.075), (3.74, ly - 0.075)],
                     fc="white", ec="black"))
ax.text(4.15, ly, "excitatory synapse", fontsize=8.5, va="center")
ax.add_patch(FancyArrowPatch((5.8, ly), (6.9, ly), arrowstyle="-", lw=1.5))
ax.add_patch(Circle((6.9, ly), 0.085, fc="black"))
ax.text(7.15, ly, "inhibitory synapse", fontsize=8.5, va="center")
ax.text(9.0, ly, "numbers = synaptic conductances (params.py: G, W_PF_MN)",
        fontsize=8.5, va="center", style="italic", color="0.3")

fig.savefig("spinal_circuit.png", dpi=150, bbox_inches="tight")
print("saved spinal_circuit.png")
