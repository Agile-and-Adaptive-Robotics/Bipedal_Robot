"""Circuit-rule figure, Di Russo Fig-3/4 style (Ben 2026-09-21).

v4 layout contract (after three visual-judge passes): panel = 10x10
axis units ~ 3.06 in tall. Title band y 9.3..10. Diagram zone
y 3.3..9.2 (all labels stay inside it). Dashed separator at y 3.05.
Caption strip: exactly 4 lines max, 10 pt Arial at 0.59-unit leading
(y 2.70, 2.11, 1.52, 0.93) - never overlaps, never crosses the frame.

Standards (AGENTS.md "Figure standards", 2026-09-21): 7.5x10 in,
Arial >=10 pt, no italics, Colors.m (Paul Tol) palette, colorblind-
safe shapes (open circle = neuron, open triangle = exc synapse,
filled dot = inh synapse, ellipse = muscle, rectangle = port).
Alt text beside the figure.
"""
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle

HERE = Path(__file__).parent
FIGS = HERE / "figures"
DISS = (HERE.parents[2] / "Documentation" / "Reports and Papers"
        / "Dissertation" / "CPG_airstepping_figs")

GOLD = "#FFD700"
ORANGE = "#FFB14E"
CORAL = "#FA8775"
PINK = "#EA5F94"
MAGENTA = "#CD34B5"
MAG2 = "#9D02D7"
INDIGO = "#0000FF"
GREY = "#9A9A9A"
LAV = "#E8E4F8"
INK = "#1A1A1A"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 10, "font.style": "normal",
    "axes.linewidth": 0.8,
})

NR = 0.46
WRAP = 44
CAP_Y = (2.70, 2.11, 1.52, 0.93)


def neuron(ax, x, y, label, color, bold=False, below=False):
    ax.add_patch(Circle((x, y), NR, facecolor="white",
                        edgecolor=color, linewidth=1.8, zorder=5))
    if below:
        ax.text(x, y - 0.70, label, ha="center", va="top",
                fontsize=10, color=INK,
                fontweight="bold" if bold else "normal",
                fontstyle="normal", zorder=8)
    else:
        ax.text(x, y + 0.54, label, ha="center", va="bottom",
                fontsize=10, color=INK,
                fontweight="bold" if bold else "normal",
                fontstyle="normal", zorder=8)


def mndot(ax, x, y, label, above=True):
    ax.add_patch(Circle((x, y), NR, facecolor="white",
                        edgecolor=INDIGO, linewidth=2.4, zorder=5))
    if above:
        ax.text(x, y + 0.60, label, ha="center", va="bottom",
                fontsize=10, color=INK, fontweight="bold",
                fontstyle="normal", zorder=8)
    else:
        ax.text(x, y - 0.68, label, ha="center", va="top",
                fontsize=10, color=INK, fontweight="bold",
                fontstyle="normal", zorder=8)


def muscle(ax, x, y, label, color=ORANGE, below=False):
    ax.add_patch(Ellipse((x, y), 1.6, 0.62, facecolor=LAV,
                         edgecolor=color, linewidth=1.8, zorder=5))
    if below:
        ax.text(x, y - 0.62, label, ha="center", va="top",
                fontsize=10, color=INK, fontstyle="normal", zorder=8)
    else:
        ax.text(x, y + 0.75, label, ha="center", va="bottom",
                fontsize=10, color=INK, fontstyle="normal", zorder=8)


def port(ax, x, y, label, w=2.1):
    ax.add_patch(Rectangle((x - w / 2, y - 0.30), w, 0.60,
                           facecolor=LAV, edgecolor=ORANGE,
                           linewidth=1.6, zorder=5))
    ax.text(x, y, label, ha="center", va="center", fontsize=10,
            color=INK, fontstyle="normal", zorder=8)


def connect(ax, p1, p2, color, kind, t=0.55, lw=1.8, ls="-"):
    x1, y1 = p1
    x2, y2 = p2
    tx, ty = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
    ax.plot([x1, x2], [y1, y2], color=color, lw=lw, ls=ls, zorder=4)
    if kind == "exc":
        ax.add_patch(Polygon([(tx, ty + 0.32),
                              (tx - 0.26, ty - 0.22),
                              (tx + 0.26, ty - 0.22)], closed=True,
                             facecolor="white", edgecolor=color,
                             linewidth=1.5, zorder=7))
    else:
        ax.add_patch(Circle((tx, ty), 0.18, facecolor=INK,
                            edgecolor=color, linewidth=1.4, zorder=7))
    ax.annotate("", xy=(x2, y2), xytext=(tx, ty),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw),
                zorder=7)


def mark(ax, x, y, num):
    ax.add_patch(Circle((x, y), 0.26, facecolor="white",
                        edgecolor=MAGENTA, linewidth=1.6, zorder=9))
    ax.text(x, y, str(num), ha="center", va="center", fontsize=10,
            color=MAGENTA, fontweight="bold", fontstyle="normal",
            zorder=10)


def caption(ax, notes):
    lines = []
    for num, txt in notes:
        for ln in textwrap.wrap(txt, WRAP):
            lines.append((num, ln))
    if len(lines) > 4:
        raise ValueError(f"caption overflow: {len(lines)} lines > 4")
    for (num, ln), y in zip(lines, CAP_Y):
        ax.text(0.25, y, ln, ha="left", va="center", fontsize=10,
                color=MAGENTA if num else INK, fontstyle="normal")


def panel(ax, title):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#777777")
    ax.axhline(3.05, color="#BBBBBB", lw=0.8, ls=(0, (4, 3)))
    ax.text(5, 9.62, title, ha="center", va="center", fontsize=11,
            fontweight="bold", fontstyle="normal", color=INK)


fig = plt.figure(figsize=(7.5, 10))
fig.suptitle("Spinal circuit rules: implementation vs literature "
             "(Ia / II / Ib with their interneurons)",
             fontsize=12, fontweight="bold", y=0.985)
gs = fig.add_gridspec(3, 2, left=0.03, right=0.985, top=0.945,
                      bottom=0.015, hspace=0.30, wspace=0.10)
AX = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(2)]

# ============================ A: architecture ================
a = AX[0]
panel(a, "A  Architecture, one side")
port(a, 1.05, 7.7, "DRIVE", w=1.7)
port(a, 1.05, 6.7, "POSTURE", w=1.9)
neuron(a, 3.7, 8.3, "RG-E", MAG2, bold=True)
neuron(a, 3.7, 6.5, "RG-F", MAG2, bold=True, below=True)
neuron(a, 2.7, 7.55, "InE", CORAL, below=True)
neuron(a, 4.7, 7.55, "InF", CORAL, below=True)
connect(a, (3.7, 7.84), (2.7, 7.85), GREY, "exc", t=0.5)
connect(a, (3.7, 6.96), (4.7, 7.25), GREY, "exc", t=0.5)
connect(a, (2.7, 7.09), (3.7, 6.04), GREY, "inh", t=0.5)
connect(a, (4.7, 7.09), (3.7, 7.76), GREY, "inh", t=0.5)
connect(a, (1.9, 7.75), (3.25, 8.22), GREY, "exc", t=0.5)
# POSTURE wire omitted for clarity (feeds all MNs; see caption)
for x, nm in ((5.9, "E1"), (6.9, "E2"), (7.9, "F1"), (8.9, "F2")):
    neuron(a, x, 7.9, nm, MAG2)
connect(a, (4.16, 8.3), (5.44, 7.95), GREY, "exc", t=0.4)
connect(a, (4.16, 6.6), (7.44, 7.76), GREY, "exc", t=0.5)
a.add_patch(Circle((8.9, 5.9), NR, facecolor="white",
                   edgecolor=INDIGO, linewidth=2.4, zorder=5))
a.text(8.2, 5.9, "MN", ha="right", va="center", fontsize=10,
       color=INK, fontweight="bold", fontstyle="normal", zorder=8)
connect(a, (8.9, 7.44), (8.9, 6.36), GREY, "exc", t=0.5)
muscle(a, 8.9, 4.4, "muscle", below=True)
connect(a, (8.9, 5.44), (8.9, 4.71), GREY, "exc", t=0.45)
for i, nm in enumerate(("Ia", "II", "Ib")):
    neuron(a, 3.7 + 1.1 * i, 4.0, nm, GOLD, bold=True)
connect(a, (8.15, 4.2), (4.85, 4.05), GOLD, "exc", t=0.55)
caption(a, [
    (0, "RG: persistent-Na half-centers, mutual"),
    (0, "inhibition via InE/InF. PF: E1 early/E2"),
    (0, "late stance, F1 early/F2 late swing."),
    (0, "Joint PF variant OFF. c1/V3 at midline."),
])

# ============================ B: Ia ===========================
a = AX[1]
panel(a, "B  Ia: monosynaptic + reciprocal")
muscle(a, 1.6, 8.0, "agonist")
neuron(a, 1.6, 6.3, "Ia", GOLD, bold=True, below=True)
connect(a, (1.6, 7.69), (1.6, 6.76), GOLD, "exc", t=0.6)
mndot(a, 5.0, 6.3, "MN ag")
connect(a, (2.06, 6.3), (4.54, 6.3), GOLD, "exc", t=0.5)
neuron(a, 2.9, 4.4, "IaIN ag", CORAL, below=True)
connect(a, (1.6, 5.84), (2.7, 4.86), GOLD, "exc", t=0.55)
neuron(a, 6.3, 4.4, "IaIN ant", CORAL, below=True)
# IaIN <-> antagonist IaIN mutual inhibition (two opposing arrows)
connect(a, (3.36, 4.4), (4.82, 4.4), CORAL, "inh", t=0.5)
connect(a, (5.84, 4.4), (4.36, 4.4), CORAL, "inh", t=0.5)
mark(a, 4.6, 3.55, 2)
muscle(a, 8.6, 8.0, "antagonist", GREY)
mndot(a, 8.6, 6.3, "MN ant", above=False)
connect(a, (8.6, 7.69), (8.6, 6.76), GREY, "exc", t=0.6)
connect(a, (3.36, 4.55), (8.35, 5.95), CORAL, "inh", t=0.62)
connect(a, (5.84, 4.55), (5.15, 5.8), CORAL, "inh", t=0.55)
neuron(a, 1.6, 4.7, "PF F1", MAG2, below=True)
connect(a, (2.06, 4.75), (2.6, 4.7), MAG2, "exc", t=0.45)
caption(a, [
    (0, "Ia to MN: monosynaptic excitation. Ia to"),
    (0, "IaIN to antagonist MN: reciprocal inh."),
    (2, "IaIN-IaIN mutual inhibition: absent."),
    (1, "ia_in = 0 (untuned): DIRECT Ia to MN ant."),
])

# ============================ C: II ===========================
a = AX[2]
panel(a, "C  II: length feedback rule")
muscle(a, 1.6, 8.0, "agonist")
neuron(a, 1.6, 6.3, "II", GOLD, bold=True, below=True)
connect(a, (1.6, 7.69), (1.6, 6.76), GOLD, "exc", t=0.6)
neuron(a, 3.5, 5.3, "", CORAL)
a.text(4.15, 5.3, "IN exc", ha="left", va="center", fontsize=10,
       color=INK, fontstyle="normal", zorder=8)
connect(a, (1.98, 6.0), (3.06, 5.6), GOLD, "exc", t=0.5)
neuron(a, 3.5, 4.4, "IN inh", CORAL, below=True)
connect(a, (1.6, 5.84), (3.06, 4.7), GOLD, "exc", t=0.5)
mndot(a, 6.0, 6.3, "MN ag")
connect(a, (3.96, 5.45), (5.7, 6.1), CORAL, "exc", t=0.55)
muscle(a, 9.0, 8.0, "antagonist", GREY)
mndot(a, 9.0, 6.3, "MN ant", above=False)
connect(a, (9.0, 7.69), (9.0, 6.76), GREY, "exc", t=0.6)
connect(a, (3.96, 4.6), (8.75, 5.9), CORAL, "inh", t=0.62)
mark(a, 6.5, 5.2, 2)
# ours: direct II -> agonist MN (dashed magenta shortcut)
connect(a, (2.06, 6.3), (5.54, 6.3), MAGENTA, "exc", t=0.42,
        ls=(0, (3, 3)))
mark(a, 3.7, 6.78, 1)
caption(a, [
    (0, "Literature rule: II drives two collateral"),
    (0, "INs - IN exc to agonist, IN inh to antag."),
    (1, "OURS: direct II to agonist MN (dashed)."),
    (2, "OURS: both IN interneurons absent."),
])


# ============================ D: Ib ===========================
a = AX[3]
panel(a, "D  Ib: force feedback + reversal")
muscle(a, 1.5, 8.1, "extensor")
neuron(a, 1.5, 6.4, "Ib", GOLD, bold=True, below=True)
connect(a, (1.5, 7.79), (1.5, 6.86), GOLD, "exc", t=0.6)
mndot(a, 4.6, 6.4, "MN ext")
connect(a, (1.96, 6.4), (4.14, 6.4), MAGENTA, "inh", t=0.38)
mark(a, 3.0, 6.88, 1)
neuron(a, 4.6, 4.3, "IBEXC", CORAL, below=True)
connect(a, (1.5, 5.94), (4.35, 4.7), GOLD, "exc", t=0.5)
connect(a, (4.9, 4.7), (4.75, 5.9), CORAL, "exc", t=0.55)
neuron(a, 8.9, 5.5, "RG-E", MAG2, bold=True)
connect(a, (8.44, 5.4), (5.05, 4.42), MAG2, "exc", t=0.42)
neuron(a, 7.0, 4.4, "Ib load IN", PINK, below=True)
connect(a, (5.0, 4.2), (6.55, 4.4), PINK, "exc", t=0.5)
connect(a, (7.45, 4.75), (8.8, 5.05), PINK, "exc", t=0.55)
caption(a, [
    (1, "DIRECT autogenic inh. (lit.: disynaptic"),
    (0, "via Ib INs, which also mutually inhibit"),
    (0, "antagonists - absent). IBEXC drives the"),
    (0, "Ib load IN (LBIN) to RG-E. Ib direct E."),
])

# ============================ E: Renshaw ======================
a = AX[4]
panel(a, "E  Renshaw: recurrent inhibition")
mndot(a, 2.7, 7.2, "MN")
neuron(a, 5.8, 7.2, "RC", CORAL, bold=True)
connect(a, (3.16, 7.2), (5.34, 7.2), CORAL, "exc", t=0.45)
connect(a, (5.3, 6.9), (3.05, 6.85), CORAL, "inh", t=0.6)
neuron(a, 8.6, 5.2, "RC2", CORAL, below=True)
connect(a, (6.22, 6.95), (8.3, 5.6), CORAL, "inh", t=0.55)
connect(a, (8.3, 5.7), (6.2, 6.88), CORAL, "inh", t=0.45)
neuron(a, 5.8, 4.5, "IaIN", CORAL, below=True)
connect(a, (5.8, 6.74), (5.62, 4.96), CORAL, "inh", t=0.55)
caption(a, [
    (0, "MN collateral excites RC; RC inhibits"),
    (0, "its own MN. RC and RC2 mutual inh. RC"),
    (0, "inhibits IaIN: disinhibition (Hultborn)."),
    (0, "Matches rule 5; needs renshaw gain > 0."),
])

# ============================ F: contact + KINH ===============
a = AX[5]
panel(a, "F  Foot sensors, KINH, crossed")
port(a, 1.15, 7.9, "HEEL_c", w=1.7)
neuron(a, 2.6, 7.9, "heel SN", GOLD, below=True)
connect(a, (2.0, 7.9), (2.14, 7.9), ORANGE, "exc", t=0.2)
neuron(a, 4.6, 7.9, "heel IN", CORAL, below=True)
connect(a, (3.06, 7.9), (4.14, 7.9), GOLD, "exc", t=0.5)
neuron(a, 8.3, 8.0, "RG-E", MAG2, bold=True)
neuron(a, 8.3, 6.7, "RG-F", MAG2, bold=True)
connect(a, (5.06, 7.95), (7.84, 7.98), CORAL, "exc", t=0.55)
connect(a, (5.06, 7.8), (7.84, 6.9), CORAL, "inh", t=0.62)
# contralateral chain: SN -> IN -> KINH
neuron(a, 2.6, 6.2, "heel SN", GOLD, below=True)
a.plot([1.25, 1.7, 2.14], [7.6, 6.6, 6.2], color=ORANGE, lw=1.8,
       ls=(0, (2, 2)), zorder=4)
neuron(a, 4.6, 5.6, "", CORAL)
a.text(5.1, 6.28, "IN (contra)", ha="left", va="center",
       fontsize=10, color=INK, fontstyle="normal", zorder=8)
connect(a, (3.06, 6.15), (4.14, 5.75), GOLD, "exc", t=0.5)
neuron(a, 6.6, 5.0, "KINH", PINK, bold=True)
connect(a, (5.06, 5.5), (6.14, 5.05), PINK, "exc", t=0.5)
mark(a, 4.2, 5.1, 3)
ax_mn = a.add_patch(Circle((8.9, 4.0), NR, facecolor="white",
                           edgecolor=INDIGO, linewidth=2.4, zorder=5))
a.text(8.9, 4.62, "knee ext +\nankle PF MN", ha="center",
       va="bottom", fontsize=10, color=INK, fontweight="bold",
       fontstyle="normal", zorder=8)
connect(a, (7.06, 4.7), (8.7, 4.25), PINK, "inh", t=0.55)
neuron(a, 1.6, 4.4, "PF F1", MAG2, below=True)
connect(a, (2.06, 4.35), (6.35, 4.85), MAG2, "exc", t=0.3)
caption(a, [
    (0, "Heel SN to IN: RG-E exc, RG-F inh (reset-"),
    (0, "to-extension). SN = runner-side encoder."),
    (3, "Crossed drive (contra_kinh): opposite"),
    (0, "strike forces swing. KINH gated by PF F1."),
])

for out in (FIGS, DISS):
    fig.savefig(out / "circuit_rules.png", dpi=300)
    fig.savefig(out / "circuit_rules.pdf")
plt.close(fig)

ALT = """Alt text - circuit_rules figure (spinal circuit rules:
implementation vs literature).

Panel A - Architecture. Signals flow left to right: descending DRIVE
and POSTURE inputs excite the rhythm generator, a pair of
persistent-sodium half-center neurons RG-E and RG-F whose mutual
inhibition is routed through interneurons InE and InF. The rhythm
generator drives four phase-window pattern-formation cells per side:
E1 early stance, E2 late stance, F1 early swing, F2 late swing, which
excite the motoneuron pool of the muscle. The muscle feeds three
somatosensory afferent neurons, Ia, II, and Ib, back into the cord. A
joint-layer pattern-formation variant (hip, knee, and ankle
half-center pairs) exists but is default-off. Commissural
interneurons c1 and V3 couple the two sides at the midline.

Panel B - Ia pathways. The agonist muscle excites its Ia afferent,
which monosynaptically excites the agonist motoneuron (matching the
literature) and also excites the agonist Ia inhibitory interneuron,
which inhibits the antagonist motoneuron (disynaptic reciprocal
inhibition, matching); symmetrically, the antagonist Ia interneuron
inhibits the agonist motoneuron. The two Ia interneurons mutually
inhibit each other in the literature rule (marker 2) - NOT
implemented in our circuit. The Ia interneuron is phase-gated by the
PF F1 swing window and disinhibited by the Renshaw cell (panel E).
Note 1: at gain ia_in = 0 (untuned configurations only) Ia connects
directly to the antagonist motoneuron.

Panel C - II pathways, literature rule. The agonist muscle excites
its II afferent, which drives TWO collateral interneurons: an
excitatory interneuron to the agonist motoneuron and an inhibitory
interneuron to the antagonist motoneuron. In our implementation
(marker 1) the II afferent excites the agonist motoneuron DIRECTLY
(shown dashed), and (marker 2) both interneurons of the rule are
absent. II afferents also project directly to flexor pattern and
rhythm cells.

Panel D - Ib pathways. The extensor muscle excites its Ib afferent.
Autogenic inhibition onto its own motoneuron is direct (marker 1);
the literature rule is disynaptic via Ib interneurons, which also
mutually inhibit antagonistic Ib interneurons (absent in ours). The
positive force reversal is interneuron-mediated: Ib excites the
stance-gated IBEXC interneuron, which excites the extensor
motoneuron and the stance-group Ib interneuron (code symbol LBIN),
which excites RG-E and prolongs stance. Ib additionally projects
directly to the extensor centers.

Panel E - Renshaw cells. The motoneuron collateral excites its
Renshaw cell, which inhibits that same motoneuron, is mutually
inhibited with the antagonist-pool Renshaw cell RC2, and inhibits the
Ia inhibitory interneuron (recurrent disinhibition, Hultborn 1971).
Matches Di Russo rule 5; the topology is conditional on the renshaw
gain.

Panel F - Foot mechanosensors, KINH, and crossed edges. The heel
contact port (a runner-side transducer/encoder standing in for the
somatosensory neuron) feeds the heel somatosensory neuron, which
excites the heel interneuron; the interneuron excites RG-E and
inhibits RG-F (reset-to-extension, Conway 1987). The contralateral
chain mirrors it: heel SN (contra) to heel IN (contra), which excites
KINH (marker 3), an inhibitory interneuron suppressing the
knee-extensor and ankle-plantarflexor motoneurons of this side: the
opposite heel strike forces this leg's stance-to-swing transition
(contra_kinh, added 2026-09-21). KINH is gated by the own-side PF F1
swing window. Runner-side mechanisms, not SNS topology:
contact_onset adds current kicks at foot loading and unloading edges;
contra_swing applies a crossed kick at the rhythm generator. The toe
chain follows the heel pattern.
"""
(FIGS / "circuit_rules_alt.txt").write_text(ALT, encoding="utf-8")
(DISS / "circuit_rules_alt.txt").write_text(ALT, encoding="utf-8")
print("saved circuit_rules.{png,pdf} + circuit_rules_alt.txt "
      "(figures + Dissertation CPG_airstepping_figs)")
