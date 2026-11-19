"""Draw Chapter 7's controlled neuromechanical experiment at its intended print size."""
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Polygon

HERE = Path(__file__).resolve().parent
CHAPTER = HERE.parents[1]
NAME = "neuromechanical_experiment_workflow"
WIDTH, HEIGHT = 396, 405  # 5.50 by 5.625 inches, in PostScript points
INK = "#233340"
LINE = "#344F5F"
INPUT = "#F5F7F8"
SHARED = "#ECF2F5"
REPEAT = "#F8F9FA"
READOUT = "#E8F1EE"

# Fail visibly if the selected font cannot be embedded at its designed size.
font_manager.findfont("Arial", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 11.0,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "Bolen-neuromechanical-experiment-v1",
})

fig = plt.figure(figsize=(WIDTH/72, HEIGHT/72), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH), ylim=(HEIGHT, 0), frameon=False)
ax.set_axis_off()


def card(x, y, w, h, fill, title, lines, title_size=12, body_size=11):
    """Center separately editable title and body lines within a printed vector frame."""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=3",
        edgecolor=LINE, facecolor=fill, linewidth=1.0, zorder=2,
    ))
    body_leading = 13.0
    title_body = 15.4
    group_height = title_body + (len(lines) - 1)*body_leading
    first = y + (h - group_height)/2
    ax.text(x+w/2, first, title, fontsize=title_size, fontweight="bold",
            ha="center", va="center", color=INK, zorder=3)
    for i, line in enumerate(lines):
        ax.text(x+w/2, first+title_body+i*body_leading, line, fontsize=body_size,
                ha="center", va="center", color=INK, zorder=3)


def stroke(points):
    x, y = zip(*points)
    ax.add_line(Line2D(x, y, linewidth=1.6, color=LINE, solid_capstyle="butt",
                       solid_joinstyle="round", zorder=1))


def down_arrow(x, start, end):
    """A filled 9 by 8 point head points into the next operation."""
    head_length = 9
    head_width = 8
    if end-start < head_length:
        raise ValueError("The shaft must leave room for its arrowhead")
    stroke([(x,start), (x,end-head_length+.5)])
    ax.add_patch(Polygon([(x,end),(x-head_width/2,end-head_length),
                          (x+head_width/2,end-head_length)],closed=True,
                         facecolor=LINE, edgecolor="none", zorder=3))


left, mid, right = 94, WIDTH/2, 302
card(4, 4, 180, 58, INPUT, "Choose or build a controller",
     ["Specify motor outputs and the", "sensory pathways to test."])
card(212, 4, 180, 58, INPUT, "Select a verified body model",
     ["Gait2392, its rerouted robot", "counterpart, or a CAD model."])
for x in (left, right):
    down_arrow(x, 62, 80)

card(4, 80, 388, 43, SHARED, "Connect muscles and sensors",
     ["Map motor commands to muscles and sensory signals to circuit inputs."])
down_arrow(mid, 123, 140)

card(4, 140, 388, 61, SHARED, "Run and save the baseline",
     ["Set the initial pose, stimulus, descending drive, and support.",
      "Use the same settings for either repeat."])
stroke([(mid,201),(mid,214)])
stroke([(left,214),(right,214)])
for x in (left, right):
    down_arrow(x, 214, 232)

card(4, 232, 180, 67, REPEAT, "Sensory-pathway repeat",
     ["Change one pathway. Keep the", "model, actuator routes, and", "other neural settings fixed."])
card(212, 232, 180, 67, REPEAT, "Actuator-route repeat",
     ["Change one attachment or route.", "Keep the neural controller and", "its sensory settings fixed."])
ax.text(mid, 265.5, "or", ha="center", va="center", fontsize=11.0,
        color=INK, fontstyle="italic")

stroke([(left,299),(left,314),(right,314),(right,299)])
down_arrow(mid, 314, 332)

card(4, 332, 388, 69, READOUT, "Compare each repeat with the baseline",
     ["Align recordings at stimulus onset or the same foot-contact event.",
      "Compare neural activity, muscle force, joint motion, and foot loads.",
      "Identify the first signal that differs from baseline."])

fig.savefig(CHAPTER/(NAME+".pdf"), metadata={
    "Title": "Controlled neuromechanical experiments",
    "Creator": "make_figure.py",
    "CreationDate": datetime(2026,11,19,tzinfo=timezone.utc),
})
fig.savefig(HERE/(NAME+".svg"), metadata={
    "Title": "Controlled neuromechanical experiments", "Date": None,
    "Description": "Flow diagram for a baseline and two independent controlled repeats.",
})
plt.close(fig)
