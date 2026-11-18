"""Generate the paired, conceptual BPA end-region and loading diagrams at print size."""
import hashlib
import io
import json
import shutil
import tempfile
import zipfile
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

notes = Path(__file__).resolve().parent
assets = notes.parents[1] / "Figures" / "70-overall"
archive = assets / "Deprecated" / "20261117"
archive.mkdir(exist_ok=True)
names = ["bpa_construction", "bpa_exploded_cutaway"]
original_hashes = {
    "bpa_construction": "5a069309f4752d04fa6c3a429792159c967e1aa57d4c1d85e6e63081cab6b7bf",
    "bpa_exploded_cutaway": "f7075e81d7f3c4fdb62959a7743718d669b3afe8d7eb7c029f41e62fc16b7da0",
}
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_changed(path, data):
    if path.exists() and path.read_bytes() == data:
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as output:
        temporary = Path(output.name)
        output.write(data)
    temporary.replace(path)

bundle = notes / "SUPERSEDED_FIGURE_SOURCES.zip"
if bundle.exists():
    with zipfile.ZipFile(bundle) as saved_sources:
        for member in saved_sources.namelist():
            p = (notes.parents[1] / member).resolve()
            assert p.is_relative_to((assets / "Deprecated").resolve())
            if not p.exists():
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_bytes(saved_sources.read(member))
for name in names:
    saved = archive / (name + ".pdf")
    if not saved.exists():
        assert digest(assets / saved.name) == original_hashes[name], "Restore the original figure bundle in the dated Notes directory."
        shutil.copy2(assets / saved.name, saved)
    assert digest(saved) == original_hashes[name], saved

mpl.rcParams.update({"font.family": "Arial", "font.size": 10.5,
                     "pdf.fonttype": 42, "lines.solid_capstyle": "round"})
ink, wall, fill = "#222222", "#225C79", "#D4E6ED"
def canvas(height, ymax):
    fig = plt.figure(figsize=(14 / 2.54, height / 2.54))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, 100), ylim=(0, ymax))
    ax.set_axis_off()
    return fig, ax

def save(fig, name, alt):
    output = io.BytesIO()
    fig.savefig(output, format="pdf", metadata={"CreationDate": None, "ModDate": None})
    write_changed(assets / (name + ".pdf"), output.getvalue())
    plt.close(fig)
    write_changed(assets / (name + ".alt.txt"), (alt + "\r\n").encode("utf-8"))

def leader(ax, text, position, point, align="center"):
    ax.annotate(text, xy=point, xytext=position, ha=align, va="center",
                arrowprops={"arrowstyle": "-", "linewidth": 0.8, "color": ink,
                            "shrinkA": 2, "shrinkB": 1}, color=ink)

fig, ax = canvas(5.5, 40)
ax.text(1, 37, "A   Axial diameter profile", fontsize=12, fontweight="bold")
# End fittings and sleeve placement follow the longitudinal construction schematic.
# The smooth transition illustrates the hypothesis and is dimensionless.
x = np.linspace(17, 83, 401)
u = np.clip(np.minimum(x - 17, 83 - x) / 15, 0, 1)
r = 2.8 + 2.5 * (1 - np.cos(np.pi * u))
center = 17.5
ax.fill_between(x, center - r, center + r, facecolor=fill, edgecolor=wall, lw=1.5)
for y in (center - 2.8, center + 2.8):
    ax.plot([17, 83], [y, y], color=ink, ls=(0, (4, 3)), lw=1)
for low in (12, 83):
    ax.add_patch(Rectangle((low, center - 4.5), 5, 9, fc="#BBBBBB", ec=ink, lw=1.1))
leader(ax, "End fitting", (10, 29), (14.5, center + 4.5))
leader(ax, "Transition", (32, 30), (26, center + 5.6))
leader(ax, "Central span", (60, 30), (59, center + 7.8))
ax.legend([Line2D([0], [0], color=wall, lw=1.6),
           Line2D([0], [0], color=ink, lw=1.1, ls=(0,(4,3)))],
          ["Pressurized profile", "Uninflated diameter"], loc="lower center", ncol=2,
          frameon=False, borderaxespad=0.2, handlelength=2.1, columnspacing=1.8)
save(fig, names[0], "Longitudinal schematic. The pressurized tube rises from its uninflated diameter beside each end fitting to a wider central span. A dashed pair of lines marks the uninflated diameter. The profiles illustrate a hypothesis.")

fig, ax = canvas(6.4, 48)
ax.text(1, 45, "B   Proposed loading conditions", fontsize=12, fontweight="bold")
for xmid, label in [(17, "Free contraction"), (50, "Suspended load"), (83, "Fixed length")]:
    ax.text(xmid, 39.5, label, ha="center", fontweight="bold")
    ax.plot([xmid - 6, xmid + 6], [36, 36], color=ink, lw=1.6)
    for offset in (-5, -2, 1, 4):
        ax.plot([xmid + offset, xmid + offset+1.8], [36, 37.6], color=ink, lw=0.8)
    ax.add_patch(Rectangle((xmid-1.8,32.5),3.6,3.2,fc="#BBBBBB",ec=ink,lw=0.8))
    y = np.linspace(18.5,32.5,101)
    u = np.clip(np.minimum(y-18.5,32.5-y)/3.5,0,1)
    width = 1.3 + 0.8*(1-np.cos(np.pi*u))
    ax.fill_betweenx(y,xmid-width,xmid+width,fc=fill,ec=wall,lw=1.1)
    ax.add_patch(Rectangle((xmid-1.8,15.5),3.6,3,fc="#BBBBBB",ec=ink,lw=0.8))
    ax.plot([xmid-3.1,xmid-3.1],[15,19.5], color=ink,lw=1)
    ax.plot([xmid+3.1,xmid+3.1],[15,19.5], color=ink,lw=1)
ax.annotate("",(23,25),(23,16),arrowprops={"arrowstyle":"<->","lw":1.2})
ax.text(17,6,"Fixture weight\ncounterbalanced",ha="center",va="center")
ax.plot([50,50],[15.5,12.3],color=ink,lw=1)
ax.add_patch(Rectangle((47.4,8.8),5.2,3.5,fc="#CCCCCC",ec=ink,lw=1))
ax.annotate("",(57,8.6),(57,15.5),arrowprops={"arrowstyle":"-|>","lw":1})
ax.text(50,3.6,"Known weight",ha="center")
ax.plot([81,85],[15.2,15.2],color=ink,lw=2)
ax.plot([83,83],[12.7,15.2],color=ink,lw=1)
ax.plot([77,89],[12.5,12.5],color=ink,lw=1.5)
for offset in (-5,-2,1,4):
    ax.plot([83+offset,84.8+offset],[10.8,12.5],color=ink,lw=0.8)
ax.text(83,5.8,"Measure axial\nreaction",ha="center",va="center")
save(fig, names[1], "All three proposed tests fix the top fitting and use the same gauge-pressure ramp. The lower fitting translates in an axial guide with its fixture weight counterbalanced in the first test and with a suspended mass in the second. Locking the lower fitting in the third test permits measurement of isometric reaction.")

sources = sorted(p for p in (assets / "Deprecated").rglob("*") if p.is_file() and
                 (p.name.lower().startswith("bpa") or p.name.lower().startswith("gen_bpa")))
record = {
    "print_width_cm": 14, "minimum_font_pt": 10.5,
    "geometry": "Fitting placement and actuator axis follow Deprecated/components/bpa_construction.tex. The normalized profiles and loading fixtures are proposed schematics, with no fitted dimensions or experimental samples.",
    "superseded": [{"original": str(assets / (n+".pdf")), "preserved": str(archive / (n+".pdf")),
                    "sha256": digest(archive / (n+".pdf"))} for n in names],
    "associated_preserved_files": [{"path": str(p.relative_to(assets)), "sha256": digest(p)} for p in sources],
    "outputs": [{"path": str(assets / (n+ext)), "sha256": digest(assets/(n+ext))}
                for n in names for ext in (".pdf",".alt.txt")]}
write_changed(notes / "FIGURE_PROVENANCE.json", (json.dumps(record,indent=2)+"\n").replace("\n","\r\n").encode("utf-8"))
archive_record={key:record[key] for key in ["superseded","associated_preserved_files"]}
write_changed(archive/"manifest.json",(json.dumps(archive_record,indent=2)+"\n").encode("utf-8"))
