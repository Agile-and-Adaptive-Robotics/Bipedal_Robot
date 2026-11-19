"""Reflow the defense toolchain diagram at dissertation column width."""
from pathlib import Path
import json
import hashlib
import shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

here=Path(__file__).resolve().parent
diss=here.parents[1]
out=diss/"Figures/95-AppendixD"
components=out/"Components/research_toolchain"
components.mkdir(exist_ok=True,parents=True)
source=Path("D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/RESEARCH_TOOLCHAIN_WORKFLOW.svg")
packed=components/source.name
if source.is_file():
    if packed.exists():
        assert source.read_bytes()==packed.read_bytes()
    else:
        shutil.copyfile(source,packed)
assert packed.is_file()
plt.rcParams.update({"font.family":"Arial","font.size":10.6,"pdf.fonttype":42})
fig,ax=plt.subplots(figsize=(6.35,3.18))
fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
ax.set(xlim=(0,6.35),ylim=(0,3.18))
ax.axis("off")
def rect(x,y,w,h,face,edge,ls="-"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.028,rounding_size=0.065",
                               facecolor=face,edgecolor=edge,lw=1.2,linestyle=ls))
xs=[.06,2.29,4.52]
w,h=1.78,.85
rows=[1.57,.29]
steps=[(0,0,"Source model","OpenSim, robot CAD,\nor an AnimatLab model"),
       (0,1,"Prepare mechanics","Check geometry, muscle\npaths, inertias, and contact"),
       (0,2,"Assemble control","Connect the neural circuit,\nactuators, and sensors"),
       (1,2,"Execute a trial","Run MuJoCo or a supported\nSimulink simulation"),
       (1,1,"Run comparisons","Tune the controller or\nperturb a saved trial"),
       (1,0,"Inspect and export","Read the neural signals,\nforces, and joint motion")]
for r,c,title,body in steps:
    x,y=xs[c],rows[r]
    rect(x,y,w,h,"#f3f6f9" if r+c else "#f5f5f5","#57758e" if r+c else "#777")
    ax.text(x+w/2,y+.60,title,ha="center",va="center",weight="bold",fontsize=11)
    ax.text(x+w/2,y+.28,body,ha="center",va="center",fontsize=10.0,linespacing=1.22)
for (x,y),(u,v) in [((xs[0]+w+.04,rows[0]+h/2),(xs[1]-.03,rows[0]+h/2)),
                      ((xs[1]+w+.04,rows[0]+h/2),(xs[2]-.03,rows[0]+h/2)),
                      ((xs[2]+w/2,rows[0]-.035),(xs[2]+w/2,rows[1]+h+.035)),
                      ((xs[2]-.04,rows[1]+h/2),(xs[1]+w+.03,rows[1]+h/2)),
                      ((xs[1]-.04,rows[1]+h/2),(xs[0]+w+.03,rows[1]+h/2))]:
    ax.add_patch(FancyArrowPatch((x,y),(u,v),arrowstyle="-|>",mutation_scale=12,
                                 color="#425869",linewidth=1.5,shrinkA=0,shrinkB=0))
rect(.06,2.66,6.24,.38,"#fff9e8","#96783e","--")
ax.text(3.18,2.85,"Proposed interface for selecting, configuring, and launching experiments",
        ha="center",va="center",fontsize=10.4)
fig.savefig(out/"research_toolchain_workflow.pdf")
plt.close(fig)
(out/"research_toolchain_workflow.alt.txt").write_text("Six connected stages lead from a source model through mechanical preparation, controller assembly, trial execution, comparison, and inspection of neural signals, force, and joint motion. A dashed band identifies the proposed interface for selecting and launching experiments across these stages. The circuit editor already supports neural wiring.\n",encoding="utf-8")
(here/"WORKFLOW_SOURCE.json").write_text(json.dumps({"original":str(source),"repository_copy":str(packed.relative_to(diss)),
    "source_sha256":hashlib.sha256(packed.read_bytes()).hexdigest(),"source_implementation":"D:/Users/Ben Bolen/ToolchainDefense_20261102/source/draw_defense_workflow.py",
    "interpretation":"Six stages from the defense SVG. Solid boxes show steps supported by separate tools and scripts. Dashed enclosure identifies the proposed cross-stage launcher. Figure 6.2 retains the two mechanical preparation routes."},indent=2)+"\n",encoding="utf-8")
