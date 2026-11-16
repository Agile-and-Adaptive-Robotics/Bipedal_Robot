"""Show the governing tibial clearance samples from the saved ankle run."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle
import numpy as np

from fr_plotDistalStages import GOLD, INDIGO
from fr_plotPilotAnatomy import mesh_data


def draw(folder, tibia_file):
    v, faces = mesh_data(tibia_file)
    triangles = 1000*v[faces]
    rows = json.loads((folder/"CLEARANCE_WITNESSES.json").read_text())
    fig,axes=plt.subplots(len(rows),2,figsize=(8.5,2.7*len(rows)),layout="constrained")
    for pair,r in zip(axes,rows):
        use_lead = "lead" in r and r["lead"]["margin_m"] < r["bone"]["margin_m"]
        d=r["lead"] if use_lead else r["bone"]
        part="free tendon lead" if use_lead else d["part"]
        pose2="MTP" if r["second"].startswith("mtp") else "subtalar"
        at=1000*np.asarray(d["closest_m"])
        segment=1000*np.asarray(d["segment_m"])
        for ax,dim,direction in zip(pair,(0,2),("Anterior +x","Right +z")):
            ax.add_collection(PolyCollection(triangles[:,:,[dim,1]],
                                             facecolor="#a8b4c1",edgecolor="none",alpha=.28))
            ax.plot(segment[:,dim],segment[:,1],"--",color=INDIGO,lw=1.25,label="Checked centerline")
            if not use_lead:
                sleeve=1000*np.asarray(d["sleeve_m"])
                ax.plot(sleeve[:,dim],sleeve[:,1],color=INDIGO,lw=2.1,label="Inflatable vessel axis")
                ax.add_patch(Circle((at[dim],at[1]),1000*d["radius_m"],
                                    fill=False,lw=1.4,color="#443827",linestyle=":"))
            color=GOLD if d["margin_m"]<0 else INDIGO
            ax.scatter(at[dim],at[1],s=55,marker="D",facecolor=color,
                       edgecolor="black",zorder=5,label="Governing sample")
            ax.set_xlim(at[dim]-52,at[dim]+52)
            ax.set_ylim(at[1]-63,at[1]+63)
            ax.set_aspect("equal",adjustable="box")
            ax.set_xlabel(f"{direction} (mm)")
            ax.grid(alpha=.18)
        pair[0].set_ylabel("Superior +y (mm)")
        pair[0].set_title(f"{r['name'].replace('_',' ')}\n{part} · {1000*d['margin_m']:+.1f} mm margin",
                          fontsize=9)
        pair[1].set_title(f"Ankle {d['ankle_deg']:+.0f}° · {pose2} {d['second_deg']:+.0f}°",
                          fontsize=9)
    fig.suptitle("Ankle packaging · governing tibial clearance samples\n"
                 "Gold diamond: shortfall · Indigo diamond: clearance passes")
    fig.savefig(folder/"ankle_clearance_witnesses.png",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    draw(Path(sys.argv[1]),Path(sys.argv[2]))
