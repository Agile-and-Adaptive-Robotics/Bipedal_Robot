"""Build defense-ready comparisons from the independently checked ankle layout."""
import csv
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image

from fr_exportPilotOpenSim import export_routes
from fr_plotAnkleClash import capsule
from fr_plotDistalStages import COLORS, GOLD, INDIGO
from fr_plotPilotAnatomy import mesh_data

COLORS_BY_GROUP = {
    "soleus": "#6b737e", "tib_post": "#394c68", "tib_ant": INDIGO,
    "per_brev": "#498788", "per_long": "#127d6d", "per_tert": "#755d9f",
    "flex_dig": "#485e44", "flex_hal": "#995c2c",
    "ext_hal": COLORS["magenta"], "ext_dig": "#563eb1",
}


def load_meshes(scene, geometry):
    model = ET.parse(scene["anatomy"])
    surfaces = []
    for body, T in scene["transforms"].items():
        T = np.asarray(T)
        for mesh in model.findall(f".//Body[@name='{body}']/attached_geometry/Mesh"):
            v, faces = mesh_data(geometry/Path(mesh.findtext("mesh_file")).name)
            v *= np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            v = 1000*(v@T[:3,:3].T+T[:3,3])
            surfaces.append(v[faces])
    return surfaces


def projections(folder, scene, surfaces, record):
    fig, axes = plt.subplots(1, 4, figsize=(13.3, 6.0), layout="constrained", sharey=True)
    for j, (key, heading) in enumerate((("before","Independent group searches"),
                                         ("after","Coordinated ankle layout"))):
        for dim, ax, xlabel in zip((0,2),axes[2*j:2*j+2],("Anterior (mm)","Right (mm)")):
            for tri in surfaces:
                ax.add_collection(PolyCollection(tri[:,:,[dim,1]],facecolor="#a9b1bc",
                                                 edgecolor="none",alpha=.3))
            for part in scene[key]:
                color = COLORS_BY_GROUP[part["group"]]
                for vessel in part["vessels"]:
                    v = np.asarray(vessel["points_mm"])[:,[dim,1]]
                    capsule(ax,*v,vessel["radius_mm"],color,.38)
                for h in part["hardware_mm"]:
                    capsule(ax,[h[dim],h[1]],[h[dim],h[1]],h[3],color,.27)
                v=np.asarray(part["wire_mm"])
                ax.plot(v[:,dim],v[:,1],color=color,lw=1.35)
                if "mount_reach_mm" in part:
                    w=np.asarray(part["mount_reach_mm"])
                    ax.plot(w[:,dim],w[:,1],":",color="#575e65",lw=.8,alpha=.6)
            ax.set(aspect="equal",xlim=(-160,225 if dim==0 else 185),ylim=(-555,55),xlabel=xlabel)
            ax.set_title(heading if dim==0 else "Frontal view",fontsize=12)
            ax.grid(alpha=.12)
    axes[0].set_ylabel("Superior (mm)")
    handles = [plt.Line2D([0],[0],color=COLORS_BY_GROUP[g],lw=3,label=g.replace("_"," "))
               for g in ("tib_ant","ext_hal","ext_dig","per_long","tib_post")]
    fig.legend(handles=handles,loc="outside lower center",ncol=5,frameon=False,fontsize=9)
    a,b=record["beforeRigidConflicts"],record["rigidConflicts"]
    fig.suptitle(f"Coordinating ankle routes reduces modeled hardware conflicts from {a} to {b}",
                 fontsize=16,fontweight="bold")
    fig.savefig(folder/"ankle_defense_before_after.png",dpi=200)
    fig.savefig(folder/"ankle_defense_before_after.pdf")
    plt.close(fig)


def triangle_tube(a,b,r):
    a,b=np.asarray(a),np.asarray(b)
    d=b-a
    if np.linalg.norm(d)<1e-8:
        return []
    d=d/np.linalg.norm(d)
    ref=np.eye(3)[np.argmin(abs(d))]
    u=np.cross(d,ref)
    u=u/np.linalg.norm(u)
    v=np.cross(d,u)
    t=np.linspace(0,2*np.pi,15)
    ring=r*(np.cos(t)[:,None]*u+np.sin(t)[:,None]*v)
    return [[a+ring[i],a+ring[i+1],b+ring[i+1],b+ring[i]] for i in range(len(t)-1)]


def render(folder,scene,surfaces):
    fig=plt.figure(figsize=(8.7,8.5),facecolor="white")
    ax=fig.add_subplot(projection="3d",proj_type="ortho")
    def reorder(x):
        return np.asarray(x)[...,[0,2,1]]
    for t in surfaces:
        ax.add_collection3d(Poly3DCollection(reorder(t),facecolor="#b9c1ca",
                                            edgecolor="none",alpha=.46))
    for p in scene["after"]:
        color=COLORS_BY_GROUP[p["group"]]
        quads=[]
        for v in p["vessels"]:
            quads+=triangle_tube(*v["points_mm"],v["radius_mm"])
        if quads:
            ax.add_collection3d(Poly3DCollection(reorder(quads),facecolor=color,
                                                edgecolor="none",alpha=.86))
        w=reorder(p["wire_mm"])
        ax.plot(*w.T,color=color,lw=1.7)
        if "mount_reach_mm" in p:
            ax.plot(*reorder(p["mount_reach_mm"]).T,color="#474f59",lw=1.15,
                    linestyle=":",alpha=.6)
        u,t=np.mgrid[0:2*np.pi:14j,0:np.pi:9j]
        for h in p["hardware_mm"]:
            # The conservative hard-envelope radius is plotted directly.
            x,z,y=h[:3]
            r=h[3]
            ax.plot_surface(x+r*np.cos(u)*np.sin(t),y+r*np.sin(u)*np.sin(t),
                            z+r*np.cos(t),color=color,alpha=.30,linewidth=0,shade=True)
    ax.set(xlim=(-160,180),ylim=(-150,150),zlim=(-550,30),box_aspect=(340,300,580))
    ax.set_axis_off()
    fig.subplots_adjust(left=0,right=1,top=.97,bottom=0)
    for name,azimuth in (("anterolateral",-55),("posterolateral",125)):
        ax.view_init(elev=10,azim=azimuth)
        fig.savefig(folder/f"ankle_{name}_solved_envelopes.png",dpi=170,facecolor="white")
    frames=[]
    fig.set_dpi(74)
    for azimuth in np.linspace(-55,125,19):
        ax.view_init(elev=10,azim=azimuth)
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert("RGB"))
    frames[0].save(folder/"ankle_180deg_orbit.gif",save_all=True,append_images=frames[1:],
                   duration=110,loop=0,optimize=True)
    plt.close(fig)


def torque_plot(folder):
    with (folder/"torque_and_force.csv").open(newline="") as stream:
        data=list(csv.DictReader(stream))
    labels=[r["group"].replace("_"," ") for r in data]
    y=np.arange(len(data))
    fig,axes=plt.subplots(1,2,figsize=(11,5.2),layout="constrained",sharey=True)
    for ax,key,label in zip(axes,("ankle_met_pct","second_met_pct"),
                            ("Ankle","MTP or subtalar")):
        for suffix,color,symbol,text in (("before_",GOLD,"o","Earlier layout"),
                                          ("",INDIGO,"s","Coordinated layout")):
            values=[float(r[suffix+key]) if r[suffix+key] else np.nan for r in data]
            ax.scatter(values,y,color=color,marker=symbol,edgecolor="#30323a",s=52,label=text)
        ax.set(xlim=(-4,104),xlabel="Sampled poses meeting signed human target (%)",title=label)
        ax.grid(axis="x",alpha=.2)
    axes[0].set_yticks(y,labels)
    axes[0].invert_yaxis()
    axes[1].legend(frameon=False,loc="lower right")
    fig.suptitle("Torque attainment varies by group after coordinated routing",fontweight="bold")
    fig.savefig(folder/"ankle_defense_signed_torque.png",dpi=190)
    plt.close(fig)


def focus_plot(folder,scene,surfaces,pairs):
    selected=[p for p in pairs if {p["group_a"],p["group_b"]}=={"tib_ant","ext_hal"}]
    worst=min((float(p["rigid_margin_mm"]) for p in selected),default=float("nan"))
    fig,axes=plt.subplots(1,4,figsize=(12,6),layout="constrained",sharey=True)
    for j,(key,title) in enumerate((("before","Earlier layout"),("after","Coordinated layout"))):
        for ax,dim,label in zip(axes[2*j:2*j+2],(0,2),("Anterior","Right")):
            for tri in surfaces:
                ax.add_collection(PolyCollection(tri[:,:,[dim,1]],facecolor="#a9b1bc",
                                                 edgecolor="none",alpha=.3))
            for part in scene[key]:
                if part["group"] not in ("tib_ant","ext_hal"):
                    continue
                color=COLORS_BY_GROUP[part["group"]]
                for vessel in part["vessels"]:
                    v=np.asarray(vessel["points_mm"])[:,[dim,1]]
                    capsule(ax,*v,vessel["radius_mm"],color,.50)
                for h in part["hardware_mm"]:
                    capsule(ax,[h[dim],h[1]],[h[dim],h[1]],h[3],color,.34)
                w=np.asarray(part["wire_mm"])
                ax.plot(w[:,dim],w[:,1],color=color,lw=1.8)
            ax.set(title=title if dim==0 else "Frontal view",aspect="equal",
                   xlim=(-125,160),ylim=(-505,-65),xlabel=f"{label} (mm)")
            ax.grid(alpha=.15)
    axes[0].set_ylabel("Superior (mm)")
    fig.suptitle(f"Tibialis anterior–EHL minimum modeled margin after routing {worst:+.2f} mm",
                 fontweight="bold")
    axes[-1].legend(handles=[
        plt.Line2D([0],[0],color=COLORS_BY_GROUP[g],lw=3,label=label)
        for g,label in (("tib_ant","Tibialis anterior"),("ext_hal","EHL"))],loc="lower left")
    fig.savefig(folder/"tib_ant_ehl_before_after.png",dpi=200)
    plt.close(fig)


def dissertation_comparison(folder):
    # The historical capability comparison is independent of this ankle run.
    fig,axes=plt.subplots(1,2,figsize=(10,4.7),layout="constrained",sharey=True)
    for ax,label,values in zip(axes,("Signed human torque target met","Active force"),
                               ((37.35,65.81),(93.00,99.31))):
        bars=ax.bar(("2019 routing\nbaseline","Selected\nconfigurations"),values,
                    color=[GOLD,INDIGO],width=.57)
        ax.bar_label(bars,labels=[f"{x:.2f}%" for x in values],padding=3,fontsize=13)
        ax.set(title=label,ylim=(0,112))
        ax.set_axisbelow(True)
        ax.grid(axis="y",alpha=.18)
    axes[0].set_ylabel("Mean angular-domain coverage (%)")
    fig.suptitle("Selected configurations improved 27-site torque and force coverage",
                 fontsize=15,fontweight="bold")
    fig.supxlabel("Dissertation capability study · 620 kPa · Chapter 4, Whole-Leg Actuator and Tackle Capability",
                  fontsize=9)
    fig.savefig(folder/"dissertation_2019_capability_comparison.png",dpi=200)
    fig.savefig(folder/"dissertation_2019_capability_comparison.pdf")
    plt.close(fig)


def main(folder,geometry):
    record=json.loads((folder/"assembly_routes.json").read_text())
    scene=json.loads((folder/"scene_clash_pose.json").read_text())
    surfaces=load_meshes(scene,geometry)
    projections(folder,scene,surfaces,record)
    render(folder,scene,surfaces)
    torque_plot(folder)
    routes=[(f"{r['group']}_{r['channel']}_independent",r,(0,.2,.85))
            for r in record["actuators"]]
    export_routes(record["model"],folder/"ankle_assembled_preview.osim",geometry,routes)
    with (folder/"pair_clearance.csv").open(newline="") as stream:
        pairs=list(csv.DictReader(stream))
    focus_plot(folder,scene,surfaces,pairs)
    dissertation_comparison(folder)
    with (folder/"torque_and_force.csv").open(newline="") as stream:
        torque=list(csv.DictReader(stream))
    active=100-np.mean([float(r["no_force_pct"]) for r in torque])
    old_active=100-np.mean([float(r["before_no_force_pct"]) for r in torque])
    max_overlap=max(float(r["sleeve_overlap_mm"]) for r in pairs)
    home=[r["minimumHomeZ_mm"] for r in record["actuators"] if "minimumHomeZ_mm" in r]
    fails=sorted((p for p in pairs if float(p["rigid_margin_mm"])<-.1),
                 key=lambda p:float(p["rigid_margin_mm"]))
    lines=["# Coordinated ankle assembly", "",
           f"All 17 nominal 20 mm BPAs were modeled with independent tendons and brackets. "
           f"The selected layout was checked at {record['poseCount']} simultaneous ankle, subtalar, and MTP poses.", "",
           f"Rigid conflicts changed from {record['beforeRigidConflicts']} to {record['rigidConflicts']} "
           f"actuator pairs. Local route failures changed from {record['beforeLocalFailures']} to "
           f"{record['localFailures']}. Mean capped signed torque coverage changed from "
           f"{record['beforeCoverage']:.1%} to {record['coverage']:.1%}. Mean active-force coverage "
           f"changed from {old_active:.1f}% to {active:.1f}%. "
           f"The largest modeled vessel overlap is {max_overlap:.2f} mm.",
           "", "| Group | Ankle target met before → after | Second target before → after | No force before → after |",
           "|---|---:|---:|---:|"]
    for r in torque:
        def change(key):
            a,b=r["before_"+key],r[key]
            if not a or not b or a.lower()=="nan" or b.lower()=="nan":
                return "—"
            return f"{float(a):.1f}% → {float(b):.1f}%"
        lines.append(f"| {r['group']} | {change('ankle_met_pct')} | {change('second_met_pct')} | "
                     f"{change('no_force_pct')} |")
    lines+=["","## Remaining modeled rigid conflicts","",
            "| Actuator A | Actuator B | Margin |","|---|---|---:|"]
    for p in fails:
        lines.append(f"| {p['group_a']} {p['channel_a']} | {p['group_b']} {p['channel_b']} | "
                     f"{float(p['rigid_margin_mm']):.2f} mm |")
    lines+=["","## Individual route blockers","",
            "| Actuator | Tibia margin | Tendon lead margin | Blockers |","|---|---:|---:|---|"]
    for r in record["actuators"]:
        notes=", ".join(r["failedConstraints"]) or "pass"
        lines.append(f"| {r['group']} {r['channel']} | {r['boneMargin_mm']:+.2f} mm | "
                     f"{r['leadMargin_mm']:+.2f} mm | {notes} |")
    if home:
        lines.extend(["", f"The smallest route-point z in the home pelvis frame is {min(home):+.2f} mm. "
                      "The side rule applies at home."])
    lines+=["","The before–after comparison uses the preceding 20 mm ankle assembly on the same pose grid. "
            "The dissertation’s 2019 comparison is a separate 27-site, area-weighted capability study. "
            "The visuals show solved vessel envelopes and provisional fitting, stop, and tackle envelopes. "
            "Dotted reach lines indicate candidate mount supports. "
            "The preview is a routing model, with CAD bracket structure and physical reeving still to be verified."]
    (folder/"ANKLE_ASSEMBLY_VERIFICATION.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


if __name__=="__main__":
    if sys.argv[2]=="--historical-only":
        dissertation_comparison(Path(sys.argv[1]))
    else:
        main(Path(sys.argv[1]),Path(sys.argv[2]))
