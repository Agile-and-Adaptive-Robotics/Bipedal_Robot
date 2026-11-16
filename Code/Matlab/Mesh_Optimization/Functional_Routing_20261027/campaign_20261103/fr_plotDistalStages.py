"""Plot distal route centerlines and signed torque margins."""
import json
import csv
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import opensim as osim

from fr_plotPilotAnatomy import mesh_data
from fr_plotGracilisPilotSweeps import human_torque

COLORS_FILE = Path(__file__).resolve().parents[3]/"Colors.m"
COLORS = {name: code for code, name in re.findall(
    r"c\{\d+\}\s*=\s*'(#[0-9A-Fa-f]{6})';\s*%\s*(\w+)",
    COLORS_FILE.read_text(encoding="utf-8"))}
GOLD, INDIGO = COLORS["gold"], COLORS["indigo"]
MARGINS = LinearSegmentedColormap.from_list("gold_indigo",[GOLD,"#ffffff",INDIGO])


def plot(folder, geometry):
    files = [p for p in folder.glob("*.json") if p.stem.endswith(("_rigid","_compliant"))]
    records = [json.loads(p.read_text()) for p in sorted(files)]
    if not records:
        raise ValueError("No branch reports")
    ref = records[0]["model"]
    osim.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry))
    model = osim.Model(ref)
    state = model.initSystem()
    model.realizePosition(state)
    bodies = model.getBodySet()
    pelvis = bodies.get("pelvis")

    def world(body, xyz):
        q = bodies.get(body).findStationLocationInAnotherFrame(state, osim.Vec3(*xyz), pelvis)
        return [q.get(i) for i in range(3)]

    meshes = {}
    xml = ET.parse(ref)
    for body in ("femur_r","tibia_r","talus_r","calcn_r","toes_r"):
        meshes[body] = []
        for mesh in xml.findall(f".//Body[@name='{body}']/attached_geometry/Mesh"):
            xyz, faces = mesh_data(geometry/Path(mesh.findtext("mesh_file")).name)
            xyz *= np.fromstring(mesh.findtext("scale_factors"),sep=" ")
            meshes[body].append((1000*np.array([world(body, p) for p in xyz]),faces))
    legacy_records=[r for r in records if r.get("legacySweep") is not False]
    lines = ["# Distal routing stage", "",
             "Isolated candidates. The recheck uses 31 poses per existing sweep and resampled human targets. "
             "Route views show centerlines at the OpenSim default pose. Gains above 1× are ideal "
             "reeving screens. Local clearance does not establish assembled packaging.", "",
             "| Group | Branch | Length fit | Gain | Signed target poses met | Local constraint failures |",
             "|---|---|---|---:|---:|---|"]
    for r in legacy_records:
        fails = [c["name"] for c in r["constraints"] if c["violated"]]
        lines.append(f"| {r['group']} | {r['branch']} | {'yes' if r['foundLengthFit'] else 'no'} | "
                     f"{r['gain']}× | {r['fractionMet']:.1%} | {', '.join(fails) or 'none'} |")
    ankle_only=[r for r in records if "noForceFraction" in r]
    if ankle_only:
        lines.extend(["","## Ankle-only force coverage","",
                      "| Group | Branch | Sampled RoM with no force |","|---|---|---:|"])
        lines.extend(f"| {r['group']} | {r['branch']} | {r['noForceFraction']:.1%} |" for r in ankle_only)
    optimization=[r for r in records if "surrogateEvaluations" in r]
    if optimization:
        lines.extend(["", "## Surrogate searches", "",
                      "| Group | Branch | Evaluations | Seed objective | Final objective |",
                      "|---|---|---:|---:|---:|"])
        for r in optimization:
            lines.append(f"| {r['group']} | {r['branch']} | {r['surrogateEvaluations']} | "
                         f"{r['initialScore']:.5g} | {r['finalScore']:.5g} |")
    if legacy_records:
        (folder/"STATUS.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    if any("coupled" in r for r in records):
        plot_coupled(folder,records,geometry)
    for branch in ("rigid","compliant"):
        rows = [r for r in records if r["branch"]==branch]
        if rows:
            legacy_rows=[r for r in rows if r.get("legacySweep") is not False]
            if legacy_rows:
                plot_margins(folder,branch,legacy_rows)
            plot_routes(folder,branch,rows,meshes,world)
            for r in rows:
                if "innerRoute" in r:
                    plot_extensor_layer(folder,branch,r,meshes,world)
    print(f"Plotted {len(records)} branch reports")


def plot_coupled(folder, records, geometry):
    ref=Path(records[0]["model"]).with_name("gait2392_simbody.osim")
    model=osim.Model(str(ref))
    state=model.initSystem()
    curve=osim.ActiveForceLengthCurve()
    compare=any("surrogateEvaluations" in r or "baselineCoupled" in r for r in records)
    if "shaft_seed" in folder.name:
        baseline=folder.with_name("ankle_edl_wrap_600")
    elif any("layer" in r for r in records):
        baseline=folder.with_name("ankle_mounts_600")
    else:
        baseline=folder.with_name("ankle_surrogate_600" if "mounts" in folder.name else "ankle_coupled")
    lines=["# Coupled ankle stage", "",
           "Each non-soleus row uses one BiPulley equilibrium per ankle × second-coordinate pose. "
           "The human targets were recomputed on matching poses. Soleus uses an ankle-only sweep. "
           "Geometry checks are local and gains above 1× are ideal topology screens.", "",
           "| Group | Branch | Gain | Ankle met | Second coordinate met | Sampled RoM with no force | Local constraints |"
           + (" Seed ankle met | Seed second coordinate met |" if compare else ""),
           "|---|---|---:|---:|---:|---:|---|" + ("---:|---:|" if compare else "")]
    retention=["# Signed torque retained", "",
               "Coverage is the sum of delivered torque in the target direction, capped at the human "
               "target at each pose, divided by the sum of human target magnitudes. "
               "Slack, stroke failure, and wrong-direction torque contribute zero.", "",
               "| Group | Branch | Ankle coverage | Second-joint coverage | Local constraints |",
               "|---|---|---:|---:|---|"]
    clearance=["# Coupled ankle clearance", "",
               "Positive margins clear the modeled requirement. Distal guide margin includes the "
               "two fitting lengths. Bundle spread is the multiplier on the resting channel spacing. "
               "The table does not establish clearance between groups or physical sheave hardware.", "",
               "| Group | Branch | Bone (mm) | Tibial lead (mm) | Distal guide (mm) | "
               "Bundle (mm) | Bundle spread | Tendon shortfall (mm) | Sampled RoM with no force |",
               "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for branch in ("rigid","compliant"):
        rows=[r for r in records if r["branch"]==branch and "coupled" in r]
        if not rows:
            continue
        panels=[]
        for r in rows:
            e=r["coupled"]
            muscle=model.getMuscles().get(r["group"]+"_r")
            ha,human_path=human_torque(model,state,muscle,curve,"ankle_angle_r",e["poses"],
                                       with_lengths=True)
            hb=human_torque(model,state,muscle,curve,e["second"],e["poses"])
            ra=np.asarray(e["ankle_Nm"],dtype=float)
            rb=np.asarray(e["second_Nm"],dtype=float)
            if "branchResults" in r:
                plot_edl_torque(folder,r,ha,hb)
            invalid=np.asarray(e["slack"])|np.asarray(e["stroke"])|np.asarray(e["invalid"])
            bad=np.asarray(e["noForce"]) if "noForce" in e else invalid
            ma=ra*np.sign(ha)-abs(ha)
            mb=rb*np.sign(hb)-abs(hb)
            a=(ma>=-1e-6)&~bad; b=(mb>=-1e-6)&~bad
            hit=a&b
            x=np.unique(e["ankle_deg"]); y=np.unique(e["second_deg"])
            panels.append((r["group"],e["second"],x,y,bad,
                           (ma,mb),(a.sum(),b.sum())))
            name=r["group"]+"_"+branch
            with (folder/f"{name}_coupled_torque.csv").open("w",newline="") as stream:
                w=csv.writer(stream)
                w.writerow(["ankle_deg",e["second"]+"_deg","human_ankle_Nm","robot_ankle_Nm",
                            "human_second_Nm","robot_second_Nm","bad_transmission","no_force","both_met",
                            "human_path_m","human_Fmax_N"])
                w.writerows(zip(e["ankle_deg"],e["second_deg"],ha,ra,hb,rb,invalid,bad,hit,
                                human_path,np.full(len(ha),muscle.getMaxIsometricForce())))
            line=(f"| {r['group']} | {branch} | {r['gain']}× | {a.mean():.1%} | {b.mean():.1%} | "
                  f"{bad.mean():.1%} | {'pass' if e['localPass'] else 'fail'} |")
            if compare:
                def joint_metrics(h,r0,old_bad):
                    met=(r0*np.sign(h)-abs(h)>=-1e-6)&~old_bad
                    return met.mean(),coverage(h,r0,old_bad)
                if "baselineCoupled" in r:
                    old=r["baselineCoupled"]
                    assert old["ankle_deg"]==e["ankle_deg"] and old["second_deg"]==e["second_deg"]
                    old_bad=np.asarray(old["noForce"])
                    old_a,old_cov_a=joint_metrics(ha,np.asarray(old["ankle_Nm"]),old_bad)
                    old_b,old_cov_b=joint_metrics(hb,np.asarray(old["second_Nm"]),old_bad)
                else:
                    with (baseline/f"{name}_coupled_torque.csv").open(newline="") as stream:
                        old=list(csv.DictReader(stream))
                    old_bad=np.array([row.get("no_force",row["bad_transmission"])=="True" for row in old])
                    def old_joint(key):
                        h=np.array([float(row[f"human_{key}_Nm"]) for row in old])
                        r0=np.array([float(row[f"robot_{key}_Nm"]) for row in old])
                        return joint_metrics(h,r0,old_bad)
                    old_a,old_cov_a=old_joint("ankle")
                    old_b,old_cov_b=old_joint("second")
                line+=f" {old_a:.1%} | {old_b:.1%} |"
                cov_a=f"{old_cov_a:.0%} → {coverage(ha,ra,bad):.0%}"
                cov_b=f"{old_cov_b:.0%} → {coverage(hb,rb,bad):.0%}"
            else:
                cov_a=f"{coverage(ha,ra,bad):.0%}"
                cov_b=f"{coverage(hb,rb,bad):.0%}"
            retention.append(f"| {r['group']} | {branch} | {cov_a} | {cov_b} | "
                             f"{'pass' if e['localPass'] else 'fail'} |")
            def margin_mm(key):
                value=e.get(key)
                return f"{1000*value:+.1f}" if isinstance(value,(int,float)) else "—"
            spread=r.get("bundleSpread")
            spread_label=f"{spread:.2f}" if spread is not None and r["count"]>1 else "—"
            shortage=max(0,1000*e["tendonRequired_m"]-r["tendon_mm"])
            clearance.append(f"| {r['group']} | {branch} | {margin_mm('boneMargin_m')} | "
                             f"{margin_mm('leadMargin_m')} | {margin_mm('distalGuideMargin_m')} | "
                             f"{margin_mm('bundleMargin_m')} | {spread_label} | {shortage:.1f} | {bad.mean():.1%} |")
            lines.append(line)
        coordinates={panel[1] for panel in panels}
        joint=("MTP" if coordinates=={"mtp_angle_r"} else
               "Subtalar" if coordinates=={"subtalar_angle_r"} else "MTP or subtalar")
        for i, (suffix,label) in enumerate((("ankle","Ankle torque"),
                                            ("second_joint",f"{joint} torque"))):
            plot_coupled_panels(folder,branch,panels,i,suffix,label)
    (folder/"COUPLED_STATUS.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    (folder/"TORQUE_RETENTION.md").write_text("\n".join(retention)+"\n",encoding="utf-8")
    (folder/"CLEARANCE.md").write_text("\n".join(clearance)+"\n",encoding="utf-8")
    layers=[r for r in records if "layer" in r]
    if layers:
        report=["# EDL around EHL", "",
                "| Branch | EDL outward margin (mm) | Fitting margin (mm) | Sleeve overlap (mm) | "
                "Bend radius (mm) | Max Xi3 length loss (mm) | Shaft tendon margin (mm) | "
                "Guide turn | Worst stop turn | Pass |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
        def mm(value):
            return f"{1000*value:+.2f}" if isinstance(value,(int,float)) else "—"
        for r in layers:
            L=r["layer"]; e=r["coupled"]
            guide=e.get("guideTurn_deg"); stop=e.get("stopTurn_deg")
            guide=f"{guide:.1f}°" if guide is not None else "—"
            stop=f"{stop:.1f}°" if stop is not None else "—"
            report.append(f"| {r['branch']} | {mm(L['orderingMargin_m'])} | {mm(L['rigidMargin_m'])} | "
                          f"{mm(L['sleeveOverlap_m'])} | {mm(e.get('bpaBendRadius_m'))} | "
                          f"{mm(e.get('bpaBendLossMax_m'))} | {mm(e.get('shaftCableMargin_m'))} | "
                          f"{guide} | {stop} | {'yes' if e['localPass'] else 'no'} |")
        (folder/"LAYER_CLEARANCE.md").write_text("\n".join(report)+"\n",encoding="utf-8")


def coverage(human,robot,bad):
    target=np.abs(human)
    delivered=np.clip(robot*np.sign(human),0,target)
    delivered[bad]=0
    return delivered.sum()/max(target.sum(),np.finfo(float).eps)


def plot_edl_torque(folder,record,ankle_target,mtp_target):
    e=record["coupled"]
    ankle=np.asarray(e["ankle_deg"])
    mtp=np.asarray(e["second_deg"])
    cases=(("ankle",ankle,mtp,"Ankle (°)","MTP",ankle_target,"ankle_Nm"),
           ("mtp",mtp,ankle,"MTP (°)","Ankle",mtp_target,"second_Nm"))
    for name,x,other,xlabel,other_label,target,key in cases:
        fig,axes=plt.subplots(1,3,figsize=(10,3.1),sharey=True,layout="constrained")
        values=np.unique(other)
        for ax,held in zip(axes,values[[0,len(values)//2,-1]]):
            ids=np.flatnonzero(other==held)
            ids=ids[np.argsort(x[ids])]
            ax.axhline(0,color="#b3bbc0",lw=.8)
            ax.plot(x[ids],target[ids],"--",color="#4c5158",lw=1.7,label="Human EDL")
            for i,(branch,style) in enumerate(zip(record["branchResults"],(":", "-.")),1):
                ax.plot(x[ids],np.asarray(branch[key])[ids],style,color=COLORS["magenta"],
                        lw=1.15,label=f"BPA {i}")
            ax.plot(x[ids],np.asarray(e[key])[ids],"-",color=INDIGO,lw=1.7,label="Two BPAs")
            ax.set(title=f"{other_label} {held:+.0f}°",xlabel=xlabel)
            ax.grid(alpha=.2)
        axes[0].set_ylabel("Signed torque (N·m)")
        axes[-1].legend(fontsize=8)
        fig.suptitle(f"Independent EDL branches · signed {name.upper() if name=='mtp' else name} torque")
        fig.savefig(folder/f"edl_{name}_signed_torque.png",dpi=170)
        plt.close(fig)


def plot_coupled_panels(folder,branch,panels,i,suffix,label):
    ncol=min(3,len(panels)); nrow=(len(panels)+ncol-1)//ncol
    fig,axes=plt.subplots(nrow,ncol,figsize=(max(5,3.35*ncol),max(4,2.7*nrow)),
                         layout="constrained",squeeze=False)
    for ax,(group,coord,x,y,bad,margins,hits) in zip(axes.flat,panels):
        margin=margins[i].copy()
        margin[bad]=np.nan
        limit=max(1.0,np.max(abs(margin[~bad]),initial=0))
        mesh=ax.pcolormesh(x,y,margin.reshape(y.size,x.size),shading="nearest",
                           cmap=MARGINS,vmin=-limit,vmax=limit)
        ax.set_facecolor("#c2c5c9")
        ax.set(title=f"{group} · {hits[i]/len(bad):.0%} of sampled range",
               xlabel="Ankle angle (°)",ylabel="MTP (°)" if coord.startswith("mtp") else "Subtalar (°)")
        fig.colorbar(mesh,ax=ax,label="Signed margin (N·m)")
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    fig.suptitle(f"{label} · {branch}\nGold: falls short   Indigo: exceeds target   Gray: no force")
    fig.savefig(folder/f"{branch}_{suffix}_heatmaps.png",dpi=160)
    plt.close(fig)


def plot_margins(folder, branch, rows):
    coords = sorted({r["coordinate"] for r in rows})
    fig, axs = plt.subplots(len(coords),1,figsize=(8,2+1.65*len(coords)),
                            layout="constrained",squeeze=False)
    for ax, coord in zip(axs[:,0],coords):
        block = [r for r in rows if r["coordinate"]==coord]
        angles = np.asarray(block[0]["angles_deg"])
        margins = []
        for r in block:
            human,robot = (np.asarray(r[key],dtype=float) for key in ("human_Nm","robot_Nm"))
            x = np.asarray(r["angles_deg"])
            value = robot*np.sign(human)-abs(human)
            order = np.argsort(x)
            margins.append(np.interp(angles,x[order],value[order]))
        limit = max(1.0,np.nanmax(np.abs(margins)))
        mesh = ax.pcolormesh(angles,np.arange(len(block)),margins,shading="nearest",
                             cmap=MARGINS,vmin=-limit,vmax=limit)
        ax.set_yticks(np.arange(len(block)),[r["group"] for r in block])
        ax.set_xlabel(coord.removesuffix("_r").replace("_"," ")+" (°)")
        fig.colorbar(mesh,ax=ax,label="Signed torque margin (N·m)")
    fig.suptitle(f"Distal stage · {branch} · primary coordinate sweeps\n"
                 "Gold: falls short   Indigo: exceeds target")
    fig.savefig(folder/f"{branch}_torque_margins.png",dpi=160)
    plt.close(fig)


def plot_routes(folder, branch, rows, meshes, world):
    ncol=5
    nrow=(len(rows)+ncol-1)//ncol
    fig, axes=plt.subplots(nrow,ncol,figsize=(13,4.6*nrow),layout="constrained",squeeze=False)
    for ax,r in zip(axes.flat,rows):
        points=1000*np.array([world(b,p) for b,p in zip(r["bodies"],r["points_m"])])
        for body,parts in meshes.items():
            if body=="femur_r" and "femur_r" not in r["bodies"]:
                continue
            for xyz,faces in parts:
                ax.add_collection(PolyCollection(xyz[faces][:,:,[0,1]],
                                                 facecolor="#a8b4c1",edgecolor="none",alpha=.25))
        ax.plot(points[:,0],points[:,1],"-o",color=INDIGO,lw=1.7,ms=3)
        if "innerRoute" in r:
            inner=r["innerRoute"]
            ehl=1000*np.array([world(b,p) for b,p in zip(inner["bodies"],inner["points_m"])])
            ax.plot(ehl[:,0],ehl[:,1],"--",color="#566472",lw=1.4,label="EHL reference")
            ax.legend(fontsize=7)
        if "mountRows" in r:
            for row,marker in zip(r["mountRows"],("^","s","D")):
                pt=points[row-1]
                ax.plot(pt[0],pt[1],marker,ms=5,mfc="white",mec="#263044",zorder=4)
        ax.set_aspect("equal",adjustable="box")
        ax.set_xlim(points[:,0].min()-60,points[:,0].max()+120)
        ax.set_ylim(points[:,1].min()-90,points[:,1].max()+50)
        if "coupled" in r:
            label=r["group"]
            clears=r["coupled"]["localPass"]
        else:
            label=f"{r['group']} · {r['fractionMet']:.0%} ankle"
            clears=r["localConstraintsPass"]
        ax.set_title(f"{label}\n"
                     f"{'local pass' if clears else 'local fail'}",fontsize=9)
        ax.set_xlabel("Anterior x (mm)",fontsize=8)
        ax.grid(alpha=.16)
    for ax in list(axes.flat)[len(rows):]:
        ax.set_visible(False)
    key="\n△ Proximal mount   □ Distal guide   ◇ Shank stop" if any("mountRows" in r for r in rows) else ""
    fig.suptitle(f"Distal stage · {branch} · default pose, pelvis axes"+key)
    fig.savefig(folder/f"{branch}_routes.png",dpi=160)
    plt.close(fig)


def plot_extensor_layer(folder,branch,r,meshes,world):
    def points(record):
        return 1000*np.array([world(b,p) for b,p in zip(record["bodies"],record["points_m"])])
    edl=points(r)
    ehl=points(r["innerRoute"])
    old=None
    previous=r.get("baselineFolder")
    if previous is None and "shaft_seed" in folder.name:
        previous="ankle_edl_wrap_600"
    if previous:
        source=folder.with_name(previous)/f"ext_dig_{branch}.json"
        if source.is_file():
            old=points(json.loads(source.read_text()))
    fig,axes=plt.subplots(1,2,figsize=(8,6),layout="constrained",sharey=True)
    for ax,dim,name in zip(axes,(0,2),("Anterior +x (mm)","Right +z (mm)")):
        for body in ("tibia_r","talus_r","calcn_r","toes_r"):
            for xyz,faces in meshes[body]:
                ax.add_collection(PolyCollection(xyz[faces][:,:,[dim,1]],facecolor="#a8b4c1",
                                                 edgecolor="none",alpha=.25))
        ax.plot(ehl[:,dim],ehl[:,1],"--",color="#566472",lw=1.8,label="Deeper EHL reference")
        if old is not None:
            ax.plot(old[:,dim],old[:,1],":",color="#86660b",lw=1.6,label="Previous EDL")
        if "routeParts" in r:
            for i,(route,style,color) in enumerate(zip(
                    r["routeParts"],("-o","-s"),(INDIGO,COLORS["magenta"])),start=1):
                branch_points=points(route)
                ax.plot(branch_points[:,dim],branch_points[:,1],style,color=color,
                        lw=1.6,ms=2,label=f"EDL {i} · separate tendon")
        else:
            ax.plot(edl[:,dim],edl[:,1],"-o",color=INDIGO,lw=1.7,ms=2,label="EDL guided wrap")
        if "guideShift_mm" in r:
            guide=r["mountRows"][1]-1
            ax.plot(edl[guide,dim],edl[guide,1],"s",mfc="white",mec=INDIGO,
                    label="Moved shank guide")
        ax.set_aspect("equal",adjustable="box")
        ax.set_xlabel(name)
        ax.autoscale()
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Superior +y (mm)")
    axes[0].set_title("Sagittal"); axes[1].set_title("Frontal")
    axes[1].legend(fontsize=8)
    L=r["layer"]
    order=L["orderingMargin_m"]; fit=L["rigidMargin_m"]
    margin=(f"Ordering {1000*order:+.1f} mm · Fittings {1000*fit:+.1f} mm"
            if isinstance(order,(int,float)) and isinstance(fit,(int,float)) else "Layer clearance unresolved")
    title=("Two independent EDL routes around EHL" if "routeParts" in r else
           "Extensor digitorum around extensor hallucis longus")
    fig.suptitle(f"{title} · {branch}\n"+margin)
    fig.savefig(folder/f"{branch}_edl_ehl_overlay.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    plot(Path(sys.argv[1]),Path(sys.argv[2]))
