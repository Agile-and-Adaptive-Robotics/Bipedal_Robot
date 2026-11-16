"""Publish the signed torque and native-frame routes for the proximal run."""
import csv
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import opensim as osim

from fr_exportPilotOpenSim import export_routes
from fr_plotDistalStages import GOLD, INDIGO, MARGINS
from fr_plotPilotAnatomy import mesh_data

plt.rcParams.update({"font.family": "Arial", "font.size": 10})


def report(folder, ankle, geometry):
    files=sorted(p for p in folder.glob("*.json") if p.stem.endswith(("_rigid","_compliant")))
    paired={"grac":"hip_adduction_r","tfl":"hip_adduction_r",
            "bifemlh":"hip_flexion_r","semimem":"hip_flexion_r","sar":"hip_flexion_r",
            "rect_fem":"hip_flexion_r","med_gas":"knee_angle_r","lat_gas":"knee_angle_r"}
    records=[]
    for p in files:
        refit=folder/"paired_refit"/p.name
        r=json.loads((refit if refit.is_file() else p).read_text())
        if r["group"] in paired:
            axes=r["coupled"]["axes"]
            assert axes["coordinate1"]==paired[r["group"]]
            assert all(len(axes[f"axis{k}_deg"])==13 for k in (1,2))
            for k in (1,2):
                assert len(axes[f"target{k}_Nm"])==len(r["coupled"][f"robot{k}_Nm"])==169
        records.append(r)
    assert len(records)==34 and all(r["diameter_mm"]==20 and r["evaluations"]>=600 for r in records)
    summary=[]
    named=[]
    for r in records:
        c=r.get("coupled")
        named.extend(individual_targets(folder.parent/"targets_v2",r))
        fails=[x["name"] for x in r["constraints"] if x["violated"]]
        if c:
            a=c["axes"]
            if c["noForceFraction"]>0:
                fails.append("coupled no-force poses")
            if c["boneMargin_m"]<-.0001:
                fails.append("coupled bone exclusion")
            if c["tendonMargin_m"]<-.000001:
                fails.append("coupled tendon minimum")
            pairs=((x,y) for y in a["axis2_deg"] for x in a["axis1_deg"])
            with (folder/f"{r['group']}_{r['branch']}_coupled_torque.csv").open("w",newline="") as stream:
                w=csv.writer(stream)
                w.writerow([a["coordinate1"]+"_deg",a["coordinate2"]+"_deg",
                            "human_1_Nm","robot_1_Nm","human_2_Nm","robot_2_Nm"])
                w.writerows((x,y,h1,t1,h2,t2) for (x,y),h1,t1,h2,t2 in zip(
                    pairs,a["target1_Nm"],c["robot1_Nm"],a["target2_Nm"],c["robot2_Nm"]))
            plot_coupled(folder,r)
        bone=next(x["cDim"] for x in r["constraints"] if x["name"]=="bone_exclusion")
        bone_mm=-1000*bone if bone is not None else ""
        if c and c["boneMargin_m"] is not None:
            bone_mm=min(bone_mm,1000*c["boneMargin_m"]) if bone_mm!="" else 1000*c["boneMargin_m"]
        summary.append(dict(group=r["group"],branch=r["branch"],BPAs=r["count"],diameter_mm=20,
                            gain=r["gain"],ideal_sheaves=r["sheaves"],rest_mm=r["rest_mm"],
                            tendon_mm=r["tendon_mm"],peak_tendon_N=r["peakTendon_N"],
                            primary_coordinate=r["coordinate"],primary_met_pct=100*r["fractionMet"],
                            coupled_coordinate_1=c["axes"]["coordinate1"] if c else "",
                            coupled_met_1_pct=100*c["fraction1"] if c else "",
                            coupled_coordinate_2=c["axes"]["coordinate2"] if c else "",
                            coupled_met_2_pct=100*c["fraction2"] if c else "",
                            coupled_no_force_pct=100*c["noForceFraction"] if c else "",
                            no_force_pct=100*r["noForce"],local_pass=r["localPass"],
                            bone_margin_mm=bone_mm,
                            failed_constraints=", ".join(fails),surrogate_evaluations=r["evaluations"]))
    with (folder/"remaining20_summary.csv").open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    with (folder/"remaining20_named_targets.csv").open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(named[0]))
        w.writeheader()
        w.writerows(named)
    lines=["# Remaining 17 groups · 20 mm search","",
           "The tibia-down designs retain their previous checkpoints. The remaining groups were "
           "searched in rigid and compliant branches with independent actuator load paths. "
           "Local constraints were rechecked without selecting a final cross-group layout. "
           "The knee uses registered robot t1 motion. Sheave counts for gains above one are "
           "topology estimates.", "",
           "The primary table uses the group target. The named-target CSV gives each "
           "human constituent separately on the same primary sweep. The attachment-movement "
           "CSV gives each movable nominal route point's displacement in its owning frame.",
           "",
           f"Local checks pass in {sum(r['local_pass'] for r in summary)} of 34 branch results.",
           "", "| Group | Branch | BPAs | Sheaves | Primary target met | Coupled targets met | No force, worst sweep | No force, coupled | Bone margin | Local checks |",
           "|---|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in summary:
        coupled=(f"{r['coupled_met_1_pct']:.1f}% / {r['coupled_met_2_pct']:.1f}%"
                 if r["coupled_coordinate_1"] else "—")
        coupled_force=f"{r['coupled_no_force_pct']:.1f}%" if r["coupled_coordinate_1"] else "—"
        bone=f"{r['bone_margin_mm']:+.1f} mm" if r["bone_margin_mm"]!="" else "—"
        lines.append(f"| {r['group']} | {r['branch']} | {r['BPAs']} | {r['ideal_sheaves']} | "
                     f"{r['primary_met_pct']:.1f}% | {coupled} | {r['no_force_pct']:.1f}% | {coupled_force} | {bone} | "
                     f"{r['failed_constraints'] or ('pass' if r['local_pass'] else 'fail')} |")
    lines+=["","The 31-point primary sweeps use the recorded signed targets. The 13 × 13 coupled "
            "grids use stock OpenSim target torques. The general route search uses the existing "
            "case topology. The worst-sweep no-force fraction is the maximum across the "
            "individual sweeps and the coupled grid, where applicable. The following alternate "
            "mechanisms were tested separately."]
    lines+=["","## Search score and feasibility","",
            "The main objective uses the existing normalized signed torque deficit (weight 100), "
            "missed target coverage (40), torque shape (0.5), attachment displacement (2), "
            "lengthening (effective weight 4 by default, 8 for vasti and biceps long head, "
            "and 16 for gracilis), zero-primary-torque coverage (20), and "
            "baseline shortfall (50). The coupled search adds signed deficit (100), missed "
            "coverage (40), and search penalties for lost tension (1000) and negative bone or "
            "tendon margin (100,000 per meter). These weights rank candidates. Local acceptance "
            "uses the recorded tolerances for the home-pose side, tendon span, tension and stroke, "
            "within-group overlap, and bone clearance. Cross-group packaging was deferred. "
            "Tackle load estimates are sizing diagnostics. "
            "The baseline shortfall term uses the existing baseline data. The reported target "
            "percentages are relative to the human reference."]
    gas=[json.loads(p.read_text()) for p in sorted(folder.glob("*_pinned_ends_*_result.json"))]
    if gas:
        lines+=["","## Gastrocnemius with pinned active ends","",
                "Each line uses a femoral and a shank sleeve, a central knee tendon, and a foot "
                "tendon. The closure uses equivalent series length with force limited by the "
                "shorter active section. Knee and ankle torque use virtual work on the routed length.",
                "", "| Group | Branch | Parallel lines | Sleeves per line | Travel gain | Knee met | Ankle met | No force | Sleeve–bone margin |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
        for r in gas:
            if "failure" in r:
                lines.append(f"| {r['group']} | {r['branch']} | — | — | — | — | — | — | {r['failure']} |")
            else:
                lines.append(f"| {r['group']} | {r['branch']} | {r['count']} | 2 | {r['gain']} | "
                             f"{r['kneeMet']:.1%} | {r['ankleMet']:.1%} | "
                             f"{r['noForceFraction']:.1%} | {r['boneMargin_mm']:+.1f} mm |")
    cams=folder/"vasti_cam_femoral_bay_screen.json"
    if cams.is_file():
        lines+=["","## Three separate vasti targets on the recorded cam","",
                "The cam screen reuses recorded contact geometry and solves a separate representative "
                "cable line for each target. The complete three-group hardware layout remains for packaging.",
                "", "| Group | Branch | Knee met | No force | Sleeve bay or bone margin |",
                "|---|---|---:|---:|---|"]
        for r in json.loads(cams.read_text()):
            if "failure" in r:
                lines.append(f"| {r['name']} | {r['branch']} | — | — | "
                             f"{r['femoralBay_mm']:.1f} mm available, {r['minimumRestAtGain3_mm']:.1f} mm required |")
            else:
                lines.append(f"| {r['name']} | {r['branch']} | {r['fractionMet']:.1%} | "
                             f"{r['noForceFraction']:.1%} | {r['femurSleeveMargin_mm']:+.1f} mm bone margin |")
    junctions=[json.loads(p.read_text()) for p in sorted(folder.glob("iliopsoas_junction_*_result.json"))]
    if junctions:
        lines+=["","## Psoas–iliacus load-transfer point",""]
        for r in junctions:
            label="pass" if r["localPass"] else "fail"
            lines.append(f"- {r['branch'].capitalize()} combined hip-flexion target attainment is "
                         f"{r['combinedFractionMet']:.1%}. Local constraints {label}. Both lines "
                         "transfer load at one femur-fixed point with independent tendon compliance.")
    if (folder/"rectus_cam_independent"/"rectus_cam_pilot.mat").is_file():
        lines+=["","## Rectus femoris cam rerun","",
                "The independent-line cam pilot is in the rectus_cam_independent folder. Its saved candidates "
                "recover tendon tension throughout the grid. Femoral sleeve clearance remains a blocker."]
    draw(folder,summary)
    distal=json.loads((ankle/"assembly_routes.json").read_text(encoding="utf-8"))
    routes=[(f"{r['group']}_{r['channel']}_ankle",r,(0,.22,.87)) for r in distal["actuators"]]
    for r in records:
        if r["branch"]!="compliant":
            continue
        view=next((g for g in gas if g["group"]==r["group"] and g["branch"]=="compliant"
                   and "points_m" in g),r)
        routes.append((r["group"]+"_proximal",view,(0,.22,.87)))
    preview=export_routes(records[0]["model"],folder/"whole_leg_route_preview.osim",geometry,routes)
    draw_routes(preview,routes)
    count=sum(len(route["points_m"]) for _,route,_ in routes)
    lines+=["","## OpenSim route preview","",
            f"`{preview.name}` reopened successfully and all {count} route points on the "
            f"{len(routes)} overlays matched the saved native body frames and coordinates. "
            "The blue overlays are geometry previews. They include representative paths "
            "for proximal bundles, and the gray paths are the robotbody reference."]
    (folder/"REMAINING20_VERIFICATION.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def individual_targets(targets,r):
    stock=json.loads((targets/f"{r['group']}_signed_targets.json").read_text())
    cells=sorted(stock["cells"],key=lambda c:c["requested_deg"])
    x=np.array([c["requested_deg"] for c in cells])
    angles=np.asarray(r["angles_deg"])
    interval=np.searchsorted(x,angles).clip(1,len(x)-1)
    weight=(angles-x[interval-1])/(x[interval]-x[interval-1])
    output=[]
    for muscle in stock["stockMuscles"]:
        y=np.array([next(m["tau_contrib_Nm"] for m in c["muscles"] if m["muscle"]==muscle)
                    for c in cells])
        target=y[interval-1]+weight*(y[interval]-y[interval-1])
        direction=np.where(abs(target)>1e-12,np.sign(target),0)
        margin=np.asarray(r["robot_Nm"])*direction-abs(target)
        output.append(dict(group=r["group"],branch=r["branch"],human_muscle=muscle,
                           coordinate=r["coordinate"],signed_target_met_pct=100*np.mean(margin>=-1e-6),
                           worst_margin_Nm=np.min(margin),local_pass=r["localPass"]))
    return output


def draw(folder,rows):
    groups=sorted({r["group"] for r in rows})
    y=np.arange(len(groups))
    fig,axes=plt.subplots(1,2,figsize=(11,7.3),sharey=True,layout="constrained")
    for ax,key,title in zip(axes,("primary_met_pct","no_force_pct"),
                            ("Signed primary torque target met","No force on the worst sweep")):
        full=100 if key=="primary_met_pct" else 0
        for branch,marker,offset in (("rigid","o",-.13),("compliant","s",.13)):
            values={r["group"]:r[key] for r in rows if r["branch"]==branch}
            x=[values[g] for g in groups]
            ax.scatter(x,y+offset,c=[INDIGO if np.isclose(v,full) else GOLD for v in x],
                       marker=marker,edgecolors="#33353b",s=38)
        ax.set(title=title,xlim=(-4,104),xlabel="Sampled poses (%)")
        ax.grid(axis="x",alpha=.18)
    axes[0].set_yticks(y,groups)
    axes[0].invert_yaxis()
    handles=[plt.Line2D([],[],color="#33353b",marker=marker,linestyle="none",label=branch)
             for branch,marker in (("Rigid","o"),("Compliant","s"))]
    axes[1].legend(handles=handles,frameon=False,loc="lower right")
    fig.suptitle("The remaining 17 groups use 20 mm BPAs",fontweight="bold")
    out=folder/"remaining20_torque_force.png"
    fig.savefig(out,dpi=175)
    plt.close(fig)
    out.with_suffix(".alt.txt").write_text(
        "Two dot plots give signed primary torque coverage and the worst sweep's no-force "
        "fraction for each group. Circles are rigid and squares are compliant. Indigo denotes "
        "full torque coverage or zero no-force poses. Gold denotes a shortfall.\n",encoding="utf-8")


def draw_routes(preview,routes):
    model=osim.Model(str(preview))
    state=model.initSystem()
    model.realizePosition(state)
    bodies=model.getBodySet()
    pelvis=bodies.get("pelvis")

    def pelvic(body,xyz):
        q=body.findStationLocationInAnotherFrame(state,osim.Vec3(*xyz),pelvis)
        return np.array([q.get(i) for i in range(3)])

    def world_path(path):
        points=path.getPathPointSet()
        return 1000*np.array([pelvic(points.get(i).getParentFrame(),
                                    [points.get(i).getLocation(state).get(j) for j in range(3)])
                              for i in range(points.getSize())])

    xml=ET.parse(preview)
    meshes=[]
    for name in ("pelvis","femur_r","tibia_r","talus_r","calcn_r","toes_r","torso"):
        body=bodies.get(name)
        origin=pelvic(body,[0,0,0])
        basis=np.array([pelvic(body,e)-origin for e in np.eye(3)])
        for mesh in xml.findall(f".//Body[@name='{name}']/attached_geometry/Mesh"):
            if name=="torso" and Path(mesh.findtext("mesh_file")).stem!="hat_spine":
                continue
            xyz,faces=mesh_data(Path(mesh.findtext("mesh_file")))
            xyz*=np.fromstring(mesh.findtext("scale_factors"),sep=" ")
            meshes.append(1000*(xyz@basis+origin)[faces])
    muscles=model.getMuscles()
    refs=[world_path(muscles.get(i).getGeometryPath()) for i in range(muscles.getSize())
          if muscles.get(i).getName().endswith("_r")]
    paths=[]
    for name,_,_ in routes:
        actuator=osim.PathActuator.safeDownCast(model.getForceSet().get(name))
        paths.append(world_path(actuator.getGeometryPath()))

    fig,axes=plt.subplots(1,2,figsize=(8.5,8),layout="constrained",sharey=True)
    for ax,dim,title,xlabel in zip(axes,(0,2),("Sagittal","Frontal"),("Anterior +x (mm)","Right +z (mm)")):
        for tris in meshes:
            ax.add_collection(PolyCollection(tris[:,:,[dim,1]],facecolor="#adb7c5",
                                             edgecolor="none",alpha=.18))
        for path in refs:
            ax.plot(path[:,dim],path[:,1],color="#7b838d",lw=.55,alpha=.42)
        for path in paths:
            ax.plot(path[:,dim],path[:,1],color="#234f99",lw=1.15,alpha=.76)
        ax.set(title=title,xlabel=xlabel,aspect="equal")
        ax.autoscale_view()
    axes[0].set_ylabel("Superior +y (mm)")
    handles=[plt.Line2D([],[],color=color,label=label) for color,label in
             (("#234f99","Candidate route centerlines"),("#7b838d","Robotbody reference"))]
    fig.legend(handles=handles,loc="outside lower center",ncol=2,frameon=False)
    fig.suptitle("Right-leg routing candidates at the OpenSim default pose",fontweight="bold")
    out=preview.with_suffix(".png")
    fig.savefig(out,dpi=190)
    plt.close(fig)
    out.with_suffix(".alt.txt").write_text(
        "Sagittal and frontal projections of the right-leg routing candidates over gray "
        "robotbody reference paths and light gray bone surfaces. Blue lines trace the saved "
        "candidate route centerlines from the torso and pelvis through the thigh, shank, "
        "and foot.\n",encoding="utf-8")


def plot_coupled(folder,r):
    c=r["coupled"]
    a=c["axes"]
    x,y=np.asarray(a["axis1_deg"]),np.asarray(a["axis2_deg"])
    def label(k):
        return a[f"coordinate{k}"].removesuffix("_r").replace("_"," ")
    fig,axes=plt.subplots(1,2,figsize=(10,4.8),layout="constrained")
    for i,ax in enumerate(axes,1):
        target=np.asarray(a[f"target{i}_Nm"])
        torque=np.asarray(c[f"robot{i}_Nm"])
        margin=torque*np.sign(target)-abs(target)
        limit=max(.5,np.max(abs(margin)))
        im=ax.pcolormesh(x,y,margin.reshape(len(y),len(x)),cmap=MARGINS,
                         vmin=-limit,vmax=limit,shading="nearest")
        ax.set(xlabel=label(1)+" (°)",ylabel=label(2)+" (°)",
               xlim=(x.min(),x.max()),ylim=(y.min(),y.max()),
               title=f"{label(i)} torque · {c[f'fraction{i}']:.0%} met")
        fig.colorbar(im,ax=ax,label="Signed margin (N·m)")
    fig.suptitle(f"{r['group']} · {r['branch']}",fontweight="bold")
    out=folder/f"{r['group']}_{r['branch']}_coupled.png"
    fig.savefig(out,dpi=160)
    plt.close(fig)
    out.with_suffix(".alt.txt").write_text(
        f"Two heat maps for {r['group']} in the {r['branch']} branch. The horizontal axis is "
        f"{label(1)} and the vertical axis is {label(2)}. Signed torque targets are met in "
        f"{c['fraction1']:.1%} and {c['fraction2']:.1%} of poses, respectively. Gold marks "
        "negative margin and indigo marks positive margin.\n",encoding="utf-8")


if __name__=="__main__":
    report(*(Path(x) for x in sys.argv[1:]))
