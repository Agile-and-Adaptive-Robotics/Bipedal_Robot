"""Plot the two-guide bypass and three separate two-BPA vasti lanes."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import numpy as np

from fr_exportPilotOpenSim import export_routes
from fr_plotAnkleClash import capsule
from fr_plotAssembledAnkle import load_meshes
from fr_plotDistalStages import MARGINS
from fr_wholeLegPack import Envelope, clearance, read, table

NAMES=("Vastus medialis","Vastus intermedius","Vastus lateralis")
COLORS=("#1b7779","#574685","#875741")
GRAC="#1d2d6a"


def route(body,points):
    return dict(bodies=body,points_m=np.asarray(points).tolist())


def routes_of(result):
    routes=[]
    for name, v, color in zip(NAMES,result["vasti"],COLORS):
        for i, offset in enumerate(v["channel_m"],1):
            off=np.array([offset[1],0,offset[0]])
            path=np.vstack([v["points_m"][0],np.asarray(v["home_cam_femur_m"])[:-1]])+off
            end=np.asarray(v["terminal_m"])+off
            r=route(["femur_r"]*len(path)+["tibia_r"],np.vstack([path,end]))
            routes.append((name.lower().replace(" ","_")+f"_{i}",r,color))
    routes.append(("gracilis_two_guides",result["grac"],GRAC))
    return routes


def geometry(pose, result, record, index, j):
    frames={b:pose["transforms"][b][index] for b in ("pelvis","femur_r","tibia_r")}
    local=np.asarray(record["points_m"])
    if j<6:
        v=result["vasti"][j//2]
        k=np.argmin(abs(np.asarray(v["knee_deg"])-pose["angles_deg"][index][3]))
        delta=v["channel_m"][j%2]
        offset=np.array([delta[1],0,delta[0]])
        contact=np.asarray(v["contactPath_m"][k])
        local=np.vstack([v["points_m"][0],contact[:-1]])+offset
        local=np.vstack([local,np.asarray(v["terminal_m"])+offset])
        frames_used=["femur_r"]*(len(local)-1)+["tibia_r"]
        sleeve,rad,gain=v["sleeve_m"][k],v["radius_m"][k],v["gain"]
    else:
        frames_used=record["bodies"]
        sleeve=record.get("sleeve",record["sleeve_m"])[index]
        rad=record.get("radius",record["radius_m"])[index]
        gain=record["gain"]
    p=np.array([(np.asarray(frames[b])@np.r_[x,1])[:3] for b,x in zip(frames_used,local)])
    u=p[1]-p[0]; u/=np.linalg.norm(u)
    a,b=p[0]+.025*u,p[0]+(.025+sleeve)*u
    return p,u,a,b,rad,gain


def group_parts(pose,result,routes,tackle_radius):
    parts=[[[] for _ in pose["angles_deg"]] for _ in range(4)]
    for i in range(len(pose["angles_deg"])):
        for j,(_,record,_) in enumerate(routes):
            p,u,a,b,rad,gain=geometry(pose,result,record,i,j)
            family=parts[j//2 if j<6 else 3][i]
            def add(x,y,radius,kind):
                family.append(np.r_[x,y,radius,kind,j+1].tolist())
            add(a,b,rad,0)
            for center in (a-.0125*u,b+.0125*u):
                add(center,center,np.hypot(rad,.0125),1)
            if gain>1:
                add(b,b+.012*u,tackle_radius,1)
            if j==6:
                for guide in p[1:3]:
                    add(guide,guide,.008,1)
            path=np.vstack([b,p[1:]])
            for x,y in zip(path[:-1],path[1:]):
                add(x,y,.001,2)
    return dict(zip((*NAMES,"Gracilis"),parts))


def check_pairs(pose,result,routes,tackle_radius):
    by_group=group_parts(pose,result,routes,tackle_radius)
    items=[Envelope(dict(parts=p)) for p in by_group.values()]
    rows=[]
    names=(*NAMES,"Gracilis")
    for i,a in enumerate(items):
        for j in range(i,len(items)):
            if i==j==3:
                continue
            gap,hit=clearance(a,items[j],same=i==j)
            rows.append(dict(branch=result["branch"],group_a=names[i],group_b=names[j],
                             sleeve_mm=gap[0],hardware_mm=gap[1],cable_mm=gap[2],
                             sleeve_pose=hit[0],hardware_pose=hit[1],cable_pose=hit[2]))
    return rows


def draw_route(out, pose, r, routes):
    fig,axes=plt.subplots(2,2,figsize=(9.5,9.2),layout="constrained")
    q=pose["angles_deg"]
    for row,index in enumerate((0,max(0,r["grac"]["governingPose"]-1))):
        frames={b:pose["transforms"][b][index] for b in ("pelvis","femur_r","tibia_r")}
        mesh=load_meshes(dict(anatomy=pose["anatomy"],transforms=frames),Path("C:/OpenSim 4.6/Geometry"))
        for ax,dim,title in zip(axes[row],(0,2),("Sagittal","Posterior")):
            for tri in mesh:
                ax.add_collection(PolyCollection(tri[:,:,[dim,1]],facecolor="#bac2cb",
                                                 edgecolor="none",alpha=.24))
            for j,(name,record,color) in enumerate(routes):
                p,_,a,b,rad,_=geometry(pose,r,record,index,j)
                capsule(ax,1000*a[[dim,1]],1000*b[[dim,1]],1000*rad,color,.50)
                ax.plot(p[:,dim]*1000,p[:,1]*1000,color=color,lw=1.3)
                if j==6:
                    ax.scatter(p[1:3,dim]*1000,p[1:3,1]*1000,color=color,marker="s",s=24,zorder=5)
            ax.set(title=f"{title} · knee {q[index][3]:.1f}°",
                   xlabel="Anterior +x (mm)" if dim==0 else "Right +z (mm)",aspect="equal")
            ax.autoscale_view()
            ax.grid(alpha=.13)
        axes[row,0].set_ylabel("Superior +y (mm)")
    fig.legend(handles=[Line2D([],[],color=c,lw=3,label=n) for n,c in
                        zip((*NAMES,"Gracilis and its guides"),(*COLORS,GRAC))],
               loc="outside lower center",ncol=2,frameon=False)
    fig.suptitle("Three vasti groups with side-by-side sleeves and a guided gracilis bypass",
                 fontweight="bold")
    fig.savefig(out/"split_vasti_gracilis_routes.png",dpi=185)
    (out/"split_vasti_gracilis_routes.alt.txt").write_text(
        "Sagittal and frontal views show two parallel inflated sleeves for each vastus "
        "and the two femur-fixed gracilis guides. The lower views show the governing "
        "sampled gracilis sleeve-and-guide clearance pose.\n",encoding="utf-8")
    plt.close(fig)


def draw_torque(out,result):
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout="constrained")
    for ax,v,name,color in zip(axes,result["vasti"],NAMES,COLORS):
        ax.plot(v["knee_deg"],v["human_Nm"],"--",color="#6b707b",label="Human target")
        ax.plot(v["knee_deg"],v["robot_Nm"],color=color,label="Two 20 mm BPAs")
        ax.set(title=f"{name}\n{v['coverage']:.0%} of poses met",xlabel="Knee angle (°)",
               xlim=(-120,0))
        ax.grid(alpha=.18)
    axes[0].set_ylabel("Signed knee torque (N·m)")
    axes[-1].legend(frameon=False)
    fig.savefig(out/"separate_vasti_signed_torque.png",dpi=185)
    (out/"separate_vasti_signed_torque.alt.txt").write_text(
        "Three panels show signed knee extension torque for the separate vastus groups "
        "against their respective human muscle torque curves.\n",encoding="utf-8")
    plt.close(fig)
    g=result["grac"]
    axes=g["axes"]
    x,y=np.array(axes["axis1_deg"]),np.array(axes["axis2_deg"])
    human,robot=np.asarray(g["human_Nm"]),np.asarray(g["robot_Nm"])
    fig,plots=plt.subplots(1,2,figsize=(10,4.6),layout="constrained")
    for i,(ax,label) in enumerate(zip(plots,("Hip adduction","Knee"))):
        margin=robot[:,i]*np.sign(human[:,i])-abs(human[:,i])
        lim=max(.1,max(abs(margin)))
        im=ax.pcolormesh(x,y,margin.reshape(len(y),len(x)),cmap=MARGINS,vmin=-lim,vmax=lim)
        ax.set(title=f"{label} · {g['coverage'][i]:.1%} met",xlabel="Hip adduction (°)",
               ylabel="Knee angle (°)",xlim=(x.min(),x.max()),ylim=(y.min(),y.max()))
        fig.colorbar(im,ax=ax,label="Signed torque margin (N·m)")
    fig.savefig(out/"gracilis_bypass_signed_torque.png",dpi=180)
    (out/"gracilis_bypass_signed_torque.alt.txt").write_text(
        "Separate hip adduction and knee torque heatmaps for gracilis. Gold is negative "
        "signed margin and indigo is positive margin.\n",encoding="utf-8")
    plt.close(fig)


def draw_section(out,result):
    fig,ax=plt.subplots(figsize=(6,4.2),layout="constrained")
    for v,name,color in zip(result["vasti"],NAMES,COLORS):
        center=np.asarray(v["points_m"][0])
        radius=1000*max(v["radius_m"])
        for j,(z,x) in enumerate(v["channel_m"]):
            point=1000*(center[[2,0]]+[z,x])
            ax.add_patch(Circle(point,radius,facecolor=color,edgecolor=color,alpha=.35,
                                label=name if j==0 else None))
            ax.scatter(*point,marker="+",s=30,color=color)
    ax.set(xlabel="Right +z in femur frame (mm)",ylabel="Anterior +x in femur frame (mm)",
           title="Vasti sleeve rows at peak modeled inflation",aspect="equal")
    ax.autoscale()
    ax.legend(frameon=False,loc="upper center",bbox_to_anchor=(.5,-.16),ncol=2)
    ax.grid(alpha=.16)
    fig.savefig(out/"vasti_sleeve_cross_section.png",dpi=185)
    plt.close(fig)
    (out/"vasti_sleeve_cross_section.alt.txt").write_text(
        "Cross-section of the six modeled vasti sleeves at each sleeve's peak inflation. "
        "Each group has two side-by-side sleeves and neighboring groups use separate rows.\n",
        encoding="utf-8")


def main(out,pack):
    plt.rcParams.update({"font.family":"Arial","font.size":10})
    pose=read(pack/"dense"/"poses.json")
    summaries,pairs=[],[]
    for branch in ("rigid","compliant"):
        r=read(out/f"split_{branch}.json")
        routes=routes_of(r)
        large=check_pairs(pose,r,routes,.03375)
        medium=check_pairs(pose,r,routes,.02105)
        for original, alternate in zip(large,medium):
            original["hardware_38mm_sheave_proxy_mm"]=alternate["hardware_mm"]
        pairs.extend(large)
        overlays=[(n,v,tuple(int(c[i:i+2],16)/255 for i in (1,3,5))) for n,v,c in routes]
        export_routes(r["model"],out/f"split_vasti_grac_{branch}.osim",
                      Path("C:/OpenSim 4.6/Geometry"),overlays)
        for name,v in zip(NAMES,r["vasti"]):
            summaries.append(dict(branch=branch,group=name,BPAs=v["count"],gain=v["gain"],
                                  hip_met_pct="",knee_met_pct=100*v["coverage"],no_force_pct=100*v["noForce"],
                                  bone_mm=1000*v["bone"],guide_or_cable_bone_mm=1000*v["cableBone"],
                                  sleeve_bay_mm=1000*v["bay"],gap_to_vasti_mm=""))
        g=r["grac"]
        summaries.append(dict(branch=branch,group="Gracilis",BPAs=g["count"],gain=g["gain"],
                              hip_met_pct=100*g["coverage"][0],knee_met_pct=100*g["coverage"][1],
                              no_force_pct=100*g["noForce"],bone_mm=1000*g["bone"],
                              guide_or_cable_bone_mm=1000*g["viaBone"],sleeve_bay_mm=1000*g["bay"],
                              gap_to_vasti_mm=1000*g["vastiGap"]))
        if branch=="compliant":
            draw_route(out,pose,r,routes)
            draw_torque(out,r)
            draw_section(out,r)
    table(out/"split_route_summary.csv",summaries)
    table(out/"split_pair_clearance.csv",pairs)
    notes=["# Gracilis bypass and three vasti groups","",
           "The three vasti groups each use two nominal 20 mm BPAs in shaft-aligned rows. "
           "Each channel has an independent tendon and brackets. The cam uses the recorded "
           "contact samples and a reduced-order tendon closure. The gracilis verification "
           "uses BiPulley for tension and virtual work on the guided cable path for signed "
           "torque. The shank attachment points are the original packed gracilis points.","",
           "The pair clearance table checks the inflated sleeves, cable centerlines with "
           "a 1 mm radius, end fittings, guides, and the original 63.5 mm tackle-envelope "
           "assumption plus 2 mm radial allowance. This tackle size is a packaging proxy. "
           "The 38.1 mm alternate is a footprint sensitivity check with the same allowance. "
           "Hardware placement, tangent geometry, and load ratings still need part-level verification.","",
           "OpenSim overlays were reopened to check native frames. The cam path points in "
           "those overlays are the home-pose preview. The saved MATLAB and JSON data carry "
           "the changing contact geometry.","",
           "| Branch | Group | Hip target met | Knee target met | No force |",
           "|---|---|---:|---:|---:|"]
    for row in summaries:
        hip=f"{row['hip_met_pct']:.1f}%" if row["hip_met_pct"]!="" else "—"
        notes.append(f"| {row['branch']} | {row['group']} | {hip} | "
                     f"{row['knee_met_pct']:.1f}% | {row['no_force_pct']:.1f}% |")
    (out/"KNEE_BYPASS_VERIFICATION.md").write_text("\n".join(notes)+"\n",encoding="utf-8")


if __name__=="__main__":
    main(*map(Path,sys.argv[1:]))
