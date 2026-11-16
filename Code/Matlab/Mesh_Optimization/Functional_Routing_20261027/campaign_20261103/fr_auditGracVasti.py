"""Audit the era deltas and replay the governing gracilis/vasti contact."""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fr_plotAnkleClash import capsule
from fr_plotAssembledAnkle import load_meshes
from fr_wholeLegPack import distances, read, table, write


def run(source, out):
    with (source/"era_comparison_by_group.csv").open() as f:
        rows={(r["group"],r["branch"],r["era"]):r for r in csv.DictReader(f)}
    delta=[]
    for group,branch,era in rows:
        if era!="2019":
            continue
        old,new=rows[group,branch,era],rows[group,branch,"current"]
        row=dict(group=group,branch=branch)
        for k in ("torque_margin_met_pct","active_rom_pct","signed_capacity_pct"):
            a,b=float(old[k]),float(new[k])
            row.update({f"old_{k}":a,f"current_{k}":b,f"delta_{k}":b-a})
        delta.append(row)
    table(out/"era_group_deltas.csv",sorted(delta,key=lambda x:(x["branch"],x["delta_torque_margin_met_pct"])))
    dense=source/"dense"
    pose=read(dense/"poses.json")
    opts=[read(dense/f"{g}_compliant.json")[-1] for g in ("grac","vas_int")]
    with (source/"compliant_packed_pairs.csv").open() as f:
        hits=next(r for r in csv.DictReader(f) if (r["group_a"],r["group_b"])==("grac","vas_int"))
    witnesses={}
    for tag,i in (("hardware",int(hits["hardware_pose"])),("sleeve",int(hits["sleeve_pose"]))):
        a,b=[np.asarray(r["parts"][i]) for r in opts]
        x,y=np.repeat(a,len(b),axis=0),np.tile(b,(len(a),1))
        keep=(x[:,7]==1)|(y[:,7]==1) if tag=="hardware" else (x[:,7]==0)&(y[:,7]==0)
        x,y=x[keep],y[keep]
        gap=distances(x,y)-x[:,6]-y[:,6]
        j=np.argmin(gap)
        np.testing.assert_allclose(1000*gap[j],float(hits[f"{tag}_gap_mm"]),atol=1e-7)
        witnesses[tag]=dict(pose_index=i,pose_deg=dict(zip(pose["keys"],pose["angles_deg"][i])),
                            gap_mm=1000*gap[j],grac_part=x[j].tolist(),vasti_part=y[j].tolist())
    write(out/"grac_vasti_collision.json",witnesses)
    draw(out,pose,opts,witnesses)


def draw(out,pose,opts,witnesses):
    from matplotlib.collections import PolyCollection
    fig,axes=plt.subplots(2,2,figsize=(9,8),layout="constrained")
    for row,tag in enumerate(("hardware","sleeve")):
        w=witnesses[tag]
        i=w["pose_index"]
        frames={k:pose["transforms"][k][i] for k in ("pelvis","femur_r","tibia_r")}
        mesh=load_meshes(dict(anatomy=pose["anatomy"],transforms=frames),Path("C:/OpenSim 4.6/Geometry"))
        local=np.array([w["grac_part"][:6],w["vasti_part"][:6]]).reshape(-1,3)*1000
        center=local.mean(axis=0)
        for ax,dim,label in zip(axes[row],(0,2),("Anterior +x (mm)","Right +z (mm)")):
            for tri in mesh:
                ax.add_collection(PolyCollection(tri[:,:,[dim,1]],facecolor="#bac2cb",edgecolor="none",alpha=.24))
            for option,key,col in zip(opts,("grac_part","vasti_part"),("#2b3082","#b68216")):
                chosen=w[key]
                for p in option["parts"][i]:
                    a,b=np.array(p[:3])*1000,np.array(p[3:6])*1000
                    if p==chosen or p[7]==0:
                        capsule(ax,a[[dim,1]],b[[dim,1]],1000*p[6],col,.64 if p==chosen else .14)
                    if p==chosen:
                        ax.plot([a[dim],b[dim]],[a[1],b[1]],"o-",color=col,lw=2.0,ms=4)
            ax.set(xlim=center[dim]+np.array([-140,140]),ylim=center[1]+np.array([-130,130]),
                   aspect="equal",xlabel=label)
            ax.grid(alpha=.18)
        d=w["pose_deg"]
        axes[row,0].set_title(f"{tag.title()} {w['gap_mm']:.1f} mm at knee {d['knee_angle_r']:.1f}°")
        axes[row,1].set_title(f"Hip flexion {d['hip_flexion_r']:.1f}°, adduction {d['hip_adduction_r']:.1f}°")
        axes[row,0].set_ylabel("Superior +y (mm)")
    fig.suptitle("Gracilis (indigo) against the previous vasti bundle (gold)",fontweight="bold")
    fig.savefig(out/"grac_vasti_collision.png",dpi=175)
    plt.close(fig)


if __name__=="__main__":
    run(*map(Path,sys.argv[1:]))
