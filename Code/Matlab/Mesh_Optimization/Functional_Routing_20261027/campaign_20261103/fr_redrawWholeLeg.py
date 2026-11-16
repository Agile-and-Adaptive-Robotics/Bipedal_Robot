"""Integrate the three vasti and guided gracilis into the full-leg preview."""
import copy
import itertools
import sys
from pathlib import Path

import numpy as np

from fr_exportPilotOpenSim import export_routes
from fr_plotSplitVastiGrac import NAMES, check_pairs, group_parts, routes_of
from fr_plotWholeLegPack import draw
from fr_wholeLegPack import Envelope, clearance, read, table, unpack, write


def shift_medialis(result, pose, shift):
    result=copy.deepcopy(result)
    if shift==0:
        return result
    fem=np.asarray(pose["transforms"]["femur_r"])[:,:3,:3]
    tib=np.asarray(pose["transforms"]["tibia_r"])[:,:3,:3]
    direction=np.einsum("nji,nj->ni",tib,fem[:,:,2])
    # A translation along the hinge axis leaves all cam cable lengths invariant.
    if np.max(np.abs(direction-direction[0]))>1e-7:
        return None
    v=result["vasti"][0]
    delta=np.array([0,0,shift])
    for key in ("points_m","home_cam_femur_m"):
        v[key]=(np.asarray(v[key])+delta).tolist()
    v["contactPath_m"]=[(np.asarray(p)+delta).tolist() for p in v["contactPath_m"]]
    v["terminal_m"]=(np.asarray(v["terminal_m"])+shift*direction[0]).tolist()
    v["bone"]-=abs(shift)
    v["cableBone"]-=abs(shift)
    return result


def choose_knee(studies,branch,pose):
    best=None
    for study in studies:
        f=study/f"split_{branch}.json"
        if not f.is_file():
            continue
        for shift in (0,.006,.012,.018):
            r=shift_medialis(read(f),pose,shift)
            if r is None:
                continue
            routes=routes_of(r)
            pair=check_pairs(pose,r,routes,.02105)
            gaps=np.array([[x[k] for k in ("sleeve_mm","hardware_mm","cable_mm")] for x in pair])
            g=r["grac"]
            local=[g["bone"],g["bay"],g["viaBone"]]
            for v in r["vasti"]:
                local.extend((v["bone"],v["bay"],v["cableBone"]))
            penalty=60*np.maximum(0,-gaps[:,0]).sum()+20*np.maximum(0,-gaps[:,1]).sum()
            penalty+=8*np.maximum(0,-gaps[:,2]).sum()+10000*np.maximum(0,-np.asarray(local)).sum()
            penalty+=1000*g["noForce"]+160*g["commonNoForce"]+70*(1-g["coverage"][1])
            penalty+=20*(1-g["signedCapacity"][1])+15*(1-g["coverage"][0])
            penalty+=100*abs(shift)
            if best is None or penalty<best[0]:
                best=(penalty,r,routes,study.name,shift)
    return best[1:]


def pair_table(options):
    shapes={g:Envelope(o[1]) for g,o in options.items()}
    rows=[]
    groups={g:{"group":g,"sleeve_mm":18.,"hardware_mm":18.,"cable_mm":18.} for g in shapes}
    for a,b in itertools.combinations_with_replacement(shapes,2):
        gap,hit=clearance(shapes[a],shapes[b],same=a==b)
        keys=("sleeve_mm","hardware_mm","cable_mm")
        rows.append(dict(group_a=a,group_b=b,**dict(zip(keys,gap)),
                         sleeve_pose=hit[0],hardware_pose=hit[1],cable_pose=hit[2]))
        for g in (a,b):
            for key,value in zip(keys,gap):
                groups[g][key]=min(groups[g][key],value)
    for g,r in groups.items():
        current=options[g][1]
        for key in ("bone_mm","foot_lead_mm","home_side_mm","sleeve_bay_mm","tendon_span_mm","foot_occupancy_mm"):
            r[key]=current.get(key)
        r["clearance_pass"]=all(v is None or v>=-1e-6 for k,v in r.items() if k.endswith("_mm"))
    return list(groups.values()),rows


def main(out,original,study):
    pose=read(out/"dense"/"poses.json")
    old_pose=read(original/"dense"/"poses.json")
    np.testing.assert_allclose(pose["angles_deg"],old_pose["angles_deg"],atol=1e-9)
    selection=read(out/"selection.json")
    notes=["# Whole-leg redraw with the separate vasti groups","",
           "Starting layout is the saved pre-coordination assembly. Updated layout uses "
           "the selected mount positions, side-by-side sleeve rows, three vasti groups, "
           "and the guided gracilis route. The posterior view places anatomical right "
           "on the right. Three vasti groups replace the single combined entry, so the "
           "redraw contains 29 functional routing entries.",""]
    studies=[study/p for p in ("flat_rows_refined","final_bypass","final_bypass_fixed_rest")]
    for branch in ("rigid","compliant"):
        options={g:{1:opt[selection[f"{g}_{branch}"]]} for g,opt in unpack(out/"dense",branch).items()
                 if g not in ("grac","vas_int")}
        os_routes=[(f"{g}_{r['channel']}",r,(0,.22,.87)) for g,opt in options.items()
                   for r in opt[1]["lines"]]
        knee,knee_routes,source,shift=choose_knee(studies,branch,pose)
        write(out/f"knee_used_{branch}.json",knee)
        parts=group_parts(pose,knee,knee_routes,.03375)
        for name,record in zip((*NAMES,"Gracilis"),(*knee["vasti"],knee["grac"])):
            options[name]={1:dict(parts=parts[name],bone_mm=1000*record["bone"],
                                   sleeve_bay_mm=1000*record["bay"])}
        os_routes.extend((n,r,tuple(int(c[i:i+2],16)/255 for i in (1,3,5))) for n,r,c in knee_routes)
        export_routes(pose["reference"],out/f"whole_leg_split_vasti_{branch}.osim",
                      Path("C:/OpenSim 4.6/Geometry"),os_routes)
        groups,pairs=pair_table(options)
        table(out/f"{branch}_redraw_group_clearance.csv",groups)
        table(out/f"{branch}_redraw_pair_clearance.csv",pairs)
        large=check_pairs(pose,knee,knee_routes,.03375)
        small=check_pairs(pose,knee,knee_routes,.02105)
        for a,b in zip(large,small):
            a["hardware_38mm_sheave_proxy_mm"]=b["hardware_mm"]
        table(out/f"{branch}_local_knee_clearance.csv",large)
        clear=[r["group"] for r in groups if r["clearance_pass"]]
        g=knee["grac"]
        notes.extend([f"## {branch.capitalize()}","",
                      f"Checked {len(pose['angles_deg'])} simultaneous samples. "
                      f"{len(clear)} of {len(groups)} routing entries clear the modeled checks. "
                      f"Gracilis meets {g['coverage'][0]:.1%} of hip targets and "
                      f"{g['coverage'][1]:.1%} of knee targets on the coupled grid. "
                      f"The local candidate came from {source}. Medialis axial shift along "
                      f"the knee rotation axis is {1000*shift:.0f} mm.",""])
        if branch=="compliant":
            before=unpack(original/"dense",branch)
            draw(out,options,dict.fromkeys(options,1),pose,groups,before,
                 ("Starting layout","Updated layout\nSplit vasti and guided gracilis"))
    notes.extend(["The torso and pelvis clearance limitations and knee/ankle CAD "
                  "registration limits from the previous report still apply. Hardware "
                  "uses packaging envelopes and the 38.1 mm alternative is a local footprint "
                  "sensitivity. OpenSim paths show the saved native frames at home."])
    (out/"WHOLE_LEG_REDRAW_REPORT.md").write_text("\n".join(notes)+"\n",encoding="utf-8")


if __name__=="__main__":
    main(*map(Path,sys.argv[1:]))
