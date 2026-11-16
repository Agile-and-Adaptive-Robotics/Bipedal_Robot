"""Stock OpenSim signed targets for paired-coordinate optimization."""
import json
import sys
from pathlib import Path

import numpy as np
import opensim as osim

from fr_plotGracilisPilotSweeps import human_torque


def run(repo, out):
    path = repo/"Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json"
    atlas = json.loads(path.read_text())
    ref = repo/"Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim"
    osim.ModelVisualizer.addDirToGeometrySearchPaths("C:/OpenSim 4.6/Geometry")
    model = osim.Model(str(ref))
    state = model.initSystem()
    coords = model.getCoordinateSet()
    home = {coords.get(i).getName(): coords.get(i).getDefaultValue()
            for i in range(coords.getSize())}
    curve = osim.ActiveForceLengthCurve()
    out.mkdir(parents=True,exist_ok=True)
    groups = {
        "bifemlh": ["bifemlh_r"],
        "semimem": ["semimem_r","semiten_r"],
        "sar": ["sar_r"],
        "rect_fem": ["rect_fem_r"],
        "grac": ["grac_r"],
        "tfl": ["tfl_r"],
        "med_gas": ["med_gas_r"],
        "lat_gas": ["lat_gas_r"],
    }
    rom = atlas["rom_table_deg"]
    for group,names in groups.items():
        calf = group.endswith("_gas")
        hip="hip_adduction_r" if group in ("grac","tfl") else "hip_flexion_r"
        a,b = ("knee_angle_r","ankle_angle_r") if calf else (hip,"knee_angle_r")
        for label,n in (("search",5),("dense",13)):
            limits = lambda k: (-120,0) if k=="knee_angle_r" else rom[k]
            first=np.linspace(*limits(a),n)
            second=np.linspace(*limits(b),n)
            poses=[dict(home,**{a:np.deg2rad(x),b:np.deg2rad(y)}) for y in second for x in first]
            torques=[]
            for coord in (a,b):
                torque=sum(human_torque(model,state,model.getMuscles().get(name),curve,coord,poses)
                           for name in names)
                torques.append(torque.tolist())
            record=dict(group=group,coordinate1=a,coordinate2=b,axis1_deg=first.tolist(),
                        axis2_deg=second.tolist(),target1_Nm=torques[0],target2_Nm=torques[1])
            (out/f"{group}_{label}.json").write_text(json.dumps(record),encoding="utf-8")
        print(f"{group}: stock signed targets on both coordinates")
    knee=np.linspace(-120,0,61)
    poses=[dict(home,knee_angle_r=np.deg2rad(angle)) for angle in knee]
    targets={}
    for name in ("vas_int_r","vas_lat_r","vas_med_r"):
        muscle=model.getMuscles().get(name)
        point=muscle.getGeometryPath().getPathPointSet().get(0).getLocation(state)
        targets[name]=dict(origin_m=[point.get(i) for i in range(3)],knee_deg=knee.tolist(),
                           target_Nm=human_torque(model,state,muscle,curve,"knee_angle_r",poses).tolist())
    (out/"vasti_individual.json").write_text(json.dumps(targets),encoding="utf-8")


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]))
