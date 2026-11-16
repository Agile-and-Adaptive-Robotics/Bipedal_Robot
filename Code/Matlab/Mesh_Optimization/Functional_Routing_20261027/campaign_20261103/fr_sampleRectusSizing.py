"""Native OpenSim path excursion and a straight sleeve bay screening proxy."""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import opensim as osim


def main(folder, geometry):
    route = json.loads((folder/"rectus_route.json").read_text())
    rom_path = Path(__file__).resolve().parents[2]/"Routing27_v2"/"actuator_map_27.json"
    limits = json.loads(rom_path.read_text())["rom_table_deg"]["hip_flexion_r"]
    hip = np.linspace(*limits, 51)
    knee = np.linspace(-120, 0, 51)
    osim.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry))
    human, robot = (osim.Model(route[key]) for key in ("stock", "robot"))
    hs, rs = human.initSystem(), robot.initSystem()
    rectus = human.getMuscles().get("rect_fem_r")
    human_knee = human.getCoordinateSet().get("knee_angle_r")
    human_hip = human.getCoordinateSet().get("hip_flexion_r")
    active_curve = osim.ActiveForceLengthCurve()
    femur = robot.getBodySet().get("femur_r")
    tibia = robot.getBodySet().get("tibia_r")
    paths = np.empty((knee.size, hip.size))
    rows = []
    for i, k in enumerate(knee):
        for j, h in enumerate(hip):
            for model, state in ((human, hs), (robot, rs)):
                coords = model.getCoordinateSet()
                coords.get("hip_flexion_r").setValue(state, np.deg2rad(h), False)
                coords.get("knee_angle_r").setValue(state, np.deg2rad(k), False)
                model.realizePosition(state)
            points = []
            for b, xyz in zip(route["bodies"], route["points_m"], strict=True):
                v = robot.getBodySet().get(b).findStationLocationInGround(rs, osim.Vec3(*xyz))
                points.append([v.get(c) for c in range(3)])
            segments = np.linalg.norm(np.diff(points, axis=0), axis=1)
            # Stop the first straight sleeve span at a knee plane normal to
            # the femur. The tibia origin is a proxy for this plane only.
            def ground(body, xyz):
                v = body.findStationLocationInGround(rs, osim.Vec3(*xyz))
                return np.array([v.get(c) for c in range(3)])
            normal = ground(femur,(0,1,0))-ground(femur,(0,0,0))
            knee_origin = ground(tibia,(0,0,0))
            a,b = np.array(points[:2])
            lead = np.dot(b-a,normal)
            fraction = np.clip(np.dot(knee_origin-a,normal)/lead,0,1) if lead < 0 else 1
            paths[i,j] = rectus.getLength(hs)
            bay = fraction*segments[0] - 2*route["endFitting_m"]
            arm = rectus.computeMomentArm(hs,human_knee)
            fiber = (paths[i,j]-rectus.getTendonSlackLength())/rectus.getOptimalFiberLength()
            force = rectus.getMaxIsometricForce()*active_curve.calcValue(float(fiber))
            rows.append([h,k,paths[i,j],segments.sum(),bay,arm,force*arm,
                         force*rectus.computeMomentArm(hs,human_hip)])
    with (folder/"rectus_pose_grid.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["hip_deg","knee_deg","human_path_m","robot_polyline_m","bay_m",
                         "human_knee_arm_m","human_knee_Nm","human_hip_Nm"])
        writer.writerows(rows)
    a = paths[-1,-1]  # hip flexed, knee extended
    b = paths[0,0]    # hip extended, knee flexed
    summary = {
        "human_fmax_N":rectus.getMaxIsometricForce(),
        "hip_deg":limits, "knee_deg":[-120,0],
        "hip_flexed_knee_extended_m":a, "hip_extended_knee_flexed_m":b,
        "human_path_min_m":float(paths.min()), "human_path_max_m":float(paths.max()),
        "robot_path_min_m":min(r[3] for r in rows), "robot_path_max_m":max(r[3] for r in rows),
        "bay_min_m":min(r[4] for r in rows), "bay_max_m":max(r[4] for r in rows),
        "bay_note":"First straight span clipped at the proxy knee plane, less two fittings; no collision qualification.",
    }
    (folder/"rectus_excursion.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    fig, ax = plt.subplots(figsize=(7.2,4.4), layout="constrained")
    mesh = ax.pcolormesh(hip,knee,1000*paths,cmap="viridis",shading="nearest")
    ax.scatter([hip[-1],hip[0]],[knee[-1],knee[0]],c="white",edgecolor="#212529",s=45)
    ax.set(title="Stock right rectus femoris path length", xlabel="Hip flexion angle (°)",
           ylabel="Knee angle (°) · flexion negative")
    fig.colorbar(mesh,ax=ax,label="Muscle–tendon path length (mm)")
    fig.savefig(folder/"rectus_path_excursion.png", dpi=200)


if __name__ == "__main__":
    main(Path(sys.argv[1]),Path(sys.argv[2]))
