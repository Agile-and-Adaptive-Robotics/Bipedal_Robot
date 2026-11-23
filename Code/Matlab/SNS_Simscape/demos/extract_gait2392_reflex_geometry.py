"""Export lengths and signed knee moment arms from the original Gait2392 paths."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import xml.etree.ElementTree as ET

import opensim as osim


def extract(source, target):
    osim.Logger.setLevelString("error")
    model = osim.Model(str(source))
    state = model.initSystem()
    coordinate = model.getCoordinateSet().get("knee_angle_r")
    muscles = ("vas_med_r", "bifemsh_r")
    data = {"source_file": source.name, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "opensim": osim.GetVersion(), "flexion_deg": list(range(96)), "muscles": {},
            "setup": "knee_angle_r set to minus plotted flexion; remaining coordinates at model defaults"}
    xml = ET.parse(source)
    for name in muscles:
        m = model.getMuscles().get(name)
        x = xml.find(f".//Thelen2003Muscle[@name='{name}']")
        rec = {"Fmax": m.getMaxIsometricForce(), "l_opt": m.getOptimalFiberLength(),
               "l_slack": m.getTendonSlackLength(),
               "alpha_deg": math.degrees(m.getPennationAngleAtOptimalFiberLength()),
               "v_max": m.getMaxContractionVelocity(), "length_m": [], "length_gradient_m_per_rad": []}
        for f in ("activation_time_constant", "deactivation_time_constant", "FmaxMuscleStrain",
                  "KshapePassive", "KshapeActive", "Af", "Flen", "FmaxTendonStrain",
                  "fv_linear_extrap_threshold", "minimum_activation"):
            element = x.findtext(f)
            thelen = osim.Thelen2003Muscle.safeDownCast(m)
            rec[f] = float(element) if element is not None else float(getattr(thelen, "get_" + f)())
        for ang in data["flexion_deg"]:
            coordinate.setValue(state, -math.radians(ang), True)
            model.realizePosition(state)
            rec["length_m"].append(m.getGeometryPath().getLength(state))
            # dL/dq_flex = moment arm for source q_extension.
            rec["length_gradient_m_per_rad"].append(m.computeMomentArm(state, coordinate))
        data["muscles"][name] = rec
        print(name, "L and dL/dq at 45 deg", rec["length_m"][45], rec["length_gradient_m_per_rad"][45])
    test, speed_checks = [], []
    for angle in (30, 45, 60):
        coordinate.setValue(state, -math.radians(angle), True)
        for activation in (0.06, 0.12, 0.24):
            for name in muscles:
                model.getMuscles().get(name).setActivation(state, activation)
            model.equilibrateMuscles(state)
            model.realizeDynamics(state)
            for name in muscles:
                m = model.getMuscles().get(name)
                test.append({"angle": angle, "activation": activation, "muscle":name,
                             "fiber_length_m": m.getFiberLength(state), "force_N": m.getTendonForce(state),
                             "fiber_velocity_m_s": m.getFiberVelocity(state)})
            for name in muscles:
                m = osim.Thelen2003Muscle.safeDownCast(model.getMuscles().get(name))
                at_rest = m.getFiberLength(state)
                for change in (-0.00015, 0.00015):
                    m.setFiberLength(state, at_rest+change)
                    model.realizeDynamics(state)
                    speed_checks.append({"angle":angle,"activation":activation,"muscle":name,
                      "fiber_length_m":m.getFiberLength(state),"force_N":m.getTendonForce(state),
                      "fiber_velocity_m_s":m.getFiberVelocity(state)})
                m.setFiberLength(state, at_rest)
    data["opensim_equilibria"] = test
    data["opensim_dynamic_states"] = speed_checks
    target.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    args = parser.parse_args()
    extract(args.source, args.target)
