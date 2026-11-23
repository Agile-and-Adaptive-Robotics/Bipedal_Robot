"""Compare a saved CAD tension request with a specified BPA at the same pose."""
import argparse
import json
from pathlib import Path

import mujoco
import numpy as np

from actuator_interface import force_at_pressure, packed_parameters


def at_time(run_file, configuration_file, route, time):
    run = json.loads(run_file.read_text())
    spec = json.loads(configuration_file.read_text())["actuators"][route]
    model = mujoco.MjModel.from_xml_path(run["model_path"])
    aid = model.actuator(route).id
    if run["actuator_names"][aid] != route:
        raise ValueError("Recorded channels differ from the model's actuator order.")
    tid = model.actuator_trnid[aid, 0]
    reference_capacity = -model.actuator_gainprm[aid, 0]
    with np.load(run_file.with_suffix(".npz")) as data:
        frames = np.flatnonzero(np.isclose(data["time"], time, rtol=0, atol=1e-8))
        if len(frames) != 1:
            raise ValueError("Choose a timestamp in the saved trajectory.")
        i = frames[0]
        ideal = -float(data["forces"][i, aid])
        if not np.isclose(ideal, data["activation"][i, aid] * reference_capacity):
            raise ValueError("This run does not use the expected ideal tension mapping.")
        request = data["desired_activation"][i, aid] * reference_capacity
        length = data["route_lengths"][i, tid]
    path = length + spec.get("constant_remote_path_offset_m", 0.)
    p = packed_parameters(spec["diameter_mm"], spec["resting_length_m"],
                          spec["maximum_contraction_fraction"], spec["nonmuscle_length_m"],
                          spec["series_stiffness_N_m"])
    p[3] *= spec["parallel_count"]
    pressure = spec["pressure_max_kpa"]
    available = force_at_pressure(path, pressure, p)
    if not np.isfinite(available):
        raise ValueError("The recorded route exceeds this actuator's force–length domain.")
    return dict(route=route, time_s=time, route_length_m=length, requested_N=request,
                filtered_ideal_N=ideal, pressure_limit_kPa=pressure, BPA_capacity_N=available)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="CAD run JSON beside its trajectory NPZ")
    parser.add_argument("configuration", type=Path)
    parser.add_argument("route")
    parser.add_argument("time", type=float)
    args = parser.parse_args()
    print(json.dumps(at_time(args.run, args.configuration, args.route, args.time), indent=2))
