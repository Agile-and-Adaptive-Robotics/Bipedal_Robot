"""Compare Appendix C vastus bundles at the recorded pose of a CAD allocator trial."""
import argparse
import json
from pathlib import Path
import sys

import mujoco
import numpy as np

from robotbody_survey import load_design

sys.path.insert(0, str(Path(__file__).resolve().parent / "robot_torso/bpa_neural_coordination"))
from actuator_interface import force_at_pressure, packed_parameters


def force_parameters(row, reference):
    """Combine the separate sleeve brackets and the single shared output cable."""
    if int(row["G_survey"]) != 1:
        raise ValueError("The documented vastus comparison requires unity transmission.")
    count = int(row["nBPA"])
    rest = float(row["resting_length_m"])
    free = float(row["kmax"])
    fixed = float(row["tendon_length_m"]) + 2 * float(row["fitting_length_m"]) + float(row["Xi0_m"])
    stiffness = 1 / (float(reference["bracket_single_m_N"]) / count + 1 / float(reference["shared_tendon_N_m"]))
    packed = packed_parameters(int(row["diameter_saved_mm"]), rest, 1 - free / rest, fixed, stiffness)
    packed[3] *= count
    return packed


def compare(record, route, time, window, previous=None):
    run = json.loads(record.read_text())
    model = mujoco.MjModel.from_xml_path(run["model_path"])
    with np.load(record.with_suffix(".npz")) as trace:
        times = trace["time"]
        select = (window[0] <= times) & (times < window[1])
        match = np.flatnonzero(select & np.isclose(times, time, rtol=0, atol=1e-8))
        if len(match) != 1:
            raise ValueError("Choose a saved timestamp within the stated comparison window.")
        i = match[0]
        aid = model.actuator(route).id
        if run["actuator_names"][aid] != route:
            raise ValueError("Recorded actuator order differs from its model.")
        reference_force = -model.actuator_gainprm[aid, 0]
        ideal = -float(trace["forces"][i, aid])
        if not np.isclose(ideal, trace["activation"][i, aid] * reference_force):
            raise ValueError("Recorded forces do not match the ideal allocator.")
        demand = trace["desired_activation"][i, aid] * reference_force
        route_length = float(trace["route_lengths"][i, model.actuator_trnid[aid, 0]])
        route_ids = [model.actuator_trnid[model.actuator("vas_med_" + side).id, 0] for side in ("l", "r")]
        lengths = trace["route_lengths"][select][:, route_ids]
    center = .5 * (lengths.min() + lengths.max())
    results = dict(route=route, time_s=time, window_s=list(window), route_m=route_length,
                   recorded_request_N=float(demand), recorded_filtered_force_N=ideal, pressure_kPa=620.)
    if previous:
        old = json.loads(previous.read_text())["actuators"][route]
        old_span = old["constant_remote_path_offset_m"]
        p = packed_parameters(old["diameter_mm"], old["resting_length_m"], old["maximum_contraction_fraction"],
                              old["nonmuscle_length_m"], old["series_stiffness_N_m"])
        sensitivities = {}
        for count, k_factor, name in ((1, 1, "one"), (4, 1, "four_shared_compliance"),
                                     (4, 4, "four_independent_compliances")):
            shape = p.copy()
            shape[3] *= count
            shape[4] *= k_factor
            sensitivities[name] = force_at_pressure(route_length + old_span, 620., shape)
        results["earlier_specimen_N"] = sensitivities
    for design in ("baseline", "selected"):
        row, ref = load_design(design)["vas_med"]
        params = force_parameters(row, ref)
        center_span = params[2] + params[0] * (1 - .5 * params[1]) - center
        spans = {"direct": 0., "centered_transfer": center_span}
        if previous:
            spans["earlier_span"] = old_span
        capacities = {}
        for name, span in spans.items():
            force = force_at_pressure(route_length + span, 620., params)
            capacities[name] = dict(extra_fixed_span_m=float(span),
                                    available_N=float(force) if np.isfinite(force) else None)
        results[design] = dict(count=int(row["nBPA"]), diameter_mm=int(row["diameter_saved_mm"]),
                               resting_m=params[0], free_m=float(row["kmax"]), maximum_force_N=params[3],
                               equivalent_stiffness_N_m=params[4], capacities=capacities)
    return results


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("run", type=Path, help="CAD run JSON beside its matching NPZ.")
    cli.add_argument("--route", choices=("vas_med_l", "vas_med_r"), default="vas_med_l")
    cli.add_argument("--time", type=float, default=21.61)
    cli.add_argument("--window", type=float, nargs=2, default=(16., 28.),
                     help="Bilateral length range used to center the additional fixed span.")
    cli.add_argument("--previous-configuration", type=Path, help="Optional earlier short-sleeve configuration JSON.")
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = compare(args.run, args.route, args.time, args.window, args.previous_configuration)
    report = json.dumps(data, indent=2, allow_nan=False)
    if args.output:
        args.output.write_text(report, encoding="utf-8")
    print(report)
