"""Check a saved pneumatic trial against the MATLAB surface and scalar pressure bank."""
import argparse
import json
from pathlib import Path

import robotbody_curriculum as curriculum
import numpy as np
from cad_interface import RoutedBPABank


def verify(prefix):
    with np.load(str(prefix) + ".npz") as trace:
        names = trace["names"].tolist()
        time, right_ground_force = trace["t"], trace["grf"][:, 0]
        pressure, target, force, length = (
            np.asarray(trace[key], dtype=float) for key in ("pressure", "command", "tension", "length")
        )
    config = json.loads(Path(str(prefix) + ".installation.json").read_text())
    summary = json.loads(Path(str(prefix) + ".json").read_text())
    if len(set(names)) != 54 or not len(pressure):
        raise ValueError("Verification requires a nonempty trial with all 54 route channels.")
    entries = [config["actuators"][name] for name in names]

    def column(key):
        return np.array([e[key] for e in entries])

    gain = column("transmission_gain")
    offset = column("bpa_side_reference_m")
    ref = column("output_reference_m")
    path = offset + (length - ref) / gain
    rest = column("resting_length_m")
    stroke = rest - column("free_length_at_620kPa_m")
    passive_length = column("nonmuscle_length_m")
    cable_stretch = force * gain / column("series_stiffness_N_m")
    relative = (rest - path + passive_length + cable_stretch) / stroke
    source = curriculum.HERE.parent / "Matlab/Functions/festo_lookup_coeffs.json"
    matlab = json.loads(source.read_text())
    dia = column("diameter_mm")
    coeff = np.array(
        [
            [matlab[{10: "f_10", 20: "f20", 40: "f40"}[d]]["coeffs"][key] for key in ("a0", "a1", "a3")]
            for d in dia
        ]
    )
    a0, a1, a3 = coeff.T
    # Values in maxBPAforce.m. The constant 40 mm force is from Festo.
    cap10 = 620.0 * 0.4895 * np.arctan(0.03068 * (rest - 0.0075) * 620.0)
    cap20 = 620.0 * 1.4877 * np.arctan(0.0248 * (rest - 0.0075) * 620.0)
    capacity = np.select([dia == 10, dia == 20], [cap10, cap20], default=6398.4)
    fraction = a0 * np.expm1(-a1 * relative) + pressure / 620.0 * np.exp(-a3 * relative**2)
    expected = column("parallel_count") * capacity * np.maximum(0.0, fraction) / gain
    maximum_force = (
        np.maximum(0.0, path - passive_length - rest + stroke) * column("series_stiffness_N_m") / gain
    )
    at_stop = (np.abs(relative - 1.0) <= 1e-5) & (force > 0.02)
    expected = np.minimum(expected, maximum_force)
    slack = path - passive_length <= rest - stroke
    expected[slack] = 0.0
    delta = np.abs(expected - force)
    bad_index = np.unravel_index(np.argmax(delta), delta.shape)
    domain = (relative < -1e-5) | ((relative > 1.0 + 1e-5) & (force > 0.02))

    # Independently advance the scalar and batched valve banks with the same
    # recorded neural commands and length changes, from identical zero pressure.
    # This tests histories with both filling and venting. Limiting the first
    # second avoids an expensive duplicate dynamics run.
    bank = RoutedBPABank(names, config, 0.002)
    actual = bank.pressure.copy()
    batch_force = np.zeros(54)
    bank_force_error = bank_pressure_error = 0.0
    for i in range(min(500, len(pressure))):
        scalar = np.minimum(bank.step(target[i] / 620.0, path[i]) / gain, maximum_force[i])
        curriculum.advance_pressures(
            target[i],
            actual,
            length[i],
            offset,
            ref,
            gain,
            bank.parameters,
            bank.fill,
            bank.vent,
            0.002,
            batch_force,
        )
        bank_force_error = max(bank_force_error, float(np.max(np.abs(batch_force - scalar))))
        bank_pressure_error = max(bank_pressure_error, float(np.max(np.abs(actual - bank.pressure))))
    result = dict(
        samples=len(pressure),
        routes=len(names),
        minimum_tension_N=float(force.min()),
        maximum_tension_N=float(force.max()),
        pressure_range_kPa=[float(pressure.min()), float(pressure.max())],
        nonzero_neural_channels=int(np.count_nonzero(target.max(axis=0) > 0.0)),
        slack_sample_fraction=float(np.mean(slack)),
        force_domain_violations=int(domain.sum()),
        shortened_stop_fraction=float(np.mean(at_stop)),
        max_independent_surface_residual_N=float(delta[bad_index]),
        worst_route=names[bad_index[1]],
        float32_force_snapshots=True,
        scalar_bank_max_pressure_error_kPa=bank_pressure_error,
        scalar_bank_max_force_error_N=bank_force_error,
        max_generalized_force_sign_error=summary["max_force_balance_error_Nm"],
    )
    result["passed"] = bool(
        np.isfinite(pressure).all()
        and np.isfinite(force).all()
        and force.min() >= 0.0
        and pressure.min() >= 0.0
        and pressure.max() <= 620.0001
        and not domain.any()
        and delta.max() < 0.05
        and bank_force_error < 1e-8
        and bank_pressure_error < 1e-8
        and summary["max_force_balance_error_Nm"] < 1e-8
    )
    stance = (time >= 5.0) & (time < 13.0) & (right_ground_force > 0.05 * summary["weight_N"])
    if stance.any():
        right = {}
        for name in ("vas_med_r", "med_gas_r"):
            j = names.index(name)
            possible = np.array(
                [
                    curriculum.force_at_pressure(x, 620.0, bank.parameters[j]) / gain[j]
                    for x in path[stance, j]
                ]
            )

            def median_and_max(values):
                return [float(np.median(values)), float(np.max(values))]

            right[name] = {
                "samples": int(stance.sum()),
                "chamber_pressure_median_max_kPa": median_and_max(pressure[stance, j]),
                "output_tension_median_max_N": median_and_max(force[stance, j]),
                "capacity_at_620_kPa_median_max_N": median_and_max(possible),
            }
        result["right_stance_capacity"] = right
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trial", type=Path, help="Output prefix, without .npz or .json.")
    args = parser.parse_args()
    result = verify(args.trial)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
