"""Replay either fixed-parameter Gait2392 controller and record actuator states."""
import argparse
import json
import os
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["s3k", "spiking"])
    parser.add_argument("--spinal", type=Path, required=True,
                        help="Repository Code/MuJoCo_SNS/spinal directory")
    parser.add_argument("--output", type=Path, required=True,
                        help="New results root used for this pair of simulations")
    args = parser.parse_args()
    spinal = args.spinal.resolve()
    dest = args.output.resolve() / ("figpair_" + args.case)
    dest.mkdir(parents=True, exist_ok=False)
    for key in list(os.environ):
        if key.startswith("AARL_"):
            del os.environ[key]
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    os.environ["AARL_NPZ"] = str(dest / "output.npz")
    os.environ["MPLBACKEND"] = "Agg"
    if args.case == "spiking":
        os.environ["AARL_NET"] = "spiking"
    os.chdir(spinal)
    sys.path.insert(0, str(spinal))

    import mujoco
    import numpy as np
    import _curriculum as curriculum
    import params
    import runner

    multiplier_path = spinal / "best_walk_params_v10.json"
    winner_path = spinal / "reports_20260923" / "s3k_trial34_full_params.json"
    multipliers = json.loads(multiplier_path.read_text())["multipliers"]
    winner = json.loads(winner_path.read_text())["params"]
    curriculum.BASE_MUL = dict(multipliers)
    curriculum.BASE_MUL["renshaw"] = 0.5
    curriculum.set_stage(4, {**multipliers, **winner, "full_rules": 0.0})
    params.G["full_rules"] = 0.0
    muscle_names = ["vas_lat_r", "semimem_r"]
    step = mujoco.mj_step
    addresses = []
    activations = []
    physics_times = []

    def record_after_step(model, data, *positional, **keyword):
        step(model, data, *positional, **keyword)
        if not addresses:
            for name in muscle_names:
                actuator = model.actuator(name).id
                assert model.actuator_dyntype[actuator] == mujoco.mjtDyn.mjDYN_MUSCLE
                addresses.append(model.actuator_actadr[actuator])
        activations.append(data.act[addresses].copy())
        physics_times.append(data.time)

    mujoco.mj_step = record_after_step
    try:
        outcome = runner.main(["--eval", "--drive", repr(winner["drive"])])
    finally:
        mujoco.mj_step = step
    with np.load(dest / "output.npz", allow_pickle=False) as samples:
        times = samples["t"]
    dt = times[1] - times[0]
    held_steps = len(times) - len(physics_times)
    assert held_steps == round(0.6 / dt)
    schedule_times = np.asarray(physics_times) + times[held_steps]
    np.testing.assert_allclose(schedule_times, times[held_steps:] + dt, atol=1e-9)
    np.savez_compressed(dest / "true_activation.npz", time_label=times[held_steps:],
                        state_time=schedule_times, activation=np.asarray(activations),
                        muscle_names=muscle_names)
    results = {key: value for key, value in outcome.items()
               if isinstance(value, (int, float, bool)) or value is None}
    results.update({"k_" + key: value for key, value in (outcome.get("kine") or {}).items()
                    if isinstance(value, (int, float, bool)) or value is None})
    (dest / "result.json").write_text(json.dumps(results, indent=2), encoding="utf8")


if __name__ == "__main__":
    main()
