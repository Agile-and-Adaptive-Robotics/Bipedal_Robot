"""Tune nearby ground-recruitment settings, then hold the motor currents fixed across support and design."""
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import json
from pathlib import Path

import robotbody_curriculum as R
from refine_robotbody_curriculum import evaluate, ordering
from tune_robotbody_curriculum import dump, run_fresh
import numpy as np
from scipy.stats import qmc


def physical_priority(row):
    if not row["complete"]:
        return (False, False, 0, row["completed_s"])
    swings = [len(leg) for leg in row["loaded_clear_swings"]]
    gates = int(row["movement"]["accepted"]) + int(row["contact"]["accepted"]) + int(min(swings) >= 3)
    return (row["accepted"], True, gates, ordering(row)[-1])


def ground_seeds(folder, design):
    rows = [json.loads(p.read_text()) for p in sorted(folder.glob(f"{design}_s4_*.json"))]
    eligible = [row for row in rows if row["complete"]]
    if not eligible:
        raise ValueError(f"No complete {design} ground trial supplies a repeatable starting configuration.")
    timed = max(
        eligible,
        key=lambda r: (
            r["contact"]["accepted"] and min(map(len, r["loaded_clear_swings"])) >= 3,
            r["contact"]["accepted"],
            *physical_priority(r),
        ),
    )
    motion = max(eligible, key=lambda r: (r["movement"]["accepted"], *physical_priority(r)))
    return [timed, motion]


def parameters(seed, u):
    p = seed["parameters"].copy()
    excessive_hip = max(seed["movement"]["hip_ranges_deg"]) > 80.0
    if excessive_hip:
        hip_flex, hip_ext = 0.46 + 0.48 * u[0], 0.46 + 0.48 * u[1]
    else:
        hip_flex, hip_ext = 0.90 + 0.24 * u[0], 0.88 + 0.23 * u[1]
    p["hip_flex_scale"] = p.get("hip_flex_scale", 1.0) * hip_flex
    p["hip_ext_scale"] = p.get("hip_ext_scale", 1.0) * hip_ext
    p["ankle_df_scale"] = p.get("ankle_df_scale", 1.0) * (0.70 + 1.30 * u[2])
    p["ankle_pf_scale"] = p.get("ankle_pf_scale", 1.0) * (0.58 + 0.62 * u[3])
    p["flex"] *= 0.94 + 0.36 * u[4]
    p["ext"] *= 0.92 + 0.20 * u[5]
    p["support4"] = float(np.clip(p["support4"] - 0.035 + 0.065 * u[6], 0.46, 0.79))
    return p


def standing_trial(design, parameters, atlas, offset, path):
    row = R.evaluate(design, 3, parameters, 0.0005, atlas, pose_offset=offset)
    row["initial_knee_offset_deg"] = offset
    dump(path, row)
    return row


def check_standing(args, atlas, folder, checks):
    """Tune postural sensory gain across two starts, then test unseen starts."""
    for design in ("baseline", "selected"):
        source = args.refined / "checks" / f"{design}_s3_rank0_frozen.json"
        p = json.loads(source.read_text())["parameters"]
        gains = (0.0, 0.025, 0.075, 0.18, 0.45, 0.9)
        rows = {g: [] for g in gains}
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {}
            for j, gain in enumerate(gains):
                for offset in (-0.5, 0.5):
                    path = folder / f"{design}_standing_{j}_{offset:+.1f}.json"
                    par = {**p, "postural_feedback_gain": gain}
                    futures[pool.submit(standing_trial, design, par, atlas, offset, path)] = gain
            for f in as_completed(futures):
                rows[futures[f]].append(f.result())
        winner = max(
            gains, key=lambda g: (min((r["accepted"], r["complete"], r["score"]) for r in rows[g]), -g)
        )
        selected = {**p, "postural_feedback_gain": winner}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            options = (
                ("frozen", 0.0005, 0.0, selected),
                ("previous_holdout", 0.0005, 1.25, selected),
                ("new_holdout", 0.0005, -1.75, selected),
                ("halfstep", 0.00025, 0.0, selected),
                ("stiffer_contact", 0.0005, 0.0, {**selected, "contact_time_s": 0.010}),
            )
            jobs = {
                pool.submit(
                    run_fresh, checks, 3, design, par, "balance_" + tag, args.atlas, dt, 0.01, offset
                ): tag
                for tag, dt, offset, par in options
            }
            final = {
                jobs[f]: {k: r[k] for k in ("accepted", "completed_s", "reason")}
                for f in as_completed(jobs)
                for r in (f.result(),)
            }
        dump(
            args.output / (design + "_standing.json"),
            {"parameters": selected, "training_trials": 2 * len(gains), "checks": final},
        )
        print("standing transfers", design, winner, final, flush=True)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--refined", type=Path, required=True)
    cli.add_argument("--atlas", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--trials", type=int, default=18, help="Additional local trials per design.")
    cli.add_argument("--workers", type=int, default=4)
    args = cli.parse_args()
    if args.trials < 2 or not 1 <= args.workers <= 4:
        cli.error("Use at least two trials and at most four workers.")
    args.output.mkdir(parents=True, exist_ok=True)
    jobs_folder, checks = args.output / "trials", args.output / "checks"
    jobs_folder.mkdir(exist_ok=True)
    checks.mkdir(exist_ok=True)
    atlas = json.loads(args.atlas.read_text())
    designs = ("baseline", "selected")
    rows = {d: [] for d in designs}
    samples = qmc.Sobol(7, scramble=True, seed=17331).random_base2(int(np.ceil(np.log2(args.trials))))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = []
        for design in designs:
            seeds = ground_seeds(args.refined, design)
            rows[design].extend(seeds)
            for i, u in enumerate(samples[: args.trials]):
                seed = seeds[i % 2]
                p = parameters(seed, u)
                if i % 3 == 1:
                    p["postural_feedback_gain"] = 0.015 + 0.16 * u[0]
                if i % 3 != 0:
                    p["load_to_vastus_uS"] = 0.20 + 3.8 * u[1]
                label = f"{seed['variant']}_local_{i:02d}"
                file = jobs_folder / f"{design}_{label}.json"
                if file.exists():
                    existing = json.loads(file.read_text())
                    if existing["parameters"] != p:
                        raise ValueError(f"Different saved trial: {file}")
                    rows[design].append(existing)
                else:
                    jobs.append(pool.submit(evaluate, (design, 4, p, 0.0005, atlas, file, label)))
        for future in as_completed(jobs):
            result = future.result()
            rows[result["design"]].append(result)
            print(
                result["design"], result["variant"], result["completed_s"], result["accepted"], flush=True
            )
    for design in designs:
        ordered = sorted(rows[design], key=physical_priority, reverse=True)
        # A rank is selected using training only, before these different starting
        # poses, timesteps, contacts, or body configurations are tested.
        row = ordered[0]
        p = row["parameters"]
        bias = row["static_neural_bias"]
        currents = bias.get("frozen_motor_current_nA")
        if currents is None:
            currents = [p["static_bias"] * x for x in bias["current_nA"]]
        frozen = {**p, "frozen_post_currents_by_route_nA": dict(zip(atlas["names"], currents))}
        other = "baseline" if design == "selected" else "selected"
        settings = [
            ("replay", design, 4, p, 0.0005, 0.0),
            ("holdout", design, 4, p, 0.0005, -1.75),
            ("prior_holdout", design, 4, p, 0.0005, 1.25),
            ("halfstep", design, 4, p, 0.00025, 0.0),
            ("stiffer_contact", design, 4, {**p, "contact_time_s": 0.010}, 0.0005, 0.0),
            ("frozen_current", design, 4, frozen, 0.0005, 0.0),
            ("same_current_other_design", other, 4, frozen, 0.0005, 0.0),
            ("less_lift", design, 4, {**frozen, "ground_lift_override": p["support4"] - 0.05}, 0.0005, 0.0),
            (
                "less_lift_toe_feedback",
                design,
                5,
                {**frozen, "support5": p["support4"] - 0.05, "layer_contacts": True},
                0.0005,
                0.0,
            ),
            ("zero_lift", design, 4, {**frozen, "ground_lift_override": 0.0}, 0.0005, 0.0),
            (
                "no_rig",
                design,
                4,
                {**frozen, "ground_lift_override": 0.0, "guidance_scale": 0.0},
                0.0005,
                0.0,
            ),
        ]
        if row["accepted"]:
            jitter = np.random.default_rng(7131)
            varied = (
                "drive",
                "tau",
                "posture",
                "ext",
                "flex",
                "hip",
                "ankle",
                "static_bias",
                "stretch",
                "speed",
                "tension",
            )
            for index in range(8):
                candidate = {**p, **{key: p[key] * (1.0 + jitter.uniform(-0.01, 0.01)) for key in varied}}
                settings.append((f"one_percent_{index}", design, 4, candidate, 0.0005, 0.0))
        reports = {}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            jobs = {
                pool.submit(
                    run_fresh, checks, stage, d, test, design + "_" + tag, args.atlas, dt, 0.01, offset
                ): tag
                for tag, d, stage, test, dt, offset in settings
            }
            for future in as_completed(jobs):
                tag = jobs[future]
                test = future.result()
                reports[tag] = {
                    k: test[k]
                    for k in (
                        "design",
                        "stage",
                        "completed_s",
                        "accepted",
                        "reason",
                        "measured_support_BW",
                        "mean_foot_load_BW",
                    )
                }
                print(design, tag, reports[tag], flush=True)
        dump(
            args.output / (design + "_winner.json"),
            {
                "source": row["variant"],
                "parameters": p,
                "fixed_currents_nA": frozen["frozen_post_currents_by_route_nA"],
                "additional_trials": args.trials,
                "checks": reports,
            },
        )
    check_standing(args, atlas, jobs_folder, checks)


if __name__ == "__main__":
    main()
