"""Repeat shortlisted robotbody controllers at a finer step and tune their poorest gates."""
import argparse
import concurrent.futures as cf
import hashlib
import json
import multiprocessing as mp
from pathlib import Path

import robotbody_curriculum as R
from tune_robotbody_curriculum import dump, run_fresh
import numpy as np
from scipy.stats import qmc


def ordering(r):
    if not r["complete"]:
        return (False, False, r["completed_s"])
    score = r["score"]
    if r["stage"] in (4, 5):
        movement = r["movement"]
        contact = r["contact"]
        score += 20 * int(contact["accepted"]) + 12 * int(movement["accepted"])
        score -= 0.8 * sum(max(0.0, x - 80.0) for x in movement["hip_ranges_deg"])
        score -= 1.0 * sum(max(0.0, x - 100.0) for x in movement["knee_ranges_deg"])
        score -= 40.0 * max(0.0, 0.8 - min(contact["coverage"])) + 60.0 * max(
            0.0, contact["phase_error"] - 0.15
        )
        score -= 8.0 * sum(max(0, 3 - len(leg)) for leg in r["loaded_clear_swings"])
    return (r["accepted"], True, score)


def evaluate(job):
    variant, stage, p, dt, atlas, path, tag = job
    result = R.evaluate(variant, stage, p, dt, atlas)
    result["variant"] = tag
    dump(path, result)
    return result


def candidates(stage, design, source, count, seed, seeds):
    opts = {
        d: sorted(
            [json.loads(p.read_text()) for p in sorted(source.glob(f"{d}_s{stage}_*.json"))],
            key=ordering,
            reverse=True,
        )
        for d in ("baseline", "selected")
    }
    other = "selected" if design == "baseline" else "baseline"
    # Always run the same unmodified settings on each design for the direct
    # transfer comparison before making a new fit.
    items = [("own_frozen", opts[design][0]["parameters"]), ("other_frozen", opts[other][0]["parameters"])]
    for j in range(1, min(5, len(opts[design]))):
        items.append((f"own_candidate_{j}", opts[design][j]["parameters"]))
    promising = opts[design][:3] + opts[other][:2]
    if stage in (4, 5):
        benchmark = sorted(opts["baseline"] + opts["selected"], key=ordering, reverse=True)
        promising = benchmark[:4] + opts[design][:1]
        for label, previous in seeds.items():
            p = previous.copy()
            if stage == 5:
                p["support5"] = max(0.0, p["support4"] - 0.05)
            items.append(("pilot_" + label, p))
        promising = [{"parameters": p} for p in seeds.values()] + promising
    if stage == 5:
        ground = sorted(
            [json.loads(path.read_text()) for path in source.glob("*_s4_*.json")],
            key=ordering,
            reverse=True,
        )
        best = ground[0]["parameters"].copy()
        best["support5"] = max(0.0, best["support4"] - 0.05)
        items.append(("frozen_ground_lower_lift", best))
        promising = ground[:3] + promising[:3]
    samples = qmc.Sobol(7, scramble=True, seed=seed + stage).random_base2(int(np.ceil(np.log2(count))))
    ankle = qmc.Sobol(2, scramble=True, seed=seed + 100 + stage).random_base2(int(np.ceil(np.log2(count))))
    for i, u in enumerate(samples[: max(0, count - len(items))]):
        p = promising[i % len(promising)]["parameters"].copy()
        if stage in (1, 2):
            p.update(
                trunk_ext_scale=0.015 + 0.43 * u[0] ** 2,
                trunk_flex_scale=0.25 + 2.4 * u[1],
                static_bias=0.35 + 0.95 * u[2],
                posture=0.02 + 0.38 * u[3],
                desc_f=1.6 + 0.8 * u[4],
                tau=0.20 + 0.27 * u[5],
                drive=2.6 + 0.85 * u[6],
            )
            p["calibrate_initial_motor_currents"] = True
            # These three optional datum choices change only the virtual
            # mounting offsets. Every source sleeve and tendon retains its size.
            if i % 4 != 0:
                p["trunk_reference_contraction"] = (0.35, 0.50, 0.65)[i % 3]
            p["trunk_voltage_mV"] = {
                "ercspn": 0.1 + 1.2 * u[0],
                "intobl": 0.35 + 1.8 * u[1],
                "extobl": 0.35 + 1.65 * u[2],
            }
        elif stage == 3:
            p.update(
                static_bias=0.7 + 0.58 * u[0],
                posture=0.01 + 0.38 * u[1],
                trunk_ext_scale=0.06 + 0.70 * u[2],
                trunk_flex_scale=0.35 + 2.2 * u[3],
                stretch=0.25 + 2.2 * u[4],
                speed=0.3 + 1.7 * u[5],
                vestibular=0.1 + 0.8 * u[6],
            )
            p["calibrate_initial_motor_currents"] = i % 3 != 0
        else:
            p.update(
                hip_flex_scale=0.20 + 0.85 * u[0],
                hip_ext_scale=0.20 + 0.85 * u[1],
                static_bias=p["static_bias"] * (0.82 + 0.38 * u[2]),
                ext=p["ext"] * (0.7 + 0.65 * u[3]),
                flex=p["flex"] * (0.75 + 0.5 * u[4]),
                drive=p["drive"] * (0.86 + 0.28 * u[5]),
            )
            p["calibrate_initial_motor_currents"] = i % 4 == 1
            if i % 2:
                p["ankle_pf_scale"] = 0.70 + 0.60 * ankle[i, 0]
                p["ankle_df_scale"] = 0.65 + 0.90 * ankle[i, 1]
            if stage == 5:
                shallow = p["support4"] - 0.03 - 0.10 * u[6]
                p["support5"] = max(
                    0.0, min(p["support4"] - 0.03, shallow if i % 3 == 0 else 0.32 + 0.32 * u[6])
                )
                p["layer_contacts"] = i % 2 == 0
            else:
                p["support4"] = np.clip(p["support4"] * (0.94 + 0.13 * u[6]), 0.45, 0.79)
        items.append((f"targeted_{i:02}", p))
    return items


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--initial", type=Path, required=True)
    cli.add_argument("--atlas", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--workers", type=int, default=4)
    cli.add_argument("--trials", type=int, default=22)
    cli.add_argument("--dt", type=float, default=0.0005)
    cli.add_argument("--stages", type=int, nargs="+", choices=R.DURATION, default=(4, 3, 1, 2, 5))
    cli.add_argument(
        "--ground-seeds",
        type=Path,
        default=Path(__file__).parent / "robotbody_curriculum_reference/ground_search_seeds.json",
    )
    args = cli.parse_args()
    if not 1 <= args.workers <= 4 or args.trials < 10:
        cli.error("Use 1–4 workers and at least ten trials.")
    args.output.mkdir(exist_ok=True, parents=True)
    checks = args.output / "checks"
    checks.mkdir(exist_ok=True)
    atlas = json.loads(args.atlas.read_text())
    seeds = json.loads(args.ground_seeds.read_text())["seeds"]
    manifest = {
        "dt": args.dt,
        "trials_per_design_stage": args.trials,
        "stages": list(args.stages),
        "seed": 7157,
        "hashes": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                Path(__file__),
                Path(R.__file__),
                args.atlas,
                args.ground_seeds,
                args.initial / "contract.json",
                args.initial / "progress.json",
            )
        },
    }
    manifest_path = args.output / "contract.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError("The output directory contains a different refinement experiment.")
    dump(manifest_path, manifest)
    progress = {}
    with cf.ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
        for s in args.stages:
            jobs = []
            results = {d: [] for d in ("baseline", "selected")}
            for d in results:
                items = candidates(s, d, args.initial / "trials", args.trials, 7157, seeds)
                for tag, p in items:
                    path = args.output / f"{d}_s{s}_{tag}.json"
                    if path.exists():
                        saved = json.loads(path.read_text())
                        if saved["parameters"] != p or saved["dt_s"] != args.dt:
                            raise ValueError(f"Existing trial used different settings: {path}")
                        results[d].append(saved)
                        continue
                    jobs.append((d, s, p, args.dt, atlas, path, tag))
            for future in cf.as_completed([pool.submit(evaluate, j) for j in jobs]):
                r = future.result()
                results[r["design"]].append(r)
                print(
                    r["design"],
                    s,
                    r["variant"],
                    r["completed_s"],
                    round(r["score"], 2),
                    r["accepted"],
                    flush=True,
                )
            for d, rows in results.items():
                rows.sort(key=ordering, reverse=True)
                # Recheck two finalists even if the highest nominal score fails.
                # Independently perturbed starts and half steps cannot enter training.
                for rank, row in enumerate(rows[:2]):
                    p = row["parameters"]
                    runs = {}
                    with cf.ThreadPoolExecutor(max_workers=args.workers) as fast:
                        futures = {}
                        for label, dt, perturb, parameters in (
                            ("frozen", args.dt, 0.0, p),
                            ("heldout", args.dt, 1.25, p),
                            ("halfstep", args.dt / 2, 0.0, p),
                            ("contact_sensitivity", args.dt, 0.0, {**p, "contact_time_s": 0.010}),
                        ):
                            tag = f"rank{rank}_{label}"
                            if s < 3 and label == "contact_sensitivity":
                                continue
                            f = fast.submit(
                                run_fresh, checks, s, d, parameters, tag, args.atlas, dt, 0.01, perturb
                            )
                            futures[f] = label
                        for f in cf.as_completed(futures):
                            runs[futures[f]] = f.result()
                    progress[f"{d}_{s}_rank{rank}"] = dict(
                        variant=row["variant"],
                        parameters=p,
                        robust=all(v["accepted"] for k, v in runs.items() if k != "contact_sensitivity"),
                        training_accepted=sum(v["accepted"] for v in rows),
                        completed=sum(v["complete"] for v in rows),
                        trials=len(rows),
                        checks={
                            k: {
                                a: r[a]
                                for a in (
                                    "complete",
                                    "accepted",
                                    "reason",
                                    "score",
                                    "completed_s",
                                    "measured_support_BW",
                                    "max_penetration_m",
                                )
                            }
                            for k, r in runs.items()
                        },
                    )
                    dump(args.output / "progress.json", progress)
                    print(
                        "verified",
                        d,
                        s,
                        rank,
                        {k: (r["accepted"], r["completed_s"]) for k, r in runs.items()},
                        flush=True,
                    )


if __name__ == "__main__":
    main()
