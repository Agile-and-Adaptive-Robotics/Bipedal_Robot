"""Tune and replay every physical stage of the survey-sized pneumatic robotbody."""
import argparse
import concurrent.futures as cf
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys

import robotbody_curriculum as R
import numpy as np
from scipy.stats import qmc


RANGES = {
    1: {
        "drive": (1.3, 4.0),
        "tau": (0.16, 0.6),
        "desc_e": (0.95, 2.2),
        "desc_f": (1.05, 2.65),
        "rg_to_pf": (1.1, 3.8),
        "posture": (0.02, 0.4),
        "ext": (0.35, 2.5),
        "flex": (0.45, 3.5),
        "hip": (0.5, 2.0),
        "static_bias": (0.1, 1.05),
    },
    2: {
        "stretch": (0.05, 3.0),
        "speed": (0.1, 2.0),
        "tension": (0.1, 3.0),
        "commissural": (0.05, 1.4),
        "static_bias": (0.1, 1.2),
        "drive": (1.5, 4.0),
    },
    3: {
        "static_bias": (0.35, 1.8),
        "posture": (0.02, 1.2),
        "vestibular": (0.05, 1.2),
        "stretch": (0.15, 3.0),
        "speed": (0.2, 2.0),
        "tension": (0.05, 1.8),
    },
    4: {
        "drive": (1.4, 4.0),
        "posture": (0.02, 0.9),
        "ext": (0.3, 2.6),
        "flex": (0.3, 3.5),
        "hip": (0.4, 2.4),
        "ankle": (0.4, 2.2),
        "heel": (0.05, 1.1),
        "stance_feedback": (0.05, 0.9),
        "static_bias": (0.2, 1.4),
        "support4": (0.48, 0.77),
    },
    5: {
        "toe": (0.015, 0.85),
        "heel": (0.03, 0.9),
        "stance_feedback": (0.02, 1.2),
        "stretch": (0.1, 2.5),
        "static_bias": (0.2, 1.5),
        "support5": (0.20, 0.46),
    },
}


def dump(path, obj):
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    temporary.replace(path)


def worker(job):
    variant, stage, p, dt, atlas, reserve, index, out = job
    path = out / f"{variant}_s{stage}_{index:03}.json"
    if path.exists():
        return json.loads(path.read_text())
    result = R.evaluate(variant, stage, p, dt, atlas, reserve)
    result["trial"] = index
    dump(path, result)
    return result


def priority(row):
    return (row["accepted"], row["complete"], row["score"], -row["trial"])


def run_fresh(directory, stage, d, p, tag, atlas, dt, reserve, heldout=0.0):
    name = f"{d}_s{stage}_{tag}"
    output = directory / name
    config = directory / (name + ".parameters.json")
    dump(config, {"parameters": p})
    args = [
        sys.executable,
        str(Path(R.__file__)),
        "--stage",
        str(stage),
        "--design",
        d,
        "--config",
        str(config),
        "--atlas",
        str(atlas),
        "--output",
        str(output),
        "--dt",
        str(dt),
        "--reserve",
        str(reserve),
        "--pose-offset",
        str(heldout),
    ]
    with open(str(output) + ".log", "w") as f:
        subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, check=True)
    return json.loads(Path(str(output) + ".json").read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--trials", type=int, default=24)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--reserve", type=float, default=0.01)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4 or args.trials < 4:
        parser.error("Use 1 to 4 workers and at least four trials per design and stage.")
    folder = args.output
    trials = folder / "trials"
    checks = folder / "checks"
    for path in (trials, checks):
        path.mkdir(exist_ok=True, parents=True)
    atlas = json.loads(args.atlas.read_text())
    variants = ("baseline", "selected")
    base = {v: dict(R.DEFAULT) for v in variants}
    contract = dict(
        stages=R.DURATION,
        ranges=RANGES,
        trials_per_design_stage=args.trials,
        dt=args.dt,
        reserve=args.reserve,
        seed=4461,
        heldout_initial_knee_asymmetry_deg=1.25,
        atlas_sha256=hashlib.sha256(args.atlas.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(R.__file__).read_bytes()).hexdigest(),
        diagnostic_continuation="After a failed prerequisite, later stages remain independent diagnostics.",
    )
    manifest = folder / "contract.json"
    if manifest.exists() and json.loads(manifest.read_text()) != json.loads(json.dumps(contract)):
        raise RuntimeError("Existing experiment has a different model, design, or search contract.")
    dump(manifest, contract)
    stages = {}
    passed = {v: True for v in variants}
    with cf.ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
        for s, ranges in RANGES.items():
            keys = list(ranges)
            bounds = np.array(list(ranges.values()))
            n = args.trials
            n0 = max(3, n // 2)
            rng = np.random.default_rng(4461 + s)
            uniform = qmc.Sobol(len(keys), scramble=True, seed=4461 + s).random_base2(
                int(np.ceil(np.log2(n0)))
            )
            prior = {v: base[v].copy() for v in variants}
            results = {v: [] for v in variants}
            for batch in (0, 1):
                jobs = []
                if batch == 0:
                    grid = bounds[:, 0] + uniform[: n0 - 1] * (bounds[:, 1] - bounds[:, 0])
                else:
                    noise = rng.normal(0, 0.18, size=(n - n0, len(keys)))
                for d in variants:
                    center = base[d] if batch == 0 else max(results[d], key=priority)["parameters"]
                    if batch == 0:
                        candidates = [center.copy(), *[{**center, **dict(zip(keys, g))} for g in grid]]
                    else:
                        middle = np.array([center[key] for key in keys])
                        points = np.clip(
                            middle + noise * (bounds[:, 1] - bounds[:, 0]), bounds[:, 0], bounds[:, 1]
                        )
                        candidates = [{**center, **dict(zip(keys, g))} for g in points]
                    for j, p in enumerate(candidates):
                        jobs.append(
                            (d, s, p, args.dt, atlas, args.reserve, j + (0 if batch == 0 else n0), trials)
                        )
                for job in cf.as_completed([pool.submit(worker, j) for j in jobs]):
                    result = job.result()
                    d = result["design"]
                    results[d].append(result)
                    print(
                        s,
                        d,
                        result["trial"],
                        round(result["completed_s"], 2),
                        result["reason"],
                        round(result["score"], 3),
                        "accepted" if result["accepted"] else "",
                        flush=True,
                    )
            winners = {d: max(results[d], key=priority) for d in variants}
            for d, row in winners.items():
                base[d] = row["parameters"].copy()
            check_jobs = {}
            # Separate processes handle exact frozen repeats, independent initial
            # poses, and step-size sensitivity. No finalist is refit on these.
            with cf.ThreadPoolExecutor(max_workers=args.workers) as fresh:
                for d in variants:
                    for tag, p, dt, perturb in (
                        ("replay", base[d], args.dt, 0.0),
                        ("heldout", base[d], args.dt, 1.25),
                        ("halfstep", base[d], args.dt / 2, 0.0),
                        ("previous_frozen", prior[d], args.dt, 0.0),
                    ):
                        task = fresh.submit(
                            run_fresh, checks, s, d, p, tag, args.atlas, dt, args.reserve, perturb
                        )
                        check_jobs[task] = (d, tag)
                snapshots = {d: {} for d in variants}
                for task in cf.as_completed(check_jobs):
                    d, tag = check_jobs[task]
                    snapshots[d][tag] = task.result()
            for d in variants:
                valid = all(snapshots[d][name]["accepted"] for name in ("replay", "heldout", "halfstep"))
                stages[f"{d}_{s}"] = dict(
                    winner=winners[d]["trial"],
                    trial_count=len(results[d]),
                    completed=sum(r["complete"] for r in results[d]),
                    passed=sum(r["accepted"] for r in results[d]),
                    accepted_in_fresh_checks=valid,
                    accepted_prerequisites=passed[d],
                    parameters=base[d],
                    checks={
                        k: {
                            a: v[a]
                            for a in ("accepted", "reason", "completed_s", "score", "measured_support_BW")
                        }
                        for k, v in snapshots[d].items()
                    },
                )
                passed[d] &= valid
            dump(folder / "progress.json", stages)


if __name__ == "__main__":
    main()
