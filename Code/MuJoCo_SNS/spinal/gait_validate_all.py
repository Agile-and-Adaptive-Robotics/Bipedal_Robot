"""Literature validation, all walkers (2026-09-30): record a WALK run
for each variant (s3k, s3kpruned combo, w2lvar, syn6) through its
authoritative config, then score it against EVERY gait_refs reference
via kine_ref.compare. Writes the underlying trajectories and scores to the
selected output directory. The runner stores angles in degrees in its NPZ.

Usage: python gait_validate_all.py [variant ...] [--output DIRECTORY]
"""
import argparse
import csv
import io
import json
import os
import sys
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import numpy as np

import kine_ref as KR
import params as P
import runner as R
from gait_lib_score_all import load_npz_ref


def _jload(fn):
    with open(HERE / fn, encoding="utf-8") as f:
        return json.load(f)


BASE = dict(_jload("best_walk_params_v10.json")["multipliers"])
PRUNE_KEYS = {"heel_rge": 0.0, "toe_rge": 0.0, "contact_onset": 0.0,
              "contra_swing": 0.0, "pm_gain": 0.0, "pm_aff": 0.0,
              "ib_rge": 0.0, "ib_e_central": 0.0, "rg_weak_exc": 0.0}


def config_variant(v):
    os.environ.pop("AARL_NET", None)
    if v == "s3k":
        import _curriculum as CS3
        W = {**BASE, **_jload(
            "reports_20260923/s3k_trial34_full_params.json")["params"]}
        CS3.set_stage(5, dict(W))
        return repr(float(W["drive"])), []
    if v == "s3kpruned":
        import _curriculum as CS3
        W = {**BASE, **_jload(
            "reports_20260923/s3k_trial34_full_params.json")["params"]}
        CS3.set_stage(5, dict(W))
        for k, vv in PRUNE_KEYS.items():
            P.G[k] = float(vv)
        return repr(float(W["drive"])), ["--no-interleg"]
    if v == "w2lvar":
        os.environ["AARL_NET"] = "w2lvar"
        import _curriculum_w2lvar as CV
        W = {**BASE, "renshaw": 0.5, **_jload(
            "curriculum_w2lvar_stage5.json")["params"]}
        CV.set_stage(5, dict(W))
        return repr(float(W["drive"])), []
    if v == "syn6":
        os.environ["AARL_NET"] = "syn6"
        import _curriculum_syn6 as CS
        W = {**BASE, "renshaw": 0.0, "syn6": 1.0, "syn6_brainstem": 0.0,
             **_jload("curriculum_syn6_stage5.json")["params"]}
        CS.set_stage(5, dict(W))
        return repr(float(W["drive"])), []


def _ensure_stdout():
    """The curriculum loaders each wrap sys.stdout; after several R.main
    calls an inner wrapper gets GC'd and closes the stream ("I/O
    operation on closed file" - killed the syn6 pass 2026-09-30).
    Re-open fd 1 directly if the current stream is dead."""
    try:
        sys.stdout.write("")
    except Exception:
        sys.stdout = io.TextIOWrapper(open(1, "wb", closefd=False),
                                      encoding="utf-8",
                                      errors="replace")


def score_trace(trajectory, references):
    """Return one bilateral result per library record using the actual run interval."""
    rows = []
    with np.load(trajectory, allow_pickle=False) as trace:
        stop = int(np.count_nonzero(trace["t"] > 0)) + 1
        t, degrees, neuro, contact = [trace[key][:stop] for key in ("t", "q", "neuro", "contact")]
    for reference in references:
        if reference.stem == "ong_selfsel_Init200":
            rows.append((reference.stem, None, "degenerate left reference cycle"))
            continue
        result = KR.compare(t, degrees, neuro, 2.0, ref=load_npz_ref(reference), contact=contact)
        rows.append((reference.stem, result["kine_score"] if result else None,
                     "scored" if result else "no qualifying model cycles"))
    return rows


def main():
    names = ("s3k", "s3kpruned", "w2lvar", "syn6")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variants", nargs="*", metavar="VARIANT")
    parser.add_argument("--output", type=Path, default=HERE / "gait_validation_results")
    args = parser.parse_args()
    if set(args.variants).difference(names):
        parser.error("Known variants are " + ", ".join(names))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    references = sorted((HERE / "gait_refs").glob("*.npz"))
    path = output / "scores.csv"
    old_rows = []
    if path.exists():
        with path.open(encoding="utf-8", newline="") as stream:
            old_rows = list(csv.DictReader(stream))
    fields = ("variant", "reference", "kine_score", "status")
    for name in args.variants or names:
        completed = {row["reference"] for row in old_rows if row["variant"] == name}
        if completed == {"kine_score", *(ref.stem for ref in references)}:
            continue
        trajectory = output / f"{name}.npz"
        os.environ["AARL_NPZ"] = str(trajectory)
        drive, extra = config_variant(name)
        metrics = R.main(["--eval", "--drive", drive, *extra]) or {}
        _ensure_stdout()
        primary = metrics.get("kine_score")
        rows = [("kine_score", primary, "scored" if primary is not None else "no qualifying model cycles")]
        rows.extend(score_trace(trajectory, references))
        old_rows = [row for row in old_rows if row["variant"] != name]
        old_rows.extend(dict(zip(fields, (name, *row))) for row in rows)
        temporary = path.with_suffix(".tmp")
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(old_rows)
        temporary.replace(path)
        print(f"[{name}] wrote {len(references)} library records to {path}", flush=True)


if __name__ == "__main__":
    main()
