"""Generalized post-run probe for a curriculum study winner (parameterizes
the _s3g/_s3j/_s3k_probe.py family). Loads the study's best trial (or a given
trial), merges curriculum_stage3.json + v10 BASE_MUL, runs the in-process
eval, and prints the full per-leg gait metric set from kine_ref.compare.

MUST run under the SNS env (walker.cmd does that):
  walker.cmd probe_study.py --study curr_s3k_nocross
  walker.cmd probe_study.py --study curr_s3j_fullrules --trial 12
  walker.cmd probe_study.py --study NAME --no-stage3 --set renshaw=0.0
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402

METRIC_KEYS = (
    "kine_score", "duty_r", "duty_l", "n_cycles_r", "n_cycles_l",
    "frozen_r", "frozen_l", "ds", "contact_frac_r", "contact_frac_l",
    "knee_min_r", "knee_min_l", "mean_hip_r", "mean_hip_l",
    "T_r", "T_l", "period_cv_r", "period_cv_l", "lag_rl",
)


def main() -> int:
    wcommon.utf8_stdio()
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--trial", type=int, default=None,
                    help="trial number (default: best by value)")
    ap.add_argument("--npz", default=None, help="npz name (default auto)")
    ap.add_argument("--out", default=None, help="output json (default <study>_probe_out.json)")
    ap.add_argument("--no-stage3", action="store_true",
                    help="skip merging curriculum_stage3.json params")
    ap.add_argument("--set", action="append", default=[],
                    help="extra param override key=float (repeatable)")
    ap.add_argument("--renshaw", type=float, default=0.5,
                    help="BASE_MUL renshaw (probes used 0.5; pass 0 to disable)")
    args = ap.parse_args()

    spinal = wcommon.resolve_spinal()
    os.chdir(spinal)
    sys.path.insert(0, str(spinal))
    wcommon.banner(f"probe study={args.study} spinal={spinal}")

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    import _curriculum as C
    import kine_ref as KR
    import numpy as np
    import runner as R

    # --- pick the trial -----------------------------------------------------
    study = optuna.load_study(study_name=args.study, storage=wcommon.db_url(spinal))
    if args.trial is not None:
        trial = next((t for t in study.trials if t.number == args.trial), None)
        if trial is None or trial.value is None:
            raise SystemExit(f"trial {args.trial} not found / has no value in {args.study}")
        src = f"given trial {trial.number}"
    else:
        trial = max((t for t in study.trials if t.value is not None),
                    key=lambda t: t.value)
        src = f"best trial {trial.number}"
    print(f"trial  : {src}  value={trial.value:.3f}")

    # --- parameter merge (s3k-style) ----------------------------------------
    mul = json.loads((spinal / "best_walk_params_v10.json").read_text())["multipliers"]
    C.BASE_MUL = dict(mul)          # module global is None until main() runs
    C.BASE_MUL["renshaw"] = args.renshaw
    params = dict(mul)
    if not args.no_stage3 and (spinal / "curriculum_stage3.json").exists():
        s3 = json.loads((spinal / "curriculum_stage3.json").read_text())
        params.update(s3.get("params", {}))
    params.update(trial.params)
    for kv in args.set:
        k, _, v = kv.partition("=")
        params[k.strip()] = float(v)

    # --- run the eval --------------------------------------------------------
    npz_name = args.npz or f"spinal_run_probe_{args.study}.npz"
    os.environ["AARL_NPZ"] = npz_name
    C.set_stage(3, params)
    t0 = time.time()
    m = R.main(["--eval", "--drive", repr(float(params["drive"]))])
    print(f"eval   : finished in {time.time() - t0:.0f}s  nan={m.get('nan')}")

    # --- gait metrics --------------------------------------------------------
    with np.load(npz_name) as z:
        k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                       ref=KR.ref_cached(), contact=z["contact"])
    out = {"study": args.study, "trial": trial.number,
           "value": trial.value, "drive": float(params["drive"])}
    out.update({kk: (round(k[kk], 3) if isinstance(k.get(kk), float) else k.get(kk))
                for kk in METRIC_KEYS})
    out["eval"] = {kk: m.get(kk) for kk in
                   ("dx", "kz", "tilt_max", "knee_min", "hip_amp", "duty", "burst_r")
                   if isinstance(m.get(kk), (int, float))}

    out_file = Path(args.out) if args.out else spinal / f"{args.study}_probe_out.json"
    out_file.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("--- metrics ---")
    for kk in METRIC_KEYS:
        print(f"  {kk:16s} {out.get(kk)}")
    print(f"wrote {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
