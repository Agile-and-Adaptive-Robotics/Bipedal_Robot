"""Kinematic-robustness sweep (Ben 2026-09-30): slightly change the
kinematics - pelvis start height +/-1, +/-2 cm (AARL_PELVIS_TY) - and
re-test WALK / STAND / PUSH. 4 walkers x 5 height levels x 4 modes.

Traps honored: stock set_stage OVERWRITES AARL_PELVIS_TY from p
["pelvis_ty"] (s3k family) -> for those the perturbation goes through
the param dict; w2lvar/syn6 set_stage POPS the env -> set it AFTER
set_stage. One jsonl line per cell -> robust_results_<variant>.jsonl.

Usage: python prune_robust.py <s3k|s3kpruned|w2lvar|syn6> <shard> <nshards>
"""
import io
import json
import os
import sys

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

VARIANT = sys.argv[1].lower()
SHARD = int(sys.argv[2]) if len(sys.argv) > 2 else 0
NSH = int(sys.argv[3]) if len(sys.argv) > 3 else 1

os.environ["AARL_NPZ"] = f"spinal_run_robust_{VARIANT}_{os.getpid()}.npz"
if VARIANT in ("w2lvar", "syn6"):
    os.environ["AARL_NET"] = VARIANT
else:
    os.environ.pop("AARL_NET", None)

import time as _time

import params as P
import runner as R


def _jload(fn):
    with open(fn, encoding="utf-8") as f:
        return json.load(f)


BASE = dict(_jload("best_walk_params_v10.json")["multipliers"])

if VARIANT == "w2lvar":
    import _curriculum_w2lvar as CV
    BASE["renshaw"] = 0.5
    WALK_P = {**BASE, **_jload("curriculum_w2lvar_stage5.json")["params"]}
    BAL_P = {**BASE, **_jload("curriculum_w2lvar_stage3.json")["params"]}
    TY_BASE = None      # no pelvis_ty param; env knob after set_stage

    def set_walk(ty):
        CV.set_stage(5, dict(WALK_P))
        if ty is not None:
            os.environ["AARL_PELVIS_TY"] = repr(ty)

    def set_bal(ty):
        CV.set_stage(3, dict(BAL_P))
        if ty is not None:
            os.environ["AARL_PELVIS_TY"] = repr(ty)
elif VARIANT == "syn6":
    import _curriculum_syn6 as CS
    BASE["renshaw"] = 0.0
    BASE["syn6"] = 1.0
    BASE["syn6_brainstem"] = 0.0
    WALK_P = {**BASE, **_jload("curriculum_syn6_stage5.json")["params"]}
    BAL_P = {**BASE, **_jload("curriculum_syn6_stage3.json")["params"]}
    TY_BASE = None      # like w2lvar: perturb via env after set_stage

    def set_walk(ty):
        CS.set_stage(5, dict(WALK_P))
        if ty is not None:
            os.environ["AARL_PELVIS_TY"] = repr(ty)

    def set_bal(ty):
        CS.set_stage(3, dict(BAL_P))
        if ty is not None:
            os.environ["AARL_PELVIS_TY"] = repr(ty)
else:
    # s3k or s3kpruned: pelvis_ty flows through the param dict
    import _curriculum as CS3
    BASE["renshaw"] = 0.5
    WALK_P = {**BASE, **_jload(
        "reports_20260923/s3k_trial34_full_params.json")["params"]}
    BAL_P = WALK_P
    TY_BASE = float(WALK_P["pelvis_ty"])
    PRUNE_KEYS = {"heel_rge": 0.0, "toe_rge": 0.0, "contact_onset": 0.0,
                  "contra_swing": 0.0, "pm_gain": 0.0, "pm_aff": 0.0,
                  "ib_rge": 0.0, "ib_e_central": 0.0,
                  "rg_weak_exc": 0.0}

    def _set(ty, stage_p):
        p = dict(stage_p)
        if ty is not None:
            p["pelvis_ty"] = float(ty)
        CS3.set_stage(5, p)
        if VARIANT == "s3kpruned":
            for k, v in PRUNE_KEYS.items():
                P.G[k] = float(v)


    def set_walk(ty):
        _set(ty, WALK_P)

    def set_bal(ty):
        _set(ty, BAL_P)

TYS = [None, 0.01, -0.01, 0.02, -0.02]     # None = nominal
MODES = [("WALK", ""), ("STAND", ""), ("PUSH", "+x"), ("PUSH", "+y")]
ARGV = ["--no-interleg"] if VARIANT == "s3kpruned" else []
OUT = f"robust_results_{VARIANT}.jsonl"


def _cells():
    out = []
    for d in TYS:
        ty = (TY_BASE + d) if (TY_BASE is not None and d is not None) \
            else (TY_BASE if d is None else d)
        for mode, ax in MODES:
            out.append((d, ty, mode, ax))
    return out


def _append(rec):
    with open(OUT, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")


def main():
    cells = _cells()
    mine = cells[SHARD::NSH]
    done = set()
    if os.path.exists(OUT):
        with open(OUT, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                    done.add((r.get("dty"), r["mode"], r.get("push", "")))
                except Exception:
                    pass
    print(f"[robust:{VARIANT}] shard {SHARD}/{NSH}: {len(mine)} cells, "
          f"{sum(1 for c in mine if (c[0], c[2], c[3]) in done)} done",
          flush=True)
    for (d, ty, mode, ax) in mine:
        if (d, mode, ax) in done:
            continue
        t0 = _time.time()
        try:
            os.environ.pop("AARL_PELVIS_TY", None)
            if mode in ("STAND", "PUSH"):
                set_bal(ty)
            else:
                set_walk(ty)
            drive = repr(float(WALK_P["drive"]))
            # CAREFUL: elif chain required - an earlier edit made the
            # STAND branch a separate `if`, so WALK cells fell through
            # to the PUSH else and every WALK row ran as a push eval
            # (caught 2026-09-30: WALK kine=-25 sentinel + nonzero
            # push-sway gave it away).
            if mode == "WALK":
                args = ["--eval", "--drive", drive] + ARGV
            elif mode == "STAND":
                args = ["--stand-eval", "8", "--rig-scale", "1.0"]
            else:
                args = ["--stand-eval", "8", "--rig-scale", "1.0",
                        "--push", "40", "--push-time", "4.0",
                        "--push-axis", ax] + ARGV
            m = R.main(args) or {}
            met = {k: v for k, v in m.items()
                   if isinstance(v, (int, float, bool))}
            kin = m.get("kine") or {}
            met["kine_components"] = {kk: vv for kk, vv in kin.items()
                                      if isinstance(vv, (int, float, bool))}
            _append({"variant": VARIANT, "dty": d, "mode": mode,
                     "push": ax, "metrics": met,
                     "wall_s": round(_time.time() - t0, 1)})
            print(f"[robust:{VARIANT}] dty={d} {mode} {ax} -> "
                  f"kine={met.get('kine_score')} fell={met.get('bal_fell')}",
                  flush=True)
        except Exception as ex:
            _append({"variant": VARIANT, "dty": d, "mode": mode,
                     "push": ax, "error": repr(ex),
                     "wall_s": round(_time.time() - t0, 1)})
            print(f"[robust:{VARIANT}] dty={d} {mode} {ax} ERR {ex!r}",
                  flush=True)
    print(f"[robust:{VARIANT}] shard {SHARD} done", flush=True)


if __name__ == "__main__":
    main()
