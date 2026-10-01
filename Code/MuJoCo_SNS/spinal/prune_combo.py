"""Round-2 prune confirmation (2026-09-28): s3k with ALL FOUR prunable
components cut together (interleg + contact + ib + rgweak - the four
PRUNE verdicts from the leave-one-out matrix). Cells: AIR, WALK, STAND,
PUSH x4, appended to prune_results_s3k.jsonl as config "combo".
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

os.environ["AARL_NPZ"] = "spinal_run_prune_s3k_combo.npz"
os.environ.pop("AARL_NET", None)

import time as _time

import numpy as np

import params as P
import runner as R
import _curriculum as CS3

OUT = "prune_results_s3k.jsonl"


def _jload(fn):
    with open(fn, encoding="utf-8") as f:
        return json.load(f)


BASE = dict(_jload("best_walk_params_v10.json")["multipliers"])
BASE["renshaw"] = 0.5
WALK_P = {**BASE, **_jload(
    "reports_20260923/s3k_trial34_full_params.json")["params"]}

COMBO_KEYS = {"heel_rge": 0.0, "toe_rge": 0.0, "contact_onset": 0.0,
              "contra_swing": 0.0, "pm_gain": 0.0, "pm_aff": 0.0,
              "ib_rge": 0.0, "ib_e_central": 0.0, "rg_weak_exc": 0.0}
COMBO_ARGV = ["--no-interleg"]

PUSHES = [("+x",), ("-x",), ("+y",), ("-y",)]


def _air_score():
    z = np.load("spinal_run_prune_s3k_combo.npz", allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    names = [str(x) for x in z["neuro_names"]]
    joints = [str(x) for x in z["key_joints"]]
    i_rge = names.index("RG_E_r")
    i_knee = joints.index("knee_angle_r")
    m = (t >= 5.0) & (t <= 17.0)
    if not np.all(np.isfinite(q[m])) or not np.all(np.isfinite(neuro[m])):
        return {"air_score": -200.0, "rises": None, "knee_min": None}
    knee = q[m, i_knee]
    rge = neuro[m, i_rge]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    if rises > 30:
        return {"air_score": -200.0, "rises": rises, "knee_min": None}
    if rises < 3 or (float(rge.max()) - float(rge.min())) < 1.0:
        score = -10.0 + 0.05 * (-float(knee.min()))
    else:
        score = 3.0 * rises + 0.5 * (-float(knee.min()))
    return {"air_score": float(score), "rises": rises,
            "knee_min": float(knee.min())}


def _append(rec):
    with open(OUT, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")


def main():
    CS3.set_stage(5, dict(WALK_P))
    snapG, snapTAU = dict(P.G), dict(P.TAU)
    drive = repr(float(WALK_P["drive"]))
    cells = [("AIR", []), ("WALK", []), ("STAND", [])] + \
            [("PUSH", ax) for ax, in PUSHES]
    for mode, ax in cells:
        t0 = _time.time()
        try:
            P.G.clear()
            P.G.update(snapG)
            P.TAU.clear()
            P.TAU.update(snapTAU)
            for k, v in COMBO_KEYS.items():
                P.G[k] = float(v)
            if mode == "AIR":
                args = ["--no-ground", "--time", "14", "--drive", drive]
            elif mode == "WALK":
                args = ["--eval", "--drive", drive]
            elif mode == "STAND":
                args = ["--stand-eval", "8", "--rig-scale", "1.0"]
            else:
                args = ["--stand-eval", "8", "--rig-scale", "1.0",
                        "--push", "40", "--push-time", "4.0",
                        "--push-axis", ax]
            args = args + COMBO_ARGV
            m = R.main(args) or {}
            met = {k: v for k, v in m.items()
                   if isinstance(v, (int, float, bool))}
            kin = m.get("kine") or {}
            met["kine_components"] = {kk: vv for kk, vv in kin.items()
                                      if isinstance(vv, (int, float, bool))}
            rec = {"variant": "s3k", "config": "combo", "mode": mode,
                   "push": ax if ax else "", "argv_extra": COMBO_ARGV,
                   "metrics": met}
            if mode == "AIR":
                rec["air"] = _air_score()
            rec["wall_s"] = round(_time.time() - t0, 1)
            _append(rec)
            print(f"[combo] {mode} {ax} done ({rec['wall_s']}s)",
                  flush=True)
        except Exception as ex:
            _append({"variant": "s3k", "config": "combo", "mode": mode,
                     "push": ax if ax else "", "error": repr(ex),
                     "wall_s": round(_time.time() - t0, 1)})
            print(f"[combo] {mode} {ax} ERROR {ex!r}", flush=True)
    print("[combo] done", flush=True)


if __name__ == "__main__":
    main()
