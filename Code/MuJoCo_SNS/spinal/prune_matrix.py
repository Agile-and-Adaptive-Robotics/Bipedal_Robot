"""Prune-matrix campaign (Ben 2026-09-28): dense-connectivity ablation.

Every variant's winner config is loaded through its OWN authoritative
curriculum set_stage (bit-faithful merge: v10 multipliers shell + stage
json), then ablation "cuts" are applied by direct params.G writes AFTER
set_stage and BEFORE R.main (a fresh net is built inside every R.main
call, so build-time conditional topology honors the cut values). A
params.G/TAU snapshot taken right after set_stage is restored before
every cell - no cut can leak into the next cell regardless of what each
set_stage happens to reset.

Configs per variant: FULL + 8 cuts (noaff / interleg / contact / ia /
ii / ib / renshaw / rgweak). Modes per config: AIR (walk cfg, no
ground), AIRDEAFF (air + --no-afferents --no-interleg), WALK (ground
--eval; for w2lvar-FULL this doubles as the t68 winner replay-verify),
STAND (balance cfg, --stand-eval 8), PUSH x4 (+x/-x/+y/-y 40 N pelvis
pulse at t=4 s via the runner --push flags).

Usage: python prune_matrix.py <s3k|w2lvar|syn6> <shard> <nshards>
One JSON line per finished cell -> prune_results_<variant>.jsonl
(crash-safe; re-running skips already-recorded cells).
"""
import io
import json
import os
import sys

# BLAS threading pinned to 1 BEFORE numpy imports: 9 concurrent shard
# processes on a 10-core box, and identical summation order inside the
# matrix (cross-shard comparability beats raw speed here).
_o = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "NUMEXPR_NUM_THREADS")
for _k in _o:
    os.environ[_k] = "1"

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

VARIANT = sys.argv[1].lower() if len(sys.argv) > 1 else "s3k"
SHARD = int(sys.argv[2]) if len(sys.argv) > 2 else 0
NSH = int(sys.argv[3]) if len(sys.argv) > 3 else 1

NPZ = f"spinal_run_prune_{VARIANT}_{SHARD}.npz"
os.environ["AARL_NPZ"] = NPZ
if VARIANT == "w2lvar":
    os.environ["AARL_NET"] = "w2lvar"
elif VARIANT == "syn6":
    os.environ["AARL_NET"] = "syn6"
else:
    os.environ.pop("AARL_NET", None)   # stock net

import time as _time

import numpy as np

import params as P
import runner as R

OUT = f"prune_results_{VARIANT}.jsonl"


def _jload(fn):
    with open(fn, encoding="utf-8") as f:
        return json.load(f)


_prev10 = _jload("best_walk_params_v10.json")
BASE = dict(_prev10["multipliers"])

if VARIANT == "w2lvar":
    import _curriculum_w2lvar as CV
    BASE["renshaw"] = 0.5
    WALK_P = {**BASE, **_jload("curriculum_w2lvar_stage5.json")["params"]}
    BAL_P = {**BASE, **_jload("curriculum_w2lvar_stage3.json")["params"]}

    def set_walk():
        CV.set_stage(5, dict(WALK_P))

    def set_bal():
        CV.set_stage(3, dict(BAL_P))

    CUTS = {
        "noaff": {"argv": ["--no-afferents"], "keys": {}},
        "interleg": {"argv": ["--no-interleg"], "keys": {}},
        "contact": {"argv": [],
                    "keys": {"contact_onset": 0.0, "contra_swing": 0.0}},
        "ia": {"argv": [],
               "keys": {"ia_to_mn": 0.0, "ia_to_antagonist": 0.0}},
        "ii": {"argv": [], "keys": {"ii_to_mn": 0.0}},
        "ib": {"argv": [],
               "keys": {"ib_to_mn_inh": 0.0, "ib_group_exc": 0.0,
                        "ib_exc_to_mn": 0.0}},
        "renshaw": {"argv": [], "keys": {"renshaw": 0.0}},
        "rgweak": {"argv": [], "keys": {"rg_weak_exc": 0.0}},
    }
elif VARIANT == "syn6":
    import _curriculum_syn6 as CS
    BASE["renshaw"] = 0.0
    BASE["syn6"] = 1.0
    BASE["syn6_brainstem"] = 0.0
    WALK_P = {**BASE, **_jload("curriculum_syn6_stage5.json")["params"]}
    BAL_P = {**BASE, **_jload("curriculum_syn6_stage3.json")["params"]}

    def set_walk():
        CS.set_stage(5, dict(WALK_P))

    def set_bal():
        CS.set_stage(3, dict(BAL_P))

    CUTS = {
        "noaff": {"argv": ["--no-afferents"], "keys": {}},
        "interleg": {"argv": ["--no-interleg"], "keys": {}},
        "contact": {"argv": [],
                    "keys": {"heel_rge": 0.0, "toe_rge": 0.0,
                             "contact_onset": 0.0, "contra_swing": 0.0}},
        "ia": {"argv": [],
               "keys": {"ia_to_mn": 0.0, "ia_to_antagonist": 0.0}},
        "ii": {"argv": [], "keys": {"ii_to_mn": 0.0}},
        "ib": {"argv": [], "keys": {"ib_rge": 0.0}},
        "renshaw": None,   # syn6 dress has NO Renshaw cells - nothing to cut
        "rgweak": {"argv": [], "keys": {"rg_weak_exc": 0.0}},
    }
else:   # s3k stock production walker
    import _curriculum as CS3
    BASE["renshaw"] = 0.5
    WALK_P = {**BASE, **_jload(
        "reports_20260923/s3k_trial34_full_params.json")["params"]}
    BAL_P = WALK_P   # s3k has no separate balance-tuned stage

    def set_walk():
        CS3.set_stage(5, dict(WALK_P))

    def set_bal():
        CS3.set_stage(5, dict(WALK_P))

    CUTS = {
        "noaff": {"argv": ["--no-afferents"], "keys": {}},
        "interleg": {"argv": ["--no-interleg"], "keys": {}},
        "contact": {"argv": [],
                    "keys": {"heel_rge": 0.0, "toe_rge": 0.0,
                             "contact_onset": 0.0, "contra_swing": 0.0,
                             "pm_gain": 0.0, "pm_aff": 0.0}},
        "ia": {"argv": [],
               "keys": {"ia_f_central": 0.0, "ia_in": 0.0,
                        "ia_f_contra_f": 0.0}},
        "ii": {"argv": [],
               "keys": {"ii_f_central": 0.0, "ii_e_central": 0.0}},
        "ib": {"argv": [], "keys": {"ib_rge": 0.0, "ib_e_central": 0.0}},
        "renshaw": {"argv": [], "keys": {"renshaw": 0.0}},
        "rgweak": {"argv": [], "keys": {"rg_weak_exc": 0.0}},
    }

PUSHES = [("+x", 40.0), ("-x", 40.0), ("+y", 40.0), ("-y", 40.0)]


def _cell_list():
    cells = []
    for cfg in ["full"] + [c for c in CUTS if CUTS[c] is not None]:
        cells += [(cfg, "AIR", ""), (cfg, "WALK", ""),
                  (cfg, "STAND", "")]
        cells += [(cfg, "PUSH", ax) for ax, _n in PUSHES]
        cells += [(cfg, "AIRDEAFF", "")]
    return cells


def _air_score():
    """Replicates the curriculum air objective verbatim (npz own-column
    identity, rises gate, deep-knee term)."""
    z = np.load(NPZ, allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    names = [str(x) for x in z["neuro_names"]]
    joints = [str(x) for x in z["key_joints"]]
    i_rge = names.index("RG_E_r")
    i_knee = joints.index("knee_angle_r")
    m = (t >= 5.0) & (t <= 17.0)
    if not np.all(np.isfinite(q[m])) or not np.all(np.isfinite(neuro[m])):
        return {"air_score": -200.0, "rises": None, "knee_min": None,
                "period": None}
    knee = q[m, i_knee]
    if not (-360.0 < float(knee.min()) < 360.0) or \
            not (-360.0 < float(knee.max()) < 360.0):
        return {"air_score": -200.0, "rises": None, "knee_min": None,
                "period": None}
    rge = neuro[m, i_rge]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    if rises > 30:
        return {"air_score": -200.0, "rises": rises, "knee_min": None,
                "period": None}
    if rises < 3 or (float(rge.max()) - float(rge.min())) < 1.0:
        score = -10.0 + 0.05 * (-float(knee.min()))
    else:
        score = 3.0 * rises + 0.5 * (-float(knee.min()))
    # burst period for the report (mean onset-to-onset over the window)
    idx = np.flatnonzero(np.diff(on.astype(int)) == 1)
    period = float(np.mean(np.diff(t[m][idx]))) if idx.size > 1 else None
    return {"air_score": float(score), "rises": rises,
            "knee_min": float(knee.min()), "period": period}


def _metrics_subset(m):
    m = m or {}
    out = {}
    for k, v in m.items():
        if isinstance(v, (int, float, bool)):
            out[k] = v
    kin = m.get("kine") or {}
    out["kine_components"] = {kk: vv for kk, vv in kin.items()
                              if isinstance(vv, (int, float, bool))}
    return out


def _append(rec):
    with open(OUT, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")


def _done_keys():
    done = set()
    if os.path.exists(OUT):
        with open(OUT, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                    done.add((r["config"], r["mode"], r.get("push", "")))
                except Exception:
                    pass
    return done


def main():
    cells = _cell_list()
    mine = cells[SHARD::NSH]
    done = _done_keys()
    print(f"[prune:{VARIANT}] shard {SHARD}/{NSH}: {len(mine)} cells, "
          f"{sum(1 for c in mine if c in done)} already done", flush=True)
    snapG, snapTAU = None, None
    cur_cfg = None
    for (cfg, mode, push) in mine:
        if (cfg, mode, push) in done:
            continue
        t0 = _time.time()
        try:
            if cfg != cur_cfg:
                # (re)build the clean config snapshot
                if mode in ("STAND", "PUSH"):
                    set_bal()
                else:
                    set_walk()
                snapG, snapTAU = dict(P.G), dict(P.TAU)
                cur_cfg = cfg
            P.G.clear()
            P.G.update(snapG)
            P.TAU.clear()
            P.TAU.update(snapTAU)
            argv_extra = []
            if cfg != "full":
                cut = CUTS[cfg]
                argv_extra += list(cut["argv"])
                for k, v in cut["keys"].items():
                    P.G[k] = float(v)
            drive = repr(float(WALK_P["drive"]))
            rig = repr(float(BAL_P.get("rig_scale", 1.0)))
            if mode == "AIR":
                args = ["--no-ground", "--time", "14", "--drive", drive]
            elif mode == "AIRDEAFF":
                args = ["--no-ground", "--no-afferents", "--no-interleg",
                        "--time", "14", "--drive", drive]
            elif mode == "WALK":
                args = ["--eval", "--drive", drive]
            elif mode == "STAND":
                args = ["--stand-eval", "8", "--rig-scale", rig]
            else:   # PUSH
                args = ["--stand-eval", "8", "--rig-scale", rig,
                        "--push", "40", "--push-time", "4.0",
                        "--push-axis", push]
            args = args + argv_extra
            m = R.main(args)
            rec = {"variant": VARIANT, "config": cfg, "mode": mode,
                   "push": push, "argv_extra": argv_extra,
                   "metrics": _metrics_subset(m)}
            if mode in ("AIR", "AIRDEAFF"):
                rec["air"] = _air_score()
            rec["wall_s"] = round(_time.time() - t0, 1)
            _append(rec)
            tag = rec["air"]["air_score"] if "air" in rec else \
                rec["metrics"].get("kine_score")
            print(f"[prune:{VARIANT}] {cfg:9s} {mode:8s} {push:2s} -> "
                  f"{tag}  ({rec['wall_s']}s)", flush=True)
        except Exception as ex:
            _append({"variant": VARIANT, "config": cfg, "mode": mode,
                     "push": push, "error": repr(ex),
                     "tb": traceback.format_exc().splitlines()[-1],
                     "wall_s": round(_time.time() - t0, 1)})
            print(f"[prune:{VARIANT}] {cfg} {mode} {push} ERROR {ex!r}",
                  flush=True)
    print(f"[prune:{VARIANT}] shard {SHARD} done", flush=True)


if __name__ == "__main__":
    main()
