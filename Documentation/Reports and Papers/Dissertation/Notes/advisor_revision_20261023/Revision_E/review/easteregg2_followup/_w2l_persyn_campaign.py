r"""W2L PER-SYNAPSE campaign (2026-10-25, egg; Ben's STOP-AND-RESCOPE ruling):
Optuna solves INDIVIDUAL conductances per synapse — no family scalars.

Search space = 59 individual strengths on bilateralrg2:
  - 24 per-synapse POOL keys (new): pf_mn_{6 pools}, rc_own_{6},
    rc_inh_{6}, rc_mut_{6} — L/R shared by symmetry (Ben's reduction).
  - 29 per-target afferent/contact keys (2026-10-25): ia/ii/ib per joint
    per target, Ia chain per joint, heel cascade per target.
  - ch_ank, ch_contra, comm + sensor encoding (ctn_amp, cref, ibnA).
Seeds: Ben's literal instance values (the 2026-10-25 re-mine: nearly all
5e-07 S -> the stack anchors express them) + Shinohara Table B.3 per-MN
pattern as literature prior on the pf_mn_* ratios. Signs/topology fixed.
Drive PINNED te=3.0/tf=4.0 (Ben's aproj defaults).

Harness: sagittal=50 + weight channel; WEAN dim = wsup (0.70 stage 1 ->
0.55/0.40/0.25/0.10 ladder, 26 trials/rung, honest stop on death).
Metric + survival gate VERBATIM (_w2l_wean_eval). Fresh-process verify per
rung. Basin/sensitivity/remove-and-retest run ONLY on the final winner
(validation, never a reason to shrink the space).

Study: curr_w2l_persyn + curr_w2l_persyn_w{NN} in optuna_w2l_persyn.db.
Usage: D:\\Anaconda\\envs\\myo\\python.exe _w2l_persyn_campaign.py [--n=120]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WEAN = os.path.normpath(os.path.join(
    HERE, "..", "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import _w2l_wean_eval as EV       # noqa: E402
import optuna  # noqa: E402  # (campaign driver import-time use only)

optuna.logging.set_verbosity(optuna.logging.WARNING)
EV.T.KV["net"] = "bilateralrg2"

DB = "sqlite:///" + os.path.join(HERE, "optuna_w2l_persyn.db").replace("\\", "/")
SEED = 20261025
DUR_TUNE, DUR_VERIFY = 20.0, 30.0
PY = sys.executable
TE_PIN, TF_PIN = 3.0, 4.0
SAG = 50.0
LADDER = [0.55, 0.40, 0.25, 0.10]

POOLS = ("hipext", "hipflx", "kneeext", "kneeflx", "ankext", "ankflx")
SPACE = dict(
    # 24 per-synapse pool conductances (anchors = former family values)
    **{f"pf_mn_{p}": (0.3, 6.0) for p in POOLS},
    **{f"rc_own_{p}": (0.15, 3.0) for p in POOLS},
    **{f"rc_inh_{p}": (0.08, 2.0) for p in POOLS},
    **{f"rc_mut_{p}": (0.08, 2.0) for p in POOLS},
    # 29 per-target afferent/contact classes
    heel_rg=(0.01, 1.0), heel_hip_pf=(0.01, 1.0), heel_knee_pf=(0.01, 1.0),
    heel_hip_mn=(0.01, 1.0), heel_knee_mn=(0.01, 1.0),
    ia_hip_rg=(0.01, 0.5), ia_hip_pf=(0.01, 0.5),
    ia_knee_rg=(0.01, 0.5), ia_knee_pf=(0.01, 0.5),
    ia_ank_rg=(0.01, 0.5), ia_ank_pf=(0.01, 0.5),
    ii_hip_rg=(0.01, 0.5), ii_hip_pf=(0.01, 0.5),
    ii_knee_rg=(0.01, 0.5), ii_knee_pf=(0.01, 0.5),
    ib_hip_rg=(0.01, 0.5), ib_hip_pf=(0.01, 0.5),
    ib_knee_rg=(0.01, 0.5), ib_knee_pf=(0.01, 0.5),
    ia_hip_sn=(0.1, 2.0), ia_hip_mn=(0.05, 1.5), ia_hip_mutual=(0.05, 1.5),
    ia_knee_sn=(0.1, 2.0), ia_knee_mn=(0.05, 1.5), ia_knee_mutual=(0.05, 1.5),
    ia_ank_sn=(0.1, 2.0), ia_ank_mn=(0.05, 1.5), ia_ank_mutual=(0.05, 1.5),
    ib_hip_mn=(0.05, 1.5), ib_knee_mn=(0.05, 1.5), ib_ank_mn=(0.05, 1.5),
    # cross-joint contact + commissural + sensor encoding
    ch_ank=(0.1, 2.0), ch_contra=(0.01, 2.0), comm=(0.2, 2.0),
    ctn_amp=(0.5, 6.0), cref=(20.0, 150.0), ibnA=(0.0, 3.0),
)
ANCHORS = dict(
    **{f"pf_mn_{p}": 2.0 for p in POOLS},
    **{f"rc_own_{p}": 1.0 for p in POOLS},
    **{f"rc_inh_{p}": 0.5 for p in POOLS},
    **{f"rc_mut_{p}": 0.5 for p in POOLS},
    heel_rg=0.08, heel_hip_pf=0.08, heel_knee_pf=0.08,
    heel_hip_mn=0.08, heel_knee_mn=0.08,
    **{k: 0.08 for k in SPACE if k.startswith(("ia_hip_r", "ia_knee_r",
        "ia_ank_r", "ia_hip_p", "ia_knee_p", "ia_ank_p", "ii_",
        "ib_hip_r", "ib_knee_r", "ib_hip_p", "ib_knee_p"))},
    ia_hip_sn=0.6, ia_knee_sn=0.6, ia_ank_sn=0.6,
    ia_hip_mn=0.4, ia_knee_mn=0.4, ia_ank_mn=0.4,
    ia_hip_mutual=0.5, ia_knee_mutual=0.5, ia_ank_mutual=0.5,
    ib_hip_mn=0.5, ib_knee_mn=0.5, ib_ank_mn=0.5,
    ch_ank=1.0, ch_contra=0.08, comm=1.0,
    ctn_amp=4.0, cref=50.0, ibnA=1.0,
)
SEEDS = [ANCHORS]   # seed 1 = Ben's literal values, every key individual


def knobs_of(wsup: float, p: dict, dur: float) -> dict:
    return dict(EV.BASE,
                dscale=1.0, te=TE_PIN, tf=TF_PIN,
                wsup=wsup, sagittal=SAG,
                **{k: float(p[k]) for k in SPACE},
                dur_walk=dur)


def fresh_verify(wsup: float, params: dict) -> dict:
    kn = knobs_of(wsup, params, DUR_VERIFY)
    tmp = os.path.join(HERE, "_finalist_knobs.json")
    json.dump(kn, open(tmp, "w"), indent=1)
    r = subprocess.run(
        [PY, os.path.join(HERE, "_verify_finalist.py"), tmp],
        capture_output=True, text=True, timeout=1800)
    for line in r.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    raise RuntimeError(f"fresh verify failed:\n{r.stdout[-800:]}\n{r.stderr[-800:]}")


def run_round(wsup: float, seed_params: dict | None, ntrials: int) -> dict:
    name = "curr_w2l_persyn" if abs(wsup - 0.70) < 1e-9 else \
        f"curr_w2l_persyn_w{int(round(wsup * 100)):02d}"
    study = optuna.create_study(
        study_name=name, storage=DB, direction="maximize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=SEED))
    have = len(study.trials)
    if have == 0:
        for sp in ([seed_params] if seed_params else []) + SEEDS:
            study.enqueue_trial(sp)
    need = max(ntrials - have, 0)
    if need == 0 and have == 0:
        need = 1

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi)
             for k, (lo, hi) in SPACE.items()}
        m = EV.evaluate(knobs_of(wsup, p, DUR_TUNE),
                        tag=f"{name}#{trial.number}",
                        csv_path=os.path.join(HERE, "persyn_trials.csv"))
        for k, v in m.items():
            if not k.startswith("_") and v is not None:
                try:
                    trial.set_user_attr(k, v)
                except (TypeError, ValueError):
                    pass
        return float(m["obj"])

    if need > 0:
        study.optimize(obj, n_trials=need)
    best = study.best_trial
    out = dict(wsup=wsup, study=name, n_trials=len(study.trials),
               params={k: round(float(best.params[k]), 5) for k in SPACE},
               obj20=float(best.value))
    m30 = fresh_verify(wsup, best.params)
    out["m30_fresh"] = {k: m30.get(k) for k in EV.METRIC_KEYS}
    out["survive30"] = bool(m30.get("survive"))
    path = os.path.join(HERE, f"round_w{int(round(wsup * 100)):02d}.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print(f"== ROUND wsup={wsup:.2f} best obj20={out['obj20']:.1f} "
          f"survive30(fresh)={out['survive30']} -> {path}", flush=True)
    return out


def main() -> int:
    n = 120
    for a in sys.argv[1:]:
        if a.startswith("--n="):
            n = int(a[4:])
    r70 = run_round(0.70, None, n)
    if not r70["survive30"]:
        print("stage-1 (wsup 0.70) did not survive fresh verify; honest stop",
              flush=True)
        return 0
    prev = r70["params"]
    for ws in LADDER:
        r = run_round(ws, prev, 26)
        if not r["survive30"]:
            print(f"WEAN DEATH at wsup={ws:.2f} (honest stop; reported)",
                  flush=True)
            return 0
        prev = r["params"]
    print("WEAN COMPLETE through wsup=0.10 — winner validation next "
          "(basin + sensitivity + remove-and-retest)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
