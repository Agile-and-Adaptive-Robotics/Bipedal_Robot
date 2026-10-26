r"""W2L PER-TARGET-CLASS campaign (2026-10-25, egg; Ben's directive after the
class-collapse ruling): "Knee and Ankle Ia feedback is a different class for
each target it reaches, same with Hip Ia, etc. etc. Run it and optimize it
so that it walks with a simulated virtual harness that takes up weight and
constrains to sagittal plane. Run iterations so that the harness takes up
less and less weight."

Net: --net=bilateralrg2 (Ben's corrected aproj). The builder now carries 29
PER-TARGET classes (build_w2l_net W2L_GAINS: ia/ii/ib per joint per target,
Ia chain per joint, Ib MN per joint, heel cascade per target) — NO family
collapsing; defaults = historical anchors = bit-identical (verified
2026-10-25: census 96/204 unchanged, override changes only its own key).

Harness: sagittal=50 (lateral leash stiffened 50x, roll+yaw locked, pitch
and fore-aft free; 0 = historical rig) + the drop-protocol weight channel.
WEAN DIMENSION = wsup (harness weight fraction), NOT drive: stage 1 tunes
at wsup=0.70; the ladder then re-optimizes at 0.55/0.40/0.25/0.10, 26
trials/rung seeded from the previous rung's winner, stopping honestly on
death. Drive PINNED te=3.0/tf=4.0 (Ben's literal aproj defaults). ch_hc is
INERT now (its edges resolve to the heel_* classes) — pinned at default.

Metric + survival gate VERBATIM (_w2l_wean_eval, unchanged). Fresh-process
verify per rung; strict basin on the final survivor (12 trials, +/-1%,
bar 11/12). Anti-lottery: per-term diagnostics per trial + csv.

Study: curr_w2l_ptclass (stage 1) + curr_w2l_ptclass_w{NN} rungs in this
folder's optuna_w2l_ptclass.db. Usage:
    D:\\Anaconda\\envs\\myo\\python.exe _w2l_ptclass_campaign.py [--n=100]
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
import optuna                     # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)
EV.T.KV["net"] = "bilateralrg2"

DB = "sqlite:///" + os.path.join(HERE, "optuna_w2l_ptclass.db").replace("\\", "/")
SEED = 20261025
DUR_TUNE, DUR_VERIFY = 20.0, 30.0
PY = sys.executable
TE_PIN, TF_PIN = 3.0, 4.0
SAG = 50.0
LADDER = [0.55, 0.40, 0.25, 0.10]

PT = dict(
    # heel cascade per target (anchor 0.08)
    heel_rg=(0.01, 1.0), heel_hip_pf=(0.01, 1.0), heel_knee_pf=(0.01, 1.0),
    heel_hip_mn=(0.01, 1.0), heel_knee_mn=(0.01, 1.0),
    # Ia -> RG / PF per joint (anchor 0.08)
    ia_hip_rg=(0.01, 0.5), ia_hip_pf=(0.01, 0.5),
    ia_knee_rg=(0.01, 0.5), ia_knee_pf=(0.01, 0.5),
    ia_ank_rg=(0.01, 0.5), ia_ank_pf=(0.01, 0.5),
    # II -> RG / PF per joint (anchor 0.08)
    ii_hip_rg=(0.01, 0.5), ii_hip_pf=(0.01, 0.5),
    ii_knee_rg=(0.01, 0.5), ii_knee_pf=(0.01, 0.5),
    # Ib -> RG / PF per joint (anchor 0.08)
    ib_hip_rg=(0.01, 0.5), ib_hip_pf=(0.01, 0.5),
    ib_knee_rg=(0.01, 0.5), ib_knee_pf=(0.01, 0.5),
    # Ia chain per joint (sn 0.6 / mn 0.4 / mutual 0.5)
    ia_hip_sn=(0.1, 2.0), ia_hip_mn=(0.05, 1.5), ia_hip_mutual=(0.05, 1.5),
    ia_knee_sn=(0.1, 2.0), ia_knee_mn=(0.05, 1.5), ia_knee_mutual=(0.05, 1.5),
    ia_ank_sn=(0.1, 2.0), ia_ank_mn=(0.05, 1.5), ia_ank_mutual=(0.05, 1.5),
    # Ib autogenic MN per joint (anchor 0.5)
    ib_hip_mn=(0.05, 1.5), ib_knee_mn=(0.05, 1.5), ib_ank_mn=(0.05, 1.5),
)
SPACE = dict(
    PT,
    ch_ank=(0.1, 2.0), ch_contra=(0.01, 2.0), comm=(0.2, 2.0),
    ctn_amp=(0.5, 6.0), cref=(20.0, 150.0), ibnA=(0.0, 3.0),
)
ANCHORS = dict(
    heel_rg=0.08, heel_hip_pf=0.08, heel_knee_pf=0.08,
    heel_hip_mn=0.08, heel_knee_mn=0.08,
    **{k: 0.08 for k in SPACE if k.startswith(("ia_hip_r", "ia_knee_r",
        "ia_ank_r", "ia_hip_p", "ia_knee_p", "ia_ank_p",
        "ii_", "ib_hip_r", "ib_knee_r", "ib_hip_p", "ib_knee_p"))},
    ia_hip_sn=0.6, ia_knee_sn=0.6, ia_ank_sn=0.6,
    ia_hip_mn=0.4, ia_knee_mn=0.4, ia_ank_mn=0.4,
    ia_hip_mutual=0.5, ia_knee_mutual=0.5, ia_ank_mutual=0.5,
    ib_hip_mn=0.5, ib_knee_mn=0.5, ib_ank_mn=0.5,
    ch_ank=1.0, ch_contra=0.08, comm=1.0,
    ctn_amp=4.0, cref=50.0, ibnA=1.0,
)
SEEDS = [ANCHORS]   # seed 1 = Ben's literal per-target defaults


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
    name = "curr_w2l_ptclass" if abs(wsup - 0.70) < 1e-9 else \
        f"curr_w2l_ptclass_w{int(round(wsup * 100)):02d}"
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
                        csv_path=os.path.join(HERE, "ptclass_trials.csv"))
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
    n = 100
    for a in sys.argv[1:]:
        if a.startswith("--n="):
            n = int(a[4:])
    r70 = run_round(0.70, None, n)
    if not r70["survive30"]:
        print("stage-1 (wsup 0.70) did not survive fresh verify; stopping",
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
    print("WEAN COMPLETE through wsup=0.10 — run basin gate next "
          "(_basin_condsolve.py pattern)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
