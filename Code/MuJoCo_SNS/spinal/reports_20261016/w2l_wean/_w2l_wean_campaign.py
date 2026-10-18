r"""W2L robot WEAN LADDER campaign driver (overnight campaign 2026-10-16,
egg; stage "feedback optimizer (robot)").

The human-side recipe applied to the W2L robot: optuna-tune the FEEDBACK
gains at the bootstrap drive (te=4/tf=5 nA, dscale 1.0), then force the
drive scale down across rounds (te=4*s, tf=5*s, ratio preserved; the
10 nA/10 ms antiphase kickoff is part of the model, kept verbatim) and
re-tune at each rung. Report the lowest dscale where ALTERNATING
PROGRESSION survives (finite, no fall, >= 6 primary stride cycles BOTH
sides, inter-leg phase 0.35-0.65, forward disp > 0.05 m -- the ask's
criterion), plus where it dies.

Search space (feedback only; protocol/body knobs pinned at the bootstrap
winner): aff_heel, aff_ib, ibnA, ctn_amp, cref, comm, and a toe-chain
multiplier toe_m (g_toe_in/toe_df/toe_df_mn scaled together, preserving
Ben's rules ratios 1 : 5 : 2.749).

Studies: curr_s6_drive_robot_r{DSC} in reports_20261016\w2l_wean\
optuna_w2l_robot.db. Every trial's diagnostics land in wean_trials.csv
(_w2l_wean_eval) AND as optuna user_attrs (anti-lottery).

Usage (nohup, egg):
    D:\Anaconda\envs\myo\python.exe _w2l_wean_campaign.py
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _w2l_wean_eval as EV          # noqa: E402

import optuna                        # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)

DB = "sqlite:///" + os.path.join(HERE, "optuna_w2l_robot.db").replace("\\", "/")
SEED = 20261016
DUR_TUNE, DUR_VERIFY = 20.0, 30.0

SPACE = dict(
    g_heel=(0.05, 2.5),     # aff_heel: heel -> InE/InF + PF-E INs
    g_ib=(0.0, 2.5),        # aff_ib: Ib grp -> RG-E/InE/PF-E
    ibnA=(0.0, 3.0),        # Ib sensor current scale (nA)
    ctn_amp=(0.5, 6.0),     # heel/toe contact-encoder amplitude (nA)
    cref=(20.0, 150.0),     # contact-encoder saturation force (N)
    comm=(0.2, 2.0),        # commissural gain scale
    toe_m=(0.0, 3.0),       # toe->TOEDF chain multiplier (ratios kept)
)
# the bootstrap winner config = Ben's rules gains (round-1 seed)
BOOT_PARAMS = dict(g_heel=0.5, g_ib=0.5, ibnA=1.0, ctn_amp=2.0,
                   cref=50.0, comm=1.0, toe_m=1.0)

LADDER = [1.0, 0.85, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20]
NTRIALS = {1.0: 36}          # default 26 otherwise


def knobs_of(dscale: float, p: dict, dur: float) -> dict:
    m = float(p.get("toe_m", 1.0))
    return dict(EV.BASE,
                dscale=dscale,
                te=EV.BOOT_TE * dscale, tf=EV.BOOT_TF * dscale,
                dur_walk=dur,
                g_heel=float(p["g_heel"]), g_ib=float(p["g_ib"]),
                ibnA=float(p["ibnA"]), ctn_amp=float(p["ctn_amp"]),
                cref=float(p["cref"]), comm=float(p["comm"]),
                g_toe_in=1.0 * m, g_toe_df=5.0 * m, g_toe_df_mn=2.749 * m)


def sname(dscale: float) -> str:
    return f"curr_s6_drive_robot_r{int(round(dscale * 100)):03d}"


def run_round(dscale: float, seed_params: dict | None) -> dict:
    study = optuna.create_study(
        study_name=sname(dscale), storage=DB, direction="maximize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=SEED))
    have = len(study.trials)
    if have == 0:
        for sp in ([seed_params] if seed_params else []) + [BOOT_PARAMS]:
            study.enqueue_trial(sp)
    need = NTRIALS.get(dscale, 26) - have

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi)
             for k, (lo, hi) in SPACE.items()}
        kn = knobs_of(dscale, p, DUR_TUNE)
        m = EV.evaluate(kn, tag=f"{sname(dscale)}#{trial.number}")
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
    out = dict(dscale=dscale, study=sname(dscale),
               n_trials=len(study.trials),
               params={k: round(float(best.params[k]), 5) for k in SPACE},
               obj20=float(best.value))
    # winner verification at the reporting duration (30 s)
    m30 = EV.evaluate(knobs_of(dscale, best.params, DUR_VERIFY),
                      tag=f"verify_r{dscale:.2f}")
    out["m30"] = {k: m30[k] for k in EV.METRIC_KEYS}
    out["survive30"] = bool(m30["survive"])
    path = os.path.join(HERE, f"round_best_r{int(round(dscale*100)):03d}.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    out["json_path"] = path
    print(f"== ROUND dscale={dscale:.2f} best obj20={out['obj20']:.1f} "
          f"survive30={out['survive30']} -> {path}", flush=True)
    return out


def main() -> int:
    t_start = time.time()
    summary_path = os.path.join(HERE, "ladder_summary.json")
    results = []
    seed = None
    fails = 0
    death = None
    for s in LADDER:
        r = run_round(s, seed)
        results.append(r)
        json.dump(dict(rounds=results, t_start=t_start),
                  open(summary_path, "w", encoding="utf-8"), indent=1)
        if r["survive30"]:
            seed = r["params"]       # wean the WINNER down, not a stranger
            fails = 0
        else:
            fails += 1
            if death is None:
                death = s
            if fails >= 2:
                break               # death confirmed by two rungs
    # one midpoint refinement between last survivor and the death rung
    last_pass = next((r["dscale"] for r in reversed(results)
                      if r["survive30"]), None)
    if death is not None and last_pass is not None:
        mid = round((last_pass + death) / 2.0, 3)
        if not any(abs(r["dscale"] - mid) < 1e-9 for r in results):
            r = run_round(mid, seed)
            results.append(r)
            json.dump(dict(rounds=results, t_start=t_start),
                      open(summary_path, "w", encoding="utf-8"), indent=1)

    # finalist = winner of the LOWEST surviving rung
    surv = [r for r in results if r["survive30"]]
    finalist = min(surv, key=lambda r: r["dscale"]) if surv else None
    json.dump(dict(rounds=results, finalist_dscale=(
        finalist["dscale"] if finalist else None),
        death_scale=death, wall_min=round((time.time() - t_start) / 60, 1)),
        open(summary_path, "w", encoding="utf-8"), indent=1)
    print(f"CAMPAIGN DONE finalist_dscale="
          f"{finalist['dscale'] if finalist else None} death={death} "
          f"({(time.time() - t_start) / 60:.1f} min)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
