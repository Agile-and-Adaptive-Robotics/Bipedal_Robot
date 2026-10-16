"""BPA-WALKER curriculum, SPIKING controller (2026-10-03 overnight
campaign). Same plant loop + score as _curriculum_bpa.py (imported, not
duplicated); the net is the SPIKING twin of the split-RG W2L controller
(build_w2l_split_net_spiking.py - the verified gait2392 spiking-mirror
cell mapping applied to the W2L split template).

Extra searched knobs vs the non-spiking twin: tau (RG SFA tau, the
analog tau_max_h role = the period knob) and thr_inc (mV/spike burst
termination) - the SPIKING_MIRROR_PLAN's "fresh short curriculum"
prediction: mapped parameters need retuning.

Usage: python _curriculum_bpa_spiking.py [n_trials]
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SPINAL = os.path.dirname(HERE)
for p in (HERE, SPINAL):
    if p not in sys.path:
        sys.path.insert(0, p)

# _curriculum_bpa wraps stdout at import (utf-8 console); do NOT wrap
# again here - a second TextIOWrapper around the same buffer closes the
# buffer when the first wrapper is garbage-collected (measured).
import _curriculum_bpa as NS  # noqa: E402  (run_bpa + score_bpa reuse)
import build_w2l_split_net_spiking as BS  # noqa: E402
import optuna  # noqa: E402

DB = NS.DB
STUDY = "curr_bpa_spk_20261003"
SEED = dict(NS.SEED, tau=0.25, thr_inc=4.0)


def _obj(trial):
    p = dict(
        te=trial.suggest_float("te", 2.0, 5.5),
        tf=trial.suggest_float("tf", 2.0, 5.5),
        tau=trial.suggest_float("tau", 0.15, 0.45),
        thr_inc=trial.suggest_float("thr_inc", 1.5, 6.0),
        cap=trial.suggest_float("cap", 0.4, 0.8),
        comm=trial.suggest_float("comm", 0.5, 1.5),
        stand_te=trial.suggest_float("stand_te", 1.0, 6.0),
        stand_tf=trial.suggest_float("stand_tf", 0.0, 4.0),
        support=trial.suggest_float("support", 0.5, 0.85),
        kmax_frac=trial.suggest_float("kmax_frac", 0.35, 0.7),
        tendon=trial.suggest_float("tendon", 0.02, 0.06),
        count=trial.suggest_int("count", 1, 3),
        dia=20.0, bias=0.02, damp_bpa=0.0,
    )
    net = BS.build(comm=p["comm"], tau=p["tau"],
                   thr_inc=p["thr_inc"],
                   template_path=os.path.join(SPINAL,
                                              "connectome_templates.json"))
    met = NS.run_bpa(net, p)
    return NS.score_bpa(met)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    st = optuna.create_study(direction="maximize", storage=DB,
                             study_name=STUDY, load_if_exists=True,
                             sampler=optuna.samplers.TPESampler(
                                 seed=32, n_startup_trials=8))
    if len(st.trials) == 0:
        st.enqueue_trial({k: v for k, v in SEED.items()
                          if k not in ("dia", "bias", "damp_bpa")})
        print("seeded", json.dumps(SEED), flush=True)
    st.optimize(_obj, n_trials=n, gc_after_trial=True)
    b = st.best_trial
    print(f"== BPA spk best {b.value:.3f} (trial {b.number})")
    print(json.dumps(b.params, indent=1))
    with open(os.path.join(SPINAL, "curriculum_bpa_spk.json"), "w",
              encoding="utf-8") as f:
        json.dump({"score": b.value, "params": b.params,
                   "trial": b.number, "study": STUDY, "db": DB,
                   "fixed": {k: SEED[k] for k in
                             ("dia", "bias", "damp_bpa")}}, f, indent=2)


if __name__ == "__main__":
    main()
