"""Stage-1 repro + ADAP-necessity test on the FIXED build_network.
(1) Re-evaluate the recorded stage-1 winner (trial 22, best 57.436) on the
fixed code - stage-1 configs are equivalent pre/post fix, so this should
reproduce bit-exact and prove the equivalence end-to-end.
(2) Same config with G['rg_adapt_inh']=0 (ADAP loops disconnected) - does
the rhythm survive without the ADAP self-adaptation loops?"""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np
import optuna

import _curriculum as C
import params as P
import runner as R

optuna.logging.set_verbosity(optuna.logging.WARNING)
st = optuna.load_study(study_name="curr_s1_air_deaff",
                       storage="sqlite:///optuna_walk.db")
trial = next(t for t in st.trials if t.number == 22)
import json
mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]
p = {**mul, **trial.params}
print(f"trial 22 params: {trial.params}")


def run_once(tag):
    C.set_stage(1, p)
    R.main(["--no-ground", "--no-afferents", "--no-interleg", "--time",
            "14", "--drive", repr(p["drive"])])
    z = np.load("spinal_run.npz", allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    m = (t >= 5.0) & (t <= 17.0)
    ok = bool(np.all(np.isfinite(q[m])) and np.all(np.isfinite(neuro[m])))
    knee = q[m, 4]
    rge = neuro[m, 2]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    score = 3.0 * rises + 0.5 * (-float(knee.min()))
    print(f"[{tag}] finite={ok} rises={rises} knee_min="
          f"{float(knee.min()):.2f} score={score:.3f}")
    return score


s1 = run_once("as-is (fixed code)")
P.G["rg_adapt_inh"] = 0.0
s2 = run_once("ADAP disconnected (rg_adapt_inh=0)")
P.G["rg_adapt_inh"] = 2.5
print(f"\nrecorded stage-1 best: 57.436 | repro: {s1:.3f} | "
      f"no-ADAP: {s2:.3f}")
