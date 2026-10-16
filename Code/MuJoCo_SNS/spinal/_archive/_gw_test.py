"""Deng 2022 G_W test, RUNNER version (Ben: 'instead of the adapt
neuron, maybe weakly self-excitatory'): full 14 s deafferented air run
per config (the exact stage-1 objective gate), stage-1 winner params.
A. ADAP on,  G_W 0     (current)   -> expect ~2 rises (known 43.9)
B. ADAP off, G_W sweep              -> G_W-only escape mode?
C. ADAP on,  G_W sweep              -> Deng-faithful combo + frequency
"""
import io
import json
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
trial = max((t for t in st.trials if t.value is not None),
            key=lambda t: t.value)
mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]
p = {**mul, **trial.params}
print(f"stage-1 winner: trial {trial.number} value {trial.value:.3f} "
      f"params {trial.params}")


def run(adap, gw):
    C.set_stage(1, p)
    P.G["rg_adapt_inh"] = adap
    P.G["rg_weak_exc"] = gw
    R.main(["--no-ground", "--no-afferents", "--no-interleg", "--time",
            "14", "--drive", repr(p["drive"])])
    z = np.load("spinal_run.npz", allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    m = (t >= 5.0) & (t <= 17.0)
    knee = q[m, 4]
    rge = neuro[m, 2]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    return rises, float(knee.min())


print(f"{'config':26s} {'rises':>5s} {'knee_min':>9s}")
for label, adap, gw in (
        ("A ADAP2.5 GW0 (current)", 2.5, 0.0),
        ("B ADAP0 GW0.5", 0.0, 0.5),
        ("B ADAP0 GW1.0", 0.0, 1.0),
        ("B ADAP0 GW2.0", 0.0, 2.0),
        ("C ADAP2.5 GW0.5", 2.5, 0.5),
        ("C ADAP2.5 GW1.0", 2.5, 1.0),
        ("C ADAP2.5 GW2.0", 2.5, 2.0)):
    rises, kmin = run(adap, gw)
    print(f"{label:26s} {rises:5d} {kmin:9.1f}")
