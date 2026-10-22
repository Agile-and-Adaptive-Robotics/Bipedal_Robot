"""Quiet-eval diagnostic for the syn6 s3q winner (trial 5): run-length,
foot-force duty over the steady window, COM height - separates true
weight-bearing from the mean-force bounce artifact."""
import io, importlib, json, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
os.environ["AARL_NPZ"] = "s3q_syn6_quietdiag_20261021.npz"
os.environ["AARL_NET"] = "syn6"
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
best = optuna.load_study(study_name="curr_syn6_s3q_perturb",
                         storage="sqlite:///optuna_syn6.db").best_trial
C = importlib.import_module("_curriculum_syn6")
import params as P
prev = json.loads(open("best_walk_params_v10.json", encoding="utf-8").read())
p = {**dict(prev["multipliers"]), **best.params}
p["renshaw"] = 0.5
C.set_stage(3, p)
P.G["coact"] = float(p.get("coact", 0.0))
if float(p.get("stance_width", 0.0)) > 0.0:
    os.environ["AARL_STANCE_W"] = repr(float(p["stance_width"]))
if float(p.get("knee_bias", 0.0)) > 0.0:
    os.environ["AARL_KNEE_BIAS"] = repr(float(p["knee_bias"]))
os.environ["AARL_PELVIS_TY"] = repr(float(p["pelvis_ty"]))
m = C.R.main(["--stand-eval", "8", "--push-time", "4.0",
              "--rig-scale", repr(float(p["rig_scale"]))])
print(f"quiet: t_end={m['t_end']:.2f}s fell={m['bal_fell']} "
      f"foot_load={m['bal_foot_load']:.3f} com_z_min={m['bal_com_z_min']:.3f}")
import numpy as np
BW = 72.0 * 9.81
z = np.load("s3q_syn6_quietdiag_20261021.npz", allow_pickle=True)
t = z["t"]; ct = z["contact"]; com = z["com"]
n = min(int(np.sum(t > 0)) + 1, len(t))
i0 = int(np.searchsorted(t[:n], 4.0))
foot = 0.5 * (ct[i0:n, 0] + ct[i0:n, 1]) / BW
cz = com[i0:n, 2]
rec = dict(t_end=float(t[n-1]), n_samples_steady=int(n - i0),
           mean_foot_bw=float(foot.mean()), std_foot_bw=float(foot.std()),
           weightbearing_duty_80_120=float(np.mean((foot > 0.8) & (foot < 1.2))),
           light_duty_lt30=float(np.mean(foot < 0.3)),
           heavy_frac_gt150=float(np.mean(foot > 1.5)),
           com_z_mean=float(cz.mean()), com_z_std=float(cz.std()),
           com_z_min=float(cz.min()))
print(json.dumps({k: round(v, 3) for k, v in rec.items()}, indent=1))
json.dump(rec, open("_s3q_syn6_quietdiag_20261021.json", "w"), indent=1)
