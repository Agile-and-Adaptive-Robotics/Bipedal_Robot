# Did the perturbation actually land? Compare saved V traces.
import numpy as np

base = np.load("deng_mujoco_base.npz")
for cond in ("p1_rg", "p2_pf", "p3_rg_long"):
    d = np.load(f"deng_mujoco_{cond}.npz")
    dvE = np.abs(d["vE"] - base["vE"]).max()
    dvF = np.abs(d["vF"] - base["vF"]).max()
    w = (base["t"] >= 1.9) & (base["t"] <= 3.0)
    dvF_w = np.abs(d["vF"][w] - base["vF"][w]).max()
    print(f"{cond}: max|dV_ext| {dvE:.3e}  max|dV_flx| {dvF:.3e}  "
          f"(near window {dvF_w:.3e})")
print("flx V in base window 1.9-3.0 s: min %.1f max %.1f mV"
      % (base["vF"][w].min(), base["vF"][w].max()))
