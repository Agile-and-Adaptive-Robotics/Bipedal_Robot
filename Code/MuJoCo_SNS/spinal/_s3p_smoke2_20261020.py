import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
os.environ["AARL_NPZ"] = "_s3p_smoke2_20261020.npz"
import runner as R
import params as P
P.G["vest_ext"] = 0.3; P.G["vest_flex_inh"] = 0.1; P.G["vest_prop"] = 0.4
P.G["coact"] = 0.6
os.environ["AARL_STANCE_W"] = "5.0"
os.environ["AARL_KNEE_BIAS"] = "8.0"
m = R.main(["--stand-eval", "3", "--push", "60", "--push-time", "1.5",
            "--push-axis", "+x"])
print("PUSH+KNOBS:", {k: m[k] for k in ("bal_fell","bal_sway","bal_tilt_max",
      "bal_push_sway","bal_effort","bal_cocon","nan")})
m = R.main(["--stand-eval", "3", "--plat-dx", "0.06", "--push-time", "1.5"])
print("PLAT-DX:", {k: m[k] for k in ("bal_fell","bal_sway","bal_push_sway",
      "bal_effort","bal_cocon","nan")})
m = R.main(["--stand-eval", "3", "--plat-tilt", "6.0", "--push-time", "1.5"])
print("PLAT-TILT:", {k: m[k] for k in ("bal_fell","bal_sway","bal_push_sway",
      "bal_effort","bal_cocon","nan")})
import numpy as np
z = np.load("_s3p_smoke2_20261020.npz", allow_pickle=True)
print("npz keys:", sorted(z.files), "ctrl shape:", z["ctrl"].shape)
