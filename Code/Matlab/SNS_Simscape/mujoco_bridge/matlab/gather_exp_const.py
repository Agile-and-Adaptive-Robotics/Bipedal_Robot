"""Constants for the SNS x MuJoCo bridge experiments (E1/E2).

Dumps vas_med_r muscle properties + tuned-network conductances to
mujoco_bridge/matlab/exp_const.mat:
  L0/Lmin/Lmax (m, keyframe rest + actuator_lengthrange), Fmax (N),
  MN tau (s), ia_to_mn / ii_to_mn / ib_to_mn_inh conductances (uS, tuned),
  Esyn +- (mV), input-port indices for Ia/II/Ib/POST_vas_med_r.
"""
import json
from pathlib import Path

import mujoco
import numpy as np
from scipy.io import savemat

HERE = Path(__file__).resolve().parent                    # mujoco_bridge\matlab
MODEL = (HERE.parents[4] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
         / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3_simbridge.xml")
SPINAL = HERE.parents[3] / "MuJoCo_SNS" / "spinal"

m = mujoco.MjModel.from_xml_path(str(MODEL))
d = mujoco.MjData(m)
mujoco.mj_resetDataKeyframe(m, d, 0)
mujoco.mj_forward(m, d)

ACT = "vas_med_r"
ia = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, ACT)
L0 = float(d.actuator_length[ia])
lr = m.actuator_lengthrange[ia]
Fmax = float(m.actuator_gainprm[ia, 2])

J = json.loads((SPINAL / "spinal_net_export.json").read_text())
syn = {(s["src"], s["dst"]): s for s in J["synapses"]}
mn = f"MN_{ACT}"
g_ia = syn[(f"Ia_{ACT}", mn)]["g_uS"]
g_ii = syn[(f"II_{ACT}", mn)]["g_uS"]
g_ib = syn[(f"Ib_{ACT}", mn)]["g_uS"]
tau_mn = next(n["tau_s"] for n in J["neurons"] if n["name"] == mn)
inp = {r["port"]: k for k, r in enumerate(J["inputs"])}
i_ia = inp[f"Ia_{ACT}"]
i_ii = inp[f"II_{ACT}"]
i_ib = inp[f"Ib_{ACT}"]
i_post = inp[f"POST_{ACT}"]

savemat(HERE / "exp_const.mat", {
    "L0": L0, "Lmin": float(lr[0]), "Lmax": float(lr[1]), "Fmax": Fmax,
    "tau_mn": tau_mn, "g_ia": g_ia, "g_ii": g_ii, "g_ib": g_ib,
    "esyn_exc": 8.0, "esyn_inh": -5.0,
    "i_ia": i_ia, "i_ii": i_ii, "i_ib": i_ib, "i_post": i_post,
    "u_len": float(inp["u_len"]) if "u_len" in inp else -1,
})
print(f"{ACT}: L0={L0:.4f} m, range [{lr[0]:.4f},{lr[1]:.4f}] m, "
      f"Fmax={Fmax:.0f} N")
print(f"MN tau={tau_mn:.3f} s; tuned g: ia={g_ia:.3f} ii={g_ii:.3f} "
      f"ib={g_ib:.3f} uS; u-port idx ia/ii/ib/post = "
      f"{i_ia}/{i_ii}/{i_ib}/{i_post}")
