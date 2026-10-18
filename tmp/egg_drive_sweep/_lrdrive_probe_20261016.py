"""L/R hand-drive topology probe at the s3k WINNER config (campaign
2026-10-16). Builds the production config (curr_s3k_nocross t34 via
set_stage(4, ...)) with the four hand-drive knobs ON (gamma_l/r,
alpha_l/r = 0.8/0.4/0.5/0.5) and dumps the GAMMA/ALPHA cells, ports,
and out-edges with conductances. Run with PYTHONIOENCODING=utf-8."""
import json
from collections import Counter

import mujoco
import build_network as bn
import params as P
import pathlib
import os

MODEL = pathlib.Path(os.environ.get(
    "AARL_MODEL",
    r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
    r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"))
model = mujoco.MjModel.from_xml_path(str(MODEL))
acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(model.nu)]

import _curriculum as C
base = dict(json.loads(
    open("best_walk_params_v10.json", encoding="utf-8").read())
    ["multipliers"])
base["renshaw"] = 0.5
s3k = json.loads(open(
    "reports_20260923/s3k_trial34_full_params.json",
    encoding="utf-8").read())["params"]
base.update(s3k)
C.set_stage(4, base)                 # production winner config
P.G["gamma_l"], P.G["gamma_r"] = 0.8, 0.4
P.G["alpha_l"], P.G["alpha_r"] = 0.5, 0.5
net = bn.build(acts, interleg=True)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print(f"winner + 4 drive knobs: counts = {counts}")
print("(winner without knobs = 878/382/8396; expect +4n/+4p/+96s)")
pop = [p["name"] for p in n.populations]
ga = [(pop[c["source"]], pop[c["destination"]],
       float(c["params"]["max_conductance"]))
      for c in n.connections if pop[c["source"]].startswith(("GAMMA_", "ALPHA_"))]
agg = Counter((s, d.split("_")[0] + "_" + d.split("_")[1]
               if d.startswith("MN_") else d, round(g, 4)) for s, d, g in ga)
for (s, d, g), v in sorted(agg.items()):
    print(f"  exc g={g:<8g} {s} -> {d} x{v}")
print("drive ports:", [p for p in net.inputs
                       if p.startswith(("GAMMA_c", "ALPHA_c"))])
print("runner feed recipe: u[GAMMA_c_s]=drive*G[gamma_s], "
      "u[ALPHA_c_s]=drive*G[alpha_s] (runner.py step loop)")
