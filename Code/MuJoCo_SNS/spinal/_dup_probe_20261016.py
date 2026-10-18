"""Duplicate-edge probe (campaign 2026-10-16): count parallel synapses on
the SAME (src,dst) pair in the s3k winner build. Read-only."""
import io
import json
import sys
from collections import Counter

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import mujoco
import build_network as bn
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
C.set_stage(4, base)
net = bn.build(acts, interleg=True)
n = net.net
pop = [p["name"] for p in n.populations]
pairs = Counter((pop[c["source"]], pop[c["destination"]])
                for c in n.connections)
dups = {k: v for k, v in pairs.items() if v > 1}
print(f"total synapses {sum(pairs.values())}, distinct pairs "
      f"{len(pairs)}, pairs with PARALLEL synapses {len(dups)}")
from collections import defaultdict
fam = defaultdict(list)
for (s, d), v in dups.items():
    fam[(s.split('_')[0], d.split('_')[0])].append((s, d, v))
for k, lst in sorted(fam.items()):
    print(f"  {k[0]}->{k[1]}: {len(lst)} duplicated pairs, "
          f"multiplicities {sorted(v for _, _, v in lst)}")
    for s, d, v in sorted(lst, key=lambda x: -x[2])[:12]:
        print(f"    x{v} {s} -> {d}")
