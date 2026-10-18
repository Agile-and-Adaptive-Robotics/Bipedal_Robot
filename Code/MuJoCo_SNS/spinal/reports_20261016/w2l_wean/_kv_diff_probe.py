r"""Decisive probe: is the evaluator-vs-testCLI trajectory difference a
CONFIG difference or process-context chaos? Print the effective T.KV
both ways and diff."""
from __future__ import annotations

import json
import os
import sys

W2L = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, W2L)

# path 1: evaluator-style (argv suppressed, dict update)
sys.argv = ["test_w2l_ground.py"]
import test_w2l_ground as T          # noqa: E402

BASE = dict(legcoll=0.0, sink=0.06, rig=0.25, walker_yield=1.5)
knobs = dict(BASE, te=4.0 * 0.6, tf=5.0 * 0.6, dur_walk=30.0,
             g_heel=0.58765, g_ib=0.11164, ibnA=1.32101,
             ctn_amp=0.71708, cref=46.37426, comm=0.41304,
             g_toe_in=1.0 * 2.21976, g_toe_df=5.0 * 2.21976,
             g_toe_df_mn=2.749 * 2.21976)
T.KV.update({k: v for k, v in knobs.items() if k in T.KV})
kv_eval = {k: T.KV[k] for k in knobs}

# path 2: CLI-style parse of the same values
cli_args = ("--phase=walk --te=2.4 --tf=3.0 --legcoll=0 --sink=0.06 "
            "--rig=0.25 --walker-yield=1.5 --ctn-amp=0.71708 "
            "--cref=46.37426 --ibnA=1.32101 --dur-walk=30 "
            "--g-heel=0.58765 --g-ib=0.11164 --comm=0.41304 "
            "--g-toe-in=2.21976 --g-toe-df=11.0988 "
            "--g-toe-df-mn=6.1021202400000005")
KV2 = dict(T.KV)
for arg in cli_args.split():
    k, v = arg[2:].split("=", 1)
    try:
        v = float(v)
    except ValueError:
        pass
    KV2[k.replace("-", "_")] = v

keys = set(kv_eval) | {"phase"}
diff = {k: (kv_eval.get(k), KV2.get(k)) for k in keys
        if kv_eval.get(k) != KV2.get(k)}
print("DIFFERING KEYS:", json.dumps(diff, indent=1, default=str)
      if diff else "NONE (configs identical -> difference is process "
      "context, not config)")
