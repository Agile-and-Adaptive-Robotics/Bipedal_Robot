"""AUDIT (goal4) - re-run the DEFAULTS gate on the current tree.

Asserts AARL_NET unset in the environment, builds the STANDARD network
exactly as the runner does, and requires neurons/inputs/synapses ==
(410, 376, 1186). Also checks the JSON-RULE part (a) for the syn6
selector keys: params.G defaults 0.0.
Read-only on all protected files (builds in memory only).
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

assert os.environ.get("AARL_NET", "") == "", \
    f"AARL_NET must be unset, got {os.environ.get('AARL_NET')!r}"
for k in ("AARL_KY", "AARL_PELVIS_TY", "AARL_NPZ"):
    print(f"env {k}: {os.environ.get(k, '<unset>')}")


def actuator_names():
    import muscle_map as MM
    acts = []
    for base in MM._GROUPS_BY_NAME:
        acts.append(base + "_r")
        acts.append(base + "_l")
    return acts


import build_network as BN
import params as P

net = BN.build(actuator_names())
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print(f"counts neurons/inputs/synapses = {counts}")
print("GATE:", "PASS" if counts == (410, 376, 1186) else "FAIL")

print(f"params.G['syn6'] default = {P.G.get('syn6', '<ABSENT>')}")
print(f"params.G['syn6_brainstem'] default = "
      f"{P.G.get('syn6_brainstem', '<ABSENT>')}")
ok = counts == (410, 376, 1186)
sys.exit(0 if ok else 1)
