r"""Standard-net default-build census — 2026-10-16 drive-implementation
verify stage. Location-parametrized twin of
reports_20260925/audit/audit_defaults_gate.py (whose SPINAL is hardcoded
to the real repo, which would defeat worktree HEAD runs).

    D:\Anaconda\envs\myo\python.exe default_census_20261016.py [spinal_dir]

spinal_dir defaults to the real repo. Asserts AARL_NET unset, builds the
standard network at params.G defaults (gamma_l/r = alpha_l/r = 0), prints
the (neurons, inputs, synapses) census, exits 0 iff == (410, 376, 1186).
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
spinal = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, spinal)
os.chdir(spinal)
assert os.environ.get("AARL_NET", "") == "", \
    f"AARL_NET must be unset for the standard census, got {os.environ.get('AARL_NET')!r}"

import muscle_map as MM
acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")

import build_network as BN
net = BN.build(acts)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print(f"spinal_dir: {spinal}")
print(f"standard default counts neurons/inputs/synapses = {counts}")
print("GATE:", "PASS" if counts == (410, 376, 1186) else "FAIL")
sys.exit(0 if counts == (410, 376, 1186) else 1)
