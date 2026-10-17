r"""Variant default-build census — 2026-10-16 drive-implementation verify
stage (dwfrun overnight walk campaign).

Prints the (neurons, inputs, synapses) census of the DEFAULT build for
AARL_NET in {w2lvar, syn6} with all params.G knobs at defaults (gamma_l/r
= alpha_l/r = 0). One process per variant:

    set AARL_NET=w2lvar&& D:\Anaconda\envs\myo\python.exe variant_census_20261016.py [spinal_dir]

spinal_dir defaults to D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal;
pass an alternative (e.g. a HEAD worktree) to census that tree. Build
recipe mirrors reports_20260925/audit/audit_defaults_gate.py.
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
spinal = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, spinal)
os.chdir(spinal)
variant = os.environ.get("AARL_NET", "").strip().lower()
assert variant in ("w2lvar", "syn6"), \
    f"set AARL_NET to w2lvar or syn6 (got {variant!r})"

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
print(f"AARL_NET={variant} default counts neurons/inputs/synapses = {counts}")
sys.exit(0)
