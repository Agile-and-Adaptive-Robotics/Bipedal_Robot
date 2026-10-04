"""Track-B gate: the default spinal network build must still produce the canonical counts
(410 neurons / 376 inputs / 1186 synapses). Run from the repo root with the myo python.

2026-10-03 REPAIR (wiring-rules audit): the original body called
``bn.build()`` with no arguments and unpacked ``net, meta = ...`` - the
current ``build_network.build(model_actuators)`` API takes a required
actuator list and returns a single SpinalNetwork, so the gate raised
TypeError and could never pass. Repaired to the same check the repo's own
defaults gate performs (reports_20260925/audit/audit_defaults_gate.py):
92-actuator list from muscle_map, counts via net.net.get_num_neurons()/
get_num_inputs_actual()/get_num_connections(). The canonical counts are
unchanged."""
import os
import sys

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())
import build_network as bn  # noqa: E402
import muscle_map as MM  # noqa: E402


def actuator_names():
    acts = []
    for base in MM._GROUPS_BY_NAME:
        acts.append(base + "_r")
        acts.append(base + "_l")
    return acts


net = bn.build(actuator_names())
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print("counts:", counts)
if counts != (410, 376, 1186):
    print("COUNTS GATE FAIL")
    sys.exit(1)
print("COUNTS GATE PASS")
