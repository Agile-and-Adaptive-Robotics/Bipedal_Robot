"""Track-B gate: the default spinal network build must still produce the canonical counts
(410 neurons / 376 inputs / 1186 synapses). Run from the repo root with the myo python."""
import os
import sys

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())
import build_network as bn  # noqa: E402

net, meta = bn.build()
counts = (len(net.neurons), len(net.inputs), len(net.synapses))
print("counts:", counts)
if counts != (410, 376, 1186):
    print("COUNTS GATE FAIL")
    sys.exit(1)
print("COUNTS GATE PASS")
