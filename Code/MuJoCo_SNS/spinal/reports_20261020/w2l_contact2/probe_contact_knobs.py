"""ON-state probe for the 2026-10-20 W2L contact-pathway knobs
(contra_flex, contact_ank_mn) across the model families.

For AARL_NET in {unset, w2lvar, syn6, spiking}: build with BOTH knobs at
0.5 (+ heel_rge/toe_rge 0.5 so the contact INs carry signal), report the
census + synapse DELTA vs the family's default census, and one finite
step. Default-state bit-identity is the separate gate
(drive_defaults_gate / variant_census); this proves the new edges appear
when the knobs are on and the step stays finite.

    D:\\Anaconda\\envs\\myo\\python.exe probe_contact_knobs.py [spinal_dir]
"""
from __future__ import annotations

import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
spinal = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, spinal)
os.chdir(spinal)

DEFAULT_CENSUS = {"": (410, 376, 1186), "w2lvar": (888, 382, 7866),
                  "syn6": (794, 382, 3033), "spiking": (422, None, 1186)}

import params as P                     # noqa: E402
import muscle_map as MM                # noqa: E402
import build_network as BN             # noqa: E402

variant = os.environ.get("AARL_NET", "").strip().lower()
P.G["contra_flex"] = 0.5
P.G["contact_ank_mn"] = 0.5
P.G["heel_rge"] = 0.5
P.G["toe_rge"] = 0.5

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")

net = BN.build(acts)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
d0 = DEFAULT_CENSUS[variant]
d_syn = counts[2] - d0[2]
d_neu = counts[0] - d0[0]
print(f"AARL_NET={variant or 'standard'} ON-knob census = {counts} "
      f"(default {d0}; delta neurons {d_neu:+d}, synapses {d_syn:+d})")
# expected delta: +2 crossed heel->RG_F edges; +2 per ankle_pf muscle MN
ankle_pf = sum(1 for a in acts if MM.classify(a).groups[0] == "ankle_pf")
print(f"ankle_pf muscles in this actuator set: {ankle_pf} "
      f"(expected synapse delta 2 + 2*{ankle_pf} = {2 + 2 * ankle_pf})")

u = net.make_inputs()
u[net.input_index("HEEL_c_r")] = 1.0
V = net.step(u)
import numpy as np                    # noqa: E402
assert np.isfinite(V).all(), "non-finite step"
print("one-step finite: OK")
