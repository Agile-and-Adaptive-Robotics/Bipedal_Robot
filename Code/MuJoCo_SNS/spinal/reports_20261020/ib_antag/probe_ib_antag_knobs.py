"""ON-state probe for the 2026-10-20 evening syn6 Ib-central +
antagonist-IN audit knobs (ib_e_central, iain_pf_gate, pf_mixed_recip).

For AARL_NET=syn6: build with all three knobs at 0.5, report the census
+ synapse DELTA vs the syn6 default census (794/382/3033), assert the
delta equals the PREDICTED edge count computed independently from the
synergy tables (ib_e_central: per extensor-stance muscle 2 fixed targets
+ its E channels; iain_pf_gate: per muscle one edge per committed F
channel; pf_mixed_recip: 2 edges per side for the mixed channel), and
take one finite step. Default-state bit-identity is the separate gate
(variant_census_20261016 / drive_defaults_gate).

    D:\\Anaconda\\envs\\myo\\python.exe probe_ib_antag_knobs.py [spinal_dir]
"""
from __future__ import annotations

import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
spinal = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, spinal)
os.chdir(spinal)

assert os.environ.get("AARL_NET", "").strip().lower() == "syn6", \
    "set AARL_NET=syn6 for this probe"

import numpy as np                     # noqa: E402
import params as P                     # noqa: E402
import muscle_map as MM                # noqa: E402
import build_network as BN             # noqa: E402
from build_network_syn6 import (N_SYN, PRUNE_MUSCLES,   # noqa: E402
                                _load_synergy_tables)

P.G["ib_e_central"] = 0.5
P.G["iain_pf_gate"] = 0.5
P.G["pf_mixed_recip"] = 0.5

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")

# --- independent prediction from the synergy tables -----------------------
tabs = _load_synergy_tables()
pred = 0
pred_breakdown = {}
for side in ("r", "l"):
    tab = tabs[side]
    fam = tab["family"]
    mixed = [tab["stance_frac"][k] and False for k in range(N_SYN)]
    # mixed flags as the builder computes them (SYN6_MIXED_BAND)
    from build_network_syn6 import SYN6_MIXED_BAND
    sf = tab["stance_frac"]
    mixed = [SYN6_MIXED_BAND[0] < float(sf[k]) < SYN6_MIXED_BAND[1]
             for k in range(N_SYN)]
    e_chans = [k for k in range(N_SYN) if fam[k] == "E"]
    f_committed = [k for k in range(N_SYN)
                   if fam[k] == "F" and not mixed[k]]
    ib_edges = 0
    for act in acts:
        mi = MM.classify(act)
        if mi.side != side or mi.groups[0] not in MM.EXTENSOR_STANCE_GROUPS:
            continue
        if act in tab["names"]:
            row = tab["names"].index(act)
            chans = [k for k in e_chans if tab["W"][row, k] > 0.0]
        else:
            chans = e_chans  # not-in-basis fallback: all E channels
        ib_edges += 2 + len(chans)
    iain_edges = sum(1 for act in acts
                     if MM.classify(act).side == side) * len(f_committed)
    mixed_edges = 2 * sum(mixed)
    pred += ib_edges + iain_edges + mixed_edges
    pred_breakdown[side] = dict(ib_e_central=ib_edges,
                                iain_pf_gate=iain_edges,
                                pf_mixed_recip=mixed_edges,
                                pruned_not_in_basis=len(
                                    PRUNE_MUSCLES & {a for a in acts
                                                     if a.endswith(side)}))

net = BN.build(acts)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
d0 = (794, 382, 3033)
d_syn = counts[2] - d0[2]
print(f"AARL_NET=syn6 ON-knob census = {counts} (default {d0}; "
      f"delta neurons {counts[0] - d0[0]:+d}, synapses {d_syn:+d})")
print(f"predicted synapse delta = {pred:+d}  breakdown={pred_breakdown}")
assert counts[0] == d0[0] and counts[1] == d0[1], "neuron/input drift"
assert d_syn == pred, f"synapse delta {d_syn} != predicted {pred}"
print("synapse delta == prediction: OK")

u = net.make_inputs()
u[net.input_index("Ib_vas_lat_r")] = 1.0
u[net.input_index("Ib_soleus_l")] = 1.0
u[net.input_index("HEEL_c_r")] = 1.0
V = net.step(u)
assert np.isfinite(V).all(), "non-finite step"
print("one-step finite: OK")
