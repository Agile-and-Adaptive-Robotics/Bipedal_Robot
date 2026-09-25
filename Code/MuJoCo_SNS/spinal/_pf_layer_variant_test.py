"""Per-PF-layer contact variant test (2026-09-24 evening; Ben: "build
it", coexists with the per-joint layering). Source drawing =
Circuit_rules_CONNECTOME_md__connectome.json (97n/109e).

gate 1: default build (all five new keys 0) -> network counts recorded
        (the new code paths are gain-gated, so this IS the pre-variant
        topology; compared against the same build with the keys forced).
gate 2: variant-on build (joint_pf + ia_in + full_rules + the five new
        gains) -> builds clean, TOEDF_r/l exist, expected edge growth.
gate 3: bit-exact ground regression (s3 stage, full_rules 0).
        BASELINE -160.23425729850192 re-recorded 2026-09-24 evening on
        PRE-VARIANT code (reverse-patched via _pfvar_toggle.py A/B):
        the variant build is bit-IDENTICAL to it. The older
        _fullrules_test.py reference -148.6878643 (09-21) is STALE --
        the s3 eval moved with the 09-23/24 physics changes, before
        this variant existed.
"""
import json
import os
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).parent
NEW_KEYS = ("heel_pf_layer", "toe_df_inh", "heel_in_f_exc",
            "ia_pf_f", "ii_pf_f")

import muscle_map as MM
import params as P
import build_network as BN


def actuator_names():
    acts = []
    for base in MM._GROUPS_BY_NAME:
        acts.append(base + "_r")
        acts.append(base + "_l")
    return acts


def counts(net):
    n = net.net
    return (n.get_num_neurons(), n.get_num_inputs_actual(),
            n.get_num_connections())


ok = True

print("=== gate 1: defaults build (new keys all 0) ===")
assert all(P.G[k] == 0.0 for k in NEW_KEYS), "new keys must default 0"
net = BN.build(actuator_names())
c_def = counts(net)
print(f"counts neurons/inputs/synapses = {c_def}")

print("=== gate 2: variant-on build ===")
# base config: the pre-existing features the variant rides on
P.G["joint_pf"] = 1.0
P.G["ia_in"] = 0.5
P.G["full_rules"] = 1.0
P.G["heel_rge"] = 0.5
P.G["toe_rge"] = 0.5
net2 = BN.build(actuator_names())
c_base = counts(net2)
# + the five variant keys
P.G["heel_pf_layer"] = 0.5
P.G["toe_df_inh"] = 5.0
P.G["heel_in_f_exc"] = 0.5
P.G["ia_pf_f"] = 0.5
P.G["ii_pf_f"] = 0.5
net3 = BN.build(actuator_names())
c_on = counts(net3)
delta = tuple(o - b for o, b in zip(c_on, c_base))
have = "TOEDF_r" in net3.idx and "TOEDF_l" in net3.idx
print(f"base counts = {c_base}")
print(f"+variant    = {c_on}  delta {delta}  TOEDF r/l: {have}")
# expected: +2 neurons (TOEDF r/l), +0 inputs, new synapses
# (heel->PF_IN_E, toe->TOEDF->ANK-F, flexor IaIN/IIX->F per side)
ok2 = have and delta[0] == 2 and delta[1] == 0 and delta[2] >= 20
print("GATE 2:", "PASS" if ok2 else "FAIL")
ok = ok and ok2
# restore defaults for gate 3
for k in NEW_KEYS:
    P.G[k] = 0.0
P.G["joint_pf"] = 0.0
P.G["ia_in"] = 0.0
P.G["full_rules"] = 0.0
P.G["heel_rge"] = 0.0
P.G["toe_rge"] = 0.0

print("=== gate 3: bit-exact ground regression (s3, full_rules 0) ===")
import _curriculum as C
import runner as R

mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
os.environ["AARL_NPZ"] = "spinal_run_pfvar.npz"
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(3, {**mul, **st3, "full_rules": 0.0})
P.G["full_rules"] = 0.0
m = R.main(["--eval", "--drive", repr(st3["drive"])])
ks = m["kine_score"]
ok3 = abs(ks - (-160.23425729850192)) < 1e-6
print(f"kine {ks!r} vs baseline -160.23425729850192")
print("GATE 3:", "PASS" if ok3 else "FAIL")
ok = ok and ok3
print("\nOVERALL:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
