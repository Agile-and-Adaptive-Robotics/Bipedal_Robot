"""goal4 gates (a) + (b) for the w2lvar variant (2026-09-25).

(a) defaults gate: env AARL_NET unset -> standard build EXACTLY
    410 neurons / 376 inputs / 1186 synapses.
(b) variant gate: env AARL_NET=w2lvar -> variant builds + compiles;
    report neuron/synapse/input counts + spot checks.
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
SPINAL = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

assert "AARL_NET" not in os.environ, "run with AARL_NET unset"


def actuator_names():
    import muscle_map as MM
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

print("=== gate (a): defaults build, AARL_NET unset ===")
import build_network as BN
net = BN.build(actuator_names())
c_def = counts(net)
print(f"counts neurons/inputs/synapses = {c_def}")
ok_a = c_def == (410, 376, 1186)
print("GATE A:", "PASS" if ok_a else "FAIL")
ok = ok and ok_a

print()
print("=== gate (b): variant build, AARL_NET=w2lvar ===")
os.environ["AARL_NET"] = "w2lvar"
import importlib
import params as P
importlib.reload(P)
importlib.reload(BN)
net2 = BN.build(actuator_names())
c_var = counts(net2)
print(f"variant counts neurons/inputs/synapses = {c_var}")
print(f"delta vs default = "
      f"{tuple(v - d for v, d in zip(c_var, c_def))}")

idx = net2.idx
need = ["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
        "PF_HIP-E_r", "PF_HIP-F_r", "PF_KNEE-E_r", "PF_KNEE-F_r",
        "PF_IN_HIP-E_r", "PF_IN_KNEE-E_r", "TOEDF_r",
        "V2a_r", "V0V_r", "V0D_r", "V3E_r", "INI_l",
        "MN_vas_lat_r", "RC_vas_lat_r", "IaIN_vas_lat_r",
        "IBEXC_knee_ext_r", "HEEL_r", "TOE_r", "LBIN_r"]
missing = [x for x in need if x not in idx]
print(f"spot-check names present: {len(need) - len(missing)}/{len(need)}"
      + (f"  MISSING: {missing}" if missing else ""))
# parent commissural block must NOT be in the variant
leak = [x for x in ("CIN_F_r", "CIN_E_r") if x in idx]
print(f"parent CIN block absent: {not leak}" + (f"  LEAKED: {leak}" if leak else ""))
# aliases point at the merged cell
alias_ok = (idx.get("PF_ANK-F_r") == idx.get("PF_KNEE-F_r")
            and idx.get("PF_ANK-E_r") == idx.get("PF_KNEE-E_r"))
print(f"watch aliases -> merged cell: {alias_ok}")
# compiled backend is the fixed-tau_h stepper
import build_network_w2lvar as WV
print(f"compiled backend class: "
      f"{type(net2.compiled).__module__}.{type(net2.compiled).__name__}")
print(f"expected  {WV.SNS_NumpyFixedTau.__module__}.{WV.SNS_NumpyFixedTau.__name__}")
ok_b = (not missing) and (not leak) and alias_ok and \
    isinstance(net2.compiled, WV.SNS_NumpyFixedTau)
print("GATE B:", "PASS" if ok_b else "FAIL")
ok = ok and ok_b

# save the variant wiring inventory for the report
edges = {}
inv = {"counts_default": c_def, "counts_variant": c_var,
       "inputs": net2.inputs}
import json
with open(os.path.join(HERE, "gate_b_inventory.json"), "w",
          encoding="utf-8") as f:
    json.dump(inv, f, indent=1)
print("\nOVERALL:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
