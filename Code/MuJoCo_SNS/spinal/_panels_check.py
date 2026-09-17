"""Smoke check: build the representative networks from spinal_layers.py
and assert the NaP-architecture wiring (2026-09-16): InE/InF laminated
RG, G_W mutual excitation, PF_IN laminated cross-reciprocal, IaIN with
afferent leg, mutual RCs, NO ADAP / PRESET / PFA anywhere."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import params as P

P.G["rg_weak_exc"] = 0.4
P.G["f1_kneext_inh"] = 0.6
P.G["renshaw"] = 0.5
P.G["ia_in"] = 0.6

from spinal_layers import (MotorColumnNetwork, PatternFormationNetwork,
                           RhythmGeneratorNetwork)


def inv(net):
    names = [p["name"] for p in net.populations]
    out = []
    for c in net.connections:
        s, d = names[c["source"]], names[c["destination"]]
        sign = "exc" if c["params"].get("reversal_potential", 0) > -1e-6 \
            else "inh"
        out.append(f"{s} -{sign}-> {d}")
    return sorted(out)


for cls, args in ((RhythmGeneratorNetwork, ("r",)),
                  (PatternFormationNetwork, ("r",)),
                  (MotorColumnNetwork, ("r",))):
    net = cls(*args)
    names = [p["name"] for p in net.populations]
    print(f"=== {cls.__name__} ({len(names)} pops, "
          f"{len(net.connections)} conns)")
    for line in inv(net):
        print("  ", line)

rgi = inv(RhythmGeneratorNetwork("r"))
assert "RG-E_r -inh-> RG-F_r" not in rgi, "direct RG inhibition present"
assert "RG-E_r -exc-> InE_r" in rgi and "InE_r -inh-> RG-F_r" in rgi
assert "RG-F_r -exc-> InF_r" in rgi and "InF_r -inh-> RG-E_r" in rgi
assert "RG-E_r -exc-> RG-F_r" in rgi and "RG-F_r -exc-> RG-E_r" in rgi, \
    "G_W weak mutual excitation missing"
assert not any("ADAP" in x or "PRESET" in x or "PREA" in x
               for x in rgi), "retired pathway present in RG"
rg_pops = [p["name"] for p in RhythmGeneratorNetwork("r").populations]
assert not any("ADAP" in n or "PRESET" in n for n in rg_pops)
pfn = PatternFormationNetwork("r")
pfi = inv(pfn)
assert "PF_E1_r -inh-> PF_F1_r" not in pfi, "direct PF edge still present"
assert "PF_E1_r -exc-> PF_IN_E_r" in pfi and \
    "PF_IN_E_r -inh-> PF_F1_r" in pfi
pf_pops = [p["name"] for p in PatternFormationNetwork("r").populations]
assert not any("PFA" in n for n in pf_pops), "PFA still present"
pf_inputs = [p["name"] for p in pfn.inputs]
assert not any("DRIVE->PF" in n for n in pf_inputs), \
    f"retired DRIVE->PF input still rendered: {pf_inputs}"
moi = inv(MotorColumnNetwork("r"))
assert "IaIN_knee_ext_r -inh-> MN_knee_flex_r" in moi
assert "Ia_knee_ext_r -exc-> IaIN_knee_ext_r" in moi
assert "RC_knee_ext_r -inh-> RC_knee_flex_r" in moi
assert "RC_knee_flex_r -inh-> RC_knee_ext_r" in moi
print("ALL NAP-ARCHITECTURE ASSERTS PASS")
