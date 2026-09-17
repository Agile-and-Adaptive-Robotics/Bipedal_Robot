"""Verify the 2026-09-16 connection fixes in build_network.py:
(1) heel pathway wired ONCE per target (was duplicated = 2x conductance);
(2) IaIN receives its Ia afferent (Ia -> IaIN exc);
(3) stage-1 equivalence: with all conditional gains 0 the conditional
    blocks (HEEL/TOE/LBIN/IaIN/RC/PRESET/PREA) are entirely absent, so
    the running stage-1 study is bit-identical under the fixed code.
"""
import io
import sys
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import params as P
import build_network as bn

ACTS = ["vas_lat_r", "semimem_r", "vas_lat_l", "semimem_l"]


def groups(gains):
    keep = {k: P.G[k] for k in gains}
    P.G.update(gains)
    net = bn.build(ACTS, interleg=True)
    P.G.update(keep)
    names = [p["name"] for p in net.net.populations]
    g = defaultdict(int)
    for c in net.net.connections:
        s, d = names[c["source"]], names[c["destination"]]
        sign = "exc" if c["params"].get("reversal_potential", 0) > -1e-6 \
            else "inh"
        g[(s, d, sign)] += 1
    return names, g


# --- (3) stage-1 gains: conditional topology absent
stage1 = dict(phase_reset_e=0.0, phase_reset_f=0.0, f1_kneext_inh=0.0,
              f1_anklepf_inh=0.0, renshaw=0.0, ia_in=0.0, heel_rge=0.0,
              toe_rge=0.0, ib_rge=0.0)
names, g = groups(stage1)
cond = [n for n in names if n.split("_")[0] in
        ("HEEL", "TOE", "LBIN", "IaIN", "RC", "PRESET", "PREA")]
assert not cond, f"conditional populations present at stage-1 gains: {cond}"
print(f"stage-1 config: {len(names)} populations, conditional blocks "
      "ABSENT (stage-1 study stays valid)")

# --- (1)+(2) tuned gains: heel once, IaIN fed, ALL new central pathways
tuned = dict(phase_reset_e=0.5, phase_reset_f=0.5, f1_kneext_inh=0.6,
             f1_anklepf_inh=0.6, renshaw=0.5, ia_in=0.6, heel_rge=0.6,
             toe_rge=0.4, ib_rge=0.6, rg_weak_exc=0.4,
             ib_e_central=0.5, ia_f_central=0.5, ii_f_central=0.3,
             ii_e_central=0.3, ia_f_contra_f=0.4, v3_to_ibexc=0.4)
names, g = groups(tuned)
heel_e = [(k, v) for k, v in g.items()
          if k[0].startswith("HEEL") and k[2] == "exc"]
heel_i = [(k, v) for k, v in g.items()
          if k[0].startswith("HEEL") and k[2] == "inh"]
assert all(v == 1 for _, v in heel_e + heel_i), \
    f"heel still duplicated: {heel_e + heel_i}"
print(f"heel edges (single per target): {heel_e + heel_i}")
iain_aff = [(k, v) for k, v in g.items()
            if k[0].startswith("Ia_") and k[1].startswith("IaIN")]
assert len(iain_aff) == 4 and all(v == 1 for _, v in iain_aff), \
    f"Ia->IaIN afferent edges wrong: {iain_aff}"
print(f"Ia->IaIN afferent edges (x4 pools, once each): {iain_aff}")
# new central-pathway groups present exactly once per pair
for want in (("Ib_vas_lat_r", "PF_E1_r", "exc"),
             ("Ib_vas_lat_r", "RG_E_r", "exc"),
             ("Ib_vas_lat_r", "InE_r", "exc"),
             ("II_vas_lat_r", "RG_E_r", "exc"),
             ("Ia_semimem_r", "PF_F1_r", "exc"),
             ("II_semimem_r", "RG_F_r", "exc"),
             ("Ia_semimem_r", "RG_F_l", "inh"),
             ("HEEL_r", "PF_E1_r", "exc"),
             ("TOE_l", "PF_E2_l", "exc"),
             ("CIN_E_r", "IBEXC_knee_ext_l", "exc")):
    assert g.get(want) == 1, f"missing/duplicated central edge {want}: " \
        f"{g.get(want)}"
print("all NEW central-pathway edges verified (Ib/II->E-centers, "
      "Ia/II->F-centers, Ia->contra-F, HEEL/TOE->PF_E, V3->contra-IBEXC)")
print("ALL FIX CHECKS PASS")
