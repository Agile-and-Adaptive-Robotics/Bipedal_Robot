"""Audit probe 5: build the STOCK network at the s3k trial-34 winner's
topology-bearing gains and dump the rule-relevant edges (heel routing,
Ib reversal chain, IaIN motif, central afferent edges, commissurals)."""
import io
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, str(SPINAL))
os.chdir(SPINAL)
assert os.environ.get("AARL_NET", "") == ""

import muscle_map as MM  # noqa: E402
import params as P  # noqa: E402
import build_network as bn  # noqa: E402

win = json.loads((SPINAL / "reports_20260923" /
                  "s3k_trial34_full_params.json").read_text(encoding="utf-8"))
wp = win["params"]
# runner-contract gain keys (the runner maps drive/rg_nap_h/desc_e/... onto
# G/TAU differently; only the G-level topology keys matter here)
G_KEYS = ("heel_rge", "toe_rge", "ib_rge", "ib_e_central", "ia_f_central",
          "ii_f_central", "ii_e_central", "c1_gain", "v3_gain",
          "v3_to_ibexc", "ia_in", "ia_f_contra_f", "f1_anklepf_inh",
          "full_rules", "contact_onset", "contra_swing", "contra_kinh")
for k in G_KEYS:
    if k in wp:
        P.G[k] = float(wp[k])
P.G["f1_kneext_inh"] = float(wp.get("f1_kneext_inh", 0.0))
print("applied winner G keys:",
      {k: P.G[k] for k in G_KEYS if k in wp})

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
net = bn.build(acts)
n = net.net
print("s3k-winner-build counts:",
      (n.get_num_neurons(), n.get_num_inputs_actual(),
       n.get_num_connections()))
pop_names = [p["name"] for p in n.populations]
edges = defaultdict(list)
for c in n.connections:
    key = (pop_names[c["source"]], pop_names[c["destination"]])
    syn = c["params"]
    edges[key].append((float(syn.get("max_conductance", np.nan)),
                       "exc" if float(syn.get(
                           "reversal_potential", 0)) > -1e-6 else "inh"))


def edge(src, dst):
    lst = edges.get((src, dst), [])
    return lst[0] if lst else None


print("\n== heel routing (full_rules branch) ==")
for side in ("r", "l"):
    other = "l" if side == "r" else "r"
    print(f"[{side}] heel->InE:", edge(f"HEEL_{side}", f"InE_{side}"),
          " heel->InF:", edge(f"HEEL_{side}", f"InF_{side}"),
          " heel->RG_E:", edge(f"HEEL_{side}", f"RG_E_{side}"),
          " heel->contra InF:", edge(f"HEEL_{side}", f"InF_{other}"),
          " heel->PF_IN_E:", edge(f"HEEL_{side}", f"PF_IN_E_{side}"))
    print(f"[{side}] toe->InE:", edge(f"TOE_{side}", f"InE_{side}"),
          " toe->RG_E:", edge(f"TOE_{side}", f"RG_E_{side}"))
    print(f"[{side}] LBIN->RG_E:", edge(f"LBIN_{side}", f"RG_E_{side}"))

print("\n== Ib reversal chain (ankle_pf sample) ==")
print("Ib_soleus_r->IBIN_soleus_r:",
      edge("Ib_soleus_r", "IBIN_soleus_r"),
      " IBIN->MN:", edge("IBIN_soleus_r", "MN_soleus_r"),
      " IBIN<->IBIN(per_long):",
      edge("IBIN_soleus_r", "IBIN_per_long_r"))
print("Ib_soleus_r->IBEXC:", edge("Ib_soleus_r", "IBEXC_ankle_pf_r"),
      " RG_E->IBEXC:", edge("RG_E_r", "IBEXC_ankle_pf_r"),
      " IBEXC->MN:", edge("IBEXC_ankle_pf_r", "MN_soleus_r"),
      " IBEXC->LBIN:", edge("IBEXC_ankle_pf_r", "LBIN_r"))

print("\n== II motif (full_rules) ==")
print("II_soleus_r->IIX:", edge("II_soleus_r", "IIX_soleus_r"),
      " IIX->MN:", edge("IIX_soleus_r", "MN_soleus_r"),
      " II->IIIN:", edge("II_soleus_r", "IIIN_soleus_r"),
      " IIIN->ant MN(tib_ant):",
      edge("IIIN_soleus_r", "MN_tib_ant_r"))
print("direct II->MN (should be None under full_rules):",
      edge("II_soleus_r", "MN_soleus_r"))
print("direct Ib->MN (should be None under full_rules):",
      edge("Ib_soleus_r", "MN_soleus_r"))

print("\n== Ia motif ==")
print("Ia_soleus_r->MN:", edge("Ia_soleus_r", "MN_soleus_r"),
      " Ia->IaIN:", edge("Ia_soleus_r", "IaIN_soleus_r"),
      " PF_F1->IaIN gate:", edge("PF_F1_r", "IaIN_soleus_r"),
      " IaIN->ant MN(tib_ant):", edge("IaIN_soleus_r", "MN_tib_ant_r"),
      " IaIN<->IaIN(tib_ant):",
      edge("IaIN_soleus_r", "IaIN_tib_ant_r"))

print("\n== central afferent edges (Deng-audit-flagged) ==")
print("Ib_soleus_r->RG_E:", edge("Ib_soleus_r", "RG_E_r"),
      " Ib->InE:", edge("Ib_soleus_r", "InE_r"),
      " Ib->PF_E1:", edge("Ib_soleus_r", "PF_E1_r"))
print("Ia_iliacus_r->RG_F:", edge("Ia_iliacus_r", "RG_F_r"),
      " Ia->InF:", edge("Ia_iliacus_r", "InF_r"),
      " Ia->PF_F1:", edge("Ia_iliacus_r", "PF_F1_r"),
      " Ia->contra RG_F:", edge("Ia_iliacus_r", "RG_F_l"))
print("II_iliacus_r->RG_F:", edge("II_iliacus_r", "RG_F_r"))
print("II_glut_max1_r->RG_E (ii_e_central):",
      edge("II_glut_max1_r", "RG_E_r"))

print("\n== commissurals ==")
print("RG_F_r->CIN_F_r:", edge("RG_F_r", "CIN_F_r"),
      " CIN_F_r->RG_F_l:", edge("CIN_F_r", "RG_F_l"))
print("RG_E_r->CIN_E_r:", edge("RG_E_r", "CIN_E_r"),
      " CIN_E_r->InE_l:", edge("CIN_E_r", "InE_l"),
      " CIN_E_r->IBEXC_ankle_pf_l:",
      edge("CIN_E_r", "IBEXC_ankle_pf_l"))

print("\n== KINH (f1_anklepf_inh only) ==")
print("PF_F1_r->KINH_r:", edge("PF_F1_r", "KINH_r"),
      " KINH_r->MN_soleus_r:", edge("KINH_r", "MN_soleus_r"),
      " KINH_r->MN_vas_lat_r:", edge("KINH_r", "MN_vas_lat_r"),
      " HEEL_l->KINH_r (contra_kinh):", edge("HEEL_l", "KINH_r"))
print("renshaw key in winner:", "renshaw" in wp,
      "-> RC cells:", bool([k for k in net.idx if k.startswith("RC_")]))
