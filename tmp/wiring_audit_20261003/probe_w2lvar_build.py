"""Audit probe 4: build w2lvar (AARL_NET=w2lvar) and verify wiring
against Ben's rule files (master rules + shevtsova + shinohara)."""
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
os.environ["AARL_NET"] = "w2lvar"

import muscle_map as MM  # noqa: E402
import build_network as bn  # noqa: E402
from build_network_w2lvar import GAINS, VARIANT_G, W2L_CELL_ROWS  # noqa: E402

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")

net = bn.build(acts)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print(f"w2lvar counts (with VARIANT_G applied) = {counts}")

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


import params as P  # noqa: E402
print("VARIANT_G applied to params.G:",
      {k: P.G[k] for k in VARIANT_G})

print("\n== PF cells + lamination ==")
pf_cells = sorted(x for x in net.idx if x.startswith("PF_")
                  and not x.startswith("PF_IN"))
print("PF cells:", pf_cells)
for side in ("r", "l"):
    for pair, (e_hc, f_hc) in (("HIP", ("HIP-E", "HIP-F")),
                               ("KNEE", ("KNEE-E", "KNEE-F"))):
        ie, iff = f"PF_IN_{pair}-E_{side}", f"PF_IN_{pair}-F_{side}"
        print(f"[{side}] {pair}: E->IN_E", edge(f"PF_{e_hc}_{side}", ie),
              " IN_E->F", edge(ie, f"PF_{f_hc}_{side}"),
              " F->IN_F", edge(f"PF_{f_hc}_{side}", iff),
              " IN_F->E", edge(iff, f"PF_{e_hc}_{side}"))
    print(f"[{side}] RG_E->HIP-E", edge(f"RG_E_{side}", f"PF_HIP-E_{side}"),
          " RG_F->HIP-F", edge(f"RG_F_{side}", f"PF_HIP-F_{side}"),
          " RG_E->KNEE-E", edge(f"RG_E_{side}", f"PF_KNEE-E_{side}"),
          " RG_F->KNEE-F", edge(f"RG_F_{side}", f"PF_KNEE-F_{side}"))

print("\n== PF->MN weights vs joint_pf_weights.json ==")
jw = json.loads((SPINAL / "joint_pf_weights.json").read_text(
    encoding="utf-8"))["w"]
print("joint_pf_weights rows:", sorted(jw.keys()))
for act in ("vas_lat_r", "soleus_r", "glut_max1_r", "iliacus_r",
            "semimem_r", "med_gas_r", "tib_ant_r", "ercspn_r"):
    base = act.rsplit("_", 1)[0]
    got = {}
    for cell in ("HIP-E", "HIP-F", "KNEE-E", "KNEE-F"):
        got[cell] = edge(f"PF_{cell}_{act[-1]}", f"MN_{act}")
    print(act, {k: (v[0] if v else None) for k, v in got.items()})

print("\n== heel / toe dress ==")
for side in ("r", "l"):
    other = "l" if side == "r" else "r"
    print(f"[{side}] heel->ipsi InE:", edge(f"HEEL_{side}", f"InE_{side}"),
          " heel->contra InF:", edge(f"HEEL_{side}", f"InF_{other}"),
          " heel->PF_IN_HIP-E:", edge(f"HEEL_{side}", f"PF_IN_HIP-E_{side}"),
          " heel->PF_IN_KNEE-E:",
          edge(f"HEEL_{side}", f"PF_IN_KNEE-E_{side}"))
    print(f"[{side}] toe->TOEDF:", edge(f"TOE_{side}", f"TOEDF_{side}"),
          " TOEDF->PF_KNEE-F:",
          edge(f"TOEDF_{side}", f"PF_KNEE-F_{side}"))
    print(f"[{side}] LBIN->RG_E:", edge(f"LBIN_{side}", f"RG_E_{side}"))

print("\n== commissurals ==")
for a, b in (("r", "l"), ("l", "r")):
    print(f"{a}->{b}: RG_F->V2a", edge(f"RG_F_{a}", f"V2a_{a}"),
          " V2a->V0V", edge(f"V2a_{a}", f"V0V_{a}"),
          " V0V->INI", edge(f"V0V_{a}", f"INI_{b}"),
          " INI->RG_F", edge(f"INI_{b}", f"RG_F_{b}"),
          " RG_F->V0D", edge(f"RG_F_{a}", f"V0D_{a}"),
          " V0D->RG_F", edge(f"V0D_{a}", f"RG_F_{b}"),
          " RG_E->V3E", edge(f"RG_E_{a}", f"V3E_{a}"),
          " V3E->RG_E", edge(f"V3E_{a}", f"RG_E_{b}"),
          " V3E->InE", edge(f"V3E_{a}", f"InE_{b}"))
print("DRIVE->V0V_r:", edge("DRIVE", "V0V_r"),
      " DRIVE->V0D_l:", edge("DRIVE", "V0D_l"))

print("\n== autogenic central projections ==")
nE = {s: sum(1 for a in acts if a.endswith("_" + s)
             and MM.classify(a).groups[0] in MM.EXTENSOR_STANCE_GROUPS)
      for s in ("r", "l")}
nF = {s: sum(1 for a in acts if a.endswith("_" + s)
             and MM.classify(a).groups[0] in ("hip_flex", "knee_flex",
                                              "ankle_df", "hip_add",
                                              "trunk_flex"))
      for s in ("r", "l")}
print("net._nE:", net._nE, "net._nF:", net._nF, "recount:", nE, nF)
for act in ("soleus_r", "glut_max1_r", "vas_lat_r"):
    s = act[-1]
    g_exp = GAINS["w2l_ib_central"] / max(nE[s], 1)
    tgts = [edge(f"Ib_{act}", t) for t in
            (f"PF_HIP-E_{s}", f"PF_KNEE-E_{s}", f"RG_E_{s}", f"InE_{s}")]
    print(act, "Ib->E-trio:", tgts, "expected g:", round(g_exp, 5))
for act in ("iliacus_r", "tib_ant_r", "semiten_r"):
    s = act[-1]
    g_exp = GAINS["w2l_ia_central"] / max(nF[s], 1)
    tgts = [edge(f"Ia_{act}", t) for t in
            (f"PF_HIP-F_{s}", f"PF_KNEE-F_{s}", f"RG_F_{s}", f"InF_{s}")]
    print(act, "Ia->F-trio:", tgts, "expected g:", round(g_exp, 5))

print("\n== Renshaw / IaIN gate / IBEXC ==")
print("MN->RC (vas_lat_r):", edge("MN_vas_lat_r", "RC_vas_lat_r"),
      " RC->MN:", edge("RC_vas_lat_r", "MN_vas_lat_r"),
      " RC->IaIN:", edge("RC_vas_lat_r", "IaIN_vas_lat_r"))
n_rc_mut = sum(1 for (s, d) in edges if s.startswith("RC_")
               and d.startswith("RC_") and s != d)
print("RC<->RC edges:", n_rc_mut)
print("IaIN gate (vas_lat_r knee_ext->KNEE-F):",
      edge("PF_KNEE-F_r", "IaIN_vas_lat_r"))
print("IaIN gate (glut_max1_r hip_ext->HIP-F):",
      edge("PF_HIP-F_r", "IaIN_glut_max1_r"))
print("IBEXC knee_ext_r: RG_E gate:",
      edge("RG_E_r", "IBEXC_knee_ext_r"),
      " Ib->IBEXC:", edge("Ib_vas_lat_r", "IBEXC_knee_ext_r"),
      " IBEXC->MN:", edge("IBEXC_knee_ext_r", "MN_vas_lat_r"))
print("KINH present?", [k for k in net.idx if k.startswith("KINH")])

print("\n== aliases (logging only) ==")
print("idx PF_ANK-F_r ->", net.idx.get("PF_ANK-F_r"),
      "== PF_KNEE-F_r", net.idx.get("PF_KNEE-F_r"))
alias_syn = [k for k in edges if k[0].startswith("PF_ANK")]
print("edges sourced on alias names (should be []):", alias_syn)
print("true neuron count check:", net._n_true, "n populations:",
      n.get_num_neurons())
