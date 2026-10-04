"""Audit probe 2: build syn6 (AARL_NET=syn6, default G) and verify the
wiring claims against the rule files + synergy basis.

Checks:
 1. counts + PF_S cell inventory + families/tau/df attributes.
 2. PF_S{k} -> MN conductance == Eq18(W[act,k]) EXACTLY for every
    covered muscle; edges only where W>0.
 3. pruned muscles: no PF_S in-edges; MN + Ia/II/Ib arc present.
 4. lamination: mixed channel excluded both directions; committed
    channels laminated with pf_recip_inh 4.0.
 5. RG drive split: mixed channel gets both RG_E*s and RG_F*(1-s);
    pure channels single-family at rg_to_pf 2.4.
 6. master-rules dress: heel -> InE/InF/PF_IN_E (ipsi), toe -> TOEDF
    -> df channel inh 2.749, toe->TOEDF exc 5.
 7. Shevtsova commissural crossed edges + gains.
 8. Shinohara autogenic motif gains per muscle.
 9. IBEXC wiring (rule gains; note absence of RG_E gate).
10. second build with ib_rge>0 to confirm the LBIN load dress is
    gain-gated (JSON-rule contract).
"""
import io
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
os.environ["AARL_NET"] = "syn6"

import build_network as bn  # noqa: E402
import muscle_map as MM
import params as P
from fsa_backsolve import analytical_conductance  # noqa: E402

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
print(f"actuator list: {len(acts)}")

net = bn.build(acts)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print(f"syn6 default counts = {counts}")

pop_names = [p["name"] for p in n.populations]
edges = defaultdict(list)
for c in n.connections:
    src, dst = pop_names[c["source"]], pop_names[c["destination"]]
    syn = c["params"]
    g = float(syn.get("max_conductance", np.nan))
    er = float(syn.get("reversal_potential", np.nan))
    edges[(src, dst)].append((g, "exc" if er > -1e-6 else "inh"))


def edge(src, dst):
    lst = edges.get((src, dst), [])
    return lst[0] if lst else None


print("\n== 1. PF cells / families ==")
print("families:", net.syn_families)
print("tau_mult:", net.syn_tau_mult)
print("df_channel:", net.syn_df_channel)
print("mixed:", net.syn_mixed)
pf_cells = [x for x in net.idx if x.startswith("PF_S")]
print("PF_S cells:", sorted(pf_cells))

print("\n== 2/3. Eq-18 PF->MN check (exact) ==")
basis = np.load(SPINAL / "synergy_basis.npz", allow_pickle=True)
from build_network_syn6 import PRUNE_MUSCLES  # noqa: E402
import runner  # noqa: E402
print("PRUNE sets equal (syn6 vs runner):",
      PRUNE_MUSCLES == runner.PRUNE_MUSCLES)
bad, ok_edges, pruned_pf, uncovered = [], 0, [], []
for side in ("r", "l"):
    names = [str(x) for x in basis[f"muscle_names_{side}"]]
    W = np.asarray(basis[f"W_{side}"], dtype=float)
    for act in [a for a in acts if a.endswith("_" + side)]:
        pf_in = [k for k in edges
                 if k[1] == f"MN_{act}" and k[0].startswith("PF_S")]
        if act in PRUNE_MUSCLES:
            if pf_in:
                pruned_pf.append(act)
            continue
        if act not in names:
            uncovered.append(act)
            if pf_in:
                pruned_pf.append(act)  # not pruned but outside basis
            continue
        row = names.index(act)
        k_vals = W[row, :]
        g_vals, valid = analytical_conductance(k_vals)
        for k in range(6):
            e = edge(f"PF_S{k + 1}_{side}", f"MN_{act}")
            if k_vals[k] > 0.0 and valid[k]:
                if e is None:
                    bad.append((act, k, "MISSING", float(g_vals[k])))
                elif not np.isclose(e[0], g_vals[k], rtol=0, atol=1e-12):
                    bad.append((act, k, "WRONG", e[0], float(g_vals[k])))
                else:
                    ok_edges += 1
            else:
                if e is not None:
                    bad.append((act, k, "EXTRA", e[0]))
print(f"Eq-18 PF->MN edges verified exact: {ok_edges}; mismatches: {bad}")
print(f"pruned/out-of-basis muscles with PF edges (should be []): "
      f"{pruned_pf}; out-of-basis non-pruned: {uncovered}")

print("\n== 4/5. lamination + RG drive ==")
fam = {s: net.syn_families[s] for s in net.syn_families}
mixed = {s: net.syn_mixed[s] for s in net.syn_mixed}
sf = {"r": None, l: None} if False else None
for side in ("r", "l"):
    for k in range(6):
        cell = f"PF_S{k + 1}_{side}"
        fam_k = fam[side][k]
        lam_in = edge(cell, f"PF_IN_{'E' if fam_k == 'E' else 'F'}_{side}")
        rg_e = edge(f"RG_E_{side}", cell)
        rg_f = edge(f"RG_F_{side}", cell)
        print(f"{cell} fam={fam_k} mixed={mixed[side][k]} "
              f"lam={'%s inh-side ok' if lam_in else 'ABSENT'} "
              f"RG_E={(rg_e[0] if rg_e else None)} "
              f"RG_F={(rg_f[0] if rg_f else None)}")
        # expected by the split rule
for side in ("r", "l"):
    st = {"r": [0.196, 0.953, 0.922, 0.246, 0.41, 0.89]}[side] if side == "r" \
        else [0.196, 0.953, 0.922, 0.246, 0.41, 0.89]
    for k in range(6):
        cell = f"PF_S{k + 1}_{side}"
        s_frac = st[k]
        e = edge(f"RG_E_{side}", cell)
        f = edge(f"RG_F_{side}", cell)
        if 0.35 < s_frac < 0.65:
            exp_e, exp_f = 2.4 * s_frac, 2.4 * (1 - s_frac)
            good = (e and f and np.isclose(e[0], exp_e, atol=1e-9)
                    and np.isclose(f[0], exp_f, atol=1e-9))
            print(f"{cell} MIXED split: got E={e[0] if e else None} "
                  f"F={f[0] if f else None} expected "
                  f"E={exp_e:.4f} F={exp_f:.4f} -> "
                  f"{'OK' if good else 'MISMATCH'}")
        else:
            src = "E" if s_frac >= 0.5 else "F"
            good = (e and not f) if src == "E" else (f and not e)
            print(f"{cell} pure {src}: "
                  f"single-family={'OK' if good else 'MISMATCH'}")

print("\n== 6. master-rules dress (heel/toe) ==")
for side in ("r", "l"):
    print(f"[{side}] heel->InE:", edge(f"HEEL_{side}", f"InE_{side}"),
          " heel->InF:", edge(f"HEEL_{side}", f"InF_{side}"),
          " heel->PF_IN_E:", edge(f"HEEL_{side}", f"PF_IN_E_{side}"))
    print(f"[{side}] toe->TOEDF:", edge(f"TOE_{side}", f"TOEDF_{side}"),
          " TOEDF->PF_S5:", edge(f"TOEDF_{side}", f"PF_S5_{side}"))
    # contra-side check: does heel_r hit the CONTRA InF? (drawing says yes)
    other = "l" if side == "r" else "r"
    print(f"[{side}] heel->contra InF ({other}):",
          edge(f"HEEL_{side}", f"InF_{other}"), " <- drawing wants exc 0.5")

print("\n== 7. commissurals ==")
for a, b in (("l", "r"), ("r", "l")):
    print(f"{a}->{b}: V0V->InE:", edge(f"V0V_{a}", f"InE_{b}"),
          " V0D->RG_F:", edge(f"V0D_{a}", f"RG_F_{b}"),
          " V3E->RG_E:", edge(f"V3E_{a}", f"RG_E_{b}"),
          " V3E->InE:", edge(f"V3E_{a}", f"InE_{b}"))
    print(f"   drive edges: RG_F->V2a:", edge(f"RG_F_{a}", f"V2a_{a}"),
          " RG_F->V0D:", edge(f"RG_F_{a}", f"V0D_{a}"),
          " RG_E->V3E:", edge(f"RG_E_{a}", f"V3E_{a}"),
          " V2a->V0V:", edge(f"V2a_{a}", f"V0V_{a}"))

print("\n== 8. Shinohara motif (sample muscles both families) ==")
for act in ("vas_lat_r", "soleus_r", "iliopsoas_r", "semimem_r",
            "tib_ant_r", "quad_fem_r"):
    mn = f"MN_{act}"
    row = dict(
        ia_mn=edge(f"Ia_{act}", mn),
        ia_iain=edge(f"Ia_{act}", f"IaIN_{act}"),
        ii_iix=edge(f"II_{act}", f"IIX_{act}"),
        iix_mn=edge(f"IIX_{act}", mn),
        ii_iiin=edge(f"II_{act}", f"IIIN_{act}"),
        ib_ibin=edge(f"Ib_{act}", f"IBIN_{act}"),
        ibin_mn=edge(f"IBIN_{act}", mn),
    )
    print(act, {k: (v if v else None) for k, v in row.items()})

print("\n== 9. IBEXC (extensor groups) ==")
for grp_side in sorted({k for k in net.idx if k.startswith("IBEXC_")}):
    pass
ibexcs = sorted({s.split("_" + s[-1])[0] + "_" + s[-1]
                 for s in net.idx if s.startswith("IBEXC_")})
for name in sorted(k for k in net.idx if k.startswith("IBEXC_")):
    side = name[-1]
    grp = name[:-2]
    src = [s for (s, d) in edges if d == name]
    print(f"{name}: in-edges {[(s, edges[(s, name)][0]) for s in src]}")
print("IBEXC gets RG_E gate?",
      bool([1 for (s, d) in edges if d.startswith("IBEXC_")
            and s.startswith("RG_E_")]))

print("\n== 10. LBIN dress gating (default ib_rge=0 -> absent) ==")
print("default: LBIN out-edges:",
      [(d, edges[(s, d)][0]) for (s, d) in edges if s.startswith("LBIN_")])
print("default: Ib->LBIN edges:",
      bool([(s, d) for (s, d) in edges if d.startswith("LBIN_")]))

print("\n== runner contract flags ==")
print("stance_fb:", net.stance_fb, "aff_loops:", net.aff_loops,
      "vest:", net.vest, "renshaw:", net.renshaw,
      "brainstem:", getattr(net, "brainstem", None))
print("G defaults: heel_rge", P.G["heel_rge"], "toe_rge", P.G["toe_rge"],
      "ib_rge", P.G["ib_rge"], "syn6_brainstem", P.G["syn6_brainstem"])
print("inputs sample:", net.inputs[:12], "...", net.inputs[-6:])

print("\n== gain-gated rebuild: ib_rge=0.5 ==")
P.G["ib_rge"] = 0.5
net2 = bn.build(acts)
n2 = net2.net
p2 = [p["name"] for p in n2.populations]
e2 = defaultdict(list)
for c in n2.connections:
    e2[(p2[c["source"]], p2[c["destination"]])].append(
        (float(c["params"].get("max_conductance", np.nan)),))
lbin_out = [(d, e2[(s, d)][0][0]) for (s, d) in e2 if s.startswith("LBIN_")]
print("ib_rge=0.5 LBIN out-edges:", lbin_out[:8], "... total",
      len(lbin_out))
ib_lbin = [(s, d) for (s, d) in e2 if d.startswith("LBIN_")]
print("Ib->LBIN edge count:", len(ib_lbin))
