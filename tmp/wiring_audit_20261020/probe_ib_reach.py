"""Audit probe 2026-10-20: where does extensor Ib reach, vs where does
HEEL contact reach, in each network family (s3k stock, syn6)?

One process per variant/overlay (params.G is module-global):

    D:/Anaconda/envs/myo/python.exe probe_ib_reach.py [overlay]

    overlay: ""            -> defaults
             "ib_rge05"    -> G["ib_rge"] = 0.5        (syn6 LBIN dress ON)
             "s3kwinner"   -> the s3k trial-34 stage-4 winner overlay

Prints census + every edge leaving {Ib-aff, HEEL_IN, TOE_IN, LBIN, IBEXC,
IBIN} grouped by (src_kind, dst_kind, sign) with counts and conductance
set, plus (syn6) the family/mixed/df-channel assignment and the
per-extensor-group channel carriage.
"""
from __future__ import annotations

import io
import os
import sys
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
spinal = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, spinal)
os.chdir(spinal)

import params as P                     # noqa: E402
import muscle_map as MM                # noqa: E402
import build_network as BN             # noqa: E402

overlay = sys.argv[1] if len(sys.argv) > 1 else ""
variant = os.environ.get("AARL_NET", "").strip().lower()

if overlay == "ib_rge05":
    P.G["ib_rge"] = 0.5
elif overlay == "ibantag05":
    P.G["ib_e_central"] = 0.5
    P.G["iain_pf_gate"] = 0.5
    P.G["pf_mixed_recip"] = 0.5
elif overlay == "s3kwinner":
    w = dict(full_rules=1.0, ia_in=0.8250460785281135,
             heel_rge=0.37530421840765804, toe_rge=0.14572281192093522,
             ib_rge=0.6741474252269108, ib_e_central=0.24354654619098992,
             ia_f_central=0.3985916558720821,
             ii_f_central=0.38389852925373125,
             ii_e_central=0.09300453229764793,
             c1_gain=0.4978970446353496, v3_gain=0.09483003906183418,
             v3_to_ibexc=0.8784035043219222,
             ia_f_contra_f=0.797446189481362,
             f1_anklepf_inh=0.9285368580077491,
             contra_kinh=1.4251105067452698,
             contact_onset=0.8138186176322883)
    P.G.update(w)

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")

net = BN.build(acts)
n = net.net
counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
          n.get_num_connections())
print(f"variant={variant or 'standard'} overlay={overlay or 'defaults'} "
      f"census = {counts}")

pop_names = [p["name"] for p in n.populations]

SRC_KINDS = ("Ib_", "HEEL_", "TOE_", "LBIN_", "IBEXC_", "IBIN",
             "II_", "Ia_", "IaIN", "IIIN")


def dst_kind(dst: str) -> str:
    for pre, k in (
            ("RG_E", "RG_E"), ("RG_F", "RG_F"), ("InE", "InE"),
            ("InF", "InF"), ("PF_IN", "PF_IN"), ("PF_S", "PF_S(p6)"),
            ("PF_ANK", "PF_ANK"), ("PF_HIP", "PF_jpf)"), ("PF_KNEE", "PF_jpf)"),
            ("PF_E", "PF_E)"), ("PF_F", "PF_F)"),
            ("HEEL", "HEEL_IN"), ("TOE", "TOE_IN"), ("LBIN", "LBIN"),
            ("IBEXC", "IBEXC"), ("IBIN", "IBIN"), ("IaIN", "IaIN"),
            ("IIIN", "IIIN"), ("IIX", "IIX"), ("KINH", "KINH"),
            ("V0", "V0"), ("V2", "V2"), ("V3", "V3"), ("CIN", "CIN"),
            ("MN_", "MN"), ("Ia_", "Ia-aff"), ("II_", "II-aff"),
            ("Ib_", "Ib-aff"), ("TOEDF", "TOEDF")):
        if dst.startswith(pre):
            return k
    return dst


groups = defaultdict(list)
for c in n.connections:
    src = pop_names[c["source"]]
    dst = pop_names[c["destination"]]
    if not src.startswith(SRC_KINDS):
        continue
    syn = c["params"]
    g = float(syn.get("max_conductance", float("nan")))
    sign = "exc" if float(syn.get("reversal_potential", 0.0)) > -1e-6 else "inh"
    groups[(src.split("_")[0], dst_kind(dst), sign)].append((src, dst, g))

print("--- edges leaving Ib / HEEL / TOE / LBIN / IBEXC / IBIN / II / Ia ---")
for (sk, dk, sign) in sorted(groups):
    edges = groups[(sk, dk, sign)]
    gs = sorted({round(e[2], 4) for e in edges})
    print(f"  {sk:>6} -> {dk:<12} {sign:>3}  n={len(edges):<5} g={gs}")

# --- syn6-specific: family/mixed/df assignment + extensor carriage -------
if variant == "syn6":
    tabs = net  # Syn6Network carries the evidence attrs
    for side in ("r", "l"):
        fam = tabs.syn_families[side]
        mixed = tabs.syn_mixed[side]
        print(f"side {side}: family={fam} mixed={mixed} "
              f"df_channel={[None if tabs.syn_df_channel[side] is None else tabs.syn_df_channel[side] + 1][0]} "
              f"tau={[round(t, 3) for t in tabs.syn_tau_mult[side]]}")
    # which E channels carry each extensor stance group (W > 0)
    from build_network_syn6 import _load_synergy_tables
    t = _load_synergy_tables()
    for side in ("r",):
        names, W, fam = t[side]["names"], t[side]["W"], t[side]["family"]
        carriage = defaultdict(set)
        for i, act in enumerate(names):
            info = MM.classify(act)
            if info is None or info.groups[0] not in MM.EXTENSOR_STANCE_GROUPS:
                continue
            for k in range(len(fam)):
                if W[i, k] > 0 and fam[k] == "E":
                    carriage[info.groups[0]].add(k + 1)
        print(f"side {side} extensor-group carriage in E channels: "
              f"{dict(sorted(carriage.items()))}")
        # antagonist pairs across families (group -> set of channels)
        chan = defaultdict(set)
        for i, act in enumerate(names):
            info = MM.classify(act)
            if info is None:
                continue
            for k in range(len(fam)):
                if W[i, k] > 0:
                    chan[info.groups[0]].add(k + 1)
        print(f"side {side} group->channels: "
              f"{dict(sorted(chan.items()))}")
print("done.")
