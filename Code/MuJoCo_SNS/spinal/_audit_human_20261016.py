"""Connectome audit driver (human-body walker nets) - 2026-10-16 campaign.

Builds the DEFAULT net (counts recorded for the gate) and the S3K WINNER
net (curr_s3k_nocross trial 34 = reports_20260923/s3k_trial34_full_params.json,
applied via _curriculum.set_stage(4, ...) with the exact merge of
_replay_best_20261003.py's s3k branch), dumps EVERY compiled connection,
and classifies the rule-relevant edge families for Ben's connectome rules
(a) per-PF-layer RG+PF feedback, (b) heel=PF-layer stance reset / toe->
TOEDF dorsiflexion-inhibition only, (c) autogenic per-muscle proprio-
ception, (d) Ib load handling (PF-E 0.5 / RG-E gate 0.1), (e) afferents
wired directly into RG half-centers. Also re-dumps Ben's three canonical
exports (ben_rules/ben_shevtsova/ben_shinohara 20260924.json) so the
cross-check evidence lives beside this run. Read-only: builds in memory.
"""
import io
import json
import sys
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import os

import mujoco

import build_network as bn
import params as P
from muscle_map import classify

MODEL = Path = None
import pathlib
_here = pathlib.Path(__file__).parent
MODEL = pathlib.Path(os.environ.get(
    "AARL_MODEL",
    r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
    r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"))


def acts_from_model():
    model = mujoco.MjModel.from_xml_path(str(MODEL))
    return [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]


def kind(name: str) -> str:
    if name.startswith(("DRIVE", "POSTURE", "BAL")):
        return "desc-in"
    for pre, k in (("GAMMA_", "gamma-in"), ("ALPHA_", "alpha-in"),
                   ("VEST_", "vest-in"),
                   ("HEEL", "HEEL_IN"), ("TOEDF", "TOEDF_IN"),
                   ("TOE", "TOE_IN"), ("LBIN", "LBIN_IN"),
                   ("IaIN", "IaIN_IN"), ("IBEXC", "IBEXC_IN"),
                   ("IIX_", "IIX_IN"), ("IIIN", "IIIN_IN"),
                   ("IBIN", "IBIN_IN"), ("KINH", "KINH_IN"),
                   ("RC_", "RC"), ("ADAP", "ADAP_IN"),
                   ("InE_", "InE"), ("InF_", "InF"),
                   ("RG_E_", "RG-E"), ("RG_F_", "RG-F"),
                   ("CIN_", "CIN"),
                   ("PF_", "PF"), ("MN_", "MN"),
                   ("Ia_", "Ia-aff"), ("II_", "II-aff"), ("Ib_", "Ib-aff")):
        if name.startswith(pre):
            return k
    return name


def dump_net(net, title):
    n = net.net
    pop = [p["name"] for p in n.populations]
    counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
              n.get_num_connections())
    print(f"\n{'='*72}\n{title}\n{'='*72}")
    print(f"counts neurons/inputs/synapses = {counts}")
    print(f"inputs ({len(net.inputs)}): {net.inputs}")
    groups = defaultdict(list)
    for c in n.connections:
        src, dst = pop[c["source"]], pop[c["destination"]]
        g = float(c["params"]["max_conductance"])
        er = float(c["params"]["reversal_potential"])
        groups[(kind(src), kind(dst), "exc" if er > -1e-6 else "inh")].append(
            (src, dst, g))
    print(f"total connections: {sum(len(v) for v in groups.values())}")
    for (sk, dk, sign), lst in sorted(groups.items()):
        gs = ", ".join(f"{g:.4g}" for g in sorted({v[2] for v in lst})[:4])
        print(f"  {sign.upper():3s} {sk:10s} -> {dk:10s} x{len(lst):4d}"
              f"  g[{gs}]")

    # ---- rule-relevant detail dumps ----
    edges = []
    for c in n.connections:
        src, dst = pop[c["source"]], pop[c["destination"]]
        edges.append((src, dst,
                      "exc" if float(c["params"]["reversal_potential"]) > -1e-6
                      else "inh",
                      float(c["params"]["max_conductance"])))

    def show(pred, header):
        sel = [e for e in edges if pred(e)]
        print(f"\n--- {header} ({len(sel)} edges) ---")
        agg = defaultdict(list)
        for src, dst, s, g in sel:
            agg[(src.split('_')[0] if src[2] != '_' else src,
                 kind(src), kind(dst), s, round(g, 4))].append((src, dst))
        for (raw, sk, dk, s, g), lst in sorted(agg.items()):
            ex = ", ".join(f"{a}->{b}" for a, b in lst[:4])
            more = f" (+{len(lst)-4} more)" if len(lst) > 4 else ""
            print(f"  {s} g={g:<9g} {sk}->{dk} x{len(lst)}: {ex}{more}")

    # (a) RG -> PF topology + PF inventory
    pfn = sorted(x for x in net.idx if x.startswith("PF_"))
    print(f"\n--- (a) PF neurons ({len(pfn)}): {pfn}")
    show(lambda e: kind(e[0]) in ("RG-E", "RG-F") and kind(e[1]) == "PF",
         "(a) RG -> PF drive edges")
    show(lambda e: kind(e[0]) == "PF" and kind(e[1]) == "MN",
         "(a) PF -> MN fan-out (per-muscle weights)")

    # (b) heel / toe / TOEDF
    show(lambda e: e[0].startswith(("HEEL_", "TOE_", "TOEDF_")),
         "(b) HEEL/TOE/TOEDF-sourced edges")
    show(lambda e: e[1].startswith("TOEDF_"), "(b) edges INTO TOEDF")

    # (c) autogenic check: afferent -> MN edges (own vs other muscle)
    own = cross = 0
    cross_ex = []
    for src, dst, s, g in edges:
        if kind(src) in ("Ia-aff", "II-aff", "Ib-aff") and kind(dst) == "MN":
            m_src = src.split("_", 1)[1]
            m_dst = dst.split("_", 1)[1]
            if m_src == m_dst:
                own += 1
            else:
                cross += 1
                cross_ex.append((src, dst, s, g))
    print(f"\n--- (c) afferent->MN direct edges: autogenic {own}, "
          f"NON-autogenic {cross} {cross_ex[:6]}")
    show(lambda e: kind(e[0]) in ("Ia-aff",) and kind(e[1]) in
         ("IaIN_IN",), "(c) Ia -> IaIN")
    show(lambda e: kind(e[0]) in ("II-aff",) and kind(e[1]) in
         ("IIX_IN", "IIIN_IN"), "(c) II -> IIX/IIIN")
    show(lambda e: kind(e[0]) in ("Ib-aff",) and kind(e[1]) in
         ("IBIN_IN", "IBEXC_IN"), "(c) Ib -> IBIN/IBEXC")

    # (d) Ib load handling: Ib afferent, IBEXC, LBIN, and RG->IBEXC gate
    show(lambda e: kind(e[0]) == "LBIN_IN", "(d) LBIN-sourced edges")
    show(lambda e: kind(e[1]) == "IBEXC_IN" and kind(e[0]) in
         ("RG-E", "RG-F", "desc-in"), "(d) RG/desc -> IBEXC gates")
    show(lambda e: kind(e[0]) == "IBEXC_IN", "(d) IBEXC-sourced edges")
    show(lambda e: kind(e[0]) == "Ib-aff" and kind(e[1]) in
         ("PF", "RG-E", "RG-F", "InE", "InF"),
         "(d)+(e) Ib -> PF/RG/In central projections")

    # (e) afferents INTO RG half-centers (direct) and into InE/InF
    show(lambda e: kind(e[1]) in ("RG-E", "RG-F") and kind(e[0]) in
         ("Ia-aff", "II-aff", "Ib-aff", "HEEL_IN", "TOE_IN", "LBIN_IN",
          "IaIN_IN", "IIX_IN", "IIIN_IN", "IBIN_IN", "IBEXC_IN"),
         "(e) afferent/IN -> RG half-centers (DIRECT)")
    show(lambda e: kind(e[1]) in ("InE", "InF") and kind(e[0]) in
         ("Ia-aff", "II-aff", "Ib-aff", "HEEL_IN", "TOE_IN", "LBIN_IN"),
         "(e/b) afferent/IN -> InE/InF lamination")

    # drive map for the L/R plan
    show(lambda e: e[0] in ("DRIVE", "POSTURE"), "drive map: DRIVE/POSTURE"
         " out-edges")
    return counts


def ben_exports():
    for f in ("ben_rules_20260924.json", "ben_shevtsova_20260924.json",
              "ben_shinohara_20260924.json"):
        d = json.load(open(f, encoding="utf-8"))
        print(f"\n{'='*72}\nBEN EXPORT {f}: {len(d['nodes'])} nodes, "
              f"{len(d['edges'])} edges\n{'='*72}")
        for e in d["edges"]:
            print(f"  {e['sign']} g={e['gain']:<6g} {e['from'][:46]:46s}"
                  f" -> {e['to'][:46]} [{e.get('tag', '')}]")


def main():
    acts = acts_from_model()
    print(f"model actuators: {len(acts)} ({acts[0]}..{acts[-1]})")
    unclass = [a for a in acts if classify(a) is None]
    print(f"unclassifiable: {unclass}")

    # ---- DEFAULT build (module-default G; canonical counts) ----
    net0 = bn.build(acts, interleg=True)
    dump_net(net0, "DEFAULT BUILD (all conditional gains at params.G "
                   "module defaults)")

    # ---- S3K WINNER build (production config) ----
    import _curriculum as C
    base = dict(json.loads(
        open("best_walk_params_v10.json", encoding="utf-8").read())
        ["multipliers"])
    base["renshaw"] = 0.5
    s3k = json.loads(open(
        "reports_20260923/s3k_trial34_full_params.json",
        encoding="utf-8").read())["params"]
    base.update(s3k)
    C.set_stage(4, base)
    print("\nwinner-config G keys of interest: "
          + ", ".join(f"{k}={P.G[k]:.4g}" for k in (
              "full_rules", "joint_pf", "heel_rge", "toe_rge", "ib_rge",
              "ib_e_central", "ia_f_central", "ii_f_central", "ii_e_central",
              "ia_f_contra_f", "c1_gain", "v3_gain", "ia_in", "renshaw",
              "contra_kinh", "heel_pf_layer", "toe_df_inh", "heel_in_f_exc",
              "ia_pf_f", "ii_pf_f", "rg_to_pf", "descend_to_rg_e",
              "descend_to_rg_f", "ib_group_exc", "ib_exc_to_mn")))
    net1 = bn.build(acts, interleg=True)
    dump_net(net1, "S3K WINNER BUILD (curr_s3k_nocross t34 via set_stage(4))"
                   " - the production human-body walker")

    ben_exports()
    print("\nAUDIT DRIVER DONE")


if __name__ == "__main__":
    main()
