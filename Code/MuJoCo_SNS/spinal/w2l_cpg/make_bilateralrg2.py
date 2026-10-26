"""Generate the "bilateralrg2" connectome template from Ben's CORRECTED
Walker_2_Layer_CPG_BilateralRG.aproj (2026-10-18, "physical contact sensing
closed (Li-style chain)" + "model updates").

WHAT CHANGED vs the "bilateralrg" template (2026-09-26, verified by this
script's delta report at generation time):
  1. Toe-contact cascade REMOVED (toe no longer drives RG ext / Hip-Knee PF
     ext / Hip-Knee MN ext). Toe now drives ONLY the ipsilateral ankle
     extensor MN (Toe_Ipsi_AnkMNext_Excite).
  2. Heel-contact cascade re-typed into the WEAK afferent-HC family
     (dedicated types Heel_Ipsi_*_Excite, gmax 1e-8 S = the build_aff.pl
     "Afferent HC Excite SynAmp 0.01" class) - contact now modulates the
     extensor layers through the same weak pathway as Ia/II/Ib, instead of
     the strong generic class.
  3. NEW crossed contact pathway: heel -> CONTRAlateral RG flx (exc,
     Heel_Contra_RGflx_Excite, weak class) - the crossed flexor-excitatory
     edge family (DESIGN.md 2026-09-22: "the W2L lesson... each heel strike
     KICKS the other leg into flexion").
  4. NEW heel -> ipsilateral Ank MN ext (Heel_Ipsi_AnkMNext_Excite, GENERIC
     5e-7 class - the only strong contact edge left).
  5. Li-style per-joint afferent chain confirmed/transcribed with dedicated
     types (Ia_HipFlx/Ia_KneeFlx -> {PF flx, RG flx}; II likewise;
     Ib_HipExt/Ib_KneeExt -> {PF ext, RG ext}; all weak class).
  6. Knee-extensor Renshaw LATERALIZED: the aproj carries TWO nodes labeled
     "L Knee MN ext RE" (GUIDs a3828afc / efa2eb22). GUID-level tracing
     shows a3828afc serves the R side (R Knee MN ext -> a3828afc ->
     R Knee MN ext / R Knee ext Ia / R Knee MN flx RE) - it is Ben's
     R-knee-ext RE with a stale L label (his 2026-10-03 editor export
     renamed it RC_8971). This generator relabels it "R Knee MN ext RE"
     (LABEL fix only; every edge keeps its GUID-resolved endpoint).
  7. One flagged edge: "R Hip MN flx RE -> R Hip MN ext RE" carries the
     V3-Commissural-EXCITE synapse type (a GUI clone-type artifact, the
     family documented in w2l_cpg README); Ben's 2026-10-03 export rules it
     inh/rc_mutual. Transcribed rc_mutual INH per Ben's export; flagged for
     his confirmation.

EFFECTIVE-STRENGTH DOCTRINE (unchanged, build_w2l_net.py): the AnimatLab
synapse TYPE's conductance is the effective strength (the .pl build chain
creates dedicated types precisely to set them; per-instance values are GUI
prefill). Template tags therefore encode the strength FAMILY; build_w2l_net
maps families to proven stack conductances.

Usage (from spinal/w2l_cpg):
    python make_bilateralrg2.py            # parse, diff, write template key
    python make_bilateralrg2.py --check    # re-run the diff only (no write)
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = Path(__file__).parent
SPINAL = HERE.parent
APROJ = (SPINAL.parents[2] / "Neuromechanical_Models" / "Walker_2_Layer_CPG_BilateralRG"
         / "Walker_2_Layer_CPG_BilateralRG.aproj").resolve()
MINE_OUT = SPINAL / "w2l_aproj_mine_20261018.json"
TEMPLATES = SPINAL / "connectome_templates.json"

# The R-knee-ext Renshaw node (GUID-prefixed key: full GUID in mine json).
RE_KNEE_EXT_R_GUID_PREFIX = "a3828afc"


def _val(el, tag):
    ch = el.find(tag)
    if ch is None:
        return None
    for k in ("Actual", "Value"):
        if k in ch.attrib:
            return float(ch.attrib[k])
    return ch.text


def _txt(el, tag, default=None):
    ch = el.find(tag)
    return (ch.text if ch is not None and ch.text else default)


def mine_aproj(path: Path) -> dict:
    """Parse nodes + synapse instances + synapse types from the aproj."""
    root = ET.parse(path).getroot()

    syn_types = {}
    for st in root.iter("SynapseTypes"):
        for lnk in st:
            syn_types[_txt(lnk, "ID")] = {
                "name": _txt(lnk, "Name"),
                "class": (_txt(lnk, "ClassName") or "").split(".")[-1],
                "eref": _val(lnk, "EquilibriumPotential"),
                "gmax": _val(lnk, "MaxSynapticConductance"),
            }

    nodes, edges = {}, []
    for node in root.iter("Node"):
        cls = (_txt(node, "ClassName") or "")
        if "NeuralModule" in cls:
            continue
        gid = _txt(node, "ID")
        nodes[gid] = {
            "guid": gid,
            "name": _txt(node, "Text") or _txt(node, "Name"),
            "class": cls.split(".")[-1],
            "tonic": _val(node, "TonicStimulus"),
        }

    for lk in root.iter("Link"):
        if not (_txt(lk, "ClassName") or "").endswith("Synapse"):
            continue
        t = syn_types.get(_txt(lk, "SynapticTypeID"), {})
        # PER-INSTANCE values (2026-10-25 ruling pt4: mines record
        # per-instance SynAmp per edge). The instance's own
        # SynapticConductance is the RUNTIME strength (Li-model-verified);
        # the type's gmax is creation-time prefill only.
        inst_sc = lk.find("SynapticConductance")
        inst_cd = lk.find("ConductionDelay")
        edges.append({
            "src_guid": _txt(lk, "OriginID"),
            "dst_guid": _txt(lk, "DestinationID"),
            "type": t.get("name"),
            "eref": t.get("eref"),
            "gmax": t.get("gmax"),
            "synamp_instance": (inst_sc.attrib.get("Actual")
                                if inst_sc is not None else None),
            "delay_ms": (inst_cd.attrib.get("Actual")
                         if inst_cd is not None else None),
        })
    return {"nodes": list(nodes.values()), "edges": edges,
            "syn_types": {v["name"]: v for v in syn_types.values()}}


# ---------------------------------------------------------------- node typing
def node_type(name: str) -> str:
    """aproj neuron label -> template node type (bilateralrg taxonomy)."""
    side, rest = name.split(" ", 1)
    if rest in ("RG ext", "RG flx"):
        return "HC-RG-E" if rest == "RG ext" else "HC-RG-F"
    if rest in ("RG ext IN", "RG flx IN"):
        return "IN-InE" if rest == "RG ext IN" else "IN-InF"
    if rest == "c1":
        return "IN-C"
    if rest == "V3":
        return "IN-V3"
    if rest.endswith("PF ext") or rest.endswith("PF flx"):
        return "HC-PF-E" if rest.endswith("PF ext") else "HC-PF-F"
    if rest.endswith("PF ext IN") or rest.endswith("PF flx IN"):
        return "IN-PF"
    if rest.endswith(" IN") and ("Ia" in rest):
        return "IN-IaIN"
    if rest.endswith("Ia"):
        return "SN-Ia"
    if rest.endswith("Ib"):
        return "SN-Ib"
    if rest.endswith("II"):
        return "SN-II"
    if rest.endswith("RE"):
        return "RC"
    if rest.endswith("MN ext") or rest.endswith("MN flx"):
        return "MN"
    if rest.endswith("heel contact"):
        return "SN-heel"
    if rest.endswith("toe contact"):
        return "SN-toe"
    raise ValueError(f"untyped aproj neuron label: {name!r}")


def edge_tag(src: str, dst: str, atype: str, stype: dict) -> tuple:
    """(tag, sign) for one aproj synapse, by endpoint semantics + aproj type."""
    exc = (stype["eref"] is not None and stype["eref"] > -0.055)
    sign = "exc" if exc else "inh"
    src_side, src_rest = src.split(" ", 1)
    dst_side, dst_rest = dst.split(" ", 1)

    # ---- contact-neuron outputs (the corrected pathway set)
    if src_rest in ("heel contact", "toe contact"):
        if dst_rest == "RG flx" and dst_side != src_side:
            return "contact_contra_flex", "exc"     # NEW crossed flexor edge
        if dst_rest == "Ank MN ext":
            return "contact_ank", "exc"             # heel+toe -> ankle MN
        if dst_rest in ("RG ext", "Hip PF ext", "Knee PF ext",
                        "Hip MN ext", "Knee MN ext"):
            return "contact_hc", "exc"              # heel cascade (weak class)
        raise ValueError(f"unknown contact edge {src}->{dst}")
    # ---- afferent -> half-centers (dedicated weak types)
    if atype.startswith(("Ia_HipFlx", "Ia_KneeFlx", "II_HipFlx", "II_KneeFlx",
                         "Ib_HipExt", "Ib_KneeExt")) and \
            (dst_rest in ("RG ext", "RG flx") or "_PF" in atype or
             "PF" in dst_rest):
        return "aff_HC_excite", "exc"
    # ---- Ib autogenic -> own MN
    if atype in ("Hip_Ib_ext_Excite", "Hip_Ib_flx_Excite") and " MN " in dst_rest:
        return "ib_auto", "exc"
    # ---- Ia family
    if atype == "Ia to MN Inhibit":
        return "ia_recip", "inh"
    if atype == "Hip_Ia_ext_Excite":            # Ia IN -> Ia (lamination)
        return "ia_recip", "exc"
    if atype == "Ia to Ia Inhibit":
        return "ia_mutual", "inh"
    # ---- PF layer
    if atype in ("PF Excite", "PF Inhibit"):
        return "pf_cross", sign
    if atype == "PF to MN Excite":
        return "pf_to_mn", "exc"
    if atype == "PF to Ia Excite":
        return "other", "exc"
    # ---- Renshaw
    if atype == "MN to RE Excite":
        return "rc_own", "exc"
    if atype == "RE to MN Inhibit":
        return "rc_own", "inh"
    if atype == "RE to RE Inhibit":
        return "rc_mutual", "inh"
    if atype == "RE to Ia Inhibit":
        return "other", "inh"
    # ---- the flagged RE->RE edge carrying a V3-clone type: Ben's 2026-10-03
    # export rules it inh/rc_mutual (README DEVIATIONS family).
    if atype == "V3 Commissural Excite":
        return "rc_mutual", "inh"
    # ---- RG layer
    if atype == "RG to PF Excite":
        return "pf_drive", "exc"
    if atype == "RG to RG Excite":
        return "rg_laminate", "exc"       # builder routes HC->HC to the weak
                                          # direct-excitation gain
    if atype == "RG Excite":
        if src_rest == "V3":
            return "v3_to_contra_InE", "exc"
        if src_rest == "c1":
            return "comm_c1", "exc"
        if dst_rest == "c1":
            return "comm_c1", "exc"
        if dst_rest == "V3":
            return "comm_v3", "exc"
        return "rg_laminate", "exc"       # RG -> own IN lamination
    if atype == "RG Inhibit":
        if src_rest == "c1":
            return "c1_SynAmp2.749", "inh"
        return "rg_laminate", "inh"       # IN -> antagonist RG half-center
    raise ValueError(f"unmapped aproj synapse type {atype!r} ({src}->{dst})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="diff only; do not write the template key")
    args = ap.parse_args()

    mine = mine_aproj(APROJ)
    Path(MINE_OUT).write_text(json.dumps(mine, indent=1), encoding="utf-8")
    print(f"mine: {len(mine['nodes'])} nodes, {len(mine['edges'])} synapses "
          f"-> {MINE_OUT.name}")

    by_guid = {n["guid"]: n for n in mine["nodes"]}
    # relabel the R-knee-ext Renshaw (label fix; GUID-traced, see docstring #6)
    label = {}
    for n in mine["nodes"]:
        if n["class"] != "NonSpiking":
            continue
        nm = n["name"]
        if nm == "L Knee MN ext RE" and n["guid"].startswith(RE_KNEE_EXT_R_GUID_PREFIX):
            nm = "R Knee MN ext RE"
        label[n["guid"]] = nm

    # ---- template nodes
    old = json.loads(TEMPLATES.read_text(encoding="utf-8"))
    old_tpl = old["bilateralrg"]
    old_pos = {n["label"]: (n.get("x", 0), n.get("y", 0)) for n in old_tpl["nodes"]}
    old_id = {n["label"]: n.get("id") for n in old_tpl["nodes"]}

    nodes, seen = [], set()
    for n in mine["nodes"]:
        if n["class"] != "NonSpiking":
            continue
        lab = label[n["guid"]]
        if lab in seen:
            raise ValueError(f"duplicate neuron label after relabel: {lab}")
        seen.add(lab)
        x, y = old_pos.get(lab, (0, 0))
        ent = {"type": node_type(lab), "label": lab, "x": x, "y": y}
        if lab in old_id and old_id[lab]:
            ent["id"] = old_id[lab]
        nodes.append(ent)

    # muscles + stimuli: the physical plant and kickoff are unchanged -
    # carried over verbatim from the bilateralrg template.
    for n in old_tpl["nodes"]:
        if n["type"] in ("MUSCLE", "PORT-load"):
            nodes.append(n)

    # ---- template edges
    stype_by_name = mine["syn_types"]
    edges = []
    for e in mine["edges"]:
        src = label.get(e["src_guid"]) or by_guid[e["src_guid"]]["name"]
        dst = label.get(e["dst_guid"]) or by_guid[e["dst_guid"]]["name"]
        st = stype_by_name[e["type"]]
        tag, sign = edge_tag(src, dst, e["type"], st)
        ent = {"from": src, "to": dst, "sign": sign, "tag": tag,
               "aproj_type": e["type"],
               "aproj_gmax": e["gmax"],
               "synamp": e.get("synamp_instance"),
               "delay_ms": e.get("delay_ms")}
        edges.append(ent)
    for e in old_tpl["edges"]:
        if e["tag"] in ("mn_to_muscle",) or e["from"].startswith("Stimulus"):
            edges.append(e)

    # ---- diff vs bilateralrg
    def keyset(tpl):
        return Counter((e["from"], e["to"], e["tag"], e["sign"])
                       for e in tpl["edges"])

    oldk, newk = keyset(old_tpl), keyset({"edges": edges})
    removed = oldk - newk
    added = newk - oldk
    print(f"\n== delta vs bilateralrg ({len(old_tpl['edges'])} -> {len(edges)} edges) ==")
    print(f"-- REMOVED ({sum(removed.values())}):")
    for (f, t, tag, s), n in sorted(removed.items()):
        print(f"   {n}x {f:<22} -> {t:<22} {tag:<20} {s}")
    print(f"++ ADDED ({sum(added.values())}):")
    for (f, t, tag, s), n in sorted(added.items()):
        print(f"   {n}x {f:<22} -> {t:<22} {tag:<20} {s}")

    tag_c = Counter(e["tag"] for e in edges)
    print("\n== new template tag census ==")
    for k, v in sorted(tag_c.items()):
        print(f"   {v:4d}  {k}")
    n_neurons = sum(1 for n in nodes if n["type"] not in ("MUSCLE", "PORT-load"))
    n_syn = sum(1 for e in edges
                if e["tag"] != "mn_to_muscle" and not e["from"].startswith("Stimulus"))
    print(f"\ncensus: {n_neurons} neurons / {n_syn} synapses / "
          f"12 muscles / 2 stimuli")

    if args.check:
        return 0

    tpl = {
        "_meta": ("bilateralrg2 - Ben's corrected W2L aproj 2026-10-18 "
                  "(physical contact sensing closed, Li-style chain; mined by "
                  "w2l_cpg/make_bilateralrg2.py - READ ITS DOCSTRING FOR THE "
                  "7-POINT DELTA). Census 96/204/12."),
        "_census": {"neurons": n_neurons, "synapses": n_syn, "outputs": 12},
        "nodes": nodes,
        "edges": edges,
    }
    old["bilateralrg2"] = tpl
    TEMPLATES.write_text(json.dumps(old, indent=1), encoding="utf-8")
    print(f"\nwrote 'bilateralrg2' -> {TEMPLATES.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
