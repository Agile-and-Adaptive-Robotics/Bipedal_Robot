r"""Tuned corrected-connectome winner -> AnimatLab realization (2026-10-20, egg).

Ben's ask (2026-10-20): port the TUNED corrected-connectome winner
(reports_20261020\w2l_contact2\round_best_r030.json, study
curr_w2l_contact2_r030, obj20 91.36 / m30 survive) back into his
AnimatLab model and run the drop protocol headless.

NEVER touches Ben's originals. Writes TWO new files into
Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\:
    Walker_2_Layer_CPG_BilateralRG_tuned20261020.aproj   (GUI vehicle)
    Walker_2_Layer_CPG_BilateralRG_tuned20261020.asim    (headless runner)
Both are PATCHED COPIES; Ben commits or deletes them.
Label: HEADLESS-RUN ONLY, GUI OPEN-CHECK PENDING BEN.

THE WINNER (round_best_r030.json):
    te 5.26781, tf 7.5507 nA; ch_hc 0.47826, ch_ank 0.91176,
    ch_contra 1.10224; ctn_amp 4.60959 nA; cref 129.28324 N;
    ibnA 0.16671; comm 0.70048.

THE MAPPING (MuJoCo knob -> AnimatLab edit; ratio-based, the two stacks'
gain scales differ - the campaign net maps FAMILIES to SNS-stack values,
the aproj carries its own class values; the winner applies as a RATIO vs
the knob's MuJoCo default, applied to the aproj type's own base gmax):
  te/tf          4 tonic TonicCurrent stimuli on L/R RG ext / RG flx
                 (the MuJoCo runner injects u[itE]=te / u[itF]=tf directly
                 into those four neurons - same semantics).
  ch_hc  x(0.47826/0.08=5.97825)  the 5 weak Heel_Ipsi_* cascade types
                 (base 0.01 uS -> 0.0597825 uS).
  ch_ank x(0.91176/1.0=0.91176)   Heel_Ipsi_AnkMNext_Excite (base 0.5 uS
                 -> 0.45588 uS) AND Toe_Ipsi_AnkMNext_Excite (base 0.01
                 -> 0.0091176 uS) - both families ride contact_ank in the
                 MuJoCo net (make_bilateralrg2 edge_tag).
  ch_contra x(1.10224/0.08=13.778) Heel_Contra_RGflx_Excite (0.01 ->
                 0.13778 uS).
  comm  x0.70048 the two commissural OUTPUT families. In the corrected
                 aproj these edges share the big "RG Inhibit" / "RG Excite"
                 types (the stock "c1 Commissural Inhibit" type is UNUSED;
                 c1->contra RG flx carries RG Inhibit, V3->contra RG ext IN
                 carries RG Excite). DEDICATED clone types
                 "RG Inhibit (c1 commissural, tuned 20261020)" and
                 "RG Excite (V3 contra, tuned 20261020)" at 2.749*comm =
                 1.92562 uS carry the 2+2 edges; the shared types stay at
                 2.749 uS (RG flx->c1 / RG ext->V3 drive edges are NOT
                 comm-scaled - the campaign's gains dict leaves comm_c1/
                 comm_v3 at 1.0).
  ctn_amp        4 ContactCount adapter Gain C -> 4.60959 nA (asim plain
                 SI 4.60959e-009; aproj Value/Scale/Actual). NOTE: the
                 corrected aproj ships C=1 (SI amp - saturated binary
                 switch, the same unit family as the old C=20 bug); the
                 tuned value makes the contact-neuron input a sane nA-scale
                 current equal to the winner's ENCODER AMPLITUDE.
  cref           NOT TRANSFERRED (documented N/A, same ruling as the
                 2026-10-16 report): the aproj adapters are ContactCount-
                 based (no force input; no force-based SourceDataType
                 exists anywhere in the repo's models). The AnimatLab port
                 therefore runs the contact encoder in its SATURATED
                 regime (any contact = full ctn_amp); the winner's
                 sub-saturation linear force scaling (cref 129 N ~ 31%
                 body weight) does not transfer. FLAGGED for Ben.
  ibnA   x0.16671 all 12 per-muscle Tension (Ib) adapter Gain C (8x
                 5.062 pA/N class, 3x 13.842, 1x 2.963 - Ben's asymmetric
                 hip values preserved; uniform ratio, the MuJoCo ibnA
                 multiplies the side Ib group port).
  L RG ext intrinsic 2 nA tonic: the corrected APROJ carries it; the stale
                 asim export has 0 -> the asim copy is SYNCED to 2e-009
                 (aproj literal wins; the MuJoCo net has no such field -
                 known aproj<->MuJoCo delta, documented).

THE ASIM GAP (why the asim needs the corrected-connectome delta): Ben's
two standalone asims predate the whole 2026-10-18 corrected connectome
(0 of the 20 dedicated types). The asim copy is built from the
battle-tested Ground_Standalone.asim + the EXACT connexion delta vs the
corrected aproj (computed at build time from w2l_aproj_mine_20261018.json
+ a live re-mine of the aproj): 39 retypes, 6 adds, 10 deletes - the
documented 7-point correction - with the TUNED conductances on top.
Neuron sets are GUID-identical (96/96, verified).

Sim window: SimEndTime 5.1 -> 30.1 s in both files; asim charts extended
to [0,30] s, CollectInterval 0.0002 -> 0.001 s (30 s of 0.2 ms rows would
be ~15 MB per chart; 1 ms resolves ~0.4 s RG periods with 400 samples).
A new "Root Track" asim chart logs Root WorldPositionX/Y/Z + hip_L/hip_R
JointRotation (drop/topple telemetry; WorldPositionX is a proven GUI
SelectedDataType - if the sim-side enum rejects it the chart writes 0
bytes but the run is unaffected; the other 16 charts carry the analysis).

Usage (egg):
    D:\Anaconda\envs\myo\python.exe build_tuned_variant_20261020.py
"""
from __future__ import annotations

import json
import os
import re
import sys
import xml.etree.ElementTree as ET

REPO = r"D:\GitHub\Bipedal_Robot"
MODEL_DIR = os.path.join(REPO, "Neuromechanical_Models",
                         "Walker_2_Layer_CPG_BilateralRG")
SRC_APROJ = os.path.join(MODEL_DIR, "Walker_2_Layer_CPG_BilateralRG.aproj")
SRC_ASIM = os.path.join(MODEL_DIR,
                        "Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim")
OUT_APROJ = os.path.join(MODEL_DIR,
                         "Walker_2_Layer_CPG_BilateralRG_tuned20261020.aproj")
OUT_ASIM = os.path.join(MODEL_DIR,
                        "Walker_2_Layer_CPG_BilateralRG_tuned20261020.asim")
MINE = os.path.join(REPO, "Code", "MuJoCo_SNS", "spinal",
                    "w2l_aproj_mine_20261018.json")
HERE = os.path.dirname(os.path.abspath(__file__))
EDITS_LOG = os.path.join(HERE, "edits_log_20261020.json")

# ---- the winner (round_best_r030.json, verbatim) -------------------------
TE, TF = 5.26781, 7.5507
CH_HC, CH_ANK, CH_CONTRA = 0.47826, 0.91176, 1.10224
CTN_AMP, CREF = 4.60959, 129.28324
IBNA, COMM = 0.16671, 0.70048

# MuJoCo knob defaults (build_w2l_net W2L_GAINS) the ratios are taken vs
DEF_HC, DEF_ANK, DEF_CONTRA = 0.08, 1.0, 0.08
R_HC, R_ANK, R_CONTRA = CH_HC / DEF_HC, CH_ANK / DEF_ANK, CH_CONTRA / DEF_CONTRA

SIM_END = 30.1
CHART_END = 30
NEWPFX = "d71a"

# aproj class bases (S)
HEEL_WEAK = 1e-8          # Heel_Ipsi cascade + Heel_Contra + Toe_Ipsi_Ank
HEEL_ANK = 5e-7           # Heel_Ipsi_AnkMNext (generic class)
RG_COMM = 2.749e-6        # c1 / V3-contra edges (shared type base)


def g(n: int) -> str:
    return f"{NEWPFX}{n:04d}-0000-4000-8000-{900000000000 + n}"


def fmt(v: float) -> str:
    return f"{v:.6g}"


LOG: list = []


def log(msg: str):
    print(msg)
    LOG.append(msg)


def rd(p: str) -> str:
    with open(p, "r", encoding="utf-8", errors="surrogateescape", newline="") as f:
        return f.read()


def wr(p: str, s: str):
    with open(p, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.write(s)


# --------------------------------------------------------------------------
# mine the corrected aproj: types + edges (GUID-level)
def mine_types(x: str) -> dict:
    """name -> {id, eref_mV, gmax_S, thresh_mV, sat_mV} from SynapseTypes."""
    s = x.find("<SynapseTypes>")
    e = x.find("</SynapseTypes>")
    out = {}
    for m in re.finditer(r"<Link>(.*?)</Link>", x[s:e], re.S):
        b = m.group(1)
        nm = re.search(r"<Name>([^<]*)</Name>", b)
        if not nm or "MaxSynapticConductance" not in b:
            continue
        gid = re.search(r"<ID>([^<]+)</ID>", b).group(1)

        def val(tag):
            mm = re.search(rf'<{tag} Value="([^"]*)" Scale="([^"]*)" Actual="([^"]*)"/>', b)
            return mm.groups() if mm else None
        out[nm.group(1)] = dict(
            id=gid, gmax=val("MaxSynapticConductance"),
            eref=val("EquilibriumPotential"), thresh=val("PreSynapticThreshold"),
            sat=val("PreSynapticSaturationLevel"))
    return out


def main() -> int:
    control = "--control" in sys.argv   # Ben's literal corrected model, no drive
    for p in (SRC_APROJ, SRC_ASIM, MINE):
        if not os.path.exists(p):
            raise RuntimeError(f"missing input: {p}")
    suffix = "_control" if control else ""
    global OUT_APROJ, OUT_ASIM
    OUT_APROJ = OUT_APROJ.replace("_tuned20261020.aproj", f"_control20261020.aproj") if control else OUT_APROJ
    OUT_ASIM = OUT_ASIM.replace("_tuned20261020.asim", f"_control20261020.asim") if control else OUT_ASIM
    xa = rd(SRC_APROJ)
    xs = rd(SRC_ASIM)
    if NEWPFX in xa or NEWPFX in xs:
        raise RuntimeError(f"GUID prefix {NEWPFX} already present - refusing")

    mine = json.load(open(MINE, encoding="utf-8"))
    nodes = {n["guid"]: (n["name"] or "?") for n in mine["nodes"]}
    ap_edges = {(e["src_guid"], e["dst_guid"]): e["type"] for e in mine["edges"]}
    atypes = mine_types(xa)
    log(f"mined aproj: {len(ap_edges)} edges, {len(atypes)} synapse types")

    # tuned type conductances (S) - in --control mode everything stays at the
    # corrected aproj's literal base values
    heel_ank_names = ["Heel_Ipsi_RGext_Excite", "Heel_Ipsi_HipMNext_Excite",
                      "Heel_Ipsi_HipPFext_Excite", "Heel_Ipsi_KneeMNext_Excite",
                      "Heel_Ipsi_KneePFext_Excite", "Heel_Ipsi_AnkMNext_Excite",
                      "Toe_Ipsi_AnkMNext_Excite", "Heel_Contra_RGflx_Excite"]
    heel_ank_base = [HEEL_WEAK, HEEL_WEAK, HEEL_WEAK, HEEL_WEAK, HEEL_WEAK,
                     HEEL_ANK, HEEL_WEAK, HEEL_WEAK]
    if control:
        tuned = dict(zip(heel_ank_names, heel_ank_base))
    else:
        tuned = {
            "Heel_Ipsi_RGext_Excite":   HEEL_WEAK * R_HC,
            "Heel_Ipsi_HipMNext_Excite": HEEL_WEAK * R_HC,
            "Heel_Ipsi_HipPFext_Excite": HEEL_WEAK * R_HC,
            "Heel_Ipsi_KneeMNext_Excite": HEEL_WEAK * R_HC,
            "Heel_Ipsi_KneePFext_Excite": HEEL_WEAK * R_HC,
            "Heel_Ipsi_AnkMNext_Excite": HEEL_ANK * R_ANK,
            "Toe_Ipsi_AnkMNext_Excite":  HEEL_WEAK * R_ANK,
            "Heel_Contra_RGflx_Excite":  HEEL_WEAK * R_CONTRA,
        }
    li_types = [k for k in atypes if k.startswith(("Ia_", "II_", "Ib_"))
                and k.endswith("_Excite")]
    for k in li_types:
        tuned[k] = 1e-8                    # weak class, unchanged
    C1_NEW_NAME = "RG Inhibit (c1 commissural, tuned 20261020)"
    V3_NEW_NAME = "RG Excite (V3 contra, tuned 20261020)"
    COMM_GMAX = RG_COMM * COMM
    C1_ID, V3_ID = g(1), g(2)

    rg_of = {name: gid for gid, name in nodes.items()}

    # ======================= APROJ =========================================
    p = xa
    # --- 1) retune the existing corrected-connectome types ----------------
    for name, newS in tuned.items():
        blk = _find_type_block(p, name)
        old = atypes[name]["gmax"]
        new = (fmt(newS * 1e6), newS)
        p = _patch_type_gmax(p, blk, name, old, new)
        log(f"aproj: type '{name}' gmax {old[2]} -> {new[1]:.6g} S")

    # --- 2) dedicated comm types (clone RG Inhibit / RG Excite) ------------
    if control:
        log("CONTROL: comm types stay shared (Ben's literal wiring)")
    else:
        p = _clone_type(p, "RG Inhibit", C1_ID, C1_NEW_NAME, COMM_GMAX, atypes)
        p = _clone_type(p, "RG Excite", V3_ID, V3_NEW_NAME, COMM_GMAX, atypes)
        log(f"aproj: +types '{C1_NEW_NAME}' / '{V3_NEW_NAME}' @ {fmt(COMM_GMAX*1e6)} uS")

        # repoint the 2+2 commissural output links
        c1_pairs = [(rg_of[s], rg_of[d]) for s, d in
                    (("L c1", "R RG flx"), ("R c1", "L RG flx"))]
        v3_pairs = [(rg_of[s], rg_of[d]) for s, d in
                    (("L V3", "R RG ext IN"), ("R V3", "L RG ext IN"))]
        p = _repoint_aproj_links(p, c1_pairs, C1_ID, COMM_GMAX, "c1->contra RG flx")
        p = _repoint_aproj_links(p, v3_pairs, V3_ID, COMM_GMAX, "V3->contra RG ext IN")

    # --- 3) adapters --------------------------------------------------------
    if not control:
        p, n_cc = _patch_aproj_contact(p, CTN_AMP)
        log(f"aproj: {n_cc} ContactCount adapters C -> {CTN_AMP} nA")
        p, n_ib = _patch_aproj_ib(p, IBNA)
        log(f"aproj: {n_ib} Tension(Ib) adapter C scaled x{IBNA}")

    # --- 4) stimuli (clone Stimulus_1) -------------------------------------
    if control:
        log("CONTROL: no gamma drive (kickoffs only, Ben's literal model)")
    else:
        p = _add_aproj_stimuli(p, rg_of, TE, TF)
        log("aproj: +4 gamma stimuli (te/tf on the four RG halves)")

    # --- 5) sim end ---------------------------------------------------------
    p, n_se = re.subn(r'<SimEndTime Value="5.1" Scale="None" Actual="5.1"/>',
                      f'<SimEndTime Value="{SIM_END}" Scale="None" Actual="{SIM_END}"/>', p)
    if n_se != 1:
        raise RuntimeError(f"aproj SimEndTime patches: {n_se}")
    log(f"aproj: SimEndTime 5.1 -> {SIM_END}")
    wr(OUT_APROJ, p)

    # ======================= ASIM ==========================================
    s = xs
    stypes = _asim_types(s)                       # name -> (id, synamp, blk span)
    # --- 1) add the corrected-connectome types (tuned, or literal base) ----
    add_blocks = []
    for name, newS in tuned.items():
        t = atypes[name]
        add_blocks.append(_asim_type_block(
            t["id"], name, newS * 1e6,
            float(t["eref"][0]), float(t["thresh"][0]), float(t["sat"][0])))
    for name in li_types:                          # li types ride `tuned`
        pass
    if not control:
        t_inh, t_exc = atypes["RG Inhibit"], atypes["RG Excite"]
        add_blocks.append(_asim_type_block(
            C1_ID, C1_NEW_NAME, COMM_GMAX * 1e6, float(t_inh["eref"][0]),
            float(t_inh["thresh"][0]), float(t_inh["sat"][0])))
        add_blocks.append(_asim_type_block(
            V3_ID, V3_NEW_NAME, COMM_GMAX * 1e6, float(t_exc["eref"][0]),
            float(t_exc["thresh"][0]), float(t_exc["sat"][0])))
    # insert after the last CHEMICAL type close (the file's overall LAST
    # </SynapseType> belongs to the ElectricalSynapses/LowCoup container -
    # appending there breaks the load: "element ... not found ... LowCoup")
    elec_i = s.find("<ElectricalSynapses>")
    anchor = s.rfind("</SynapseType>", 0, elec_i) + len("</SynapseType>")
    s = s[:anchor] + "\n" + "\n".join(add_blocks) + s[anchor:]
    log(f"asim: +{len(add_blocks)} SynapseTypes "
        f"({'corrected literal base' if control else 'incl. 2 comm dedicated'})")

    # --- 2) apply the connexion delta vs the corrected aproj ---------------
    name_to_id = {k: v[0] for k, v in stypes.items()}
    name_to_id.update({k: atypes[k]["id"] for k in list(tuned) + li_types})
    name_to_id[C1_NEW_NAME], name_to_id[V3_NEW_NAME] = C1_ID, V3_ID
    retarget = {}
    asim_edges = _asim_connexions(s)              # (src,dst) -> typeid
    for key, atype in ap_edges.items():
        if key not in asim_edges:
            continue
        cur = _typeid_name(asim_edges[key], stypes)
        want = atype
        # c1 / V3 commissural outputs -> dedicated tuned types (shared types
        # in --control mode: Ben's literal corrected wiring keeps them there)
        if (nodes[key[0]], nodes[key[1]]) in (("L c1", "R RG flx"), ("R c1", "L RG flx")):
            want = C1_NEW_NAME if not control else "RG Inhibit"
        elif (nodes[key[0]], nodes[key[1]]) in (("L V3", "R RG ext IN"),
                                                ("R V3", "L RG ext IN")):
            want = V3_NEW_NAME if not control else "RG Excite"
        if cur != want:
            retarget[key] = name_to_id[want]
    n_ret = 0
    for key, tid in retarget.items():
        s, ok = _asim_retype(s, key, tid)
        n_ret += ok
    log(f"asim: retyped {n_ret}/{len(retarget)} connexions to corrected types")

    # deletes (toe cascade, 10) + adds (ankle + crossed, 6)
    del_keys = [k for k in asim_edges if k not in ap_edges]
    for key in del_keys:
        s, ok = _asim_del_connexion(s, key)
        if not ok:
            raise RuntimeError(f"asim delete failed {key}")
    log(f"asim: deleted {len(del_keys)} connexions (toe cascade)")

    add_edges = [k for k in ap_edges if k not in asim_edges]
    tmpl = _asim_connexion_template(s)
    for n, key in enumerate(sorted(add_edges)):
        tid = name_to_id[ap_edges[key]]
        blk = tmpl.replace(f"<ID>{re.search(r'<ID>([^<]+)</ID>', tmpl).group(1)}</ID>",
                           f"<ID>{g(30 + n)}</ID>", 1)
        blk = re.sub(r"<SourceID>[^<]*</SourceID>",
                     f"<SourceID>{key[0]}</SourceID>", blk, count=1)
        blk = re.sub(r"<TargetID>[^<]*</TargetID>",
                     f"<TargetID>{key[1]}</TargetID>", blk, count=1)
        blk = re.sub(r"<SynapseTypeID>[^<]*</SynapseTypeID>",
                     f"<SynapseTypeID>{tid}</SynapseTypeID>", blk, count=1)
        anchor = s.rfind("</Connexion>") + len("</Connexion>")
        s = s[:anchor] + "\n" + blk + s[anchor:]
    log(f"asim: added {len(add_edges)} connexions (ankle + crossed flexor)")

    # --- 3) adapters --------------------------------------------------------
    if control:
        # Ben's literal C=1 (SI amps, the saturated binary contact switch)
        s, n_cc2 = _patch_asim_contact(s, None, raw_c="1")
        log(f"asim: CONTROL {n_cc2} ContactCount adapters C 20 -> 1 (Ben's literal)")
    else:
        s, n_cc2 = _patch_asim_contact(s, CTN_AMP)
        log(f"asim: {n_cc2} ContactCount adapters C -> {CTN_AMP}e-009 A")
        s, n_ib2 = _patch_asim_ib(s, IBNA, [gid for gid, nm in nodes.items()
                                            if nm.endswith("Ib")])
        log(f"asim: {n_ib2} Tension(Ib) adapter C scaled x{IBNA}")

    # --- 4) L RG ext intrinsic tonic sync (aproj literal 2 nA) -------------
    i = s.find("<ID>17fdb1f6-")
    blk_s = s.rfind("<Neuron>", 0, i)
    blk_e = s.find("</Neuron>", i)
    blk = s[blk_s:blk_e]
    nb = blk.replace("<TonicStimulus>0</TonicStimulus>",
                     "<TonicStimulus>2e-009</TonicStimulus>", 1)
    if nb == blk:
        raise RuntimeError("asim L RG ext tonic patch failed")
    s = s[:blk_s] + nb + s[blk_e:]
    log("asim: L RG ext intrinsic TonicStimulus 0 -> 2e-009 (aproj sync)")

    # --- 5) stimuli ---------------------------------------------------------
    if not control:
        s = _add_asim_stimuli(s, rg_of, TE, TF)
        log("asim: +4 gamma stimuli (te/tf on the four RG halves)")

    # --- 6) sim window + charts --------------------------------------------
    s, n_se2 = re.subn(r"<SimEndTime>5.1</SimEndTime>",
                       f"<SimEndTime>{SIM_END}</SimEndTime>", s)
    if n_se2 != 1:
        raise RuntimeError(f"asim SimEndTime patches: {n_se2}")
    s, n_ce = re.subn(r"<EndTime>5</EndTime>", f"<EndTime>{CHART_END}</EndTime>", s)
    s, n_cs = re.subn(r"<StartTime>0.5</StartTime>", "<StartTime>0</StartTime>", s)
    s, n_ci = re.subn(r"<CollectInterval>0.0002</CollectInterval>",
                      "<CollectInterval>0.001</CollectInterval>", s)
    log(f"asim: charts EndTime 5->{CHART_END} x{n_ce}, StartTime 0.5->0 x{n_cs}, "
        f"CollectInterval 0.0002->0.001 x{n_ci}; SimEndTime -> {SIM_END}")

    # Root Track chart (drop telemetry)
    joints = _asim_joint_guids(s)
    i_root = s.find("<Name>Root</Name>")
    root_id = re.search(r"<ID>([^<]+)</ID>", s[i_root:i_root + 200]).group(1)
    cols = [("RootX", "WorldPositionX", root_id),
            ("RootY", "WorldPositionY", root_id),
            ("RootZ", "WorldPositionZ", root_id),
            ("hip_L", "JointRotation", joints["hip_L"]),
            ("hip_R", "JointRotation", joints["hip_R"])]
    chart = [f"<DataChart>\n<ID>{g(41)}</ID>\n<Type>FileChart</Type>\n"
             f"<OutputFilename>Root Track.txt</OutputFilename>\n"
             f"<Name>Root Track</Name>\n<SetStartEndTime>True</SetStartEndTime>\n"
             f"<CollectTimeWindow>-1</CollectTimeWindow>\n"
             f"<AlwaysActive>False</AlwaysActive>\n<StartTime>0</StartTime>\n"
             f"<EndTime>{CHART_END}</EndTime>\n<CollectInterval>0.001</CollectInterval>\n"
             f"<DataColumns>"]
    for n, (cn, dt, tid) in enumerate(cols):
        chart.append(f"<DataColumn>\n<ModuleName/>\n<Type>DataColumn</Type>\n"
                     f"<ID>{g(42 + n)}</ID>\n<ColumnName>{cn}</ColumnName>\n"
                     f"<DataType>{dt}</DataType>\n<TargetID>{tid}</TargetID>\n"
                     f"</DataColumn>")
    chart.append("</DataColumns>\n</DataChart>")
    a2 = s.rfind("</DataChart>") + len("</DataChart>")
    s = s[:a2] + "\n" + "\n".join(chart) + s[a2:]
    log("asim: +Root Track chart (WorldPosition XYZ + hip rotations)")

    wr(OUT_ASIM, s)

    # ======================= VERIFY ========================================
    verify(OUT_APROJ, aproj=True)
    verify(OUT_ASIM, aproj=False)
    n_syn = len(_asim_connexions(rd(OUT_ASIM)))
    log(f"VERIFY: asim connexions = {n_syn} (aproj has {len(ap_edges)})")
    if n_syn != len(ap_edges):
        raise RuntimeError(f"connexion count mismatch {n_syn} vs {len(ap_edges)}")
    json.dump({"log": LOG, "winner": dict(te=TE, tf=TF, ch_hc=CH_HC, ch_ank=CH_ANK,
               ch_contra=CH_CONTRA, ctn_amp=CTN_AMP, cref=CREF, ibnA=IBNA,
               comm=COMM),
               "tuned_types_S": {k: v for k, v in tuned.items()},
               "comm_dedicated_S": COMM_GMAX},
              open(EDITS_LOG, "w", encoding="utf-8"), indent=2)
    print("BUILT + VERIFIED:", OUT_APROJ, "and", OUT_ASIM, sep="\n  ")
    return 0


# ---------------------------------------------------------------- helpers
def _find_type_block(x: str, name: str):
    i = x.find(f"<Name>{name}</Name>", x.find("<SynapseTypes>"))
    if i < 0:
        raise RuntimeError(f"type {name} not found in SynapseTypes")
    s = x.rfind("<Link>", 0, i)
    e = x.find("</Link>", i)
    return s, e, x[s:e]


def _patch_type_gmax(x, blk, name, old, new):
    s, e, b = blk
    oldv = (f'<MaxSynapticConductance Value="{old[0]}" Scale="{old[1]}" '
            f'Actual="{old[2]}"/>')
    newv = (f'<MaxSynapticConductance Value="{new[0]}" Scale="micro" '
            f'Actual="{fmt(float(new[1]))}"/>')
    if oldv not in b:
        raise RuntimeError(f"{name}: gmax triplet not found: {oldv}")
    return x[:s] + b.replace(oldv, newv, 1) + x[e:]


def _clone_type(x, base_name, new_id, new_name, gmax_S, atypes):
    s, e, b = _find_type_block(x, base_name)
    full = x[s:e + len("</Link>")]          # block WITH its closer
    old_id = re.search(r"<ID>([^<]+)</ID>", full).group(1)
    nb = full.replace(f"<ID>{old_id}</ID>", f"<ID>{new_id}</ID>", 1)
    nb = re.sub(r"<Name>[^<]*</Name>", f"<Name>{new_name}</Name>", nb, count=1)
    old = atypes[base_name]["gmax"]
    newv = (f'<MaxSynapticConductance Value="{fmt(gmax_S*1e6)}" Scale="micro" '
            f'Actual="{fmt(gmax_S)}"/>')
    oldv = (f'<MaxSynapticConductance Value="{old[0]}" Scale="{old[1]}" '
            f'Actual="{old[2]}"/>')
    nb = nb.replace(oldv, newv, 1)
    anchor = e + len("</Link>")             # sibling AFTER the base block
    return x[:anchor] + "\n" + nb + x[anchor:]


def _repoint_aproj_links(x, pairs, new_tid, gmax_S, what):
    want = set(pairs)
    n = 0
    out, pos = [], 0
    for m in re.finditer(r"<Link>(.*?)</Link>", x, re.S):
        b = m.group(1)
        if "Behavior.Synapse" not in b:
            continue
        o = re.search(r"<OriginID>([^<]*)</OriginID>", b)
        d = re.search(r"<DestinationID>([^<]*)</DestinationID>", b)
        if not (o and d and (o.group(1), d.group(1)) in want):
            continue
        nb = re.sub(r"<SynapticTypeID>[^<]*</SynapticTypeID>",
                    f"<SynapticTypeID>{new_tid}</SynapticTypeID>", b, count=1)
        nb = re.sub(r'<SynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>',
                    f'<SynapticConductance Value="{fmt(gmax_S*1e6)}" Scale="micro" '
                    f'Actual="{fmt(gmax_S)}"/>', nb, count=1)
        out.append(x[pos:m.start()])
        out.append(f"<Link>{nb}</Link>")
        pos = m.end()
        n += 1
    out.append(x[pos:])
    log(f"aproj: repointed {n}/{len(pairs)} links ({what})")
    if n != len(pairs):
        raise RuntimeError(f"aproj repoint {what}: {n} != {len(pairs)}")
    return "".join(out)


def _patch_aproj_contact(x, ctn_na):
    out, pos, n = [], 0, 0
    for m in re.finditer(r"<Node>(.*?)</Node>", x, re.S):
        b = m.group(1)
        if "PhysicalToNodeAdapter" not in b or "<DataTypeID>ContactCount" not in b:
            continue
        C = re.search(r'<C Value="([^"]*)" Scale="([^"]*)" Actual="([^"]*)"/>', b)
        if not C or C.group(2) != "None":
            continue
        newc = f'<C Value="{ctn_na:.6g}" Scale="nano" Actual="{ctn_na:.6g}e-009"/>'
        nb = b.replace(C.group(0), newc, 1)
        out.append(x[pos:m.start()])
        out.append(f"<Node>{nb}</Node>")
        pos = m.end()
        n += 1
    out.append(x[pos:])
    return "".join(out), n


def _patch_aproj_ib(x, ratio):
    out, pos, n = [], 0, 0
    for m in re.finditer(r"<Node>(.*?)</Node>", x, re.S):
        b = m.group(1)
        if "PhysicalToNodeAdapter" not in b or "<DataTypeID>Tension" not in b:
            continue
        C = re.search(r'<C Value="([^"]*)" Scale="([^"]*)" Actual="([^"]*)"/>', b)
        if not C:
            continue
        val = float(C.group(1)) * ratio
        newa = val * {"pico": 1e-12, "nano": 1e-9, "micro": 1e-6,
                      "None": 1.0, "milli": 1e-3}[C.group(2)]
        newc = (f'<C Value="{val:.6g}" Scale="{C.group(2)}" '
                f'Actual="{fmt(newa)}"/>')
        nb = b.replace(C.group(0), newc, 1)
        out.append(x[pos:m.start()])
        out.append(f"<Node>{nb}</Node>")
        pos = m.end()
        n += 1
    out.append(x[pos:])
    return "".join(out), n


def _add_aproj_stimuli(x, rg_of, te, tf):
    i = x.find("<Name>Stimulus_1</Name>", x.find("<Stimuli>"))
    s = x.rfind("<Stimulus>", 0, i)
    e = x.find("</Stimulus>", i) + len("</Stimulus>")
    tmpl = x[s:e]
    specs = [("Gamma drive L RG ext (te=5.26781nA, tuned 20261020)", "L RG ext", te),
             ("Gamma drive L RG flx (tf=7.5507nA, tuned 20261020)", "L RG flx", tf),
             ("Gamma drive R RG ext (te=5.26781nA, tuned 20261020)", "R RG ext", te),
             ("Gamma drive R RG flx (tf=7.5507nA, tuned 20261020)", "R RG flx", tf)]
    ins = ""
    for k, (nm, key, na) in enumerate(specs):
        blk = tmpl
        old_id = re.search(r"<ID>([^<]+)</ID>", blk).group(1)
        blk = blk.replace(f"<ID>{old_id}</ID>", f"<ID>{g(20 + k)}</ID>", 1)
        blk = re.sub(r"<Name>[^<]*</Name>", f"<Name>{nm}</Name>", blk, count=1)
        blk = re.sub(r"<NodeID>[^<]*</NodeID>", f"<NodeID>{rg_of[key]}</NodeID>",
                     blk, count=1)
        blk = re.sub(r'<CurrentOn Value="[^"]*" Scale="nano" Actual="[^"]*"/>',
                     f'<CurrentOn Value="{na:.6g}" Scale="nano" '
                     f'Actual="{na:.6g}e-009"/>', blk, count=1)
        blk = re.sub(r'<EndTime Value="[^"]*" Scale="None" Actual="[^"]*"/>',
                     f'<EndTime Value="{SIM_END}" Scale="None" Actual="{SIM_END}"/>',
                     blk, count=1)
        for tag in ("CycleOnDuration", "BurstOnDuration"):
            blk = re.sub(rf'<{tag} Value="[^"]*" Scale="milli" Actual="[^"]*"/>',
                         f'<{tag} Value="{SIM_END * 1000:.6g}" Scale="milli" '
                         f'Actual="{SIM_END}"/>', blk, count=1)
        ins += "\n" + blk
    return x[:e] + ins + x[e:]


def _asim_types(x):
    out = {}
    for m in re.finditer(r"<SynapseType>(.*?)</SynapseType>", x, re.S):
        b = m.group(1)
        nm = re.search(r"<Name>([^<]*)</Name>", b)
        gid = re.search(r"<ID>([^<]+)</ID>", b)
        sa = re.search(r"<SynAmp>([^<]*)</SynAmp>", b)
        out[nm.group(1)] = (gid.group(1), sa.group(1) if sa else None)
    return out


def _typeid_name(tid, stypes):
    for nm, (gid, _) in stypes.items():
        if gid == tid:
            return nm
    return "?"


def _asim_connexions(x):
    out = {}
    for m in re.finditer(r"<Connexion>(.*?)</Connexion>", x, re.S):
        b = m.group(1)
        s = re.search(r"<SourceID>([^<]*)</SourceID>", b)
        t = re.search(r"<TargetID>([^<]*)</TargetID>", b)
        ty = re.search(r"<SynapseTypeID>([^<]*)</SynapseTypeID>", b)
        out[(s.group(1), t.group(1))] = ty.group(1)
    return out


def _asim_connexion_template(x):
    j = x.find("<Connexion>")
    return x[j:x.find("</Connexion>", j) + len("</Connexion>")]


def _asim_retype(x, key, tid):
    n = 0
    out, pos = [], 0
    for m in re.finditer(r"<Connexion>(.*?)</Connexion>", x, re.S):
        b = m.group(1)
        s = re.search(r"<SourceID>([^<]*)</SourceID>", b)
        t = re.search(r"<TargetID>([^<]*)</TargetID>", b)
        if (s.group(1), t.group(1)) == key:
            nb = re.sub(r"<SynapseTypeID>[^<]*</SynapseTypeID>",
                        f"<SynapseTypeID>{tid}</SynapseTypeID>", b, count=1)
            out.append(x[pos:m.start()])
            out.append(f"<Connexion>{nb}</Connexion>")
            pos = m.end()
            n += 1
    out.append(x[pos:])
    if n != 1:
        raise RuntimeError(f"retype {key}: matched {n}")
    return "".join(out), n


def _asim_del_connexion(x, key):
    n = 0
    out, pos = [], 0
    for m in re.finditer(r"<Connexion>(.*?)</Connexion>", x, re.S):
        b = m.group(1)
        s = re.search(r"<SourceID>([^<]*)</SourceID>", b)
        t = re.search(r"<TargetID>([^<]*)</TargetID>", b)
        if (s.group(1), t.group(1)) == key:
            out.append(x[pos:m.start()])
            pos = m.end()
            n += 1
    out.append(x[pos:])
    if n != 1:
        raise RuntimeError(f"delete {key}: matched {n}")
    return "".join(out), n


def _asim_type_block(gid, name, synamp_uS, equil_mV, thresh_mV, sat_mV):
    return (f"<SynapseType>\n<Name>{name}</Name>\n<ID>{gid}</ID>\n"
            f"<Type>NonSpikingChemical</Type>\n<Equil>{equil_mV:.6g}</Equil>\n"
            f"<SynAmp>{synamp_uS:.6g}</SynAmp>\n"
            f"<ThreshV>{thresh_mV:.6g}</ThreshV>\n"
            f"<SaturateV>{sat_mV:.6g}</SaturateV>\n</SynapseType>")


def _patch_asim_contact(x, ctn_na, raw_c=None):
    new_c = raw_c if raw_c is not None else f"{ctn_na:.6g}e-009"
    out, pos, n = [], 0, 0
    for m in re.finditer(r"<Adapter>(.*?)</Adapter>", x, re.S):
        b = m.group(1)
        if "<SourceDataType>ContactCount" not in b:
            continue
        C = re.search(r"<C>([^<]*)</C>", b)
        D = re.search(r"<D>([^<]*)</D>", b)
        nb = b.replace(f"<C>{C.group(1)}</C>", f"<C>{new_c}</C>", 1)
        if D:
            nb = nb.replace(f"<D>{D.group(1)}</D>", "<D>0</D>", 1)
        out.append(x[pos:m.start()])
        out.append(f"<Adapter>{nb}</Adapter>")
        pos = m.end()
        n += 1
    out.append(x[pos:])
    return "".join(out), n


def _patch_asim_ib(x, ratio, ib_guids):
    want = set(ib_guids)
    out, pos, n = [], 0, 0
    for m in re.finditer(r"<Adapter>(.*?)</Adapter>", x, re.S):
        b = m.group(1)
        if "<SourceDataType>Tension" not in b:
            continue
        t = re.search(r"<TargetID>([^<]*)</TargetID>", b)
        if not t or t.group(1) not in want:
            continue
        C = re.search(r"<C>([^<]*)</C>", b)
        val = float(C.group(1)) * ratio
        nb = b.replace(f"<C>{C.group(1)}</C>", f"<C>{fmt(val)}</C>", 1)
        out.append(x[pos:m.start()])
        out.append(f"<Adapter>{nb}</Adapter>")
        pos = m.end()
        n += 1
    out.append(x[pos:])
    return "".join(out), n


def _add_asim_stimuli(x, rg_of, te, tf):
    i = x.find("<Name>Stimulus_1</Name>")
    s = x.rfind("<Stimulus>", 0, i)
    e = x.find("</Stimulus>", i) + len("</Stimulus>")
    tmpl = x[s:e]
    specs = [("Gamma drive L RG ext (te=5.26781nA, tuned 20261020)", "L RG ext", te),
             ("Gamma drive L RG flx (tf=7.5507nA, tuned 20261020)", "L RG flx", tf),
             ("Gamma drive R RG ext (te=5.26781nA, tuned 20261020)", "R RG ext", te),
             ("Gamma drive R RG flx (tf=7.5507nA, tuned 20261020)", "R RG flx", tf)]
    ins = ""
    for k, (nm, key, na) in enumerate(specs):
        blk = tmpl
        old_id = re.search(r"<ID>([^<]+)</ID>", blk).group(1)
        blk = blk.replace(f"<ID>{old_id}</ID>", f"<ID>{g(10 + k)}</ID>", 1)
        blk = re.sub(r"<Name>[^<]*</Name>", f"<Name>{nm}</Name>", blk, count=1)
        blk = re.sub(r"<TargetNodeID>[^<]*</TargetNodeID>",
                     f"<TargetNodeID>{rg_of[key]}</TargetNodeID>", blk, count=1)
        blk = re.sub(r"<CurrentOn>[^<]*</CurrentOn>",
                     f"<CurrentOn>{na:.6g}e-009</CurrentOn>", blk, count=1)
        for tag in ("EndTime", "CycleOnDuration", "BurstOnDuration"):
            blk = re.sub(rf"<{tag}>[^<]*</{tag}>", f"<{tag}>{SIM_END}</{tag}>",
                         blk, count=1)
        ins += "\n" + blk
    close = x.find("</ExternalStimuli>", e)
    return x[:close] + ins + x[close:]


def _asim_joint_guids(x):
    out = {}
    for m in re.finditer(r"<Joint>(.*?)</Joint>", x, re.S):
        b = m.group(1)
        nm = re.search(r"<Name>([^<]*)</Name>", b)
        gid = re.search(r"<ID>([^<]+)</ID>", b)
        if nm and gid:
            out[nm.group(1)] = gid.group(1)
    return out


def verify(path, aproj):
    x = rd(path)
    ET.fromstring(x.encode("utf-8", "surrogateescape"))
    pages = 0
    for m in re.finditer(r"<DiagramXml><!\[CDATA\[(.*?)\]\]></DiagramXml>", x, re.S):
        cd = m.group(1)
        ET.fromstring(("<r>" + cd + "</r>").encode("utf-8", "surrogateescape"))
        hdr = re.search(r'<AddFlow Nodes="(\d+)" Links="(\d+)"', cd)
        nn = len(re.findall(r"<Node\s+Left=", cd))
        ln = len(re.findall(r"<Link\s+Org=", cd))
        if hdr and (int(hdr.group(1)) != nn or int(hdr.group(2)) != ln):
            raise RuntimeError(f"{path}: page CDATA count mismatch")
        pages += 1
    print(f"verify {os.path.basename(path)}: XML OK, {pages} pages CDATA OK")


if __name__ == "__main__":
    sys.exit(main())
