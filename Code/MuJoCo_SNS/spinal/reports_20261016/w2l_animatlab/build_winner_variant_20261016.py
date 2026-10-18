r"""AnimatLab realization of the 2026-10-16 W2L robot campaign winner
(curr_s6_drive_robot_r060) on Walker_2_Layer_CPG_BilateralRG.

Overnight campaign stage "AnimatLab realizer" (2026-10-16). NO commits, NO
pushes; this writes ONLY into the campaign reports folder, never into
Ben's model files.

WHAT THE WINNER IS (from reports_20261016\w2l_wean\W2L_ROBOT_WEANLADDER_20261016.md
+ w2l_wean_best_params.json, study curr_s6_drive_robot_r060):
    te=2.4, tf=3.0 nA   (gamma = tonic drive into the four TONIC {L,R} RG
                          {ext,flx} ports; 60% of the bootstrap 4/5)
    alpha = 0             (alpha MN drive: present in the campaign as knobs,
                          the winner uses none)
    g_heel=0.58765        (heel afferent gain; Ben default 0.5 -> x1.1753)
    ctn_amp=0.71708       (contact force->current encoder amplitude; the
                          bootstrap used 2.0 -> x0.35854, saturating at
                          cref=46.37 N)
    comm=0.41304          (commissural scale: c1 inhibit + V3 weak excite +
                          their RG->c1/V3 drive edges, ALL x0.41304)
    (g_ib/ibnA NOT transferred: the campaign's own sensitivity probe at the
     winner shows the Ib pathway is NOT load-bearing there -- removal delta
     0 on the objective -- and in this aproj the Ib afferents share the
     "Afferent HC Excite" type with the Ia/II flexor afferents, which have
     no campaign analog; toe->TOEDF chain: no DF path exists in this
     AnimatLab net, same documented reduction family as the MuJoCo port.)

MAPPING (documented, ratio-based; the two stacks' gain scales differ):
    gamma te/tf  -> 4 tonic Stimuli on L/R RG ext (2.4 nA) + L/R RG flx
                    (3.0 nA). The stock 10 nA/10 ms kickoff stimuli stay
                    (the wean report: the kickoff is "part of the model,
                    unscaled").
    alpha        -> 4 Stimuli on L/R Hip MN ext + L/R Knee MN ext at 0 nA
                    (winner alpha=0; GUI-visible knobs for Ben).
    comm         -> SynAmp edits on "c1 Commissural Inhibit" (2.749 ->
                    1.1354) and "V3 Commissural Excite" (0.1 -> 0.0413), +
                    a NEW dedicated type "RG Excite (comm drive, campaign)"
                    (SynAmp 2.749x0.41304 = 1.1354) repointing the 4
                    RG->{c1,V3} input connexions (the shared "RG Excite"
                    type itself is untouched -- the RG half-centers keep
                    their stock wiring).
    ctn_amp      -> the 4 ContactCount adapters' Gain C. UNIT BUG FOUND:
                    the shipped asim has <C>20</C> in plain SI (AMPS per
                    contact) while every other adapter in the file uses
                    amp-scale values (e.g. 1.741e-009), and the aproj mirror
                    says Value="20" Scale="nano" Actual="1.741e-009"
                    (Value/Actual inconsistent). build_contact.pl evidently
                    intended 20 nA. We apply the campaign ratio to the
                    INTENDED value: 20 nA x 0.35854 = 7.17 nA ->
                    asim <C>7.17e-009</C>, aproj Value="7.17" nano
                    Actual="7.17e-009". V0 (control) keeps the shipped C=20.
    g_heel       -> NEW type "Afferent HC Excite (Contact, campaign)"
                    SynAmp 0.01 x 1.1753 = 0.01175, repointing ONLY the 20
                    contact-neuron->extensor-layer connexions (the shared
                    "Afferent HC Excite" stays for the 24 Ia/II/Ib afferent
                    links).

STRUCTURAL DELTA (not fixable tonight, documented for Ben): the AnimatLab
contact wiring is ipsilateral-extensor excitation only; Ben's rules crossed
heel->contralateral InF reset exists only in the MuJoCo M5 build.

OUTPUTS (all under --out):
    V0/<name>_V0_control.asim      untouched copy (faithful baseline)
    V1/<name>_V1_drive.asim        + gamma/alpha stimulators only
    V2/<name>_V2_winner.asim       + stimulators + comm + contact retune
    V3/<name>_V3_feedback.asim     comm + contact retune, NO tonic gamma
    V2/<name>_V2_winner.aproj      aproj mirror of V2 (GUI vehicle; NOT
                                   GUI-verified tonight -- labeled
                                   headless-run only, GUI open-check pending)
    edits_log.json                 every edit made, machine-readable

Usage (egg):
    D:\Anaconda\envs\myo\python.exe build_winner_variant_20261016.py ^
        --src D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG ^
        --out D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

SRC_ASIM = "Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim"
SRC_APROJ = "Walker_2_Layer_CPG_BilateralRG.aproj"

# winner values (w2l_wean_best_params.json / wean report)
TE, TF, ALPHA = 2.4, 3.0, 0.0
COMM = 0.41304
CTN_RATIO = 0.71708 / 2.0          # winner / bootstrap encoder amplitude
HEEL_RATIO = 0.58765 / 0.5         # winner / Ben-default heel afferent gain
CTN_DESIGN_NA = 20.0               # build_contact.pl's intended nA/count
SIM_END = "5.1"                    # keep the stock sim length

# neuron ID prefixes (8-hex, resolved per file; from build_rg.pl + audit)
PFX = dict(LRGext="17fdb1f6", LRGflx="7bd02f2b", RRGext="cafe0001",
           RRGflx="cafe0003", LHipMNExt="64083a38", LKneeMNExt="c0986b70",
           RHipMNExt="56e26621", RKneeMNExt="cabac9d1",
           Lc1="cafe0005", Rc1="cafe0006", LV3="cafe0007", RV3="cafe0008",
           Lheel="cafe0201", Ltoe="cafe0202", Rheel="cafe0203",
           Rtoe="cafe0204")
# contact-adapter Gain child IDs (the 4 ContactCount adapters, asim audit;
# the ADAPTER nodes are cafe0211-214, their Gain children cafe1308-1311)
ADAPTER_IDS = ["cafe1308", "cafe1309", "cafe1310", "cafe1311"]

C1_SYNAMP_OLD, C1_SYNAMP_NEW = "2.749", f"{2.749*COMM:.4f}"
V3_SYNAMP_OLD, V3_SYNAMP_NEW = "0.1", f"{0.1*COMM:.6f}"
COMMDRIVE_SYNAMP = f"{2.749*COMM:.4f}"
CONTACT_SYNAMP = f"{0.01*HEEL_RATIO:.6f}"
CTN_NA_NEW = f"{CTN_DESIGN_NA*CTN_RATIO:.4f}"           # 7.1708 nA
CTN_AMP_NEW = f"{CTN_DESIGN_NA*CTN_RATIO*1e-9:.6g}"     # 7.1708e-009 SI amps

# deterministic fresh GUID prefix (collision-checked at runtime)
NEWPFX = "d51a"


def g(n: int) -> str:
    return f"{NEWPFX}{n:04d}-0000-4000-8000-{900000000000+n}"


LOG: list = []


def log(msg: str):
    print(msg)
    LOG.append(msg)


def rd(p: str) -> str:
    with open(p, "r", encoding="utf-8", errors="surrogateescape", newline="") as f:
        return f.read()


def wr(p: str, s: str):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.write(s)


def resolve(x: str, pfx: str) -> str:
    i = x.find(f"<ID>{pfx}")
    if i < 0:
        raise RuntimeError(f"prefix {pfx} not found")
    j = x.find("</ID>", i + 4)
    return x[i + 4:j]


# ---------------------------------------------------------------- stimulators
def fmt_amp(na: float) -> str:
    """0 -> '0'; else nA in SI amps, '2.4e-009'-style like the stock file."""
    if na == 0:
        return "0"
    return f"{na:.6g}e-009"


def make_stim_asim(tmpl: str, sid: str, name: str, target: str, na: float) -> str:
    s = tmpl
    s = s.replace(f"<ID>{re.search(r'<ID>([^<]+)</ID>', tmpl).group(1)}</ID>",
                  f"<ID>{sid}</ID>", 1)
    s = re.sub(r"<Name>[^<]*</Name>", f"<Name>{name}</Name>", s, count=1)
    s = re.sub(r"<TargetNodeID>[^<]*</TargetNodeID>",
               f"<TargetNodeID>{target}</TargetNodeID>", s, count=1)
    cur = fmt_amp(na)
    s = re.sub(r"<CurrentOn>[^<]*</CurrentOn>", f"<CurrentOn>{cur}</CurrentOn>", s, count=1)
    for tag in ("EndTime", "CycleOnDuration", "BurstOnDuration"):
        s = re.sub(f"<{tag}>[^<]*</{tag}>", f"<{tag}>{SIM_END}</{tag}>", s, count=1)
    return s


def add_stimuli_asim(x: str, idmap: dict) -> str:
    i = x.find("<Name>Stimulus_1</Name>")
    if i < 0:
        raise RuntimeError("Stimulus_1 not found")
    s = x.rfind("<Stimulus>", 0, i)
    e = x.find("</Stimulus>", i) + len("</Stimulus>")
    tmpl = x[s:e]
    specs = [
        ("gamma_L_ext", "Gamma drive L RG ext (te=2.4nA)", "LRGext", TE),
        ("gamma_L_flx", "Gamma drive L RG flx (tf=3.0nA)", "LRGflx", TF),
        ("gamma_R_ext", "Gamma drive R RG ext (te=2.4nA)", "RRGext", TE),
        ("gamma_R_flx", "Gamma drive R RG flx (tf=3.0nA)", "RRGflx", TF),
        ("alpha_L_hip", "Alpha drive L Hip MN ext (0nA)", "LHipMNExt", ALPHA),
        ("alpha_L_knee", "Alpha drive L Knee MN ext (0nA)", "LKneeMNExt", ALPHA),
        ("alpha_R_hip", "Alpha drive R Hip MN ext (0nA)", "RHipMNExt", ALPHA),
        ("alpha_R_knee", "Alpha drive R Knee MN ext (0nA)", "RKneeMNExt", ALPHA),
    ]
    ins = ""
    for k, (_gk, nm, key, na) in enumerate(specs):
        blk = make_stim_asim(tmpl, g(100 + k), nm, idmap[key], na)
        ins += "\n" + blk
        log(f"asim: +Stimulus {nm} -> {key} @ {na} nA")
    close = x.find("</ExternalStimuli>", e)
    return x[:close] + ins + x[close:]


# --------------------------------------------------------------------- types
def clone_type_asim(x: str, type_name: str, new_id: str, new_name: str,
                    new_synamp: str) -> str:
    """Append a cloned SynapseType with new name/id/SynAmp; returns (x, blk)."""
    i = x.find(f"<Name>{type_name}</Name>")
    if i < 0:
        raise RuntimeError(f"type {type_name} not found")
    s = x.rfind("<SynapseType>", 0, i)
    e = x.find("</SynapseType>", i) + len("</SynapseType>")
    blk = x[s:e]
    old_id = re.search(r"<ID>([^<]+)</ID>", blk).group(1)
    nb = blk.replace(f"<ID>{old_id}</ID>", f"<ID>{new_id}</ID>", 1)
    nb = re.sub(r"<Name>[^<]*</Name>", f"<Name>{new_name}</Name>", nb, count=1)
    nb = re.sub(r"<SynAmp>[^<]*</SynAmp>", f"<SynAmp>{new_synamp}</SynAmp>", nb, count=1)
    return x[:e] + "\n" + nb + x[e:], nb


def set_type_synamp(x: str, type_name: str, old: str, new: str) -> str:
    i = x.find(f"<Name>{type_name}</Name>")
    if i < 0:
        raise RuntimeError(f"type {type_name} not found")
    s = x.rfind("<SynapseType>", 0, i)
    e = x.find("</SynapseType>", i)
    blk = x[s:e]
    if f"<SynAmp>{old}</SynAmp>" not in blk:
        raise RuntimeError(f"{type_name}: SynAmp {old} not found in block")
    nb = blk.replace(f"<SynAmp>{old}</SynAmp>", f"<SynAmp>{new}</SynAmp>", 1)
    log(f"asim: type '{type_name}' SynAmp {old} -> {new}")
    return x[:s] + nb + x[s + len(blk):]


def repoint_connexions(x: str, pairs: list, new_type_id: str, what: str) -> str:
    """pairs = [(src_full_id, dst_full_id)]; sets SynapseTypeID on each."""
    n = 0
    out = []
    pos = 0
    for m in re.finditer(r"<Connexion>(.*?)</Connexion>", x, re.S):
        b = m.group(1)
        src = re.search(r"<SourceID>([^<]*)</SourceID>", b)
        dst = re.search(r"<TargetID>([^<]*)</TargetID>", b)
        if src and dst and (src.group(1), dst.group(1)) in pairs:
            nb = re.sub(r"<SynapseTypeID>[^<]*</SynapseTypeID>",
                        f"<SynapseTypeID>{new_type_id}</SynapseTypeID>", b, count=1)
            out.append(x[pos:m.start()])
            out.append(f"<Connexion>{nb}</Connexion>")
            pos = m.end()
            n += 1
    out.append(x[pos:])
    log(f"asim: repointed {n}/{len(pairs)} connexions ({what})")
    if n != len(pairs):
        raise RuntimeError(f"expected {len(pairs)} repoints, got {n}")
    return "".join(out)


def set_adapter_gain_c(x: str, adapter_ids: list, new_c: str) -> str:
    for aid in adapter_ids:
        i = x.find(f"<ID>{aid}")
        if i < 0:
            raise RuntimeError(f"adapter {aid} not found")
        s = x.rfind("<Adapter>", 0, i)
        e = x.find("</Adapter>", i)
        blk = x[s:e]
        gb = re.search(r"<Gain>.*?</Gain>", blk, re.S).group(0)
        if "<C>20</C>" not in gb:
            raise RuntimeError(f"adapter {aid}: <C>20</C> not found")
        ngb = gb.replace("<C>20</C>", f"<C>{new_c}</C>", 1)
        x = x[:s] + blk.replace(gb, ngb, 1) + x[e:]
        log(f"asim: adapter {aid} Gain C 20 -> {new_c}")
    return x


# ------------------------------------------------------------------ asim main
def build_asim(src_path: str, out_path: str, want_stim: bool, want_fb: bool):
    x = rd(src_path)
    if NEWPFX in x:
        raise RuntimeError(f"GUID prefix {NEWPFX} already present -- refusing")
    idmap = {k: resolve(x, p) for k, p in PFX.items()}

    if want_stim:
        x = add_stimuli_asim(x, idmap)

    if want_fb:
        # 1) scale the two dedicated commissural output types
        x = set_type_synamp(x, "c1 Commissural Inhibit", C1_SYNAMP_OLD, C1_SYNAMP_NEW)
        x = set_type_synamp(x, "V3 Commissural Excite", V3_SYNAMP_OLD, V3_SYNAMP_NEW)
        # 2) new type for the RG->{c1,V3} drive edges, repoint those 4 edges
        x, _ = clone_type_asim(x, "RG Excite", g(1),
                               "RG Excite (comm drive, campaign)", COMMDRIVE_SYNAMP)
        comm_pairs = [(idmap["LRGflx"], idmap["Lc1"]), (idmap["RRGflx"], idmap["Rc1"]),
                      (idmap["LRGext"], idmap["LV3"]), (idmap["RRGext"], idmap["RV3"])]
        x = repoint_connexions(x, comm_pairs, g(1), "RG->c1/V3 comm drive")
        # 3) new type for the 20 contact edges, repoint them
        x, _ = clone_type_asim(x, "Afferent HC Excite", g(2),
                               "Afferent HC Excite (Contact, campaign)", CONTACT_SYNAMP)
        contact_srcs = {idmap[k] for k in ("Lheel", "Ltoe", "Rheel", "Rtoe")}
        contact_pairs = []
        for m in re.finditer(r"<Connexion>(.*?)</Connexion>", x, re.S):
            b = m.group(1)
            src = re.search(r"<SourceID>([^<]*)</SourceID>", b)
            if src and src.group(1) in contact_srcs:
                dst = re.search(r"<TargetID>([^<]*)</TargetID>", b).group(1)
                contact_pairs.append((src.group(1), dst))
        x = repoint_connexions(x, contact_pairs, g(2), "contact->extensor layers")
        # 4) contact encoder amplitude (design-intent nA; see module docstring)
        x = set_adapter_gain_c(x, ADAPTER_IDS, CTN_AMP_NEW)

    wr(out_path, x)
    return x


# ----------------------------------------------------------------- aproj side
def build_aproj(src_path: str, out_path: str):
    p = rd(src_path)
    if NEWPFX in p:
        raise RuntimeError(f"GUID prefix {NEWPFX} already present -- refusing")
    idmap = {k: resolve(p, pfx) for k, pfx in PFX.items()}

    # ---- 1) stimuli (clone the aproj Stimulus_1 block) ----
    i = p.find("<Name>Stimulus_1</Name>")
    s = p.rfind("<Stimulus>", 0, i)
    e = p.find("</Stimulus>", i) + len("</Stimulus>")
    tmpl = p[s:e]
    specs = [
        ("Gamma drive L RG ext (te=2.4nA)", "LRGext", TE),
        ("Gamma drive L RG flx (tf=3.0nA)", "LRGflx", TF),
        ("Gamma drive R RG ext (te=2.4nA)", "RRGext", TE),
        ("Gamma drive R RG flx (tf=3.0nA)", "RRGflx", TF),
        ("Alpha drive L Hip MN ext (0nA)", "LHipMNExt", ALPHA),
        ("Alpha drive L Knee MN ext (0nA)", "LKneeMNExt", ALPHA),
        ("Alpha drive R Hip MN ext (0nA)", "RHipMNExt", ALPHA),
        ("Alpha drive R Knee MN ext (0nA)", "RKneeMNExt", ALPHA),
    ]
    ins = ""
    for k, (nm, key, na) in enumerate(specs):
        blk = tmpl
        old_id = re.search(r"<ID>([^<]+)</ID>", blk).group(1)
        blk = blk.replace(f"<ID>{old_id}</ID>", f"<ID>{g(200+k)}</ID>", 1)
        blk = re.sub(r"<Name>[^<]*</Name>", f"<Name>{nm}</Name>", blk, count=1)
        blk = re.sub(r"<NodeID>[^<]*</NodeID>", f"<NodeID>{idmap[key]}</NodeID>",
                     blk, count=1)
        actual = f"{na:.6g}e-009" if na else "0"
        blk = re.sub(r'<CurrentOn Value="[^"]*" Scale="nano" Actual="[^"]*"/>',
                     f'<CurrentOn Value="{na:.6g}" Scale="nano" Actual="{actual}"/>',
                     blk, count=1)
        blk = re.sub(r'<EndTime Value="[^"]*" Scale="None" Actual="[^"]*"/>',
                     f'<EndTime Value="{SIM_END}" Scale="None" Actual="{SIM_END}"/>',
                     blk, count=1)
        for tag in ("CycleOnDuration", "BurstOnDuration"):
            blk = re.sub(rf'<{tag} Value="[^"]*" Scale="milli" Actual="[^"]*"/>',
                         f'<{tag} Value="{float(SIM_END)*1000:.6g}" Scale="milli" '
                         f'Actual="{SIM_END}"/>', blk, count=1)
        ins += "\n" + blk
        log(f"aproj: +Stimulus {nm} -> {key} @ {na} nA")
    # insert right after the template Stimulus_1 block (order is irrelevant;
    # avoids guessing the aproj stimuli container element name)
    p = p[:e] + ins + p[e:]

    # ---- 2) comm output type conductances ----
    def patch_type_link(p, name, new_val, new_actual):
        for m in re.finditer(r"<Link>(.*?)</Link>", p, re.S):
            b = m.group(1)
            if "SynapseTypes." not in b:
                continue
            nm = re.search(r"<Name>([^<]*)</Name>", b)
            if nm and nm.group(1) == name:
                nb = re.sub(r'<MaxSynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>',
                            f'<MaxSynapticConductance Value="{new_val}" Scale="micro" '
                            f'Actual="{new_actual}"/>', b, count=1)
                if nb == b:
                    raise RuntimeError(f"aproj type {name}: conductance patch failed")
                return p[:m.start()] + f"<Link>{nb}</Link>" + p[m.end():]
        raise RuntimeError(f"aproj type link {name} not found")

    p = patch_type_link(p, "c1 Commissural Inhibit", C1_SYNAMP_NEW, f"{2.749*COMM:.4g}e-007")
    p = patch_type_link(p, "V3 Commissural Excite", V3_SYNAMP_NEW, f"{0.1*COMM:.6g}e-008")
    log("aproj: c1/V3 type MaxSynapticConductance scaled by comm")

    # ---- 3) new aproj type links + repoint 4 comm links + 20 contact links ----
    def clone_type_link(p, name, new_id, new_name, new_val, new_actual):
        for m in re.finditer(r"<Link>(.*?)</Link>", p, re.S):
            b = m.group(1)
            if "SynapseTypes." not in b:
                continue
            nm = re.search(r"<Name>([^<]*)</Name>", b)
            if nm and nm.group(1) == name:
                old_id = re.search(r"<ID>([^<]+)</ID>", b).group(1)
                nb = b.replace(f"<ID>{old_id}</ID>", f"<ID>{new_id}</ID>", 1)
                nb = re.sub(r"<Name>[^<]*</Name>", f"<Name>{new_name}</Name>", nb, count=1)
                nb = re.sub(r'<MaxSynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>',
                            f'<MaxSynapticConductance Value="{new_val}" Scale="micro" '
                            f'Actual="{new_actual}"/>', nb, count=1)
                at = m.end()
                return p[:at] + f"\n<Link>{nb}</Link>" + p[at:]
        raise RuntimeError(f"aproj type link {name} not found for clone")

    p = clone_type_link(p, "RG Excite", g(1), "RG Excite (comm drive, campaign)",
                        COMMDRIVE_SYNAMP, f"{2.749*COMM:.4g}e-007")
    p = clone_type_link(p, "Afferent HC Excite", g(2),
                        "Afferent HC Excite (Contact, campaign)",
                        CONTACT_SYNAMP, f"{0.01*HEEL_RATIO:.6g}e-008")

    def repoint_links(p, pairs, new_tid, new_cond, new_actual, what):
        n = 0
        out, pos = [], 0
        for m in re.finditer(r"<Link>(.*?)</Link>", p, re.S):
            b = m.group(1)
            if "Behavior.Synapse" not in b:
                continue
            o = re.search(r"<OriginID>([^<]*)</OriginID>", b)
            d = re.search(r"<DestinationID>([^<]*)</DestinationID>", b)
            if not (o and d):
                continue
            if (o.group(1), d.group(1)) in pairs:
                nb = re.sub(r"<SynapticTypeID>[^<]*</SynapticTypeID>",
                            f"<SynapticTypeID>{new_tid}</SynapticTypeID>", b, count=1)
                nb = re.sub(r'<SynapticConductance Value="[^"]*" Scale="micro" Actual="[^"]*"/>',
                            f'<SynapticConductance Value="{new_cond}" Scale="micro" '
                            f'Actual="{new_actual}"/>', nb, count=1)
                out.append(p[pos:m.start()])
                out.append(f"<Link>{nb}</Link>")
                pos = m.end()
                n += 1
        out.append(p[pos:])
        log(f"aproj: repointed {n}/{len(pairs)} links ({what})")
        if n != len(pairs):
            raise RuntimeError(f"aproj: expected {len(pairs)} repoints, got {n}")
        return "".join(out)

    comm_pairs = [(idmap["LRGflx"], idmap["Lc1"]), (idmap["RRGflx"], idmap["Rc1"]),
                  (idmap["LRGext"], idmap["LV3"]), (idmap["RRGext"], idmap["RV3"])]
    p = repoint_links(p, comm_pairs, g(1), COMMDRIVE_SYNAMP, f"{2.749*COMM:.4g}e-007",
                      "RG->c1/V3 comm drive")
    contact_srcs = {idmap[k] for k in ("Lheel", "Ltoe", "Rheel", "Rtoe")}
    contact_pairs = []
    for m in re.finditer(r"<Link>(.*?)</Link>", p, re.S):
        b = m.group(1)
        if "Behavior.Synapse" not in b:
            continue
        o = re.search(r"<OriginID>([^<]*)</OriginID>", b)
        if o and o.group(1) in contact_srcs:
            d = re.search(r"<DestinationID>([^<]*)</DestinationID>", b)
            contact_pairs.append((o.group(1), d.group(1)))
    p = repoint_links(p, contact_pairs, g(2), CONTACT_SYNAMP,
                      f"{0.01*HEEL_RATIO:.6g}e-008", "contact->extensor layers")

    # ---- 4) adapter Gain C (design-intent nA); SINGLE-PASS rebuild
    #      (mutating p inside a finditer loop invalidates later offsets --
    #       that corruption is exactly what XML verify is here to catch)
    outp_, pos_, n_ad = [], 0, 0
    for m in re.finditer(r"<Node>(.*?)</Node>", p, re.S):
        b = m.group(1)
        if "PhysicalToNodeAdapter" not in b or "<DataTypeID>ContactCount" not in b:
            continue
        gb = re.search(r"<Gain>.*?</Gain>", b, re.S)
        if not gb or 'Value="20"' not in gb.group(0):
            continue
        ngb = gb.group(0).replace(
            '<C Value="20" Scale="nano" Actual="1.741e-009"/>',
            f'<C Value="{CTN_NA_NEW}" Scale="nano" Actual="{CTN_AMP_NEW}"/>', 1)
        if ngb == gb.group(0):
            raise RuntimeError("aproj contact adapter C patch failed (format drift)")
        outp_.append(p[pos_:m.start()])
        outp_.append(f"<Node>{b.replace(gb.group(0), ngb, 1)}</Node>")
        pos_ = m.end()
        n_ad += 1
    outp_.append(p[pos_:])
    p = "".join(outp_)
    log(f"aproj: contact adapter Gain C patched on {n_ad} nodes")
    if n_ad != 4:
        raise RuntimeError(f"expected 4 contact adapters, patched {n_ad}")

    wr(out_path, p)
    return p


# --------------------------------------------------------------- verification
def verify_file(path: str, is_aproj: bool):
    x = rd(path)
    ET.fromstring(x.encode("utf-8", "surrogateescape"))
    # per-page CDATA sanity (untouched, but confirm not corrupted)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    base = "Walker_2_Layer_CPG_BilateralRG"
    variants = [
        ("V0", f"{base}_V0_control.asim", False, False),
        ("V1", f"{base}_V1_drive.asim", True, False),
        ("V2", f"{base}_V2_winner.asim", True, True),
        ("V3", f"{base}_V3_feedback.asim", False, True),
    ]
    for tag, fname, stim, fb in variants:
        log(f"===== variant {tag} (stim={stim}, feedback={fb}) =====")
        outp = os.path.join(a.out, tag, fname)
        build_asim(os.path.join(a.src, SRC_ASIM), outp, stim, fb)
        verify_file(outp, False)
        log(f"wrote {outp}")

    log("===== variant V2 aproj mirror =====")
    aproj_out = os.path.join(a.out, "V2", f"{base}_V2_winner.aproj")
    build_aproj(os.path.join(a.src, SRC_APROJ), aproj_out)
    verify_file(aproj_out, True)
    log(f"wrote {aproj_out}")

    with open(os.path.join(a.out, "edits_log.json"), "w", encoding="utf-8") as f:
        json.dump({"log": LOG,
                   "knobs": dict(te=TE, tf=TF, alpha=ALPHA, comm=COMM,
                                 ctn_ratio=CTN_RATIO, heel_ratio=HEEL_RATIO,
                                 ctn_design_nA=float(CTN_NA_NEW),
                                 c1_synamp=C1_SYNAMP_NEW, v3_synamp=V3_SYNAMP_NEW,
                                 commdrive_synamp=COMMDRIVE_SYNAMP,
                                 contact_synamp=CONTACT_SYNAMP)}, f, indent=2)
    print("ALL VARIANTS BUILT + VERIFIED (XML + pages)")


if __name__ == "__main__":
    sys.exit(main())
