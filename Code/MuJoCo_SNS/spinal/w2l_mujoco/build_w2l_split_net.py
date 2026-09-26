"""Milestone 4: the SPLIT-RG W2L architecture — Ben's ask, verbatim:

  "DISCONNECT the LH RG from directly controlling the RH PF layers. Create a
   copy of the LH RG for the RH side and use it to power the RH PF layers;
   then couple the two RGs with the commissural wiring of the documented
   bilateralrg build."

This is exactly SESSION_NOTES_20260916.md build-chain steps 1-3
(Neuromechanical_Models\\Walker_2_Layer_CPG_BilateralRG\\tools\\), applied to
the M3 2023-original transcription (build_w2l_orig_net.py, w2laproj template):

  step 1 (build_rg.pl)    remove the 4 CROSSED pf_drive edges (L RG ext ->
                          R Knee/R Hip PF *flx*, L RG flx -> R Hip/R Knee PF
                          *ext*) = the L-RG->R-PF disconnect; add the R RG
                          half-center block mirroring the L block exactly
                          (SESSION_NOTES line 85) and 4 IPSILATERAL R pf_drive
                          edges; Stimulus_2 (10 nA, 0-10 ms, into R RG flx =
                          antiphase kickoff) becomes an input port.
  steps 2-3 (build_comm.pl + patch_comm_types.pl) the Shinohara commissurals
                          between the two RGs: S RG flx -> S c1 -(inh)- contra
                          RG flx; S RG ext -> S V3 -(exc, weak)- contra RG ext
                          AND contra RG ext IN (both V3 out-edges, as in the
                          bilateralrg template). Built ONLY when the coupling
                          scale > 0 (conditional-topology lab law: zero-gain
                          synapses alone perturb results via BLAS summation
                          order) — so --comm=0 is a true ablation, not a
                          zero-gain one.

GAINS: "the gains that build used" = the w2l_cpg realization of this same
documented build (spinal\\w2l_cpg\\build_w2l_net.py W2L_GAINS), which re-express
the AnimatLab SynAmps (c1 2.749 inh / V3 0.1 weak) on this stack's 0..5 mV
scale with the measured anti-latch fix: comm_c1=1.0, c1_inh=3.0 (locks the two
sides' periods; 2.0 lets them drift, measured 2026-09-24), comm_v3=1.0,
v3_weak=0.08 (full-strength V3 E-latches both RGs, measured). R RG mirror and
ipsilateral pf_drive use the same L-side gains (rg_direct_exc 0.5 /
rg_laminate 4.0 / pf_drive 2.4) — the mirror is a COPY of the L RG.

CAUSAL DISCONNECT ASSERTED: after the edit, no edge into any "R Hip PF*" /
"R Knee PF*" node has an L-side source; every cross-side edge in the coupled
net is RG<->RG commissural (c1/V3) or the shared knee-extensor Renshaw quirk
transcribed from the 2023 original (MN/RC level, not RG->PF).

CENSUS (asserted):
  coupled  : 87 neurons / 166 synapses / 6 inputs / 12 outputs
  --comm=0 : 83 neurons / 156 synapses / 6 inputs / 12 outputs
  (M3 baseline: 79 / 150 / 3 / 12; delta = +4 R RG +4 comm neurons,
   -4 crossed pf_drive +6 R rg_laminate +4 R pf_drive +10 comm synapses)

USAGE
    net = build(comm=1.0, tau_rg_nap_h=0.25)
    u = net.make_inputs()                       # len 6
    u[net.input_index("Stimulus_1")] = 10.0     # nA, t in [0, 0.01] (L RG ext)
    u[net.input_index("Stimulus_2")] = 10.0     # nA, t in [0, 0.01] (R RG flx)
    u[net.input_index("TONIC L RG ext")] = 3.0  # nA (te)
    ... V = net.step(u)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)

HERE = Path(__file__).parent
SPINAL = HERE.parent
if str(SPINAL) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(SPINAL))
from build_network import SNS_NumpyFixedTau  # noqa: E402  (read-only import)
from w2l_cpg.build_w2l_net import DT, NAP, TAU, W2L_GAINS, _neu, _syn  # noqa: E402
import build_w2l_orig_net as ORIG  # noqa: E402  (M3 transcription we split)

E_HI = ORIG.E_HI
MUSCLE_MAP = ORIG.MUSCLE_MAP

# ---------------------------------------------------------------- edit spec
# The 4 CROSSED pf_drive edges removed from the w2laproj template (the
# L-RG->R-PF disconnect). (from, to) verbatim from the template.
CROSSED_REMOVE = [
    ("L RG ext", "R Knee PF flx"),
    ("L RG ext", "R Hip PF flx"),
    ("L RG flx", "R Hip PF ext"),
    ("L RG flx", "R Knee PF ext"),
]

# R RG half-center block: exact mirror of the L RG block's 6 rg_laminate-family
# edges (SESSION_NOTES line 85 "mirrors the L RG block exactly"; the bilateralrg
# template ships only the 4 laminated ones — the 2 direct HC<->HC escape edges
# are the documented MIRROR FIX, w2l_cpg/build_w2l_net.py:321-335).
R_MIRROR_EDGES = [
    ("R RG flx", "R RG ext", "exc", "rg_laminate"),   # direct HC<->HC escape
    ("R RG ext", "R RG flx", "exc", "rg_laminate"),   # direct HC<->HC escape
    ("R RG ext", "R RG ext IN", "exc", "rg_laminate"),
    ("R RG flx", "R RG flx IN", "exc", "rg_laminate"),
    ("R RG flx IN", "R RG ext", "inh", "rg_laminate"),
    ("R RG ext IN", "R RG flx", "inh", "rg_laminate"),
]

# IPSILATERAL R PF drive: R RG now powers the R PF layers (the ask).
R_PF_DRIVE_EDGES = [
    ("R RG ext", "R Hip PF ext", "exc", "pf_drive"),
    ("R RG ext", "R Knee PF ext", "exc", "pf_drive"),
    ("R RG flx", "R Hip PF flx", "exc", "pf_drive"),
    ("R RG flx", "R Knee PF flx", "exc", "pf_drive"),
]

# Shinohara commissurals, verbatim (from, to, sign, tag) from the bilateralrg
# template (steps 2-3). Both V3 out-edges per side (contra RG ext + contra
# RG ext IN) as the template carries them.
COMM_EDGES = [
    ("L RG flx", "c1_L", "exc", "comm_c1"),
    ("c1_L", "R RG flx", "inh", "c1_SynAmp2.749"),
    ("R RG flx", "c1_R", "exc", "comm_c1"),
    ("c1_R", "L RG flx", "inh", "c1_SynAmp2.749"),
    ("L RG ext", "V3_L", "exc", "comm_v3"),
    ("V3_L", "R RG ext", "exc", "v3_SynAmp0.1_weak"),
    ("V3_L", "R RG ext IN", "exc", "v3_to_contra_InE"),
    ("R RG ext", "V3_R", "exc", "comm_v3"),
    ("V3_R", "L RG ext", "exc", "v3_SynAmp0.1_weak"),
    ("V3_R", "L RG ext IN", "exc", "v3_to_contra_InE"),
]

COMM_NODES = [("c1_L", "IN-C"), ("c1_R", "IN-C"),
              ("V3_L", "IN-V3"), ("V3_R", "IN-V3")]

RG_MIRROR_NODES = [("R RG ext", "HC-RG-E"), ("R RG flx", "HC-RG-F"),
                   ("R RG ext IN", "IN-InE"), ("R RG flx IN", "IN-InF")]

TONIC_PORTS = ("L RG ext", "L RG flx", "R RG ext", "R RG flx")


def make_split_template(tpl: dict, comm: float = 1.0) -> dict:
    """w2laproj template -> split-RG template (remove crossed pf_drive, add
    R RG mirror + ipsilateral R pf_drive + optional commissurals).

    comm <= 0 -> the commissural neurons AND synapses are NOT built
    (conditional topology: a zero-gain synapse still changes summation order)."""
    nodes = [dict(n) for n in tpl["nodes"]]
    labels = {n["label"] for n in nodes}
    edges = [dict(e) for e in tpl["edges"]]

    # ---- 1. remove the 4 crossed pf_drive edges
    removed = [e for e in edges if (e["from"], e["to"]) in
               [(a, b) for a, b in CROSSED_REMOVE]]
    assert len(removed) == 4, f"expected 4 crossed pf_drive edges, found {len(removed)}"
    for e in removed:
        assert e["tag"] == "pf_drive" and e["sign"] == "exc", \
            f"removal spec mismatch: {e}"
    edges = [e for e in edges if e not in removed]

    # ---- 2. add R RG mirror + ipsilateral R pf_drive + the antiphase
    # kickoff port (SESSION_NOTES step 1: Stimulus_2, 10 nA 0-10 ms, into
    # R RG flx -> "legs start neurally antiphase"; the waveform lives in the
    # caller, as for Stimulus_1).
    for lab, typ in RG_MIRROR_NODES:
        assert lab not in labels, f"{lab} already in template"
        nodes.append({"label": lab, "type": typ})
        labels.add(lab)
    edges += [{"from": f, "to": t, "sign": s, "tag": g}
              for f, t, s, g in R_MIRROR_EDGES + R_PF_DRIVE_EDGES]
    nodes.append({"label": "Stimulus_2", "type": "PORT-load"})
    labels.add("Stimulus_2")
    edges.append({"from": "Stimulus_2", "to": "R RG flx", "sign": "exc",
                  "tag": "kickoff_antiphase"})

    # ---- 3. commissurals (coupled only)
    if comm > 0:
        for lab, typ in COMM_NODES:
            assert lab not in labels, f"{lab} already in template"
            nodes.append({"label": lab, "type": typ})
            labels.add(lab)
        edges += [{"from": f, "to": t, "sign": s, "tag": g}
                  for f, t, s, g in COMM_EDGES]

    # ---- CAUSAL DISCONNECT ASSERT: no edge into an R-side PF half-center has
    # an L-side source. (R PF INs are fed by R PFs; R PFs are fed by R RG and
    # their own side's INs only.)
    for e in edges:
        if e["to"].startswith(("R Hip PF", "R Knee PF")) and \
                not e["to"].endswith(" IN"):
            assert not e["from"].startswith("L "), \
                f"L->R PF edge survived the disconnect: {e}"

    return {"nodes": nodes, "edges": edges}


@dataclass
class W2LSplitNet:
    """Compiled wrapper over the split-RG network (same interface as
    W2LOrigNet: idx / inputs / make_inputs / step / muscle_ctrl)."""
    template: dict
    comm: float = 1.0
    gains: dict = field(default_factory=dict)   # merged over W2L_GAINS
    net: object = field(init=False, default=None)
    compiled: object = field(init=False, default=None)
    idx: dict = field(default_factory=dict)
    inputs: list = field(default_factory=list)
    port_targets: dict = field(default_factory=dict)
    muscle_outputs: dict = field(default_factory=dict)
    synapse_counts: dict = field(default_factory=dict)
    n_synapses: int = field(init=False, default=0)
    removed_edges: list = field(default_factory=list)

    def __post_init__(self):
        self.G = dict(W2L_GAINS)
        self.G.update(self.gains or {})
        tpl = self.template
        nodes = {n["label"]: n for n in tpl["nodes"]}
        edges = tpl["edges"]
        self.net = ORIG.Network(name="W2L split-RG (M4, commissural bilateral)")

        for e in edges:
            for k in ("from", "to"):
                if e[k] not in nodes:
                    raise ValueError(f"template edge endpoint unknown: {e}")

        # ---- neurons (identical logic to M3; + IN-C/IN-V3 -> rg tau)
        for nd in tpl["nodes"]:
            lab, typ = nd["label"], nd["type"]
            if typ in ("MUSCLE", "PORT-load"):
                continue
            tau_key = _TYPE_TAU[typ]
            if typ in ("HC-RG-E", "HC-RG-F"):
                nap = NonSpikingNeuronWithPersistentSodiumChannel(
                    membrane_capacitance=TAU["rg"], membrane_conductance=1.0,
                    resting_potential=0.0, bias=0.0,
                    g_ion=np.array([NAP["g_ion"]]),
                    e_ion=np.array([NAP["e_ion"]]),
                    k_m=np.array([NAP["k_m"]]),
                    slope_m=np.array([NAP["slope_m"]]),
                    e_m=np.array([NAP["e_m"]]),
                    k_h=np.array([NAP["k_h"]]),
                    slope_h=np.array([NAP["slope_h"]]),
                    e_h=np.array([NAP["e_h"]]),
                    tau_max_h=np.array([NAP["tau_max_h"]]))
                self.net.add_neuron(nap, name=lab)
            else:
                self.net.add_neuron(_neu(TAU[tau_key]), name=lab)
            self.idx[lab] = len(self.idx)

        # ---- input ports: template PORT-load stimuli (Stimulus_1 -> L RG ext,
        # Stimulus_2 -> R RG flx) + TONIC on all four RG half-centers
        outgoing = {}
        for e in edges:
            outgoing.setdefault(e["from"], []).append(e)
        for nd in tpl["nodes"]:
            lab, typ = nd["label"], nd["type"]
            if typ == "PORT-load":
                outs = outgoing.get(lab, [])
                if len(outs) != 1:
                    raise ValueError(f"{lab}: expected 1 outgoing edge")
                tgt = outs[0]["to"]
                self.net.add_input(tgt)
                self.inputs.append(lab)
                self.port_targets[lab] = tgt
        for lab in TONIC_PORTS:
            port = "TONIC " + lab
            self.net.add_input(lab)
            self.inputs.append(port)
            self.port_targets[port] = lab

        # ---- muscles: output ports of their MN (pass-through, g=1)
        for e in edges:
            if e["tag"] == "mn_to_muscle":
                self.net.add_output(e["from"], name=e["to"])
                self.muscle_outputs[MUSCLE_MAP[e["to"]]] = e["from"]

        # ---- synapses: every remaining chemical edge, per-tag conductance
        n_types = {lab: nd["type"] for lab, nd in nodes.items()}
        for e in edges:
            tag, sign = e["tag"], e["sign"]
            if tag == "mn_to_muscle":
                continue
            if n_types[e["from"]] == "PORT-load":
                continue
            g = self._gain(tag, sign, e, n_types)
            self.net.add_connection(_syn(g, exc=(sign == "exc")),
                                    e["from"], e["to"])
            self.synapse_counts[tag] = self.synapse_counts.get(tag, 0) + 1
            self.n_synapses += 1

        self._assert_complete(edges, n_types)

    def _gain(self, tag: str, sign: str, e: dict, n_types: dict) -> float:
        G = self.G
        if tag == "rg_laminate":
            if n_types[e["from"]] in ("HC-RG-E", "HC-RG-F") and \
                    n_types[e["to"]] in ("HC-RG-E", "HC-RG-F"):
                return G["rg_direct_exc"]
            return G["rg_laminate"]
        if tag == "pf_drive":
            return G["pf_drive"]
        if tag == "pf_cross":
            return G["pf_cross"]
        if tag == "pf_to_mn":
            return G["pf_to_mn"]
        if tag == "rc_own":
            return G["rc_own_exc"] if sign == "exc" else G["rc_inh"]
        if tag == "rc_mutual":
            return G["rc_inh"]
        if tag == "ia_recip":
            return G["ia_sn_exc"] if sign == "exc" else G["ia_inh"]
        if tag == "ia_mutual":
            return G["ia_mutual"]
        if tag == "ib_auto":
            return G["ib_auto"]
        if tag == "other":
            return G["other"]
        # ---- commissural tags (bilateralrg verbatim, w2l_cpg gain names),
        # scaled by self.comm (ablation knob; comm=0 never reaches here —
        # make_split_template dropped the nodes/edges entirely)
        if tag == "comm_c1":
            return G["comm_c1"] * self.comm
        if tag == "c1_SynAmp2.749":
            return G["c1_inh"] * self.comm
        if tag == "comm_v3":
            return G["comm_v3"] * self.comm
        if tag in ("v3_SynAmp0.1_weak", "v3_to_contra_InE"):
            return G["v3_weak"] * self.comm
        raise ValueError(f"unknown template tag: {tag}")

    def _assert_complete(self, edges, n_types):
        expected = {}
        for e in edges:
            if e["tag"] == "mn_to_muscle":
                continue
            if n_types[e["from"]] == "PORT-load":
                continue
            expected[e["tag"]] = expected.get(e["tag"], 0) + 1
        comm = self.comm > 0
        n_neu, n_syn = (87, 166) if comm else (83, 156)
        assert len(self.idx) == n_neu, \
            f"neuron count {len(self.idx)} != {n_neu}"
        assert self.n_synapses == sum(expected.values()) == n_syn, \
            f"synapse transcription mismatch: {self.n_synapses} vs {n_syn}"
        assert dict(self.synapse_counts) == expected, "per-tag mismatch"
        assert len(self.inputs) == 6, f"port count {len(self.inputs)} != 6"
        assert len(self.muscle_outputs) == 12, "muscle outputs != 12"
        # the 4 crossed edges must be gone and the mirror/comm families present
        tags = expected
        assert tags.get("pf_drive") == 8, tags.get("pf_drive")  # 4 L + 4 R ipsi
        assert tags.get("rg_laminate") == 12, tags.get("rg_laminate")  # 6+6
        if comm:
            assert tags.get("comm_c1") == 2 and tags.get("c1_SynAmp2.749") == 2
            assert tags.get("comm_v3") == 2 and \
                tags.get("v3_SynAmp0.1_weak") == 2 and \
                tags.get("v3_to_contra_InE") == 2
        else:
            assert not any(k.startswith(("comm_", "c1_", "v3_"))
                           for k in tags)

    # ------------------------------------------------------------------ run
    def compile(self, dt: float = DT):
        self.compiled = self.net.compile(dt=dt, backend="numpy")
        self.compiled.__class__ = SNS_NumpyFixedTau
        assert self.compiled.V.shape[0] == len(self.idx), "index mismatch"
        return self.compiled

    def step(self, input_currents: np.ndarray) -> np.ndarray:
        self.compiled.forward(list(input_currents))
        return self.compiled.V

    def input_index(self, port: str) -> int:
        return self.inputs.index(port)

    def make_inputs(self) -> np.ndarray:
        return np.zeros(len(self.inputs))

    def muscle_ctrl(self, V: np.ndarray) -> dict:
        return {a: float(min(max(V[self.idx[mn]] / E_HI, 0.0), 1.0))
                for a, mn in self.muscle_outputs.items()}


# template node type -> membrane tau key (M3 table + the commissural INs,
# which live in the RG layer: w2l_cpg TAU maps IN-C/IN-V3 -> "rg")
_TYPE_TAU = dict(ORIG._TYPE_TAU, **{"IN-C": "rg", "IN-V3": "rg"})


def build(template_path=None, dt: float = DT, gains: dict | None = None,
          tau_rg_nap_h: float | None = None, comm: float = 1.0) -> W2LSplitNet:
    """Build + compile the split-RG network.

    comm: commissural coupling scale. 1.0 = the documented bilateralrg gains
    (W2L_GAINS comm_c1/c1_inh/comm_v3/v3_weak); 0 = ABLATION (the 4 c1/V3
    neurons and all 10 commissural synapses are not built — conditional
    topology, not zero-gain edges). tau_rg_nap_h = the fixed h-gate tau
    (period knob), applied to BOTH RGs."""
    if template_path is None:
        template_path = SPINAL / "connectome_templates.json"
    tpl = json.loads(Path(template_path).read_text(encoding="utf-8"))
    if "w2laproj" not in tpl:
        raise KeyError(f"{template_path}: no 'w2laproj' key")
    if tau_rg_nap_h is not None:
        NAP["tau_max_h"] = float(tau_rg_nap_h)
    split = make_split_template(tpl["w2laproj"], comm=comm)
    net = W2LSplitNet(template=split, comm=comm, gains=gains or {})
    net.compile(dt=dt)
    return net


if __name__ == "__main__":
    import io
    import sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    for comm in (1.0, 0.0):
        w = build(comm=comm)
        print(f"comm={comm}: neurons {len(w.idx)}  synapses {w.n_synapses}  "
              f"inputs {len(w.inputs)}  outputs {len(w.muscle_outputs)}")
        print("  ports:", w.inputs)
        print("  synapses per tag:", w.synapse_counts)
