"""Milestone 3: the 2023 ORIGINAL Walker_2_Layer_CPG neural architecture in
sns-toolbox, wired 1:1 to the M1 MuJoCo body's 12 muscles.

WIRING SOURCE OF TRUTH (this module designs NO topology - it transcribes):
  spinal/connectome_templates.json, key "w2laproj"
  (92 nodes / 163 edges, mined from Neuromechanical_Models\\Walker_2_Layer_CPG\\
  Walker_2_Layer_CPG.aproj). THE 2023 ORIGINAL: a SINGLE LEFT RG (HC-RG-E/F
  persistent-Na half-centers + laminated InE/InF + direct HC<->HC escape
  excitation), endogenous after the Stimulus_1 kickoff (10 nA, 10 ms, -> L RG
  ext). The RIGHT side has NO RG: it is driven ANTI-PHASically through the 4
  CROSSED pf_drive edges (L RG ext -> R Hip/R Knee PF *flx*, L RG flx ->
  R Hip/R Knee PF *ext*), while the 4 ipsilateral pf_drive edges drive the L
  PFs in phase. Ankle MNs ride the Knee PFs (pf_to_mn, the biarticular-gas
  routing mined from the aproj). Renshaw pool = 11 RCs with the shared
  knee-extensor RE (L Knee MN ext RE is excited by BOTH sides' knee-ext MNs
  and inhibits both). Ia (SN-Ia -> IN-IaIN -> antagonist MN, mutual) and Ib
  (SN-Ib -> homonymous MN exc) chains are built and asserted but receive ZERO
  input at this milestone (the Deng afferent doctrine: nothing sensory
  reaches RG/PF; in AIR the afferents are silent altogether).

CONventions identical to ..\\w2l_cpg\\build_w2l_net.py (which transcribed the
LATER bilateralrg variant): 0..5 mV scale, E_REV_EXC +8 / E_REV_INH -5,
NaP half-centers stepped by SNS_NumpyFixedTau (stock tau_h(V) quenches them;
AnimatLab treats tau_h.max as fixed - same semantics), TAU table and
W2L_GAINS calibration table verbatim from that module. The original template
has NO contact/c1/V3/aff_HC tags, so those gain keys are simply unused here.

MUSCLE MAP (template MUSCLE node -> M1 actuator; by anatomy, not by name -
the aproj's R-side hip muscle NAMES are crossed vs their attachment sites,
M1 report section 4c, kinematically confirmed 2026-09-25 by tendon-length
finite differences on w2l_mjcf_fixed.xml: flexion = -qpos on EVERY joint):

    L Hip flx -> hip_L_flx      L Hip ext -> hip_L_ext
    R Hip flx -> hip_R_ext      R Hip ext -> hip_R_flx   (crossed)
    L/R knee flexion -> knee_*_flx    L/R knee extension -> knee_*_ext
    L/R ankle plantarflexion -> ankle_*_ext (PF)
    L/R ankle dorsiflexion  -> ankle_*_flx (DF)

MN -> MuJoCo ctrl uses the stack activation convention a = clip(V/E_HI,0,1).

USAGE
    net = build()                        # compile, dt = 2 ms default
    u = net.make_inputs()                # zeros, len = 3
    u[net.input_index("Stimulus_1")] = 10.0    # nA, t in [0, 0.01]
    u[net.input_index("TONIC L RG ext")] = 2.0 # nA, the drive regime
    V = net.step(u)                      # V[net.idx["L Knee MN ext"]] etc.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)
from sns_toolbox.networks import Network

HERE = Path(__file__).parent
SPINAL = HERE.parent
if str(SPINAL) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(SPINAL))
from build_network import SNS_NumpyFixedTau  # noqa: E402  (read-only import)
from w2l_cpg.build_w2l_net import DT, NAP, TAU, W2L_GAINS, _neu, _syn  # noqa: E402

E_HI = 5.0   # mV, = w2l_cpg SYN_E_HI / stack activation scale

# template MUSCLE label -> M1 actuator name (see module docstring)
MUSCLE_MAP = {
    "L Hip ext": "hip_L_ext", "L Hip flx": "hip_L_flx",
    "R Hip ext": "hip_R_flx", "R Hip flx": "hip_R_ext",   # aproj name crossing
    "L_knee extension": "knee_L_ext", "L_knee flexion": "knee_L_flx",
    "R_knee extension": "knee_R_ext", "R_knee flexion": "knee_R_flx",
    "L_ankle plantarflexion": "ankle_L_ext",
    "L_ankle dorsiflexion": "ankle_L_flx",
    "R_ankle plantarflexion": "ankle_R_ext",
    "R_ankle dorsiflexion": "ankle_R_flx",
}


@dataclass
class W2LOrigNet:
    """Compiled wrapper over the transcribed 2023-original W2L network."""
    template: dict
    gains: dict = field(default_factory=dict)   # merged over W2L_GAINS
    net: Network = field(init=False)
    compiled: object = field(init=False, default=None)
    idx: dict = field(default_factory=dict)
    inputs: list = field(default_factory=list)   # ordered port names
    port_targets: dict = field(default_factory=dict)
    muscle_outputs: dict = field(default_factory=dict)  # actuator -> MN label
    synapse_counts: dict = field(default_factory=dict)
    n_synapses: int = field(init=False, default=0)

    def __post_init__(self):
        self.G = dict(W2L_GAINS)
        self.G.update(self.gains or {})
        tpl = self.template
        nodes = {n["label"]: n for n in tpl["nodes"]}
        edges = tpl["edges"]
        self.net = Network(name="W2L 2023 original (single LH RG)")

        for e in edges:
            for k in ("from", "to"):
                if e[k] not in nodes:
                    raise ValueError(f"template edge endpoint unknown: {e}")

        # ---- neurons (template order)
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

        # ---- input ports: the template PORT-load stimulus + TONIC drive on
        # the two RG half-centers (the smoke's measured drive regime: F tonic
        # is the escape engine and must exceed E tonic, w2l_cpg README DEv #2).
        # Input ports add no synapses; zero current leaves the graph unchanged.
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
        for lab in ("L RG ext", "L RG flx"):
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
            # the two direct HC<->HC edges (HC -> HC, EXC) are the escape
            # assist at the template's raw 0.5; the four laminated edges
            # (E->InE, F->InF exc; InE->F, InF->E inh) carry the stack's
            # proven rg_mutual_inh. Same split as w2l_cpg's transcription.
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
        raise ValueError(f"unknown template tag: {tag}")

    def _assert_complete(self, edges, n_types):
        expected = {}
        for e in edges:
            if e["tag"] == "mn_to_muscle":
                continue
            if n_types[e["from"]] == "PORT-load":
                continue
            expected[e["tag"]] = expected.get(e["tag"], 0) + 1
        assert self.n_synapses == sum(expected.values()) == 150, \
            f"synapse transcription mismatch: {self.n_synapses} vs 150"
        assert dict(self.synapse_counts) == expected, "per-tag mismatch"
        assert len(self.idx) == 79, f"neuron count {len(self.idx)} != 79"
        assert len(self.inputs) == 3, f"port count {len(self.inputs)} != 3"
        assert len(self.muscle_outputs) == 12, "muscle outputs != 12"

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
        """MN membrane potential [mV] -> MuJoCo ctrl (stack convention
        a = clip(V/E_HI, 0, 1)), keyed by actuator name."""
        return {a: float(min(max(V[self.idx[mn]] / E_HI, 0.0), 1.0))
                for a, mn in self.muscle_outputs.items()}


# template node type -> membrane tau key (subset of w2l_cpg's table)
_TYPE_TAU = {
    "HC-RG-E": "rg", "HC-RG-F": "rg",
    "IN-InE": "rg", "IN-InF": "rg",
    "HC-PF-E": "pf", "HC-PF-F": "pf", "IN-PF": "pf",
    "IN-IaIN": "afferent", "SN-Ia": "afferent", "SN-Ib": "afferent",
    "MN": "mn", "RC": "rc",
}


def build(template_path=None, dt: float = DT, gains: dict | None = None,
          tau_rg_nap_h: float | None = None) -> W2LOrigNet:
    """Load the w2laproj template and build + compile the network.

    gains: per-key overrides merged over W2L_GAINS.
    tau_rg_nap_h: override the fixed h-gate tau = the PERIOD knob
    (w2l_cpg TAU['rg_nap_h'], default 0.25 s ~ 1 s intrinsic period)."""
    if template_path is None:
        template_path = SPINAL / "connectome_templates.json"
    tpl = json.loads(Path(template_path).read_text(encoding="utf-8"))
    if "w2laproj" not in tpl:
        raise KeyError(f"{template_path}: no 'w2laproj' key")
    if tau_rg_nap_h is not None:
        NAP["tau_max_h"] = float(tau_rg_nap_h)
    net = W2LOrigNet(template=tpl["w2laproj"], gains=gains or {})
    net.compile(dt=dt)
    return net


if __name__ == "__main__":
    import io
    import sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    w = build()
    print(f"neurons: {len(w.idx)}  synapses: {w.n_synapses}  "
          f"inputs: {len(w.inputs)}  outputs: {len(w.muscle_outputs)}")
    print("ports:", w.inputs)
    print("ports inject into:", {p: w.port_targets[p] for p in w.inputs})
    print("synapses per tag:", w.synapse_counts)
    print("muscle map:", w.muscle_outputs)
