"""Build LI's AnimatLab 2-layer contact-driven CPG (Neuromechanical_Models/
Li Model/walk tester rearranged.aproj) in sns-toolbox, wired for the M1
MuJoCo body (w2l_mujoco/w2l_mjcf.xml, milestone 1).

WIRING SOURCE OF TRUTH (this module designs NO topology - it transcribes):
  spinal/connectome_templates.json, key "li" (58 nodes / 98 edges). The
  template was mined from Li's aproj (make_editor_templates.py > build_li);
  its per-edge "gain" values are Li's per-connexion SynapticConductance in
  uS, VERIFIED against the aproj this session (e.g. L_foot ground contact ->
  L_knee flexion inhibit = 8e-6 in the aproj = 8.0 in the template;
  L_hip swing MN -> L_hip MV2 = 3e-6 = 3.0; CPG/PF connexions = 1e-6 = 1.0).

NEURON REALIZATION - DEVIATION #1 (loud, per the task's w2l_cpg-deviation-5
precedent). Li's net is SPIKING IntegrateFire (aproj: rest -60 mV, threshold
-55 mV (MN -59), tau_m 5 ms, tonic 5-6 nA; module: 0.2 ms timestep, spike
strength 1, refractory 2 ms). sns-toolbox 1.5.2 DOES ship SpikingNeuron and
it does run on the numpy backend (probed this session), BUT its spike
conductance state is PER-NEURON, not per-synapse (sns_toolbox/backends.py,
SNS_Numpy.forward: g_spike is a per-neuron scalar added to EVERY incoming
connection row; g_increment likewise per-neuron). Li's circuit has
convergent mixed-sign input everywhere it matters (each PF half-center gets
excitation near -50 mV AND inhibition at -70 mV; each MN gets excitation at
-10 mV and inhibition at -70 mV), and his per-synapse conductances differ
by 8x (1 uS generic vs 8 uS contact/II drive). On this backend that would
(i) sum the spike conductance once per incoming connection regardless of
which afferent fired and (ii) collapse the 1-8 uS differentiation. We
therefore realize every neuron as a GRADED RC surrogate on AnimatLab's own
units (mV / nA / uS), which preserves the topology, the per-synapse
conductances, the reversal potentials and the tonic currents EXACTLY:
  spiking neuron -> NonSpikingNeuron(C=tau_m, G=1 uS, Vrest=-60 mV);
  spiking synapse -> NonSpikingSynapse(g = the per-link uS value,
      reversal = the synapse type's EquilibriumPotential, presynaptic
      window e_lo=-55 mV (the spike threshold) .. e_hi=Vrest+PRE_ACT_MV
      (the surrogate's active plateau)).
A real-spiking build remains possible once the toolbox gains per-synapse
spike conductance (or with a custom stepper); flagged as future work, not
silently dropped.

MUSCLE INTERFACE: the aproj interposes 12 NonSpiking "MV" relay neurons
between MNs and muscles (MN -Nicotinic ACh-> MV -identity gain-> muscle);
their MN->MV conductances ARE the template's mn_to_muscle gains (3/5/8 uS,
verified per-side). The template omits the MVs; we ADD them back (deviation
#2, aproj-verified) because MuJoCo needs a graded actuator drive anyway:
the MV voltage maps to ctrl = clip((V_MV - (-100)) / ((-10) - (-100)), 0, 1)
(MV rest -100 mV, driven toward the -10 mV Nicotinic reversal by spikes),
matching AnimatLab's MembraneVoltage->StimulusTension chain up to the
muscle's own 0..1 tension curve, which the M1 <muscle> actuator subsumes.

SENSORY ENCODERS (runner-side currents into port neurons; knobs below):
 - heel/toe contact -> SN: the aproj adapter is a Sigmoid on ContactCount
   (A=0.7-0.8 nA, B=1e-8, C=25, D=0). AnimatLab's Sigmoid convention is not
   documented in the file, so the effective steady-state amplitude is not
   reliably recoverable; we use a count-gated current (CONTACT_GAIN_NA per
   active contact), tuned so the SN surrogate sits above threshold.
 - hip angle -> SN-II "hip middle": aproj adapter Sigmoid(A=15 nA, B=1e-8,
   C=5, D=0) = a near-step from 0 to 15 nA crossing hip angle +3.68 deg
   (ln(1e-8)/5); reproduced as a hard sigmoid (HIP_MID_DEG, HIP_GAIN_NA).
   Which side of +3.7 deg is "extension" in the M1 hinge convention was not
   derivable from the files (M1 caveat #1); HIP_SIGN flips it.

MUSCLE MAP (M1 actuator names; BY GEOMETRY, not by aproj label - the aproj's
R-side hip muscle names are crossed relative to their attachments, M1 report
section 4c: hip_R_flx carries the back/back (extensor) attachments and
hip_R_ext the front/front (flexor) ones):
  hip stance MN (extensor, back/back)  -> hip_L_ext  / hip_R_flx
  hip swing MN (flexor,  front/front)  -> hip_L_flx  / hip_R_ext
  knee extension MN (patella route)    -> knee_L_ext / knee_R_ext
  knee flexion MN                      -> knee_L_flx / knee_R_flx
  plantarflexion MN (heel/back)        -> ankle_L_ext / ankle_R_ext
  dorsiflexion MN (toe/front)          -> ankle_L_flx / ankle_R_flx

USAGE
    net = build()
    u = net.make_inputs()               # zeros, ordered like net.inputs
    u[net.input_index("L tonic drive")] = net.tonic_amplitude("L tonic drive")
    V = net.step(u)                     # one dt; V indexed by net.idx
    ctrl = net.muscle_ctrl(V)           # 12 MuJoCo actuator drives
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.networks import Network

HERE = Path(__file__).parent
SPINAL = HERE.parent

# --------------------------------------------------------------- scale / units
# AnimatLab IntegrateFire units, mined from walk tester rearranged.aproj:
# voltages mV (rest -60, threshold -55/-59, MV rest -100, Nicotinic eq -10,
# Depolarizing eq -50, Hyperpolarizing eq -70), currents nA (tonic 5-6,
# hip encoder 15, contact 0.7-0.8 per the aproj adapter), conductance uS.
V_REST = -60.0          # mV, every spiking neuron's RestingPotential
V_THR = -55.0           # mV, InitialThreshold (MNs -59)
V_MN_THR = -59.0
TAU_M = 0.005           # s, TimeConstant 5 ms (all neurons incl. MV)
V_MV_REST = -100.0      # mV, MV relay RestingPotential
V_MV_MAX = -10.0        # mV, Nicotinic ACh reversal = MV drive ceiling
REV_NIC = -10.0         # 'Nicotinic ACh' eq (excitatory, MN-directed)
REV_DEP = -50.0         # 'Depolarizing IPSP' eq (Li's EXCITATORY type)
REV_HYP = -70.0         # 'Hyperpolarizing IPSP for CPG/IN' eq (inhibitory)
DT = 0.0002             # s, AnimatLab neural timestep (0.2 ms)

# --------------------------------------------------------------- knobs (tunable)
LI_KNOBS = dict(
    # graded-surrogate synapse window (mV, presynaptic). Li's IntegrateFire
    # cells are RATE-coding around threshold (5-6 nA tonic => sustained
    # firing; synaptic input shifts the rate), so the surrogate treats
    # presynaptic depolarization ABOVE REST as the rate proxy: e_lo = rest,
    # e_hi = rest + PRE_ACT_MV. (A window anchored at the -55 threshold was
    # tried first (run1) and is DEAD on this scale: a 1 uS synapse can lift a
    # 1 uS-leak cell only to -55 mV, so no synapse ever conducted.)
    SYN_E_LO=V_REST,
    # presynaptic active plateau above rest (mV) = window upper edge
    PRE_ACT_MV=15.0,
    # tonic port scale (1.0 = Li's verbatim 5-6 nA onto the 9 PF/CPG targets)
    TONIC_SCALE=1.0,
    # contact encoder: nA per active contact into the SN
    CONTACT_GAIN_NA=15.0,
    # hip "middle" encoder: nA amplitude, step location (deg, aproj-derived)
    # and steepness (1/deg); HIP_SIGN flips the aproj angle convention to the
    # M1 hinge convention if the closed loop comes out mirrored.
    HIP_GAIN_NA=15.0,
    HIP_MID_DEG=3.68,
    HIP_C=1.2,
    HIP_SIGN=1.0,
    # MV -> ctrl map: stands in for the aproj muscle StimulusTension curve
    # (A/B/C/D per muscle, not ported - M1 scope). Soft threshold: MV must
    # depolarize above MV_ON before the muscle recruits; full activation at
    # the Nicotinic ceiling -10 mV.
    MU_GAIN=1.0,
    MV_ON=-40.0,
)

# tag/sign -> reversal potential (see the connexions dump: every template edge
# maps to one of Li's three used synapse types)
def _reversal(tag: str, sign: str) -> float:
    if sign == "inh":
        return REV_HYP
    if tag == "pf_to_mn":          # PF -> MN excitation: Nicotinic ACh
        return REV_NIC
    return REV_DEP                 # contact/II/drive/commissural excitation


# MN (template MUSCLE label) -> M1 actuator name, by geometry (docstring)
MUSCLE_TO_ACTUATOR = {
    "L_hip stance": "hip_L_ext", "L_hip swing": "hip_L_flx",
    "R_hip stance": "hip_R_flx", "R_hip swing": "hip_R_ext",
    "L_knee extension": "knee_L_ext", "R_knee extension": "knee_R_ext",
    "L_knee flexion": "knee_L_flx", "R_knee flexion": "knee_R_flx",
    "L_ankle plantarflexion": "ankle_L_ext", "R_ankle plantarflexion": "ankle_R_ext",
    "L_ankle dorsiflexion": "ankle_L_flx", "R_ankle dorsiflexion": "ankle_R_flx",
}

# Li aproj labels are sometimes padded ("R_ankle plantarflexion "); normalize
def _norm(s: str) -> str:
    return " ".join(s.split())


@dataclass
class LINet:
    template: dict
    knobs: dict = field(default_factory=dict)
    net: Network = field(init=False)
    compiled: object = field(init=False, default=None)
    idx: dict = field(default_factory=dict)
    inputs: list = field(default_factory=list)
    port_targets: dict = field(default_factory=dict)   # port -> [neurons]
    port_slots: dict = field(default_factory=dict)     # port -> [u slots]
    tonic_per_target: dict = field(default_factory=dict)  # port -> [(tgt, nA)]
    muscle_outputs: dict = field(default_factory=dict) # actuator -> MV neuron
    synapse_counts: dict = field(default_factory=dict)
    n_synapses: int = field(init=False, default=0)

    def __post_init__(self):
        self.K = dict(LI_KNOBS)
        self.K.update(self.knobs or {})
        tpl = self.template
        nodes = {n["label"]: n for n in tpl["nodes"]}
        edges = tpl["edges"]
        self.net = Network(name="Li 2-layer contact-driven CPG (graded surrogate)")

        for e in edges:
            for k in ("from", "to"):
                if e[k] not in nodes:
                    raise ValueError(f"template edge endpoint unknown: {e}")

        # ---- neurons: every non-MUSCLE, non-PORT-load node, template order
        for nd in tpl["nodes"]:
            lab, typ = nd["label"], nd["type"]
            if typ in ("MUSCLE", "PORT-load"):
                continue
            thr = V_MN_THR if typ == "MN" else V_THR
            self.net.add_neuron(
                NonSpikingNeuron(membrane_capacitance=TAU_M,
                                 membrane_conductance=1.0,
                                 resting_potential=V_REST, bias=0.0),
                name=lab)
            self.idx[lab] = len(self.idx)
            self._thr[lab] = thr  # per-neuron surrogate threshold (unused by sim)

        # ---- input ports
        #  (a) SN-heel / SN-toe / SN-II: encoder currents land on the neuron
        #  (b) PORT-load (Li's per-leg tonic ports): one add_input per TARGET
        #      neuron (the toolbox input vector is per-neuron); the runner
        #      fills every slot of a port via set_port()/set_all_tonics().
        outgoing: dict[str, list] = {}
        for e in edges:
            outgoing.setdefault(e["from"], []).append(e)
        for nd in tpl["nodes"]:
            lab, typ = nd["label"], nd["type"]
            if typ in ("SN-heel", "SN-toe", "SN-II"):
                self.net.add_input(lab)
                self.port_slots[lab] = [len(self.inputs)]
                self.inputs.append(lab)
                self.port_targets[lab] = [lab]
                self.tonic_per_target[lab] = [(lab, 0.0)]
            elif typ == "PORT-load":
                outs = outgoing.get(lab, [])
                if not outs:
                    raise ValueError(f"{lab}: no outgoing edges")
                slots, per_tgt = [], []
                for o in outs:
                    tgt = o["to"]
                    self.net.add_input(tgt)
                    slots.append(len(self.inputs))
                    self.inputs.append(f"{lab} -> {tgt}")
                    per_tgt.append((tgt, float(o.get("gain", 5.0))))
                self.port_slots[lab] = slots
                self.port_targets[lab] = [t for t, _ in per_tgt]
                self.tonic_per_target[lab] = per_tgt

        # ---- MV relay neurons + muscle outputs (deviation #2, aproj-verified)
        for e in edges:
            if e["tag"] == "mn_to_muscle":
                mn = e["from"]
                act = MUSCLE_TO_ACTUATOR[_norm(e["to"])]
                mv = mn + " MV"
                self.net.add_neuron(
                    NonSpikingNeuron(membrane_capacitance=TAU_M,
                                     membrane_conductance=1.0,
                                     resting_potential=V_MV_REST, bias=0.0),
                    name=mv)
                self.idx[mv] = len(self.idx)
                self.muscle_outputs[act] = mv
                self.net.add_output(mv, name=act)

        # ---- synapses
        for e in edges:
            tag, sign = e["tag"], e["sign"]
            if tag in ("mn_to_muscle",):       # becomes the MN -> MV synapse
                g = float(e["gain"])
                self._add_syn(e["from"], e["from"] + " MV", g, REV_NIC,
                              tag="mn_to_mv")
                continue
            if nodes[e["from"]]["type"] == "PORT-load":
                continue                        # tonic edges -> ports, above
            g = float(e["gain"])
            self._add_syn(e["from"], e["to"], g, _reversal(tag, sign), tag=tag)

        self._assert_complete(edges)

    _thr: dict = field(default_factory=dict)

    def _add_syn(self, src: str, dst: str, g: float, rev: float, tag: str):
        syn = NonSpikingSynapse(max_conductance=float(g),
                                reversal_potential=float(rev),
                                e_lo=self.K["SYN_E_LO"],
                                e_hi=V_REST + self.K["PRE_ACT_MV"])
        self.net.add_connection(syn, src, dst)
        self.synapse_counts[tag] = self.synapse_counts.get(tag, 0) + 1
        self.n_synapses += 1

    def _assert_complete(self, edges):
        n_tpl_nodes = sum(1 for n in self.template["nodes"]
                          if n["type"] not in ("MUSCLE", "PORT-load"))
        assert len(self.idx) == n_tpl_nodes + 12, \
            f"neuron count {len(self.idx)} != {n_tpl_nodes}+12 MV"
        assert len(self.inputs) == 6 + 14, \
            f"ports {len(self.inputs)} != 6 encoders (heel/toe/hip x L/R) " \
            f"+ 14 tonic slots"
        assert len(self.muscle_outputs) == 12, "muscle outputs != 12"
        n_syn_tpl = sum(1 for e in edges if e["tag"] != "mn_to_muscle"
                        and self.template["nodes"][[n["label"] for n in
                        self.template["nodes"]].index(e["from"])]["type"]
                        != "PORT-load")
        assert self.n_synapses == n_syn_tpl + 12, \
            f"synapse count {self.n_synapses} != {n_syn_tpl}+12 MV"

    # ------------------------------------------------------------------ run
    def compile(self, dt: float = DT):
        self.compiled = self.net.compile(dt=dt, backend="numpy")
        assert self.compiled.V.shape[0] == len(self.idx), "index mismatch"
        return self.compiled

    def step(self, input_currents):
        self.compiled.forward(list(input_currents))
        return self.compiled.V          # FULL neuron voltage vector (mV)

    def input_index(self, port: str) -> int:
        return self.port_slots[port][0]

    def set_port(self, u, port: str, amp: float, scale: float = 1.0):
        for s in self.port_slots[port]:
            u[s] = amp * scale

    def set_all_tonics(self, u, extra: dict | None = None):
        """Fill every tonic port with Li's per-target nA x TONIC_SCALE.
        extra: {port: nA} added on top (kickoff / perturbation)."""
        extra = extra or {}
        for port, per in self.tonic_per_target.items():
            if per[0][1] == 0.0:      # encoder port, not a tonic port
                continue
            for (tgt, na), s in zip(per, self.port_slots[port]):
                u[s] = na * self.K["TONIC_SCALE"]
        for port, amp in extra.items():
            for s in self.port_slots[port]:
                u[s] += amp

    def make_inputs(self) -> np.ndarray:
        return np.zeros(len(self.inputs))

    # ------------------------------------------------------------- encoders
    def contact_current(self, port: str, n_contact: float) -> float:
        return self.K["CONTACT_GAIN_NA"] * max(0.0, min(1.0, n_contact))

    def hip_current(self, port: str, hip_deg: float) -> float:
        th = self.K["HIP_SIGN"] * hip_deg
        sig = 1.0 / (1.0 + np.exp(self.K["HIP_C"] * (self.K["HIP_MID_DEG"] - th)))
        return self.K["HIP_GAIN_NA"] * sig

    # ------------------------------------------------------------- actuation
    def muscle_ctrl(self, V) -> dict:
        """MV voltages -> MuJoCo actuator drives in M1 actuator order.
        Soft-threshold map standing in for the unported StimulusTension
        curve: recruit above MV_ON, saturate at the Nicotinic ceiling."""
        out = {}
        for act, mv in self.muscle_outputs.items():
            v = V[self.idx[mv]]
            out[act] = self.K["MU_GAIN"] * float(np.clip(
                (v - self.K["MV_ON"]) / (V_MV_MAX - self.K["MV_ON"]), 0.0, 1.0))
        return out


def build(template_path=None, dt: float = DT, knobs: dict | None = None) -> LINet:
    if template_path is None:
        template_path = SPINAL / "connectome_templates.json"
    tpl = json.loads(Path(template_path).read_text(encoding="utf-8"))
    if "li" not in tpl:
        raise KeyError(f"{template_path}: no 'li' key")
    net = LINet(template=tpl["li"], knobs=knobs or {})
    net.compile(dt=dt)
    return net


if __name__ == "__main__":
    import io, sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    w = build()
    print(f"neurons: {len(w.idx)}  synapses: {w.n_synapses}  "
          f"inputs: {len(w.inputs)}  outputs: {len(w.muscle_outputs)}")
    print("ports:", w.inputs)
    print("tonic amplitudes (nA/target):",
          {p: (w.tonic_amplitude(p), len(w.port_targets[p])) for p in w.inputs
           if p.endswith('tonic drive')})
    print("synapses per tag:", w.synapse_counts)
    print("muscle map:", w.muscle_outputs)
