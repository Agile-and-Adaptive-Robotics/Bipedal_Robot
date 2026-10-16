"""SPIKING twin of the split-RG W2L network (2026-10-03 overnight
campaign, Ben's ask: the BPA walker driven by a SPIKING controller).

Same cell-class mapping as the verified gait2392 spiking mirror
(spinal/build_network_spiking.py, goal-1 2026-10-02: topology mirror
exact, basin 12/12) applied to the split-RG W2L template
(build_w2l_split_net.make_split_template):

  HC-RG-E/F   -> adapting LIF ('rg'; threshold-increment SFA, tau_theta
                 = the tau knob = the analog tau_max_h role)
  HC-PF-E/F   -> adapting LIF ('pf', mild SFA)
  IN-* (IaIN, InE, InF, PF, C, V3) -> plain LIF ('in')
  RC          -> LIF ('rc' thresholds)
  SN-Ia/SN-Ib -> LIF ('aff' thresholds)
  MN          -> NonSpikingNeuron, UNCHANGED 0..5 mV frame; spiking
                 sources reach MNs via calibrated spike synapses
                 (_ssyn onto_analog=True); the graded MN->RC edge is the
                 mirror's rate-encoder pattern (_ns2s).

Topology is TRANSFORMED from the same split template (identical edge
multiset; only cell types and synapse models change). Readout taps
RD_<RG> (calibrated non-spiking sinks) repoint idx so the driver's
V[net.idx['L RG ext']] logging reads an analog level.

API identical to build_w2l_split_net.W2LSplitNet: idx / ridx / inputs /
make_inputs / input_index / step / muscle_ctrl / muscle_outputs.

Usage:
    import build_w2l_split_net_spiking as BS
    net = BS.build(comm=1.0, tau=0.25, thr_inc=1.0)
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SPINAL = HERE.parent
ROOT = SPINAL.parent
for p in (str(HERE), str(SPINAL), str(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from sns_toolbox.connections import SpikingSynapse  # noqa: E402
from sns_toolbox.neurons import NonSpikingNeuron  # noqa: E402

import params as P  # noqa: E402  (TAU["rg_nap_h"] read by _spk('rg'))
import build_network_spiking as MS  # noqa: E402  (SPIKE/CAL/_spk/_ssyn)
from sns_toolbox.networks import Network  # noqa: E402
import build_w2l_split_net as B  # noqa: E402 (template + gains)
import build_w2l_orig_net as ORIG  # noqa: E402 (MUSCLE_MAP)

E_HI = ORIG.E_HI
MUSCLE_MAP = ORIG.MUSCLE_MAP
_TYPE_TAU = B._TYPE_TAU
TONIC_PORTS = B.TONIC_PORTS

# template node type -> spiking cell kind (None = stays analog)
TYPE_KIND = {
    "HC-RG-E": "rg", "HC-RG-F": "rg",
    "HC-PF-E": "pf", "HC-PF-F": "pf",
    "IN-IaIN": "in", "IN-InE": "in", "IN-InF": "in", "IN-PF": "in",
    "IN-C": "in", "IN-V3": "in",
    "RC": "rc",
    "SN-Ia": "aff", "SN-Ib": "aff",
    "MN": None,
}


@dataclass
class W2LSplitNetSpiking:
    """Spiking twin over the split template (W2LSplitNet interface)."""
    template: dict
    comm: float = 1.0
    gains: dict = field(default_factory=dict)
    tau_rg: float = 0.25          # seconds, RG SFA tau (the period knob)
    thr_inc: float = 4.0          # mV/spike RG adaptation (CAL winner)
    net: object = field(init=False, default=None)
    compiled: object = field(init=False, default=None)
    idx: dict = field(default_factory=dict)
    ridx: dict = field(default_factory=dict)
    inputs: list = field(default_factory=list)
    port_targets: dict = field(default_factory=dict)
    muscle_outputs: dict = field(default_factory=dict)
    synapse_counts: dict = field(default_factory=dict)
    n_synapses: int = field(init=False, default=0)
    n_sub: int = 4                # 2 ms call / 0.5 ms net sub-step
    spike_counts: object = field(init=False, default=None)

    def __post_init__(self):
        self.G = dict(B.W2L_GAINS)
        self.G.update(self.gains or {})
        # the big mirror's _spk('rg') reads params.TAU at CALL time
        P.TAU["rg_nap_h"] = float(self.tau_rg)
        MS.CAL["rg_thr_inc"] = float(self.thr_inc)
        tpl = self.template
        nodes = {n["label"]: n for n in tpl["nodes"]}
        edges = tpl["edges"]
        self.net = Network(name="W2L split-RG (M4, SPIKING twin)")
        for e in edges:
            for k in ("from", "to"):
                if e[k] not in nodes:
                    raise ValueError(f"template edge endpoint unknown: {e}")

        # ---- neurons (cell mapping table; MN analog)
        for nd in tpl["nodes"]:
            lab, typ = nd["label"], nd["type"]
            if typ in ("MUSCLE", "PORT-load"):
                continue
            kind = TYPE_KIND[typ]
            tau = B.TAU[_TYPE_TAU[typ]]
            if kind is None:                      # MN: analog 0..5 mV
                self.net.add_neuron(MS._neu(tau), name=lab)
            else:
                self.net.add_neuron(MS._spk(tau, kind), name=lab)
            self.idx[lab] = len(self.idx)

        # ---- input ports. The analog twin injects port current
        # DIRECTLY into its cells; here every port first charges an
        # ANALOG relay neuron (the big mirror's DRIVE pattern: analog
        # brainstem surrogate -> calibrated graded synapse _ns2s onto
        # the spiking RG, the k_ns rate map). Direct nA injection into
        # the real-mV LIFs leaves them subthreshold (bias 16 + te 3 =
        # V_eq -51 mV < -50 thr) - measured in the first smoke.
        outgoing = {}
        for e in edges:
            outgoing.setdefault(e["from"], []).append(e)
        DESC = {"L RG ext": 1.7, "R RG ext": 1.7,      # descend_to_rg_e
                "L RG flx": 1.4, "R RG flx": 1.4}      # descend_to_rg_f
        for nd in tpl["nodes"]:
            lab, typ = nd["label"], nd["type"]
            if typ == "PORT-load":
                outs = outgoing.get(lab, [])
                if len(outs) != 1:
                    raise ValueError(f"{lab}: expected 1 outgoing edge")
                tgt = outs[0]["to"]
                self._add_port(lab, tgt, 1.7)
        for lab in TONIC_PORTS:
            self._add_port("TONIC " + lab, lab, DESC[lab])
        self._wire_remaining(tpl, nodes, edges)

    def _add_port(self, port, target, g_desc):
        """port current -> analog relay -> graded _ns2s onto target."""
        relay = f"RELAY_{port.replace(' ', '_')}"
        self.net.add_neuron(MS._neu(0.05), name=relay)
        self.idx[relay] = len(self.idx)
        self.net.add_input(relay)
        self.inputs.append(port)
        self.port_targets[port] = target
        self.net.add_connection(MS._ns2s(g_desc, exc=True), relay, target)

    def _wire_remaining(self, tpl, nodes, edges):
        """muscle outputs + chemical synapses + readout taps."""
        # ---- muscle outputs (MN pass-through)
        for e in edges:
            if e["tag"] == "mn_to_muscle":
                self.net.add_output(e["from"], name=e["to"])
                self.muscle_outputs[MUSCLE_MAP[e["to"]]] = e["from"]

        # ---- synapses: per-tag conductance (same _gain as the analog
        # twin), synapse MODEL chosen by source/target cell class
        n_types = {lab: nd["type"] for lab, nd in nodes.items()}
        analog = {lab: TYPE_KIND.get(n_types[lab]) is None
                  for lab in nodes}
        for e in edges:
            tag, sign = e["tag"], e["sign"]
            if tag == "mn_to_muscle":
                continue
            if n_types[e["from"]] == "PORT-load":
                continue
            g = B.W2LSplitNet._gain(self, tag, sign, e, n_types)
            exc = sign == "exc"
            if analog[e["to"]]:            # spiking -> MN (calibrated)
                syn = MS._ssyn(g, exc=exc, onto_analog=True)
            elif analog[e["from"]]:        # MN -> RC rate encoder
                syn = MS._ns2s(g, exc=exc)
            else:                          # spiking -> spiking
                syn = MS._ssyn(g, exc=exc)
            self.net.add_connection(syn, e["from"], e["to"])
            self.synapse_counts[tag] = self.synapse_counts.get(tag, 0) + 1
            self.n_synapses += 1

        # ---- readout taps for the 4 RG cells (driver logging)
        for lab in ("L RG ext", "L RG flx", "R RG ext", "R RG flx"):
            rd = f"RD_{lab.replace(' ', '_')}"
            self.net.add_neuron(MS._neu(0.05), name=rd)
            g_inc = float(MS.CAL["g_inc_readout"])
            self.net.add_connection(
                SpikingSynapse(max_conductance=4.0 * g_inc,
                               reversal_potential=MS.E_REV_EXC,
                               time_constant=MS.TAU_SYN["readout"],
                               transmission_delay=0,
                               conductance_increment=g_inc), lab, rd)
            self.ridx[lab] = self.idx[lab]
            self.idx[rd] = len(self.idx)       # the RD neuron itself
            self.idx[lab] = self.idx[rd]       # repoint at the readout

    # ------------------------------------------------------------------ run
    def compile(self, dt: float = 0.0005):
        self.compiled = self.net.compile(dt=dt, backend="numpy")
        _ns_mask = np.array(
            [not self.net.populations[i]["type"].params.get("spiking",
                                                            False)
             for i in range(len(self.net.populations))])
        for attr in ("theta_0", "theta", "theta_last"):
            arr = getattr(self.compiled, attr, None)
            if arr is not None:
                arr[_ns_mask] = 1.0e6
        self.spike_counts = np.zeros(len(self.idx))
        return self.compiled

    def step(self, input_currents: np.ndarray) -> np.ndarray:
        """One 2 ms DRIVER step = n_sub network sub-steps (constant u).
        RG entries of the returned V are the READOUT analog levels."""
        u = list(input_currents)
        for _ in range(self.n_sub):
            self.compiled.forward(u)
            self.spike_counts += (self.compiled.spikes == -1)
        return self.compiled.V

    def input_index(self, port: str) -> int:
        return self.inputs.index(port)

    def make_inputs(self) -> np.ndarray:
        return np.zeros(len(self.inputs))

    def muscle_ctrl(self, V: np.ndarray) -> dict:
        return {a: float(min(max(V[self.idx[mn]] / E_HI, 0.0), 1.0))
                for a, mn in self.muscle_outputs.items()}


def build(template_path=None, dt: float = 0.0005, gains: dict | None = None,
          tau: float = 0.25, thr_inc: float = 4.0,
          comm: float = 1.0) -> W2LSplitNetSpiking:
    """Build + compile the spiking split-RG twin. Same knobs as the
    analog twin plus tau (RG SFA tau, seconds) and thr_inc (mV/spike)."""
    if template_path is None:
        template_path = SPINAL / "connectome_templates.json"
    tpl = json.loads(Path(template_path).read_text(encoding="utf-8"))
    if "w2laproj" not in tpl:
        raise KeyError(f"{template_path}: no 'w2laproj' key")
    split = B.make_split_template(tpl["w2laproj"], comm=comm)
    net = W2LSplitNetSpiking(template=split, comm=comm,
                             gains=gains or {}, tau_rg=float(tau),
                             thr_inc=float(thr_inc))
    net.compile(dt=dt)
    return net


if __name__ == "__main__":
    import io as _io
    sys.stdout = _io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                   errors="replace")
    for comm in (1.0, 0.0):
        w = build(comm=comm)
        print(f"comm={comm}: neurons {len(w.idx)} (+{len(w.ridx)} raw) "
              f"synapses {w.n_synapses}+4 inputs {len(w.inputs)} "
              f"outputs {len(w.muscle_outputs)}")
        # 1 s smoke: tonic drive both sides, count spikes
        u = w.make_inputs()
        iE, iF = w.input_index("TONIC L RG ext"), w.input_index(
            "TONIC L RG flx")
        iEr, iFr = w.input_index("TONIC R RG ext"), w.input_index(
            "TONIC R RG flx")
        iS = w.input_index("Stimulus_1")
        for i in range(500):                      # 1 s at 2 ms calls
            u[iE], u[iEr] = 3.0, 3.0
            u[iF], u[iFr] = 4.0, 4.0
            u[iS] = 10.0 if i == 0 else 0.0
            V = w.step(u)
        sc = w.spike_counts
        raw = {lab: int(sc[w.ridx[lab]]) for lab in
               ("L RG ext", "L RG flx", "R RG ext", "R RG flx")}
        rd = {lab: round(float(V[w.idx[lab]]), 2) for lab in
              ("L RG ext", "R RG ext")}
        print(f"  1 s @ te=3/tf=4: RG spikes {raw}, readouts {rd} mV")
