"""Shevtsova 2026 lumbar circuit - SNS-TOOLBOX build (Ben directive
2026-10-20: "Repeat her results in MuJoCo is one of your many /goals").

Paper: Shevtsova, Lockhart, Rybak, Magnuson & Danner, eLife 14:RP107480
(DOI 10.7554/eLife.107480). LUMBAR-ONLY biped reduction transcribed from
replication/shevtsova_rules.json (Table 1 weights, every edge cited
there). The calibrated parameter set is the one that makes the
standalone numpy ODE PASS (shevtsova_standalone.py + verify_shevtsova.py
- read those first): division-form Boltzmann gates, tau 150/400 ms bell,
k_tau 6.67 mV, drives dE 0.1 (RG-F) / bE 0.1 (RG-E) / dI 0.75 (V0D) /
0.15 (V0V), dt 1 ms.

This module re-expresses THAT circuit with sns_toolbox primitives
(Network / NonSpikingNeuron / NonSpikingNeuronWithPersistentSodiumChannel
/ NonSpikingSynapse + the repo's SNS_NumpyFixedTau backend patch imported
from ../build_network.py), so the circuit definition is MuJoCo-stack
facing: input ports carry alpha, stepping is the production
compiled.forward(port_currents) pattern.

EXACTNESS OF THE MAPPING. After per-neuron normalization (divide each
population's equation by its own gL: g_m = 1), the toolbox update
    V += dt/c_m * ( -(V-Vrest) + sum_j g_j f(V_pre_j)(E_rev_j - V)
                    + i_app + g_ion m_inf h (E_Na - V) )
with c_m = C/gL, g_j = gSyn*w_j/gL_post, e_lo/e_hi = Vthr/Vmax is the
SAME Euler map as the paper's eq 1-13 term by term (verified symbol by
symbol against shevtsova_standalone.simulate). Her Boltzmann gates ARE
expressible exactly in the toolbox by slope-sign choice:
    m_inf: k=1, e_a=-40 mV, slope_a=+1/6 per-mV   (paper V12,m -40, k_m -6)
    h_inf: k=1, e_b=-45 mV, slope_b=-0.25 per-mV  (paper V12,h -45, k_h +4)

DEVIATION FLAGS (every departure from the standalone's dynamics):
  S-1  tau_h BELL -> the repo's PROVEN fixed-tau patch was tried FIRST
      (build_network.py SNS_NumpyFixedTau, Deng/Animatlab semantics;
      stock voltage-dependent tau_h measured 10x-too-fast, 2026-09-16),
      sweep over her whole stated 150-400 ms band (tau_mode="fixed"):
        tau 0.15/0.20/0.25/0.30/0.35/0.40 s at alpha=0.55 ->
        7.41/6.04/5.11/4.45/3.95/3.55 Hz (0.30 s matches the standalone's
        4.51 Hz; 0.25-0.40 hold her 3-6 Hz trot band) - but at alpha=1.05
        EVERY fixed value LATCHES at the co-active fixed point (RG-F/RG-E
        both pinned ~-42.3 mV, f(V_RG-F) ~ 0.15 > 0.1, no cycles). Her
      bell tau_h(V) = tau0 + (taumax-tau0)/cosh((V-V12)/k_tau) is slow
      through the burst band (-45..-25 mV: 256-400 ms) and fast in
      silence (-70..-55 mV: 152-175 ms); that asymmetry is what
      destabilizes the co-active fixed point at high drive. DEFAULT
      therefore restores HER EXACT eq-8 bell through the same
      backend-patch mechanism (class SNS_NumpyBellTau below: single
      tau_b line of the repo's patched forward overridden with her
      published tau_h constants - her stated form, no free knobs);
      tau_mode="fixed" keeps the repo-pattern reduction for the record
      (holds C1, latches C2 - documented above).
  S-2  h-gate INIT CONVENTION. sns_toolbox 1.5.2's compiler initializes
      b_gate with a sign-flipped Boltzmann vs its own forward()
      (compilers.py uses exp(slope*(V-e)), forward uses exp(slope*(e-V)));
      compile_circuit() overwrites b_gate/b_gate_last/b_gate_0 with the
      forward-convention h_inf(V_init) after compiling.
  S-3  DRIVES AS PORT CURRENTS. The paper's tonic drive is a CONDUCTANCE
      gSyn*(d*alpha+b) toward ESynE/ESynI. The toolbox input port is a
      current source only, so each driven population gets its own input
      port and the stepping loop supplies the EXACT paper drive current
      evaluated at the backend's own V (i_app = -(gSyn*(dE*alpha+bE)*
      (V-ESynE) + gSyn*(dI*alpha+bI)*(V-ESynI))/gL) - the production
      runner pattern of computing port currents from state. Numerically
      exact at Euler order; the only architectural mapping.
  S-4  UPDATE ORDER. The toolbox advances the h gate BEFORE the membrane
      step, so i_ion uses the NEW h; the standalone computes iNap with
      the OLD h. Half-step phase difference: dt/tau_h <= 1/150 = 0.7%.
  S-5  No membrane clipping (standalone clips V to [-100,50] mV) and no
      OU noise (sigma 0 = the verification protocol, standalone F-E).

Census: 22 populations (4 persistent-Na RG centers + 18 non-spiking
interneurons), 32 synapses (Table 1, cited per edge in the rules json),
8 input ports (RG-F_l/r, RG-E_l/r, V0D_l/r, V0V_l/r).

Run:      D:/Anaconda/envs/myo/python.exe shevtsova_sns.py
Verify:   D:/Anaconda/envs/myo/python.exe verify_shevtsova_sns.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SPINAL = HERE.parent
RULES = HERE / "shevtsova_rules.json"

sys.path.insert(0, str(SPINAL))       # SNS_NumpyFixedTau (repo-proven)
sys.path.insert(0, str(HERE))         # shevtsova_standalone (analysis)

from sns_toolbox.connections import NonSpikingSynapse          # noqa: E402
from sns_toolbox.neurons import NonSpikingNeuron               # noqa: E402
from sns_toolbox.neurons import \
    NonSpikingNeuronWithPersistentSodiumChannel                # noqa: E402
from sns_toolbox.networks import Network                       # noqa: E402

from build_network import SNS_NumpyFixedTau                    # noqa: E402


class SNS_NumpyBellTau(SNS_NumpyFixedTau):
    """SNS_NumpyFixedTau with Shevtsova 2026's OWN eq-8 tau_h bell
    restored on the single tau_b line (flag S-1): the fixed-tau repo
    reduction holds her trot band but latches at alpha=1.05 across her
    whole stated 150-400 ms band, so the default build patches her
    published tau_h(V) = tau0 + (taumax-tau0)/cosh((V-V12)/k_tau)
    (constants TAU_BELL, no free knobs) through the same one-line
    mechanism the repo uses for its fixed-tau patch. Verbatim forward()
    copy of SNS_NumpyFixedTau (itself the verbatim sns_toolbox 1.5.2
    copy); the single changed line is marked. Her bell applies to every
    neuron's h gate; non-RG neurons carry g_ion = 0 so theirs is inert.
    """

    def forward(self, x=None):
        self.V_last = np.copy(self.V)
        if x is None:
            i_app = 0
        else:
            i_app = np.matmul(self.input_connectivity, x)
        g_syn = np.maximum(0, np.minimum(
            self.g_max_non * ((self.V_last - self.e_lo) /
                              (self.e_hi - self.e_lo)), self.g_max_non))
        if self.spiking:
            self.theta_last = np.copy(self.theta)
            self.g_spike = self.g_spike * (1 - self.time_factor_synapse)
            g_syn += self.g_spike
        i_syn = np.sum(g_syn * self.del_e, axis=1) - \
            self.V_last * np.sum(g_syn, axis=1)
        if self.electrical:
            i_syn += (np.sum(self.g_electrical * self.V_last, axis=1) -
                      self.V_last * np.sum(self.g_electrical, axis=1))
        if self.electrical_rectified:
            mask = np.subtract.outer(self.V_last, self.V_last).T > 0
            masked_g = mask * self.g_rectified
            diag_masked = masked_g + masked_g.T - \
                np.diag(masked_g.diagonal())
            i_syn += np.sum(diag_masked * self.V_last, axis=1) - \
                self.V_last * np.sum(diag_masked, axis=1)
        if self.gated:
            a_inf = 1 / (1 + self.k_a * np.exp(
                self.slope_a * (self.e_a - self.V_last)))
            b_inf = 1 / (1 + self.k_b * np.exp(
                self.slope_b * (self.e_b - self.V_last)))
            c_inf = 1 / (1 + self.k_c * np.exp(
                self.slope_c * (self.e_c - self.V_last)))
            # *** THE PATCH: her eq-8 bell tau_h (was constant tau_max_b) ***
            tau_b = self.tau0_h + (self.taumax_h - self.tau0_h) / np.cosh(
                (self.V_last - self.e_tau_h) / self.k_tau_h)
            tau_c = self.tau_max_c * c_inf * np.sqrt(
                self.k_c * np.exp(self.slope_c * (self.e_c - self.V_last)))
            self.b_gate_last = np.copy(self.b_gate)
            self.c_gate_last = np.copy(self.c_gate)
            self.b_gate = self.b_gate_last + self.dt * (
                (b_inf - self.b_gate_last) / tau_b)
            self.c_gate = self.c_gate_last + self.dt * (
                (c_inf - self.c_gate_last) / tau_c)
            i_ion = self.g_ion * (a_inf ** self.pow_a) * \
                (self.b_gate ** self.pow_b) * \
                (self.c_gate ** self.pow_c) * (self.e_ion - self.V_last)
            i_gated = np.sum(i_ion, axis=0)
            self.V = self.V_last + self.time_factor_membrane * (
                -self.g_m * (self.V_last - self.V_rest) + self.i_b +
                i_syn + i_app + i_gated)
        else:
            self.V = self.V_last + self.time_factor_membrane * (
                -self.g_m * (self.V_last - self.V_rest) + self.i_b +
                i_syn + i_app)
        if self.spiking:
            self.theta = self.theta_last + self.time_factor_threshold * (
                self.theta_leak * (self.theta_0 - self.theta_last) +
                self.m * (self.V_last - self.V_rest))
            self.spikes = np.sign(np.minimum(0, self.theta - self.V))
            if self.delay:
                self.spike_buffer = np.roll(self.spike_buffer, 1, axis=0)
                self.spike_buffer[0, :] = self.spikes
                self.delayed_spikes[self.spike_rows, self.spike_cols] = \
                    self.spike_buffer[self.buffer_steps, self.buffer_nrns]
                self.g_spike += np.minimum(
                    (-self.delayed_spikes * self.g_increment),
                    (-self.delayed_spikes) *
                    (self.g_max_spike - self.g_spike))
            else:
                self.g_spike += np.minimum(
                    (-self.spikes * self.g_increment),
                    (-self.spikes) * (self.g_max_spike - self.g_spike))
            self.V = ((self.V - self.V_reset) * (self.spikes + 1)) + \
                self.V_reset
            self.theta = np.maximum(
                self.theta_increment,
                self.theta_floor - self.theta) * (-self.spikes) + self.theta
        self.outputs = np.matmul(self.output_voltage_connectivity, self.V)
        if self.spiking:
            self.outputs += np.matmul(self.output_spike_connectivity,
                                      -self.spikes)
        return self.outputs

# ------------------------------------------------- paper constants (pF/nS/mV)
# Identical to shevtsova_standalone.py (the calibrated set); see its
# flags F-A..F-E for provenance.
C_PF, GSYN = 10.0, 10.0            # membrane capacitance; gSynE = gSynI
E_SYN_E, E_SYN_I = -10.0, -75.0    # synaptic reversals
E_NA, G_NAP = 50.0, 4.5            # NaP reversal + conductance (RG only)
V_THR, V_MAX = -50.0, 0.0          # synaptic threshold/saturation (f(V))
GL_RG, EL_RG = 4.5, -62.5          # leak, RG centers (standalone F-B)
GL_IN, EL_IN = 2.8, -60.0          # leak, all INs
V_INIT = -60.0                     # mV, standalone initial condition
# Boltzmann gates, division form -> toolbox slopes (exact; see docstring)
SLOPE_M, E_M = 1.0 / 6.0, -40.0    # m_inf: paper k_m -6 mV, V12,m -40
SLOPE_H, E_H = -0.25, -45.0        # h_inf: paper k_h +4 mV, V12,h -45
# tau_h (flag S-1): the paper's bell (Methods eq 8), constants = the
# standalone's calibrated set (its flag F-A); "fixed" = the repo's
# SNS_NumpyFixedTau reduction (retune band = her stated tau_h range).
TAU_BELL = dict(tau0_s=0.150, taumax_s=0.400, v12_mV=-35.0, k_mV=6.67)
TAU_H_DEFAULT_S = 0.300              # fixed-mode stand-in (band-matched)
TAU_H_BAND_S = (0.150, 0.400)
DRIVE = dict(  # (dE, bE, dI, bI) per population-substring; standalone F-D
    RG=("RG-F", 0.1, 0.0, 0.0, 0.0), RGE=("RG-E", 0.0, 0.1, 0.0, 0.0),
    V0D=("V0D", 0.0, 0.0, 0.75, 0.0), V0V=("V0V", 0.0, 0.0, 0.15, 0.0),
)


def _drives(label: str):
    for _, (sub, dE, bE, dI, bI) in DRIVE.items():
        if sub in label:
            return dE, bE, dI, bI
    return 0.0, 0.0, 0.0, 0.0


def _gl_of(label: str) -> float:
    return GL_RG if ("RG-F" in label or "RG-E" in label) else GL_IN


def build_circuit(tau_h_s: float = TAU_H_DEFAULT_S,
                  tau_mode: str = "bell"):
    """Build the 22-population lumbar circuit in sns_toolbox.

    tau_mode "bell" (default): her eq-8 tau_h bell in SNS_NumpyBellTau
    (flag S-1). tau_mode "fixed": the repo's SNS_NumpyFixedTau reduction
    with constant tau_h_s (band 150-400 ms enforced).

    Returns a dict with the Network plus the bookkeeping the simulate
    loop and shevtsova_standalone.analyze need (labels/idx in build
    order = compiled.V order).
    """
    if tau_mode not in ("bell", "fixed"):
        raise ValueError(f"tau_mode must be 'bell' or 'fixed', got {tau_mode}")
    if tau_mode == "fixed" and not (
            TAU_H_BAND_S[0] <= tau_h_s <= TAU_H_BAND_S[1]):
        raise ValueError(
            f"tau_h_s={tau_h_s} outside the paper's stated tau_h range "
            f"{TAU_H_BAND_S} s (flag S-1: retune only inside her band)")
    rules = json.loads(RULES.read_text(encoding="utf-8"))
    labels = [n["label"] for n in rules["nodes"] if n["type"] != "DRIVE"]
    idx = {l: i for i, l in enumerate(labels)}
    net = Network(name="shevtsova2026 lumbar sns")
    ports = []                       # (label, dE, bE, dI, bI, gL), port order
    for l in labels:
        gl = _gl_of(l)
        if "RG-F" in l or "RG-E" in l:
            neu = NonSpikingNeuronWithPersistentSodiumChannel(
                membrane_capacitance=C_PF / gl / 1000.0,  # C/gL in seconds
                membrane_conductance=1.0,
                resting_potential=EL_RG, bias=0.0,
                g_ion=np.array([G_NAP / gl]),   # normalized NaP conductance
                e_ion=np.array([E_NA]),
                k_m=np.array([1.0]), slope_m=np.array([SLOPE_M]),
                e_m=np.array([E_M]),
                k_h=np.array([1.0]), slope_h=np.array([SLOPE_H]),
                e_h=np.array([E_H]),
                tau_max_h=np.array([tau_h_s]))  # flag S-1
        else:
            neu = NonSpikingNeuron(
                membrane_capacitance=C_PF / gl / 1000.0,
                membrane_conductance=1.0,
                resting_potential=EL_IN, bias=0.0)
        net.add_neuron(neu, name=l, initial_value=V_INIT)
        dE, bE, dI, bI = _drives(l)
        if (dE, bE, dI, bI) != (0.0, 0.0, 0.0, 0.0):
            net.add_input(l)                       # port named after its cell
            ports.append((l, dE, bE, dI, bI, gl))
    for s in rules["synapses"]:
        w = float(s["gain"])
        exc = s["sign"] == "exc"
        gl_dst = _gl_of(s["to"])
        syn = NonSpikingSynapse(
            max_conductance=GSYN * w / gl_dst,    # gSyn*w normalized by gL
            reversal_potential=E_SYN_E if exc else E_SYN_I,
            e_lo=V_THR, e_hi=V_MAX)               # paper f(V) exactly
        net.add_connection(syn, s["from"], s["to"])
    return dict(labels=labels, idx=idx, net=net, ports=ports,
                n_pop=len(labels),
                n_syn=len(rules["synapses"]),
                tau_h_s=tau_h_s, tau_mode=tau_mode)


def compile_circuit(c: dict, dt_ms: float = 1.0):
    """Compile with the tau_h backend patch + h-init fix (flags S-1/S-2)."""
    m = c["net"].compile(dt=dt_ms / 1000.0, backend="numpy")
    if c["tau_mode"] == "bell":
        m.__class__ = SNS_NumpyBellTau
        m.tau0_h = TAU_BELL["tau0_s"]
        m.taumax_h = TAU_BELL["taumax_s"]
        m.e_tau_h = TAU_BELL["v12_mV"]
        m.k_tau_h = TAU_BELL["k_mV"]
    else:
        m.__class__ = SNS_NumpyFixedTau
    h_inf_init = 1.0 / (1.0 + np.exp(SLOPE_H * (E_H - V_INIT)))
    m.b_gate[:] = h_inf_init
    m.b_gate_last[:] = h_inf_init
    m.b_gate_0[:] = h_inf_init
    return m


def simulate(alpha: float, t_end_s: float = 20.0, dt_ms: float = 1.0,
             tau_h_s: float = TAU_H_DEFAULT_S, tau_mode: str = "bell",
             circuit: dict = None, model=None):
    """Step the compiled toolbox circuit at the paper's dt = 1 ms.

    Fresh initial condition per call (V = -60 mV, h = h_inf(-60)),
    mirroring shevtsova_standalone.simulate. alpha enters through the
    drive input ports (flag S-3).
    """
    c = circuit or build_circuit(tau_h_s, tau_mode)
    m = model or compile_circuit(c, dt_ms)
    n = c["n_pop"]
    nsteps = int(round(t_end_s * 1000.0 / dt_ms))
    out = np.zeros((nsteps + 1, n))
    out[0] = m.V
    ports = c["ports"]
    x = np.zeros(len(ports))
    for k in range(nsteps):
        V = m.V                        # becomes V_last inside forward()
        for j, (lab, dE, bE, dI, bI, gl) in enumerate(ports):
            vi = V[c["idx"][lab]]
            x[j] = -(GSYN * (dE * alpha + bE) * (vi - E_SYN_E)
                     + GSYN * (dI * alpha + bI) * (vi - E_SYN_I)) / gl
        m.forward(list(x))
        out[k + 1] = m.V
    t = np.arange(nsteps + 1) * dt_ms / 1000.0
    return t, out, c


if __name__ == "__main__":
    import io
    import sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    from shevtsova_standalone import analyze
    c = build_circuit()
    print(f"census: {c['n_pop']} populations, {c['n_syn']} synapses, "
          f"{len(c['ports'])} input ports, tau_mode {c['tau_mode']} "
          f"(fixed-tau stand-in {c['tau_h_s']:.3f} s; flag S-1)")
    for a in (0.55, 1.05):
        t, out, _ = simulate(a, t_end_s=20.0, circuit=c)
        r = analyze(t, out, c)
        print(f"alpha={a}: freq L {r['l']['freq']:.2f} R {r['r']['freq']:.2f} Hz, "
              f"LR phase {r['phase_med']:.3f} (n={r['phase_n']}), "
              f"ext_frac L {r['l']['ext_frac']:.2f}")
