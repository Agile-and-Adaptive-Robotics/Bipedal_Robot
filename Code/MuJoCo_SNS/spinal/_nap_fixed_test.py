"""NaP half-center with FIXED tau_h — v2 (correct override).

v1 failed silently: the installed sns_toolbox (1.5.2) implements
everything inside SNS_Numpy.forward (no __forward_pass__), so the
mangled-name override never ran. v2 subclasses SNS_Numpy and overrides
FORWARD with a verbatim copy, changing exactly one line:
    tau_b = self.tau_max_b          # was: tau_max_b*b_inf*sqrt(...)
(Deng/Animatlab semantics; audit row #2's prescribed patch.)
"""
import io
import sys

import numpy as np
from scipy.signal import find_peaks

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from sns_toolbox import backends as _B
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)
from sns_toolbox.networks import Network


class SNS_NumpyFixedTau(_B.SNS_Numpy):
    """SNS_Numpy with a CONSTANT h-gate time constant (Deng semantics).
    Verbatim forward() copy; the single changed line is marked."""

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
            # *** THE PATCH: constant tau_b (was voltage-dependent) ***
            tau_b = self.tau_max_b
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


E_HI, DT, T_END = 5.0, 0.002, 30.0


def _neu(tau):
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _syn(g, exc):
    return NonSpikingSynapse(max_conductance=float(g),
                             reversal_potential=E_HI if exc else -E_HI,
                             e_lo=0.0, e_hi=E_HI)


def run(g_na, e_m, s_h, e_h, tau_h, e_ion=8.0, g_inh=4.0, g_w=0.4,
        drive=1.0, cm=0.05, s_m=0.8):
    n = Network()
    nap = lambda: NonSpikingNeuronWithPersistentSodiumChannel(
        membrane_capacitance=cm, membrane_conductance=1.0,
        resting_potential=0.0, bias=0.0,
        g_ion=np.array([g_na]), e_ion=np.array([e_ion]),
        k_m=np.array([1.0]), slope_m=np.array([s_m]),
        e_m=np.array([e_m]),
        k_h=np.array([1.0]), slope_h=np.array([s_h]),
        e_h=np.array([e_h]), tau_max_h=np.array([tau_h]))
    n.add_neuron(nap(), name="RG_E")
    n.add_neuron(nap(), name="RG_F")
    n.add_neuron(_neu(0.05), name="InE")
    n.add_neuron(_neu(0.05), name="InF")
    n.add_input("RG_E", name="DRIVE_E")
    n.add_input("RG_F", name="DRIVE_F")
    n.add_connection(_syn(g_inh, True), "RG_E", "InE")
    n.add_connection(_syn(g_inh, False), "InE", "RG_F")
    n.add_connection(_syn(g_inh, True), "RG_F", "InF")
    n.add_connection(_syn(g_inh, False), "InF", "RG_E")
    n.add_connection(_syn(g_w, True), "RG_E", "RG_F")
    n.add_connection(_syn(g_w, True), "RG_F", "RG_E")
    net = n.compile(backend="numpy", dt=DT)
    net.__class__ = SNS_NumpyFixedTau      # swap in the fixed-tau stepper
    n_steps = int(round(T_END / DT))
    kick = int(round(0.2 / DT))
    v = np.zeros(n_steps + 1)
    vfull = np.zeros((n_steps + 1, 2))
    for k in range(n_steps):
        net.forward([drive + (2.0 if k < kick else 0.0), drive])
        V = net.V
        v[k + 1] = V[0] - V[1]
        vfull[k + 1] = (V[0], V[1])
    m = v[int(5.0 / DT):]
    pk, _ = find_peaks(m, prominence=0.3)
    cyc = len(pk)
    per = float(np.median(np.diff(pk))) * DT if cyc >= 2 else float("nan")
    return cyc, per, float(m.max() - m.min()), float(vfull.min()), \
        float(vfull.max())


print("FIXED-tau_h NaP — e_ion sweep (plateau target 3..6 mV)")
print(f"{'g_na':>5s} {'e_ion':>5s} {'e_m':>4s} {'s_h':>5s} {'e_h':>4s} "
      f"{'tau_h':>6s} {'cyc':>4s} {'period':>7s} {'amp':>6s} "
      f"{'Vmin':>6s} {'Vmax':>6s}")
for g_na, e_ion, e_m, s_h, e_h, tau_h in (
        (6.0, 8.0, 1.5, -2.0, 3.0, 0.35),
        (6.0, 8.0, 2.0, -2.0, 3.5, 0.35),
        (6.0, 12.0, 2.0, -2.0, 3.5, 0.35),
        (10.0, 8.0, 2.0, -2.0, 3.5, 0.35),
        (10.0, 12.0, 2.0, -2.0, 3.5, 0.35),
        (10.0, 12.0, 2.0, -2.0, 3.5, 0.70),
        (15.0, 8.0, 2.0, -2.0, 3.5, 0.35),
        (15.0, 12.0, 2.5, -2.0, 3.5, 0.35),
        (15.0, 12.0, 2.5, -2.0, 4.0, 0.70),
        (20.0, 8.0, 2.0, -2.0, 3.5, 0.35),
        (20.0, 12.0, 2.5, -2.0, 4.0, 0.35),
        (20.0, 12.0, 2.5, -2.5, 4.0, 0.70)):
    cyc, per, amp, vmin, vmax = run(g_na, e_m, s_h, e_h, tau_h,
                                    e_ion=e_ion)
    print(f"{g_na:5.2f} {e_ion:5.1f} {e_m:4.1f} {s_h:5.1f} {e_h:4.1f} "
          f"{tau_h:6.2f} {cyc:4d} {per:7.2f} {amp:6.2f} "
          f"{vmin:6.2f} {vmax:6.2f}")
