"""Sparse conductance computation retaining SNS NumPy reduction order.

Only actual edges are calculated by Numba. Preallocated dense workspaces
retain the reference np.sum axis/order, which is necessary for trajectory
identity in the closed loop. Physics and channel equations are unchanged.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _edges(v, post, pre, gmax, elo, erange, reversal, g, ge):
    for k in range(len(post)):
        p, q = post[k], pre[k]
        value = max(0., min(gmax[k] * ((v[q] - elo[k]) / erange[k]), gmax[k]))
        g[p, q] = value
        ge[p, q] = value * reversal[k]


class ExactFastForward:
    def __init__(self, compiled):
        assert not (compiled.spiking or compiled.electrical or
                    compiled.electrical_rectified), 'Unsupported channel type'
        self.c = compiled
        post, pre = np.nonzero(compiled.g_max_non)
        self.post, self.pre = post, pre
        self.gmax = compiled.g_max_non[post, pre].copy()
        self.elo = compiled.e_lo[post, pre].copy()
        self.erange = (compiled.e_hi - compiled.e_lo)[post, pre].copy()
        self.reversal = compiled.del_e[post, pre].copy()
        self.g = np.zeros_like(compiled.g_max_non)
        self.ge = np.zeros_like(self.g)

    def __getattr__(self, name):
        return getattr(self.c, name)

    def forward(self, x=None):
        c = self.c
        c.V_last = np.copy(c.V)
        i_app = 0 if x is None else np.matmul(c.input_connectivity, x)
        _edges(c.V_last, self.post, self.pre, self.gmax, self.elo,
               self.erange, self.reversal, self.g, self.ge)
        i_syn = np.sum(self.ge, axis=1) - c.V_last * np.sum(self.g, axis=1)
        if c.gated:
            a_inf = 1 / (1 + c.k_a * np.exp(c.slope_a * (c.e_a - c.V_last)))
            b_inf = 1 / (1 + c.k_b * np.exp(c.slope_b * (c.e_b - c.V_last)))
            c_inf = 1 / (1 + c.k_c * np.exp(c.slope_c * (c.e_c - c.V_last)))
            tau_c = c.tau_max_c * c_inf * np.sqrt(
                c.k_c * np.exp(c.slope_c * (c.e_c - c.V_last)))
            c.b_gate_last, c.c_gate_last = np.copy(c.b_gate), np.copy(c.c_gate)
            c.b_gate = c.b_gate_last + c.dt * ((b_inf - c.b_gate_last) / c.tau_max_b)
            c.c_gate = c.c_gate_last + c.dt * ((c_inf - c.c_gate_last) / tau_c)
            i_ion = c.g_ion * (a_inf ** c.pow_a) * (c.b_gate ** c.pow_b) * (
                c.c_gate ** c.pow_c) * (c.e_ion - c.V_last)
            i_gated = np.sum(i_ion, axis=0)
            c.V = c.V_last + c.time_factor_membrane * (
                -c.g_m * (c.V_last - c.V_rest) + c.i_b + i_syn + i_app + i_gated)
        else:
            c.V = c.V_last + c.time_factor_membrane * (
                -c.g_m * (c.V_last - c.V_rest) + c.i_b + i_syn + i_app)
        c.outputs = np.matmul(c.output_voltage_connectivity, c.V)
        return c.outputs
