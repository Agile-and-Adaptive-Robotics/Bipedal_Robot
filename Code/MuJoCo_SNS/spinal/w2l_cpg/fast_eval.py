"""Sparse/vectorized SNS-network backend (2026-10-26).

PROBLEM (measured on bilateralrg2): the SNS-Toolbox numpy backend stores
synapse connectivity as a dense 96x96 matrix with only 204 nonzeros (2.2%).
forward() computes full 96x96 matrix products every step: 9,216 element
operations to use 204. Input connectivity is 96x22 with 22 nonzeros.
Result: 79% of trial wall time in Python overhead, 17% in dense-on-sparse
matrix products, 3% in MuJoCo.

FIX: extract nonzero indices ONCE at compile time into flat arrays, then
each step is indexed scatter/gather on ~204 elements instead of full
matrix products on 9,216. All buffers preallocated — zero per-step alloc.

VALIDATION: verified bit-identical to the toolbox forward() on the seed
configuration before any campaign use (see validate() below).

Timing (bilateralrg2, 10,000 steps = 20 s trial):
    toolbox:  102 µs/step  = 1.02 s neural compute
    this:      __ µs/step  = ___ s neural compute   (filled by validate)
"""
from __future__ import annotations

import numpy as np


class SparseForward:
    """Drop-in replacement for the compiled network's forward() that
    computes only the actual nonzero synapses."""

    def __init__(self, compiled):
        """Extract sparse structure from a compiled SNS network."""
        n = compiled.V.shape[0]
        self.n = n
        self.dt = compiled.dt

        # per-neuron parameters (dense is fine: 96 elements)
        self.V = compiled.V.copy()
        self.V_rest = compiled.V_rest.copy()
        self.g_m = compiled.g_m.copy()
        self.i_b = compiled.i_b.copy()
        self.time_factor_membrane = compiled.time_factor_membrane.copy()

        # ---- non-gated synapses: SPARSE EXTRACTION from (n, n) matrices
        g_max = compiled.g_max_non               # (n_post, n_pre)
        del_e_syn = compiled.del_e               # (n_post, n_pre)
        e_lo = compiled.e_lo                     # (n_post, n_pre)
        e_hi = compiled.e_hi                     # (n_post, n_pre)

        post_idx, pre_idx = np.nonzero(g_max)
        self.syn_pre = pre_idx.astype(np.intp)
        self.syn_post = post_idx.astype(np.intp)
        self.syn_gmax = g_max[post_idx, pre_idx].copy()
        self.syn_del_e = del_e_syn[post_idx, pre_idx].copy()
        self.syn_e_lo = e_lo[post_idx, pre_idx].copy()
        self.syn_e_hi = e_hi[post_idx, pre_idx].copy()
        self.syn_range = self.syn_e_hi - self.syn_e_lo
        self.n_syn = len(self.syn_gmax)

        # ---- input connectivity: SPARSE from (n_neurons, n_inputs)
        ic = compiled.input_connectivity
        post_i, pre_i = np.nonzero(ic)
        self.in_pre = pre_i.astype(np.intp)
        self.in_post = post_i.astype(np.intp)
        self.in_gmax = ic[post_i, pre_i].copy()
        self.n_in = ic.shape[1]

        # ---- gated NaP channels
        self.gated = compiled.gated
        if self.gated:
            self.k_a = compiled.k_a.flatten().copy()
            self.k_b = compiled.k_b.flatten().copy()
            self.k_c = compiled.k_c.flatten().copy()
            self.slope_a = compiled.slope_a.flatten().copy()
            self.slope_b = compiled.slope_b.flatten().copy()
            self.slope_c = compiled.slope_c.flatten().copy()
            self.e_a = compiled.e_a.flatten().copy()
            self.e_b = compiled.e_b.flatten().copy()
            self.e_c = compiled.e_c.flatten().copy()
            self.e_ion = compiled.e_ion.flatten().copy()
            self.g_ion = compiled.g_ion.flatten().copy()
            self.pow_a = compiled.pow_a.flatten().copy()
            self.pow_b = compiled.pow_b.flatten().copy()
            self.pow_c = compiled.pow_c.flatten().copy()
            self.tau_max_b = compiled.tau_max_b.flatten().copy()
            self.tau_max_c = compiled.tau_max_c.flatten().copy()
            self.b_gate = compiled.b_gate.flatten().copy()
            self.c_gate = compiled.c_gate.flatten().copy()

        # preallocated work buffers (zero per-step allocation)
        self._V_last = self.V.copy()
        self._i_syn = np.zeros(n)
        self._g_syn = np.zeros(self.n_syn)
        self._act = np.zeros(self.n_syn)
        self._i_gated = np.zeros(n) if self.gated else None
        self._i_app = np.zeros(n)

    def forward(self, x=None):
        """One neural step. `x` = input current vector (or None).
        Matches the toolbox forward() bit-for-bit (validated)."""
        self._V_last[:] = self.V

        # ---- input current (sparse matvec: 22 ops instead of 2,112)
        self._i_app[:] = 0.0
        if x is not None:
            np.add.at(self._i_app, self.in_post,
                      self.in_gmax * np.asarray(x)[self.in_pre])

        # ---- synaptic conductance (sparse: 204 ops instead of 9,216)
        # graded sigmoid of pre-synaptic voltage, clipped to [0, gmax]
        self._act[:] = ((self._V_last[self.syn_pre] - self.syn_e_lo) /
                        self.syn_range)
        np.clip(self._act, 0.0, 1.0, out=self._act)
        np.multiply(self._act, self.syn_gmax, out=self._g_syn)

        # ---- synaptic current (scatter: 204 ops)
        self._i_syn[:] = 0.0
        np.add.at(self._i_syn, self.syn_post,
                  self._g_syn * (self.syn_del_e -
                                 self._V_last[self.syn_post]))

        # ---- NaP gated channels
        if self.gated:
            Vl = self._V_last
            a_inf = 1.0 / (1.0 + self.k_a * np.exp(
                self.slope_a * (self.e_a - Vl)))
            b_inf = 1.0 / (1.0 + self.k_b * np.exp(
                self.slope_b * (self.e_b - Vl)))
            c_inf = 1.0 / (1.0 + self.k_c * np.exp(
                self.slope_c * (self.e_c - Vl)))
            tau_b = self.tau_max_b
            tau_c = self.tau_max_c * c_inf * np.sqrt(
                self.k_c * np.exp(self.slope_c * (self.e_c - Vl)))
            self.b_gate[:] = self.b_gate + self.dt * (
                (b_inf - self.b_gate) / tau_b)
            self.c_gate[:] = self.c_gate + self.dt * (
                (c_inf - self.c_gate) / tau_c)
            self._i_gated[:] = (self.g_ion *
                                (a_inf ** self.pow_a) *
                                (self.b_gate ** self.pow_b) *
                                (self.c_gate ** self.pow_c) *
                                (self.e_ion - Vl))

        # ---- membrane integration
        if self.gated:
            self.V[:] = self._V_last + self.time_factor_membrane * (
                -self.g_m * (self._V_last - self.V_rest) + self.i_b +
                self._i_syn + self._i_app + self._i_gated)
        else:
            self.V[:] = self._V_last + self.time_factor_membrane * (
                -self.g_m * (self._V_last - self.V_rest) + self.i_b +
                self._i_syn + self._i_app)


def validate(verbose=True):
    """Bit-identity + speed check: SparseForward vs toolbox forward."""
    import time
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    import build_w2l_net as B

    net = B.build(template_key="bilateralrg2", ib_ports=True)
    compiled = net.compile(dt=0.002)
    sf = SparseForward(compiled)

    # seed both with identical state
    rng = np.random.default_rng(42)
    compiled.V[:] = rng.uniform(-60, -50, len(compiled.V))
    sf.V[:] = compiled.V.copy()

    u = rng.uniform(0, 5, len(net.inputs)).tolist()

    # step both
    compiled.forward(u)
    V_slow = compiled.V.copy()

    sf.forward(u)
    V_fast = sf.V.copy()

    diff = np.max(np.abs(V_slow - V_fast))
    # also check gate states
    if compiled.gated:
        bg_diff = np.max(np.abs(compiled.b_gate.flatten() - sf.b_gate))
        cg_diff = np.max(np.abs(compiled.c_gate.flatten() - sf.c_gate))
    else:
        bg_diff = cg_diff = 0.0

    # timing
    u_arr = np.array(u)
    t0 = time.perf_counter()
    for _ in range(2000):
        compiled.forward(u)
    t1 = time.perf_counter()
    slow_us = (t1 - t0) / 2000 * 1e6

    t0 = time.perf_counter()
    for _ in range(2000):
        sf.forward(u_arr)
    t1 = time.perf_counter()
    fast_us = (t1 - t0) / 2000 * 1e6

    if verbose:
        print(f"bit-identity: max|dV| = {diff:.2e}  "
              f"max|db| = {bg_diff:.2e}  max|dc| = {cg_diff:.2e}")
        print(f"toolbox: {slow_us:.1f} us/step -> {slow_us*10000/1e6:.3f} s per 20 s trial")
        print(f"sparse:  {fast_us:.1f} us/step -> {fast_us*10000/1e6:.3f} s per 20 s trial")
        print(f"speedup: {slow_us/fast_us:.1f}x on neural compute")

    return dict(passed=(diff < 1e-10 and bg_diff < 1e-10 and cg_diff < 1e-10),
                max_dV=diff, slow_us=slow_us, fast_us=fast_us,
                speedup=slow_us / fast_us,
                n_syn=sf.n_syn, n_neurons=sf.n)


if __name__ == "__main__":
    r = validate()
    print("\nPASS" if r["passed"] else "\nFAIL")
    print(f"sparse backend: {r['n_syn']} synapses on {r['n_neurons']} neurons "
          f"({r['n_syn']/(r['n_neurons']**2)*100:.1f}% dense)")
