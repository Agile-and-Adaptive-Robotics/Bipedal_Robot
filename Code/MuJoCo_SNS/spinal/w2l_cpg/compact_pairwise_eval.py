"""Opt-in fixed-tau graded SNS workspace, with reference-order synaptic sums.

Keeps the existing SparsePairwiseForward reduction tree and mutable edge gains.
Reuses state/current buffers, evaluates voltage-dependent gate functions only
where their slopes are nonzero, and computes ionic force only on nonzero ionic
channels. Inactive gate states still advance. No fastmath, timestep, topology,
or neuron-law change. Returned state/output arrays are borrowed buffers.

Call refresh_parameters() after changing intrinsic parameters. Construction is
lazy so a subsequently installed FSA profile is consumed on the first step.
Spiking/electrical channels remain explicitly unsupported, as in the parent.
"""
import numpy as np
from numba import njit
from sparse_pairwise_eval import SparsePairwiseForward


@njit(cache=True)
def compact_sums(v, post, pre, gmax, elo, erange, reversal, group, lane,
                 tail, group_post, group_leaf, counts, joins, root,
                 active_posts, acc, values, out):
    """Same eight-lane leaves/tree; materialize only nonempty row/leaf groups."""
    acc[:] = 0.
    values[:] = 0.
    out[:] = 0.
    for e in range(len(post)):
        f = (v[pre[e]] - elo[e]) / erange[e]
        if f < 0.:
            f = 0.
        elif f > 1.:
            f = 1.
        g = gmax[e] * f
        ge = g * reversal[e]
        b, k = group[e], lane[e]
        if tail[e]:
            acc[b, 8 + k, 0] = g
            acc[b, 8 + k, 1] = ge
        else:
            acc[b, k, 0] += g
            acc[b, k, 1] += ge
    for b in range(len(group_post)):
        p, leaf = group_post[b], group_leaf[b]
        for kind in range(2):
            r = acc[b, :, kind]
            value = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))
            for k in range(counts[leaf] % 8):
                value += r[8 + k]
            values[p, leaf, kind] = value
    offset = len(counts)
    for p in active_posts:
        for j in range(len(joins)):
            for kind in range(2):
                values[p, offset + j, kind] = values[p, joins[j, 0], kind] + values[p, joins[j, 1], kind]
        out[p, 0] = values[p, root, 0]
        out[p, 1] = values[p, root, 1]


class CompactPairwiseForward(SparsePairwiseForward):
    def __init__(self, compiled):
        super().__init__(compiled)
        # These parent dense workspaces are unused by sparse_sums.
        del self.g, self.ge
        count = len(self.counts)
        keys, self.edge_group = np.unique(self.post * count + self.leaf, return_inverse=True)
        self.group_post, self.group_leaf = keys // count, keys % count
        self.active_posts = np.unique(self.post)
        self.acc = np.zeros((len(keys), 16, 2))
        c = self.c
        self.next_voltage = np.empty_like(c.V)
        self.chemical = np.empty_like(c.V)
        self.current = np.empty_like(c.V)
        self.ionic = np.empty_like(c.V)
        self.negative_leak = np.empty_like(c.g_m)
        self.ready = False

    def refresh_parameters(self):
        """Rebuild constant/index tables without resetting neural states."""
        c = self.c
        self.gate_tables = {}
        if c.gated:
            for label in ('a', 'b', 'c'):
                k = getattr(c, 'k_' + label)
                slope = getattr(c, 'slope_' + label)
                e = getattr(c, 'e_' + label)
                indices = np.nonzero(slope)
                constant = 1. / (1. + k)
                self.gate_tables[label] = {
                    'indices': indices, 'k': k[indices].copy(),
                    'slope': slope[indices].copy(), 'e': e[indices].copy(),
                    'inf': constant.copy(), 'work': np.empty(len(indices[0])),
                    'voltage': np.empty(len(indices[0])),
                }
            shape = c.g_ion.shape
            self.b_last = np.empty_like(c.b_gate)
            self.c_last = np.empty_like(c.c_gate)
            self.b_next = np.empty_like(c.b_gate)
            self.c_next = np.empty_like(c.c_gate)
            self.gate_work = np.empty(shape)
            self.tau_c = c.tau_max_c * self.gate_tables['c']['inf'] * np.sqrt(c.k_c)
            self.ion_indices = np.nonzero(c.g_ion)
            count = len(self.ion_indices[0])
            self.ion_work = np.empty(count)
            self.power_work = np.empty(count)
            self.ion_voltage = np.empty(count)
            self.ion_matrix = np.zeros(shape)
        self.ready = True

    def _gate(self, label):
        table = self.gate_tables[label]
        ids = table['indices']
        if len(ids[0]):
            np.take(self.c.V_last, ids[1], out=table['voltage'])
            work = table['work']
            np.subtract(table['e'], table['voltage'], out=work)
            np.multiply(table['slope'], work, out=work)
            np.exp(work, out=work)
            np.multiply(table['k'], work, out=work)
            if label == 'c':
                np.sqrt(work, out=table['voltage'])
                self.tau_c[ids] = table['voltage']
            np.add(1., work, out=work)
            np.reciprocal(work, out=work)
            table['inf'][ids] = work
            if label == 'c':
                self.tau_c[ids] *= self.c.tau_max_c[ids] * work
        return table['inf']

    def forward(self, x=None):
        c = self.c
        if not self.ready:
            self.refresh_parameters()
        if not np.isfinite(c.V).all():
            raise FloatingPointError('Nonfinite SNS state outside compact kernel contract')
        np.copyto(self.last_voltage, c.V)
        c.V_last = self.last_voltage
        if x is None:
            self.i_app.fill(0.)
        elif self.input_onehot:
            np.take(x, self.iq, out=self.input_values)
            self.i_app[self.ip] = self.input_values
        else:
            np.matmul(c.input_connectivity, x, out=self.i_app)
        compact_sums(c.V_last, self.post, self.pre, self.gmax, self.elo,
                     self.erange, self.reversal, self.edge_group, self.lane, self.tail,
                     self.group_post, self.group_leaf, self.counts, self.joins, self.root,
                     self.active_posts, self.acc, self.values, self.sums)
        np.multiply(c.V_last, self.sums[:, 0], out=self.chemical)
        np.subtract(self.sums[:, 1], self.chemical, out=self.chemical)
        if c.gated:
            a_inf, b_inf, c_inf = self._gate('a'), self._gate('b'), self._gate('c')
            np.copyto(self.b_last, c.b_gate)
            np.copyto(self.c_last, c.c_gate)
            c.b_gate_last, c.c_gate_last = self.b_last, self.c_last
            np.subtract(b_inf, self.b_last, out=self.gate_work)
            np.divide(self.gate_work, c.tau_max_b, out=self.gate_work)
            np.multiply(c.dt, self.gate_work, out=self.gate_work)
            np.add(self.b_last, self.gate_work, out=self.b_next)
            np.subtract(c_inf, self.c_last, out=self.gate_work)
            np.divide(self.gate_work, self.tau_c, out=self.gate_work)
            np.multiply(c.dt, self.gate_work, out=self.gate_work)
            np.add(self.c_last, self.gate_work, out=self.c_next)
            c.b_gate, c.c_gate = self.b_next, self.c_next
            ids = self.ion_indices
            np.power(a_inf[ids], c.pow_a[ids], out=self.ion_work)
            np.multiply(c.g_ion[ids], self.ion_work, out=self.ion_work)
            np.power(c.b_gate[ids], c.pow_b[ids], out=self.power_work)
            np.multiply(self.ion_work, self.power_work, out=self.ion_work)
            np.power(c.c_gate[ids], c.pow_c[ids], out=self.power_work)
            np.multiply(self.ion_work, self.power_work, out=self.ion_work)
            np.take(c.V_last, ids[1], out=self.ion_voltage)
            np.subtract(c.e_ion[ids], self.ion_voltage, out=self.power_work)
            np.multiply(self.ion_work, self.power_work, out=self.ion_work)
            self.ion_matrix[ids] = self.ion_work
            np.sum(self.ion_matrix, axis=0, out=self.ionic)
        np.subtract(c.V_last, c.V_rest, out=self.current)
        np.negative(c.g_m, out=self.negative_leak)
        np.multiply(self.negative_leak, self.current, out=self.current)
        np.add(self.current, c.i_b, out=self.current)
        np.add(self.current, self.chemical, out=self.current)
        np.add(self.current, self.i_app, out=self.current)
        if c.gated:
            np.add(self.current, self.ionic, out=self.current)
        np.multiply(c.time_factor_membrane, self.current, out=self.current)
        np.add(c.V_last, self.current, out=self.next_voltage)
        c.V = self.next_voltage
        if self.output_onehot:
            np.take(c.V, self.oq, out=self.output_values)
        else:
            np.matmul(c.output_voltage_connectivity, c.V, out=self.output_values)
        c.outputs = self.output_values
        return c.outputs
