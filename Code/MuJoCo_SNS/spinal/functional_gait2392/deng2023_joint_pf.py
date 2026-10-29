"""Experimental human adaptation of Deng 2023 Fig. 1E, not a replication.

Three single premotor integrators per leg combine the joint PFs corresponding
to HAM/BFP-like, rectus, and gastrocnemius actions. Rat input-conductance ratios
are normalized to a configurable human-network budget, not copied as microS.
The W2L joint PFs remain leaky cells, unlike the paper's NaP joint oscillators.
Other biarticular families retain the prior experimental implementation.
HAM/RECT/GAS retain base-W2L reflex paths, without the prototype's additional
direct Ia/II-to-biarticular-PF inputs. Use a dedicated process; production
selectors are untouched.
"""
import json
from dataclasses import dataclass, field
import network as F
import build_network as BN
import build_network_w2lvar as W
from params import TAU, G

# Fig. 1E and author-manuscript Table 2. These are relative source weights.
JOINT_INPUTS = {
    'HAM': (('PF_HIP-E_{side}', .967), ('PF_KNEE-F_{side}', .365)),
    'RECT': (('PF_HIP-F_{side}', .236), ('PF_KNEE-E_{side}', .872)),
    'GAS': (('PF_KNEE-F_{side}', 1.), ('PF_FUNC_ANKLE_Z_P_{side}', .3)),
}

@dataclass
class JointIntegrationNetwork(F.FunctionalNetwork):
    bi_integrators: dict = field(default_factory=dict)
    joint_integration_edges: list = field(default_factory=list)

    def _build_pf_w2l(self, n, side):
        W.W2LVarnet._build_pf_w2l(self, n, side)
        for label in ('HIP_X', 'HIP_Y', 'FOOT_X'):
            self._pair(n, label, side)
        if side == self.sides[0]:
            for label in ('LUM_X', 'LUM_Y', 'LUM_Z'):
                self._pair(n, label, 'mid')
        self._pair(n, 'ANKLE_Z', side, rg=True)
        present = {info['family'] for act, info in self.routing['actuator_mapping'].items()
                   if act.endswith('_' + side)}
        for family in F.BI_FAMILIES:
            if family not in present:
                continue
            if family not in JOINT_INPUTS:
                self._pair(n, 'BI_' + family, side, rg=True)
                continue
            cell = 'PF_BI_JOINT_' + family + '_' + side
            self._add(cell, TAU['pf'], n)
            self.bi_integrators[family + '_' + side] = cell
            sources = JOINT_INPUTS[family]
            total = sum(weight for _, weight in sources)
            for template, weight in sources:
                source = template.format(side=side)
                gain = self.options['joint_mix_gain'] * W.GAINS['w2l_rg_to_pf'] * weight / total
                n.add_connection(BN._syn(gain, True), source, cell)
                self.joint_integration_edges.append(dict(source=source, target=cell,
                    raw_rat_weight=weight, normalized_fraction=weight / total,
                    human_network_conductance=gain))
        names = [p['name'] for p in n.populations]
        for connection in n.connections:
            if (names[connection['source']] == 'TOEDF_' + side and
                    names[connection['destination']] == 'PF_KNEE-F_' + side):
                connection['destination'] = names.index(self.extra_pairs['ANKLE_Z_' + side][1])

    def _wire_muscle_w2l(self, n, act, mi):
        family = self.routing['actuator_mapping'][act]['family']
        if family not in JOINT_INPUTS:
            return F.FunctionalNetwork._wire_muscle_w2l(self, n, act, mi)
        W.W2LVarnet._wire_muscle_w2l(self, n, act, mi)
        names = [p['name'] for p in n.populations]
        mn = 'MN_' + act
        old = [c for c in n.connections if names[c['destination']] == mn and
               names[c['source']].startswith('PF_') and
               not names[c['source']].startswith('PF_IN')]
        budget = sum(float(c['params']['max_conductance']) for c in old)
        self.legacy_pf_budget[act] = budget
        for connection in old:
            n.connections.remove(connection)
        cell = self.bi_integrators[family + '_' + mi.side]
        self._edge(n, cell, mn, self.options['bi_gain'] * budget)
        # Preserve the prior model-aware axis recruitment as a separate input.
        leg = 0 if mi.side == 'r' else 7
        for key, row in [('HIP_X_' + mi.side, leg + 1),
                         ('HIP_Y_' + mi.side, leg + 2),
                         ('FOOT_X_' + mi.side, leg + 5)]:
            for sign, source in zip((1., -1.), self.extra_pairs[key]):
                self._edge(n, source, mn, 1e-12, record=(act, row, sign, key))

def build(acts, dt=.002, interleg=True, options=None):
    options = {**F.DEFAULTS, 'joint_mix_gain': 1., **(options or {})}
    if not options['axes'] or not options['biarticular'] or options['primitive']:
        raise ValueError('Joint-integration variant requires axes/biarticular on, primitive off.')
    if options['joint_mix_gain'] < 0:
        raise ValueError('joint_mix_gain must be nonnegative.')
    for key, value in W.VARIANT_G.items():
        G[key] = float(value)
    muscles = {act: W.classify(act) for act in acts}
    if not all(muscles.values()):
        raise ValueError('Unclassified actuator.')
    routing = json.loads((F.OUT / 'routing_audit.json').read_text())
    net = JointIntegrationNetwork(muscles=muscles,
        sides=tuple(sorted({m.side for m in muscles.values()})), interleg=interleg,
        options=options, routing=routing)
    net.compile(dt)
    return net
