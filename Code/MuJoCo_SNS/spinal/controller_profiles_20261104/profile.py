"""Compose a declared circuit revision before compilation; never run a trial.

Archived builders/checkpoints are unchanged. This is a graded human adaptation,
not a complete Rybak replication. Gains below are initial functional budgets for
retraining, not fitted results or literal cat-model conductances.
"""
from collections import defaultdict
from copy import deepcopy
import json
from pathlib import Path
import sys

import numpy as np

S = Path(__file__).resolve().parent.parent
for p in (S, S / 'hunt_training', S / 'functional_gait2392', S / 'w2l_cpg'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

FAMILIES = ('hunt', 's3k', 'w2lvar', 'syn6')
DISPLAY_NAMES = dict(hunt='BiW2Lvar+', s3k='s3k', w2lvar='w2lvar', syn6='syn6')
DEFAULTS = dict(
    schema='gait_four_commands_central_afferents_v1',
    command_voltage_fullscale_mV=5., command_current_fullscale_nA=5., command_tau_s=.05,
    flexor_reference_current_nA=1.2, extensor_reference_current_nA=1.2,
    contralateral_flexor_fraction=1. / 3.,
    alpha_presynaptic_k=2.,
    rg_type_budget_nA=.5, pf_type_budget_nA=.5,
    heel_reference_conductance_uS=.1,
    heel_on=.05, heel_off=.02, heel_rearm_s=.08, heel_decay_s=.10,
    command=dict(F_l=0., E_l=0., F_r=0., E_r=0.),
    pattern_enable=1.,
)


def nodes(net):
    return [p['name'] for p in net.net.populations]


def connections(net):
    names = nodes(net)
    for c in net.net.connections:
        yield c, names[c['source']], names[c['destination']]


def central(name):
    return name.startswith(('RG_', 'InE_', 'InF_', 'PF_'))


def pf_role(net, name):
    """Phase class for contact auditing; a mixed primitive is explicit."""
    if name.startswith('PF_S'):
        side = name[-1]
        k = int(name.split('_')[1][1:]) - 1
        return ('PF', side, 'phase', net.syn_families[side][k])
    if name.startswith(('PF_IN_E_', 'PF_E1_', 'PF_E2_')):
        return ('PF', name[-1], 'phase', 'E')
    if name.startswith(('PF_IN_F_', 'PF_F1_', 'PF_F2_')):
        return ('PF', name[-1], 'phase', 'F')
    for joint in ('HIP', 'KNEE', 'ANK'):
        for phase in ('E', 'F'):
            if name.startswith((f'PF_{joint}-{phase}_', f'PF_IN_{joint}-{phase}_')):
                return ('PF', name[-1], joint, phase)
    for key, pair in getattr(net, 'extra_pairs', {}).items():
        if name in pair:
            if key.startswith('ANKLE_Z_'):
                return ('PF', key[-1], 'ANKLE', 'E' if name == pair[0] else 'F')
            return ('PF', key[-1], key, 'P' if name == pair[0] else 'N')
    return None


def contact_flags(net):
    """Trace contacts to first central entry, including noncentral relays.

    Do not traverse intrinsic RG/PF circuitry: excitation of InE followed by
    inhibition of F is the intended complementary action, not a second entry.
    """
    adjacency = defaultdict(list)
    roots = []
    for c, a, b in connections(net):
        g = float(c['params']['max_conductance'])
        if g > 0:
            adjacency[a].append((b, 1 if c['params']['reversal_potential'] > 0 else -1))
    roots = [n for n in nodes(net) if n.startswith(('HEEL_', 'TOE_', 'LBIN_'))]
    flags, entries = [], []
    for root in roots:
        queue = [(root, 1, [root])]
        for current, sign, path in queue:
            if len(path) > 5:
                continue
            for target, edge_sign in adjacency[current]:
                if target in path:
                    continue
                total = sign * edge_sign
                # Inspect the sensory entry itself, before intrinsic reciprocal
                # inhibition. Ben's baseline allows no inhibitory contact entry,
                # and no toe-specific entry at RG level.
                if total < 0:
                    flags.append(dict(source=root, target=target, sign=total,
                                      layer='PF' if target.startswith('PF_') else 'RG' if central(target) else 'local',
                                      path=path + [target],
                                      issue='inhibitory_contact_entry'))
                if root.startswith(('TOE_', 'TOEDF_')) and target.startswith(('RG_', 'InE_', 'InF_')):
                    flags.append(dict(source=root, target=target, sign=total,
                                      path=path + [target],
                                      issue='toe_contact_at_RG'))
                role = None
                if target.startswith(('RG_E_', 'InE_')):
                    role = ('RG', target[-1], 'phase', 'E')
                elif target.startswith(('RG_F_', 'InF_')):
                    role = ('RG', target[-1], 'phase', 'F')
                elif target.startswith('PF_'):
                    role = pf_role(net, target)
                if central(target):
                    if role:
                        entries.append(dict(source=root, target=target, sign=total,
                                            role=role, path=path + [target]))
                    continue
                if target.startswith(('TOEDF_', 'HEEL_', 'TOE_', 'LBIN_', 'CONTACT_', 'CT_')):
                    queue.append((target, total, path + [target]))
        grouped = defaultdict(set)
        for e in entries:
            if e['source'] == root:
                grouped[(*e['role'][:3], e['sign'])].add(e['role'][3])
        for key, phases in grouped.items():
            if {'E', 'F'} <= phases or {'P', 'N'} <= phases:
                flags.append(dict(source=root, layer=key[0], side=key[1],
                                  subsystem=key[2], sign=key[3],
                                  issue='same_sign_contact_to_antagonistic_channels'))
    return dict(flags=flags, entries=entries)


def compose(net, family, config=None):
    import build_network as B
    cfg = deepcopy(DEFAULTS)
    cfg.update(config or {})
    if family not in FAMILIES or hasattr(net, 'central_revision'):
        raise ValueError('Invalid family or circuit composed twice')
    if set(cfg['command']) != {'F_l', 'E_l', 'F_r', 'E_r'}:
        raise ValueError('Exactly four declared commands are required')
    if cfg.get('rg_design') is not None:
        if family=='hunt':raise ValueError('Hunt already loads its source-locked intrinsic design')
        from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
        spec=cfg['rg_design'];nv=spec['nap'];tau=float(spec['tau_h']);mutual=float(spec['mutual_inhibition']);c1=float(spec['c1_gain'])
        if set(nv)!=set(B.NAP) or not 0<tau<=3 or not 0<mutual<=12 or not 0<c1<=1:raise ValueError('Invalid RG seed')
        # Replace only the RG neurons. PF intrinsic parameters stay specific to
        # each family, and no global params dictionary is modified.
        for population in net.net.populations:
            if population['name'].startswith(('RG_E_','RG_F_')):
                population['type']=NonSpikingNeuronWithPersistentSodiumChannel(
                    membrane_capacitance=B.TAU['rg'],membrane_conductance=1.,resting_potential=0.,bias=0.,
                    g_ion=np.array([nv['g_ion']]),e_ion=np.array([nv['e_ion']]),
                    k_m=np.array([nv['k_m']]),slope_m=np.array([nv['slope_m']]),e_m=np.array([nv['e_m']]),
                    k_h=np.array([nv['k_h']]),slope_h=np.array([nv['slope_h']]),e_h=np.array([nv['e_h']]),tau_max_h=np.array([tau]))
        for connection,source,target in connections(net):
            if (source.startswith('RG_') and target.startswith(('InE_','InF_')))or(source.startswith(('InE_','InF_'))and target.startswith('RG_')):
                connection['params']['max_conductance']=mutual
            if source.startswith('CIN_F_') or (source.startswith('RG_F_')and target.startswith('CIN_F_')):
                connection['params']['max_conductance']=c1*mutual
    before = contact_flags(net)
    removed = []
    synergy_mirror = None
    if family == 'syn6':
        # Original left/right NMF columns do not have identical membership.
        # A new symmetric candidate uses the right PF circuit/readout on both
        # sides. The source NPZ and every historical candidate are untouched.
        names = nodes(net)
        index = {n: i for i, n in enumerate(names)}
        def opposite(n):
            return n[:-1] + ('l' if n[-1] == 'r' else 'r') if n.endswith(('_r', '_l')) else n
        def is_pf(n):
            return n.startswith('PF_')
        originals = list(connections(net))
        copied = []
        for c, a, b in originals:
            if (is_pf(b) and b.endswith('_l')) or (is_pf(a) and a.endswith('_l') and b.startswith(('MN_', 'IaIN_'))):
                net.net.connections.remove(c)
            if (is_pf(b) and b.endswith('_r')) or (is_pf(a) and a.endswith('_r') and b.startswith(('MN_', 'IaIN_'))):
                replacement = deepcopy(c)
                replacement['source'] = index[opposite(a)]
                replacement['destination'] = index[opposite(b)]
                copied.append(replacement)
        net.net.connections.extend(copied)
        for n in names:
            if is_pf(n) and n.endswith('_r'):
                net.net.populations[index[opposite(n)]]['type'] = deepcopy(net.net.populations[index[n]]['type'])
        for attr in ('syn_families', 'syn_tau_mult', 'syn_df_channel', 'syn_mixed'):
            getattr(net, attr)['l'] = deepcopy(getattr(net, attr)['r'])
        synergy_mirror = dict(reference_side='r', copied_edges=len(copied),
                              policy='Right PF dynamics, incoming circuit, PF readout and IaIN gates mirrored to left; archived NMF files unchanged')
    # SNS-Toolbox 1.5.2 assigns repeated nonspiking endpoint matrix entries
    # in order (last wins); preserve that exact choice, not their sum.
    seen, deduplicated = set(), []
    for c, a, b in reversed(list(connections(net))):
        key = (a, b)
        if key in seen:
            net.net.connections.remove(c)
            deduplicated.append(dict(source=a, target=b, policy='preserve_compiler_last_assignment'))
        seen.add(key)
    # Remove whole source edges, not zero-gain placeholders. Local MN/reflex
    # arcs and locomotor RG/PF interneuron wiring remain untouched.
    for c, a, b in list(connections(net)):
        descending = a in ('DRIVE', 'BS_ALPHA', 'BS_GAMMA') or a.startswith(('ALPHA_', 'GAMMA_'))
        old_afferent = a.startswith(('Ia_', 'II_', 'Ib_', 'IaIN_', 'IIIN_', 'IIX_', 'IBIN_', 'IBEXC_', 'AFF_', 'HUNT_AF_', 'HEEL_', 'TOE_', 'TOEDF_', 'LBIN_', 'PRESET_'))
        # Toe-dependent local inhibition is not in the selected baseline either.
        # Toe->ankle PF remains a separately named future hypothesis.
        toe_path = a.startswith(('TOE_', 'TOEDF_'))
        if (descending and (central(b) or b.startswith('MN_'))) or (a == 'POSTURE' and b.startswith('RG_')) or (old_afferent and central(b)) or toe_path:
            removed.append(dict(source=a, target=b, gmax_uS=float(c['params']['max_conductance'])))
            net.net.connections.remove(c)
    # Aliases are not actual additional populations.
    aliases = getattr(net, '_aliases', {})
    for a in aliases:
        net.idx.pop(a, None)
    added = []
    def edge(a, b, g, purpose, **extra):
        if g <= 0 or not np.isfinite(g):
            raise ValueError('New edge conductances must be positive and finite')
        if any(x == a and y == b for _, x, y in connections(net)):
            raise ValueError('Duplicate endpoint ' + a + ' -> ' + b)
        net.net.add_connection(B._syn(g, True), a, b)
        added.append(dict(source=a, target=b, gmax_uS=g, purpose=purpose, **extra))
    for side, other in (('l', 'r'), ('r', 'l')):
        for phase in ('F', 'E'):
            cell = f'DESC_RG_{phase}_{side}'
            net._add(cell, cfg['command_tau_s'], net.net)
            net.net.add_input(cell)
            net.inputs.append(f'RG_{phase}_drive_{side}')
            current = cfg['flexor_reference_current_nA' if phase == 'F' else 'extensor_reference_current_nA']
            # At source saturation and Vtarget=0, g*(Erev-Vtarget)=budget.
            edge(cell, f'RG_{phase}_{side}', current / B.E_REV_EXC, 'descending')
            if phase == 'F':
                edge(cell, f'RG_F_{other}', current * cfg['contralateral_flexor_fraction'] / B.E_REV_EXC, 'descending')
                edge(cell, f'InF_{side}', current * (.5 / .45) / B.E_REV_EXC, 'descending')
                for target in (f'V3_E_{side}', f'V3_F_{side}'):
                    if target in net.idx:
                        edge(cell, target, current * (.1 / .45) / B.E_REV_EXC, 'descending')
    # PF feedback is routed only to channels that actually recruit that muscle.
    # Tiny geometric reservations and the tonic bank are not phase primitives.
    carrying = defaultdict(dict)
    for c, a, b in connections(net):
        g = float(c['params']['max_conductance'])
        if b.startswith('MN_') and a.startswith('PF_') and not a.startswith(('PF_IN', 'PF_HUNT_TONIC')) and g > 1e-10 and c['params']['reversal_potential'] > 0:
            carrying[b[3:]][a] = g
    candidates = defaultdict(dict)
    unmapped = []
    for act, muscle in net.muscles.items():
        side, base = muscle.side, muscle.base
        if side not in ('l', 'r'):
            continue
        # Explicit monoarticular timing sources avoid assigning a biarticular
        # muscle one global flexor/extensor identity from its first label.
        rg_phase = ('F' if base in ('iliacus', 'psoas') else
                    'E' if base in ('glut_max1', 'glut_max2', 'glut_max3', 'vas_med', 'vas_int', 'vas_lat', 'soleus') else None)
        stance_force = base in ('glut_max1', 'glut_max2', 'glut_max3', 'vas_med', 'vas_int', 'vas_lat', 'soleus')
        if not carrying[act]:
            unmapped.append(act)
        for aff in ('Ia', 'II', 'Ib'):
            source = f'{aff}_{act}'
            if source not in net.idx:
                raise ValueError('Missing afferent ' + source)
            if rg_phase and (aff != 'Ib' or stance_force):
                candidates[(aff, f'RG_{rg_phase}_{side}')][source] = 1.
            if aff == 'Ib' and not stance_force:
                continue
            for target, weight in carrying[act].items():
                role = pf_role(net, target)
                if aff == 'Ib' and (role is None or role[3] != 'E'):
                    continue
                candidates[(aff, target)][source] = weight
    # Equalize homologous budgets before fan-in normalization. Differences in
    # legacy side-specific W rows must not create independently fitted gains.
    canonical = defaultdict(dict)
    for (aff, target), sources in candidates.items():
        for source, weight in sources.items():
            key = (aff, target[:-1] if target.endswith(('_l', '_r')) else target)
            source_key = source[:-1]
            canonical[key][source_key] = max(canonical[key].get(source_key, 0.), weight)
    for (aff, target), sources in sorted(candidates.items()):
        key = (aff, target[:-1] if target.endswith(('_l', '_r')) else target)
        weights = canonical[key]
        total = sum(weights.values())
        budget = cfg['rg_type_budget_nA' if target.startswith('RG_') else 'pf_type_budget_nA']
        for source in sorted(sources):
            edge(source, target, budget * weights[source[:-1]] / total / B.E_REV_EXC,
                 'afferent', afferent=aff, target_layer='RG' if target.startswith('RG_') else 'PF',
                 source_side=source[-1], stance_gated=aff == 'Ib')
    for side, other in (('l', 'r'), ('r', 'l')):
        if f'HEEL_{side}' not in net.idx:
            raise ValueError('Heel input missing')
        edge(f'HEEL_{side}', f'InE_{side}', cfg['heel_reference_conductance_uS'], 'heel_event')
        edge(f'HEEL_{side}', f'RG_F_{other}', cfg['heel_reference_conductance_uS'], 'heel_event')
        # Contact excites only the extensor-side PF inhibitory relay(s),
        # thereby suppressing opposing flexor PF cells through lamination.
        relays = [n for n in nodes(net) if n.endswith('_' + side) and
                  (n == f'PF_IN_E_{side}' or n.startswith(('PF_IN_HIP-E_', 'PF_IN_KNEE-E_')))]
        for target in relays:
            edge(f'HEEL_{side}', target, cfg['heel_reference_conductance_uS'] / len(relays), 'heel_event')
    net._n_true = len(net.idx)
    net.idx.update(aliases)
    after = contact_flags(net)
    if after['flags']:
        raise ValueError('Unsafe contact fanout: ' + json.dumps(after['flags']))
    net.central_revision = dict(family=family, config=cfg, removed=removed,
                                deduplicated_legacy=deduplicated,
                                synergy_mirror=synergy_mirror,
                                added=added, contact_before=before, contact_after=after,
                                pf_unmapped_muscles=unmapped)


class Runtime:
    """Preallocated command/feedback gate applied to the actual sparse slots."""
    def __init__(self, net, dt):
        self.net, self.dt = net, dt
        self.config = net.central_revision['config']
        self.commands = deepcopy(self.config['command'])
        self.pattern_enable = self.config['pattern_enable']
        self.ports = {n: i for i, n in enumerate(net.inputs)}
        c = net.compiled
        names = nodes(net)
        lookup = {(names[q], names[p]): i for i, (p, q) in enumerate(zip(c.post, c.pre))}
        self.slots, self.base, self.side, self.ib = [], [], [], []
        for e in net.central_revision['added']:
            i = lookup[e['source'], e['target']]
            c.gmax[i] = e['gmax_uS']  # Restore after Hunt install masks.
            if e['purpose'] == 'afferent':
                self.slots.append(i); self.base.append(e['gmax_uS'])
                self.side.append(0 if e['source_side'] == 'l' else 1)
                self.ib.append(e['stance_gated'])
        self.slots = np.array(self.slots, dtype=int)
        self.base = np.array(self.base)
        self.side = np.array(self.side, dtype=int)
        self.ib = np.array(self.ib, dtype=bool)
        self.gains = np.ones(len(self.base))
        self.alpha = np.zeros(2); self.stance = np.zeros(2)
        self.unloaded = np.zeros(2); self.armed = np.zeros(2, dtype=bool)
        self.loaded = np.zeros(2, dtype=bool); self.heel_trace = np.zeros(2)
        self.filtered_alpha = np.zeros(2)
        self.shared_gains = dict(Ia_RG=1., II_RG=1., Ib_RG=1., Ia_PF=1., II_PF=1., Ib_PF=1.)
        self.groups = [e['afferent'] + '_' + e['target_layer'] for e in net.central_revision['added'] if e['purpose'] == 'afferent']
        self.group_names = tuple(self.shared_gains)
        self.group_indices = np.array([self.group_names.index(n) for n in self.groups], dtype=int)
        self.group_values = np.ones(len(self.group_names))
        self.disabled_ports = np.array([i for name, i in self.ports.items() if name == 'DRIVE' or name.startswith(('HUNT_INIT_', 'ALPHA_c_', 'GAMMA_c_'))], dtype=int)
        self.sensory_ports = np.array([i for name, i in self.ports.items() if name.startswith(('Ia_', 'II_', 'Ib_', 'HUNT_AF_', 'HEEL_c_', 'TOE_c_', 'LOAD_c_'))], dtype=int)
        self.sensory_enabled = True
        self.contact_enabled = True
        self.original_step = net.step
        net.step = self.step

    def prepare_inputs(self, u):
        for value in self.commands.values():
            if not 0 <= value <= 1 or not np.isfinite(value):
                raise ValueError('Commands must be normalized values in [0,1]')
        self.filtered_alpha[:] = [np.clip(self.net.compiled.V[self.net.idx[f'DESC_RG_F_{s}']] / self.config['command_voltage_fullscale_mV'], 0., 1.) for s in ('l', 'r')]
        u[self.disabled_ports] = 0.
        if not self.sensory_enabled:
            u[self.sensory_ports] = 0.
        for s in ('l', 'r'):
            for phase in ('F', 'E'):
                u[self.ports[f'RG_{phase}_drive_{s}']] = self.config['command_current_fullscale_nA'] * self.commands[f'{phase}_{s}']
        for i, s in enumerate(('l', 'r')):
            port = self.ports['HEEL_c_' + s]
            raw = max(float(u[port]), 0.)
            if raw <= self.config['heel_off']:
                self.loaded[i] = False; self.unloaded[i] += self.dt
                if self.unloaded[i] >= self.config['heel_rearm_s']:
                    self.armed[i] = True
            elif raw >= self.config['heel_on']:
                if not self.loaded[i] and self.armed[i]:
                    self.heel_trace[i] = 1.
                self.loaded[i] = True; self.armed[i] = False; self.unloaded[i] = 0.
            self.heel_trace[i] *= np.exp(-self.dt / self.config['heel_decay_s'])
            u[port] = self.heel_trace[i] if self.contact_enabled and self.sensory_enabled else 0.
            # Measured LOAD port, not the heel pulse, gates stance Ib.
            load_port = self.ports.get('LOAD_c_' + s)
            self.stance[i] = np.clip(float(u[load_port]) if load_port is not None else raw, 0., 1.)
        self.gains[:] = self.base
        self.gains /= 1. + self.config['alpha_presynaptic_k'] * self.filtered_alpha[self.side]
        self.gains[self.ib] *= self.stance[self.side[self.ib]]
        for i, group in enumerate(self.group_names):
            gain = float(self.shared_gains[group])
            if gain < 0 or not np.isfinite(gain):
                raise ValueError('Invalid shared afferent gain')
            self.group_values[i] = gain
        self.gains *= self.group_values[self.group_indices]
        if not self.sensory_enabled:
            self.gains.fill(0.)
        self.net.compiled.gmax[self.slots] = self.gains
        if hasattr(self.net, 'training_gate'):
            self.net.training_gate = float(self.pattern_enable)
        if hasattr(self.net, 'sync'):
            self.net.sync()
        return u

    def step(self, u):
        self.prepare_inputs(u)
        return self.original_step(u)


def build(family, acts, dt=.002, options=None, config=None, backend='compact'):
    """Build/compile only. No neuron or mechanical step is executed."""
    if family in ('BiW2Lvar+', 'BiW2LvarX'):family='hunt'
    import build_network as B
    import build_network_w2lvar as W
    import build_network_syn6 as Y
    from muscle_map import classify
    from compact_pairwise_eval import CompactPairwiseForward
    from sparse_pairwise_eval import SparsePairwiseForward
    cls = {'s3k': B.SpinalNetwork, 'w2lvar': W.W2LVarnet, 'syn6': Y.Syn6Network}
    if backend not in ('compact', 'sparse'):
        raise ValueError('Only verified graded backends are supported')
    # These runner knobs also scale the measured contact input encoders. Their
    # inherited central projections are removed by compose(), not stacked.
    for key in ('heel_rge', 'toe_rge', 'ib_rge'):
        B.G[key] = 1.
    if family == 'hunt':
        import controller as C
        original, original_backend = C.HuntNetwork, C.SparsePairwiseForward
        class PreparedHunt(original):
            def compile(self, *args, **kwargs):
                compose(self, family, config)
                return super().compile(*args, **kwargs)
        # Dedicated-process factory; never use global builder patching in threads.
        C.HuntNetwork = PreparedHunt
        C.SparsePairwiseForward = CompactPairwiseForward if backend == 'compact' else SparsePairwiseForward
        try:
            net = C.build(acts, dt=dt, options=options)
        finally:
            C.HuntNetwork = original
            C.SparsePairwiseForward = original_backend
    else:
        if family not in cls:
            raise ValueError('Unknown controller family')
        muscles = {a: classify(a) for a in acts}
        if any(m is None for m in muscles.values()):
            raise ValueError('Unknown actuator')
        net = cls[family](muscles=muscles, sides=('l', 'r'), interleg=True)
        compose(net, family, config)
        net.compile(dt)
        net.compiled = (CompactPairwiseForward if backend == 'compact' else SparsePairwiseForward)(net.compiled)
    for key in ('heel_rge', 'toe_rge', 'ib_rge'):
        B.G[key] = 1.
    # Retain Hunt's exact pre-install sparse slot map, including fitted zeros.
    net.command_runtime = Runtime(net, dt)
    net.controller_display_name = DISPLAY_NAMES[family]
    return net
