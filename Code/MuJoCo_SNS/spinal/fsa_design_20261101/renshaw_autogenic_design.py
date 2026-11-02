"""Opt-in autogenic Renshaw candidate; no historical source/runtime edits.

The contract is an isolated, closed two-cell motor-collateral operator, not
afferent feedback or a whole-network physiological-gain estimate.
"""
import hashlib
import json
from pathlib import Path
from types import MethodType
import numpy as np

SCHEMA = 'aarl_hunt_autogenic_renshaw_v1'


def _alpha(value, low, span):
    return float(np.clip((value-low)/span, 0., 1.))


def describe(net, delta=-.1):
    if not -.5 <= delta < 0:
        raise ValueError('Declared isolated target must lie in [-0.5,0) mV')
    c = net.compiled
    rows = []
    for act in net.targets_data['act_names'].tolist():
        mn, rc = net.idx['MN_'+act], net.idx['RC_'+act]
        forward = np.flatnonzero((c.pre == mn) & (c.post == rc))
        feedback = np.flatnonzero((c.pre == rc) & (c.post == mn))
        if len(forward) != 1 or len(feedback) != 1:
            raise ValueError('Expected one actual motor collateral and feedback slot: '+act)
        a, b = int(forward[0]), int(feedback[0])
        if c.gmax[a] <= 0 or c.reversal[a] <= 2.5 or c.reversal[b] >= 2.5:
            raise ValueError('Invalid actual recurrent-inhibition operator: '+act)
        if hasattr(c.c, 'g_ion') and np.any(c.c.g_ion[:, [mn, rc]]):
            raise ValueError('Isolated passive derivation does not cover ionic MN/RC')
        target = 2.5+delta
        drive = c.gmax[a]*_alpha(target, c.elo[a], c.erange[a])
        rc_target = (c.c.g_m[rc]*c.c.V_rest[rc]+c.c.i_b[rc]+drive*c.reversal[a])/(c.c.g_m[rc]+drive)
        fraction = _alpha(rc_target, c.elo[b], c.erange[b])
        if fraction <= 0:
            raise ValueError('Actual collateral cannot activate recurrent output: '+act)
        gain = -c.c.g_m[mn]*delta/(fraction*(target-c.reversal[b]))
        row = dict(actuator=act, MN='MN_'+act, RC='RC_'+act,
                   collateral_g_uS=float(c.gmax[a]), collateral_E_mV=float(c.reversal[a]),
                   collateral_low_mV=float(c.elo[a]), collateral_span_mV=float(c.erange[a]),
                   feedback_E_mV=float(c.reversal[b]), feedback_low_mV=float(c.elo[b]),
                   feedback_span_mV=float(c.erange[b]),
                   MN_G_uS=float(c.c.g_m[mn]), MN_C_microF=float(c.c.c_m[mn]),
                   MN_rest_mV=float(c.c.V_rest[mn]), MN_bias_nA=float(c.c.i_b[mn]),
                   RC_G_uS=float(c.c.g_m[rc]), RC_C_microF=float(c.c.c_m[rc]),
                   RC_rest_mV=float(c.c.V_rest[rc]), RC_bias_nA=float(c.c.i_b[rc]),
                   reference_MN_mV=2.5, target_MN_mV=float(target),
                   derived_RC_mV=float(rc_target), derived_feedback_g_uS=float(gain))
        rows.append(row)
    if len(rows) != 92:
        raise ValueError('This named candidate requires the actual 92-muscle Hunt family')
    by_act = {r['actuator']: r for r in rows}
    for name, row in by_act.items():
        other = name[:-2]+('_l' if name.endswith('_r') else '_r')
        if other not in by_act:
            raise ValueError('Unpaired actual muscle: '+name)
        fields = set(row)-{'actuator', 'MN', 'RC'}
        if any(row[k] != by_act[other][k] for k in fields):
            raise ValueError('Source bilateral RC parameters differ: '+name)
    return rows


def make_profile(net, delta=-.1):
    return {'schema': SCHEMA, 'target_delta_mV': float(delta),
            'source_and_derived_operators': describe(net, delta),
            'scope': 'Isolated motor-collateral recurrent inhibition at MN2.5mV; '
                     'other incoming circuits silent for unit derivation, retained in production. '
                     'Not an afferent; retained in deafferented protocol. No body pass inherited.'}


def load_profile(profile):
    p = json.loads(Path(profile).read_text()) if isinstance(profile, (str, Path)) else profile
    if set(p) != {'schema', 'target_delta_mV', 'source_and_derived_operators', 'scope'} or p['schema'] != SCHEMA:
        raise ValueError('Unknown Renshaw profile schema/key')
    return p


def install(net, profile=None):
    if profile is None:
        return {'applied': False}
    p = load_profile(profile)
    if hasattr(net, 'renshaw_autogenic_design') or 'apply_weights' in net.__dict__:
        raise ValueError('Refuse duplicate/unknown command wrapper')
    rows = describe(net, p['target_delta_mV'])
    if rows != p['source_and_derived_operators']:
        raise ValueError('Actual source RC operator differs from profile')
    c = net.compiled
    ids = np.array([np.flatnonzero((c.pre == net.idx[r['RC']]) &
                                  (c.post == net.idx[r['MN']]))[0] for r in rows], dtype=int)
    if np.any(c.gmax[ids]):
        raise ValueError('Candidate specifically repairs the known masked RC outputs')
    native = net.apply_weights.__func__
    net.renshaw_autogenic_gain = np.array([r['derived_feedback_g_uS'] for r in rows])
    net.renshaw_autogenic_design = {'profile': p, 'indices': ids,
                                  'native_apply_function': native}

    def apply_weights(self, command):
        native(self, command)
        self.compiled.gmax[ids] = self.renshaw_autogenic_gain
        self.compiled.c.g_max_non[self.compiled.post[ids], self.compiled.pre[ids]] = self.renshaw_autogenic_gain
        verify(self)

    net.apply_weights = MethodType(apply_weights, net)
    c.gmax[ids] = net.renshaw_autogenic_gain
    c.c.g_max_non[c.post[ids], c.pre[ids]] = net.renshaw_autogenic_gain
    return verify(net)


def verify(net):
    state = net.renshaw_autogenic_design
    ids = state['indices']
    expected = np.array([r['derived_feedback_g_uS'] for r in state['profile']['source_and_derived_operators']])
    c = net.compiled
    if not np.array_equal(net.renshaw_autogenic_gain, expected):
        raise ValueError('Authoritative recurrent-feedback gain overwritten')
    if not np.array_equal(c.gmax[ids], expected) or not np.array_equal(c.c.g_max_non[c.post[ids], c.pre[ids]], expected):
        raise ValueError('Actual RC conductance/cache differs from its contract')
    return {'applied': True, 'autogenic_RC_outputs': len(ids),
            'afferent_protocol_mask_independent': True,
            'profile_sha256': hashlib.sha256(json.dumps(state['profile'], sort_keys=True).encode()).hexdigest()}
