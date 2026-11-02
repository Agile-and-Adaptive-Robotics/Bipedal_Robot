"""Opt-in, endpoint-identified four-coefficient Hunt ankle output checkpoint.

No fitting, model stepping, or changes to historical weights. Conductances are uS.
The integrated local/timing rules are installed before changing motor buffers.
"""
import argparse
import copy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import integrated_hunt_design as I

BUFFERS = ('pf_gain', 'ff_gain', 'reflex_gain', 'bi_gain')
BOUNDS = {'soleus': (0., .15), 'tib_ant': (.3, 1.2)}


def base_from_best(best):
    best = Path(best).resolve()
    b = json.loads(best.read_text())
    cfg_path = Path(b['trial']['diagnostics']['configuration'])
    cfg = json.loads(cfg_path.read_text())
    if cfg.get('motor_scale', 1.) != 1 or cfg.get('reflex_scale', 1.) != 1 or cfg.get('motor_scales'):
        raise ValueError('Base has additional caller scaling; explicitly compose before using this API')
    p = I.load_profile(cfg['integrated_fsa_profile'])
    acts = [q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')]
    n = I.build(acts, options=I.condition_options(cfg), profile=p, protocol_stage=cfg['stage'])
    return n, p, cfg, {'best': str(best), 'best_sha256': I.sha(best),
                      'base_configuration': str(cfg_path), 'base_configuration_sha256': I.sha(cfg_path),
                      'base_profile': cfg['integrated_fsa_profile'],
                      'base_profile_sha256': I.sha(cfg['integrated_fsa_profile']),
                      'trial': b['trial']['number']}


def endpoints(net):
    c = net.compiled
    result = {}
    for muscle, sign in (('soleus', 'P'), ('tib_ant', 'N')):
        for side in ('r', 'l'):
            pre, post = f'PF_FUNC_ANKLE_Z_{sign}_{side}', f'MN_{muscle}_{side}'
            ids = [i for i, (a, b) in enumerate(zip(c.pre, c.post))
                   if net.names[int(a)] == pre and net.names[int(b)] == post]
            if len(ids) != 1 or ids[0] not in net.pf_edges:
                raise ValueError('Expected exactly one actual PF output edge: ' + pre + ' -> ' + post)
            result[muscle + '_' + side] = dict(index=ids[0], source=pre, target=post,
                                              base_uS=float(net.ff_gain[ids[0]]))
    for muscle in BOUNDS:
        if result[muscle + '_r']['base_uS'] != result[muscle + '_l']['base_uS']:
            raise ValueError('Base output coefficients are not tied')
    return result


def set_pair(net, soleus_uS, tib_ant_uS):
    values = dict(soleus=float(soleus_uS), tib_ant=float(tib_ant_uS))
    for muscle, value in values.items():
        lo, hi = BOUNDS[muscle]
        if not np.isfinite(value) or not lo <= value <= hi:
            raise ValueError(f'{muscle} conductance must be within [{lo},{hi}] uS')
    before = {q: getattr(net, q).copy() for q in BUFFERS}
    edges = endpoints(net)
    for name, row in edges.items():
        row['requested_uS'] = values[name.rsplit('_', 1)[0]]
        net.ff_gain[row['index']] = row['requested_uS']
    allowed = {v['index'] for v in edges.values()}
    changed = {}
    for field in BUFFERS:
        ids = np.flatnonzero(getattr(net, field) != before[field]).tolist()
        if (field != 'ff_gain' and ids) or not set(ids) <= allowed:
            raise AssertionError('Unexpected authoritative coefficient changed')
        changed[field] = ids
    return edges, changed


def command_proof(net, edges):
    rows = []
    with I.unbound_electrical(net):
        for command in (0., .5, 1., 0.):
            u = np.zeros(len(net.inputs))
            u[net.ports['DRIVE']] = 2.5 * command
            net.step(u)
            I.verify(net, full=True)
            for name, edge in edges.items():
                i = edge['index']; c = net.compiled
                expected = command * edge['requested_uS']
                actual = float(c.gmax[i])
                if actual != expected or float(c.c.g_max_non[c.post[i], c.pre[i]]) != expected:
                    raise AssertionError('Actual command/current cache did not consume authoritative ff_gain')
                on, off = I.clone(c), I.clone(c)
                a, b = int(c.pre[i]), int(c.post[i])
                on.c.V[a] = off.c.V[a] = c.elo[i] + c.erange[i]
                on.c.V[b] = off.c.V[b] = 2.5
                off.gmax[i] = 0.; off.c.g_max_non[b, a] = 0.
                on.forward(np.zeros(len(net.inputs))); off.forward(np.zeros(len(net.inputs)))
                current = float((on.V[b] - off.V[b]) / on.c.time_factor_membrane[b])
                error = abs(current - expected * (c.reversal[i] - 2.5))
                if error > 1e-10:
                    raise AssertionError('Isolated motor-current operator mismatch')
                rows.append(dict(command=command, edge=name, gmax_uS=actual,
                                 current_error_nA=error, exact_zero=(actual == 0.)))
    return rows


def generate(best, out, soleus_uS, tib_ant_uS):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    if any((out / name).exists() for name in ('PROFILE.json', 'CONFIGURATION.json', 'PARAMETERS.npz', 'PROOF.json')):
        raise FileExistsError('Fresh output directory required; no checkpoint overwrites')
    n, p, cfg, provenance = base_from_best(best)
    before = {q: getattr(n, q).copy() for q in BUFFERS}
    edges, changed = set_pair(n, soleus_uS, tib_ant_uS)
    p['source_sha256'][str(Path(__file__).resolve())] = I.sha(__file__)
    p = I.snapshot_parameters(n, out / 'PARAMETERS.npz', p)
    (out / 'PROFILE.json').write_text(json.dumps(p, indent=2))
    cfg['integrated_fsa_profile'] = str(out / 'PROFILE.json')
    (out / 'CONFIGURATION.json').write_text(json.dumps(cfg, indent=2))
    proof = command_proof(n, edges)
    acts = [q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')]
    fresh = I.build(acts, options=I.condition_options(cfg), profile=out / 'PROFILE.json', protocol_stage=cfg['stage'])
    for field in BUFFERS:
        if not np.array_equal(getattr(n, field), getattr(fresh, field)):
            raise AssertionError('Endpoint checkpoint reload differs: ' + field)
    expected_changes = sum(row['base_uS'] != row['requested_uS'] for row in edges.values())
    if len(changed['ff_gain']) != expected_changes:
        raise AssertionError('Unexpected coefficient delta count')
    result = dict(provenance=provenance, generator_source=str(Path(__file__).resolve()),
                  generator_sha256=I.sha(__file__), bounds_uS=BOUNDS, edges=edges,
                  changed_authoritative_indices=changed, expected_changed_count=expected_changes,
                  command_proof=proof, fresh_reload_command_proof=command_proof(fresh, edges),
                  unchanged_buffer_elements={q:int(np.sum(getattr(n,q)==before[q])) for q in BUFFERS},
                  source_physics_unchanged=True, inputs_unchanged=True, body_steps=0, fitting=False,
                  full_body_acceptance=False, checkpoint_sha256=p['parameter_checkpoint']['sha256'])
    (out / 'PROOF.json').write_text(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--best', required=True); parser.add_argument('--out', required=True)
    parser.add_argument('--soleus', type=float, required=True); parser.add_argument('--tib-ant', type=float, required=True)
    args = parser.parse_args()
    result = generate(args.best, args.out, args.soleus, args.tib_ant)
    print(json.dumps(dict(changed=result['expected_changed_count'], body_steps=0, proof=str(Path(args.out)/'PROOF.json'))))
