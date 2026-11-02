"""Read-only selected Hunt Renshaw routing/buffer audit; no body integration."""
import argparse
import ast
import json
import os
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
os.environ['AARL_NET'] = 'w2lvar'
import integrated_hunt_design as I


def audit(condition):
    cfg = json.loads(Path(condition).read_text())
    profile = I.load_profile(cfg['integrated_fsa_profile'])
    acts = [q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')]
    caller = I.H / 'physical_eval.py'
    tree = ast.parse(caller.read_text())
    node = next(q for q in tree.body if isinstance(q, ast.FunctionDef) and q.name == 'build')
    env = {'H': I.H, 'C': I.C, 'np': np, 'cfg': cfg, 'stage': cfg['stage'],
           'capture': {}, 'trace_cells': [], 'event_amplitude': 0.}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(caller), 'exec'), env)
    n = env['build'](acts, .002, True)
    c = n.compiled
    mechanics = Path(condition).with_suffix('.mechanics.npz')
    with np.load(mechanics, allow_pickle=True) as saved_mechanics:
        names_exact = np.array_equal(saved_mechanics['names'], n.names[:len(c.c.V)])
    if not names_exact:
        raise AssertionError('Actual compiled name order differs from saved selected mechanics')
    rc_out = np.array([i for i, pre in enumerate(c.pre)
                       if n.names[pre].startswith('RC_')], dtype=int)
    rc_mn = np.array([i for i in rc_out
                      if n.names[c.post[i]].startswith('MN_')], dtype=int)
    groups = sorted({n.names[c.post[i]].split('_', 1)[0] for i in rc_out})
    routes = {}
    for group in groups:
        ids = np.array([i for i in rc_out
                        if n.names[c.post[i]].split('_', 1)[0] == group])
        routes[group] = {'allocated_pairs': len(ids),
                         'original_nonzero': int(np.count_nonzero(n.original_g[ids])),
                         'live_nonzero': int(np.count_nonzero(c.gmax[ids])),
                         'live_gmax_range_uS': [float(c.gmax[ids].min()),
                                                float(c.gmax[ids].max())]}
    commands = []
    with I.unbound_electrical(n):
        for command in (0., .5, 1., 0.):
            n.apply_weights(command)
            commands.append({'command': command,
                             'RC_to_MN_live_nonzero': int(np.count_nonzero(c.gmax[rc_mn])),
                             'RC_to_MN_dense_nonzero': int(np.count_nonzero(
                                 c.c.g_max_non[c.post[rc_mn], c.pre[rc_mn]])),
                             'RC_outgoing_active_targets': dict(Counter(
                                 n.names[c.post[i]].split('_', 1)[0]
                                 for i in rc_out if c.gmax[i] != 0))})
    saved = {}
    for filename in ('pf_weights.json', 'bi_weights.json', 'mn_weights.json',
                     'reflex_weights.json'):
        path = Path(cfg['weights_directory']) / filename
        edges = json.loads(path.read_text()).get('edges', {})
        rc = {key: value for key, value in edges.items()
              if key.startswith('RC_') and ' -> MN_' in key}
        saved[filename] = {'RC_to_MN_keys': len(rc),
                           'nonzero_keys': sum(v != 0 for v in rc.values()),
                           'value_range': ([min(rc.values()), max(rc.values())]
                                           if rc else None)}
    return {'condition': str(Path(condition).resolve()),
            'profile': cfg['integrated_fsa_profile'],
            'actual_physical_eval_neural_constructor_used': True,
            'physical_eval_sha256': I.sha(caller),
            'saved_mechanics': str(mechanics),
            'all_actual_compiled_names_exact_to_saved_mechanics': names_exact,
            'scope': 'Neural construction/command application only; PT not body-bound. '
                     'MN-only POSTURE mask cannot alter any RC-source slot.',
            'neurons': len(c.c.V), 'name_entries_including_logging_aliases': len(n.names),
            'RC_neurons': sum(s.startswith('RC_') for s in n.names),
            'RC_outgoing': routes, 'command_checks': commands,
            'RC_to_MN_in_reflex_edge_set': len(set(rc_mn) & set(n.reflex_edges)),
            'RC_to_MN_original_nonzero': int(np.count_nonzero(n.original_g[rc_mn])),
            'RC_to_MN_reflex_buffer_nonzero': int(np.count_nonzero(n.reflex_gain[rc_mn])),
            'saved_weight_keys': saved,
            'interpretation': 'Direct recurrent RC-to-MN inhibition is suppressed by '
                              'Hunt MN clearing/restoration policy, although other '
                              'Renshaw outputs may remain active. Cell presence and '
                              'local-pool FSA pass do not certify complete lamination.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('condition')
    parser.add_argument('--out')
    args = parser.parse_args()
    result = audit(args.condition)
    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
