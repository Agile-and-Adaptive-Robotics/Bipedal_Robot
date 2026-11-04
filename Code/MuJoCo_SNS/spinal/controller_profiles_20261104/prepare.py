"""Construct and audit one revised circuit; no forward, physics, or optimizer."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from controller_profiles_20261104 import profile as P


def audit(net):
    c = net.compiled
    names = P.nodes(net)
    lookup = {(names[pre], names[post]): float(g) for pre, post, g in zip(c.pre, c.post, c.gmax)}
    change = net.central_revision
    missing = [e for e in change['added'] if lookup.get((e['source'], e['target'])) != e['gmax_uS']]
    counts = Counter(e['afferent'] + '_' + e['target_layer'] for e in change['added'] if e['purpose'] == 'afferent')
    expected = {f'{a}_{layer}' for a in ('Ia', 'II', 'Ib') for layer in ('RG', 'PF')}
    load_sources = [f'Ib_{act}' for act, info in net.muscles.items() if info.base in ('glut_max1', 'glut_max2', 'glut_max3', 'vas_med', 'vas_int', 'vas_lat', 'soleus')]
    missing_load_paths = [dict(source=source, layer=layer) for source in load_sources for layer in ('RG', 'PF')
                          if not any(e['source'] == source and e.get('target_layer') == layer for e in change['added'])]
    def mirror(name):
        return name[:-1] + ('r' if name[-1] == 'l' else 'l') if name.endswith(('_l', '_r')) else name
    added = {(e['source'], e['target']): e['gmax_uS'] for e in change['added']}
    mirror_errors = []
    for key, g in added.items():
        other = tuple(mirror(n) for n in key)
        if other != key and (other not in added or abs(added[other] - g) > 1e-14):
            mirror_errors.append(dict(edge=key, mirror=other, g=g, other_g=added.get(other)))
    flag = P.contact_flags(net)
    # Detect accidental duplicates before SNS compilation can merge them.
    endpoints = [(a, b) for _, a, b in P.connections(net)]
    duplicates = [key for key, count in Counter(endpoints).items() if count > 1]
    # Fault-injection of graph entries only: no neural/mechanical integration.
    import build_network as B
    original = next(c for c, source, target in P.connections(net) if source == 'HEEL_r' and target == 'InE_r')
    net.net.add_connection(B._syn(.1, True), 'HEEL_r', 'InF_r')
    injected = net.net.connections[-1]
    try:
        detects_exc = bool(P.contact_flags(net)['flags'])
        old_reversal = original['params']['reversal_potential']
        original['params']['reversal_potential'] = -8.
        injected['params']['reversal_potential'] = -8.
        detects_inh = any(f['sign'] == -1 for f in P.contact_flags(net)['flags'])
    finally:
        original['params']['reversal_potential'] = old_reversal
        net.net.connections.remove(injected)
    pf_original = next(c for c, source, target in P.connections(net) if source == 'HEEL_r' and target.startswith('PF_IN_'))
    pf_name = names[pf_original['destination']]
    antagonist = pf_name.replace('PF_IN_E_', 'PF_IN_F_').replace('HIP-E_', 'HIP-F_').replace('KNEE-E_', 'KNEE-F_')
    if antagonist == pf_name or antagonist not in net.idx:
        raise ValueError('No named antagonistic PF relay for contact guard check')
    net.net.add_connection(B._syn(.1, True), 'HEEL_r', antagonist)
    injected = net.net.connections[-1]
    old_reversal = pf_original['params']['reversal_potential']
    try:
        detects_pf_exc = any(f['layer'] == 'PF' and f['sign'] == 1 for f in P.contact_flags(net)['flags'])
        pf_original['params']['reversal_potential'] = -8.
        injected['params']['reversal_potential'] = -8.
        detects_pf_inh = any(f['layer'] == 'PF' and f['sign'] == -1 for f in P.contact_flags(net)['flags'])
    finally:
        pf_original['params']['reversal_potential'] = old_reversal
        net.net.connections.remove(injected)
    return dict(constructed=True, neural_steps=0, physics_steps=0, optimizer_trials=0,
                neurons=len(names), inputs=len(net.inputs), chemical_slots=len(c.pre),
                new_edges_match_compiled=not missing,
                central_coverage=dict(counts), all_afferent_types_reach_RG_and_PF=expected <= set(counts),
                new_parameters_mirrored=not mirror_errors, mirror_errors=mirror_errors,
                missing_extensor_Ib_RG_or_PF=missing_load_paths,
                same_sign_contact_flags=flag['flags'], duplicate_endpoints=duplicates,
                contact_guard_detects_injected_excitation=detects_exc,
                contact_guard_detects_injected_inhibition=detects_inh,
                contact_guard_detects_PF_excitation=detects_pf_exc,
                contact_guard_detects_PF_inhibition=detects_pf_inh,
                static_pass=not (missing or missing_load_paths or mirror_errors or flag['flags'] or duplicates) and expected <= set(counts) and detects_exc and detects_inh and detects_pf_exc and detects_pf_inh,
                scope='Static graded construction only. No dynamical, training, stage, or spiking acceptance.',
                stage3_postural_premotor_replacement='still_pending_not_created_by_four_command_profile')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', choices=P.FAMILIES, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--hunt-config', type=Path)
    a = p.parse_args()
    model = P.S.parents[2] / 'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    acts = [n.get('name') for n in ET.parse(model).getroot().find('actuator')]
    options = None
    if a.family == 'hunt' and a.hunt_config:
        old = json.loads(a.hunt_config.read_text())
        keys = ('weights_directory', 'timing_scale', 'biarticular_pair_gain', 'learn_reciprocal', 'bilateral_midline_drive', 'hip_abduction_goal')
        options = {k: old[k] for k in keys if k in old}
        options.update(balance=None, use_timing_afferents=False)
    net = P.build(a.family, acts, options=options)
    result = audit(net)
    a.output.mkdir(parents=True, exist_ok=True)
    (a.output / 'WIRING_CHANGES.json').write_text(json.dumps(net.central_revision, indent=2), encoding='utf-8')
    (a.output / 'PROFILE.json').write_text(json.dumps(net.central_revision['config'], indent=2), encoding='utf-8')
    if options is not None:
        (a.output / 'HUNT_OPTIONS.json').write_text(json.dumps(options, indent=2), encoding='utf-8')
    (a.output / 'STATIC_CHECK.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    sources = list(Path(__file__).parent.glob('*.py'))
    sources += [P.S / name for name in ('params.py', 'muscle_map.py', 'build_network.py', 'build_network_w2lvar.py', 'build_network_syn6.py', 'synergy_basis.npz',
                'hunt_training/controller.py', 'functional_gait2392/network.py', 'functional_gait2392/deng2023_joint_pf.py',
                'w2l_cpg/sparse_pairwise_eval.py', 'w2l_cpg/compact_pairwise_eval.py')]
    sources.append(model)
    if options and 'weights_directory' in options:
        sources += list(Path(options['weights_directory']).glob('*.json'))
    (a.output / 'SOURCE_HASHES.json').write_text(json.dumps({str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}, indent=2))
    import params
    (a.output / 'BUILD_SETTINGS.json').write_text(json.dumps(dict(G=params.G, NAP=params.NAP, TAU=params.TAU), indent=2, default=lambda value: value.tolist()))
    print(json.dumps({k: v for k, v in result.items() if k not in ('mirror_errors', 'duplicate_endpoints')}, indent=2))
    if not result['static_pass']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
