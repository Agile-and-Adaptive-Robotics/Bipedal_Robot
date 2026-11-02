"""Prepare, never apply, the optional caller patch and bounded electrical proofs."""
import os
os.environ['AARL_NET'] = 'w2lvar'
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import ast
import copy
import difflib
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import brentq
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
import integrated_hunt_design as I
import renshaw_autogenic_design as RC

HOOK = '''    if cfg.get('renshaw_autogenic_profile'):
        if not cfg.get('integrated_fsa_profile'):
            raise ValueError('The named RC candidate requires its reviewed integrated Hunt source')
        from fsa_design_20261101 import renshaw_autogenic_design as renshaw
        rc_path=Path(cfg['renshaw_autogenic_profile']).resolve()
        rc_proof=renshaw.install(n,renshaw.load_profile(rc_path))
        rc_helper=Path(renshaw.__file__).resolve()
        source_hashes['renshaw_autogenic_design.py']=hashlib.sha256(rc_helper.read_bytes()).hexdigest()
        source_hashes['renshaw_autogenic_profile']=hashlib.sha256(rc_path.read_bytes()).hexdigest()
        capture['renshaw_autogenic_design']={**rc_proof,'profile_file':str(rc_path),'profile_file_sha256':source_hashes['renshaw_autogenic_profile'],'helper_file':str(rc_helper),'helper_file_sha256':source_hashes['renshaw_autogenic_design.py']}
'''


def construct(source, cfg, stage):
    node = next(q for q in ast.parse(source).body if isinstance(q, ast.FunctionDef) and q.name == 'build')
    env = {'H': I.H, 'C': I.C, 'np': np, 'cfg': {**cfg, 'stage': stage},
           'stage': stage, 'capture': {}, 'trace_cells': [], 'event_amplitude': 0.,
           'Path': Path, 'hashlib': hashlib, 'source_hashes': {}}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(I.H/'physical_eval.py'), 'exec'), env)
    acts = [q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')]
    return env['build'](acts, .002, stage != 1), env


def coactive_check(n):
    """Complete retained MN/RC subsystem; all other source inputs silent."""
    c = n.compiled
    mn_ids = [n.idx['MN_'+a] for a in n.targets_data['act_names'].tolist()]
    rc_ids = [n.idx['RC_'+a] for a in n.targets_data['act_names'].tolist()]
    keep = mn_ids+rc_ids;mapping = {old: new for new, old in enumerate(keep)}
    rc_set = set(rc_ids);mn_set = set(mn_ids)
    peer_rows = []
    for post in rc_ids:
        ids = [i for i, (a, b) in enumerate(zip(c.pre, c.post)) if b == post and a in rc_set]
        expected = {q for q in rc_ids if q != post and n.names[q][-1] == n.names[post][-1]}
        if {int(c.pre[i]) for i in ids} != expected or len(ids) != 45:
            raise AssertionError('RC fan-in is not the actual 45 same-side peers')
        if any(c.gmax[i] != .5 or c.reversal[i] != -5 or c.elo[i] != 0 or c.erange[i] != 5 for i in ids):
            raise AssertionError('Unexpected retained peer operator')
        peer_rows.append({'RC': n.names[post], 'peer_inputs': len(ids),
                          'same_side_only': True, 'self_inhibition_present': False})
    row = RC.describe(n)[0]
    g = row['derived_feedback_g_uS']
    u0 = 2.5
    def motor(v):
        b = RC._alpha(v, 0., 5.)
        return (u0+g*b*(-5))/(1+g*b)
    def balance(v):
        u = motor(v)
        return -v+RC._alpha(u, 0., 5.)*(8-v)+45*.5*RC._alpha(v, 0., 5.)*(-5-v)
    vstar = brentq(balance, 0., 8., xtol=1e-14);ustar = motor(vstar)
    pair = Network()
    for old in keep:
        pair.add_neuron(NonSpikingNeuron(membrane_capacitance=float(c.c.c_m[old]),
            membrane_conductance=float(c.c.g_m[old]), resting_potential=float(c.c.V_rest[old]),
            bias=float(c.c.i_b[old])), name=n.names[old])
    active = []
    for i, (pre, post) in enumerate(zip(c.pre, c.post)):
        if int(pre) in mapping and int(post) in mapping and c.gmax[i] != 0:
            # Exactly actual motor collateral, peer inhibition, and new autogenic output.
            active.append(i)
            pair.add_connection(NonSpikingSynapse(max_conductance=float(c.gmax[i]),
                reversal_potential=float(c.reversal[i]), e_lo=float(c.elo[i]),
                e_hi=float(c.elo[i]+c.erange[i])), mapping[int(pre)], mapping[int(post)])
    if len(active) != 4324:
        raise AssertionError(('Unexpected retained subsystem edges', len(active)))
    for old in mn_ids:
        pair.add_input(mapping[old]);pair.add_output(mapping[old])
    for old in rc_ids:
        pair.add_output(mapping[old])
    dense = pair.compile(dt=.002, backend='numpy')
    sparse = I.SparsePairwiseForward(copy.deepcopy(dense))
    current = np.full(92, 2.5);error = 0.;low = np.full(184, np.inf);high = -low
    for _ in range(1500):
        dense.forward(current);sparse.forward(current)
        values = sparse.c.V
        # Abort guards, never clipping or a control clamp.
        if not np.isfinite(values).all() or np.any(values[:92] < -5-1e-10) or np.any(values[:92] > 2.5+1e-10) or np.any(values[92:] < -5-1e-10) or np.any(values[92:] > 8+1e-10):
            raise AssertionError('Coactive fixture left its declared continuous voltage bounds')
        low = np.minimum(low, values);high = np.maximum(high, values)
        error = max(error, float(np.max(abs(dense.V-values))))
    target = np.r_[np.full(92, ustar), np.full(92, vstar)]
    residual = float(np.max(abs(sparse.c.V-target)))
    if error != 0. or residual > 1e-10:
        raise AssertionError(('Actual coactive subsystem mismatch', error, residual))
    # Complete184-state operator, not an isolated-pair stability inference.
    jac = np.diag(-np.asarray(dense.g_m)/np.asarray(dense.c_m))
    for edge in active:
        post, pre = mapping[int(c.post[edge])], mapping[int(c.pre[edge])]
        fraction = RC._alpha(target[pre], c.elo[edge], c.erange[edge])
        if not 0 < fraction < 1:
            raise AssertionError('Coactive Jacobian requires an interior source window')
        jac[post, post] -= c.gmax[edge]*fraction/dense.c_m[post]
        jac[post, pre] += c.gmax[edge]*(c.reversal[edge]-target[post])/(c.erange[edge]*dense.c_m[post])
    sampled = np.eye(184)+.002*jac
    fd = np.zeros_like(sampled);epsilon = 1e-5
    for column in range(184):
        dense.V[:] = target;dense.V[column] += epsilon;dense.forward(current);plus = dense.V.copy()
        dense.V[:] = target;dense.V[column] -= epsilon;dense.forward(current);minus = dense.V.copy()
        fd[:, column] = (plus-minus)/(2*epsilon)
    fd_error = float(np.max(abs(fd-sampled)))
    poles = np.linalg.eigvals(jac);discrete_poles = np.linalg.eigvals(sampled)
    if fd_error > 1e-7:
        raise AssertionError(('Full184-state sampled derivative mismatch', fd_error))
    linearization = {'continuous_poles_real_imag': [[float(z.real), float(z.imag)] for z in poles],
                     'actual_2ms_Euler_poles_real_imag': [[float(z.real), float(z.imag)] for z in discrete_poles],
                     'maximum_continuous_real_per_s': float(poles.real.max()),
                     'minimum_continuous_real_per_s': float(poles.real.min()),
                     'actual_2ms_spectral_radius': float(abs(discrete_poles).max()),
                     'actual_SNS_one_step_FD_max_error': fd_error,
                     'continuous_locally_stable': bool(np.all(poles.real < 0)),
                     'actual_2ms_locally_stable': bool(np.all(abs(discrete_poles) < 1)),
                     'scope': 'Full184-state MN/RC subsystem at declared coactive equilibrium, other controller inputs silent; not full1033-state or body stability.'}
    return {'actual_RC_peer_membership': peer_rows, 'peer_inputs_per_RC': 45,
            'other_same_side_MNs_per_group': 46, 'cross_side_peer_edges': 0,
            'actual_subsystem_neurons': 184, 'actual_subsystem_active_pairs': len(active),
            'coactive_MN_equilibrium_mV': ustar, 'coactive_RC_equilibrium_mV': vstar,
            'MN_delta_from_2p5_mV': ustar-2.5,
            'autogenic_current_nA': g*(vstar/5)*(-5-ustar),
            'total_peer_current_into_each_RC_nA': 45*.5*(vstar/5)*(-5-vstar),
            'RC_layer_voltage_min_mV': float(low[92:].min()),
            'RC_layer_voltage_max_mV': float(high[92:].max()),
            'MN_voltage_min_mV': float(low[:92].min()), 'MN_voltage_max_mV': float(high[:92].max()),
            'dense_sparse_max_error_mV': error, 'equilibrium_error_mV': residual,
            'full184_state_linearization': linearization,
            '_linearization_arrays': {'A_per_s': jac, 'A_Euler_2ms': sampled,
                                      'actual_SNS_one_step_Jacobian_FD': fd,
                                      'continuous_poles_per_s': poles, 'Euler_poles': discrete_poles,
                                      'equilibrium_mV': target, 'names': np.array([n.names[i] for i in keep])},
            'continuous_invariant_bounds_mV': {'MN': [-5., 2.5], 'RC': [-5., 8.]},
            'bounds_scope': 'This constant2.5nA-per-MN fixture only; other full-controller currents are not bounded by this result. Guards abort and never clip.',
            'worst_saturated_RC_Euler_conductance_factor': .002*(1+1+45*.5)/.03,
            'global_monotone_Euler_bound_s': .03/(1+1+45*.5),
            'discrete_bound_caveat': '2ms is not globally a convex update over all saturated peer states; observed guarded trajectory is verified, not a universal discrete-invariance theorem.',
            'scope': 'Complete184-cell MN/RC operator with actual4140 peer paths and92+92 collaterals/outputs. Other local/central/controller inputs silent for this fixture, retained by the production candidate. Not a whole-controller current claim.'}


def main(condition, profile, out, baseline=None):
    out = Path(out);out.mkdir(parents=True, exist_ok=True)
    caller = I.H/'physical_eval.py';actual = caller.read_text();digest = I.sha(caller)
    original = Path(baseline).read_text() if baseline else actual
    marker = '    native_step=n.step\n'
    if original.count(marker) != 1 or 'renshaw_autogenic_profile' in original:
        raise ValueError('Caller insertion precondition no longer matches')
    proposed = original.replace(marker, HOOK+marker)
    if baseline and proposed != actual:
        raise AssertionError('Installed runtime differs from exactly the reviewed patch')
    patch = ''.join(difflib.unified_diff(original.splitlines(True), proposed.splitlines(True),
        fromfile='a/Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py',
        tofile='b/Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py'))
    patch_path = I.HERE/'physical_eval_autogenic_rc_v1.patch';patch_path.write_text(patch, encoding='utf-8')
    cfg = json.loads(Path(condition).read_text());rows = []
    for stage in (1, 2, 3, 4):
        old, old_env = construct(original, cfg, stage)
        unchanged, none_env = construct(proposed, {**cfg, 'renshaw_autogenic_profile': None}, stage)
        if I.GP.field_digest(old.compiled.c) != I.GP.field_digest(unchanged.compiled.c) or old_env['capture'].keys() != none_env['capture'].keys() or none_env['source_hashes']:
            raise AssertionError('DefaultNone constructor/provenance changed')
        for field in ('V', 'V_last', 'b_gate', 'c_gate'):
            if hasattr(old.compiled.c, field) and not np.array_equal(getattr(old.compiled.c, field), getattr(unchanged.compiled.c, field)):
                raise AssertionError('DefaultNone initial state changed: '+field)
        default_error = 0.
        x = np.zeros(len(old.inputs))
        for _ in range(128):
            old.compiled.forward(x);unchanged.compiled.forward(x)
            default_error = max(default_error, float(np.max(abs(old.compiled.c.V-unchanged.compiled.c.V))))
        if default_error != 0.:
            raise AssertionError('DefaultNone electrical trajectory changed')
        new, env = construct(proposed, {**cfg, 'renshaw_autogenic_profile': str(Path(profile).resolve())}, stage)
        for key in ('original_g', 'pf_gain', 'ff_gain', 'bi_gain', 'reflex_gain'):
            if not np.array_equal(getattr(old, key), getattr(new, key)):
                raise AssertionError('Existing gain buffer changed')
        with I.unbound_electrical(new):
            for command in (0., .5, 1., 0.):
                new.apply_weights(command);RC.verify(new);I.verify(new, full=True)
        rows.append({'stage': stage, 'default_None_exact': True,
                     'default128step_max_error_mV': default_error,
                     'existing_gain_buffers_exact': True, 'RC_outputs': RC.verify(new),
                     'captured_RC': env['capture']['renshaw_autogenic_design'],
                     'captured_source_hashes': env['source_hashes']})
    if I.sha(caller) != digest:
        raise AssertionError('Current runtime changed during preparation')
    proof = {'actual_runtime_sha256_unchanged': digest, 'patch_applied': bool(baseline),
             'actual_on_disk_constructor_matches_proposed_source': bool(baseline),
             'baseline_source': str(baseline) if baseline else str(caller),
             'baseline_source_sha256': I.sha(Path(baseline)) if baseline else digest,
             'patch_file': str(patch_path), 'patch_sha256': I.sha(patch_path),
             'application_requires_root_confirmation_after_ground_v2_study': not bool(baseline),
             'source_precondition': 'Exact recorded baseline hash; unique native_step insertion point; reviewed integrated Hunt source; actual dormant92 RC outputs.',
             'actual_constructor_proofs': rows,
             'complete_coactive_operator': coactive_check(new),
             'physics_or_body_steps': 0}
    arrays = proof['complete_coactive_operator'].pop('_linearization_arrays')
    np.savez_compressed(out/'COACTIVE_LINEARIZATION.npz', **arrays)
    proof['preparation_source_sha256'] = I.sha(Path(__file__))
    (out/'OPTIONAL_PATCH_PROOF.json').write_text(json.dumps(proof, indent=2, sort_keys=True), encoding='utf-8')
    print(json.dumps({'patch': str(patch_path), 'proof': str(out/'OPTIONAL_PATCH_PROOF.json'),
                      'coactive_delta_mV': proof['complete_coactive_operator']['MN_delta_from_2p5_mV']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('condition');parser.add_argument('profile');parser.add_argument('out');parser.add_argument('--baseline')
    args = parser.parse_args();main(args.condition, args.profile, args.out, args.baseline)
