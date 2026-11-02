"""Bounded electrical proofs for the separate opt-in motor-collateral design."""
import os
os.environ['AARL_NET'] = 'w2lvar'
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import ast
import copy
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
import integrated_hunt_design as I
import renshaw_autogenic_design as R


def selected(cfg, stage):
    caller = I.H/'physical_eval.py'
    tree = ast.parse(caller.read_text())
    node = next(q for q in tree.body if isinstance(q, ast.FunctionDef) and q.name == 'build')
    env = {'H': I.H, 'C': I.C, 'np': np, 'cfg': {**cfg, 'stage': stage},
           'stage': stage, 'capture': {}, 'trace_cells': [], 'event_amplitude': 0.}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(caller), 'exec'), env)
    acts = [q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')]
    return env['build'](acts, .002, stage != 1)


def pair_check(row):
    pair = Network()
    for name in ('MN', 'RC'):
        pair.add_neuron(NonSpikingNeuron(
            membrane_capacitance=row[name+'_C_microF'],
            membrane_conductance=row[name+'_G_uS'],
            resting_potential=row[name+'_rest_mV'], bias=row[name+'_bias_nA']), name=name)
    pair.add_connection(NonSpikingSynapse(max_conductance=row['collateral_g_uS'],
        reversal_potential=row['collateral_E_mV'], e_lo=row['collateral_low_mV'],
        e_hi=row['collateral_low_mV']+row['collateral_span_mV']), 'MN', 'RC')
    pair.add_connection(NonSpikingSynapse(max_conductance=row['derived_feedback_g_uS'],
        reversal_potential=row['feedback_E_mV'], e_lo=row['feedback_low_mV'],
        e_hi=row['feedback_low_mV']+row['feedback_span_mV']), 'RC', 'MN')
    pair.add_input('MN');pair.add_output('MN');pair.add_output('RC')
    dense = pair.compile(dt=.002, backend='numpy')
    sparse = I.SparsePairwiseForward(copy.deepcopy(dense))
    current = np.array([row['MN_G_uS']*(2.5-row['MN_rest_mV'])-row['MN_bias_nA']])
    error = 0.
    for _ in range(1500):
        dense.forward(current);sparse.forward(current)
        error = max(error, float(np.max(abs(dense.V-sparse.c.V))))
    target = np.array([row['target_MN_mV'], row['derived_RC_mV']])
    residual = float(np.max(abs(dense.V-target)))
    if error != 0. or residual > 1e-10:
        raise AssertionError(('Closed recurrent pair mismatch', error, residual))
    final = dense.V.copy()
    u, v = target
    gm, gr = row['MN_G_uS'], row['RC_G_uS']
    cm, cr = row['MN_C_microF'], row['RC_C_microF']
    g, h = row['derived_feedback_g_uS'], row['collateral_g_uS']
    ef, ec = row['feedback_E_mV'], row['collateral_E_mV']
    am = R._alpha(u, row['collateral_low_mV'], row['collateral_span_mV'])
    ar = R._alpha(v, row['feedback_low_mV'], row['feedback_span_mV'])
    dm, dr = 1/row['collateral_span_mV'], 1/row['feedback_span_mV']
    if not 0 < am < 1 or not 0 < ar < 1:
        raise AssertionError('Interior small-signal derivation does not cover this branch')
    jac = np.array([[-(gm+g*ar)/cm, g*dr*(ef-u)/cm],
                    [h*dm*(ec-v)/cr, -(gr+h*am)/cr]])
    discrete = np.eye(2)+.002*jac
    b = np.array([1/cm, 0.])
    step_fd = np.zeros((2, 2));epsilon = 1e-5
    for j in range(2):
        delta = np.eye(2)[j]*epsilon
        dense.V[:] = target+delta;dense.forward(current);plus = dense.V.copy()
        dense.V[:] = target-delta;dense.forward(current);minus = dense.V.copy()
        step_fd[:, j] = (plus-minus)/(2*epsilon)
    derivative_error = float(np.max(abs(step_fd-discrete)))
    dense.V[:] = target;linear = np.zeros(2);pulse_error = 0.
    for k in range(150):
        pulse = .001 if k < 50 else 0.
        dense.forward(current+pulse);linear = discrete@linear+.002*b*pulse
        pulse_error = max(pulse_error, float(np.max(abs(dense.V-(target+linear)))))
    poles = np.linalg.eigvals(jac);sampled_poles = np.linalg.eigvals(discrete)
    if derivative_error > 1e-8 or pulse_error > 1e-6 or np.any(poles.real >= 0) or np.any(abs(sampled_poles) >= 1):
        raise AssertionError(('Small-signal proof failed', derivative_error, pulse_error))
    transfer = {'A_per_s': jac.tolist(), 'B_mV_per_s_per_nA': b.tolist(),
                'C': [[1., 0.]], 'D': 0.,
                'numerator_descending_s': [float(b[0]), float(-b[0]*jac[1, 1])],
                'denominator_descending_s': [1., float(-np.trace(jac)), float(np.linalg.det(jac))],
                'DC_gain_mV_per_nA': float(-np.linalg.solve(jac, b)[0]),
                'continuous_poles_per_s': [[float(z.real), float(z.imag)] for z in poles],
                'actual_2ms_Euler_poles': [[float(z.real), float(z.imag)] for z in sampled_poles],
                'actual_SNS_one_step_Jacobian_FD_error': derivative_error,
                'pulse_current_nA': .001, 'pulse_duration_s': .1,
                'nonlinear_SNS_vs_linear_max_error_mV': pulse_error,
                'scope': 'Isolated recurrent pair only, interior saturation branch; '
                         'not the full SNS, 105-coordinate model, or hybrid gait stability.'}
    return {'source_actuator': row['actuator'], 'actual_source_parameters': row,
            'external_MN_baseline_current_nA': float(current[0]),
            'MN_RC_final_mV': final.tolist(), 'target_error_mV': residual,
            'dense_sparse_max_error_mV': error,
            'local_linearization': transfer,
            'unit_fixture': 'Actual two-cell leak/collateral/feedback law, constant baseline input; '
                            'both voltages free, no voltage clamp. Other production circuits retained.'}


def run(condition, out):
    cfg = json.loads(Path(condition).read_text());out = Path(out);out.mkdir(parents=True, exist_ok=True)
    n = selected(cfg, 3);c = n.compiled
    digest = I.GP.field_digest(c.c)
    R.install(n, None)
    if digest != I.GP.field_digest(c.c) or 'apply_weights' in n.__dict__:
        raise AssertionError('Default None changed actual network')
    profile = R.make_profile(n)
    before = {key: np.asarray(getattr(c.c, key)).copy()
              for key in I.GP.FIELDS+I.GP.ION_FIELDS if hasattr(c.c, key)}
    gains = {key: getattr(n, key).copy() for key in ('original_g', 'pf_gain', 'ff_gain', 'bi_gain', 'reflex_gain')}
    prior = c.gmax.copy();R.install(n, profile)
    ids = n.renshaw_autogenic_design['indices'];other = np.setdiff1d(np.arange(len(c.gmax)), ids)
    if not np.array_equal(np.flatnonzero(c.gmax != prior), np.sort(ids)):
        raise AssertionError('Changes extend beyond the92 RC output slots')
    for key, value in before.items():
        if key != 'g_max_non' and not np.array_equal(value, getattr(c.c, key)):
            raise AssertionError('Non-conductance field changed: '+key)
    for key, value in gains.items():
        if not np.array_equal(value, getattr(n, key)):
            raise AssertionError('Existing authoritative gain buffer changed: '+key)
    rows = []
    baseline = selected(cfg, 3)
    with I.unbound_electrical(n), I.unbound_electrical(baseline):
        for command in (0., .5, 1., 0.):
            baseline.apply_weights(command);n.apply_weights(command)
            if not np.array_equal(c.gmax[other], baseline.compiled.gmax[other]):
                raise AssertionError('Other command-dependent routes changed')
            rows.append({'command': command, 'new_RC': R.verify(n), 'old_integrated': I.verify(n, full=True)})
    rejected = {}
    edge = ids[0];post, pre = c.post[edge], c.pre[edge]
    for name, array, index in (('sparse', c.gmax, edge),
                               ('dense', c.c.g_max_non, (post, pre)),
                               ('authoritative_buffer', n.renshaw_autogenic_gain, 0)):
        value = array[index];array[index] += .001
        try:
            R.verify(n);rejected[name] = False
        except ValueError:
            rejected[name] = True
        finally:
            array[index] = value
        if not rejected[name]:
            raise AssertionError('Corrupted RC value accepted: '+name)
    dense = copy.deepcopy(c.c)
    x = np.zeros(len(n.inputs));max_error = 0.
    for _ in range(128):
        dense.forward(x);c.forward(x)
        max_error = max(max_error, float(np.max(abs(dense.V-c.c.V))))
    if max_error != 0.:
        raise AssertionError(('Actual whole neural dense/sparse differs', max_error))
    deaff = selected(cfg, 1);R.install(deaff, profile)
    for command in (0., .5, 1., 0.):
        deaff.apply_weights(command);R.verify(deaff);I.verify(deaff, full=True)
    unique = {}
    for row in profile['source_and_derived_operators']:
        key = tuple((k, v) for k, v in row.items() if k not in ('actuator', 'MN', 'RC'))
        unique.setdefault(key, row)
    paired = [pair_check(row) for row in unique.values()]
    result = {'schema': 'aarl_hunt_autogenic_renshaw_proof_v1', 'condition': str(condition),
              'helper_sha256': I.sha(Path(R.__file__)), 'proof_source_sha256': I.sha(Path(__file__)),
              'physical_eval_sha256': I.sha(I.H/'physical_eval.py'),
              'actual_neurons': len(c.c.V), 'only92_RC_gains_changed': True,
              'all_other_parameters_and_gain_buffers_exact': True,
              'other_routes_exact_after_all_commands': True, 'default_None_exact': True,
              'actual128step_dense_sparse_max_error_mV': max_error,
              'tamper_rejected': rejected,
              'command_checks': rows, 'deafferented_collateral_retained': True,
              'unique_actual_operator_pair_checks': paired,
              'scope': 'Separate opt-in neural wiring candidate; not applied to main runtime. '
                       'No body, training, physics, or accepted-controller changes. '
                       'Whole-circuit recurrent interaction/phase/body impact not established.'}
    # Exclude environment-dependent elapsed times; a fresh process must reproduce these bytes.
    (out/'PROFILE.json').write_text(json.dumps(profile, indent=2, sort_keys=True), encoding='utf-8')
    (out/'PROOF.json').write_text(json.dumps(result, indent=2, sort_keys=True), encoding='utf-8')
    print(json.dumps({'profile': str(out/'PROFILE.json'), 'proof': str(out/'PROOF.json'),
                      'feedback_g_uS': paired[0]['actual_source_parameters']['derived_feedback_g_uS'],
                      'unique_operators': len(paired), 'dense_sparse_error': max_error}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('condition');parser.add_argument('out')
    args = parser.parse_args();run(args.condition, args.out)
