"""Repeat an archived Gait2392 air or supported-contact condition.

Pass a ControllerTraining study directory and a condition relative to it.
Keep the original circuit or select a new fitted D directory with --checkpoint.
The output records the original settings, applied neural and physics sources,
and muscle forces. An existing completed experiment cannot be overwritten.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys
from unittest.mock import patch

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'

from coupled_muscle_context import SPINAL, use_four_source_muscles


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('condition', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--arm', choices=('native', 'four'), default='four')
    parser.add_argument('--checkpoint', type=Path)
    args = parser.parse_args()
    study, output = args.study.resolve(), args.output.resolve()
    original = (study / args.condition).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Choose an empty output directory')
    output.mkdir(parents=True, exist_ok=True)
    source = read(original / 'SOURCE.json')
    config = read(original / 'body.json')
    old_fit = Path(source['reflex_fit'])
    fit = args.checkpoint.resolve() if args.checkpoint else old_fit
    old_gain = read(old_fit / 'FINAL.json')['parameters'].get('local_II', 1.)
    gain = read(fit / 'FINAL.json')['parameters'].get('local_II', 1.)
    if fit != old_fit:
        if not old_gain:
            raise ValueError('Cannot transfer a relative sensory gain from a zero baseline')
        config['local_II_feedback_scale'] = config.get('local_II_feedback_scale', old_gain) * gain / old_gain
    body = output / 'body'
    config['verification_capture'] = str(body / 'capture')
    family = config.get('controller_family') or read(Path(source['peak_fit']) / 'MANIFEST.json')['configuration']['family']
    profiles = SPINAL / 'controller_profiles_20261104'
    history = study / 'source_snapshot_final'
    adapter = history / 'controller_profiles_20261104/body_probe.py' if family == 'hunt' else profiles / 'body_probe.py'
    physical = history / 'existing_dependencies/hunt_training__physical_eval.py'
    metadata = original / 'capture/ALIGNED_TRACE_METADATA.json'
    expected_physics = read(metadata)['provenance']['source_hashes']['physical_eval.py'] if metadata.exists() else 'ecda3d37012e593b54ece1166cb9df440c9f961acbc172689db258535bf6ac33'
    if digest(physical) != expected_physics:
        raise ValueError('The physical evaluator differs from the source used by this archived condition')
    arguments = ['--reflex-fit', str(fit), '--stage', str(config['stage']), '--expression',
                 str(config.get('motor_scale', source['air_expression_scale'])), '--output', str(body)]
    flags = [('central_RG_feedback_scale', '--feedback-scale', 1.),
             ('local_II_feedback_scale', '--local-ii-scale', gain),
             ('ankle_crossing_local_Ia_scale', '--ankle-ia-scale', 1.),
             ('ankle_crossing_local_II_scale', '--ankle-ii-scale', 1.),
             ('hip_PF_to_MN_expression', '--hip-expression', 1.)]
    for key, flag, default in flags:
        value = config.get(key, default)
        if value != default:
            arguments.extend((flag, str(value / gain if key == 'local_II_feedback_scale' else value)))
    if config.get('preserve_knee_crossing_afferents'):
        arguments.append('--preserve-knee-crossing')
    noise = config.get('neural_noise')
    if noise:
        arguments.extend(('--noise-rg-sigma', str(noise['stationary_sigma_nA']), '--noise-seed', str(noise['seed'])))
    stimulus = config.get('neural_stimulus')
    if stimulus:
        arguments.extend(('--stim-cell', stimulus['cell'], '--stim-nA', str(stimulus['current_nA']),
                          '--stim-start', str(stimulus['start_s']), '--stim-duration', str(stimulus['duration_s'])))

    # The historical module resolves dependencies relative to its original location.
    physical_path = SPINAL / 'hunt_training/physical_eval.py'
    original_run = runpy.run_path

    def load_physics(path, init_globals=None, run_name=None):
        if Path(path).resolve() != physical_path:
            return original_run(path, init_globals=init_globals, run_name=run_name)
        config_path = Path(sys.argv[1])
        generated = read(config_path)
        for key, _, _ in flags:
            if key in generated and key in config:
                if abs(generated[key] - config[key]) > 1e-12:
                    raise ValueError(f'Mismatch in {key}: {generated[key]} vs {config[key]}')
        config_path.write_text(json.dumps(config, indent=2), encoding='utf-8')
        scope = {**(init_globals or {}), '__name__': run_name or '<run_path>', '__file__': str(physical_path)}
        exec(compile(physical.read_bytes(), str(physical_path), 'exec'), scope)
        return scope

    for key in list(os.environ):
        if key.startswith('AARL_'):
            os.environ.pop(key)
    if config.get('native_pf_spec'):
        os.environ['AARL_SYN6_NATIVE_PF'] = config['native_pf_spec']
    os.environ['AARL_AIR_START_KNEE_DEG'] = str(config['start_knee_deg'])
    os.environ.update(CONDA_PREFIX=str(Path(sys.executable).parent), MPLBACKEND='Agg')
    sys.path[:0] = [str(profiles), str(physical_path.parent), str(SPINAL),
                    str(SPINAL / 'functional_gait2392'), str(SPINAL / 'w2l_cpg')]
    os.chdir(SPINAL)
    adapter_location = str(profiles / 'body_probe.py')
    meta = {'original_condition': str(original), 'original_fit': str(old_fit),
            'applied_fit': str(fit), 'fitted_parameters_sha256': digest(fit / 'FINAL.json'),
            'physical_evaluator_sha256': digest(physical), 'adapter_sha256': digest(adapter),
            'original_condition_sha256': digest(original / 'body.json'), 'arm': args.arm}
    sys.argv = [adapter_location, *arguments]
    with use_four_source_muscles(output, args.arm), patch.object(runpy, 'run_path', load_physics):
        try:
            exec(compile(adapter.read_bytes(), adapter_location, 'exec'),
                 {'__name__': '__main__', '__file__': adapter_location})
        except RuntimeError:
            stop = body / 'STOP.json'
            if not stop.exists() or 'Predictive knee guard' not in read(stop).get('error', ''):
                raise
            meta['knee_guard'] = read(stop)['knee_guard']
    report = read(output / 'force_audit.json')
    if len(report['records']) != 1 or not report['records'][0]['steps'] or report['force_failure']:
        raise RuntimeError('The coupled force audit did not complete')
    meta['coupled_steps'] = report['records'][0]['steps']
    (output / 'condition.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
