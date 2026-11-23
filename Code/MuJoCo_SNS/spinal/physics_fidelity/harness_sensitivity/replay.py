"""Reproduce the human Gait2392 pelvic-harness muscle-law comparison.

Run with the spinal stack's MuJoCo 2.3.7 environment on easteregg2. The
archived runtime and its side-by-side SOURCE_MANIFEST.json, plus the
saved experiment directory containing repeat.json and capture/, are inputs.
The submitted config also references its original learned weights, profiles
and mechanical assets. These must be accessible at the recorded paths.

python replay.py --bundle <runtime.zip> --reference <other_models/gait_hunt_rc> --out <scratch-directory>

The two evaluations run in separate processes. The NPZ/MJB recordings and
logs stay in the output directory. The compact results go to result.json.
In the selective run the Python adapter advances the four fiber states and
updates applied tensions at every physics step.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

import numpy as np


HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent.parent
EXPECTED = json.loads((HERE / 'reference_result.json').read_text(encoding='utf-8'))


def digest(path):
    with path.open('rb') as stream:
        checksum = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def require_digest(path, expected):
    actual = digest(path)
    if actual != expected:
        raise ValueError(f'Source hash mismatch at {path}: {actual} != {expected}')


def restore_snapshot(bundle, target):
    hashes = EXPECTED['sources']
    require_digest(bundle, hashes['runtime_zip_sha256'])
    manifest = json.loads((bundle.parent / 'SOURCE_MANIFEST.json').read_bytes())
    if len(manifest['files']) != hashes['validated_source_entries']:
        raise ValueError('Source manifest differs from the selected runtime')
    sources = {}
    target.mkdir()
    with ZipFile(bundle) as archive:
        for entry in manifest['files']:
            name = entry['zip_relative_path']
            dest = (target / name).resolve()
            if not dest.is_relative_to(target):
                raise ValueError(f'Unexpected member path {name}')
            data = archive.read(name)
            if hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise ValueError(f'Archived source hash mismatch: {name}')
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            sources[Path(entry['original_path'])] = dest

    frozen_spinal = target / 'repo/Code/MuJoCo_SNS/spinal'
    gate = Path('fsa_design_20261101/CAMPAIGN_GATE.json')
    require_digest(SPINAL / gate, hashes['campaign_gate_sha256'])
    shutil.copyfile(SPINAL / gate, frozen_spinal / gate)
    return frozen_spinal, sources


def bind_profile(config, sources):
    profile = json.loads(Path(config['integrated_fsa_profile']).read_bytes())
    checked_sources = {}
    for name, expected in profile['source_sha256'].items():
        original = Path(name)
        path = sources.get(original, original)
        require_digest(path, expected)
        checked_sources[str(path)] = expected
    checkpoint = profile['parameter_checkpoint']
    require_digest(Path(checkpoint['path']), checkpoint['sha256'])
    return {**profile, 'source_sha256': checked_sources}


def summarize(directory):
    metrics = json.loads((directory / 'repeat.result.json').read_bytes())
    capture = directory / 'capture'
    meta = json.loads((capture / 'ALIGNED_TRACE_METADATA.json').read_bytes())
    trace = capture / 'aligned_trace.npz'
    with np.load(trace) as data:
        time = data['force_interval_end']
        selected = (time >= 5.0 - 1e-10) & (time <= 13.0 + 1e-10)
        harness = data['harness_wrench_anatomical'][selected]
        explicit = data['explicit_support_anatomical_N'][selected, 1]
        bw = meta['body_weight_N']
        support = np.mean(harness[:, 1] + explicit) * 100 / bw
        horizontal = np.sqrt(np.mean(harness[:, 0]**2 + harness[:, 2]**2))
    feet = metrics['foot_clearance']['feet']
    return {
        'coupled_steps': len(time),
        'scored_force_intervals': int(np.count_nonzero(selected)),
        'body_weight_N': bw,
        'trace_sha256': digest(trace),
        'finite': metrics['assessment']['finite'],
        'complete': metrics['assessment']['complete'],
        'right_left_measured_contact_cycles': metrics['ground_contacts']['cycles'],
        'right_left_knee_cycles': metrics['assessment']['knee_cycles'],
        'right_left_clearance_episodes': [foot['sustained_swings'] for foot in feet],
        'right_left_peak_clearance_mm': [1000 * foot['maximum_clearance_m'] for foot in feet],
        'mean_net_superior_percent_BW': float(support),
        'prescribed_superior_N': float(np.mean(explicit)),
        'horizontal_restraint_rms_N': float(horizontal),
    }


def run_arm(name, directory, config, profile, evaluator, corrected_start=False):
    directory.mkdir()
    prof_file = directory / 'bound_input_profile.json'
    prof_file.write_text(json.dumps(profile), encoding='utf-8')
    request = {**config, 'verification_capture': str(directory / 'capture'),
               'integrated_fsa_profile': str(prof_file)}
    run_file = directory / 'repeat.json'
    run_file.write_text(json.dumps(request), encoding='utf-8')
    env = {key: value for key, value in os.environ.items() if not key.startswith('AARL_')}
    env.update(CONDA_PREFIX=sys.prefix, MPLBACKEND='Agg', PYTHONIOENCODING='utf-8',
               PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    fallback = ('', 'hunt_training', 'functional_gait2392', 'fsa_design_20261101', 'w2l_cpg')
    env['PYTHONPATH'] = os.pathsep.join(str(SPINAL / path) for path in fallback)
    cmd = [sys.executable, str(HERE / 'selective_force.py'), name,
           str(evaluator), str(run_file)]
    if corrected_start:
        cmd.append('--corrected-start')
    with (directory / 'run.log').open('wb') as log:
        subprocess.run(cmd, cwd=evaluator.parent.parent, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    return summarize(directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--families', action='store_true',
                        help='Also test the two individual bilateral muscle groups with corrected startup.')
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error('--out must be new or empty. Full traces are written there.')
    out.mkdir(parents=True, exist_ok=True)
    reference = args.reference.resolve()
    config = json.loads((reference / 'repeat.json').read_bytes())
    require_digest(reference / 'capture/aligned_trace.npz', EXPECTED['native']['trace_sha256'])
    require_digest(SPINAL.parent.parent.parent /
                   'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim',
                   EXPECTED['sources']['source_osim_sha256'])
    require_digest(HERE.parent / 'thelen_dynamic.py', EXPECTED['sources']['thelen_dynamic_sha256'])
    frozen_spinal, sources = restore_snapshot(args.bundle, out / 'snapshot')
    profile = bind_profile(config, sources)
    evaluator = frozen_spinal / 'hunt_training/physical_eval.py'
    require_digest(evaluator, EXPECTED['sources']['physical_eval_sha256'])
    results = {}
    arms = ('native', 'rect_fem', 'bifemsh', 'selective') if args.families else ('native', 'selective')
    for name in arms:
        results[name] = run_arm(name, out / name, config, profile, evaluator,
                                corrected_start=args.families)
        if name in EXPECTED and results[name]['trace_sha256'] != EXPECTED[name]['trace_sha256']:
            raise RuntimeError(f'{name} trajectory differs from the recorded trajectory. '
                               f'Inspect {out / name}.')
    native_start = json.loads((out / 'native/controlled_start.json').read_bytes())
    other_start = json.loads((out / 'selective/controlled_start.json').read_bytes())
    if native_start != other_start:
        raise RuntimeError('Arms did not start from identical model, state, input and harness.')
    if args.families:
        for name in ('rect_fem', 'bifemsh'):
            start = json.loads((out / name / 'controlled_start.json').read_bytes())
            for field in ('qpos', 'qvel', 'act', 'ctrl', 'xfrc_applied'):
                if start[field] != native_start[field]:
                    raise RuntimeError(f'{name} started with a different {field}.')
    (out / 'result.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
