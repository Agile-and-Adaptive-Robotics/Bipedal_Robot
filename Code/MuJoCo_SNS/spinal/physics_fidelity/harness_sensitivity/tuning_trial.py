"""Evaluate one fixed-support, bilateral neural recruitment change in myo."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from replay import HERE, SPINAL, digest, require_digest, run_arm, EXPECTED


# Source pool, target muscle(s). One multiplier applies on both sides.
OUTPUTS = {
    'knee_flexor': ('PF_KNEE-F', ('bifemsh',)),
    'rectus_femoris': ('PF_BI_JOINT_RECT', ('rect_fem',)),
    'hip_flexors': ('PF_HIP-F', ('iliacus', 'psoas')),
    'dorsiflexor': ('PF_FUNC_ANKLE_Z_N', ('tib_ant',)),
}
BOUNDS = {'knee_flexor': (0.7, 2.4), 'rectus_femoris': (0.45, 1.6),
          'hip_flexors': (0.85, 2.15), 'dorsiflexor': (0.80, 1.36)}


def adjusted_profile(profile, values, directory):
    if set(values) != set(BOUNDS):
        raise ValueError('Expected all four recruitment multipliers')
    checkpoint = profile['parameter_checkpoint']
    source = Path(checkpoint['path'])
    require_digest(source, checkpoint['sha256'])
    with np.load(source) as data:
        arrays = {key: data[key] for key in data.files}
    original = arrays['ff_gain'].copy()
    names = arrays['names']
    sources, targets = names[arrays['pre']], names[arrays['post']]
    edges = {}
    for group, (pool, muscles) in OUTPUTS.items():
        factor = values[group]
        low, high = BOUNDS[group]
        if not low <= factor <= high:
            raise ValueError(f'{group} outside [{low}, {high}]')
        for muscle in muscles:
            baseline_pair = []
            for side in ('r', 'l'):
                pre, post = f'{pool}_{side}', f'MN_{muscle}_{side}'
                indices = np.flatnonzero((sources == pre) & (targets == post))
                if len(indices) != 1 or original[indices[0]] <= 0:
                    raise ValueError(f'Expected one active edge {pre} to {post}')
                i = int(indices[0])
                new = float(original[i] * factor)
                arrays['ff_gain'][i] = new
                edges[str(i)] = dict(source=pre, target=post, before=float(original[i]), after=new)
                baseline_pair.append(original[i])
            if not np.isclose(baseline_pair[0], baseline_pair[1], rtol=1e-10, atol=1e-12):
                raise ValueError(f'Original paired muscle coefficients differ for {muscle}')
    actual_changes = set(np.flatnonzero(arrays['ff_gain'] != original))
    if not actual_changes.issubset({int(i) for i in edges}):
        raise ValueError('An unintended synapse changed')
    path = directory / 'parameters.npz'
    np.savez_compressed(path, **arrays)
    return {**profile, 'parameter_checkpoint': {'path': str(path), 'sha256': digest(path)}}, edges


def score_trial(folder, frozen):
    for subdirectory in ('', 'fsa_design_20261101', 'hunt_training', 'functional_gait2392'):
        sys.path.insert(0, str(frozen / subdirectory))
    # Score with the exact physical gate and optimization objective used in the
    # original campaign, including both fall fields and the 5–13 s window.
    from ground_motor_refinement_v2_adapter import measure
    result = json.loads((folder / 'repeat.result.json').read_bytes())
    objective = measure(result, folder / 'repeat.assistance.npz')
    gates = {key: value for key, value in result['assessment'].items()
             if key.endswith('_gate') or key in ('finite', 'complete', 'accepted', 'ground_load_verified')}
    return dict(objective=objective, gates=gates,
                fall=bool(result['metrics'].get('fall', False)),
                balance_fall=bool(result['metrics'].get('bal_fell', False)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--params', help='JSON object containing the four bilateral multipliers.')
    choice.add_argument('--recorded-recovery', choices=('first', 'matched-support'),
                        nargs='?', const='first', help='Use a recorded five-right, four-left solution.')
    parser.add_argument('--arm', choices=('selective', 'native'), default='selective')
    parser.add_argument('--keep-trace', action='store_true')
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    frozen = args.replay.resolve() / 'snapshot/repo/Code/MuJoCo_SNS/spinal'
    evaluator = frozen / 'hunt_training/physical_eval.py'
    require_digest(evaluator, EXPECTED['sources']['physical_eval_sha256'])
    require_digest(HERE.parent / 'thelen_dynamic.py', EXPECTED['sources']['thelen_dynamic_sha256'])
    source_model = SPINAL.parents[2] / 'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim'
    require_digest(source_model, EXPECTED['sources']['source_osim_sha256'])
    reference = args.reference.resolve()
    require_digest(reference / 'capture/aligned_trace.npz', EXPECTED['native']['trace_sha256'])
    profile = json.loads((args.replay / 'native/bound_input_profile.json').read_bytes())
    if args.recorded_recovery:
        reference_result = json.loads((HERE / 'tuning_result.json').read_bytes())
        key = 'first_recovery' if args.recorded_recovery == 'first' else 'adaptive_matched_support_recovery'
        parameters = reference_result[key]['group_multipliers']
    else:
        parameters = json.loads(args.params)
    tuned, changes = adjusted_profile(profile, parameters, out)
    config = json.loads((reference / 'repeat.json').read_bytes())
    # On Windows a log viewer can briefly hold the handle after the result is
    # complete. The compact result survives a deferred temporary-file cleanup.
    with tempfile.TemporaryDirectory(dir=out, ignore_cleanup_errors=True) as temporary:
        temp = Path(temporary).resolve()
        if temp.parent != out:
            raise ValueError('Unexpected temporary trial folder')
        folder = (out if args.keep_trace else temp) / 'forward'
        try:
            summary = run_arm(args.arm, folder, config, tuned, evaluator, corrected_start=True)
            result = score_trial(folder, frozen)
        except subprocess.CalledProcessError:
            log_tail = (folder / 'run.log').read_bytes()[-5000:].decode('utf-8', errors='replace')
            (out / 'failed_forward.log').write_text(log_tail, encoding='utf-8')
            reasons = ('Fiber equilibrium or length limit at', 'Invalid muscle or body state at')
            reason = next((line for line in log_tail.splitlines()
                           if line.startswith('RuntimeError:') and any(s in line for s in reasons)), None)
            if reason is None:
                raise
            # A partial trajectory has no clearance count or full-window assistance.
            summary = dict(complete=False, trace_sha256=None, right_left_clearance_episodes=None,
                           right_left_measured_contact_cycles=None, mean_net_superior_percent_BW=None)
            result = dict(objective=dict(value=1e6, strict_accepted=False, unsafe=True, reason=reason),
                          gates={'complete': False}, fall=None, balance_fall=None)
        result.update(arm=args.arm, parameters=parameters, checkpoint=tuned['parameter_checkpoint']['sha256'],
                      edges=changes, trace=summary,
                      projection_sha256=digest(SPINAL / 'joint_projection.py'))
        (out / 'trial.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        print(json.dumps({'value': result['objective']['value'], 'summary': summary}))


if __name__ == '__main__':
    main()
