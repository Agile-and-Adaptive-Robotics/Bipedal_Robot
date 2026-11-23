"""Bounded local Rustuna search for bilateral clearance with four Thelen muscles.

Run with the 3.11 Rustuna environment. Supply the myo 3.10 Python executable,
the original reference directory, and the output of replay.py. Each candidate
runs the complete closed-loop experiment in a fresh myo process. Records are
written in numbered order so each batch gives Rustuna a deterministic history.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import importlib.metadata
import os
from pathlib import Path
import subprocess
import sys

from replay import HERE, SPINAL, digest
from tuning_trial import BOUNDS

sys.path.insert(0, str(SPINAL / 'virtual_walker_20261031'))
from rustuna_runtime import load_rustuna


SEED = 823419
INITIAL = [
    (1.0, 1.0, 1.0, 1.0),
    (1.55, 0.75, 1.4, 1.12),
    (2.15, 0.55, 1.7, 1.24),
    (1.75, 1.15, 1.35, 1.0),
    (1.0, 1.0, 1.65, 1.16),
    (1.3, 0.7, 1.45, 0.91),
]


def call_trial(args, label, parameters, keep_trace=False, arm='selective'):
    output = args.out / label
    result = output / 'trial.json'
    if result.is_file():
        previous = json.loads(result.read_bytes())
        if previous['parameters'] != parameters or previous.get('arm', 'selective') != arm:
            raise ValueError(f'Incompatible saved parameters in {output}')
        return previous
    command = [str(args.myo), str(HERE / 'tuning_trial.py'), '--reference', str(args.reference),
               '--replay', str(args.replay), '--out', str(output),
               '--params', json.dumps(parameters), '--arm', arm]
    if keep_trace:
        command.append('--keep-trace')
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'OPENBLAS_NUM_THREADS': '1',
           'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1'}
    # If a stopped trial left an incomplete folder, preserve its evidence and
    # relaunch under a new directory before assigning the requested result.
    if output.exists():
        attempt = 1
        while output.with_name(f'{label}_incomplete_{attempt}').exists():
            attempt += 1
        output.rename(output.with_name(f'{label}_incomplete_{attempt}'))
    with (args.out / f'{label}.log').open('w', encoding='utf-8') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    return json.loads(result.read_bytes())


def record_progress(args, rows):
    ledger = [{'number': i, 'params': row['parameters'], 'value': row['objective']['value'],
               'strict_accepted': row['objective']['strict_accepted'], 'trace': row['trace'],
               'gates': row['gates'], 'fall': row['fall'], 'balance_fall': row['balance_fall']}
              for i, row in enumerate(rows)]
    (args.out / 'trials.json').write_text(json.dumps(ledger, indent=2), encoding='utf-8')


def record_protocol(args):
    frozen = args.replay / 'snapshot/repo/Code/MuJoCo_SNS/spinal'
    protocol = {
        'optimizer': 'Rustuna TPE', 'version': importlib.metadata.version('rustuna'),
        'seed': SEED, 'workers_and_batch_width': args.workers, 'startup_trials': len(INITIAL),
        'bounds': {name: list(pair) for name, pair in BOUNDS.items()},
        'initial_candidates': [dict(zip(BOUNDS, row)) for row in INITIAL],
        'reference_configuration_sha256': digest(args.reference / 'repeat.json'),
        'frozen_evaluator_sha256': digest(frozen / 'hunt_training/physical_eval.py'),
        'original_objective_sha256': digest(frozen / 'fsa_design_20261101/ground_motor_refinement_v2_adapter.py'),
        'projection_sha256': digest(SPINAL / 'joint_projection.py'),
        'four_actuator_law_sha256': digest(HERE.parent / 'thelen_dynamic.py'),
    }
    path = args.out / 'search_protocol.json'
    if path.is_file() and json.loads(path.read_bytes()) != protocol:
        raise ValueError('Existing output uses a different protocol or simulator configuration')
    path.write_text(json.dumps(protocol, indent=2), encoding='utf-8')


def verify_winner(args):
    ledger = json.loads((args.out / 'trials.json').read_bytes())
    best = min(ledger, key=lambda row: row['value'])
    if args.verify_trial is not None:
        best = next((row for row in ledger if row['number'] == args.verify_trial), None)
        if best is None:
            raise ValueError('Selected trial is absent from the recorded evaluations.')
    if not best['strict_accepted']:
        raise RuntimeError('No candidate passed the original physical gates.')
    targets = [('nominal', best['params'], 'selective'), ('native_transfer', best['params'], 'native')]
    for name in BOUNDS:
        for factor, direction in ((0.99, 'lower'), (1.01, 'higher')):
            changed = {**best['params'], name: best['params'][name] * factor}
            if not BOUNDS[name][0] <= changed[name] <= BOUNDS[name][1]:
                raise ValueError(f'Validation {direction} for {name} lies outside its search interval.')
            targets.append((f'{name}_{direction}', changed, 'selective'))
    args.out = args.out / f'verify_{best["number"]:02}'
    args.out.mkdir(exist_ok=True)
    checks = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [(label, pool.submit(call_trial, args, label, params,
                                       label in ('nominal', 'native_transfer'), arm))
                   for label, params, arm in targets]
        for label, future in futures:
            result = future.result()
            checks[label] = {'params': result['parameters'], 'trace': result['trace'],
                             'strict_accepted': result['objective']['strict_accepted'],
                             'value': result['objective']['value']}
            (args.out / 'checks.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
            print(label, checks[label]['strict_accepted'], result['trace']['right_left_clearance_episodes'],
                  flush=True)
    if checks['nominal']['trace']['trace_sha256'] != best['trace']['trace_sha256']:
        raise RuntimeError('Fresh nominal trajectory differs from selected trial.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('myo', 'reference', 'replay', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--trials', type=int, default=12)
    parser.add_argument('--workers', type=int, default=3)
    verification = parser.add_mutually_exclusive_group()
    verification.add_argument('--verify-winner', action='store_true',
                        help='Freshly test the recorded winner and each of its four parameters at ±1 percent.')
    verification.add_argument('--verify-trial', type=int,
                              help='Apply the same repeat, native transfer, and ±1 percent tests to a recorded trial.')
    args = parser.parse_args()
    if not 1 <= args.workers <= 8 or not 1 <= args.trials <= 32:
        parser.error('Choose 1–8 workers and 1–32 trials.')
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.verify_winner or args.verify_trial is not None:
        verify_winner(args)
        return
    backend = load_rustuna()
    record_protocol(args)
    study = backend.create_study(direction='minimize',
        sampler=backend.samplers.TPESampler(seed=SEED, n_startup_trials=len(INITIAL), multivariate=False))
    for candidate in INITIAL:
        study.enqueue_trial(dict(zip(BOUNDS, candidate)))
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        while len(rows) < args.trials:
            batch = []
            for _ in range(min(args.workers, args.trials - len(rows))):
                trial = study.ask()
                params = {name: trial.suggest_float(name, low, high)
                          for name, (low, high) in BOUNDS.items()}
                future = pool.submit(call_trial, args, f'trial_{trial.number:02}', params)
                batch.append((trial.number, future))
            for number, future in batch:
                row = future.result()
                study.tell(number, row['objective']['value'])
                rows.append(row)
                record_progress(args, rows)
                print(number, row['trace']['right_left_clearance_episodes'],
                      row['objective']['strict_accepted'], row['objective']['value'], flush=True)
    best = min(rows, key=lambda row: row['objective']['value'])
    print('Best', best['parameters'], best['objective'], flush=True)


if __name__ == '__main__':
    main()
