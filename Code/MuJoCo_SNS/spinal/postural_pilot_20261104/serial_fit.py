"""Bounded Rustuna tuning and fresh holdouts for each coupled pendulum."""
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'virtual_walker_20261031'))
sys.path.insert(0, r'D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311')
from rustuna_transfer import Worker
import rustuna


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--single-fit', type=Path, required=True)
    parser.add_argument('--links', type=int, choices=(2, 3), required=True)
    parser.add_argument('--axes', type=int, choices=(1, 2), default=1)
    parser.add_argument('--seed-fit', type=Path)
    parser.add_argument('--model-predictor', action='store_true')
    parser.add_argument('--disturbance-observer', action='store_true')
    parser.add_argument('--joint-shaping', action='store_true')
    parser.add_argument('--trials', type=int, default=32)
    parser.add_argument('--budget-s', type=float, default=720.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.trials <= 128 or not 45 <= args.budget_s <= 2400:
        raise ValueError('Use a bounded serial fit')
    if args.disturbance_observer and not args.model_predictor:
        raise ValueError('The disturbance observer requires the model predictor')
    if args.joint_shaping and not args.disturbance_observer:
        raise ValueError('Joint shaping requires the disturbance-observer fit')
    args.output.mkdir(parents=True, exist_ok=False)
    preflight = dict(op='preflight', single_fit=str(args.single_fit), links=args.links, axes=args.axes,
                     model_predictor=args.model_predictor, disturbance_observer=args.disturbance_observer,
                     joint_shaping=args.joint_shaping)
    (args.output / 'MANIFEST.json').write_text(json.dumps(
        dict(preflight, seed=61106, trials=args.trials, budget_s=args.budget_s,
             scope='Fixed-base serial pendulum; frozen antagonistic SNS modules; tune lag-aware task gains and kinematic prediction'),
        indent=2))
    study = rustuna.create_study(direction='minimize', sampler=rustuna.samplers.TPESampler(
        seed=61106, n_startup_trials=8))
    if args.seed_fit:
        prior = json.loads(args.seed_fit.read_text())
        if 'best' in prior:
            prior = prior['best']['parameters']
        for prediction in (0., .75, 1.25):
            candidate=dict(effort=prior['effort'],
                           neural_tau_s=prior.get('neural_tau_s', prior.get('modeled_neural_tau_s')),
                           delay_prediction=prediction)
            if args.disturbance_observer:
                candidate['disturbance_gain'] = prior.get('disturbance_gain', .65)
            if args.joint_shaping:
                candidate.update(distal_effort=1., velocity_weight=1.)
            study.enqueue_trial(candidate)
    deadline = time.monotonic() + args.budget_s
    worker = Worker(HERE / 'serial_fit_adapter.py', args.output)
    best = None
    try:
        worker.request(preflight, 45.)
        for number in range(args.trials):
            if time.monotonic() + 25. >= deadline:
                break
            trial = study.ask()
            effort_bounds = (.0003, .1) if args.joint_shaping else (.003, 3.)
            params = dict(effort=trial.suggest_float('effort', *effort_bounds, log=True),
                          neural_tau_s=trial.suggest_float('neural_tau_s', .005, .05),
                          delay_prediction=trial.suggest_float('delay_prediction', 0., 1.5))
            if args.disturbance_observer:
                params['disturbance_gain'] = trial.suggest_float('disturbance_gain', 0., 1.)
            if args.joint_shaping:
                params['distal_effort'] = trial.suggest_float('distal_effort', .2, 3., log=True)
                params['velocity_weight'] = trial.suggest_float('velocity_weight', .45, 3., log=True)
            result = worker.request(dict(op='evaluate', params=params), 60.)
            study.tell(trial.number, result['value'])
            row = dict(number=number, parameters=params, **result)
            with (args.output / 'TRIALS.jsonl').open('a', encoding='utf-8') as stream:
                stream.write(json.dumps(row) + '\n')
            if best is None or result['value'] < best['value']:
                best = row
                (args.output / 'BEST.json').write_text(json.dumps(best, indent=2))
    finally:
        worker.close()
    if best is None:
        raise RuntimeError('No completed fit trials')
    fresh = Worker(HERE / 'serial_fit_adapter.py', args.output)
    try:
        fresh.request(preflight, 45.)
        repeat = fresh.request(dict(op='evaluate', params=best['parameters']), 60.)
        unseen = fresh.request(dict(op='held_out', params=best['parameters']), 90.)
    finally:
        fresh.close()
    final = dict(best=best, repeat=repeat, holdout=unseen,
                 fresh_process_matched=abs(repeat['value'] - best['value']) < 1e-10)
    final['accepted'] = bool(final['fresh_process_matched'] and best['passed'] == best['total']
                             and unseen['passed'] == unseen['total'])
    (args.output / 'FINAL.json').write_text(json.dumps(final, indent=2))


if __name__ == '__main__':
    main()
