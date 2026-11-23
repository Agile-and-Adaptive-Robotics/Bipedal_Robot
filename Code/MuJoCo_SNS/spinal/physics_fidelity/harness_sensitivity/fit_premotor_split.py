"""Refine measurable Pareto parents by separating the existing swing recruitment of antagonist muscle groups."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

from pareto_feedback import HERE, np, one_job, optuna, read


def make_candidate(base, parameters):
    p = {**base, 'force_convention': 'physical_50_N'}
    p['premotor_group_gains'] = {
        'knee_flex': parameters['swing_knee'],
        'ankle_df': parameters['swing_ankle'],
        'hip_flex': parameters['swing_hip'],
        'ankle_pf': parameters['plantar_release'],
    }
    p['right_drive'] *= parameters['right_drive']
    p['left_drive'] *= parameters['left_drive']
    p['pattern_weights'] = dict(base.get('pattern_weights', {}))
    p['pattern_weights']['E2/hip_ext'] = (
        p['pattern_weights'].get('E2/hip_ext', 1.0) * parameters['stance_hip_ext']
    )
    p['regional_synapses'] = dict(base.get('regional_synapses', {}))
    for group, gain in [
        ('hip_abd', parameters['right_abductors']),
        ('hip_add', 1.0 / parameters['right_abductors']),
    ]:
        key = group + '/r/ii_to_mn'
        p['regional_synapses'][key] = p['regional_synapses'].get(key, 1.0) * gain
    return p


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--parents', type=Path, required=True)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--anchors', type=Path, default=HERE / 'pareto_reference/anchors.json')
    cli.add_argument('--budget', type=int, default=96)
    cli.add_argument('--jobs', type=int, default=6)
    cli.add_argument('--concurrency', type=Path)
    cli.add_argument('--seed', type=int, default=238411)
    args = cli.parse_args()
    parents = read(args.parents.resolve())
    anchors = read(args.anchors.resolve())
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sampler = optuna.samplers.TPESampler(
        seed=args.seed,
        multivariate=True,
        n_startup_trials=25,
        constraints_func=lambda t: t.user_attrs['constraints'],
    )
    study = optuna.create_study(
        study_name='split_existing_premotor_recruitment',
        directions=['minimize'] * 3,
        sampler=sampler,
        storage='sqlite:///' + (out / 'optimizer.sqlite').as_posix(),
        load_if_exists=True,
    )
    bounds = {
        'swing_knee': (0.36, 3.5),
        'swing_ankle': (0.07, 2.7),
        'swing_hip': (0.22, 2.3),
        'plantar_release': (0.24, 2.45),
        'stance_hip_ext': (0.5, 5.0),
        'right_abductors': (0.40, 2.7),
        'right_drive': (0.71, 1.34),
        'left_drive': (0.71, 1.34),
    }
    if not study.trials:
        base = {k: 1.0 for k in bounds}
        for adjustment in (
            {},
            {'swing_knee': 1.6, 'swing_ankle': 0.45, 'swing_hip': 0.75, 'stance_hip_ext': 1.75},
            {'swing_knee': 2.5, 'swing_ankle': 0.22, 'swing_hip': 0.50, 'stance_hip_ext': 2.6},
            {'right_abductors': 1.38, 'right_drive': 1.07},
            {'right_abductors': 0.75, 'right_drive': 0.92},
            {'swing_knee': 1.22, 'swing_ankle': 0.75, 'swing_hip': 0.90, 'right_abductors': 1.15},
        ):
            for name in parents:
                study.enqueue_trial({**base, **adjustment, 'parent': name})
    count = sum(t.state == optuna.trial.TrialState.COMPLETE for t in study.trials)
    running = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        while count < args.budget:
            available = (
                min(args.jobs, int(args.concurrency.read_text().strip())) if args.concurrency else args.jobs
            )
            while len(running) < min(available, args.budget - count):
                t = study.ask()
                label = t.suggest_categorical('parent', sorted(parents))
                fit = {k: t.suggest_float(k, *v, log=True) for k, v in bounds.items()}
                p = make_candidate(parents[label], fit)
                running[
                    pool.submit(one_job, (p, anchors, out / f'trial_{t.number:04d}', False, None, False))
                ] = t
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for f in done:
                t = running.pop(f)
                try:
                    r = f.result()
                except Exception:
                    study.tell(t, state=optuna.trial.TrialState.FAIL)
                    raise
                t.set_user_attr('constraints', r['constraints'])
                values = [v if v is not None else 180.0 for v in r['objectives']]
                study.tell(t, values)
                print(t.number, values, r['cycles'], r['meets_3pp'], flush=True)
                count += 1


if __name__ == '__main__':
    main()
