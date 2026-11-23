"""Explore existing premotor gains and phase-specific motor-neuron weights on the same three objectives."""
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import argparse
from pathlib import Path

from pareto_feedback import BOUNDS, HERE, one_job, optuna, read


def search(output, anchors, budget, jobs, seed):
    sampler = optuna.samplers.NSGAIISampler(
        seed=seed, population_size=16, constraints_func=lambda t: t.user_attrs['constraints']
    )
    study = optuna.create_study(
        directions=['minimize'] * 3,
        sampler=sampler,
        storage='sqlite:///' + (output / 'optimizer.sqlite').as_posix(),
        study_name='premotor_and_existing_synapses',
        load_if_exists=True,
    )
    gains = {
        'swing_suppression': 'pm_gain',
        'swing_recruitment': 'pm_add',
        'swing_sensory_gain': 'pm_aff',
        'hip_abductors': 'pm_ws',
    }
    synapses = {
        'late_dorsiflexor': 'F2/ankle_df',
        'swing_hip_flexor': 'F1/hip_flex',
        'early_knee_extensor': 'E1/knee_ext',
        'late_knee_extensor': 'E2/knee_ext',
    }
    bounds = {
        'swing_suppression': (0.20, 2.7),
        'swing_recruitment': (0.10, 2.5),
        'swing_sensory_gain': (0.35, 2.5),
        'hip_abductors': (0.18, 2.8),
        'late_dorsiflexor': (0.12, 3.0),
        'swing_hip_flexor': (0.32, 2.5),
        'early_knee_extensor': (0.30, 3.0),
        'late_knee_extensor': (0.30, 3.0),
        'plantar_flexion': (0.34, 3.0),
        'premotor_period': (0.68, 1.17),
    }
    base = {k: 1.0 for k in bounds}
    if not study.trials:
        for mix in (0.0, 1.0, 0.55):
            for tuning in (
                {},
                {'swing_recruitment': 0.35, 'swing_suppression': 1.5},
                {'swing_recruitment': 1.8, 'swing_suppression': 2.0},
                {'late_dorsiflexor': 0.30, 'plantar_flexion': 2.0, 'early_knee_extensor': 1.6},
                {
                    'swing_recruitment': 0.35,
                    'late_dorsiflexor': 0.3,
                    'plantar_flexion': 2.0,
                    'swing_hip_flexor': 1.7,
                },
            ):
                study.enqueue_trial({**base, 'mix': mix, **tuning})
    if any(x.state == optuna.trial.TrialState.RUNNING for x in study.trials):
        raise RuntimeError('Complete or explicitly recover any interrupted trials before resuming')
    count = sum(x.state == optuna.trial.TrialState.COMPLETE for x in study.trials)
    running = {}
    with ProcessPoolExecutor(max_workers=jobs) as workers:
        while count < budget:
            while len(running) < min(jobs, budget - count):
                trial = study.ask()
                mix = trial.suggest_float('mix', 0.0, 1.0)
                fitted = {k: trial.suggest_float(k, *interval, log=True) for k, interval in bounds.items()}
                spec = {
                    **{k: 1.0 for k in BOUNDS},
                    'anchor': 'stance',
                    'mix': mix,
                    'plantar_flexion': fitted['plantar_flexion'],
                    'premotor_period': fitted['premotor_period'],
                    'premotor_factors': {label: fitted[k] for k, label in gains.items()},
                    'pattern_weights': {label: fitted[k] for k, label in synapses.items()},
                }
                running[
                    workers.submit(
                        one_job, (spec, anchors, output / f'trial_{trial.number:04d}', False, None, False)
                    )
                ] = trial
            ready, _ = wait(running, return_when=FIRST_COMPLETED)
            for task in ready:
                trial = running.pop(task)
                try:
                    row = task.result()
                except Exception:
                    study.tell(trial, state=optuna.trial.TrialState.FAIL)
                    raise
                trial.set_user_attr('constraints', row['constraints'])
                values = [v if v is not None else 180.0 for v in row['objectives']]
                study.tell(trial, values)
                print(trial.number, values, row['cycles'], row['meets_3pp'], flush=True)
                count += 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--anchors', default=HERE / 'pareto_reference/anchors.json', type=Path)
    p.add_argument('--budget', type=int, default=120)
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--seed', type=int, default=94862)
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    search(out, read(args.anchors.resolve()), args.budget, args.jobs, args.seed)


if __name__ == '__main__':
    main()
