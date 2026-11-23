"""Fit contact-reset swing timing and existing hip/knee motor projections to single-counted foot loading and joint angles."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

from pareto_feedback import BOUNDS, HERE, one_job, optuna, read


def candidate(d):
    return {
        **{k: 1.0 for k in BOUNDS},
        'anchor': d['anchor'],
        'mix': d['mix'],
        'right_drive': d['right_drive'],
        'left_drive': d['left_drive'],
        'force_convention': 'physical_50_N',
        'swing_window': [d['onset'], min(d['offset'], 1.0 - d['edge']), d['edge']],
        'pattern_weights': {
            'E2/hip_ext': d['hip_stance_ext'],
            'F1/knee_flex': d['swing_knee'],
            'F2/ankle_df': d['late_dorsi'],
        },
        'premotor_factors': {'pm_add': d['swing_drive']},
        'regional_synapses': {'hip_flex/both/posture_to_mn': d['hip_flex_tone']},
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--anchors', type=Path, default=HERE / 'pareto_reference/anchors.json')
    cli.add_argument('--jobs', type=int, default=4)
    cli.add_argument('--budget', type=int, default=112)
    cli.add_argument('--seed', type=int, default=138714)
    args = cli.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    anchors = read(args.anchors.resolve())
    sampler = optuna.samplers.TPESampler(
        seed=args.seed,
        multivariate=True,
        n_startup_trials=25,
        constraints_func=lambda t: t.user_attrs['constraints'],
    )
    study = optuna.create_study(
        directions=['minimize'] * 3,
        study_name='physical_swing_and_angle',
        sampler=sampler,
        storage='sqlite:///' + (out / 'optimizer.sqlite').as_posix(),
        load_if_exists=True,
    )
    ranges = {
        'onset': (0.46, 0.74),
        'offset': (0.79, 0.966),
        'edge': (0.025, 0.09),
        'right_drive': (0.55, 1.55),
        'left_drive': (0.55, 1.55),
        'hip_stance_ext': (0.35, 7.0),
        'swing_knee': (0.40, 5.0),
        'late_dorsi': (0.16, 2.3),
        'swing_drive': (0.2, 2.25),
        'hip_flex_tone': (0.12, 1.6),
    }
    if not study.trials:
        default = {k: 1.0 for k in ranges}
        default.update(onset=0.62, offset=0.95, edge=0.05)
        for adjustment in (
            {},
            {'onset': 0.54, 'offset': 0.91},
            {'onset': 0.50, 'offset': 0.88},
            {'hip_stance_ext': 2.4, 'swing_knee': 1.85, 'hip_flex_tone': 0.55},
            {'hip_stance_ext': 4.0, 'swing_knee': 2.8, 'hip_flex_tone': 0.30, 'late_dorsi': 0.48},
            {'onset': 0.56, 'offset': 0.91, 'hip_stance_ext': 3.0, 'swing_knee': 2.3, 'hip_flex_tone': 0.45},
            {'onset': 0.68, 'offset': 0.965, 'edge': 0.032, 'hip_stance_ext': 2.5, 'swing_knee': 2.5},
        ):
            for anchor, mix in (('stance', 0.0), ('stance', 0.30), ('r6t10', 1.0), ('r6t10', 0.50)):
                study.enqueue_trial({**default, **adjustment, 'anchor': anchor, 'mix': mix})
    done_count = sum(t.state == optuna.trial.TrialState.COMPLETE for t in study.trials)
    running = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as workers:
        while done_count < args.budget:
            while len(running) < min(args.jobs, args.budget - done_count):
                t = study.ask()
                p = {
                    'anchor': t.suggest_categorical('anchor', ['stance', 'r6t10', 'symmetric']),
                    'mix': t.suggest_float('mix', 0.0, 1.0),
                }
                for k, bounds in ranges.items():
                    p[k] = t.suggest_float(k, *bounds, log=k not in ('onset', 'offset', 'edge'))
                running[
                    workers.submit(
                        one_job, (candidate(p), anchors, out / f'trial_{t.number:04d}', False, None, False)
                    )
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
                v = [k if k is not None else 180.0 for k in r['objectives']]
                study.tell(t, v)
                print(f"{t.number} {v} {r['cycles']} {r['meets_3pp']}", flush=True)
                done_count += 1


if __name__ == '__main__':
    main()
