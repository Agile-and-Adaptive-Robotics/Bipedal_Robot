"""Fit absolute joint-angle error and the two stance fractions with four source Thelen muscles."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import hashlib
import json
import os
from pathlib import Path
import sys

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import optuna

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parents[1]
sys.path[:0] = [str(HERE), str(SPINAL)]
from reference_angle_errors import compare_angles

JOINTS = {
    'plantar_flexion': 'e2_pf',
    'dorsiflexion': 'f1_df',
    'knee_flexion': 'f1_kf',
    'knee_extension': 'post_kneext',
    'hip_extension': 'post_hipext',
    'ankle_tone': 'ankle_post_walk_trim',
    'premotor_period': 'pm_T',
}
BOUNDS = {
    'mix': (0.0, 1.0),
    'right_drive': (0.40, 1.70),
    'left_drive': (0.40, 1.70),
    'plantar_flexion': (0.35, 4.0),
    'dorsiflexion': (0.20, 2.70),
    'knee_flexion': (0.4, 2.1),
    'knee_extension': (0.55, 2.1),
    'hip_extension': (0.55, 1.9),
    'ankle_tone': (0.18, 3.0),
    'premotor_period': (0.66, 1.17),
}
OPTUNA_SEED = 240178
ANGLE_LIMIT = 17.78
FROZEN_PHYSICAL_RMSE = 17.777920236828237
INITIALIZED = False
PRODUCTION_SCHEDULE = None


@contextmanager
def regional_feedback(factors):
    """Scale established homonymous afferent or tonic projections during circuit construction."""
    import build_network
    import params as neural

    original = build_network.SpinalNetwork._wire_muscle
    checked = {}
    entries = {}
    for label, factor in factors.items():
        group, side, gain = label.split('/')
        if (
            group
            not in (
                'ankle_pf',
                'ankle_df',
                'knee_ext',
                'knee_flex',
                'hip_flex',
                'hip_ext',
                'hip_abd',
                'hip_add',
            )
            or side not in ('r', 'l', 'both')
            or gain not in ('ia_to_mn', 'ii_to_mn', 'ib_to_mn_inh', 'posture_to_mn')
            or not np.isfinite(factor)
            or factor <= 0
        ):
            raise ValueError(f'Unsupported muscle feedback factor {label}: {factor}')
        for target_side in ('r', 'l') if side == 'both' else (side,):
            entries.setdefault((group, target_side), []).append((gain, factor, label))

    def write_muscle(self, network, actuator, muscle):
        factors_here = entries.get((muscle.groups[0], muscle.side), ())
        prior = {gain: neural.G[gain] for gain, _, _ in factors_here}
        try:
            for gain, factor, label in factors_here:
                neural.G[gain] *= factor
                checked.setdefault(label, set()).add(actuator)
            return original(self, network, actuator, muscle)
        finally:
            neural.G.update(prior)

    if factors:
        build_network.SpinalNetwork._wire_muscle = write_muscle
    try:
        yield checked
    finally:
        build_network.SpinalNetwork._wire_muscle = original


@contextmanager
def swing_timing(runner, interval):
    original = runner._pm_window
    if interval:
        lo, hi, edge = interval
        if not (0 < edge < lo < hi < 1 and hi + edge <= 1.0000001):
            raise ValueError(f'Invalid raised-cosine swing interval {interval}')

        def configured_window(phase):
            return original(phase, lo, hi, edge)

        runner._pm_window = configured_window
    try:
        yield
    finally:
        runner._pm_window = original


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf8')


def fitted_config(parameters, anchors, study):
    parent = anchors[parameters['anchor']]
    base = anchors['frozen']
    fraction = parameters['mix']

    def blended(left, right):
        return (1.0 - fraction) * left + fraction * right

    knobs = {
        k: blended(base['knobs'][k] * base['scale'], parent['knobs'][k] * parent['scale'])
        * parameters['left_drive' if k.endswith('_l') else 'right_drive']
        for k in ('gamma_l', 'alpha_l', 'gamma_r', 'alpha_r')
    }
    knobs['alpha_to_mn'] = blended(base['knobs']['alpha_to_mn'], parent['knobs']['alpha_to_mn'])
    feedback = {k: blended(base['params'][k], parent['params'][k]) for k in study.FB_KEYS}
    joint = {
        k: blended(base['params'][k], parent['params'][k]) * parameters[label] for label, k in JOINTS.items()
    }
    for name, factor in parameters.get('premotor_factors', {}).items():
        if name not in ('pm_gain', 'pm_add', 'pm_aff', 'pm_ws'):
            raise ValueError(f'Unknown premotor controller parameter {name}')
        joint[name] = min(1.0, blended(base['params'][name], parent['params'][name]) * factor)
    return {
        'knobs': knobs,
        'effective_feedback_and_joint_parameters': {**feedback, **joint},
        'shared_drive_nA': study.BASE_DRIVE,
    }


def one_job(job):
    global INITIALIZED, PRODUCTION_SCHEDULE
    params, anchors, destination, production, basin_seed, retain = job
    output = Path(destination)
    output.mkdir(parents=True, exist_ok=True)
    os.chdir(SPINAL)
    with (output / 'process.log').open('w', encoding='utf8') as log, redirect_stdout(log), redirect_stderr(
        log
    ):
        from coupled_muscle_context import use_four_source_muscles
        from fit_reference_feedback import use_compact_nonspiking_network
        import wean_tune_20261016 as study

        if not INITIALIZED:
            use_compact_nonspiking_network()
            import matplotlib.figure

            matplotlib.figure.Figure.savefig = lambda *args, **kwargs: None
            INITIALIZED = True
        spec = fitted_config(params, anchors, study)
        study.apply_config(spec['knobs'], 1.0, fb=spec['effective_feedback_and_joint_parameters'])
        groups = params.get('premotor_group_gains')
        if groups and any(
            k not in ('knee_ext', 'ankle_pf', 'knee_flex', 'ankle_df', 'hip_flex')
            or not np.isfinite(v)
            or v < 0.0
            for k, v in groups.items()
        ):
            raise ValueError(
                'Swing motor gains must be nonnegative and refer to an existing recruited muscle group'
            )
        study.MODS['R'].PREMOTOR_GROUP_GAINS = groups
        schedule = study.MODS['R'].SCHEDULE
        if PRODUCTION_SCHEDULE is None:
            PRODUCTION_SCHEDULE = dict(schedule)
        schedule.update(PRODUCTION_SCHEDULE)
        if production == 'holdout':
            schedule.update(walk=(5.0, 21.0), ramp_down=(21.0, 23.0), stand2=(23.0, 28.0))
        # These are existing pattern-formation-to-motor-neuron synapses. The
        # original fitted table is reloaded by apply_config for each trial.
        if 'pattern_weights' in params:
            import params as neural

            applied_pattern = {}
            for label, multiplier in params['pattern_weights'].items():
                phase, muscle_group = label.split('/')
                if phase not in neural.W_PF_MN or muscle_group not in neural.W_PF_MN[phase]:
                    raise ValueError(f'Unknown existing neural projection {label}')
                neural.W_PF_MN[phase][muscle_group] *= multiplier
                applied_pattern[label] = neural.W_PF_MN[phase][muscle_group]
            spec['pattern_formation_to_motoneuron_conductances'] = applied_pattern
        if basin_seed is not None:
            import basin_gate

            basin_gate.perturb(0.01, np.random.default_rng(basin_seed))
        trace = output / 'full.npz'
        os.environ['AARL_NPZ'] = str(trace)
        physical_load_file = output / 'foot_ground.npz'
        os.environ['AARL_FOOT_FORCE_NPZ'] = str(physical_load_file)
        foot_inputs = params.get('foot_sensor_mode', 'historical')
        os.environ['AARL_FOOT_SENSOR_MODE'] = foot_inputs
        arguments = ['--drive', repr(study.BASE_DRIVE)]
        if not production:
            arguments.insert(0, '--eval')
        regions = params.get('regional_synapses', {})
        with swing_timing(study.MODS['R'], params.get('swing_window')), regional_feedback(
            regions
        ) as reached, use_four_source_muscles(output, 'four', keep_trace=retain):
            study.MODS['R'].main(arguments)
        assert set(reached) == set(regions), 'A requested muscle group received no weights'
        if regions:
            spec['regional_synapse_targets'] = {k: sorted(v) for k, v in reached.items()}
        window = 5.0 if production else 2.0
        reductions = {
            'historical_stored_50': compare_angles(trace, window, recorded_force_scale=1.0),
            'physical_50_N': compare_angles(trace, window, physical_loads=physical_load_file),
        }
        score_mode = params.get('force_convention', 'historical_stored_50')
        angles = reductions[score_mode]
        with np.load(trace, allow_pickle=False) as data:
            stop = np.count_nonzero(data['t'] > 0) + 1
            time, q, com = [data[k][:stop] for k in ('t', 'q', 'com')]
            finite = bool(np.isfinite(q).all() and np.isfinite(com).all())
            minimum_height = float(np.min(com[:, 2]))
            maximum_tilt = float(np.max(np.abs(q[:, 0])))
            duration = float(time[-1])
            with np.load(physical_load_file, allow_pickle=False) as loads:
                recorded_half = data['contact'][:stop] * 0.5
                load_error = float(np.max(np.abs(recorded_half - loads['upward_N'][:stop])))
                normal_error = float(np.max(np.abs(loads['normal_N'][:stop] - loads['upward_N'][:stop])))
            if retain:
                np.savez_compressed(
                    output / 'trajectory.npz',
                    **{k: data[k][:stop] for k in ('t', 'q', 'neuro', 'contact', 'com', 'footz')},
                    **{k: data[k] for k in ('key_joints', 'neuro_names')},
                )
        force_check = read(output / 'force_audit.json')
        assert force_check['force_failure'] is None and len(force_check['records']) == 1
        cycles = [angles['sides'][s]['cycles'] for s in ('r', 'l')]
        errors = [
            abs(angles['sides'][s]['stance_fraction'] - angles['sides'][s]['reference_stance_fraction'])
            for s in ('r', 'l')
        ]
        expected_end = 27.998 if production == 'holdout' else (21.998 if production else 15.998)
        physically_valid = bool(
            finite and duration >= expected_end - 1e-9 and minimum_height >= 0.62 and maximum_tilt <= 40.0
        )
        constraints = [6 - cycles[0], 6 - cycles[1], 0 if physically_valid else 1]
        eligible = max(constraints) <= 0
        mean_error = angles['mean_of_six_joint_RMSE_deg']
        measured = {
            'parameters': params,
            'applied': spec,
            'production': bool(production),
            'schedule': '28s_extended_walking_holdout'
            if production == 'holdout'
            else ('22s_production' if production else '16s_evaluation'),
            'schedule_intervals_s': schedule.copy(),
            'foot_sensor_mode': foot_inputs,
            'basin_seed': basin_seed,
            'constraints': constraints,
            'eligible': eligible,
            'cycles': cycles,
            'stance': [angles['sides'][s]['stance_fraction'] for s in ('r', 'l')],
            'objectives': [mean_error, *errors],
            'absolute_angles': angles,
            'force_convention': score_mode,
            'threshold_comparisons': reductions,
            'physical': {
                'finite': finite,
                'end_s': duration,
                'minimum_COM_m': minimum_height,
                'maximum_absolute_pelvis_tilt_deg': maximum_tilt,
                'active_steps': force_check['records'][0]['steps'],
                'half_stored_minus_upward_max_N': load_error,
                'normal_minus_upward_max_N': normal_error,
            },
            'meets_3pp': bool(eligible and mean_error < ANGLE_LIMIT and max(errors) <= 0.030000000001),
            'meets_2pp': bool(eligible and mean_error < ANGLE_LIMIT and max(errors) <= 0.020000000001),
            'lower_than_frozen_physical_RMSE': bool(
                mean_error is not None and mean_error < FROZEN_PHYSICAL_RMSE
            ),
        }
        if retain:
            measured['constant_walking_drive'] = compare_angles(
                trace, window, physical_loads=physical_load_file, end=schedule['walk'][1]
            )
            measured['full_trace_sha256'] = hashlib.sha256(trace.read_bytes()).hexdigest()
            measured['compact_trace_sha256'] = hashlib.sha256(
                (output / 'trajectory.npz').read_bytes()
            ).hexdigest()
        write(output / 'result.json', measured)
        if not retain:
            trace.unlink()
        return measured


def seeds():
    base = {k: 1.0 for k in BOUNDS}
    for anchor, mix in (
        ('stance', 0.0),
        ('stance', 1.0),
        ('production', 1.0),
        ('symmetric', 1.0),
        ('stance', 0.5),
    ):
        yield {**base, 'anchor': anchor, 'mix': mix}
    # Screen the ankle recruiting weights that differed most in the trial-38 angles.
    for anchor, mix in (('stance', 0.0), ('stance', 0.45), ('stance', 1.0), ('production', 1.0)):
        for plantar, dorsi, tone in ((1.6, 0.7, 1.5), (2.5, 0.45, 2.0), (0.55, 1.8, 0.65)):
            yield {
                **base,
                'anchor': anchor,
                'mix': mix,
                'plantar_flexion': plantar,
                'dorsiflexion': dorsi,
                'ankle_tone': tone,
            }
    for mix in (0.0, 0.35, 0.70, 1.0):
        for right, left in ((0.8, 1.15), (1.15, 0.80)):
            yield {**base, 'anchor': 'stance', 'mix': mix, 'right_drive': right, 'left_drive': left}


def run_search(output, anchors, budget, jobs, population, convention, sampler_name):
    def constraints(trial):
        return trial.user_attrs['constraints']

    if sampler_name == 'tpe':
        sampler = optuna.samplers.TPESampler(
            seed=OPTUNA_SEED + 827, multivariate=True, n_startup_trials=24, constraints_func=constraints
        )
    else:
        sampler = optuna.samplers.NSGAIISampler(
            seed=OPTUNA_SEED, population_size=population, constraints_func=constraints
        )
    store = 'sqlite:///' + (output / 'optuna.sqlite').as_posix()
    study = optuna.create_study(
        directions=['minimize'] * 3,
        sampler=sampler,
        storage=store,
        study_name='absolute_joint_angles_and_stance'
        if convention == 'historical_stored_50'
        else 'physical_50N_and_absolute_joint_angles',
        load_if_exists=True,
    )
    if not study.trials:
        for candidate in seeds():
            study.enqueue_trial(candidate)
        if convention == 'physical_50_N':
            base = {k: 1.0 for k in BOUNDS}
            for anchor in ('r6t10', 'r5t54'):
                if anchor in anchors:
                    for mix in (1.0, 0.55):
                        study.enqueue_trial({**base, 'anchor': anchor, 'mix': mix})
            for mix in (0.0, 0.10, 0.25, 0.40):
                for right, left in ((0.85, 0.85), (0.85, 1.05), (1.05, 0.85), (1.08, 1.08), (0.63, 1.18)):
                    study.enqueue_trial(
                        {**base, 'anchor': 'stance', 'mix': mix, 'right_drive': right, 'left_drive': left}
                    )
    all_trials = study.get_trials(deepcopy=False)
    for previous in all_trials:
        if previous.state == optuna.trial.TrialState.RUNNING:
            stored = output / f'trial_{previous.number:04d}' / 'result.json'
            if stored.exists():
                completed = read(stored)
                trial = optuna.Trial(study, previous._trial_id)
                trial.set_user_attr('constraints', completed['constraints'])
                vals = [v if v is not None else 180.0 for v in completed['objectives']]
                study.tell(trial, vals)
            else:
                study.tell(previous.number, state=optuna.trial.TrialState.FAIL)
    finished = sum(t.state == optuna.trial.TrialState.COMPLETE for t in study.trials)
    queue = {}
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        while finished < budget:
            while len(queue) < min(jobs, budget - finished):
                trial = study.ask()
                options = ['stance', 'symmetric', 'production']
                if convention == 'physical_50_N':
                    options.extend(k for k in ('r6t10', 'r5t54') if k in anchors)
                params = {'anchor': trial.suggest_categorical('anchor', options)}
                params['force_convention'] = convention
                for key, (low, high) in BOUNDS.items():
                    params[key] = trial.suggest_float(key, low, high, log=key != 'mix')
                directory = output / f'trial_{trial.number:04d}'
                queue[pool.submit(one_job, (params, anchors, directory, False, None, False))] = trial
            ready, _ = wait(queue, return_when=FIRST_COMPLETED)
            for future in ready:
                trial = queue.pop(future)
                try:
                    row = future.result()
                except Exception:
                    study.tell(trial, state=optuna.trial.TrialState.FAIL)
                    raise
                trial.set_user_attr('constraints', row['constraints'])
                values = [v if v is not None else 180.0 for v in row['objectives']]
                study.tell(trial, values)
                print(f"{trial.number}: {values} cycles={row['cycles']} pass={row['meets_3pp']}", flush=True)
                finished += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('fit', 'replay'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--anchors', type=Path, default=HERE / 'pareto_reference/anchors.json')
    parser.add_argument('--budget', type=int, default=192)
    parser.add_argument('--jobs', type=int, default=6)
    parser.add_argument('--population', type=int, default=20)
    parser.add_argument('--sampler', choices=('nsga2', 'tpe'), default='nsga2')
    parser.add_argument(
        '--force-convention',
        choices=('physical_50_N', 'historical_stored_50'),
        default='physical_50_N',
        help='Objective used in a newly fitted study; the replay reads its convention from the saved candidate',
    )
    parser.add_argument('--parameters', type=Path)
    parser.add_argument('--production', action='store_true')
    parser.add_argument(
        '--foot-sensor-mode',
        choices=('historical', 'anatomical'),
        help='Override the saved neural input adapter in a separate replay only',
    )
    parser.add_argument(
        '--extended-holdout',
        action='store_true',
        help='Extend the original production walking interval by six seconds, giving a 28 s independent schedule',
    )
    parser.add_argument('--basin-seed', type=int)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    anchors = read(args.anchors.resolve())
    if args.action == 'fit':
        run_search(
            output, anchors, args.budget, args.jobs, args.population, args.force_convention, args.sampler
        )
    elif args.parameters is not None:
        candidate = read(args.parameters.resolve())
        params = candidate.get('parameters', candidate)
        if args.foot_sensor_mode:
            params = {**params, 'foot_sensor_mode': args.foot_sensor_mode}
        if (output / 'result.json').exists():
            parser.error('A replay result already exists at that destination')
        result = one_job(
            (
                params,
                anchors,
                output,
                'holdout' if args.extended_holdout else args.production,
                args.basin_seed,
                True,
            )
        )
        print(json.dumps(result), flush=True)
    else:
        parser.error('A replay requires --parameters')


if __name__ == '__main__':
    main()
