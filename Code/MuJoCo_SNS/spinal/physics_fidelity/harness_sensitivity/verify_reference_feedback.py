"""Fresh coupled finalist or basin-perturbation replay of the saved feedback weights."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np

from coupled_muscle_context import SPINAL, use_four_source_muscles
from fit_reference_feedback import use_compact_nonspiking_network


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('winner', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--arm', choices=['four', 'native'], default='four')
    parser.add_argument('--production', action='store_true', help='Use the 22 s schedule and score its walking window')
    parser.add_argument('--basin-seed', type=int, help='Perturb scalar network parameters by ±1%% using the original basin gate')
    parser.add_argument('--scale', type=float, help='Override only the four supplementary left/right drive gains')
    args = parser.parse_args()
    winner = args.winner.resolve()
    out = args.output.resolve()
    if (out / 'result.json').exists():
        parser.error('The output already contains a completed result')
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(SPINAL)
    sys.path.insert(0, str(SPINAL))
    use_compact_nonspiking_network()
    import wean_tune_20261016 as study

    settings = json.loads(winner.read_text())
    feedback = {key: settings['params'][key] for key in study.FB_KEYS}
    scale = settings.get('scale', 1.) if args.scale is None else args.scale
    if not (0. <= scale <= 1.):
        parser.error('Supplementary drive scale must lie between zero and one')
    study.apply_config(settings['knobs'], scale, fb=feedback)
    if args.basin_seed is not None:
        import basin_gate
        basin_gate.perturb(.01, np.random.default_rng(args.basin_seed))
    output = out / 'trajectory.npz'
    os.environ['AARL_NPZ'] = str(output)
    runner = study.MODS['R']
    if args.production:
        # The inherited production runner also writes a default PNG into the repository.
        # This task already records the complete numeric trace in its isolated destination.
        import matplotlib.figure
        matplotlib.figure.Figure.savefig = lambda *args, **kwargs: None
    with use_four_source_muscles(out, args.arm):
        inputs = ['--drive', repr(study.BASE_DRIVE)]
        if not args.production:
            inputs.insert(0, '--eval')
        result = runner.main(inputs)
    if args.production:
        import kine_ref
        with np.load(output) as arrays:
            samples = int(np.sum(arrays['t'] > 0.0)) + 1
            t = arrays['t'][:samples]
            measurements = kine_ref.compare(t, arrays['q'][:samples], arrays['neuro'][:samples], 5.0,
                                            ref=kine_ref.ref_cached(), contact=arrays['contact'][:samples]) or {}
            result = {'nan': not np.isfinite(arrays['q'][:samples]).all(), 't_end': float(t[-1]),
                      'kz': float(np.min(arrays['com'][:samples, 2])), 'kine': measurements,
                      'kine_score': measurements.get('kine_score')}
    metrics = study.flat_metrics(result)
    import kine_ref
    with np.load(output) as arrays:
        stop = int(np.count_nonzero(arrays['t'] > 0.)) + 1
        t, degrees, neuro, contacts = [arrays[key][:stop] for key in ('t', 'q', 'neuro', 'contact')]
    window = 5. if args.production else 2.
    detected = {side: kine_ref.sim_side(t, degrees, neuro, window, side, contacts) for side in ('r', 'l')}
    metrics.update(winner_sha256=hashlib.sha256(winner.read_bytes()).hexdigest(),
                   trace_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                   arm=args.arm, basin_seed=args.basin_seed, supplementary_scale=scale,
                   shared_drive_nA=study.BASE_DRIVE,
                   detected_cycles={side: int(leg[5]) for side, leg in detected.items()},
                   kinematic_score_available=any(leg[0] is not None for leg in detected.values()),
                   onset_source={side: 'foot loading' if leg[4] else 'RG extensor threshold'
                                 for side, leg in detected.items()},
                   bilateral_basin_gate=bool(not result['nan'] and
                                             (metrics.get('k_n_cycles_r') or 0) >= 3 and
                                             (metrics.get('k_n_cycles_l') or 0) >= 3))
    (out / 'result.json').write_text(json.dumps(metrics, indent=2), encoding='utf-8')
    print(f"Score {metrics.get('kine_score')}, right/left cycles {metrics.get('k_n_cycles_r')}/{metrics.get('k_n_cycles_l')}")


if __name__ == '__main__':
    main()
