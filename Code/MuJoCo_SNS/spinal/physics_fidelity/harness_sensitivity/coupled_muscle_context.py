"""Use the source Thelen muscles in an existing Gait2392 forward simulation.

Examples (from spinal, in the myo environment):

    python physics_fidelity/harness_sensitivity/coupled_muscle_context.py --output D:/temp/selective_eval runner.py --eval --drive 2.0773198435110425
    python physics_fidelity/harness_sensitivity/coupled_muscle_context.py --output D:/temp/air controller_profiles_20261104/body_probe.py --reflex-fit FIT --stage 2 --expression .1 --output D:/temp/air/capture

The second output argument belongs to body_probe. The first receives the muscle
forces, equilibrium checks, and a cryptographic signature of the common start.
The native arm is the same experiment with ordinary MuJoCo actuator mechanics.
"""
from contextlib import contextmanager
import argparse
import json
import os
from pathlib import Path
import runpy
import sys

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
FIDELITY = HERE.parent
SPINAL = FIDELITY.parent
sys.path[:0] = [str(HERE), str(FIDELITY)]
from selective_force import MUSCLES, SelectiveMuscleStep, record_start


def is_coupled_step(frame):
    """Exclude auxiliary MuJoCo steps made during pose and support setup."""
    while frame is not None:
        name = frame.f_code.co_name
        if name == 'main' and frame.f_code.co_filename.endswith('runner.py'):
            return 'k' in frame.f_locals and 't' in frame.f_locals
        if name == 'pilot' and frame.f_code.co_filename.endswith('paper_checked_variant.py'):
            return 'plog' in frame.f_locals
        if name == 'simulate' and frame.f_code.co_filename.endswith('dirusso_recovery.py'):
            return 'log' in frame.f_locals and 'i' in frame.f_locals
        frame = frame.f_back
    return False


@contextmanager
def use_four_source_muscles(out, arm='four', keep_trace=True, record_activation=False):
    """Advance each selected fiber once per forward step in a single process."""
    if arm not in ('native', 'four'):
        raise ValueError('Choose native or four')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    original = mujoco.mj_step
    runs = {}
    failure = None

    def advance(model, data, *args, **kwargs):
        nonlocal failure
        if not is_coupled_step(sys._getframe(1)):
            return original(model, data, *args, **kwargs)
        identifier = id(model)
        if identifier not in runs:
            ids = np.array([model.actuator(name).id for name in MUSCLES], dtype=int)
            if (model.na != 92 or model.nu != 92 or
                    not np.array_equal(model.actuator_actadr[ids], ids) or
                    np.any(model.actuator_gaintype[ids] != mujoco.mjtGain.mjGAIN_MUSCLE) or
                    np.any(model.actuator_biastype[ids] != mujoco.mjtBias.mjBIAS_MUSCLE) or
                    np.any(model.actuator_gainprm[ids, 2] <= 0.0)):
                raise ValueError('The four named active muscle paths must belong to the native 92-actuator Gait2392 body')
            record_start(model, data, ids, MUSCLES, out / f'start_{len(runs)}.json')
            runs[identifier] = {
                'model': model, 'steps': 0, 'time': None, 'activation': [],
                'readout': [model.actuator_actadr[model.actuator(name).id]
                            for name in ('vas_lat_r', 'semimem_r')],
                'muscles': SelectiveMuscleStep(model, data, ids, MUSCLES, original) if arm == 'four' else None,
            }
        current = runs[identifier]
        muscles = current['muscles']
        try:
            if muscles is None:
                original(model, data, *args, **kwargs)
            else:
                muscles.advance(model, data, *args, **kwargs)
        except Exception as error:
            failure = {'reason': repr(error), 'time': float(data.time), 'steps': current['steps']}
            if muscles is not None:
                failure.update(max_equilibrium_error_N=float(np.max(np.abs(muscles.bank.error))),
                               bounds=muscles.bank.bounds.tolist())
            raise
        current['steps'] += 1
        current['time'] = float(data.time)
        if record_activation:
            current['activation'].append(np.r_[data.time, data.act[current['readout']]])
        if muscles is not None and not keep_trace:
            muscles.samples.clear()

    mujoco.mj_step = advance
    try:
        yield
    finally:
        mujoco.mj_step = original
        records = []
        for index, current in enumerate(runs.values()):
            muscles = current['muscles']
            records.append({'steps': current['steps'], 'last_time': current['time'],
                            'sampled': len(muscles.samples) if muscles is not None else 0})
            if muscles is not None and keep_trace:
                np.savez_compressed(out / f'muscles_{index}.npz', samples=np.array(muscles.samples), names=MUSCLES)
            if record_activation:
                np.savez_compressed(out / f'activation_{index}.npz', samples=np.array(current['activation']))
        report = {'arm': arm, 'paths': MUSCLES, 'records': records, 'force_failure': failure}
        (out / 'force_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--arm', choices=('native', 'four'), default='four')
    parser.add_argument('--activation', action='store_true', help='Save vas_lat_r and semimem_r actuator activation every step')
    parser.add_argument('script', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    output = args.output.resolve()
    script = args.script.resolve()
    if (output / 'force_audit.json').exists():
        parser.error(f'This output already contains a completed or stopped trial: {output}')
    sys.path[:0] = [str(script.parent), str(SPINAL), str(SPINAL / 'functional_gait2392')]
    os.chdir(SPINAL)
    os.environ['AARL_NPZ'] = str(output / 'trajectory.npz')
    sys.argv = [str(script), *args.arguments]
    with use_four_source_muscles(output, args.arm, record_activation=args.activation):
        runpy.run_path(str(script), run_name='__main__')
    if not json.loads((output / 'force_audit.json').read_text())['records']:
        raise RuntimeError('No forward Gait2392 steps were observed in the requested evaluator')


if __name__ == '__main__':
    main()
