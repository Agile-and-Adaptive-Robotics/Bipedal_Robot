"""Score absolute hip, knee, and ankle angle errors using the published 50 N cycle reduction."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

SPINAL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SPINAL))
import kine_ref


def compare_angles(recording, start=2.0, recorded_force_scale=1.0, physical_loads=None, end=None):
    """Reduce a trace using calibrated force samples or explicitly scaled historical channels.

    The unit-scale positional default retains the earlier Python API. For an
    archived Gait2392 trace without the force sidecar, pass 0.5 for physical
    newtons. New trials should pass their calibrated sidecar via physical_loads.
    """
    if not np.isfinite(recorded_force_scale) or recorded_force_scale <= 0.0:
        raise ValueError('The recorded-to-physical contact force scale must be positive')
    reference = kine_ref.load_reference()
    sides = {}
    with np.load(recording, allow_pickle=False) as data:
        count = np.count_nonzero(data['t'] > 0) + 1
        t, degrees, neural, contact = [data[key][:count] for key in ('t', 'q', 'neuro', 'contact')]
    if physical_loads is None:
        contact = contact * recorded_force_scale
    else:
        if recorded_force_scale != 1.0:
            raise ValueError('A calibrated foot-force file is already in newtons')
        with np.load(physical_loads, allow_pickle=False) as loads:
            np.testing.assert_array_equal(loads['t'][:count], t)
            contact = loads['upward_N'][:count]
    if end is not None:
        if end <= start:
            raise ValueError('The analysis end must follow the start')
        selected = t < end
        t, degrees, neural, contact = [values[selected] for values in (t, degrees, neural, contact)]
    for side in ('r', 'l'):
        curve, periods, duty, _, used_contact, cycles, _ = kine_ref.sim_side(
            t, degrees, neural, start, side, contact
        )
        valid = bool(used_contact and cycles)
        errors = (
            {
                joint: float(np.sqrt(np.mean((curve[joint] - reference[side][joint]) ** 2)))
                for joint in ('hip', 'knee', 'ankle')
            }
            if valid
            else None
        )
        sides[side] = {
            'cycles': cycles if valid else 0,
            'period_s': float(np.mean(periods)) if valid else None,
            'reference_period_s': reference['T_' + side],
            'stance_fraction': duty,
            'reference_stance_fraction': reference['duty_' + side],
            'joint_RMSE_deg': errors,
        }
    errors = [value for side in sides.values() for value in (side['joint_RMSE_deg'] or {}).values()]
    return {
        'mean_of_six_joint_RMSE_deg': float(np.mean(errors)) if len(errors) == 6 else None,
        'sides': sides,
        'reference': 'subject01_walk1_ik.mot, subject01_walk1_grf.mot',
        'window_start_s': start,
        'window_end_exclusive_s': end,
        'load_threshold_after_scaling': kine_ref.LOAD_N,
        'recorded_force_scale': recorded_force_scale,
        'force_source': 'single-counted upward ground reaction (N)'
        if physical_loads is not None
        else 'stored contact channels times recorded_force_scale',
        'cycle_interval_s': [0.35, 2.5],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('recording', type=Path)
    parser.add_argument('--start', type=float, default=2.0)
    parser.add_argument('--end', type=float, help='Exclusive end of an optional constant-drive window')
    convention = parser.add_mutually_exclusive_group()
    convention.add_argument(
        '--force-convention',
        choices=('physical', 'historical'),
        help='On an original spinal runner trace, physical uses half the stored heel-plus-toe sum; historical reproduces the original stored-50 threshold',
    )
    convention.add_argument(
        '--foot-loads', type=Path, help='Calibrated AARL_FOOT_FORCE_NPZ file from this identical run'
    )
    convention.add_argument(
        '--recorded-force-scale',
        type=float,
        help='Use 0.5 for a physical 50 N threshold on historical runner heel-plus-toe sums, or 1 for calibrated forces or the historical scoring convention',
    )
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    scale = (
        args.recorded_force_scale
        if args.recorded_force_scale is not None
        else (1.0 if args.foot_loads or args.force_convention == 'historical' else 0.5)
    )
    result = compare_angles(args.recording, args.start, scale, args.foot_loads, args.end)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
    print(json.dumps(result))
