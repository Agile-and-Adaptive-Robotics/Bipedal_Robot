"""Collect completed controller evaluations, their settings, and the non-dominated physical results."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


def fingerprint(parameters):
    return hashlib.sha256(json.dumps(parameters, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def gather(label, source):
    if source.is_dir():
        entries = (
            (file.parent.name, json.loads(file.read_text()))
            for file in sorted(source.glob('trial_*/result.json'))
        )
    else:
        entries = json.loads(source.read_text()).items()
    for name, result in entries:
        angle, right, left = result['objectives']
        mode = result.get(
            'force_convention', result['parameters'].get('force_convention', 'historical_stored_50')
        )
        eligible = bool(result['eligible'] and angle is not None)
        spec = result['parameters']
        legs = result['absolute_angles']['sides']
        yield {
            'study': label,
            'trial': int(name.removeprefix('trial_')),
            'convention': mode,
            'eligible': eligible,
            'objectives': result['objectives'],
            'stance': result['stance'],
            'cycles': result['cycles'],
            'joint_RMSE_deg': {s: legs[s]['joint_RMSE_deg'] for s in 'rl'},
            'period_s': [legs[s]['period_s'] for s in 'rl'],
            'three_pp': bool(eligible and max(right, left) <= 0.030000000001),
            'two_pp': bool(eligible and max(right, left) <= 0.020000000001),
            'target': bool(eligible and angle < 17.78 and max(right, left) <= 0.030000000001),
            'target_two_pp': bool(eligible and angle < 17.78 and max(right, left) <= 0.020000000001),
            'physical': result['physical'],
            'parameters': spec,
            'parameter_hash': fingerprint(spec),
        }


def front(rows):
    keep = []
    values = np.asarray([r['objectives'] for r in rows])
    seen = set()
    for i, row in enumerate(rows):
        triple = tuple(np.round(values[i], 10))
        if triple in seen:
            continue
        if not np.any(np.all(values <= values[i], axis=1) & np.any(values < values[i], axis=1)):
            keep.append(row)
            seen.add(triple)
    return sorted(keep, key=lambda r: (r['objectives'][0], max(r['objectives'][1:])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--source', nargs=2, action='append', metavar=('LABEL', 'DIRECTORY_OR_INDEX'), required=True
    )
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    records = [r for label, source in args.source for r in gather(label, Path(source))]
    with (out / 'trials.json.gz').open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as zipped:
            zipped.write(json.dumps(records, separators=(',', ':'), allow_nan=False).encode())
    by_study = {}
    for label, _ in args.source:
        rows = [r for r in records if r['study'] == label]
        by_study[label] = {
            'evaluations': len(rows),
            'convention': rows[0]['convention'] if rows else None,
            'eligible': sum(r['eligible'] for r in rows),
            'angle_and_three_pp': sum(r['target'] for r in rows),
            'both_stance_within_three_pp': sum(r['three_pp'] for r in rows),
            'both_stance_within_two_pp': sum(r['two_pp'] for r in rows),
        }
    physical = [r for r in records if r['eligible'] and r['convention'] == 'physical_50_N']
    frontier = front(physical)
    if not physical:
        raise ValueError('There are no eligible single-counted physical trials')
    pp3 = [r for r in physical if r['three_pp']]
    pp2 = [r for r in physical if r['two_pp']]
    angle_order = lambda r: (r['objectives'][0], max(r['objectives'][1:]), r['study'], r['trial'])

    def gap(r):
        e, a, b = r['objectives']
        return (
            (max(0.0, e - 17.78) / 3.0) ** 2
            + (max(0.0, a - 0.03) / 0.07) ** 2
            + (max(0.0, b - 0.03) / 0.07) ** 2,
            e,
        )

    selected = {'minimum_angle': min(physical, key=angle_order), 'balanced': min(physical, key=gap)}
    if pp3:
        selected['best_angle_at_three_pp'] = min(pp3, key=angle_order)
    if pp2:
        selected['best_angle_at_two_pp'] = min(pp2, key=angle_order)
    over = [r for r in physical if r['target']]
    tight = [r for r in physical if r['target_two_pp']]
    if over:
        selected['joint_target'] = min(tight or over, key=angle_order)
    summary = {
        'studies': by_study,
        'physical_evaluations': sum(r['convention'] == 'physical_50_N' for r in records),
        'physical_eligible': len(physical),
        'non_dominated_unique_objective_vectors': len(frontier),
        'physical_target_hits': [{'study': r['study'], 'trial': r['trial']} for r in over],
        'selection': 'Rank completed physical 16 s evaluations only. Prioritize the lowest angle within the two-point stance interval when there is a combined hit. The balanced example minimizes normalized shortfall using 3 degrees and 7 percentage points. Longer schedules are reserved for verification.',
        'selected': {k: {a: v for a, v in row.items() if a != 'parameters'} for k, row in selected.items()},
    }
    (out / 'outcomes.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    selected_dir = out / 'selected'
    selected_dir.mkdir(exist_ok=True)
    for label, row in selected.items():
        (selected_dir / (label + '.json')).write_text(json.dumps(row['parameters'], indent=2) + '\n')
    with (out / 'front.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                'study',
                'trial',
                'mean_joint_RMSE_deg',
                'right_stance',
                'left_stance',
                'right_stance_error',
                'left_stance_error',
                'right_cycles',
                'left_cycles',
                'three_pp',
                'two_pp',
                'combined_target',
            ]
        )
        for r in frontier:
            writer.writerow(
                [
                    r['study'],
                    r['trial'],
                    r['objectives'][0],
                    *r['stance'],
                    *r['objectives'][1:],
                    *r['cycles'],
                    r['three_pp'],
                    r['two_pp'],
                    r['target'],
                ]
            )
    print(
        json.dumps(
            {
                'studies': by_study,
                'selected': {k: [v['study'], v['trial'], v['objectives']] for k, v in selected.items()},
            }
        )
    )


if __name__ == '__main__':
    main()
