"""Bracket existing right hip spindle asymmetry and bilateral controller drive about one fixed fit."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

from pareto_feedback import HERE, one_job, read, write
from scipy.stats import qmc


def bracket():
    factors = list(product((0.98, 1.0, 1.02), repeat=3))
    sample = qmc.Sobol(d=3, scramble=True, seed=280913).random_base2(4)
    factors.extend(tuple(row) for row in (0.97 + 0.06 * sample))
    return factors


def candidate(base, right_drive, left_drive, right_hip):
    regions = dict(base['regional_synapses'])
    regions['hip_abd/r/ii_to_mn'] *= right_hip
    regions['hip_add/r/ii_to_mn'] /= right_hip
    return {
        **base,
        'right_drive': base['right_drive'] * right_drive,
        'left_drive': base['left_drive'] * left_drive,
        'regional_synapses': regions,
        'force_convention': 'physical_50_N',
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--parent', type=Path, required=True)
    cli.add_argument('--anchors', type=Path, default=HERE / 'pareto_reference/anchors.json')
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--jobs', type=int, default=4)
    args = cli.parse_args()
    parent = read(args.parent)
    parent = parent.get('parameters', parent)
    anchors = read(args.anchors)
    out = args.output.resolve()
    out.mkdir(exist_ok=True, parents=True)
    jobs = []
    changes = []
    for number, scales in enumerate(bracket()):
        params = candidate(parent, *scales)
        folder = out / f'trial_{number:04d}'
        changes.append(
            {'trial': number, 'right_drive': scales[0], 'left_drive': scales[1], 'right_hip': scales[2]}
        )
        if (folder / 'result.json').exists():
            existing = read(folder / 'result.json')
            if existing['parameters'] != params:
                raise ValueError(f'Different controller already occupies {folder}')
        else:
            jobs.append((params, anchors, folder, False, None, False))
    write(out / 'design.json', changes)
    with ProcessPoolExecutor(max_workers=args.jobs) as workers:
        for result in workers.map(one_job, jobs):
            print(result['objectives'], result['stance'], result['cycles'], result['meets_3pp'], flush=True)
    results = {path.parent.name: read(path) for path in sorted(out.glob('trial_*/result.json'))}
    write(out / 'measured.json', results)


if __name__ == '__main__':
    main()
