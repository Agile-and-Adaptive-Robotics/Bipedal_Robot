"""Local Pareto refinement of anatomically grouped existing afferent projections."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

from pareto_feedback import ANGLE_LIMIT, BOUNDS, HERE, np, one_job, read, write


def nondominated(rows):
    return [
        row
        for row in rows
        if not any(
            all(a <= b for a, b in zip(other['objectives'], row['objectives']))
            and any(a < b for a, b in zip(other['objectives'], row['objectives']))
            for other in rows
        )
    ]


def gap(row):
    angle, right, left = row['objectives']
    return np.linalg.norm(
        [
            max(0.0, (angle - ANGLE_LIMIT) / 4.0),
            max(0.0, (right - 0.03) / 0.08),
            max(0.0, (left - 0.03) / 0.08),
        ]
    )


def regional(pf, asymmetry, knee, load, frontal):
    factors = {}
    for side, ankle in [('r', pf / asymmetry), ('l', pf * asymmetry)]:
        for group, strength in [
            ('ankle_pf', ankle),
            ('ankle_df', 1.0 / ankle),
            ('knee_ext', knee),
            ('knee_flex', 1.0 / knee),
        ]:
            for gain, value in [
                ('ii_to_mn', strength),
                ('ia_to_mn', np.sqrt(strength)),
                ('ib_to_mn_inh', strength ** (-load)),
                ('posture_to_mn', strength),
            ]:
                factors[f'{group}/{side}/{gain}'] = float(np.clip(value, 0.07, 14.0))
    for group in ('hip_abd', 'hip_add'):
        factors[f'{group}/both/ii_to_mn'] = frontal if group == 'hip_abd' else 1.0 / frontal
    return factors


def initial():
    default = {**{k: 1.0 for k in BOUNDS}, 'anchor': 'stance'}
    tests = [
        (1.0, 1.0, 1.0, 0.6, 1.0),
        (2.2, 1.0, 1.0, 0.6, 1.0),
        (4.0, 1.0, 1.35, 0.7, 1.0),
        (0.45, 1.0, 1.0, 0.6, 1.0),
        (2.8, 1.35, 1.8, 0.8, 2.0),
        (1.0, 1.0, 1.0, 0.6, 3.0),
        (1.0, 1.0, 1.0, 0.6, 0.35),
    ]
    for factors in tests:
        for mix in (1.0, 0.0, 0.55):
            yield {**default, 'mix': mix, 'regional_synapses': regional(*factors)}, list(factors)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--anchors', type=Path, default=HERE / 'pareto_reference/anchors.json')
    cli.add_argument('--seed-search', type=Path, required=True)
    cli.add_argument('--budget', type=int, default=100)
    cli.add_argument('--jobs', type=int, default=8)
    cli.add_argument('--concurrency', type=Path)
    cli.add_argument('--proposals', type=Path)
    cli.add_argument('--seed', type=int, default=96173)
    cli.add_argument(
        '--force-convention', choices=('physical_50_N', 'historical_stored_50'), default='physical_50_N'
    )
    args = cli.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    anchors = read(args.anchors.resolve())
    initial_rows = list(initial())
    completed = {}
    parent_cache = {}
    running = {}
    submitted = set()
    for p in out.glob('trial_*/result.json'):
        n = int(p.parent.name[6:])
        completed[n] = read(p)
        submitted.add(n)

    def read_external():
        for p in args.seed_search.glob('trial_*/result.json'):
            if p not in parent_cache:
                try:
                    parent_cache[p] = read(p)
                except (OSError, ValueError):
                    pass
        return [
            row
            for row in parent_cache.values()
            if row.get('force_convention', 'historical_stored_50') == args.force_convention
        ]

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        while len(completed) < args.budget:
            active = int(args.concurrency.read_text().strip()) if args.concurrency else args.jobs
            while len(running) < min(active, args.budget - len(completed)):
                n = next(i for i in range(args.budget) if i not in submitted)
                submitted.add(n)
                rng = np.random.default_rng(args.seed + n * 997)
                case = out / f'trial_{n:04d}'
                case.mkdir(exist_ok=True)
                supplied = (
                    []
                    if not args.proposals
                    else [
                        p
                        for p in sorted(args.proposals.glob('*.json'))
                        if not (out / (p.stem + '.used')).exists()
                    ]
                )
                if supplied:
                    inp = supplied[0]
                    spec = read(inp)
                    p = spec.get('parameters', spec)
                    factors = spec.get('regional_factors', [1.0, 1.0, 1.0, 0.6, 1.0])
                    lineage = {'external_nomination': inp.name}
                    (out / (inp.stem + '.used')).write_text(str(n))
                elif n < len(initial_rows):
                    p, factors = initial_rows[n]
                    lineage = {'initial_design': n}
                else:
                    pool_rows = [r for r in read_external() + list(completed.values()) if r['eligible']]
                    frontier = nondominated(pool_rows)
                    best = sorted(frontier, key=gap)
                    if n % 4 != 0:
                        best = best[: max(2, len(best) // 3)]
                    parent = best[int(rng.integers(len(best)))]
                    p = parent['parameters'].copy()
                    existing = parent.get('regional_factors', [1.0, 1.0, 1.0, 0.6, 1.0])
                    factors = list(existing)
                    step = rng.choice([0.12, 0.28, 0.52, 0.78])
                    selected = rng.choice(5, int(rng.integers(1, 4)), replace=False)
                    for i in selected:
                        factors[i] = float(
                            np.clip(
                                factors[i] * np.exp(rng.normal(0.0, step)),
                                [0.13, 0.48, 0.25, 0.12, 0.15][i],
                                [8.0, 2.3, 5.0, 1.9, 5.0][i],
                            )
                        )
                    for key in rng.choice(
                        [
                            'mix',
                            'right_drive',
                            'left_drive',
                            'premotor_period',
                            'knee_extension',
                            'plantar_flexion',
                        ],
                        int(rng.integers(1, 4)),
                        replace=False,
                    ):
                        lo, hi = BOUNDS[key]
                        p[key] = float(
                            np.clip(
                                (p[key] + rng.normal(0, 0.075))
                                if key == 'mix'
                                else p[key] * np.exp(rng.normal(0, 0.12)),
                                lo,
                                hi,
                            )
                        )
                    p['regional_synapses'] = regional(*factors)
                    lineage = {
                        'parent_objectives': parent['objectives'],
                        'choice': 'non-dominated parents with emphasis on the joint target',
                    }
                p = {**p, 'force_convention': args.force_convention}
                write(case / 'proposal.json', {'parameters': p, 'regional_factors': factors, **lineage})
                running[pool.submit(one_job, (p, anchors, case, False, None, False))] = (n, factors)
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for f in done:
                n, factors = running.pop(f)
                r = f.result()
                r['regional_factors'] = factors
                write(out / f'trial_{n:04d}' / 'result.json', r)
                completed[n] = r
                print(n, r['objectives'], r['cycles'], r['meets_3pp'], flush=True)


if __name__ == '__main__':
    main()
