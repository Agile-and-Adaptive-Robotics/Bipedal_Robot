"""Repeat, basin, and frozen-tuning support ladder; independent processes."""
import concurrent.futures
import json
import os
from pathlib import Path
import numpy as np
from recover_w2l import HERE, run_child


def main():
    recovered = json.loads((HERE / 'recovery_results.json').read_text())
    source = [x for x in recovered if x['job']['knobs']['dur_walk'] == 30.]
    jobs = []
    for item in source:
        for support in (0.55, 0.40, 0.25, 0.10, 0.):
            knobs = dict(item['job']['knobs'], wsup=support)
            jobs.append(dict(name=item['job']['name'] + f'_ws{round(100*support)}',
                             knobs=knobs, kind='frozen_support_ladder'))
    best = source[0]
    for i in range(3):
        jobs.append(dict(name=f'repeat_{i}', knobs=best['job']['knobs'],
                         kind='repeat'))
    rng = np.random.default_rng(20261027)
    for i in range(12):
        knobs = dict(best['job']['knobs'])
        # Every searched scalar is perturbed, with independent +/-1% draws.
        for key in best['expected_gains']:
            knobs[key] *= rng.uniform(.99, 1.01)
        for key in ('ch_ank', 'ch_contra', 'comm', 'ctn_amp', 'cref', 'ibnA'):
            if key not in best['expected_gains']:
                knobs[key] *= rng.uniform(.99, 1.01)
        jobs.append(dict(name=f'basin_{i:02}', knobs=knobs, kind='basin'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(run_child, jobs))
    (HERE / 'validation_results.json').write_text(json.dumps(results, indent=2),
                                                 encoding='utf-8')


if __name__ == '__main__':
    main()
