"""Fresh-process perturbations and independent schedules for a fixed Pareto candidate."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('candidate', type=Path)
    cli.add_argument('--output', required=True, type=Path)
    cli.add_argument('--jobs', type=int, default=3)
    cli.add_argument('--nominal-only', action='store_true')
    cli.add_argument('--extended-holdout', action='store_true')
    cli.add_argument('--anchors', type=Path)
    args = cli.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    source = json.loads(args.candidate.read_text())
    parameters = {**source.get('parameters', source), 'force_convention': 'physical_50_N'}
    spec = out / 'parameters.json'
    content = json.dumps(parameters, sort_keys=True, indent=2) + '\n'
    if spec.exists() and spec.read_text() != content:
        raise ValueError('Existing verification belongs to a different controller')
    spec.write_text(content)
    runner = Path(__file__).with_name('pareto_feedback.py')
    cases = [('nominal_16s', []), ('nominal_22s', ['--production'])]
    if args.extended_holdout:
        cases.append(('holdout_28s', ['--extended-holdout']))
    if not args.nominal_only:
        cases.extend((f'perturbation_{n}', ['--basin-seed', str(n)]) for n in range(3000, 3008))

    def check(pair):
        name, opts = pair
        dest = out / name
        dest.mkdir(exist_ok=True)
        if not (dest / 'result.json').exists():
            call = [
                sys.executable,
                '-X',
                'utf8',
                str(runner),
                'replay',
                '--parameters',
                str(spec),
                '--output',
                str(dest),
                *opts,
            ]
            if args.anchors:
                call += ['--anchors', str(args.anchors.resolve())]
            with (dest / 'launch.log').open('w') as log:
                subprocess.run(call, stdout=log, stderr=subprocess.STDOUT, check=True)
        res = json.loads((dest / 'result.json').read_text())
        assert res['force_convention'] == 'physical_50_N' and res['parameters'] == parameters
        keys = (
            'schedule',
            'basin_seed',
            'cycles',
            'stance',
            'objectives',
            'eligible',
            'meets_3pp',
            'meets_2pp',
            'physical',
            'threshold_comparisons',
            'constant_walking_drive',
            'compact_trace_sha256',
            'full_trace_sha256',
        )
        return name, {k: res[k] for k in keys}

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = dict(pool.map(check, cases))
    bundled = {'parameters_sha256': hashlib.sha256(spec.read_bytes()).hexdigest(), 'cases': rows}
    (out / 'summary.json').write_text(json.dumps(bundled, indent=2) + '\n')
    print(json.dumps({k: [v['objectives'], v['cycles'], v['meets_3pp']] for k, v in rows.items()}))


if __name__ == '__main__':
    main()
