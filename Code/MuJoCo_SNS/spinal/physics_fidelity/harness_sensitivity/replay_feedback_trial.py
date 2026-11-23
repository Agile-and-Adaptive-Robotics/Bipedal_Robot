"""Replay one archived optimizer evaluation, independent of sampler scheduling."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=(5, 6, 7, 8))
    parser.add_argument('trial', type=int)
    parser.add_argument('--reference', type=Path, default=here / 'reference_feedback')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--arm', choices=('four', 'native'), default='four')
    args = parser.parse_args()
    study = json.loads((args.reference / 'study.json').read_text(encoding='utf8'))['rounds']
    if str(args.round) not in study:
        parser.error('The reference package does not contain that round')
    with (args.reference / 'trials.csv').open(newline='', encoding='utf8') as stream:
        row = next((r for r in csv.DictReader(stream)
                    if int(r['round']) == args.round and int(r['trial']) == args.trial), None)
    if row is None:
        parser.error('The requested trial is absent from the selected round')
    settings = study[str(args.round)]['winners'][0]['fit']
    params = settings['params'].copy()
    for key in params:
        if row.get(key):
            params[key] = float(row[key])
    if args.round in (5, 7):
        knobs = {k: float(row[k]) for k in settings['knobs']}
    else:
        parent = 5 if args.round == 6 else 6
        knobs = study[str(parent)]['winners'][0]['fit']['knobs']
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        parser.error('Choose an empty output directory')
    spec = out / 'settings.json'
    spec.write_text(json.dumps({'round': args.round, 'trial': args.trial, 'knobs': knobs,
                                'scale': float(row.get('drive_scale') or 1.), 'params': params}, indent=2),
                    encoding='utf8')
    call = [sys.executable, '-X', 'utf8', str(here / 'verify_reference_feedback.py'), str(spec),
            '--output', str(out), '--arm', args.arm]
    if args.round == 7:
        call.append('--production')
    subprocess.run(call, check=True)


if __name__ == '__main__':
    main()
