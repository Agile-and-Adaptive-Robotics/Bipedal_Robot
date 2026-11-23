"""Run a complete four-muscle feedback curriculum in an isolated results directory.

From spinal, run `python physics_fidelity/harness_sensitivity/fit_reference_feedback.py
--output D:/temp/gait_feedback` using the myo environment. Rounds retain the original
parameter ranges, score, bilateral eligibility, seeding, and four worker seeds.
Each completed worker resumes from its JSON-lines record. Use a separate directory
for each experiment. Pass --arm native to obtain a matched native-muscle study.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parents[1]
sys.path[:0] = [str(HERE), str(SPINAL), str(SPINAL / 'w2l_cpg')]
from coupled_muscle_context import use_four_source_muscles

TRIALS_PER_WORKER = {5: 22, 6: 20, 7: 16, 8: 12}


def use_compact_nonspiking_network():
    import build_network
    from compact_pairwise_eval import CompactPairwiseForward

    original = build_network.build

    def compact_build(*args, **kwargs):
        network = original(*args, **kwargs)
        solver = network.compiled
        if solver.spiking or solver.electrical or solver.electrical_rectified:
            raise ValueError('The paired-order calculation applies only to this nonspiking network')
        solver.forward = CompactPairwiseForward(solver).forward
        return network

    build_network.build = compact_build


def study_worker(args):
    import kine_ref
    import runner
    import wean_tune_20261016 as study

    use_compact_nonspiking_network()
    original_comparison = kine_ref.compare

    def registered_comparison(*call_args, **call_kwargs):
        result = original_comparison(*call_args, **call_kwargs)
        # The original search assigned its infeasible score to left-only runs.
        # Preserve its sampling and ranking; fresh reports retain real leg counts.
        if result and not result.get('n_cycles_r') and result.get('n_cycles_l'):
            return None
        return result

    kine_ref.compare = registered_comparison
    study.OUT = str(args.output)
    study.DB = 'sqlite:///' + (args.output / 'optuna.sqlite3').as_posix()
    if args.select:
        sys.argv = [str(__file__), 'best', str(args.round), str(args.select)]
        study.cmd_best()
        return
    name = f'round{args.round}_w{args.worker}'
    audits = args.output / 'forces' / name
    audits.mkdir(parents=True, exist_ok=True)
    original = runner.main

    def measured(*run_args, **run_kwargs):
        trial_path = audits / f'trial_{len(list(audits.iterdir())):03d}'
        with use_four_source_muscles(trial_path, args.arm, keep_trace=False):
            metrics = original(*run_args, **run_kwargs)
        audit = json.loads((trial_path / 'force_audit.json').read_text())
        if len(audit['records']) != 1 or not audit['records'][0]['steps'] or audit['force_failure']:
            raise RuntimeError(f'No verified four-muscle trial in {trial_path}')
        return metrics

    runner.main = measured
    log = args.output / (name + '.jsonl')
    completed = len(log.read_text().splitlines()) if log.exists() else 0
    remaining = TRIALS_PER_WORKER[args.round] - completed
    if remaining:
        sys.argv = [str(__file__), 'round', str(args.round), str(remaining), str(args.worker)]
        study.cmd_round()


def study_schedule(args):
    common = [sys.executable, '-X', 'utf8', str(Path(__file__).resolve()),
              '--output', str(args.output), '--arm', args.arm]
    for stage in TRIALS_PER_WORKER:
        if (args.output / f'round{stage}_best.json').exists():
            continue
        processes = {}
        for worker in range(4):
            if worker == 1:
                first = args.output / f'round{stage}_w0.jsonl'
                while not first.exists() and processes[0].poll() is None:
                    time.sleep(2)
                if processes[0].poll() not in (None, 0):
                    raise RuntimeError(f'The first worker exited before it seeded round {stage}')
            while sum(p.poll() is None for p in processes.values()) >= args.jobs:
                time.sleep(2)
            with (args.output / f'round{stage}_w{worker}.log').open('a', encoding='utf-8') as log:
                cmd = [*common, '--round', str(stage), '--worker', str(worker)]
                processes[worker] = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        for worker, process in processes.items():
            if process.wait():
                raise RuntimeError(f'Round {stage}, worker {worker} failed. Details are in its log.')
        # A round without an eligible bilateral finalist cannot seed a dependent round.
        for rank in (1, 2, 3):
            selection = subprocess.run([*common, '--round', str(stage), '--select', str(rank)])
            if selection.returncode:
                if rank == 1:
                    selection.check_returncode()
                break


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--arm', choices=('four', 'native'), default='four')
    parser.add_argument('--jobs', type=int, choices=(1, 2, 3, 4), default=4)
    parser.add_argument('--round', type=int, choices=TRIALS_PER_WORKER)
    parser.add_argument('--worker', type=int, choices=range(4))
    parser.add_argument('--select', type=int, choices=(1, 2, 3))
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    os.chdir(SPINAL)
    if args.round is None:
        record = args.output / 'study_config.json'
        config = {'arm': args.arm, 'trials_per_worker': TRIALS_PER_WORKER,
                  'source_study': 'wean_tune_20261016.py', 'workers': 4,
                  'left_only_search_trials': 'original infeasible objective; recount both legs in independent replays'}
        encoded = json.dumps(config, indent=2)
        if record.exists() and record.read_text() != encoded:
            parser.error('The existing study configuration differs from this invocation')
        record.write_text(encoded, encoding='utf-8')
        study_schedule(args)
    elif (args.worker is None) == (args.select is None):
        parser.error('A stage requires either a worker number or a finalist rank')
    else:
        study_worker(args)


if __name__ == '__main__':
    main()
