"""Bounded fresh-process pilots for the explicitly authorized Nov 2 campaign."""
import argparse
from collections import deque
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

from campaign_budget import require_time, VERIFICATION_STOP

ROOT = Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h').resolve()
SPINAL = Path(__file__).resolve().parents[1]
MYO = Path(r'D:\Anaconda\envs\myo\python.exe')
sys.path.insert(0, str(SPINAL / 'virtual_walker_20261031'))
from owned_process import OwnedProcessJob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('configuration', type=Path)
    parser.add_argument('--label', required=True)
    parser.add_argument('--max-wall-s', type=float, default=900.)
    args = parser.parse_args()
    config = args.configuration.resolve()
    if not config.is_relative_to(ROOT):
        raise ValueError('Only this campaign owns the supplied configuration')
    if not args.label.replace('_', '').replace('-', '').isalnum():
        raise ValueError('Use a simple unique run label')
    require_time(args.max_wall_s)
    folder = ROOT / 'timings' / args.label
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / 'RESULT.json').exists() or (folder / 'RUNNING.json').exists():
        raise FileExistsError('Existing timed run is preserved; use a new label')
    command = [str(MYO), str(SPINAL / 'hunt_training' / 'physical_eval.py'), str(config)]
    env = os.environ.copy()
    env.update({key: '1' for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
    env['PYTHONIOENCODING'] = 'utf-8'
    started = datetime.now(timezone.utc)
    tick = time.perf_counter()
    process = subprocess.Popen(command, cwd=SPINAL, env=env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, encoding='utf-8',
                               errors='replace', creationflags=subprocess.CREATE_NO_WINDOW)
    process_job = OwnedProcessJob(process)
    provenance = dict(configuration=str(config), configuration_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
                      command=command, pid=process.pid, start_utc=started.isoformat(),
                      max_wall_s=args.max_wall_s, verification_stop_utc=VERIFICATION_STOP.isoformat())
    (folder / 'RUNNING.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    tail = deque(maxlen=120)

    def consume():
        for line in process.stdout:
            tail.append(line.rstrip('\n')[:2000])

    reader = threading.Thread(target=consume, daemon=True)
    reader.start()
    timed_out = False
    try:
        try:
            code = process.wait(timeout=args.max_wall_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            process_job.close()
            code = process.wait(timeout=5.)
    finally:
        process_job.close()
    reader.join(timeout=2.)
    result = dict(**provenance, end_utc=datetime.now(timezone.utc).isoformat(),
                  total_process_wall_s=time.perf_counter()-tick, exit_code=code,
                  timed_out=timed_out, success=code == 0 and not timed_out,
                  timing_scope='Fresh interpreter startup, compilation, full body trial, assessment, and saved outputs. Plotting is separate.')
    (folder / 'PROCESS_TAIL.log').write_text('\n'.join(tail)+'\n', encoding='utf-8')
    (folder / 'RESULT.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
    raise SystemExit(0 if result['success'] else 1)


if __name__ == '__main__':
    main()
