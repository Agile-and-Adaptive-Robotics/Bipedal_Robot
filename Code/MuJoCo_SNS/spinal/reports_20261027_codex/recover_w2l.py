"""Fresh-process, full-precision W2L recovery using the authoritative runner."""
from __future__ import annotations

import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent
OLD = SPINAL / 'reports_20261025_persyn'
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)


def run_child(job):
    path = HERE / (job['name'] + '_input.json')
    path.write_text(json.dumps(job, indent=2), encoding='utf-8')
    started = time.perf_counter()
    result = subprocess.run([sys.executable, str(Path(__file__)), '--one', str(path)],
                            capture_output=True, text=True, timeout=300)
    (HERE / (job['name'] + '.log')).write_text(result.stdout + result.stderr,
                                             encoding='utf-8')
    if result.returncode:
        raise RuntimeError(job['name'] + ': ' + result.stderr[-2000:])
    out = json.loads((HERE / (job['name'] + '.json')).read_text())
    print(job['name'], 'survive', out['metrics']['survive'],
          'cycles', out['metrics']['cycL'], out['metrics']['cycR'],
          'phase', out['metrics']['phase_med'],
          'wall', round(time.perf_counter() - started, 2), flush=True)
    return out


def one(path):
    import numpy as np
    sys.path.insert(0, str(SPINAL / 'reports_20261016/w2l_wean'))
    sys.path.insert(0, str(SPINAL / 'w2l_cpg'))
    import _w2l_wean_eval as EV
    import build_w2l_net as B
    job = json.loads(Path(path).read_text())
    os.environ['AARL_W2L_BACKEND'] = ('exact_fast' if job.get('backend') == 'exact_fast'
                                    else 'reference')
    EV.T.KV['net'] = 'bilateralrg2'
    EV.T.KV.update(job.get('text_knobs', {}))
    delivered = {}
    original = B.build

    def checked_build(*args, **kwargs):
        delivered.update(kwargs.get('gains', {}))
        net = original(*args, **kwargs)
        if job.get('backend') == 'sparse':
            from fast_eval import SparseForward
            net.compiled = SparseForward(net.compiled)
        if job.get('backend') == 'exact_fast':
            from exact_fast_eval import ExactFastForward
            if not isinstance(net.compiled, ExactFastForward):
                net.compiled = ExactFastForward(net.compiled)
        return net

    B.build = checked_build
    knobs = job['knobs']
    result = EV.evaluate(knobs, tag=job['name'],
                         csv_path=str(HERE / (job['name'] + '.csv')),
                         want_log=True)
    expected = {k: v for k, v in knobs.items() if k in B.W2L_GAINS}
    assert all(k in delivered and delivered[k] == v
               for k, v in expected.items()), 'Conductance forwarding mismatch'
    log = result.pop('_log')
    np.savez_compressed(HERE / (job['name'] + '.npz'), **log)
    output = dict(job=job, metrics=result, expected_gains=expected,
                  delivered_gains=delivered, python=sys.version,
                  numpy=np.__version__)
    (HERE / (job['name'] + '.json')).write_text(json.dumps(output, indent=2),
                                              encoding='utf-8')


def main():
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.load_study(study_name='curr_w2l_persyn_v2',
        storage='sqlite:///' + (OLD / 'optuna_w2l_persyn_v2.db').as_posix())
    top = sorted([t for t in study.trials if t.value is not None],
                 key=lambda t: t.value, reverse=True)[:8]
    jobs = []
    for trial in top:
        for dur in (20., 30.):
            knobs = dict(legcoll=0., sink=0.06, rig=0.25, walker_yield=1.5,
                         te=3., tf=4., wsup=0.7, sagittal=50.,
                         dur_walk=dur, **trial.params)
            jobs.append(dict(name=f'trial{trial.number}_d{int(dur)}',
                             knobs=knobs, original_value=trial.value,
                             original_metrics=trial.user_attrs))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(run_child, jobs))
    (HERE / 'recovery_results.json').write_text(json.dumps(results, indent=2),
                                               encoding='utf-8')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--one':
        one(sys.argv[2])
    else:
        main()
