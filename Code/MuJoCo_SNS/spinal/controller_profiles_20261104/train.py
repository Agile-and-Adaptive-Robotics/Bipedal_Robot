"""Bounded, explicit Rustuna launch for prepared offline Hunt A-D blocks.

Use D:\Anaconda\python.exe (3.11 optimizer); scientific worker remains myo3.10.
Without --execute --computer-free this prints a plan and starts nothing.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
S = HERE.parent
SPACES = dict(
    A=dict(tau_scale=(.5, 3.), g_ion_scale=(.7, 1.3)),
    B={k + '_RG': (.05, 128.) for k in ('Ia', 'II', 'Ib')},
    C={'ff_' + k: (.2, 3.) for k in ('hip_ext', 'hip_flex', 'knee_ext', 'knee_flex', 'ankle_pf', 'ankle_df')},
    D={'local_' + k: (.2, 2.) for k in ('Ia', 'II', 'Ib')},
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', choices=('hunt', 's3k', 'w2lvar', 'syn6'), required=True)
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--phase', choices=tuple(SPACES), required=True)
    p.add_argument('--trials', type=int, default=16)
    p.add_argument('--budget-s', type=float, default=900.)
    p.add_argument('--execute', action='store_true')
    p.add_argument('--computer-free', action='store_true')
    p.add_argument('--runtime', type=Path, help='Existing isolated Python-3.11 Rustuna package directory')
    a = p.parse_args()
    space=dict(SPACES[a.phase])
    if a.phase=='A' and a.family in ('w2lvar','syn6'):
        space.update(commissural_F=(.001,1.),commissural_E=(.001,1.))
    if not 1 <= a.trials <= 64 or not 0 < a.budget_s <= 3600:
        p.error('Use bounded blocks: 1..64 trials and at most one hour')
    print(json.dumps(dict(family=a.family, phase=a.phase, trials=a.trials,
                         wall_budget_s=a.budget_s, execution_requested=a.execute and a.computer_free,
                         scope='Offline neural A-D fitting; corrected physical curriculum follows separately'), indent=2))
    if not (a.execute and a.computer_free):
        return
    static = json.loads((a.prepared / 'STATIC_CHECK.json').read_text())
    if not static['static_pass']:
        raise ValueError('Prepared static gate failed')
    hashes = json.loads((a.prepared / 'SOURCE_HASHES.json').read_text())
    for path, expected in hashes.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError('Prepared source changed; repeat static preparation')
    profile = json.loads((a.prepared / 'PROFILE.json').read_text())
    fixed = {}
    phase_index = tuple(SPACES).index(a.phase)
    for phase in tuple(SPACES)[:phase_index]:
        final = json.loads((a.prepared / 'training' / phase / 'FINAL.json').read_text())
        if final['accepted'] is not True:
            raise ValueError('Earlier Hunt block did not pass: ' + phase)
        fixed.update(final['parameters'])
    folder = a.prepared / 'training' / a.phase
    if folder.exists():
        raise FileExistsError('Keep the completed/interrupted study; choose a separate prepared candidate for a new run')
    folder.mkdir(parents=True)
    target = S / 'hunt_training/runs/targets/targets.npz'
    options = json.loads((a.prepared / 'HUNT_OPTIONS.json').read_text()) if (a.prepared / 'HUNT_OPTIONS.json').exists() else None
    cfg = dict(family=a.family, phase=a.phase, targets=str(target),
               targets_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
               target_start_s=.617, target_period_s=1.233,
               command_F=.15, command_E=.15,
               controller_profile=profile, hunt_options=options, fixed_parameters=fixed)
    cfg['source_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (HERE / 'profile.py', HERE / 'replay_adapter.py', Path(__file__))}
    (folder / 'MANIFEST.json').write_text(json.dumps(cfg, indent=2))
    sys.path.insert(0, str(S / 'virtual_walker_20261031'))
    from rustuna_transfer import Worker
    runtime = a.runtime or a.prepared.parent / 'rustuna_runtime_py311'
    if sys.version_info[:2] != (3, 11):
        raise RuntimeError('Use D:\\Anaconda\\python.exe (Python 3.11) for the isolated Rustuna wheel; myo remains the worker')
    if not runtime.is_dir():
        raise FileNotFoundError('Persistent isolated Rustuna runtime missing: ' + str(runtime))
    sys.path.insert(0, str(runtime))
    import rustuna as backend
    study = backend.create_study(direction='minimize', sampler=backend.samplers.TPESampler(seed=61104, n_startup_trials=min(8, a.trials), multivariate=False))
    study.enqueue_trial({key: 1. for key in space})
    if a.phase=='A' and a.family in ('w2lvar','syn6'):
        for f,e in ((.03,.001),(.1,.001),(.3,.001),(.03,.03)):
            study.enqueue_trial(dict(tau_scale=1.,g_ion_scale=1.,commissural_F=f,commissural_E=e))
    if a.phase == 'B':
        # Separate ablation probes entrained Hunt at these budgets. Preserve
        # that evidence as initial candidates, not as predetermined winners.
        for values in ((96., .05, .05), (.05, 96., .05), (32., 32., 32.)):
            study.enqueue_trial(dict(zip(SPACES['B'], values)))
    worker = Worker(HERE / 'replay_adapter.py', folder)
    deadline = time.time() + a.budget_s
    best = None
    try:
        ready = worker.request(dict(op='preflight', manifest=cfg, deadline_epoch_s=min(deadline, time.time() + 180.)), min(180., a.budget_s))
        (folder / 'PREFLIGHT.json').write_text(json.dumps(ready, indent=2))
        if not ready.get('ready_for_optimization'):
            raise RuntimeError('Neural replay preflight failed')
        for i in range(a.trials):
            if time.time() >= deadline:
                break
            trial = study.ask()
            params = {**fixed, **{k: trial.suggest_float(k, *bounds, log=a.phase == 'B'or k.startswith('commissural_')) for k, bounds in space.items()}}
            result = worker.request(dict(op='evaluate', params=params, deadline_epoch_s=min(deadline, time.time() + 180.)), min(180., deadline - time.time()))
            if 'value' not in result:
                raise RuntimeError('Scientific worker rejected trial: ' + json.dumps(result))
            study.tell(trial.number, float(result['value']))
            row = dict(trial=i, parameters=params, **result)
            with (folder / 'TRIALS.jsonl').open('a', encoding='utf-8') as stream:
                stream.write(json.dumps(row, allow_nan=False) + '\n')
            if best is None or (row['diagnostics']['accepted'] and not best['diagnostics']['accepted']) or (row['diagnostics']['accepted'] == best['diagnostics']['accepted'] and row['value'] < best['value']):
                best = row
                (folder / 'BEST.json').write_text(json.dumps(best, indent=2))
    finally:
        worker.close()
    if best is not None:
        # New process, same source-locked manifest, no state carried from TPE.
        fresh = Worker(HERE / 'replay_adapter.py', folder)
        try:
            fresh.request(dict(op='preflight', manifest=cfg, deadline_epoch_s=time.time()+180.),180.)
            repeat=fresh.request(dict(op='capture',params=best['parameters'],path=str(folder/'BEST_TRACE.npz'),deadline_epoch_s=time.time()+180.),180.)
        finally:
            fresh.close()
        delta=abs(float(repeat['value'])-float(best['value']))
        matched=delta<1e-10
        (folder/'FRESH_REPEAT.json').write_text(json.dumps(dict(matched=matched,objective_abs_difference=delta,**repeat),indent=2))
        (folder / 'FINAL.json').write_text(json.dumps(dict(accepted=bool(best['diagnostics']['accepted'] and repeat['diagnostics']['accepted'] and matched), parameters=best['parameters'], best=best,
            fresh_process_matched=matched,scope='Offline neural block only; body curriculum not yet passed'), indent=2))


if __name__ == '__main__':
    main()
