"""Small mirrored Hunt fits through Rustuna, with fresh physical interpreters.

The dispatch worker persists. The mutable Hunt/runner modules deliberately do
not persist across physical trials. This avoids accumulating monkeypatches or
global neural parameters, while retaining the verified sparse neural backend.
"""
from collections import deque
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent
ROOT = Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h').resolve()
TRAINING_STOP = datetime(2026, 11, 2, 17, 43, 22, tzinfo=timezone.utc).timestamp()
VARIABLES = {
    'Hunt_B_timing': {'e_budget_nA', 'f_budget_nA'},
    'Hunt_D_local_reflex': {'ia_delta_mV', 'ii_delta_mV', 'ib_abs_delta_mV', 'reciprocal_abs_delta_mV'},
}


def owned(path):
    path = Path(path).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('Trial writes must stay in the owned nine-hour campaign')
    return path


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, value):
    owned(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def preflight(manifest):
    if int(manifest['max_trials']) > 24:
        raise ValueError('Only small runs of at most 24 trials are authorized here')
    gate = read(ROOT / 'CAMPAIGN_MANIFEST.json')
    for name in gate['training_prerequisites']:
        item = gate['prerequisites_passed'].get(name, {})
        if item.get('passed') is not True or not item.get('evidence'):
            raise ValueError('Unmet selected-campaign prerequisite: '+name)
        for path in item['evidence']:
            if not Path(path).is_file():
                raise FileNotFoundError(path)
    cfg = read(owned(manifest['baseline_configuration']))
    if cfg['stage'] != 2 or cfg.get('backend') != 'sparse':
        raise ValueError('This bounded adapter fits measured stage-2 motion with the sparse backend')
    if cfg.get('motor_teacher'):
        raise ValueError('A motor teacher is excluded from this neural fit')
    if cfg.get('knee_training_limits_deg', [-120, 0])[1] > 0:
        raise ValueError('Unprotected knees')
    profile = read(manifest['base_profile'])
    if cfg.get('integrated_fsa_profile') != manifest['base_profile']:
        raise ValueError('Baseline does not consume the exact selected profile')
    for path, digest in profile['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('Profile source changed; review and explicitly regenerate it: '+path)
    kind = manifest.get('study_kind', 'Hunt_B_timing')
    if kind not in VARIABLES or set(manifest['search_space']) != VARIABLES[kind] or manifest.get('fixed_parameters'):
        raise ValueError('Use the explicitly declared shared functional-target variables')
    if kind == 'Hunt_D_local_reflex':
        proof = read(owned(manifest['functional_target_api_proof']))
        if proof.get('local_targets_consumed') is not True or proof.get('other_gain_buffers_exact') is not True:
            raise ValueError('Local-target API needs consumed-target and frozen-other-buffer evidence')
    baseline = read(owned(manifest['baseline_result']))
    if baseline.get('configuration') != cfg:
        raise ValueError('Measured baseline configuration differs from the selected baseline')
    consumed = baseline.get('compiled', {}).get('integrated_fsa_design', {})
    if consumed.get('profile') != profile:
        raise ValueError('The measured baseline did not consume this exact integrated profile')
    if not baseline['assessment'].get('finite') or not baseline['assessment'].get('complete'):
        raise ValueError('A complete finite measured baseline is required before training')
    if baseline['assessment'].get('knee_extension_gate') is not True:
        raise ValueError('Baseline violated protected knee limits')
    return cfg, profile


def evaluate(request):
    manifest = request['manifest']
    cfg, profile = preflight(manifest)
    params = request['params']
    kind = manifest.get('study_kind', 'Hunt_B_timing')
    if set(params) != VARIABLES[kind]:
        raise ValueError('Unexpected fit variables')
    for key, value in params.items():
        bounds = manifest['search_space'][key]
        if not math.isfinite(value) or not bounds['low'] <= value <= bounds['high']:
            raise ValueError('Timing budget outside the declared study domain')
    deadline = min(TRAINING_STOP, float(request['deadline_epoch_s']))
    timeout = deadline-time.time()-3.
    if timeout <= 1.:
        raise TimeoutError('Insufficient time remains to launch a physical trial')
    evaluation_id = request['evaluation_id']
    if not evaluation_id.replace('_', '').isalnum():
        raise ValueError('Unsafe evaluation identifier')
    folder = owned(manifest['output_dir']) / 'trials' / evaluation_id
    folder.mkdir(parents=True, exist_ok=False)
    if kind == 'Hunt_B_timing':
        profile['timing_profile']['mode'] = 'custom_budgets'
        profile['timing_profile']['per_RG_budget_nA'] = dict(E=params['e_budget_nA'], F=params['f_budget_nA'])
    else:
        profile['local_profile']['contracts'] = dict(Ia_homonymous=params['ia_delta_mV'],
            II_homonymous=params['ii_delta_mV'], Ib_autogenic=-params['ib_abs_delta_mV'],
            Ia_reciprocal_pool=-params['reciprocal_abs_delta_mV'])
        # Re-derive the signed-pool conductances from their new targets. A saved
        # old reflex_gain checkpoint would otherwise overwrite these targets.
        # The API proof confirms that PF, feedforward, biarticular, and timing
        # parameters rebuild exactly while only these local targets change.
        profile['parameter_checkpoint'] = None
    profile_path = folder / 'integrated_profile.json'
    write(profile_path, profile)
    cfg['integrated_fsa_profile'] = str(profile_path)
    cfg.pop('verification_capture', None)
    if request['trial_number'] in manifest.get('capture_trial_numbers', []):
        cfg['verification_capture'] = str(folder/'verification')
    cfg.pop('export_graph', None)
    cfg_path = folder / 'stage2.json'
    write(cfg_path, cfg)
    env = os.environ.copy()
    env.update({key:'1' for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    env['PYTHONIOENCODING'] = 'utf-8'
    tick = time.perf_counter()
    try:
        child = subprocess.run([sys.executable, str(SPINAL/'hunt_training/physical_eval.py'), str(cfg_path)],
                               cwd=SPINAL, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, encoding='utf-8', errors='replace', timeout=timeout,
                               creationflags=subprocess.CREATE_NO_WINDOW)
    except subprocess.TimeoutExpired as exc:
        write(folder/'TIMEOUT.json', dict(timeout_s=timeout, deadline_epoch_s=deadline, wall_s=time.perf_counter()-tick))
        return dict(ok=True, value=1.e6, diagnostics=dict(reason='timed_out', configuration=str(cfg_path)))
    tail = deque(child.stdout.splitlines(), maxlen=60)
    (folder/'PROCESS_TAIL.log').write_text('\n'.join(line[:2000] for line in tail)+'\n', encoding='utf-8')
    if child.returncode:
        raise RuntimeError('Fresh physical trial failed; do not continue a structurally invalid sweep: '+str(folder))
    result_path = folder/'stage2.result.json'
    result = read(result_path)
    assessment = result['assessment']
    if not assessment.get('finite') or not assessment.get('complete'):
        value = 1.e5 + 100.*max(0., 16.-float(result['metrics'].get('t_end', 0.)))
    elif not assessment.get('knee_extension_gate', False):
        value = 1.e5
    else:
        value = -float(assessment['score']) + (0. if assessment.get('frontal_gate') else 100.)
    diagnostics = dict(result_file=str(result_path), profile_file=str(profile_path), configuration=str(cfg_path),
                       assessment=assessment, measured=result['measured'], fresh_process_wall_s=time.perf_counter()-tick,
                       objective='Negative existing continuous stage-2 motion score, plus 100 for frontal-gate failure; incomplete/nonfinite or hyperextension receives >=100000. No standing/walking claim follows from optimizer rank.')
    write(folder/'TRIAL_SUMMARY.json', dict(parameters=params, value=value, **diagnostics))
    return dict(ok=True, value=value, diagnostics=diagnostics)


def handle_request(request):
    if time.time() >= TRAINING_STOP:
        raise TimeoutError('Training reserve has started; no new fit evaluations')
    if request['op'] == 'preflight':
        cfg, profile = preflight(request['manifest'])
        return dict(ok=True, ready_for_optimization=True, stage=cfg['stage'],
                    study_kind=request['manifest'].get('study_kind', 'Hunt_B_timing'),
                    decision_variables=sorted(request['manifest']['search_space']), symmetry='Each declared functional target is shared across both sides; left and right are not fitted independently.',
                    dispatcher='Persistent myo dispatcher; every physical trial uses a fresh interpreter.',
                    physics='Exact baseline model, muscle force law, knee protection, and rig parameters are retained.')
    if request['op'] == 'evaluate':
        return evaluate(request)
    raise ValueError('Unsupported request operation')
