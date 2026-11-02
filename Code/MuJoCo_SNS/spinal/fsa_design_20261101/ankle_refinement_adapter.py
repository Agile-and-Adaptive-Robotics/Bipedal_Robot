"""Bounded two-output Hunt refinement. No body trial runs during preflight.

The existing Rustuna dispatcher invokes evaluate only after explicit launch.
Every neural preparation and body evaluation has its own fresh interpreter.
"""
import argparse
import ast
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
import numpy as np
import ankle_motor_checkpoint as A

I = A.I
ROOT = Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h').resolve()
STOP = datetime(2026, 11, 2, 17, 43, 22, tzinfo=timezone.utc).timestamp()
sys.path.insert(0, str(I.S / 'virtual_walker_20261031'))
from owned_process import OwnedProcessJob
import kine_ref as K
VARIABLES = {'soleus_uS': (0., .15), 'tib_ant_uS': (.3, 1.2)}


def read(path): return json.loads(Path(path).read_text(encoding='utf-8'))


def owned(path):
    p = Path(path).resolve()
    if not p.is_relative_to(ROOT): raise ValueError('Output outside owned campaign')
    return p


def write(path, data):
    owned(path).write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')


def target():
    text = K.IK_MOT.read_text()
    if 'inDegrees=yes' not in text.split('endheader')[0]:
        raise ValueError('Reference must explicitly declare degree units')
    r = K.load_reference(); values = np.asarray(r['r']['ankle'])
    return dict(source=str(K.IK_MOT), source_sha256=I.sha(K.IK_MOT),
                force_reference=str(K.GRF_MOT), force_reference_sha256=I.sha(K.GRF_MOT),
                extraction_source=str(Path(K.__file__).resolve()), extraction_sha256=I.sha(K.__file__),
                canonical_side='right', units='degrees', native_ankle_sign=1.,
                mean_deg=float(values.mean()), range_deg=float(np.ptp(values)),
                cycle_period_s=float(r['T_r']), samples=len(values),
                provenance='First complete right GRF-loading cycle fully inside subject01 IK support; original 100-point kine_ref cycle interpolation. Same right target on both sides; no independently fitted left target.')


def ankle_metrics(npz, reference):
    tree = ast.parse((I.S/'runner.py').read_text())
    assignment = next(n for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'KEY_JOINTS' for t in n.targets))
    names = ast.literal_eval(assignment.value)
    with np.load(npz) as a:
        t = a['t']; window = (t >= 5.) & (t < 13.)
        if window.sum() < 3000 or t[window][-1] < 12.99:
            raise ValueError('Complete predeclared late full-command window is unavailable')
        if not np.all(a['neuro'][window, 0] == 2.5):
            raise ValueError('Ankle window contains a ramp/standing command')
        q = a['q'][window]
        out = {}
        for side in ('r', 'l'):
            values = q[:, names.index('ankle_angle_'+side)]
            if not np.isfinite(values).all(): raise ValueError('Nonfinite ankle coordinate')
            mean, span = float(values.mean()), float(np.ptp(values))
            out[side] = dict(mean_deg=mean, range_deg=span, min_deg=float(values.min()), max_deg=float(values.max()),
                             mean_error_deg=mean-reference['mean_deg'], range_error_deg=span-reference['range_deg'])
        loss = float(np.mean([(row[key]/reference['range_deg'])**2
                              for row in out.values() for key in ('mean_error_deg', 'range_error_deg')]))
        return dict(window_s=[5.,13.], end_exclusive=True, verified_drive_nA=2.5,
                    recorded_units='runner q is degrees', sides=out, normalized_squared_error=loss)


def preflight(m):
    campaign = read(ROOT/'CAMPAIGN_MANIFEST.json')
    for name in campaign['training_prerequisites']:
        item = campaign['prerequisites_passed'].get(name,{})
        if item.get('passed') is not True or not item.get('evidence'):
            raise ValueError('Unmet selected campaign prerequisite: '+name)
        for path in item['evidence']:
            if not Path(path).is_file(): raise FileNotFoundError(path)
    if int(m['max_trials']) > 16 or set(m['search_space']) != set(VARIABLES) or m.get('fixed_parameters'):
        raise ValueError('Only this bounded two-parameter study is supported')
    for key, bounds in VARIABLES.items():
        if (m['search_space'][key]['low'], m['search_space'][key]['high']) != bounds:
            raise ValueError('Changed endpoint domain')
    for path, digest in m['source_sha256'].items():
        if I.sha(path) != digest: raise ValueError('Source changed: '+path)
    cfg = read(m['baseline_configuration']); origin = read(m['posture_configuration'])
    for key in set(cfg)|set(origin):
        if key not in ('stage','verification_capture','diagnostic_scope') and cfg.get(key) != origin.get(key):
            raise ValueError('Posture/neural/physics option mismatch: '+key)
    if cfg['stage'] != 2 or cfg.get('motor_teacher') or not cfg.get('suppress_inherited_posture_mn'):
        raise ValueError('Wrong stage or posture policy')
    if cfg.get('motor_scale',1.) != 1. or cfg.get('reflex_scale',1.) != 1. or cfg.get('motor_scales'):
        raise ValueError('Conflicting broad gain scales')
    p = I.load_profile(cfg['integrated_fsa_profile'])
    b = read(m['baseline_result'])
    if b['configuration'] != cfg or b['compiled']['integrated_fsa_design']['profile'] != p:
        raise ValueError('Baseline did not consume this exact selected recipe')
    gates = ('accepted','finite','complete','knee_extension_gate','frontal_gate')
    if not all(b['assessment'].get(k) is True for k in gates):
        raise ValueError('Selected full inherited S2 baseline has not passed: '+str({k:b['assessment'].get(k) for k in gates}))
    t = target()
    if t != m['reference']: raise ValueError('Reference extraction changed')
    return cfg, p, b, ankle_metrics(Path(m['baseline_configuration']).with_suffix('.npz'),t)


def prepare(m, params, folder):
    cfg, p, _, _ = preflight(m)
    if set(params) != set(VARIABLES): raise ValueError('Unexpected motor parameter')
    acts = [x.attrib['name'] for x in A.ET.parse(I.MODEL).getroot().find('actuator')]
    n = I.C.build(acts, options=I.condition_options(cfg))
    ids = [i for i in n.stand_edges if n.names[n.compiled.pre[i]] == 'POSTURE' and n.names[n.compiled.post[i]].startswith('MN_')]
    n.original_g[ids] = 0.; n.compiled.gmax[ids] = 0.; n.sync()
    I.apply_to_network(n, p, protocol_stage=2)
    edges, changed = A.set_pair(n, params['soleus_uS'], params['tib_ant_uS'])
    p['source_sha256'][str(Path(__file__).resolve())] = I.sha(__file__)
    p['source_sha256'][str(Path(A.__file__).resolve())] = I.sha(A.__file__)
    folder = owned(folder)
    p = I.snapshot_parameters(n, folder/'PARAMETERS.npz', p)
    write(folder/'integrated_profile.json', p)
    number = int(read(folder/'REQUEST.json')['trial_number']) if (folder/'REQUEST.json').exists() else -1
    cfg['integrated_fsa_profile'] = str(folder/'integrated_profile.json')
    cfg.pop('verification_capture',None); cfg.pop('export_graph',None)
    if number in m['capture_trial_numbers']: cfg['verification_capture'] = str(folder/'verification')
    else: cfg['diagnose_swing'] = False
    write(folder/'stage2.json',cfg)
    proof = A.command_proof(n,edges)
    write(folder/'MOTOR_CHECKPOINT_PROOF.json',dict(edges=edges,changed_authoritative_indices=changed,
          command_proof=proof,local_targets_installed_first=p['local_profile']['contracts'],
          posture_MN_mask_count=len(ids),body_steps=0,checkpoint=p['parameter_checkpoint']))


def child(command, deadline):
    env=os.environ.copy();env['PYTHONIOENCODING']='utf-8'
    for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    p=subprocess.Popen(command,cwd=I.S,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                       text=True,encoding='utf-8',errors='replace',creationflags=subprocess.CREATE_NO_WINDOW)
    job=OwnedProcessJob(p)
    try:
        output,_=p.communicate(timeout=max(.1,deadline-time.time()))
        if p.returncode:raise RuntimeError('Fresh child failed: '+output[-6000:])
        return output
    finally:
        job.close()
        if p.poll() is None:p.wait(timeout=5.)


def evaluate(request):
    m=request['manifest'];preflight(m)
    deadline=min(STOP,float(request['deadline_epoch_s']),time.time()+float(m['trial_timeout_s']))
    if deadline-time.time()<30:raise TimeoutError('Insufficient reserved time for a new trial')
    name=request['evaluation_id']
    if not name.replace('_','').isalnum():raise ValueError('Unsafe evaluation ID')
    folder=owned(m['output_dir'])/'trials'/name;folder.mkdir(parents=True,exist_ok=False)
    write(folder/'REQUEST.json',request)
    start=time.perf_counter()
    try:
        output=child([sys.executable,str(Path(__file__).resolve()),'--prepare',m['manifest_path'],
                      '--folder',str(folder),'--params',json.dumps(request['params'])],deadline)
        preparation_s=time.perf_counter()-start
        output+=child([sys.executable,str(I.H/'physical_eval.py'),str(folder/'stage2.json')],deadline)
    except subprocess.TimeoutExpired:
        diagnostics=dict(reason='owned_child_timeout',deadline_epoch_s=deadline,
                         total_wall_s=time.perf_counter()-start,folder=str(folder))
        write(folder/'TIMEOUT.json',diagnostics)
        return dict(ok=True,value=1.e6,diagnostics=diagnostics)
    (folder/'PROCESS_TAIL.log').write_text('\n'.join(deque(output.splitlines(),maxlen=35)),encoding='utf-8')
    result=read(folder/'stage2.result.json');a=result['assessment']
    strict=all(a.get(k) is True for k in ('accepted','finite','complete','knee_extension_gate','frontal_gate'))
    metrics=ankle_metrics(folder/'stage2.npz',m['reference']) if a.get('complete') and a.get('finite') else None
    value=metrics['normalized_squared_error'] if strict and metrics else 1.e5
    diagnostics=dict(configuration=str(folder/'stage2.json'),result_file=str(folder/'stage2.result.json'),
                     strict_S2_pass=strict,assessment=a,ankle=metrics,measured=result['measured'],
                     preparation_wall_s=preparation_s,total_wall_s=time.perf_counter()-start,
                     objective='Minimize mean of squared ankle mean and full-range errors, each normalized by actual canonical-right reference RoM, across both sides; strict inherited S2/knee/frontal failure receives100000. Air only, no ground-walking claim.')
    write(folder/'TRIAL_SUMMARY.json',dict(params=request['params'],value=value,**diagnostics))
    if request['trial_number'] not in m['capture_trial_numbers'] and strict:
        # Preserve compact actual ankle/drive trace; keep all recipes and result/gate evidence.
        with np.load(folder/'stage2.npz') as z:
            window=(z['t']>=5.)&(z['t']<13.)
            np.savez_compressed(folder/'ANKLE_WINDOW.npz',t=z['t'][window],q_degrees=z['q'][window][:,[5,10]],drive_nA=z['neuro'][window,0])
        for path in folder.glob('stage2*.npz'):path.unlink()
    return dict(ok=True,value=float(value),diagnostics=diagnostics)


def handle_request(request):
    if time.time()>=STOP:raise TimeoutError('Absolute training stop reached')
    if request['op']=='preflight':
        _,p,_,metrics=preflight(request['manifest'])
        return dict(ok=True,ready_for_optimization=True,body_steps=0,baseline_ankle=metrics,
                    local_targets=p['local_profile']['contracts'],decision_variables=VARIABLES,
                    scope='Prepared two-output study; preflight performs no body trial or fit.')
    if request['op']=='evaluate':return evaluate(request)
    raise ValueError('Unsupported operation')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',required=True)
    parser.add_argument('--folder',required=True);parser.add_argument('--params',required=True)
    a=parser.parse_args();prepare(read(a.prepare),json.loads(a.params),a.folder)
