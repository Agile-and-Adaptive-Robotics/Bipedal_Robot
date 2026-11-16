"""Read saved PF experiments; no controller, neural, or mechanical simulation runs."""
import ast
import hashlib
import json
import importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(r'D:\GitHub\Bipedal_Robot')
RAW = Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
OUT = Path(__file__).resolve().parent
PROFILE = ROOT / 'Code/MuJoCo_SNS/spinal/controller_profiles_20261104'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(path):
    return json.loads(path.read_text(encoding='utf-8'))

def inspect():
    for relative in ['syn6_PF_relays_20261105/combined/RESULT.json',
                     'syn6_native_20261105_v7/final_verified/TIMING.json',
                     'syn6_adaptive_20261105/replay/TIMING.json',
                     'syn6_adaptive_20261105/replay_live_range/TIMING.json',
                     'syn6_adaptive_20261105/live_full_native/TIMING.json',
                     'syn6_flexor_handover_20261105/predict150_live_timing/TIMING.json']:
        obj = load(RAW / relative)
        print(relative, 'KEYS', list(obj))
        if 'offline' in obj:
            print('HEAD', {k:v for k,v in obj.items() if k!='offline'})
            for factor, row in obj['offline'].items():
                print('FACTOR',factor,'KEYS', list(row),
                      'SUMMARY',{side:{k:v for k,v in val.items() if k!='channels'}
                                 for side,val in row['checks'].items()})
        elif 'closed_loop' in obj:
            print('HEAD', {k:v for k,v in obj.items() if k!='closed_loop'})
            print('CHECKS',{side:{k:v for k,v in val.items() if k!='channels'}
                           for side,val in obj['closed_loop'].items()})
        else:
            print('SCALARS',{k:v for k,v in obj.items() if not isinstance(v,(dict,list))})
    for relative in ['syn6_native_20261105_v7/final_verified/phase_1.0.npz',
                     'syn6_flexor_handover_20261105/predict150_stage2/capture/aligned_trace.npz']:
        with np.load(RAW/relative) as z:
            print(relative, {k:z[k].shape for k in z.files})

def extracted(path, names, namespace):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    selected = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name in names]
    for node in selected:
        node.decorator_list = []
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), namespace)

def array_repeat(first, second):
    with np.load(first) as a, np.load(second) as b:
        return dict(same_keys=set(a.files)==set(b.files),
                    equal={key:bool(np.array_equal(a[key],b[key],equal_nan=True)
                                   if np.issubdtype(a[key].dtype,np.number)
                                   else np.array_equal(a[key],b[key]))
                           for key in a.files if key in b.files})

def config_chain():
    path=RAW/'waveform_edge_D_20261104/syn6'
    chain=[]
    for _ in range(8):
        manifest=load(path/'MANIFEST.json')
        chain.append(str(path/'MANIFEST.json'))
        if 'configuration' in manifest:
            return manifest['configuration'], chain
        if 'fit' not in manifest:
            raise ValueError(manifest.keys())
        path=Path(manifest['fit'])
    raise ValueError('no configuration in manifest chain')

def audit():
    spec=importlib.util.spec_from_file_location('pf_audit_kine_ref',PROFILE.parent/'kine_ref.py')
    kine=importlib.util.module_from_spec(spec);spec.loader.exec_module(kine)
    offset=float(kine.load_reference()['lag_rl'])
    ns=dict(np=np,Path=Path,__file__=str(PROFILE/'native_pf_check.py'),left_phase_offset=lambda:offset)
    extracted(PROFILE/'pf_timing_audit.py',{'features','primitive_target'},ns)
    extracted(PROFILE/'syn6_s5_timing.py',{'measure'},ns)
    extracted(PROFILE/'native_pf_check.py',{'rises','assess'},ns)
    manifest,chain=config_chain()
    result=dict(scope='Saved-array statistics only; no neural/mechanical simulation, training, or optimization.',
                source_files={}, target_configuration=manifest,manifest_chain=chain,left_phase_offset=offset,
                offline=[],body_timing=[],body_results=[],repeat_checks=[])
    paths=[PROFILE/'pf_timing_audit.py',PROFILE/'syn6_s5_timing.py',PROFILE/'native_pf_check.py',
           PROFILE/'native_pf.py',PROFILE/'TRAINING_RECORD.md',PROFILE.parent/'fsa_results/fsa_backsolve.npz',
           kine.IK_MOT,kine.GRF_MOT]+[Path(p) for p in chain]
    for folder in ['syn6_native_20261105_v7/final_verified','syn6_adaptive_20261105/replay',
                   'syn6_adaptive_20261105/replay_live_range']:
        report=load(RAW/folder/'TIMING.json')
        paths.extend([RAW/folder/'TIMING.json',RAW/folder/'SPEC.json'])
        for factor,row in report['offline'].items():
            path=RAW/folder/f'phase_{factor}.npz'
            with np.load(path) as z:
                check=ns['assess'](z['time'],z['volts'],z['names'].tolist(),manifest,
                                   report['reference_rg_extensor_phase'])
                windows={side:list(zip(ns['rises'](z['time'],z['volts'][:,z['names'].tolist().index('RG_E_'+side)])[:-1],
                                       ns['rises'](z['time'],z['volts'][:,z['names'].tolist().index('RG_E_'+side)])[1:]))[-4:]
                         for side in ('l','r')}
            summary={side:{k:v for k,v in val.items() if k!='channels'} for side,val in check.items()}
            agrees=all(check[side]['passed']==row['checks'][side]['passed'] and
                       check[side]['checks']==row['checks'][side]['checks'] for side in ('l','r'))
            result['offline'].append(dict(folder=folder,factor=factor,summary=summary,windows=windows,
                                          saved_verdict_reproduced=agrees))
            paths.append(path)
    for folder,body in [('syn6_adaptive_20261105/live_full_native','syn6_adaptive_20261105/body_slow_handover'),
                        ('syn6_flexor_handover_20261105/predict150_live_timing','syn6_flexor_handover_20261105/predict150_stage2')]:
        report=load(RAW/folder/'TIMING.json')
        path=RAW/body/'capture/aligned_trace.npz'
        with np.load(path) as z:
            keep=z['t_neural']>=report['full_native_since_neural_s']
            check=ns['assess'](z['t_neural'][keep],z['neural_voltage'][keep],z['neural_names'].tolist(),
                               manifest,report['reference_rg_extensor_phase'])
        result['body_timing'].append(dict(folder=folder,body=body,
                   full_native_since_neural_s=report['full_native_since_neural_s'],
                   capture_duration_s=report['capture_duration_s'],
                   summary={side:{k:v for k,v in val.items() if k!='channels'} for side,val in check.items()},
                   saved_verdict_reproduced=all(check[s]['passed']==report['closed_loop'][s]['passed'] for s in ('l','r'))))
        paths.extend([path,RAW/folder/'TIMING.json',RAW/body/'capture/ALIGNED_TRACE_METADATA.json',RAW/body/'capture/model.mjb'])
    for body in ['syn6_flexor_handover_20261105/predict150_stage2','syn6_flexor_handover_20261105/predict150_stage3',
                 'syn6_flexor_handover_20261105/predict150_stage4']:
        report=load(RAW/body/'body.result.json')
        with np.load(RAW/body/'capture/aligned_trace.npz') as z:
            meta=load(RAW/body/'capture/ALIGNED_TRACE_METADATA.json')
            q=z['qpos'][:,z['joint_qpos_addresses']]*180/np.pi
            names=z['joint_names'].tolist()
            maximum_knee=[float(q[:,names.index('knee_angle_'+s)].max()) for s in ('r','l')]
            finite_states={key:bool(np.isfinite(z[key]).all()) for key in
                           ['qpos','qvel','activation','control','neural_voltage']}
            actual_hash=digest(RAW/body/'capture/model.mjb')
            compiled_hash_matches=actual_hash==meta['model_sha256']
        result['body_results'].append(dict(body=body,assessment=report['assessment'],maximum_knee_deg_RL=maximum_knee,
                     finite_saved_states=finite_states,
                     compiled_hash_matches=compiled_hash_matches,compiled_model_sha256=actual_hash,
                     wall_s=report['wall_s'],simulated_s=report.get('simulated_s')))
        paths.extend([RAW/body/'body.result.json',RAW/body/'body.json',RAW/body/'MOTOR_SHA256.txt',
                      RAW/body/'capture/ALIGNED_TRACE_METADATA.json',RAW/body/'capture/aligned_trace.npz',RAW/body/'capture/model.mjb'])
    for folder in ['syn6_native_20261105_v7','syn6_flexor_handover_20261105','syn6_s3_single_20261105']:
        for path in (RAW/folder).glob('*/STOP.json'):
            result.setdefault('failed_body_stops',[]).append(dict(path=str(path),stop=load(path)))
            paths.append(path)
    for first,second in [('syn6_native_20261105_v7/final_verified/phase_1.0.npz','syn6_native_20261105_v7/fresh_repeat/phase_1.0.npz'),
                         ('syn6_native_20261105_v7/final_verified/phase_0.85.npz','syn6_native_20261105_v7/fresh_repeat/phase_0.85.npz'),
                         ('syn6_native_20261105_v7/final_verified/phase_1.15.npz','syn6_native_20261105_v7/fresh_repeat/phase_1.15.npz'),
                         ('syn6_adaptive_20261105/body_slow_handover/capture/aligned_trace.npz','syn6_adaptive_20261105/body_slow_repeat/capture/aligned_trace.npz'),
                         ('syn6_flexor_handover_20261105/predict150_stage2/capture/aligned_trace.npz','syn6_flexor_handover_20261105/predict150_stage2_repeat/capture/aligned_trace.npz')]:
        result['repeat_checks'].append(dict(first=first,second=second,**array_repeat(RAW/first,RAW/second)))
        paths.extend([RAW/first,RAW/second])
    relay=load(RAW/'syn6_PF_relays_20261105/combined/RESULT.json')
    result['hybrid_relay_record']={k:v for k,v in relay.items() if k!='channels' and k!='designs'}
    result['hybrid_relay_record']['channel_summary']={k:{key:(sum(r['passed'] for r in val) if key=='mean' else val)
                                                      for key,val in rows.items() if key in ['mean','mn_weights_exact','rg_trace_exact']}
                                                    for k,rows in relay['channels'].items()}
    paths.extend([RAW/'syn6_PF_relays_20261105/combined/RESULT.json',RAW/'syn6_PF_relays_20261105/combined/FRESH_REPEAT.json'])
    paths.extend([RAW/'syn6_flexor_handover_20261105/SPEC_150ms.json',
                  RAW/'syn6_adaptive_20261105/SPEC_slow_handover.json'])
    paths.extend((RAW/'syn6_s3_single_20261105').glob('SPEC*.json'))
    result['source_files']={str(p):digest(p) for p in dict.fromkeys(paths) if p.exists()}
    (OUT/'PF_AUDIT_RESULTS.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('Wrote',OUT/'PF_AUDIT_RESULTS.json')
    print('Offline',[(r['factor'],r['summary'],r['saved_verdict_reproduced']) for r in result['offline']])
    print('Body timing',result['body_timing'])
    print('Body acceptance',[(r['body'],r['assessment'].get('accepted'),r['compiled_hash_matches']) for r in result['body_results']])
    print('Repeat equality',[(r['first'],all(r['equal'].values())) for r in result['repeat_checks']])

if __name__=='__main__':
    audit()
