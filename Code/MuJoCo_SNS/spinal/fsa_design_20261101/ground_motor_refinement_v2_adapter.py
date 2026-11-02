"""Opt-in v2 ground motor study; both Hunt and generic fall fields are checked."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
import ankle_refinement_adapter as B

A, I = B.A, B.I
VARIABLES = {'soleus_uS':(0.,.15), 'tib_ant_uS':(.3,.6), 'knee_flexor_multiplier':(1.,2.5)}


def finite_json(value):
    if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
    if isinstance(value,list):return [finite_json(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):return None
    return value


def set_parameters(net, params):
    if set(params) != set(VARIABLES): raise ValueError('Wrong decision variables')
    for key,value in params.items():
        if not math.isfinite(value) or not VARIABLES[key][0] <= value <= VARIABLES[key][1]:
            raise ValueError('Parameter outside registered bounds: '+key)
    before = {q:getattr(net,q).copy() for q in A.BUFFERS}
    edges,_ = A.set_pair(net,params['soleus_uS'],params['tib_ant_uS'])
    c=net.compiled
    for side in ('r','l'):
        source,target='PF_KNEE-F_'+side,'MN_bifemsh_'+side
        ids=[i for i,(a,b) in enumerate(zip(c.pre,c.post)) if net.names[int(a)]==source and net.names[int(b)]==target]
        if len(ids)!=1 or ids[0] not in net.pf_edges:raise ValueError('Missing exact knee-flexor PF output: '+source)
        i=ids[0];base=float(before['ff_gain'][i])
        if base<=0:raise ValueError('Source knee-flexor output is not positive')
        net.ff_gain[i]=base*params['knee_flexor_multiplier']
        edges['bifemsh_'+side]=dict(index=i,source=source,target=target,base_uS=base,requested_uS=float(net.ff_gain[i]))
    if edges['bifemsh_r']['base_uS']!=edges['bifemsh_l']['base_uS']:raise ValueError('Untied baseline knee outputs')
    allowed={e['index'] for e in edges.values()};changed={}
    for key in A.BUFFERS:
        ids=np.flatnonzero(getattr(net,key)!=before[key]).tolist()
        if (key!='ff_gain' and ids) or not set(ids)<=allowed:raise AssertionError('Non-target coefficient changed')
        changed[key]=ids
    return edges,changed


def measure(result, assistance):
    assessment=result['assessment'];bw=result['measured']['body_weight_N']
    if not assessment.get('finite') or not assessment.get('complete'):
        return dict(value=1.e6,unsafe=True,strict_accepted=False,
                    reason='Incomplete/nonfinite physical evaluation; no late-window metric invented')
    with np.load(assistance) as z:
        x=z['samples'];mask=(x[:,0]>=5)&(x[:,0]<13)
        if mask.sum()<3000:raise ValueError('Incomplete objective measurement window')
        spring=x[mask,1]/bw;prescribed=x[mask,2]/bw
        up=float(np.maximum(spring,0).mean()+np.maximum(prescribed,0).mean())
        down=float(np.maximum(-spring,0).mean()+np.maximum(-prescribed,0).mean())
        horizontal=float(np.sqrt(np.mean(np.sum(x[mask,9:11]**2,axis=1)))/bw)
        effort=float(x[mask,3].mean())
    clearance=result.get('foot_clearance',{});counts=[int(f['sustained_swings']) for f in clearance.get('feet',[])]
    if len(counts)!=2:raise ValueError('Both whole-mesh clearance measurements required')
    missing=sum(max(0,3-n) for n in counts)
    contacts=result.get('ground_contacts',{});phase=float(contacts.get('phase_error',.5))
    unsafe=not all(assessment.get(k) is True for k in ('finite','complete','frontal_gate','knee_extension_gate','ground_tilt_gate')) or bool(result['metrics'].get('fall',False)) or bool(result['metrics'].get('bal_fell',False))
    failures=[key for key in ('movement_gate','measured_contact_gate','ground_load_verified') if assessment.get(key) is not True]
    terms=dict(missing_clearance=1000.*missing,movement_contact=5000.*len(failures),phase=100.*phase,
               upward_assistance=10.*up,horizontal_force=2.*horizontal,activation_effort=effort)
    value=1.e6 if unsafe else sum(terms.values())
    if not np.isfinite(value):raise ValueError('Nonfinite objective')
    return dict(value=float(value),terms=terms,unsafe=unsafe,strict_accepted=(assessment.get('accepted') is True and not unsafe),
                clearance_swings=counts,missing_swings=missing,movement_contact_failures=failures,
                phase_error=phase,upward_support_bw=up,downward_support_bw=down,net_support_bw=up-down,
                horizontal_harness_rms_bw=horizontal,activation_squared=effort,window_s=[5.,13.],end_exclusive=True)


def preflight(m):
    if int(m['max_trials'])>11 or set(m['search_space'])!=set(VARIABLES) or m.get('fixed_parameters'):raise ValueError('Only the registered small three-pair study')
    for key,bounds in VARIABLES.items():
        if (m['search_space'][key]['low'],m['search_space'][key]['high'])!=bounds:raise ValueError('Changed study domain')
    for path,digest in m['source_sha256'].items():
        if I.sha(path)!=digest:raise ValueError('Source or selected evidence changed: '+path)
    campaign=B.read(B.ROOT/'CAMPAIGN_MANIFEST.json')
    for name in campaign['training_prerequisites']:
        item=campaign['prerequisites_passed'].get(name,{})
        if item.get('passed') is not True or not item.get('evidence'):raise ValueError('Unmet prerequisite: '+name)
        for path in item['evidence']:
            if not Path(path).is_file():raise FileNotFoundError(path)
    cfg=B.read(m['baseline_configuration']);p=I.load_profile(cfg['integrated_fsa_profile'])
    if p['local_profile']['contracts']!=m['local_targets']:
        raise ValueError('Original winner23 local targets changed')
    if cfg['stage']!=4 or cfg['support_fraction']!=.4 or cfg['balance']['kp']!=5. or cfg['balance']['kd']!=2.:
        raise ValueError('Wrong fixed ground/support/balance condition')
    if not cfg.get('suppress_inherited_posture_mn') or cfg.get('motor_teacher') or p['parameter_checkpoint'] is not None:
        raise ValueError('Expected original winner23 local rules and source motor baseline')
    if cfg.get('motor_scales') or cfg.get('motor_scale',1.)!=1. or cfg.get('reflex_scale',1.)!=1.:
        raise ValueError('Conflicting broad gain scales')
    if cfg['knee_training_limits_deg'][1]>0:raise ValueError('Knees unprotected')
    r=B.read(m['baseline_result'])
    if r['configuration']!=cfg or r['compiled']['integrated_fsa_design']['profile']!=p:raise ValueError('Baseline consumption mismatch')
    if not all(r['assessment'].get(k) is True for k in ('finite','complete','frontal_gate','knee_extension_gate','ground_tilt_gate')):
        raise ValueError('Unsafe baseline is not a fitting fixture')
    binding=B.read(m['mechanical_binding_proof'])
    if not binding['bilateral_bifemsh_is_knee_flexor'] or not binding['hip_moment_abs_max_m']<1e-8:
        raise ValueError('Unverified functional knee-flexor routing')
    return cfg,p,measure(r,Path(m['baseline_configuration']).with_suffix('.assistance.npz'))


def prepare(m,params,folder):
    cfg,p,_=preflight(m);acts=[x.attrib['name'] for x in A.ET.parse(I.MODEL).getroot().find('actuator')]
    net=I.C.build(acts,options=I.condition_options(cfg));c=net.compiled
    mask=[i for i in net.stand_edges if net.names[c.pre[i]]=='POSTURE' and net.names[c.post[i]].startswith('MN_')]
    net.original_g[mask]=0.;c.gmax[mask]=0.;net.sync();I.apply_to_network(net,p,protocol_stage=4)
    edges,changed=set_parameters(net,params)
    for module in (Path(__file__).resolve(),Path(A.__file__).resolve(),Path(B.__file__).resolve()):p['source_sha256'][str(module)]=I.sha(module)
    folder=B.owned(folder);p=I.snapshot_parameters(net,folder/'PARAMETERS.npz',p);B.write(folder/'integrated_profile.json',p)
    cfg['integrated_fsa_profile']=str(folder/'integrated_profile.json');cfg.pop('verification_capture',None);cfg.pop('export_graph',None)
    number=B.read(folder/'REQUEST.json')['trial_number'] if (folder/'REQUEST.json').exists() else -1
    if number in m['capture_trial_numbers']:cfg['verification_capture']=str(folder/'verification')
    else:cfg['diagnose_swing']=False
    B.write(folder/'stage4.json',cfg)
    proof=A.command_proof(net,edges)
    B.write(folder/'MOTOR_PROOF.json',dict(edges=edges,changed_authoritative_indices=changed,
          current_cache_command_proof=proof,posture_mask_count=len(mask),body_steps=0,
          local_targets=p['local_profile']['contracts'],checkpoint=p['parameter_checkpoint']))


def evaluate(request):
    m=request['manifest'];preflight(m)
    deadline=min(B.STOP,float(request['deadline_epoch_s']),time.time()+float(m['trial_timeout_s']))
    if deadline-time.time()<30:raise TimeoutError('Insufficient reserved trial time')
    name=request['evaluation_id']
    if not name.replace('_','').isalnum():raise ValueError('Unsafe evaluation name')
    folder=B.owned(m['output_dir'])/'trials'/name;folder.mkdir(parents=True,exist_ok=False);B.write(folder/'REQUEST.json',request)
    start=time.perf_counter()
    try:
        output=B.child([sys.executable,str(Path(__file__).resolve()),'--prepare',m['manifest_path'],
                        '--folder',str(folder),'--params',json.dumps(request['params'])],deadline)
        preparation=time.perf_counter()-start
        output+=B.child([sys.executable,str(I.H/'physical_eval.py'),str(folder/'stage4.json')],deadline)
    except subprocess.TimeoutExpired:
        d=dict(reason='owned_process_timeout',folder=str(folder),wall_s=time.perf_counter()-start)
        B.write(folder/'TIMEOUT.json',d);return dict(ok=True,value=1.e6,diagnostics=d)
    (folder/'PROCESS_TAIL.log').write_text('\n'.join(line[:2000] for line in output.splitlines()[-35:]),encoding='utf-8')
    result=B.read(folder/'stage4.result.json');scored=measure(result,folder/'stage4.assistance.npz')
    diagnostics=finite_json(dict(configuration=str(folder/'stage4.json'),result_file=str(folder/'stage4.result.json'),
                     assessment=result['assessment'],objective=scored,preparation_wall_s=preparation,total_wall_s=time.perf_counter()-start)
                     )
    B.write(folder/'TRIAL_SUMMARY.json',dict(params=request['params'],**diagnostics))
    if request['trial_number'] not in m['capture_trial_numbers'] and scored['strict_accepted']:
        with np.load(folder/'stage4.assistance.npz') as a:
            x=a['samples'];np.savez_compressed(folder/'LOAD_WINDOW.npz',samples=x[(x[:,0]>=5)&(x[:,0]<13)])
        for path in folder.glob('stage4*.npz'):path.unlink()
    return dict(ok=True,value=scored['value'],diagnostics=diagnostics)


def handle_request(request):
    if time.time()>=B.STOP:raise TimeoutError('Absolute training stop reached')
    if request['op']=='preflight':
        _,p,baseline=preflight(request['manifest'])
        return dict(ok=True,ready_for_optimization=True,baseline=baseline,local_targets=p['local_profile']['contracts'],
                    body_steps=0,decision_variables=VARIABLES,finalist_requires_fresh_S1_to_S4=True)
    if request['op']=='evaluate':return evaluate(request)
    raise ValueError('Unsupported operation')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',required=True)
    parser.add_argument('--folder',required=True);parser.add_argument('--params',required=True)
    args=parser.parse_args();prepare(B.read(args.prepare),json.loads(args.params),args.folder)

