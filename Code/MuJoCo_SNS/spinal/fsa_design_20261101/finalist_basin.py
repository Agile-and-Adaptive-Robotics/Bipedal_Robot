"""Seeded, mirrored finalist probes. No optimizer; explicit run command required.

Primary scope is selected10, not all coefficients or whole-model robustness.
The optional matrix probe is electrical-only until separately reviewed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import sys
import time
import numpy as np
import ankle_motor_checkpoint as A
import ground_motor_refinement_v2_adapter as G
import campaign_budget

I,B=A.I,G.B
SELECTED=('soleus_uS','tib_ant_uS','knee_flexor_multiplier','Ia_homonymous',
          'II_homonymous','Ib_autogenic','Ia_reciprocal_pool','rg_tonic_current','PT_kp','PT_kd')
DEADLINE=campaign_budget.VERIFICATION_STOP.timestamp()


def factors(seed):
    if seed is None:return {key:1. for key in SELECTED}
    rng=np.random.default_rng(int(seed))
    return {key:float(rng.uniform(.99,1.01)) for key in SELECTED}


def checkpoint_values(net,profile):
    item=profile['parameter_checkpoint']
    if item is None:raise ValueError('Finalist requires its exact fitted gain checkpoint')
    with np.load(item['path']) as a:
        names=a['names'].tolist();lookup={(names[int(x)],names[int(y)]):i for i,(x,y) in enumerate(zip(a['pre'],a['post']))}
        if len(lookup)!=len(a['pre']):raise ValueError('Duplicate endpoint checkpoint')
        take=[lookup[(net.names[int(x)],net.names[int(y)])] for x,y in zip(net.compiled.pre,net.compiled.post)]
        return {q:a[q][take].copy() for q in A.BUFFERS}


def matrix_probe(net,seed,excluded):
    """Qualified fitted pairs only; geometry/unpaired/asymmetric cases disclosed."""
    source_files={'pf_gain':'pf_weights.json','ff_gain':'mn_weights.json','bi_gain':'bi_weights.json'}
    masks={'pf_gain':np.union1d(net.pfbank_edges,net.pf_internal_edges),
           'ff_gain':np.union1d(net.pf_edges,net.geom_edges),'bi_gain':net.bi_edges}
    c=net.compiled;lookup={(net.names[int(a)],net.names[int(b)]):i for i,(a,b) in enumerate(zip(c.pre,c.post))}
    explicit={}
    for a,b in net.midline_drive_pairs:explicit[int(a)]=int(b);explicit[int(b)]=int(a)
    mirror=lambda name:re.sub(r'_([rl])(?=$|[^a-zA-Z0-9])',lambda m:'_l' if m[1]=='r' else '_r',name)
    rng=np.random.default_rng((0 if seed is None else int(seed))+1000003);ledger=[]
    for field,file in source_files.items():
        fitted=B.read(Path(net.options['weights_directory'])/file)['edges'];selected=set(map(int,masks[field]));seen=set()
        for i in sorted(selected):
            if i in seen:continue
            pre,post=net.names[c.pre[i]],net.names[c.post[i]];edge=pre+' -> '+post
            if edge not in fitted or getattr(net,field)[i]==0:continue
            mapped=(mirror(pre),mirror(post));j=explicit.get(i,lookup.get(mapped));row=dict(buffer=field,source=pre,target=post,index=i)
            if i in excluded:row['decision']='motor-owned; selected10 handles this';seen.add(i)
            elif i in set(net.geom_edges.tolist()):row['decision']='frozen: anatomical geometry reflection not verified';seen.add(i)
            elif j is None or j not in selected:row['decision']='frozen: unpaired fitted nonzero';seen.add(i)
            elif (net.names[c.pre[j]]+' -> '+net.names[c.post[j]]) not in fitted:
                row['decision']='frozen: counterpart not in fitted file';seen.add(i)
            elif getattr(net,field)[i]!=getattr(net,field)[j]:row['decision']='frozen: unequal source coefficients';seen.update((i,j))
            else:
                # Right postsynaptic label is canonical; shared unsided slots use one factor.
                if post.endswith('_l') or (not post.endswith('_r') and pre.endswith('_l')):i,j=j,i
                factor=1. if seed is None else float(rng.uniform(.99,1.01));base=float(getattr(net,field)[i]);getattr(net,field)[i]=base*factor;getattr(net,field)[j]=base*factor
                row.update(decision='perturbed qualified pair',canonical_source=net.names[c.pre[i]],canonical_target=net.names[c.post[i]],
                           mirror_source=net.names[c.pre[j]],mirror_target=net.names[c.post[j]],factor=factor,base_uS=base)
                seen.update((i,j))
            ledger.append(row)
    return ledger


def make_case(base_config,motor_proof,out,seed=None,forced=None,mode='selected10'):
    base_config=Path(base_config).resolve();cfg=B.read(base_config);original=I.load_profile(cfg['integrated_fsa_profile'])
    if cfg['stage']!=4:raise ValueError('Fresh S4 basin only')
    if cfg.get('motor_scale',1.)!=1 or cfg.get('reflex_scale',1.)!=1 or cfg.get('motor_scales'):raise ValueError('Ambiguous broad gain scales')
    if not cfg.get('suppress_inherited_posture_mn'):raise ValueError('Expected selected masked-posture architecture')
    motor_proof=Path(motor_proof).resolve();proof=B.read(motor_proof);motor={}
    for muscle,key in [('soleus','soleus_uS'),('tib_ant','tib_ant_uS')]:
        r,l=(proof['edges'][muscle+'_'+s]['requested_uS'] for s in ('r','l'))
        if r!=l:raise ValueError('Untied motor source')
        motor[key]=r
    # Read the actual checkpoint and inherited source coefficient, also allowing
    # older two-ankle-output ledgers that predate the separate knee multiplier.
    learned=B.read(Path(cfg['weights_directory'])/'mn_weights.json')['edges']
    with np.load(original['parameter_checkpoint']['path']) as saved:
        names=saved['names'].tolist();edge_lookup={(names[int(a)],names[int(b)]):i for i,(a,b) in enumerate(zip(saved['pre'],saved['post']))}
        knee={s:dict(base_uS=learned[f'PF_KNEE-F_{s} -> MN_bifemsh_{s}'],requested_uS=float(saved['ff_gain'][edge_lookup[(f'PF_KNEE-F_{s}',f'MN_bifemsh_{s}')]])) for s in ('r','l')}
    kr,kl=(knee[s] for s in ('r','l'))
    if kr['base_uS']!=kl['base_uS'] or kr['requested_uS']!=kl['requested_uS']:raise ValueError('Untied knee source')
    motor['knee_flexor_multiplier']=kr['requested_uS']/kr['base_uS']
    values={**motor,**original['local_profile']['contracts'],'rg_tonic_current':cfg['rg_tonic_current'],
            'PT_kp':cfg['balance']['kp'],'PT_kd':cfg['balance']['kd']}
    scale=factors(seed) if forced is None else dict(forced)
    if set(scale)!=set(SELECTED) or any(not .99<=v<=1.01 for v in scale.values()):raise ValueError('Exactly ten positive ±1% factors required')
    perturbed={key:float(values[key]*scale[key]) for key in SELECTED}
    if any(values[k]==0 and perturbed[k]!=0 for k in SELECTED):raise AssertionError('Zero reactivated')
    p=copy.deepcopy(original);p['parameter_checkpoint']=None
    p['local_profile']['contracts']={key:perturbed[key] for key in original['local_profile']['contracts']}
    cfg['rg_tonic_current']=perturbed['rg_tonic_current'];cfg['balance']['kp']=perturbed['PT_kp'];cfg['balance']['kd']=perturbed['PT_kd']
    acts=[x.attrib['name'] for x in A.ET.parse(I.MODEL).getroot().find('actuator')]
    n=I.C.build(acts,options=I.condition_options(cfg));c=n.compiled
    mask=[i for i in n.stand_edges if n.names[c.pre[i]]=='POSTURE' and n.names[c.post[i]].startswith('MN_')]
    n.original_g[mask]=0.;c.gmax[mask]=0.;n.sync();I.apply_to_network(n,p,protocol_stage=4)
    old=checkpoint_values(n,original);protected=set(map(int,n.hunt_graded_fsa['protected_indices']))
    for field in A.BUFFERS:
        if field=='reflex_gain':
            ids=np.array([i for i in range(len(c.pre)) if i not in protected],dtype=int);n.reflex_gain[ids]=old[field][ids]
        else:getattr(n,field)[:]=old[field]
    # No fit-domain clipping: a boundary value may move1% beyond optimizer bounds.
    edges=A.endpoints(n)
    for name,row in edges.items():
        row['requested_uS']=perturbed['soleus_uS' if name.startswith('soleus') else 'tib_ant_uS'];n.ff_gain[row['index']]=row['requested_uS']
    for side in ('r','l'):
        source,target='PF_KNEE-F_'+side,'MN_bifemsh_'+side
        ids=[i for i,(a,b) in enumerate(zip(c.pre,c.post)) if n.names[a]==source and n.names[b]==target]
        if len(ids)!=1 or ids[0] not in n.pf_edges:raise ValueError('Missing knee output')
        i=ids[0];value=float(kr['base_uS']*perturbed['knee_flexor_multiplier']);n.ff_gain[i]=value
        edges['bifemsh_'+side]=dict(index=i,source=source,target=target,requested_uS=value)
    matrix_ledger=[]
    if mode=='matrix_probe':matrix_ledger=matrix_probe(n,seed,{x['index'] for x in edges.values()})
    elif mode!='selected10':raise ValueError('Unknown scope')
    rc_sources={};rc_proof={'applied':False}
    if cfg.get('renshaw_autogenic_profile'):
        import renshaw_autogenic_design as RC
        rc_path=Path(cfg['renshaw_autogenic_profile']).resolve()
        rc_proof=RC.install(n,RC.load_profile(rc_path))
        rc_sources={str(rc_path):I.sha(rc_path),str(Path(RC.__file__).resolve()):I.sha(RC.__file__)}
    with I.unbound_electrical(n):n.apply_weights(1.);n.sync();pairing=I.IB.audit_hunt_gain_pairing(n);I.verify(n,full=True)
    if any(q['mismatched_pairs'] for q in pairing['buffers'].values()):raise ValueError('Pairing changed')
    out=B.owned(out);out.mkdir(parents=True,exist_ok=False)
    cfg['integrated_fsa_profile']=str(out/'PROFILE.json');cfg.pop('verification_capture',None);cfg.pop('export_graph',None)
    cfg['diagnose_swing']=False
    if seed is None:cfg['verification_capture']=str(out/'verification')
    B.write(out/'CONFIGURATION.json',cfg);p['condition_file']=str(out/'CONFIGURATION.json')
    p['source_sha256'][str(Path(__file__).resolve())]=I.sha(__file__);p['source_sha256'][str(out/'CONFIGURATION.json')]=I.sha(out/'CONFIGURATION.json')
    p['source_sha256'][str(Path(campaign_budget.__file__).resolve())]=I.sha(campaign_budget.__file__)
    for dependency in (I.H/'physical_eval.py',I.S/'runner.py',I.H/'controller.py',I.H/'balance.py',I.MODEL,Path(cfg['force_law_parameters_file'])):
        p['source_sha256'][str(dependency.resolve())]=I.sha(dependency)
    p['source_sha256'].update(rc_sources)
    p=I.snapshot_parameters(n,out/'PARAMETERS.npz',p);B.write(out/'PROFILE.json',p)
    command=A.command_proof(n,edges)
    changed={field:int(np.count_nonzero(getattr(n,field)!=old[field])) for field in A.BUFFERS}
    base_topology=hashlib.sha256(np.c_[c.pre,c.post].astype(np.int64).tobytes()).hexdigest()
    if seed is None and forced is None and any(changed.values()):raise AssertionError('Default is not exact')
    record=dict(seed=seed,mode=mode,factors=scale,base_values=values,perturbed_values=perturbed,
                changed_buffer_counts=changed,local_contracts_rederived=True,protected_reflex_indices=len(protected),
                raw_reflex_scaling=False,source_configuration=str(base_config),source_configuration_sha256=I.sha(base_config),
                source_profile_sha256=I.sha(B.read(base_config)['integrated_fsa_profile']),motor_proof=str(motor_proof),motor_proof_sha256=I.sha(motor_proof),
                topology_sha256=base_topology,neurons=len(c.V),logging_name_count=len(n.names),edges=len(c.pre),pairing=pairing,
                fixed_RC_sources=rc_sources,fixed_RC_operator=rc_proof,
                matrix_ledger=matrix_ledger,actual_command_current_checks=command,body_steps=0,
                scope='Ten selected shared scalars; not all fitted coefficients, physical parameters, or whole-model robustness' if mode=='selected10' else 'Separate electrical matrix proposal; geometry/unpaired/unfitted/frozen entries explicitly excluded; body review required')
    B.write(out/'CASE.json',record);return record


def prepare(base_config,motor_proof,out,mode='selected10',static_fixture=False):
    out=B.owned(out);out.mkdir(parents=True,exist_ok=False)
    result=Path(base_config).with_suffix('.result.json');r=B.read(result)
    if not static_fixture and not r['assessment']['accepted']:raise ValueError('Final base has not passed strict S4')
    if not static_fixture and (r['metrics'].get('bal_fell') or r['metrics'].get('fall')):raise ValueError('Fallen baseline')
    default=make_case(base_config,motor_proof,out/'default',None,mode=mode)
    for seed in range(110301,110309):make_case(base_config,motor_proof,out/f'seed_{seed}',seed,mode=mode)
    plan=dict(source=str(Path(__file__).resolve()),source_sha256=I.sha(__file__),base_configuration=str(Path(base_config).resolve()),
              base_result=str(result),base_result_sha256=I.sha(result),mode=mode,body_ready=not static_fixture and mode=='selected10',
              seeds=list(range(110301,110309)),max_processes=3,default_processes=1,absolute_deadline_epoch_s=DEADLINE,per_case_timeout_s=240.,
              verification_stop_utc=campaign_budget.VERIFICATION_STOP.isoformat(),budget_source=str(Path(campaign_budget.__file__).resolve()),budget_source_sha256=I.sha(campaign_budget.__file__),
              physical_runner=str(I.H/'physical_eval.py'),physical_runner_sha256=I.sha(I.H/'physical_eval.py'),
              optional_RC_sources=default['fixed_RC_sources'],
              scope=default['scope'],case_paths=[str(out/f'seed_{s}') for s in range(110301,110309)],body_trials_started=0)
    B.write(out/'PLAN.json',plan);return plan


def run_case(folder,deadline,timeout,first):
    deadline=min(DEADLINE,deadline)
    folder=Path(folder);done=folder/'BODY_RESULT.json';active=folder/'ACTIVE.json'
    case_sha=I.sha(folder/'CASE.json')
    if done.exists():
        saved=B.read(done)
        if saved['case_sha256']!=case_sha:raise ValueError('Completed case provenance changed')
        return saved
    if deadline-time.time()<30:return dict(skipped_deadline=True,folder=str(folder))
    try:
        with active.open('x',encoding='utf-8') as f:json.dump(dict(host=socket.gethostname(),owner_pid=os.getpid(),started_epoch_s=time.time()),f)
    except FileExistsError:return dict(skipped_active=True,folder=str(folder),reason='No automatic stale-active takeover; verify process ownership before recovery')
    try:
        cfg=B.read(folder/'CONFIGURATION.json');cfg['diagnose_swing']=bool(first)
        cfg['verification_capture']=str(folder/'body_verification')
        # Avoid modifying source-hashed CONFIGURATION used by its profile.
        body_cfg=folder/'BODY_CONFIGURATION.json';B.write(body_cfg,cfg)
        try:output=B.child([sys.executable,str(I.H/'physical_eval.py'),str(body_cfg)],min(deadline,time.time()+timeout))
        except __import__('subprocess').TimeoutExpired:
            failed=dict(folder=str(folder),case_sha256=case_sha,timeout=True,assessment={'accepted':False},objective={'unsafe':True})
            B.write(done,failed);return failed
        (folder/'BODY_TAIL.log').write_text('\n'.join(line[:2000] for line in output.splitlines()[-30:]))
        r=B.read(folder/'BODY_CONFIGURATION.result.json');scored=G.measure(r,folder/'BODY_CONFIGURATION.assistance.npz')
        d=dict(folder=str(folder),case_sha256=case_sha,result=str(folder/'BODY_CONFIGURATION.result.json'),assessment=r['assessment'],objective=scored)
        if scored['strict_accepted'] and not first:
            capture=(folder/'body_verification').resolve()
            if not capture.is_relative_to(B.ROOT) or capture.parent!=folder.resolve():raise ValueError('Invalid owned capture cleanup target')
            removed=[]
            for file in capture.iterdir():
                if file.is_file():removed.append(dict(path=str(file),sha256=I.sha(file)));file.unlink()
            d['passing_capture_pruned']=removed
        B.write(done,G.finite_json(d));return d
    finally:active.unlink(missing_ok=True)


def run(plan_path,workers=1):
    p=B.read(plan_path)
    if not p['body_ready'] or p['mode']!='selected10':raise ValueError('Static fixture/matrix proposal not authorized for bodies')
    if I.sha(p['source'])!=p['source_sha256'] or I.sha(p['base_result'])!=p['base_result_sha256']:raise ValueError('Basin source changed')
    if I.sha(p['budget_source'])!=p['budget_source_sha256']:raise ValueError('Campaign budget changed; review/rebind')
    if I.sha(p['physical_runner'])!=p['physical_runner_sha256']:raise ValueError('Physical runtime changed; review/rebind')
    for path,digest in p.get('optional_RC_sources',{}).items():
        if I.sha(path)!=digest:raise ValueError('Optional RC source changed; review/rebind')
    if not 1<=workers<=p['max_processes']<=3:raise ValueError('At most three; coordinate shared worker cap externally')
    deadline=min(DEADLINE,p['absolute_deadline_epoch_s'])
    if deadline-time.time()<30:raise TimeoutError('Insufficient verification reserve')
    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(run_case,q,deadline,p['per_case_timeout_s'],i==0) for i,q in enumerate(p['case_paths'])]
        rows=[j.result() for j in as_completed(jobs)]
    B.write(Path(plan_path).parent/'BASIN_RESULT.json',G.finite_json(dict(rows=rows,workers=workers,scope=p['scope'],all_completed=len(rows)==8 and all('assessment' in x for x in rows),
            strict_pass_count=sum(bool(x.get('assessment',{}).get('accepted')) and not x.get('objective',{}).get('unsafe',True) for x in rows))))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--base-config');parser.add_argument('--motor-proof');parser.add_argument('--out')
    parser.add_argument('--mode',choices=('selected10','matrix_probe'),default='selected10');parser.add_argument('--static-fixture',action='store_true')
    parser.add_argument('--run');parser.add_argument('--workers',type=int,default=1);a=parser.parse_args()
    if a.run:run(a.run,a.workers)
    elif a.base_config and a.motor_proof and a.out:print(json.dumps(prepare(a.base_config,a.motor_proof,a.out,a.mode,a.static_fixture)))
    else:parser.error('Provide preparation arguments or an explicit --run plan')
