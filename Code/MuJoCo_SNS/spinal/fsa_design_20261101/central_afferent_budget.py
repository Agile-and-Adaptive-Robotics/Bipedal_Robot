"""Opt-in central afferent CURRENT budgets on existing compiled Gait edges.

No new topology, sensor encoding, NaP tuning, body binding, or Hunt adapter.
The source w2lvar ALREADY divides its central gains by muscle-family counts.
This candidate is a new current-domain contract, not a missing-count fix.
"""
from __future__ import annotations
import argparse,ast,copy,hashlib,inspect,json,os,sys
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np

HERE=Path(__file__).resolve().parent;SPINAL=HERE.parent;REPO=SPINAL.parents[2]
EXTERNAL=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')
OUT=EXTERNAL/'fsa_design/central_afferent_budget'
CASES={'s3k':('other_family_pressure/s3k_bio_prepared_v1','rustuna_s3k_v1/BEST.json'),
       'w2lvar':('normalized_prepared_transfer/w2lvar_bio_graded_v1','group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json'),
       'syn6':('pinned_preparation/syn6_prepared_bio_nominal_v1','rustuna_syn6_v1/BEST.json')}
SCHEMA='aarl_central_afferent_current_budget_v1'


def default_profile():
    return dict(schema=SCHEMA,operating_post_mV=2.5,absolute_current_per_post_type_sign_nA=.5,
                runner_conditions=dict(gamma_coact=0.,vest_prop=0.,muscle_afferent_scale=1.),
                source_current_upper_nA=dict(Ia=2.88,II=2.92,Ib=1.125),
                interpretation='Engineering incremental current at declared post Vm, independently per post/type/sign; not a NaP equilibrium voltage or global current clamp.')


def core(net):return getattr(net.compiled,'c',net.compiled)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def mirror(name):return name[:-2]+('_l' if name.endswith('_r') else '_r') if name.endswith(('_r','_l')) else None
def clone(forward):
    item=object.__new__(type(forward));item.__dict__.update(copy.deepcopy(forward.__dict__));return item


def build_saved(family,projection_callback=None):
    """Validate saved baseline first; optionally rebuild with a NEW projection.

    The callback replaces only the constructor's local tie_weights binding.
    Candidate output never inherits the legacy averaged-graph hash claim.
    """
    if projection_callback is not None and not callable(projection_callback):raise ValueError('Explicit projection callback must be callable')
    os.environ['AARL_NET']=family
    sys.path[:0]=[str(SPINAL),str(SPINAL/'w2l_cpg'),str(SPINAL.parent/'robot_torso/bpa_neural_coordination')]
    import gait_standing_family as F
    root=EXTERNAL/'gait_bpa';stem,config=CASES[family];prefix=root/stem
    saved=json.loads(prefix.with_suffix('.result.json').read_text());cfg=json.loads((root/config).read_text())
    settings=cfg['global_parameters'] if family=='w2lvar' else cfg['trial']['params']
    checked={};changed=[]
    for name in ('runner.py','params.py','build_network.py','build_network_w2lvar.py','build_network_syn6.py'):
        checked[name]=sha(SPINAL/name)
        if checked[name]!=saved['source_sha256'][name]:changed.append(name)
    if any(name not in ('build_network.py','build_network_w2lvar.py','build_network_syn6.py') for name in changed):raise ValueError('Non-builder saved source changed: '+str(changed))
    if saved['failure'] is not None or saved['family']!=family:raise ValueError('Not a saved accepted biological standing case')
    F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.
    m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));acts=[m.actuator(i).name for i in range(m.nu)]
    capture={}
    class Captured(BaseException):pass
    def snapshot(inner,names):
        capture.update(net=inspect.currentframe().f_back.f_locals['n'],names=names)
        raise Captured()
    tree=ast.parse(Path(F.__file__).read_text());outer=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run');node=next(x for x in outer.body if isinstance(x,ast.FunctionDef) and x.name=='build')
    module=ast.Module(body=[node],type_ignores=[]);ast.fix_missing_locations(module);env=dict(F.__dict__)
    env.update(oldbuild=F.B.build,capture=capture,gain=settings['pf_mn_gain'],ib_group_gain=settings['ib_group_gain'],pathway_group_scales=cfg.get('pathway_group_scales'),backend='sparse',compiled_transform=snapshot,fast_motor_io=True,pinned_actuator_settling=True)
    exec(compile(module,str(Path(F.__file__))+'::design_build','exec'),env)
    try:env['build'](acts,F.P.DT,True)
    except Captured:pass
    else:raise AssertionError('Saved build not intercepted')
    if json.loads((prefix.parent/(prefix.name+'_projection.json')).read_text())!=capture['projection']:raise ValueError('Saved bilateral projection changed')
    ledger_path=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'/(family+'_graph_ledger.json')
    ledger=json.loads(ledger_path.read_text(encoding='utf-8'));n=capture['net'];c=core(n)
    matrix_checks={key:hashlib.sha256(np.asarray(getattr(c,key)).tobytes()).hexdigest()==value for key,value in ledger['compiled_matrix_sha256'].items()}
    if not all(matrix_checks.values()) or c.dt!=ledger['dt_s']:raise ValueError('Archived compiled matrix/timestep mismatch after source-hook change')
    node_checks=[]
    for i,(pop,node) in enumerate(zip(n.net.populations,ledger['nodes'])):
        actual=pop.get('params',getattr(pop.get('type'),'params',{}));actual=json.loads(json.dumps(actual,default=lambda v:v.tolist() if hasattr(v,'tolist') else v))
        if pop['name']!=node['name'] or actual!=node['source_neuron_parameters']:raise ValueError('Archived source neuron parameters differ: '+node['name'])
        for key,param in (('c_m','membrane_capacitance'),('g_m','membrane_conductance'),('V_rest','resting_potential'),('i_b','bias')):
            if getattr(c,key)[i]!=actual[param]:raise ValueError('Actual membrane field differs from archived source: '+node['name'])
        if actual['gated']:
            for field,param in (('g_ion','Gion'),('e_ion','Eion')):
                if not np.array_equal(getattr(c,field)[:,i],actual[param]):raise ValueError('Actual ionic field differs: '+field)
            for suffix,component in (('a','paramsA'),('b','paramsB'),('c','paramsC')):
                for field,param in ((f'pow_{suffix}','pow'),(f'k_{suffix}','k'),(f'slope_{suffix}','slope'),(f'e_{suffix}','reversal')):
                    if not np.array_equal(getattr(c,field)[:,i],actual[component][param]):raise ValueError('Actual ionic kinetic field differs: '+field)
                if 'TauMax' in actual[component] and not np.array_equal(getattr(c,'tau_max_'+suffix)[:,i],actual[component]['TauMax']):raise ValueError('Actual ionic time constant differs')
        node_checks.append(True)
    adapters=[dict(input=n.inputs[port],target=capture['names'][post],weight=float(c.input_connectivity[post,port]),units='nA') for post,port in zip(*np.nonzero(c.input_connectivity))]
    if adapters!=ledger['input_adapters'] or len(node_checks)!=len(c.V):raise ValueError('Archived input mapping/census differs')
    if any(F.P.G.get(k,0.)!=0 for k in ('gamma_coact','vest_prop')):raise ValueError('Canonical source-current bounds do not cover current gamma/vest settings')
    bounds=dict(Ia=F.P.AFF['ia_gain']*(1+F.P.MOD['ia']),II=F.P.AFF['i0_ii']+F.P.AFF['ii_gain']*(1+F.P.MOD['ii']),Ib=.5*F.P.AFF['ib_gain']*1.5)
    if any(not np.isclose(bounds[k],v,rtol=0,atol=1e-14) for k,v in default_profile()['source_current_upper_nA'].items()):raise ValueError('Actual AFF/MOD values changed; canonical bounds must be rederived')
    anchor=dict(family=family,saved_result=str(prefix.with_suffix('.result.json')),saved_config=str(root/config),saved_constructor_settings=copy.deepcopy(settings),saved_pathway_group_scales=copy.deepcopy(cfg.get('pathway_group_scales')),source_sha256=checked,default_off_builder_source_changes=changed,compiled_graph_anchor=str(ledger_path),compiled_graph_anchor_sha256=sha(ledger_path),archived_matrix_hashes_exact=matrix_checks,archived_neuron_parameters_and_compiled_cap_ionic_fields_exact=all(node_checks),archived_input_bindings_exact=True,actual_AFF=copy.deepcopy(F.P.AFF),actual_MOD=copy.deepcopy(F.P.MOD),derived_source_upper_nA=bounds,derivation='Runner clips length/velocity to[-1,1],force/Fmax to[0,1.5],drive/stance to[0,1]. Ia=1.6*(1+.8),II=1+1.2*(1+.6),Ib=.5*1.5*1.5. gamma0/vest_prop0/afferent scale1. This is a conservative input envelope, not simultaneous physical attainability.',runner_sha256=sha(F.__file__),neural_only=True,body_binding_calls=0)
    if projection_callback is None:return n,anchor,F.P
    # Reuse the identical extracted saved constructor,not a hand-copied build.
    capture={};env.update(capture=capture,tie_weights=projection_callback)
    try:env['build'](acts,F.P.DT,True)
    except Captured:pass
    else:raise AssertionError('Projected candidate build not intercepted')
    candidate=capture['net']
    return candidate,dict(family=family,baseline_saved_validation=anchor,saved_constructor_settings=copy.deepcopy(settings),saved_pathway_group_scales=copy.deepcopy(cfg.get('pathway_group_scales')),actual_AFF=copy.deepcopy(F.P.AFF),actual_MOD=copy.deepcopy(F.P.MOD),derived_source_upper_nA=bounds,new_projection_policy=capture['projection'],candidate_is_new_condition=True,candidate_claims_legacy_graph_hash=False,candidate_compiled_sha256={key:hashlib.sha256(np.asarray(getattr(core(candidate),key)).tobytes()).hexdigest() for key in ledger['compiled_matrix_sha256']},body_binding_calls=0),F.P


def propose(net,profile):
    p=copy.deepcopy(profile)
    if p.get('schema')!=SCHEMA or p['runner_conditions']!=dict(gamma_coact=0.,vest_prop=0.,muscle_afferent_scale=1.) or p['source_current_upper_nA']!=dict(Ia=2.88,II=2.92,Ib=1.125):raise ValueError('Unsupported input-domain contract; derive changed runner bounds explicitly')
    if p['absolute_current_per_post_type_sign_nA']<=0 or not np.isfinite([p['operating_post_mV'],p['absolute_current_per_post_type_sign_nA']]).all():raise ValueError('Invalid budget')
    c=core(net);names=[pop['name'] for pop in net.net.populations];ids={name:i for i,name in enumerate(names)}
    if len(names)!=len(c.V):raise ValueError('Population/index mismatch')
    groups={};excluded=[];composites=[]
    def homologous_cells(a,b):
        for field in ('c_m','g_m','V_rest','i_b'):
            if getattr(c,field)[a]!=getattr(c,field)[b]:return False
        for field in ('g_ion','e_ion','pow_a','k_a','slope_a','e_a','pow_b','k_b','slope_b','e_b','tau_max_b','pow_c','k_c','slope_c','e_c','tau_max_c'):
            if not np.array_equal(getattr(c,field)[:,a],getattr(c,field)[:,b]):return False
        return True
    for post,pre in zip(*np.nonzero(c.g_max_non)):
        a,b=names[pre],names[post]
        if not b.startswith(('RG_','PF_','InE_','InF_')):continue
        kind=next((k for k in ('Ia','II','Ib') if a.startswith(k+'_')),None)
        if kind is None:
            if a.startswith(('LBIN_','IaIN_','IIX_','IIIN_','IBIN_','IBEXC_')):
                composites.append(dict(source=a,target=b,gmax_uS=float(c.g_max_non[post,pre]),status='Retained composite, no automatic direct-afferent bound',reason='LBIN combines LOAD current and recurrent/group paths; other INs include central feedback. Separate qualified bound required.'))
            continue
        incoming=np.flatnonzero(c.g_max_non[pre]);ports=np.flatnonzero(c.input_connectivity[pre])
        if len(incoming) or c.g_ion[:,pre].any() or len(ports)!=1 or c.V_rest[pre]!=0 or c.i_b[pre]!=0 or c.g_m[pre]<=0:
            excluded.append(dict(source=a,target=b,reason='Not an isolated passive zero-rest direct input encoder'));continue
        gain=float(c.input_connectivity[pre,ports[0]]);upper=p['source_current_upper_nA'][kind]*gain/float(c.g_m[pre]);lo=float(c.e_lo[post,pre]);hi=float(c.e_hi[post,pre]);reverse=float(c.del_e[post,pre]);drive=reverse-p['operating_post_mV']
        if gain<=0 or hi<=lo or c.dt*c.g_m[pre]/c.c_m[pre]>1 or drive==0:excluded.append(dict(source=a,target=b,reason='Unqualified gain/window/Euler pole/operating drive'));continue
        fraction=float(np.clip((upper-lo)/(hi-lo),0.,1.));sign=1 if drive>0 else -1;key=(b,kind,sign)
        groups.setdefault(key,[]).append(dict(source=a,target=b,post=int(post),pre=int(pre),kind=kind,sign=sign,old_gmax_uS=float(c.g_max_non[post,pre]),source_upper_mV=upper,source_upper_nA=p['source_current_upper_nA'][kind],threshold_mV=lo,saturation_mV=hi,maximum_fraction=fraction,reversal_mV=reverse))
    ready=[];unmatched=[];seen=set()
    for key,rows in sorted(groups.items()):
        if key in seen:continue
        target,kind,sign=key;other=(mirror(target),kind,sign);twins=groups.get(other)
        mapped={r['source']:r for r in twins or []}
        comparable=bool(twins) and len(twins)==len(rows)
        for row in rows:
            twin=mapped.get(mirror(row['source']))
            comparable=comparable and twin is not None and all(row[k]==twin[k] for k in ('old_gmax_uS','source_upper_mV','threshold_mV','saturation_mV','reversal_mV')) and homologous_cells(row['pre'],twin['pre']) and homologous_cells(row['post'],twin['post'])
        if not comparable:
            unmatched.append(dict(target=target,kind=kind,sign=sign,reason='Complete identical mirrored graph/fields not verified; group retained unchanged'));seen.add(key);continue
        # One half determines both halves; no independent left/right fitting.
        reference=rows if target.endswith('_r') else twins;mirror_rows=twins if target.endswith('_r') else rows
        original=sum(r['old_gmax_uS']*r['maximum_fraction']*abs(r['reversal_mV']-p['operating_post_mV']) for r in reference)
        if original<=0:unmatched.append(dict(target=target,kind=kind,reason='Source domain cannot activate these existing edges'));seen.update((key,other));continue
        factor=p['absolute_current_per_post_type_sign_nA']/original
        for half in (reference,mirror_rows):
            for row in half:row['new_gmax_uS']=row['old_gmax_uS']*factor
            ready.append(dict(target=half[0]['target'],kind=kind,sign=sign,old_maximum_reference_current_nA=sign*original,new_maximum_reference_current_nA=sign*p['absolute_current_per_post_type_sign_nA'],shared_mirrored_scale=factor,edges=half))
        seen.update((key,other))
    return dict(profile=p,groups=ready,unmatched_groups=unmatched,unqualified_direct_edges=excluded,composite_inventory=composites,
                source_count_normalization_already_present='w2lvar central source divides by _nE/_nF; this is new current-domain budgeting, not a missing-count correction',
                additional_scope='HUNT_AF excluded. Composite paths not silently renamed. Budget is per post/type/sign; totals across types may exceed0.5nA. NaP current contract, not voltage equilibrium.',
                training_allowed=False,Hunt_runtime_integration=False)


def install(net,profile=None):
    if profile is None:return dict(applied=False,default_noop=True)
    if hasattr(net,'apply_weights') or hasattr(net,'hunt_timing_normalized_design'):raise ValueError('Hunt dynamic masks require a separate verified policy adapter; no integration claim')
    if hasattr(net,'central_afferent_budget_design'):raise ValueError('Already installed')
    if isinstance(profile,(str,Path)):profile=json.loads(Path(profile).read_text())
    design=propose(net,profile);c=core(net);forward=net.compiled
    cache={(int(post),int(pre)):i for i,(post,pre) in enumerate(zip(forward.post,forward.pre))} if hasattr(forward,'gmax') else None
    for group in design['groups']:
        for row in group['edges']:
            if cache is not None and (row['post'],row['pre']) not in cache:raise ValueError('Missing actual cached source edge')
    for group in design['groups']:
        for row in group['edges']:
            c.g_max_non[row['post'],row['pre']]=row['new_gmax_uS']
            if cache is not None:forward.gmax[cache[row['post'],row['pre']]]=row['new_gmax_uS']
    net.central_afferent_budget_design=design;verify_runtime(net);return design


def verify_runtime(net):
    c=core(net);forward=net.compiled;cache={(int(p),int(q)):i for i,(p,q) in enumerate(zip(forward.post,forward.pre))} if hasattr(forward,'gmax') else None
    for group in net.central_afferent_budget_design['groups']:
        for row in group['edges']:
            if c.g_max_non[row['post'],row['pre']]!=row['new_gmax_uS'] or (cache is not None and forward.gmax[cache[row['post'],row['pre']]]!=row['new_gmax_uS']):raise AssertionError('Compiled central budget overwritten; caller policy must reapply explicitly')
    return True


def run(family,out):
    from graded_profile import field_digest
    net,anchor,parameters=build_saved(family);c=core(net);before=field_digest(c);install(net,None);default_exact=field_digest(c)==before
    all_arrays={k:v.copy() for k,v in c.__dict__.items() if isinstance(v,np.ndarray) and k!='g_max_non'};old=c.g_max_non.copy()
    design=install(net,default_profile());selected=[r for g in design['groups'] for r in g['edges']];mask=np.ones_like(old,dtype=bool)
    for row in selected:mask[row['post'],row['pre']]=False
    other_weights_exact=np.array_equal(old[mask],c.g_max_non[mask]);all_other_fields_exact=all(np.array_equal(value,getattr(c,k)) for k,value in all_arrays.items())
    names=[p['name'] for p in net.net.populations];targets=sorted({r['post'] for r in selected});cache={(int(p),int(q)):i for i,(p,q) in enumerate(zip(net.compiled.post,net.compiled.pre))};probes=[]
    for kind in ('Ia','II','Ib','ALL'):
        chosen=[r for r in selected if kind=='ALL' or r['kind']==kind]
        if not chosen:continue
        for fraction in (0.,.25,.5,1.):
            on=clone(net.compiled);off=clone(net.compiled)
            for row in selected:
                on.c.V[row['pre']]=row['source_upper_mV']*fraction;off.c.V[row['pre']]=on.c.V[row['pre']]
            on.c.V[targets]=2.5;off.c.V[targets]=2.5;expected={post:0. for post in targets}
            for row in chosen:
                off.gmax[cache[row['post'],row['pre']]]=0.
                activation=np.clip((row['source_upper_mV']*fraction-row['threshold_mV'])/(row['saturation_mV']-row['threshold_mV']),0,1)
                expected[row['post']]+=row['new_gmax_uS']*activation*(row['reversal_mV']-2.5)
            u=net.make_inputs();on.forward(u);off.forward(u);measured=(on.V[targets]-off.V[targets])/on.c.time_factor_membrane[targets]
            error=max(abs(float(v)-expected[post]) for post,v in zip(targets,measured))
            probes.append(dict(kind=kind,source_bound_fraction=fraction,maximum_native_incremental_current_error_nA=error))
    # Actual net.step and pure-input transient, without a mechanical body.
    source_ids={row['pre']:row for row in selected};bound_violation=0.;finite=True
    for step in range(120):
        u=net.make_inputs();amplitude=1. if 40<=step<80 else 0.
        for row in source_ids.values():u[net.input_index(row['source'])]=amplitude*row['source_upper_nA']
        net.step(u);verify_runtime(net);finite=finite and bool(np.isfinite(c.V).all())
        for i,row in source_ids.items():bound_violation=max(bound_violation,-float(c.V[i]),float(c.V[i])-row['source_upper_mV'])
    passed=default_exact and other_weights_exact and all_other_fields_exact and finite and bound_violation<1e-10 and all(p['maximum_native_incremental_current_error_nA']<1e-10 for p in probes)
    result=dict(schema=SCHEMA,family=family,anchor=anchor,neurons=len(c.V),input_count=len(net.inputs),selected_existing_edges=len(selected),selected_post_type_sign_groups=len(design['groups']),default_None_exact=default_exact,all_other_numeric_fields_exact_at_install=all_other_fields_exact,nonselected_weights_exact=other_weights_exact,native_current_probes=probes,source_bound_transient_violation_mV=bound_violation,actual_base_net_step_finite=finite,actual_step_policy_checked=120,
                electrical_passed=passed,applicability='existing selected groups' if selected else 'No qualifying existing active direct-central paths: inventory only; no invented edges',global_family_design_passed=False,Hunt_integration=False,physics_steps=0,training_steps=0,helper_sha256=sha(__file__))
    out.mkdir(parents=True,exist_ok=True)
    for name,payload in (('PROFILE.json',design['profile']),('APPLIED_LEDGER.json',design),('RESULTS.json',result)):
        path=out/name
        if path.exists():raise FileExistsError('Preserve earlier evidence: '+str(path))
        path.write_text(json.dumps(payload,indent=2))
    print(json.dumps({k:result[k] for k in ('family','selected_existing_edges','selected_post_type_sign_groups','electrical_passed','applicability')}))
    if not passed:raise AssertionError('Central fan-in electrical check failed; no promotion')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--family',choices=tuple(CASES),required=True);parser.add_argument('--out',type=Path);args=parser.parse_args()
    sys.path[:0]=[str(HERE),str(SPINAL),str(SPINAL/'w2l_cpg')]
    run(args.family,args.out or OUT/args.family)
