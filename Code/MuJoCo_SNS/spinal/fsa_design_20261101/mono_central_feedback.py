"""NEW explicit sagittal-monoarticular central afferent completion.

Opt-in compiled overlay, never a source-reproduction or all-muscle/3D claim.
No body integration, optimization, Hunt, or curriculum operation is provided.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,os,sys
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import central_afferent_budget as CB
from graded_profile import assert_campaign_allowed,field_digest,population_names

SCHEMA='aarl_explicit_sagittal_mono_central_feedback_v1'
ROLES={'iliacus':'F','psoas':'F','bifemsh':'F','tib_ant':'F',
       'glut_max1':'E','glut_max2':'E','glut_max3':'E',
       'vas_med':'E','vas_int':'E','vas_lat':'E','soleus':'E'}
TYPES={'F':('Ia','II'),'E':('Ib','II')}

def default_profile(family):
    if family not in ('s3k','syn6'):raise ValueError('This new completion is scoped to s3k/syn6')
    return dict(schema=SCHEMA,family=family,role_by_muscle_base=ROLES.copy(),
        source_current_max_nA=dict(Ia=2.88,II=2.92,Ib=1.125),
        caller_conditions=dict(gamma_coact=0.,vest_prop=0.,muscle_afferent_scale=1.),
        reference_post_mV=2.5,per_type_post_current_nA=.5,reversal_mV=8.,activation_window_mV=[0.,5.],
        relative_RG_IN_weights='equal_per_eligible_muscle',relative_PF_weights='positive_actual_PF_to_MN_conductance_normalized_within_muscle_same_phase',
        existing_selected_policy='reject_nonzero_not_add_or_sum',bilateral_policy='derive_right_then_exact_named_mirror',
        mixed_synergy_policy='exclude_and_report',new_contralateral_paths=False)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2),encoding='utf-8')
def mirror(n):return CB.mirror(n)
def context(P):
    upper=dict(Ia=float(P.AFF['ia_gain']*(1+P.MOD['ia'])),II=float(P.AFF['i0_ii']+P.AFF['ii_gain']*(1+P.MOD['ii'])),Ib=float(.5*P.AFF['ib_gain']*1.5))
    if P.G.get('gamma_coact',0.)!=0 or P.G.get('vest_prop',0.)!=0:raise ValueError('Caller context outside the declared source bound')
    if not np.allclose([upper[k] for k in ('Ia','II','Ib')],[2.88,2.92,1.125],rtol=0,atol=1e-14):raise ValueError('Actual AFF/MOD changed; rederive source bounds')
    return dict(upper_nA=upper,AFF=copy.deepcopy(P.AFF),MOD=copy.deepcopy(P.MOD),gamma_coact=0.,vest_prop=0.,caller_muscle_afferent_scale=1.)

def phase(net,name):
    if not name.startswith('PF_') or name.startswith('PF_IN_'):return None
    side=name[-1];base=name[:-2]
    if base.startswith('PF_S'):
        k=int(base[4:])-1
        if bool(net.syn_mixed[side][k]):return 'mixed'
        return str(net.syn_families[side][k])
    if base.startswith('PF_E') or base.endswith('-E'):return 'E'
    if base.startswith('PF_F') or base.endswith('-F'):return 'F'
    return None

def validate_input_currents(net,x,*,P):
    ledger=net.mono_central_feedback_ledger
    if context(P)!=ledger['source_input_context']:raise ValueError('Caller source settings changed')
    x=np.asarray(x,float)
    if x.shape!=(len(net.inputs),) or not np.isfinite(x).all():raise ValueError('Invalid neural input vector')
    checked=set()
    for group in ledger['groups']:
        for r in group['edges']:
            if r['source'] in checked:continue
            checked.add(r['source']);value=x[net.inputs.index(r['input_port'])]*r['input_weight']
            if not 0<=value<=r['source_max_current_nA']*(1+1e-12):raise ValueError('Physical current outside new central-source domain: '+r['source'])
    return x

def propose(net,p,P):
    assert_campaign_allowed('optional_profile_implementation')
    if p!=default_profile(p.get('family')):raise ValueError('Unknown/missing/unused profile fields; register a different design explicitly')
    if getattr(net,'graded_fsa_family',None)!=p['family']:raise ValueError('Actual family does not match profile')
    if hasattr(net,'apply_weights') or getattr(net.compiled,'inner',None) is not None:raise ValueError('Dynamic/opaque neural wrappers need their own integration; apply before wrapping')
    c=CB.core(net);names=population_names(net);ix={n:i for i,n in enumerate(names)};ctx=context(P)
    groups={};routing=[];unresolved=[];support=[]
    # Both sides are enumerated, then actual fields and membership must mirror.
    for act,info in sorted(net.muscles.items()):
        role=ROLES.get(info.base)
        if role is None:
            reserved=info.base in ('quad_fem','gem','peri')
            unresolved.append(dict(muscle=act,groups=list(info.groups),status='source_force_pruned_reserved' if reserved else 'unreviewed_active_route',reason='Source force-pruned route retained, not silently assigned' if reserved else 'Outside explicitly reviewed sagittal-mono set; biarticular/secondary-axis/distal/trunk or phase-role decision not supplied'))
            continue
        mn=ix['MN_'+act];pf=[]
        for i,name in enumerate(names):
            if name[-2:]!='_'+info.side or not name.startswith('PF_'):continue
            g=float(c.g_max_non[mn,i])
            if g<=0 or c.del_e[mn,i]<=c.V_rest[mn]:continue
            family=phase(net,name);support.append(dict(muscle=act,PF=name,phase_family=family,actual_PF_MN_g_uS=g,eligible=family==role))
            if family==role:pf.append((name,g))
        total=sum(g for _,g in pf)
        if not pf:unresolved.append(dict(muscle=act,reason='No positive same-phase PF→MN support; no all-PF fallback invented'))
        destinations=[('RG_'+role+'_'+info.side,1.),('In'+role+'_'+info.side,1.)]+[(n,g/total) for n,g in pf]
        for typ in TYPES[role]:
            source=typ+'_'+act
            if source not in ix:raise ValueError('Missing actual sensory node '+source)
            pre=ix[source];ports=np.flatnonzero(c.input_connectivity[pre]);incoming=np.flatnonzero(c.g_max_non[pre])
            if len(ports)!=1 or len(incoming) or c.i_b[pre]!=0 or c.V_rest[pre]!=0 or c.g_ion[:,pre].any():raise ValueError('Unqualified physical input encoder '+source)
            gain=float(c.input_connectivity[pre,ports[0]])
            if gain<=0 or c.g_m[pre]<=0:raise ValueError('Unqualified source electrical gain')
            vmax=ctx['upper_nA'][typ]*gain/float(c.g_m[pre]);fraction=float(np.clip(vmax/5.,0,1))
            for target,relative in destinations:
                if target not in ix:raise ValueError('No actual central target '+target)
                post=ix[target]
                if c.g_max_non[post,pre]!=0:raise ValueError('Existing selected conditional connection: explicit rejection avoids double addition '+source+'→'+target)
                row=dict(source=source,target=target,pre=pre,post=post,type=typ,phase=role,relative_weight=relative,source_max_current_nA=ctx['upper_nA'][typ],source_max_voltage_mV=vmax,source_max_activation_fraction=fraction,input_port=net.inputs[int(ports[0])],input_weight=gain)
                groups.setdefault((target,typ),[]).append(row);routing.append(row)
    allpairs=[(r['post'],r['pre']) for r in routing]
    if len(allpairs)!=len(set(allpairs)):raise ValueError('Duplicate planned edge')
    complete=[]
    for (target,typ),rows in sorted(groups.items()):
        if not target.endswith('_r'):continue
        other=groups.get((mirror(target),typ));twins={r['source']:r for r in other or []}
        if len(twins)!=len(rows):raise ValueError('Incomplete mirrored target fanin '+target)
        denominator=sum(r['relative_weight']*r['source_max_activation_fraction']*(p['reversal_mV']-p['reference_post_mV']) for r in rows)
        if denominator<=0:raise ValueError('Inactive/invalid reference current function')
        for r in rows:
            twin=twins.get(mirror(r['source']))
            if twin is None:raise ValueError('Missing mirrored source')
            for key in ('relative_weight','source_max_voltage_mV','input_weight','source_max_activation_fraction'):
                if not np.isclose(r[key],twin[key],rtol=0,atol=1e-14):raise ValueError('Nonhomologous source/PF support '+key)
            for field in ('c_m','g_m','V_rest','i_b'):
                if getattr(c,field)[r['pre']]!=getattr(c,field)[twin['pre']] or getattr(c,field)[r['post']]!=getattr(c,field)[twin['post']]:raise ValueError('Nonhomologous membrane parameters; canonical policy required')
            for field in ('g_ion','e_ion','pow_a','k_a','slope_a','e_a','pow_b','k_b','slope_b','e_b','tau_max_b','pow_c','k_c','slope_c','e_c','tau_max_c'):
                if not np.array_equal(getattr(c,field)[:,r['post']],getattr(c,field)[:,twin['post']]):raise ValueError('Nonhomologous postsynaptic ionic parameters')
            r['gmax_uS']=p['per_type_post_current_nA']*r['relative_weight']/denominator;twin['gmax_uS']=r['gmax_uS']
        for half in (rows,other):complete.append(dict(target=half[0]['target'],type=typ,reference_current_nA=p['per_type_post_current_nA'],full_fanin=len(half),edges=half))
    if sum(len(g['edges']) for g in complete)!=len(routing):raise ValueError('Unmatched source/target groups')
    rules=HERE.parent/'replication/shinohara_rules.json';rulings=HERE.parent/'WIRING_RULINGS_20261003.md';skill=Path(r'C:\Users\Ben Bolen\.agents\skills\neuroscience\SKILL.md')
    inherited=dict(files={str(f):sha(f) for f in (rules,rulings,skill)},F_Ia_II='Local replication shinohara_rules.json shin2025_eq10 target set RG-F/IN-F/PF-F supports positive flexor velocity/length feedback. Actual bounded SNSIa/II encoding differs from its v^0.6 and90%-length equations; no copied paper gain or exact replication.',E_Ib='Local replication shin2025_eq11 target set RG-E/IN-E/PF-E supports positive extensor force feedback. Human mono roles are an explicit engineering adaptation.',E_II='Local neuroscience same-group typeII rule plus Ben-requested EII central extension; NOT attributed to Shinoharaeq11 force formula.',implementation='New engineering gains from independent current budgets, actual PF support and declared roles. No extensorIa universal sign, reflex-state reversal, biarticular integration or historical-dissertation repair claimed.')
    return dict(profile=p,source_input_context=ctx,groups=complete,PF_support=support,unresolved_muscles=unresolved,new_edges=len(routing),new_neurons=0,source_numeric_sha256=field_digest(c),scope='NEW sagittal-mono subset only. All group currents are independent per type/postsynaptic cell; budgets across types add. Noncentral/local reflexes unchanged. No contralateral route, extensorIa universal sign, biarticular integration, or independent3D PF claim.',rule_inheritance=inherited,training_allowed=False)

def install(net,p=None,*,P=None):
    if p is None:return dict(applied=False,default_noop=True)
    if P is None:raise ValueError('Actual caller parameter module required')
    plan=propose(net,p,P);c=CB.core(net);old={k:getattr(c,k).copy() for k in ('g_max_non','del_e','e_lo','e_hi')};states={k:getattr(c,k).copy() for k in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last')}
    selected=[]
    for group in plan['groups']:
        for r in group['edges']:
            t,s=r['post'],r['pre'];selected.append((t,s));c.g_max_non[t,s]=r['gmax_uS'];c.del_e[t,s]=p['reversal_mV'];c.e_lo[t,s],c.e_hi[t,s]=p['activation_window_mV']
    mask=np.ones(c.g_max_non.shape,dtype=bool)
    for t,s in selected:mask[t,s]=False
    for key,a in old.items():assert np.array_equal(getattr(c,key)[mask],a[mask])
    for key,a in states.items():assert np.array_equal(getattr(c,key),a)
    from sparse_pairwise_eval import SparsePairwiseForward
    if hasattr(net.compiled,'c'):net.compiled=SparsePairwiseForward(c)
    net.mono_central_feedback_profile=copy.deepcopy(p);net.mono_central_feedback_ledger=plan
    plan.update(compiled_after_sha256=field_digest(c),preserved_numerical_states=True,nonselected_matrix_entries_exact=True,compiled_graph_policy='Actual dense/sparse runtime is authority; original builder drawing is a baseline and must not be relabelled as the augmented graph')
    verify(net,P=P);return plan

def verify(net,*,P):
    p=net.mono_central_feedback_profile;ledger=net.mono_central_feedback_ledger;c=CB.core(net)
    if context(P)!=ledger['source_input_context']:raise ValueError('Caller context changed')
    if field_digest(c)!=ledger['compiled_after_sha256']:raise ValueError('Profile invalidated by later parameter/topology edit')
    for group in ledger['groups']:
        for r in group['edges']:
            t,s=r['post'],r['pre']
            if c.g_max_non[t,s]!=r['gmax_uS'] or c.del_e[t,s]!=p['reversal_mV'] or [c.e_lo[t,s],c.e_hi[t,s]]!=p['activation_window_mV']:raise ValueError('Unconsumed connection')
    if hasattr(net.compiled,'c'):
        w=net.compiled
        if not np.array_equal(w.gmax,c.g_max_non[w.post,w.pre]):raise ValueError('Stale sparse cache')
        expected={(r['post'],r['pre']) for g in ledger['groups'] for r in g['edges']}
        if not expected.issubset(set(zip(w.post.tolist(),w.pre.tolist()))):raise ValueError('New graph edge absent from accelerator')
    return dict(actual_edges_consumed=ledger['new_edges'],all_new_cached_edges_present=True,parameter_context_verified=True)

def compiled_graph(net,out):
    """Export final runtime truth, not source builder labels after an overlay."""
    c=CB.core(net);w=net.compiled;names=population_names(net);plan=net.mono_central_feedback_ledger
    new={(r['post'],r['pre']) for g in plan['groups'] for r in g['edges']};nodes=[];edges=[]
    for i,name in enumerate(names):
        ion=[]
        for channel in range(c.g_ion.shape[0]):
            if c.g_ion[channel,i]==0:continue
            ion.append({k:float(getattr(c,k)[channel,i]) for k in ('g_ion','e_ion','pow_a','k_a','slope_a','e_a','pow_b','k_b','slope_b','e_b','tau_max_b','pow_c','k_c','slope_c','e_c','tau_max_c')})
        nodes.append(dict(index=i,name=name,electrical_model='graded_ionic' if ion else 'graded_passive',C_microF=float(c.c_m[i]),G_uS=float(c.g_m[i]),rest_mV=float(c.V_rest[i]),bias_nA=float(c.i_b[i]),ionic_channels=ion,PF_phase_role=phase(net,name)))
    for post,pre in zip(*np.nonzero(c.g_max_non)):
        edges.append(dict(source_index=int(pre),target_index=int(post),source=names[pre],target=names[post],gmax_uS=float(c.g_max_non[post,pre]),reversal_mV=float(c.del_e[post,pre]),threshold_mV=float(c.e_lo[post,pre]),saturation_mV=float(c.e_hi[post,pre]),origin='NEW_MONO_CENTRAL_FEEDBACK' if (post,pre) in new else 'RETAINED_ACTUAL_COMPILED_EDGE'))
    rows=np.array([[e['target_index'],e['source_index'],e['gmax_uS'],e['reversal_mV'],e['threshold_mV'],e['saturation_mV']] for e in edges],dtype=np.float64)
    cached=np.c_[w.post,w.pre,w.gmax,w.reversal,w.elo,w.elo+w.erange].astype(np.float64)
    if not np.array_equal(rows,cached):raise ValueError('Exported graph/operator differs from actual sparse cache')
    digest=hashlib.sha256(rows.tobytes()).hexdigest()
    inputs=[dict(input=net.inputs[j],target=names[i],target_index=int(i),weight=float(c.input_connectivity[i,j]),units='nA') for i,j in zip(*np.nonzero(c.input_connectivity))]
    motors=[dict(actuator=act,MN=net.mn_names[act],node_index=int(net.idx[net.mn_names[act]]),conversion='clip(V_MN/5mV,0,1)') for act in net.muscles]
    report=dict(schema='aarl_final_compiled_neural_graph_v1',machine='easteregg2',condition='NEW full canonical-template sagittal-mono central-feedback design',body_run=False,baseline_result_inherited=False,training_allowed=False,draw_from_this_graph_not_source_builder=True,compiled_parameter_sha256=field_digest(c),edge_operator_sha256=digest,cache_operator_sha256=hashlib.sha256(cached.tobytes()).hexdigest(),cache_exact=True,census=dict(nodes=len(nodes),inputs=len(net.inputs),active_chemical_edges=len(edges),new_central_edges=len(new),new_neurons=0),units=dict(time='s',voltage='mV',conductance='uS',current='nA',numeric_capacitance='microF'),nodes=nodes,chemical_edges=edges,input_adapters=inputs,motor_readouts=motors,coverage=dict(reviewed_role_muscles=ROLES,PF_support=plan['PF_support'],omitted_or_unresolved=plan['unresolved_muscles'],scope='22reviewed sagittal-mono muscles;64active unreviewed routes and6force-pruned reserved routes omitted from NEW feedback. Existing base neural arcs retained. No biarticular or independent3D completion.'),source_drawing_policy='Original n.net/source JSON remains a baseline. These numeric fields come from final compiled matrices; no stale source synapse type name is attached.')
    write(out/'FINAL_COMPILED_GRAPH.json',report);return dict(path=str(out/'FINAL_COMPILED_GRAPH.json'),sha256=sha(out/'FINAL_COMPILED_GRAPH.json'),edge_operator_sha256=digest,cache_operator_sha256=digest,counts=report['census'])

def probe(family,out,source_condition='full_template'):
    from sparse_pairwise_eval import SparsePairwiseForward
    from mirror_canonical import apply_canonical_right,default_policy
    out.mkdir(parents=True,exist_ok=True)
    callback=(lambda n:apply_canonical_right(n,default_policy())) if family=='syn6' else None
    os.environ['AARL_NET']=family
    sys.path.insert(0,str(HERE.parent.parent/'robot_torso/bpa_neural_coordination'))
    import gait_standing_family as F
    defaults={k:F.P.G[k] for k in ('ii_set_pose','posture_to_mn')}
    saved,anchor,P=CB.build_saved(family,projection_callback=callback)
    if source_condition=='full_template':
        for k,v in defaults.items():P.G[k]=v
        model=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));net=F.B.build([model.actuator(i).name for i in range(model.nu)],dt=P.DT,interleg=True)
        from gait_design_composition import ANCHORS
        baseline=copy.deepcopy(net.compiled);rng=np.random.default_rng(1101);h=hashlib.sha256()
        for _ in range(128):
            x=rng.uniform(0,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;baseline.forward(x);h.update(baseline.V.tobytes())
        if h.hexdigest()!=ANCHORS[family]:raise ValueError('Unmasked default template anchor differs')
        projection=apply_canonical_right(net,default_policy()) if family=='syn6' else F.tie_weights(net)
        anchor=dict(condition='NEW_full_unmasked_canonical_template',saved_standing_comparison_anchor=anchor,standing_fitted_gains_not_applied=True,source_default_anchor_exact=True,projection=projection)
    else:net=saved
    source=CB.core(net);before=field_digest(source);states={k:getattr(source,k).copy() for k in ('V','V_last','b_gate','c_gate')}
    noop=install(net);assert noop['default_noop'] and field_digest(source)==before
    baseline=SparsePairwiseForward(copy.deepcopy(source));unchanged=SparsePairwiseForward(copy.deepcopy(source));u=np.zeros(len(net.inputs));u[net.input_index('DRIVE')]=2.5
    for _ in range(64):baseline.forward(u);unchanged.forward(u);assert np.array_equal(baseline.V,unchanged.V)
    profile=default_profile(family);plan=propose(net,profile,P);write(out/'RULE_MATRIX.json',plan);write(out/'PROFILE.json',profile)
    net.compiled=SparsePairwiseForward(source);plan=install(net,profile,P=P);c=CB.core(net);rows=[];maximum=0.
    for group in plan['groups']:
        for fraction in (0.,.25,.5,1.):
            v=c.V.copy()
            for e in group['edges']:v[e['pre']]=e['source_max_voltage_mV']*fraction
            # Only one extracted afferent pool is active in this reference current measurement.
            current=sum(c.g_max_non[e['post'],e['pre']]*np.clip(v[e['pre']]/5,0,1)*(profile['reversal_mV']-profile['reference_post_mV']) for e in group['edges'])
            error=abs(current-.5*fraction);maximum=max(maximum,float(error));assert error<1e-12
            # Actual sparse operator comparison via matched per-connection zeroing.
            on=SparsePairwiseForward(copy.deepcopy(c));offc=copy.deepcopy(c)
            for e in group['edges']:offc.g_max_non[e['post'],e['pre']]=0
            off=SparsePairwiseForward(offc);post=group['edges'][0]['post'];on.c.V[:]=off.c.V[:]=v;on.c.V[post]=off.c.V[post]=profile['reference_post_mV'];on.forward(np.zeros(len(net.inputs)));off.forward(np.zeros(len(net.inputs)))
            actual=(on.V[post]-off.V[post])/c.time_factor_membrane[post];actual_error=abs(actual-current);maximum=max(maximum,float(actual_error));assert actual_error<1e-10
            rows.append(dict(target=group['target'],type=group['type'],input_fraction=fraction,expected_nA=.5*fraction,actual_increment_current_nA=float(actual),error_nA=float(actual_error)))
    rejected=False
    try:install(net,profile,P=P)
    except ValueError:rejected=True
    assert rejected
    # Policies are consumed strictly; corrupted caches and changed domains fail.
    cache=net.compiled;first=plan['groups'][0]['edges'][0];slot=np.flatnonzero((cache.post==first['post'])&(cache.pre==first['pre']))[0];old=cache.gmax[slot];cache.gmax[slot]*=1.01
    try:verify(net,P=P);cache_rejected=False
    except ValueError:cache_rejected=True
    finally:cache.gmax[slot]=old
    assert cache_rejected
    invalid=copy.deepcopy(profile);invalid['unused_typo']=1
    try:propose(net,invalid,P);unused_rejected=False
    except ValueError:unused_rejected=True
    assert unused_rejected
    zero=net.make_inputs();assert np.array_equal(validate_input_currents(net,zero,P=P),zero)
    zero[net.inputs.index(first['input_port'])]=first['source_max_current_nA']/first['input_weight']+.01
    before_invalid=c.V.copy()
    try:validate_input_currents(net,zero,P=P);domain_rejected=False
    except ValueError:domain_rejected=True
    assert domain_rejected and np.array_equal(before_invalid,c.V)
    # Node/input/ionic/time fields remain the baseline; only selected compiled chemical edges changed.
    assert len(c.V)==len(source.V)
    for key,a in states.items():assert np.array_equal(getattr(c,key),a)
    verify(net,P=P)
    graph=compiled_graph(net,out)
    files=[__file__,CB.__file__,HERE.parent/'synergy_basis.npz',HERE.parent/'fsa_results/fsa_backsolve.npz']
    result=dict(machine='easteregg2',family=family,source_condition=source_condition,source_anchor=anchor,default_None_exact=True,default64_V_exact=True,new_edges=plan['new_edges'],new_groups=len(plan['groups']),group_endpoint_max_error_nA=maximum,isolated_actual_sparse_current_checks=rows,duplicate_install_rejected=True,corrupted_cache_rejected=cache_rejected,unused_policy_rejected=unused_rejected,out_of_domain_current_rejected_without_mutation=domain_rejected,nonselected_entries_exact=True,neural_states_preserved=True,runtime=verify(net,P=P),final_compiled_graph=graph,PF_support_and_unresolved_in_rule_matrix=True,source_hashes={str(p):sha(p) for p in files},source_models_modified=False,full_all_muscle_or3D_verification=False,training_allowed=False,body_trials=0,optimization_trials=0)
    write(out/'RESULTS.json',result);print(json.dumps(dict(family=family,edges=plan['new_edges'],groups=len(plan['groups']),maximum_error_nA=maximum,unresolved=len(plan['unresolved_muscles']),record=str(out/'RESULTS.json'))),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--family',choices=['s3k','syn6'],required=True);a.add_argument('--source-condition',choices=['full_template','saved_standing'],default='full_template');a.add_argument('--output',type=Path,required=True);args=a.parse_args();probe(args.family,args.output,args.source_condition)
