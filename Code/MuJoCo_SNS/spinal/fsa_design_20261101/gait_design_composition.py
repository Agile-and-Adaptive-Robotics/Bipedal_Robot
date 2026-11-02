"""Gait neural pathway composition reviewer: electrical design only.

No body binding, MuJoCo stepping, tendon/compliance/damping application,
optimization, Hunt fitting, or curriculum. Saved standing masks are preserved;
separate canonical-expression fixtures exercise their dormant pathways.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,copy,hashlib,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import graded_profile as GP,central_afferent_budget as CB,ibexc_coincidence as IB,mirror_canonical as MC

ANCHORS={'s3k':'765e35efff4a8d0eb3c4dc5fd94e1c0d15d36a3c5ee21a6afa8976e8ea841d45','w2lvar':'960fad929ae58737d452797d23501067608377ac97d8947eb979acddac5abc75','syn6':'0283aae4d70531aad60a8816611af7c00dc76d7419a646cc8947ccddb12456a4'}


def states(c):
    return {k:v.copy() for k,v in c.__dict__.items() if isinstance(v,np.ndarray) and k in MC.DYNAMIC_FIELDS}


def verify(n,gate):
    result={'graded':GP.verify_profile(n),'central':CB.verify_runtime(n)}
    if gate:result['coincidence']=IB.verify_coincidence(n)
    if hasattr(n,'canonical_right_ledger'):result['canonical_right']=MC.verify_canonical(n)
    return result


def current_probe(n,slots):
    """Isolate delivered synaptic current in actual compiled one-step subtraction."""
    c=n.compiled;core=CB.core(n);targets=sorted(set(c.post[slots].tolist()));tests=[]
    if not len(slots):return {'active_slots':0,'status':'Masked/inventory only; no path reactivated'}
    for fraction in (0.,.5,1.,1.5):
        on=CB.clone(c);off=CB.clone(c);values={};expected={j:0. for j in targets}
        for i in slots:
            pre,post=int(c.pre[i]),int(c.post[i]);voltage=float(core.e_lo[post,pre]+fraction*(core.e_hi[post,pre]-core.e_lo[post,pre]))
            if pre in values and abs(values[pre]-voltage)>1e-12:raise ValueError('A common source has incompatible probe voltage windows')
            values[pre]=voltage;expected[post]+=float(c.gmax[i]*np.clip(fraction,0,1)*(c.reversal[i]-2.5));off.gmax[i]=0.;off.c.g_max_non[post,pre]=0.
        for pre,v in values.items():on.c.V[pre]=v;off.c.V[pre]=v
        on.c.V[targets]=2.5;off.c.V[targets]=2.5
        u=n.make_inputs();on.forward(u);off.forward(u)
        actual=(on.V[targets]-off.V[targets])/on.c.time_factor_membrane[targets];error=max(abs(float(x)-expected[j]) for j,x in zip(targets,actual))
        if error>1e-10:raise AssertionError('Compiled pathway current differs from declared conductance law')
        tests.append({'source_fraction':fraction,'maximum_native_current_error_nA':error})
    return {'active_slots':len(slots),'status':'Actual compiled current verified; source voltages clamped,not an intrinsic PF waveform fit','tests':tests}


def central_probe(n):
    design=n.central_afferent_budget_design;c=n.compiled;selected=[r for g in design['groups'] for r in g['edges']]
    if not selected:return {'groups':0,'status':'No qualifying existing direct-central group; inventory only'}
    cache={(int(p),int(q)):i for i,(p,q) in enumerate(zip(c.post,c.pre))};targets=sorted({r['post'] for r in selected});tests=[]
    for kind in ('Ia','II','Ib','ALL'):
        chosen=[r for r in selected if kind=='ALL' or r['kind']==kind]
        if not chosen:continue
        for fraction in (0.,.5,1.):
            on=CB.clone(c);off=CB.clone(c);expected={j:0. for j in targets}
            for row in selected:on.c.V[row['pre']]=row['source_upper_mV']*fraction;off.c.V[row['pre']]=on.c.V[row['pre']]
            on.c.V[targets]=2.5;off.c.V[targets]=2.5
            for row in chosen:
                i=cache[(row['post'],row['pre'])];off.gmax[i]=0.;off.c.g_max_non[row['post'],row['pre']]=0.
                activation=np.clip((fraction*row['source_upper_mV']-row['threshold_mV'])/(row['saturation_mV']-row['threshold_mV']),0,1)
                expected[row['post']]+=row['new_gmax_uS']*activation*(row['reversal_mV']-2.5)
            u=n.make_inputs();on.forward(u);off.forward(u);actual=(on.V[targets]-off.V[targets])/on.c.time_factor_membrane[targets]
            error=max(abs(float(x)-expected[j]) for j,x in zip(targets,actual))
            if error>1e-10:raise AssertionError('Central current budget overwritten or misconsumed')
            tests.append({'kind':kind,'bound_fraction':fraction,'maximum_native_current_error_nA':error})
    return {'groups':len(design['groups']),'active_edges':len(selected),'tests':tests,'scope':'Per post/type/sign current contract at2.5mV; not globally additive voltage or NaP equilibrium'}


def domain_and_backend_probe(n,gate):
    bounds=n.graded_fsa_profile['physical_input_domains_nA'];core=CB.core(n);names=GP.population_names(n);largest=0.;finite=True
    upper={n.inputs.index(name):interval[1] for name,interval in bounds.items()};idx={name:names.index(name) for name in bounds}
    for step in range(120):
        u=n.make_inputs();amplitude=0. if step<30 or step>=90 else .5 if step<60 else 1.
        for port,high in upper.items():u[port]=amplitude*high
        u,_=GP.adapt_inputs(n,u,enforce_domain=True);n.step(u);verify(n,gate);finite=finite and bool(np.isfinite(core.V).all())
        for name,i in idx.items():
            port=n.inputs.index(name);high=bounds[name][1]*core.input_connectivity[i,port]/core.g_m[i]
            largest=max(largest,-float(core.V[i]),float(core.V[i])-high)
    if not finite or largest>1e-10:raise AssertionError('Actual neural input-envelope/transient check failed')
    invalid=n.make_inputs();name=next(iter(bounds));invalid[n.inputs.index(name)]=bounds[name][1]+1.
    try:GP.adapt_inputs(n,invalid,enforce_domain=True);rejected=False
    except ValueError:rejected=True
    if not rejected:raise AssertionError('Out-of-domain current was silently accepted')
    from sparse_pairwise_eval import SparsePairwiseForward
    dense=copy.deepcopy(core);fast=SparsePairwiseForward(copy.deepcopy(core));rng=np.random.default_rng(1107);h=hashlib.sha256()
    for step in range(128):
        u=np.zeros(len(n.inputs))
        for port,high in upper.items():u[port]=rng.uniform(0,high)
        dense.forward(u);fast.forward(u)
        if not np.array_equal(dense.V,fast.V):raise AssertionError('Final composed dense/sparse output mismatch')
        h.update(dense.V.tobytes())
    return {'actual_base_net_step_checked':120,'finite':finite,'maximum_encoder_envelope_violation_mV':largest,'out_of_domain_rejected':rejected,'conditional_runner_domains':bounds,
            'domain_scope':'Analytically clipped runner envelope with gamma0,vest_prop0,afferent scale1; other ports zero in this test. Bounds need not be simultaneously physically attainable. No body-range measurement.',
            'dense_sparse128_exact':True,'voltage_trace_sha256':h.hexdigest()}


def coverage(n,condition):
    c=n.compiled;names=GP.population_names(n);active_pf=[i for i,(post,pre) in enumerate(zip(c.post,c.pre)) if names[post].startswith('MN_') and names[pre].startswith('PF_') and c.gmax[i]>0]
    rows=[]
    for joint,mechanics in [('hip_flexion','sagittal reference'),('knee_angle','sagittal reference; protective extension stops remain external setup'),('ankle_angle','sagittal reference plus declared foot restraint'),('hip_adduction','distributed muscle drive plus frontal restraint surrogate'),('hip_rotation','distributed muscle drive plus artificial rotation RIG'),('subtalar_angle','distributed muscle drive plus foot restraint surrogate'),('mtp_angle','no independent canonical MTP PF target; muscle-route coupling plus foot restraint')]:
        rows.append({'coordinates':[joint+'_'+s for s in ('r','l')],'neural_status':'PF-to-MN expression masked' if not active_pf else 'Distributed existing PF/muscle/afferent paths; independent joint target/control not certified','mechanical_contributions':mechanics})
    rows.append({'coordinates':['lumbar_extension','lumbar_bending','lumbar_rotation'],'neural_status':'Trunk MN/local afferents and external balance inputs; no separately certified3-axis PF loop','mechanical_contributions':'Existing lumbar RIG springs,passive forces,and source damping; not free balance'})
    rows.append({'coordinates':['pelvis_'+k for k in ('tx','ty','tz','tilt','list','rotation')],'neural_status':'External balance inputs may modulate muscles; floating-base control not proved','mechanical_contributions':'Existing horizontal/rotational/vertical apparatus remains a distinct assistance source'})
    return {'base_ledger':GP.gait_coverage_ledger(n),'condition':condition,'active_PF_to_MN_slots':len(active_pf),'coordinate_status':rows,'physical_control_rank_not_inferred_from_paths':True,'free3Dwalking_is_not_a_design_prerequisite':True,'restraint_passive_or_unsupported_status_explicit':True}


def compose(n,family,condition,anchor):
    from sparse_pairwise_eval import SparsePairwiseForward
    n.compiled=CB.core(n);c=n.compiled;before=states(c);names=GP.population_names(n);p=GP.default_profile(family)
    p['physical_input_domains_nA']={name:[0.,anchor['derived_source_upper_nA'][name.split('_')[0]]] for name in n.inputs if name.startswith(('Ia_','II_','Ib_'))}
    GP.apply_profile(n,p,family=family);central=CB.install(n,CB.default_profile())
    ids=[i for i,name in enumerate(names) if name.startswith('IBEXC_')];active=[i for i in ids if np.any(c.g_max_non[:,i])];gate=None
    if len(active)==len(ids):gate=IB.apply_coincidence(n,IB.make_policy(family,{name:.225 for name in IB.eligible_load_sources(n,family)},'Canonical full-stance runner AFF1.5*0.3*0.5 at F/Fmax1; engineering reference,not body-validated stance force'))
    gate_inventory={'allocated':len(ids),'active_output_gates':[names[i] for i in active],'status':'applied' if gate else 'Inactive/partial outputs preserved; coincidence inventory only','masked_gates':[names[i] for i in ids if i not in active]}
    n.compiled=SparsePairwiseForward(c);GP.apply_profile(n,p,family=family);final=verify(n,gate)
    if any(not np.array_equal(v,getattr(c,k)) for k,v in before.items()):raise AssertionError('Design composition changed initial dynamic state')
    gp_edges={(r['post_index'],r['pre_index']) for r in n.graded_fsa_ledger['modified_edges']};cb_edges={(r['post'],r['pre']) for g in central['groups'] for r in g['edges']}
    if gp_edges&cb_edges:raise AssertionError('Terminal and central budget selections overlap')
    pools=GP.probe_pools(n);encoder=GP.probe_encoder(n)
    pf=np.array([i for i,(post,pre) in enumerate(zip(n.compiled.post,n.compiled.pre)) if names[post].startswith('MN_') and names[pre].startswith('PF_') and n.compiled.gmax[i]>0],dtype=int)
    pfcheck=current_probe(n,pf);cbcheck=central_probe(n);gatecheck=None
    if gate:
        row=max(gate['rows'],key=lambda r:len(r['candidate']['loads']));gatecheck={'representative':row['name'],'truth':[IB.isolated(row['candidate'],phase,load) for phase,load in ((False,0),(True,0),(False,1),(True,1))],'release':[IB.isolated(row['candidate'],True,r,True) for r in (1.,1.5)]}
        if [r['final_fraction'] for r in gatecheck['truth']]!=[0.,0.,0.,1.]:raise AssertionError('Composed coincidence truth failed')
    runtime=domain_and_backend_probe(n,gate)
    return {'condition':condition,'family':family,'anchor':anchor,'neurons':len(c.V),'inputs':len(n.inputs),'initial_states_preserved':True,'final_consumption':verify(n,gate),
            'graded_profile':p,'graded_static_pool_checks':{'pools':pools['pools_checked'],'maximum_error_mV':pools['all_pool_static_error_max_mV'],'representatives':pools['representatives']},
            'encoder_phase_checks':encoder['results'],'central_design':central,'central_current_checks':cbcheck,'coincidence_inventory':gate_inventory,'coincidence_design':gate,'coincidence_checks':gatecheck,
            'PF_current_check':pfcheck,'actual_runtime_and_domains':runtime,'control_coverage':coverage(n,condition),'selected_GP_CB_conductance_intersection':0,'mechanics_changed':False,
            'role':'Gait neural pathway composition reviewer','body_bindings':0,'physics_steps':0,'training':False,'Hunt':False,'curriculum':False,'damping_repair_enabled':False,'tendon_compliance_enabled':False,'user_training_approval_required':True}


def run(family,out):
    files=[Path(__file__),Path(GP.__file__),Path(CB.__file__),Path(IB.__file__),Path(MC.__file__),HERE.parent/'synergy_basis.npz',HERE.parent/'runner.py']
    frozen={str(p):CB.sha(p) for p in files};callback=(lambda n:MC.apply_canonical_right(n,MC.default_policy())) if family=='syn6' else None
    os.environ['AARL_NET']=family;sys.path.insert(0,str(HERE.parent.parent/'robot_torso/bpa_neural_coordination'))
    import gait_standing_family as F
    default_flags={k:F.P.G[k] for k in ('ii_set_pose','posture_to_mn')}
    n,anchor,P=CB.build_saved(family,projection_callback=callback);standing=compose(n,family,'saved_standing_neural_design_candidate',anchor)
    import mujoco
    # The standing constructor intentionally overrides these globals. Restore
    # source defaults only for this separately named canonical-expression fixture.
    for key,value in default_flags.items():P.G[key]=value
    model=mujoco.MjModel.from_xml_path(str(F.R.MODEL));acts=[model.actuator(i).name for i in range(model.nu)];template=F.B.build(acts,dt=P.DT,interleg=True)
    # Verify original builder anchor BEFORE any new projection/design changes.
    clone=copy.deepcopy(template.compiled);rng=np.random.default_rng(1101);h=hashlib.sha256()
    for _ in range(128):
        u=rng.uniform(0.,.2,len(template.inputs));u[template.input_index('DRIVE')]=1.2;clone.forward(u);h.update(clone.V.tobytes())
    if h.hexdigest()!=ANCHORS[family]:raise AssertionError('Canonical default neural anchor changed')
    projection=MC.apply_canonical_right(template,MC.default_policy()) if family=='syn6' else F.tie_weights(template)
    expression=compose(template,family,'NEW_canonical_template_neural_expression_fixture',{'derived_source_upper_nA':anchor['derived_source_upper_nA'],'source_saved_anchor_validated_separately':True,'saved_gain_suppression_not_applied':True,'projection':projection,'default_neural_anchor_exact':True})
    if frozen!={str(p):CB.sha(p) for p in files}:raise RuntimeError('Neural composition source changed during audit')
    result={'role':'Gait neural pathway composition reviewer','family':family,'source_sha256':frozen,'saved_standing':standing,'expression_fixture':expression,'complete':'Scoped electrical composition,unit/current/encoding/backend checks only','remaining':'No full-neural/body control performance or training result. Ben explicit approval required before any training.','global_release':False,'mechanics_unchanged':True,'no_damping_or_tendon_application':True}
    out.mkdir(parents=True,exist_ok=True);target=out/(family+'_COMPOSITION.json')
    if target.exists():raise FileExistsError('Preserve previous composition record: '+str(target))
    target.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({'family':family,'standing_gates':standing['coincidence_inventory']['status'],'standing_PF_slots':standing['PF_current_check']['active_slots'],'expression_PF_slots':expression['PF_current_check']['active_slots'],'central_groups':expression['central_current_checks']['groups'],'dense_sparse_exact':True,'physics_steps':0,'output':str(target)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=('s3k','w2lvar','syn6'),required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.family,a.out)
