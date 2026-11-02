"""Explicit Hunt1033 neural composition; historical defaults are unchanged.

Use apply_to_network after a caller's gain transforms and before body binding.
The standalone builder retains the real PT object; unbound electrical tests
temporarily detach it and never claim a full PT/body validation.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ.setdefault(k,'1')
import argparse,ast,copy,hashlib,json,sys
from pathlib import Path
from contextlib import contextmanager
import xml.etree.ElementTree as ET
import numpy as np
HERE=Path(__file__).resolve().parent;S=HERE.parent;H=S/'hunt_training'
sys.path[:0]=[str(HERE),str(H),str(S),str(S/'w2l_cpg'),str(S/'functional_gait2392')]
import controller as C,graded_profile as GP,hunt_graded_profile as LOCAL,hunt_timing_profile as TIMING,ibexc_coincidence as IB
from central_afferent_budget import clone
from sparse_pairwise_eval import SparsePairwiseForward
SCHEMA='aarl_integrated_hunt1033_v1'
DEFAULT_CONDITION=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/aligned_support_ladder_v1/packaged_configuration.json'
MODEL=S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def condition_options(cfg):
    source=HERE/'hunt_pathway_audit.py';tree=ast.parse(source.read_text());node=next(q for q in tree.body if isinstance(q,ast.FunctionDef) and q.name=='options_from_physical_condition');env={'H':H}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),env)
    return env['options_from_physical_condition'](cfg)


@contextmanager
def unbound_electrical(net):
    balance=net.__dict__.pop('balance',None)
    try:yield
    finally:
        if balance is not None:net.balance=balance


def make_profile(condition=DEFAULT_CONDITION,timing_mode='function_equivalent',budgets=None):
    condition=Path(condition).resolve();cfg=json.loads(condition.read_text());w=Path(cfg['weights_directory'])
    files=[condition,MODEL,Path(__file__),HERE/'hunt_pathway_audit.py',HERE/'hunt_graded_profile.py',HERE/'hunt_timing_profile.py',HERE/'ibexc_coincidence.py',HERE/'graded_profile.py',H/'controller.py',H/'balance.py',H/'runs/step_A/winner.json',H/'runs/step_B/winner.json',H/'runs/step_B_left/winner.json',H/'runs/targets/targets.npz',H/'runs/targets/metadata.json',C.F.OUT/'routing_audit.json',C.F.OUT/'routing_samples.npz']
    if (S/'joint_pf_weights.json').exists():files.append(S/'joint_pf_weights.json')
    files+=list(w.glob('*.json'))
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and Path(path).suffix=='.py' and Path(path).resolve().is_relative_to(S):files.append(Path(path).resolve())
    files=list(dict.fromkeys(files))
    local=GP.default_profile('w2lvar');bounds={'Ia':C.P.AFF['ia_gain']*(1+C.P.MOD['ia']),'II':C.P.AFF['i0_ii']+C.P.AFF['ii_gain']*(1+C.P.MOD['ii']),'Ib':.5*C.P.AFF['ib_gain']*1.5}
    acts=[q.attrib['name'] for q in ET.parse(MODEL).getroot().find('actuator')]
    local['physical_input_domains_nA']={f'{kind}_{act}':[0.,upper] for kind,upper in bounds.items() for act in acts}
    timing=TIMING.default_profile(mode=timing_mode,relay_timing='fast_role_budget')
    if budgets is not None:
        if timing_mode!='custom_budgets' or set(budgets)!={'E','F'}:raise ValueError('Two explicit shared E/F budgets require custom_budgets mode')
        timing['per_RG_budget_nA']={k:float(v) for k,v in budgets.items()}
    return {'schema':SCHEMA,'condition_file':str(condition),'source_sha256':{str(q):sha(q) for q in files},'local_profile':local,'timing_profile':timing,'Ib_reference_current_nA':.225,'heel_hypothesis':'Preserve actual source contact paths; B requires active inherited heel-IN routes and no stacked direct reset','parameter_checkpoint':None,'approval_scope':'Neural composition/electrical verification. Body and training require the separately authorized bounded campaign.'}


def load_profile(value):
    p=json.loads(Path(value).read_text()) if isinstance(value,(str,Path)) else copy.deepcopy(value)
    GP.exact_keys(p,('schema','condition_file','source_sha256','local_profile','timing_profile','Ib_reference_current_nA','heel_hypothesis','parameter_checkpoint','approval_scope'),'integrated Hunt profile')
    if p['schema']!=SCHEMA:raise ValueError('Wrong integrated schema')
    for path,digest in p['source_sha256'].items():
        if sha(path)!=digest:raise ValueError('Integrated source changed; explicitly rebind/review profile: '+path)
    if p['parameter_checkpoint'] is not None:
        if set(p['parameter_checkpoint'])!={'path','sha256'} or sha(p['parameter_checkpoint']['path'])!=p['parameter_checkpoint']['sha256']:raise ValueError('Gain checkpoint changed')
    return p


def apply_to_network(net,profile=None,protocol_stage=None):
    if profile is None:return {'applied':False,'default_noop':True}
    if hasattr(net,'integrated_fsa_design'):raise ValueError('Composition already installed; do not reinitialize fitted gains every tick')
    p=load_profile(profile);stage=int(protocol_stage or json.loads(Path(p['condition_file']).read_text())['stage'])
    if type(net) is not C.HuntNetwork:raise ValueError('Actual current HuntNetwork required')
    if stage==1:
        if len(net.compiled.V)!=1025 or hasattr(net,'balance') or net.options['use_timing_afferents']:raise ValueError('Expected explicit deafferented1025/noPT protocol projection')
    else:
        IB.hunt_identity(net)
        if len(net.compiled.V)!=1033:raise ValueError('Expected full1033 architecture outside deafferented protocol')
    parameter_mapping=None;gate=None
    with unbound_electrical(net):
        net.apply_weights(1.)
        if stage!=1:
            LOCAL.install(net,p['local_profile']);TIMING.install(net,p['timing_profile'])
            policy=IB.make_policy('Hunt_PT',{q:p['Ib_reference_current_nA'] for q in IB.eligible_load_sources(net,'Hunt_PT')},'Canonical AFF.ib_gain1.5*full-stance0.3*0.5*F/Fmax1=.225nA; engineering reference,not a measured stance input')
            gate=IB.apply_coincidence(net,policy)
        if p['parameter_checkpoint'] is not None:
            with np.load(p['parameter_checkpoint']['path']) as a:
                oldnames=a['names'].tolist();lookup={(oldnames[int(pre)],oldnames[int(post)]):i for i,(pre,post) in enumerate(zip(a['pre'],a['post']))};current=[(net.names[int(pre)],net.names[int(post)]) for pre,post in zip(net.compiled.pre,net.compiled.post)]
                if any(k not in lookup for k in current):raise ValueError('Unexpected new endpoint in protocol projection')
                omitted=set(lookup)-set(current)
                if omitted and (stage!=1 or any(not any(name.startswith(('PT_BAL_','PT_IN_')) for name in pair) for pair in omitted)):raise ValueError('Only explicit omitted PT endpoints allowed')
                take=np.array([lookup[k] for k in current]);parameter_mapping={'mapped_current_slots':len(current),'omitted_PT_slots':len(omitted),'binding':'Exact source/target names,not assumed equal sparse indices'}
                for field in ('pf_gain','ff_gain','reflex_gain','bi_gain'):
                    if not np.isfinite(a[field]).all() or np.any(a[field]<0):raise ValueError('Invalid fitted gain checkpoint')
                    if stage!=1 or field!='reflex_gain':getattr(net,field)[:]=a[field][take]
        if stage==1:net.reflex_gain[:]=0.
        net.apply_weights(1.);net.sync()
    c=net.compiled;names=net.names;heel=[]
    for i,(post,pre) in enumerate(zip(c.post,c.pre)):
        if names[pre].startswith('HEEL_') and names[post].startswith(('RG_','InE_','InF_')):heel.append({'source':names[pre],'target':names[post],'conductance_uS':float(c.gmax[i]),'reversal_mV':float(c.reversal[i])})
    direct=[q for q in heel if q['target'].startswith('RG_') and q['conductance_uS']>0];masked=[q for q in heel if q['target'].startswith(('InE_','InF_')) and q['conductance_uS']==0]
    # Preserve existing contact policy: no inherited masked path is reactivated.
    net.integrated_fsa_design={'profile':p,'protocol_stage':stage,'local':LOCAL.serializable_ledger(net) if stage!=1 else None,'timing':getattr(net,'hunt_timing_normalized_design',None),'coincidence':gate,'parameter_checkpoint_mapping':parameter_mapping,'protocol_omissions':['Local afferent outputs remain zero','Timing afferents disabled','Active Ib motor-gate checks omitted','PT omitted by source protocol'] if stage==1 else [],'heel':{'B_applied':False,'reason':'Existing direct heel-RG reset would stack and inherited heel-IN routes are masked' if direct or masked else 'No eligible active source heel-IN contract','actual_paths':heel},'PT_object_retained':hasattr(net,'balance'),'PT_bound':False,'no_physics_changes':True}
    verify(net,full=True)
    return net.integrated_fsa_design


def verify(net,full=False,require_bound_PT=False):
    if not hasattr(net,'integrated_fsa_design'):return {'applied':False}
    c=net.compiled;stage=net.integrated_fsa_design['protocol_stage'];ids=net.hunt_graded_fsa['protected_indices'] if stage!=1 else net.reflex_edges
    if not np.array_equal(c.gmax[ids],net.reflex_gain[ids]) or not np.array_equal(c.gmax,c.c.g_max_non[c.post,c.pre]):raise ValueError('Actual authoritative gains/cache diverged')
    result={'applied':True,'authoritative_cache_match':True}
    if stage==1:
        if np.any(net.reflex_gain) or np.any(c.gmax[net.reflex_edges]) or net.options['use_timing_afferents'] or hasattr(net,'balance'):raise ValueError('Deafferented protocol mask was reactivated')
        result['protocol_omissions']=net.integrated_fsa_design['protocol_omissions'];result['deafferented_masks_verified']=True
    elif full:
        result.update(local=LOCAL.verify_runtime(net),timing=TIMING.verify_runtime(net),coincidence=IB.verify_coincidence(net))
        pair=IB.audit_hunt_gain_pairing(net)
        if any(v['mismatched_pairs'] for v in pair['buffers'].values()):raise ValueError('Authoritative bilateral gain pairing changed')
        result['bilateral_gain_pairing']=True
    else:
        TIMING.verify_runtime(net);IB.verify_coincidence(net)
    bound=hasattr(net,'balance') and hasattr(net.balance,'calibration') and hasattr(net.balance,'edges')
    if require_bound_PT and not bound:raise ValueError('Actual PT mechanical calibration is not bound')
    result['PT_bound']=bound
    if bound:
        b=net.balance
        if not np.isfinite(b.basis).all() or not np.isfinite(b.response).all() or not np.isfinite(b.gain).all():raise ValueError('Nonfinite actual PT calibration')
        if b.options.get('bilateral_recruitment',False):
            names=net.act_order;mirror=np.array([names.index(a[:-2]+('_l' if a.endswith('_r') else '_r')) if a.endswith(('_r','_l')) else i for i,a in enumerate(names)])
            if not all(np.array_equal(b.basis[k],b.basis[k][mirror]) for k in (0,1)) or not np.array_equal(b.basis[2],b.basis[3][mirror]):raise ValueError('Actual PT recruitment is not tied')
        result['PT_active_directions']=np.any(abs(b.basis)>1e-12,axis=1).tolist();result['PT_calibration']=b.calibration
        net.integrated_fsa_design['PT_bound']=True;net.integrated_fsa_design['PT_calibration']=b.calibration
    return result


def build(acts,dt=.002,interleg=False,options=None,*,profile=None,protocol_stage=None):
    net=C.build(acts,dt=dt,interleg=interleg,options=options);apply_to_network(net,profile,protocol_stage);return net


def fresh_check(profile,out):
    p=load_profile(profile);cfg=json.loads(Path(p['condition_file']).read_text());acts=[q.attrib['name'] for q in ET.parse(MODEL).getroot().find('actuator')];n=build(acts,options=condition_options(cfg),profile=p);c=n.compiled
    with np.load(Path(profile).parent/'COMPILED_GRAPH.npz') as saved:
        for key in saved.files:
            actual=getattr(c,key) if key in ('pre','post','gmax','reversal','elo','erange') else getattr(c.c,key)
            if not np.array_equal(actual,saved[key]):raise AssertionError('Fresh-process final graph mismatch: '+key)
    result={'fresh_process_parameter_graph_exact':True,'verification':verify(n,full=True),'neurons':len(c.V),'body_steps':0,'PT_bound':False,'training':False};Path(out).write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))


def snapshot_parameters(net,out,profile):
    """Explicit checkpoint of fitted authoritative gains; never reset per step."""
    out=Path(out);p=copy.deepcopy(profile);np.savez_compressed(out,names=np.array(net.names),pre=net.compiled.pre,post=net.compiled.post,**{q:getattr(net,q) for q in ('pf_gain','ff_gain','reflex_gain','bi_gain')});p['parameter_checkpoint']={'path':str(out.resolve()),'sha256':sha(out)};return p


def electrical_checks(net):
    checks=[];rng=np.random.default_rng(110209)
    with unbound_electrical(net):
        for command in (0.,.5,1.,0.):
            for j in range(32):
                u=np.zeros(len(net.inputs));u[net.ports['DRIVE']]=2.5*command
                for name,index in net.ports.items():
                    if name.startswith('HUNT_AF_'):u[index]=float((j%8)/8)
                net.step(u);checks.append(verify(net,full=(j==31)))
        net.apply_weights(1.);net.sync();source=clone(net.compiled);zero=clone(source)
        positive=next(i for i in net.reflex_edges if zero.gmax[i]>0);zero.gmax[positive]=0.;zero.c.g_max_non[zero.post[positive],zero.pre[positive]]=0.
        rebuilt=SparsePairwiseForward(copy.deepcopy(zero.c));dense=copy.deepcopy(zero.c);max_rebuild=0.;max_dense=0.
        for j in range(64):
            u=rng.uniform(0.,.1,len(net.inputs));zero.forward(u);rebuilt.forward(u);dense.forward(u)
            max_rebuild=max(max_rebuild,float(np.max(abs(zero.V-rebuilt.V))));max_dense=max(max_dense,float(np.max(abs(zero.V-dense.V))))
        if max_rebuild!=0 or max_dense!=0:raise AssertionError('Zero/delete/rebuilt/dense electrical mismatch')
        on=clone(source);off=clone(source);dst=int(on.post[positive]);pre=int(on.pre[positive]);on.c.V[pre]=on.elo[positive]+on.erange[positive];off.c.V[pre]=on.c.V[pre];on.c.V[dst]=off.c.V[dst]=2.5;off.gmax[positive]=0.;off.c.g_max_non[dst,pre]=0.
        expected=float(on.gmax[positive]*(on.reversal[positive]-2.5));u=np.zeros(len(net.inputs));on.forward(u);off.forward(u);actual=float((on.V[dst]-off.V[dst])/on.c.time_factor_membrane[dst]);error=abs(actual-expected)
        if error>1e-10:raise AssertionError('Actual motor current mismatch')
    return {'actual_command_runtime_steps':128,'all_commands_verified':[0.,.5,1.,0.],'finite':bool(np.isfinite(net.compiled.V).all()),'zero_delete_rebuild_steps':64,'zero_delete_rebuild_max_error_mV':max_rebuild,'dense_reference_max_error_mV':max_dense,'zero_edge_scope':'Finite electrical inputs; fixed finalized gains. Rebuilt wrapper removes zero slots. Live Hunt keeps allocated slots for legitimate future command changes.','one_synapse_actual_motor_current_error_nA':error,'deleted_edge':{'source':net.names[pre],'target':net.names[dst]},'PT_scope':'PT object retained in returned candidate; temporarily detached only for these unbound electrical checks. Actual bound PT verification must follow.','body_steps':0,'training':False}


def default_endhook_parity(cfg):
    """Execute only physical_eval's neural constructor,never its module/main."""
    tree=ast.parse((H/'physical_eval.py').read_text());source=next(q for q in tree.body if isinstance(q,ast.FunctionDef) and q.name=='build')
    class WithoutHook(ast.NodeTransformer):
        def visit_Assign(self,node):
            return None if any(isinstance(t,ast.Name) and t.id=='integrated_check' for t in node.targets) else self.generic_visit(node)
        def visit_If(self,node):
            marker=ast.dump(node.test)
            return None if 'integrated_fsa_profile' in marker or "id='integrated_check'" in marker else self.generic_visit(node)
    source_old=WithoutHook().visit(copy.deepcopy(source));ast.fix_missing_locations(source_old);made=[];acts=[q.attrib['name'] for q in ET.parse(MODEL).getroot().find('actuator')]
    for node in (source_old,source):
        env={'H':H,'C':C,'np':np,'cfg':{k:v for k,v in cfg.items() if k!='integrated_fsa_profile'},'stage':cfg['stage'],'capture':{},'trace_cells':[],'event_amplitude':0.}
        exec(compile(ast.Module(body=[copy.deepcopy(node)],type_ignores=[]),str(H/'physical_eval.py'),'exec'),env);made.append(env['build'](acts,.002,True))
    if GP.field_digest(made[0].compiled.c)!=GP.field_digest(made[1].compiled.c) or any(not np.array_equal(getattr(made[0],k),getattr(made[1],k)) for k in ('pf_gain','ff_gain','reflex_gain','bi_gain')):raise AssertionError('Default physical_eval constructor changed')
    return {'default_neural_constructor_exact_to_pre_hook_AST':True,'body_steps':0,'scope':'Actual unchanged-condition build function with/without only the new opt-in blocks; no whole-body trajectory replay'}


def protocol_setup_checks(cfg,profile):
    tree=ast.parse((H/'physical_eval.py').read_text());node=next(q for q in tree.body if isinstance(q,ast.FunctionDef) and q.name=='build');acts=[q.attrib['name'] for q in ET.parse(MODEL).getroot().find('actuator')];rows=[]
    for stage in (1,2,3):
        env={'H':H,'C':C,'np':np,'cfg':{**cfg,'stage':stage,'integrated_fsa_profile':profile},'stage':stage,'capture':{},'trace_cells':[],'event_amplitude':0.}
        exec(compile(ast.Module(body=[copy.deepcopy(node)],type_ignores=[]),str(H/'physical_eval.py'),'exec'),env);n=env['build'](acts,.002,stage!=1);rows.append({'stage':stage,'neurons':len(n.compiled.V),'verification':verify(n,full=True),'PT_object_present':hasattr(n,'balance'),'timing_enabled':n.options['use_timing_afferents'],'protocol_omissions':n.integrated_fsa_design['protocol_omissions']})
    return {'actual_physical_eval_builder_stages':rows,'neural_steps':0,'body_steps':0}


def run(out,condition=DEFAULT_CONDITION,timing_mode='function_equivalent',budgets=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);p=make_profile(condition,timing_mode,budgets);cfg=json.loads(Path(condition).read_text());acts=[q.attrib['name'] for q in ET.parse(MODEL).getroot().find('actuator')]
    baseline=C.build(acts,options=condition_options(cfg))
    with unbound_electrical(baseline):baseline.apply_weights(1.)
    baseline_hash=GP.field_digest(baseline.compiled.c);apply_to_network(baseline,None)
    if GP.field_digest(baseline.compiled.c)!=baseline_hash:raise AssertionError('Default None changed historical network')
    net=build(acts,options=condition_options(cfg),profile=p);checks=electrical_checks(net)
    baseline_again=C.build(acts,options=condition_options(cfg))
    with unbound_electrical(baseline_again):baseline_again.apply_weights(1.)
    if GP.field_digest(baseline_again.compiled.c)!=baseline_hash:raise AssertionError('Candidate leaked into default builder')
    checks['None_default_parameter_exact']=True;checks['fresh_default_after_candidate_parameter_exact']=True
    checks['physical_eval_default_hook_parity']=default_endhook_parity(cfg)
    p=snapshot_parameters(net,out/'GAIN_PARAMETERS.npz',p);(out/'PROFILE.json').write_text(json.dumps(p,indent=2),encoding='utf-8');checks['protocol_stage_setup']=protocol_setup_checks(cfg,p)
    fresh=build(acts,options=condition_options(cfg),profile=p)
    with unbound_electrical(net),unbound_electrical(fresh):net.apply_weights(1.);fresh.apply_weights(1.)
    if GP.field_digest(net.compiled.c)!=GP.field_digest(fresh.compiled.c):raise AssertionError('Serialized candidate rebuild differs')
    c=fresh.compiled;arrays={q:np.asarray(getattr(c.c,q)) for q in GP.FIELDS+GP.ION_FIELDS if hasattr(c.c,q)};arrays.update(pre=c.pre,post=c.post,gmax=c.gmax,reversal=c.reversal,elo=c.elo,erange=c.erange);np.savez_compressed(out/'COMPILED_GRAPH.npz',**arrays)
    record={'schema':SCHEMA,'machine':'easteregg2','condition':str(condition),'neurons':len(c.V),'inputs':len(fresh.inputs),'names':fresh.names,'input_order':fresh.inputs,'checks':checks,'serialized_parameter_rebuild_exact':True,'final_verification':verify(fresh,full=True),'composition':fresh.integrated_fsa_design,'PT_full_validation_pending_static_bind':True,'physics_changed':False,'body_trials':0,'training':False,'source_sha256':p['source_sha256']};(out/'INTEGRATED_REVIEW.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps({'neurons':len(c.V),'serialized_rebuild_exact':True,'zero_delete_rebuild_error':checks['zero_delete_rebuild_max_error_mV'],'profile':str(out/'PROFILE.json'),'PT_bound':False}))


def static_bind_preflight(profile,out):
    """Actual physical_eval initialization,stopped after bind; isolated CLI only."""
    if __name__!='__main__':raise ValueError('Run this preflight in its own process; it intercepts initialization')
    import runpy,mujoco
    p=load_profile(profile);out=Path(out);cfg=json.loads(Path(p['condition_file']).read_text());cfg['integrated_fsa_profile']=str(Path(profile).resolve());cfg['export_graph']=False
    config_path=out.parent/'STATIC_BIND_CONFIG.json';config_path.write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    original_bind=C.HuntNetwork.bind_mujoco;original_step=mujoco.mj_step;original_argv=sys.argv[:]
    class Finished(Exception):pass
    def no_body_step(*args,**kwargs):raise AssertionError('Static preflight attempted a body step')
    def bind(net,m,d):
        original_bind(net,m,d);net.apply_weights(1.);net.sync();verified=verify(net,full=True,require_bound_PT=True);c=net.compiled
        chosen=np.array([i for i,(pre,post) in enumerate(zip(c.pre,c.post)) if net.names[pre].startswith(('PT_BAL_','PT_IN_')) and net.names[post].startswith('MN_') and c.gmax[i]>0],dtype=int)
        if not len(chosen):raise AssertionError('Actual bound PT has no nonzero motor connections')
        targets=np.unique(c.post[chosen]);pres=np.unique(c.pre[chosen]);expected=np.zeros(len(c.V))
        for i in chosen:expected[c.post[i]]+=c.gmax[i]*np.clip((3.-c.elo[i])/c.erange[i],0.,1.)*(c.reversal[i]-2.5)
        # A caller may install an instance forward closure. Copying __dict__
        # does not rebind that closure. Exercise the actual bound dispatcher,
        # resetting neural state/weights between counterfactuals instead.
        saved_state=copy.deepcopy(c.c.__dict__);saved_g=c.gmax.copy();voltages=[]
        for enabled in (True,False):
            c.c.__dict__.clear();c.c.__dict__.update(copy.deepcopy(saved_state));c.gmax[:]=saved_g
            if not enabled:c.gmax[chosen]=0.
            net.sync();c.c.V[pres]=3.;c.c.V[targets]=2.5;c.forward(np.zeros(len(net.inputs)));voltages.append(c.c.V[targets].copy())
        actual=(voltages[0]-voltages[1])/c.c.time_factor_membrane[targets];error=float(np.max(abs(actual-expected[targets])))
        c.c.__dict__.clear();c.c.__dict__.update(saved_state);c.gmax[:]=saved_g;net.sync();verify(net,full=True,require_bound_PT=True)
        if error>1e-10:
            failure={'maximum_error_nA':error,'actual':actual.tolist(),'expected':expected[targets].tolist(),'targets':[net.names[i] for i in targets],'backend':type(c).__name__,'configured_backend':cfg.get('backend','sparse'),'active_PT_edges':len(chosen),'actual_forward_module':getattr(c.forward,'__module__',None),'actual_forward_qualname':getattr(c.forward,'__qualname__',None)}
            out.with_name('STATIC_BIND_CURRENT_FAILURE.json').write_text(json.dumps(failure,indent=2),encoding='utf-8');raise AssertionError('Bound PT compiled-current mismatch: '+str(error))
        knees={side:float(np.degrees(d.qpos[m.joint('knee_angle_'+side).qposadr[0]])) for side in ('r','l')}
        if max(knees.values())>0.:raise AssertionError('Initial pose violates zero-hyperextension requirement')
        physical_hashes={k:hashlib.sha256(np.ascontiguousarray(getattr(m,k)).tobytes()).hexdigest() for k in ('body_mass','body_inertia','dof_armature','dof_damping','jnt_stiffness','jnt_axis','eq_data','actuator_gainprm','actuator_biasprm')}
        result={'actual_physical_eval_source':str(H/'physical_eval.py'),'source_sha256':sha(H/'physical_eval.py'),'profile':str(profile),'configuration':cfg,'verification':verified,'PT_motor_connections_active':len(chosen),'actual_PT_motor_current_error_nA':error,'actual_forward_module':getattr(c.forward,'__module__',None),'actual_forward_qualname':getattr(c.forward,'__qualname__',None),'neural_state_restored_between_counterfactuals':True,'initial_knee_angles_degrees':knees,'compiled_physics_parameter_sha256':physical_hashes,'physics_policy':'Exact existing physical_eval configuration/setup; no new tendon,damping,contact,or dimension change','PT_scope':'Actual original bind_mujoco and fixed-pose PT mechanical calibration,followed by real compiled electrical current subtraction; not a standing/gait/stability trial','body_steps':0,'neural_current_probe_steps':2,'training':False}
        np.savez_compressed(out.parent/'COMPILED_BOUND_GRAPH.npz',pre=c.pre,post=c.post,gmax=c.gmax,reversal=c.reversal,elo=c.elo,erange=c.erange,PT_basis=net.balance.basis,PT_response=net.balance.response)
        out.write_text(json.dumps(result,indent=2),encoding='utf-8');raise Finished
    C.HuntNetwork.bind_mujoco=bind;mujoco.mj_step=no_body_step;sys.argv=[str(H/'physical_eval.py'),str(config_path)]
    try:
        try:runpy.run_path(str(H/'physical_eval.py'),run_name='__main__')
        except Finished:pass
        else:raise AssertionError('Actual bind preflight was not intercepted')
    finally:C.HuntNetwork.bind_mujoco=original_bind;mujoco.mj_step=original_step;sys.argv=original_argv
    print(json.dumps({'static_PT_bind_passed':True,'body_steps':0,'output':str(out)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--condition',type=Path,default=DEFAULT_CONDITION);parser.add_argument('--verify-profile',type=Path);parser.add_argument('--static-bind-profile',type=Path);parser.add_argument('--timing-mode',choices=('function_equivalent','canonical_right_role_budget','custom_budgets'),default='function_equivalent');parser.add_argument('--E-budget',type=float);parser.add_argument('--F-budget',type=float);a=parser.parse_args();budgets=None if a.E_budget is None and a.F_budget is None else {'E':a.E_budget,'F':a.F_budget}
    if a.static_bind_profile:static_bind_preflight(a.static_bind_profile,a.out)
    elif a.verify_profile:fresh_check(a.verify_profile,a.out)
    else:run(a.out,a.condition,a.timing_mode,budgets)
