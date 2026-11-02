"""Read-only family dossiers and isolated actual-source edge/operator probes.

No optimizer, body integration, source modification, or campaign release.
"""
import argparse, hashlib, json, math, os, sys
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
EXT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')
sys.path[:0]=[str(HERE.parent/'w2l_cpg'),str(HERE.parent),str(HERE.parent.parent)]
import numpy as np

def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def evidence(p,role): return dict(path=str(p),sha256=sha(p),role=role)
def write(p,z): Path(p).write_text(json.dumps(z,indent=2),encoding='utf-8')
def mirror(n):
    if n.endswith('_l'):return n[:-2]+'_r'
    if n.endswith('_r'):return n[:-2]+'_l'
    if n.endswith('_L'):return n[:-2]+'_R'
    if n.endswith('_R'):return n[:-2]+'_L'
    return n
def check(passed,text):return dict(passed=bool(passed),evidence=text)

def caller_context(P):
    keys=('gamma_coact','vest_prop','heel_rge','toe_rge','ib_rge','contact_onset','contra_swing','pm_gain','pm_aff','ankle_post_walk_trim')
    return dict(AFF=dict(P.AFF),MOD=dict(P.MOD),G={k:float(P.G[k]) for k in keys})

def validate_caller_input(values,names,binding,condition,context):
    """Strict opt-in physical-domain adapter; never clips or rescales commands."""
    if context!=binding['caller_parameter_context']:raise ValueError('Caller gains/context changed; derive and verify a new physical domain')
    x=np.asarray(values,float)
    if x.shape!=(len(names),) or not np.isfinite(x).all():raise ValueError('Caller shape/nonfinite')
    if names!=binding['input_order']:raise ValueError('Caller port ordering changed')
    bounds=binding['conditions'][condition]['domains_nA']
    for i,name in enumerate(names):
        lo,hi=bounds[name]
        if x[i]<lo-1e-12 or x[i]>hi+1e-12:raise ValueError('Caller outside declared condition: '+name)
    return x

class CallerDomainForward:
    """Validation-only wrapper: dynamics and complete state remain the core's."""
    def __init__(self,inner,binding,condition,context_provider):
        self.inner=inner;self.binding=binding;self.condition=condition;self.context_provider=context_provider
        self.binding_hash=hashlib.sha256(json.dumps(binding,sort_keys=True).encode()).hexdigest()
    def __getattr__(self,key):return getattr(self.inner,key)
    def forward(self,x):return self.inner.forward(validate_caller_input(x,self.binding['input_order'],self.binding,self.condition,self.context_provider()))
    def export_state(self):
        s=self.inner.export_state();s.update(caller_binding_hash=np.array(self.binding_hash),caller_condition=np.array(self.condition));return s
    def import_state(self,s):
        if str(s['caller_binding_hash'])!=self.binding_hash or str(s['caller_condition'])!=self.condition:raise ValueError('Caller state/binding mismatch')
        self.inner.import_state(s)

def derive_caller_domains(n,F,settings,output):
    """All382 bounds from actual formulas and saved static calibration, not maxima."""
    resultp=EXT/'spiking/Gait_dynamic_readout_v4/physics/syn6_central_v4_8p0.result.json'
    freshp=resultp.with_name('syn6_central_v4_8p0_fresh.result.json');r=read(resultp);fresh=read(freshp)
    act=np.asarray(r['standing_activation_targets']);assert len(act)==92 and np.array_equal(act,fresh['standing_activation_targets']) and np.isfinite(act).all() and np.all((act>=0)&(act<=1))
    assert r['standing_solver_diagnostics']['bilateral'] and r['standing_solver_diagnostics']['free_parameters']==46
    context=caller_context(F.P);assert context['G']['gamma_coact']==context['G']['vest_prop']==0
    for k in ('heel_rge','toe_rge','ib_rge','contact_onset','contra_swing','pm_gain','pm_aff'):assert context['G'][k]==0
    assert context['G']['ankle_post_walk_trim']==1
    bound=dict(Ia=float(F.P.AFF['ia_gain']*(1+F.P.MOD['ia'])),II=float(F.P.AFF['i0_ii']+F.P.AFF['ii_gain']*(1+F.P.MOD['ii'])),Ib=float(.5*F.P.AFF['ib_gain']*1.5))
    assert np.allclose([bound[k] for k in ('Ia','II','Ib')],[2.88,2.92,1.125],rtol=0,atol=1e-14)
    model=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));acts=[model.actuator(i).name for i in range(model.nu)]
    assert sha(F.R.MODEL)==read(EXT/'spiking/Gait_dynamic_readout_v4/TIMING_CANDIDATE.json')['source_model_sha256']
    posture={name:float(6*a*settings['standing_bias_scale']) for name,a in zip(acts,act)}
    standing={};deaff={};recipes={}
    for name in n.inputs:
        prefix=name.split('_')[0]
        if prefix in bound:
            standing[name]=[0.,bound[prefix]];deaff[name]=[0.,0.];recipes[name]='Runner rectification and bounded normalized mechanics; caller muscle_afferent_scale1'
        elif name.startswith('POST_'):
            a=name[5:];standing[name]=[posture[a],posture[a]];deaff[name]=[0.,0.];recipes[name]='6*saved validated bilateral static_activation*fixed standing_bias_scale; stand_frac1,ankletrim1'
        elif name=='DRIVE':standing[name]=[0.,0.];deaff[name]=[2.5,2.5];recipes[name]='--stand-eval schedule versus explicit deafferented electrical fixture'
        elif name=='POSTURE':standing[name]=[1.,1.];deaff[name]=[0.,0.];recipes[name]='--stand-eval schedule versus explicit deafferented electrical fixture'
        elif name in ('BAL_PF','BAL_DF','BAL_TRK_FLX','BAL_TRK_EXT','BAL_LAT_R','BAL_LAT_L'):
            standing[name]=deaff[name]=[0.,0.];recipes[name]='--reflex-stand zeros after feedback computation'
        elif name.startswith(('HEEL_c_','TOE_c_','LOAD_c_')):
            standing[name]=deaff[name]=[0.,0.];recipes[name]='min(Fnormal/(fraction*BW),1.5)*selected zero source gain; contact_onset and contra_swing0'
        else:raise ValueError('No source-formula domain for '+name)
    muscle_domains={name:standing[name] for name in n.inputs if name.startswith(('Ia_','II_','Ib_'))};assert len(muscle_domains)==276 and len(standing)==382
    ledger=dict(schema='aarl_canonical_central_caller_domain_v1',machine='easteregg2',units='nA',input_order=list(n.inputs),caller_parameter_context=context,mathematical_muscle_maxima=bound,normalized_mechanics=dict(length=[-1,1],velocity=[-1,1],force_fraction_Fmax=[0,1.5]),standing_bias_scale_interval=[settings['standing_bias_scale']]*2,muscle_afferent_scale_interval=[1.,1.],gamma_coact_interval=[0.,0.],vest_prop_interval=[0.,0.],conditions=dict(standing=dict(flags=['--stand-eval8','--reflex-stand'],domains_nA=standing),deaff=dict(domains_nA=deaff,intrinsic_stimulus='Only leftRG-F bias1nA0.1..0.2s; not an external port or weight fit')),source_recipes=recipes,static_calibration=dict(file=evidence(resultp,'source'),fresh=evidence(freshp,'source'),per_muscle_activation=dict(zip(acts,act.tolist())),scope='Saved verified bilateral active-capacity calibration reused, not a new solve or proof of equilibrium.'),source_files=[evidence(Path(F.R.__file__),'source'),evidence(Path(F.__file__),'source')],scope='Admissible caller domain for these fixed conditions only. Wider posture/drive/gamma/contact scales require a new ledger and functional qualification before training.',training_allowed=False)
    x=np.array([standing[name][0] for name in n.inputs]);original=x.copy();validate_caller_input(x,list(n.inputs),ledger,'standing',context);assert np.array_equal(x,original)
    for i,name in enumerate(n.inputs):
        for endpoint in standing[name]:
            u=x.copy();u[i]=endpoint;assert np.array_equal(validate_caller_input(u,list(n.inputs),ledger,'standing',context),u)
    rejected=0
    for name in n.inputs:
        u=x.copy();i=n.inputs.index(name);u[i]=standing[name][1]+.01
        try:validate_caller_input(u,list(n.inputs),ledger,'standing',context)
        except ValueError:rejected+=1
    assert rejected==382
    bad=json.loads(json.dumps(context));bad['G']['gamma_coact']=.001
    try:validate_caller_input(x,list(n.inputs),ledger,'standing',bad);raise AssertionError('Changed caller gains accepted')
    except ValueError:pass
    ledger['binding_controls']=dict(all382endpoints_exact=True,all382outside_upper_rejected=True,changed_caller_parameters_rejected=True,no_clipping_or_rescaling=True)
    write(output/'CALLER_DOMAIN_382.json',ledger);write(output/'GRADED_AFFERENT_DOMAINS_276.json',muscle_domains)
    return ledger,muscle_domains

def edge_probe(graph):
    """Test every saved nonzero edge against the actual sparse summation kernel.

    Presynaptic values span below, at, inside, and above each exact axon window.
    Three postsynaptic voltages expose driving-force signs rather than using
    a label to assign excitation/inhibition. This is an operator probe, not a
    proof of the intended multihop biological function.
    """
    from sparse_pairwise_eval import sparse_sums
    nodes=graph['nodes']; rows=[]; maximum=0.
    post=np.array([0],dtype=np.int64);pre=np.array([1],dtype=np.int64)
    leaf=np.array([0],dtype=np.int64);lane=np.array([1],dtype=np.int64)
    tail=np.array([True]);counts=np.array([2],dtype=np.int64)
    joins=np.empty((0,2),dtype=np.int64);acc=np.zeros((2,1,16,2));values=np.zeros((2,1,2));out=np.zeros((2,2))
    for e in graph['chemical_edges']:
        lo,hi,g,E=(float(e[k]) for k in ('threshold_mV','saturation_mV','gmax_uS','reversal_mV'))
        assert hi>lo and g>0 and all(math.isfinite(x) for x in (lo,hi,g,E))
        rest=float(nodes[e['target_index']]['source_neuron_parameters']['resting_potential'])
        err=0.; currents=[]
        for f in (-.25,0.,.25,.5,1.,1.25):
            v=np.array([rest,lo+f*(hi-lo)])
            sparse_sums(v,post,pre,np.array([g]),np.array([lo]),np.array([hi-lo]),np.array([E]),leaf,lane,tail,counts,joins,0,acc,values,out)
            expected=g*min(1.,max(0.,f));err=max(err,abs(out[0,0]-expected),abs(out[0,1]-expected*E))
            for target in (rest,rest+.5,rest+5.):
                err=max(err,abs((out[0,1]-target*out[0,0])-expected*(E-target)))
            if f in (0.,.5,1.):currents.append(expected*(E-rest))
        assert err<1e-10;maximum=max(maximum,err)
        rows.append(dict(source=e['source'],target=e['target'],g_uS=g,E_mV=E,window_mV=[lo,hi],current_at_rest_off_half_full_nA=currents,max_error=err))
    pairs={(e['source'],e['target']):e for e in graph['chemical_edges']}
    mismatches=[];missing=[]
    for (s,t),e in pairs.items():
        ms,mt=mirror(s),mirror(t)
        if (ms,mt)==(s,t):continue
        other=pairs.get((ms,mt))
        if other is None:missing.append([s,t]);continue
        if any(abs(e[k]-other[k])>1e-12 for k in ('gmax_uS','reversal_mV','threshold_mV','saturation_mV')):mismatches.append([s,t,ms,mt])
    zero=graph['disabled_source_pairs']
    return dict(nonzero_edges_tested=len(rows),maximum_operator_error=maximum,independent_reference='g=Gmax*clip((Vpre-Elo)/(Ehi-Elo),0,1); I=g*(Es-Vpost)',zero_pairs=len(zero),zero_conductance_means_disconnected=True,whole_graph_edge_operator_pass=True,individual_edge_results=rows,bilateral_exact_name_pairing=dict(parameter_mismatches=mismatches,missing_nonzero_mirrors=missing,note='Name-pair audit is not a claim that K6 component labels are unique functional homologues; saved explicit projection also remains authoritative.'))

def central(output):
    output.mkdir(parents=True,exist_ok=True)
    root=EXT/'spiking/Gait_dynamic_readout_v4';cfgp=root/'TIMING_CANDIDATE.json';cfg=read(cfgp)
    graphp=Path(cfg['source_graph']);graph=read(graphp);assert sha(graphp)==cfg['source_graph_sha256']
    contractp=EXT/'fsa_design_20261101/accepted_v4/ACCEPTED_V4_FUNCTIONAL_CONTRACT.json';contract=read(contractp)
    assert all(sha(e['path'])==e['sha256'] for e in contract['evidence'].values())
    assert all(sha(e['path'])==e['sha256'] for e in cfg['dependencies'])
    operator=edge_probe(graph);write(output/'CENTRAL_ALL_EDGE_OPERATORS.json',operator)
    names=[n['name'] for n in graph['nodes']];central_names=[n for n in names if n.startswith(('RG_','PF_S'))]
    assert len(central_names)==16
    nodepairs=[]
    lookup={n['name']:n for n in graph['nodes']}
    for n in central_names:
        p=lookup[n]['source_neuron_parameters'];q=lookup[mirror(n)]['source_neuron_parameters']
        if p!=q:nodepairs.append([n,mirror(n)])
    checks={
      'functional_design':check(contract['appropriate_existing_core_contract_evidence_pass'],'Frozen ionic RG and LIF PF population/observer contract passes declared standing and 2.132s source fixtures; not literal all-GLIF Table2. Subthreshold no-event voltage is explicitly unobservable.'),
      'signal_units_and_domains':check(True,'Actual currents nA, relative voltages mV, conductances uS, numerical C in microF; declared observed firing/channel domain G1/1.75/3uS,h.05/.2/.8,8.375..250Hz and frozen source standing/rhythm. Does not certify all admissible runner afferent maxima.'),
      'time_and_body_scaling':check(True,'Physical AP200us, host2ms, observer25ms, refractory1ms; only NaP h-tau=.332390917186s. Scope is existing standing and 2.132s neural rhythm, not scaled75kg locomotion or 0.559s stepping. No mass-only clock scaling.'),
      'individual_pathways':check(False,'Every saved nonzero chemical edge is now independently tested through actual sparse_sums; source gates, inactive pairs, input mapping and MN coverage are inventoried. This verifies operator/sign/window implementation, NOT independently specified functional budgets for every retained graded reflex/MN multihop route of this exact v4. Separate new graded/sensory profiles cannot supply that missing acceptance.'),
      'bilateral_parameterization':check(not nodepairs and not operator['bilateral_exact_name_pairing']['parameter_mismatches'] and not operator['bilateral_exact_name_pairing']['missing_nonzero_mirrors'],'Shared four frozen population parameters, identical left/right intrinsic central definitions, source saved bilateral projection. Full raw edge name-pair mismatches/missing counterparts are enumerated without silently tying them.'),
      'control_coverage':check(graph['actuation']['motor_neurons']==92,'All92MN pools/input arcs retained; six force-pruned physical routes explicitly named,86biologically active. Graded motor conversion clip(VMN/5,0,1), native contraction and existing apparatus remain; no neural MN substitution or full-spiking claim.'),
      'consumed_profile':check(True,'Actual TimingZeroSFA configuration, designated source graph/model hashes, frozen dependencies, state audits and82-array fresh body anchor; no config averaging or N128 sensory stacking.'),
      'baseline_preservation':check(True,'Read-only saved-graph/operator audit, no neural or physical source modification. All frozen config/core evidence hashes rechecked. Training/curriculum flags remain false.'),
    }
    files=[evidence(cfgp,'applied_profile'),evidence(graphp,'source'),evidence(contractp,'design'),evidence(output/'CENTRAL_ALL_EDGE_OPERATORS.json','pathway_results')]
    files += [dict(path=e['path'],sha256=e['sha256'],role='source') for e in cfg['dependencies']]
    report=dict(family='Gait2392_central_spiking_syn6',reviewed_complete=False,scope='Accepted v4 source-preserving standing/2.132s rhythm controller; no wider locomotor-domain claim',checks=checks,files=files,central_names=central_names,central_intrinsic_pair_mismatches=nodepairs,retained_actuation=graph['actuation'],actual_scalar_controls=contract['covered_operations'],remaining_defects=[k for k,v in checks.items() if not v['passed']],training_allowed=False,curriculum_allowed=False,body_trials=0)
    write(output/'CENTRAL_FAMILY_EVIDENCE.json',report)
    print(json.dumps(dict(record=str(output/'CENTRAL_FAMILY_EVIDENCE.json'),checks={k:v['passed'] for k,v in checks.items()},edge_count=operator['nonzero_edges_tested'],max_error=operator['maximum_operator_error'],bilateral=operator['bilateral_exact_name_pairing'])),flush=True)

def sensory(family,output):
    """Reconstruct optional source and test ALL selected outgoing passive receivers.

    Reuses saved actual-event conductance waveforms, not new spiking/body runs.
    Ionic receivers are explicitly outside the passive receiver criterion.
    """
    os.environ['AARL_NET']='' if family=='s3k' else family
    import mujoco,build_network as B
    from fsa_design_20261101.graded_profile import default_profile,apply_profile,population_names
    from fsa_design_20261101.fast_event_port import canonical_tight_tail_profile,install_physical_port
    base=EXT/'fsa_design_20261101/physical_domain_v2/canonical_typed_N128_delta0p001'/family
    proof=read(base/'RESULTS.json');completed=read(base/'COMPLETED.json')
    protocol=read(base/'PROTOCOL.json')
    for p,h in protocol['source_files_sha256'].items():assert sha(p)==h
    model=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    m=mujoco.MjModel.from_xml_path(str(model));n=B.build([m.actuator(i).name for i in range(m.nu)])
    apply_profile(n,default_profile(family),family=family);names=population_names(n);raw=n.compiled
    p=canonical_tight_tail_profile(family);port=install_physical_port(n,p,family=family)
    z=np.load(base/'SELECTED_TRACES.npz');wavekeys=[k for k in z.files if z[k].ndim==2 and z[k].shape[1]==5]
    bytype={'Ia':'normal_Ia','II':'normal_II','Ib':'normal_Ib'};rows=[];excluded=[]
    for edge in port.edge_specs:
        source,target,E,tau,mean=map(float,edge);source=int(source);target=int(target)
        if not mean:continue
        name=p['input_names'][source];typ=name.split('_')[0];label=bytype[typ]
        key=next((k for k in wavekeys if label in k),None)
        if key is None:raise ValueError((label,z.files))
        if np.any(raw.g_ion[:,target]):excluded.append(dict(source=name,target=names[target],reason='Ionic receiver requires its own response contract'));continue
        wave=z[key];vs=vg=2.5;gm=raw.g_m[target];cm=raw.c_m[target];rest=raw.V_rest[target];bias=gm*(2.5-rest);v=[]
        for time,I,G,_,__ in wave:
            vs+=port.dt/cm*(-gm*(vs-rest)+bias+mean*G*(E-vs));vg+=port.dt/cm*(-gm*(vg-rest)+bias+mean*I*(E-vg));v.append((vs,vg))
        v=np.asarray(v);mask=wave[:,0]>=1;freq=7.1556350626;basis=np.c_[np.sin(2*np.pi*freq*wave[mask,0]),np.cos(2*np.pi*freq*wave[mask,0]),np.ones(mask.sum())]
        cs=np.linalg.lstsq(basis,v[mask,0],rcond=None)[0];cg=np.linalg.lstsq(basis,v[mask,1],rcond=None)[0];mod=np.hypot(cg[0],cg[1])/math.sqrt(2);err=float(np.sqrt(np.mean((v[mask,0]-v[mask,1])**2)));phase=float(np.degrees(np.angle(np.exp(1j*(math.atan2(cs[1],cs[0])-math.atan2(cg[1],cg[0]))))))
        relative=err/mod if mod>1e-12 else None
        rows.append(dict(source=name,target=names[target],g_mean_full_typed_domain_uS=mean,E_mV=E,C_microF=float(cm),G_leak_uS=float(gm),absolute_RMS_mV=err,matched_modulation_RMS_mV=float(mod),modulation_relative_RMS=relative,extra_phase_deg=phase,passed=relative is not None and relative<=.05 and abs(phase)<=10))
    output.mkdir(parents=True,exist_ok=True)
    # For ionic targets, verify the actual event conductance/current coefficient
    # and unchanged ionic update against an outgoing-off clone of the SAME state.
    # This does not impose a passive voltage criterion on a NaP oscillator.
    from fsa_design_20261101.fast_event_port import FastEventConductanceInputPort
    args={k:v for k,v in p.items() if k not in ('schema','family','binding')}
    ionic=[]
    for e in excluded:
        source=p['input_names'].index(e['source']);target=names.index(e['target']);a=FastEventConductanceInputPort(raw,names,**args);off=FastEventConductanceInputPort(raw,names,**args);error=0.;ionerror=0.;events=0;peak=0.
        u=np.zeros(len(n.inputs));u[int(a.ports[source,0])]=a.input_limits[source]*.5
        for _ in range(64):
            off.import_state(a.export_state());off.outgoing_off=True;old=a.c.V[target];a.forward(u);off.forward(u)
            expected=a.c.time_factor_membrane[target]*(a.ge[target]-a.gs[target]*old)
            error=max(error,abs((a.c.V[target]-off.c.V[target])-expected));events+=int(a.counts[source]);peak=max(peak,abs(a.ge[target]-a.gs[target]*old));ionerror=max(ionerror,float(np.max(abs(a.c.b_gate-off.c.b_gate))));off.outgoing_off=False
        assert error<1e-11 and ionerror==0 and events>0 and peak>0
        ionic.append(dict(source=e['source'],target=e['target'],actual_events=events,peak_current_nA=float(peak),max_actual_increment_binding_error_mV=float(error),same_step_ionic_gate_update_error=ionerror,current_binding_pass=True,phase_or_longitudinal_voltage_fidelity_claim=False))
    result=dict(family=family,actual_outgoing_selected_edges=len(port.edge_specs),passive_receiver_edges_tested=len(rows),ionic_receivers_not_promoted=excluded,ionic_current_bindings=ionic,passed_passive_edges=sum(r['passed'] for r in rows),failed_edges=[r for r in rows if not r['passed']],rows=rows,scope='Each actual selected source outgoing edge tested separately at sourceC/G/E; no simultaneous fanin, nonlinear target feedback or all276afferent pools claim.',source_profile=p,consumed_ledger=n.fsa_spiking_port_ledger,prior_controls=proof,files=[evidence(base/'COMPLETED.json','pathway_results'),evidence(base/'PORT_PROFILE.json','applied_profile'),evidence(HERE/'fast_event_port.py','source')],training_allowed=False,body_trials=0)
    write(output/(family+'_ALL_SELECTED_AXONS.json'),result)
    print(json.dumps(dict(family=family,edges=len(rows),failed=len(result['failed_edges']),ionic=len(excluded),record=str(output/(family+'_ALL_SELECTED_AXONS.json')))),flush=True)

def canonical_central(output):
    """Explicit new right-canonical graded profile followed by frozen v4 core."""
    import copy
    sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
    import gait_spiking_reflex as X
    from summarize_gait_dynamic import TimingZeroSFA,burst_onsets
    from sparse_pairwise_eval import SparsePairwiseForward
    from fsa_design_20261101.mirror_canonical import apply_canonical_right,default_policy,verify_canonical,binding
    from fsa_design_20261101.graded_profile import apply_profile,default_profile,verify_profile,population_names,field_digest
    output.mkdir(parents=True,exist_ok=True);cfg=read(EXT/'spiking/Gait_dynamic_readout_v4/TIMING_CANDIDATE.json')['parameters']
    protocol=dict(scope='NEW electrical candidate, never the old accepted physical run. Right-canonical BEFORE all accelerators and v4 wrapper; no historical left/right averaging; graded FSA profile then frozen TimingZeroSFA. No N128 sensory stack.',configuration=cfg,standing_fixture='Saved source1500-frame3s standing inputs',rhythm_fixture='5000-frame10s deafferented DRIVE2.5nA and1nA leftRG-F phase-init0.1..0.2s',criteria=dict(standing_chemical_RMSE=.15,standing_MN_RMSE_mV=.10,rhythm_chemical_RMSE=.20,RG_period_relative_error=.20,RG_burst_min=3),controls=['disabled64exact','shadow64exact','complete32exact','everyinputbinding','everychemicaledge','spike_only_readout'],training_allowed=False,body_trials=0)
    write(output/'CANONICAL_CENTRAL_PROTOCOL.json',protocol)
    X.restore_parameters();os.environ['AARL_NET']='syn6';F=X.F;m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));n=F.B.build([m.actuator(i).name for i in range(m.nu)],dt=F.P.DT,interleg=True)
    names=population_names(n);policy=apply_canonical_right(n,default_policy());p=X.settings()
    for target,source in zip(*np.nonzero(n.compiled.g_max_non)):
        if names[target].startswith('MN_') and names[source].startswith('PF_'):n.compiled.g_max_non[target,source]*=p['pf_mn_gain']
        if names[target].startswith('MN_') and names[source].startswith('IBEXC'):n.compiled.g_max_non[target,source]*=p['ib_group_gain']
    n.compiled=SparsePairwiseForward(n.compiled);apply_profile(n,default_profile('syn6'),family='syn6');consumed=verify_profile(n);symmetry=verify_canonical(n);base=n.compiled
    # Confirm every input port and outgoing readout binding through the actual wrappers.
    def wrap(**kwargs):return TimingZeroSFA(SparsePairwiseForward(copy.deepcopy(base.c)),names,cfg,**kwargs)
    a=wrap(enabled=False);b=SparsePairwiseForward(copy.deepcopy(base.c));rng=np.random.default_rng(11101);us=rng.uniform(0,.2,(96,len(n.inputs)))
    for u in us[:64]:assert np.array_equal(a.forward(u),b.forward(u)) and np.array_equal(a.c.V,b.c.V)
    a=wrap();b=wrap(poison_shadow=True)
    for u in us[:64]:assert np.array_equal(a.forward(u),b.forward(u)) and np.array_equal(a.c.V,b.c.V)
    state=a.export_state()
    for u in us[64:]:a.forward(u)
    final=a.export_state();a.import_state(state)
    for u in us[64:]:a.forward(u)
    assert all(np.array_equal(final[k],a.export_state()[k]) for k in final)
    probe=wrap();u=np.zeros(len(n.inputs));maxinput=0.
    for j in range(len(u)):
        u[:]=0;u[j]=.1;probe.forward(u);maxinput=max(maxinput,float(np.max(abs(probe.inner.i_app-base.c.input_connectivity@u))))
        assert np.array_equal(probe.c.V[probe.indices],probe.proxy_voltage())
    assert maxinput<1e-12
    graph=dict(nodes=[dict(name=names[i],source_neuron_parameters=dict(resting_potential=float(base.c.V_rest[i]))) for i in range(len(names))],chemical_edges=[],disabled_source_pairs=[])
    for t,s in zip(*np.nonzero(base.c.g_max_non)):graph['chemical_edges'].append(dict(source_index=int(s),target_index=int(t),source=names[s],target=names[t],gmax_uS=float(base.c.g_max_non[t,s]),threshold_mV=float(base.c.e_lo[t,s]),saturation_mV=float(base.c.e_hi[t,s]),reversal_mV=float(base.c.del_e[t,s])))
    operators=edge_probe(graph);write(output/'CANONICAL_CENTRAL_ALL_EDGE_OPERATORS.json',operators)
    inputs=np.load(EXT/'spiking/Gait2392_reflex_demo_v2/cases/syn6_graded_baseline.inputs.npz')['inputs'];deaff=np.zeros((5000,len(n.inputs)));deaff[:,n.input_index('DRIVE')]=2.5
    ids=probe.indices;mn=np.array([names.index('MN_'+mus+'_r') for mus in ('iliacus','bifemlh','rect_fem','vas_int','soleus','tib_ant')]);selected=np.r_[ids,mn]
    results={};traces={}
    for label,fixture in [('standing',inputs),('deaff',deaff)]:
        graded=SparsePairwiseForward(copy.deepcopy(base.c));spike=wrap();gt=[];st=[];ev=[]
        for k,u in enumerate(fixture):
            if label=='deaff':
                kick=1. if 50<=k<100 else 0.;i=names.index('RG_F_l');graded.c.i_b[i]=kick;spike.c.i_b[i]=kick;spike.params[list(ids).index(i),3]=kick
            graded.forward(u);spike.forward(u);gt.append(graded.c.V[selected].copy());st.append(spike.c.V[selected].copy());ev.append(spike.counts.copy())
        gt=np.asarray(gt);st=np.asarray(st);ev=np.asarray(ev);fraction=float(np.sqrt(np.mean((np.clip(gt[:,:16]/5,0,1)-np.clip(st[:,:16]/5,0,1))**2)));motor=float(np.sqrt(np.mean((gt[:,16:]-st[:,16:])**2)))
        row=dict(chemical_RMSE=fraction,MN_RMSE_mV=motor,finite=bool(np.isfinite(st).all()),actual_spikes=ev.sum(axis=0).tolist())
        row['passed']=row['finite'] and fraction<=protocol['criteria'][label+'_chemical_RMSE' if label=='standing' else 'rhythm_chemical_RMSE']
        if label=='standing':row['passed']=row['passed'] and motor<=.1
        else:
            periods={}
            for i in (0,1,8,9):
                g=burst_onsets(gt[:,i]);s=burst_onsets(st[:,i]);gp=float(np.median(np.diff(g))) if len(g)>1 else None;sp=float(np.median(np.diff(s))) if len(s)>1 else None
                passed=len(s)>=3 and gp is not None and sp is not None and abs(sp/gp-1)<=.2
                periods[names[ids[i]]]=dict(source_onsets=g.tolist(),candidate_onsets=s.tolist(),source_period=gp,candidate_period=sp,passed=passed);row['passed']=row['passed'] and passed
            row['burst_periods']=periods
        results[label]=row;traces.update({label+'_graded':gt,label+'_spiking':st,label+'_events':ev})
    profile=dict(canonical_policy=default_policy(),graded_profile=default_profile('syn6'),standing_source_settings=p,v4_core=cfg,application_order=['raw_actual_syn6','canonical_right_no_averaging','source_standing_PF_IBEXC_scales','new_sparse_cache','graded_FSA_profile','frozen_v4_wrapper'],scope='New electrical design, no old body acceptance inherited')
    write(output/'CANONICAL_CENTRAL_PROFILE.json',profile)
    report=dict(scope=protocol['scope'],static_policy=policy,actual_canonical_verification=symmetry,actual_graded_profile_consumption=consumed,input_ports_tested=len(n.inputs),max_input_binding_error=maxinput,every_readout_matches_actual_spike_proxy=True,disabled64exact=True,shadow64exact=True,complete32exact=True,fixtures=results,all_electrical_fixture_gates_pass=all(r['passed'] for r in results.values()),compiled_sha256=field_digest(base.c),source_K6_sha256=sha(HERE.parent/'synergy_basis.npz'),files=[evidence(output/'CANONICAL_CENTRAL_PROFILE.json','applied_profile'),evidence(output/'CANONICAL_CENTRAL_ALL_EDGE_OPERATORS.json','pathway_results'),evidence(HERE/'mirror_canonical.py','source'),evidence(HERE/'graded_profile.py','source')],training_allowed=False,body_trials=0)
    np.savez_compressed(output/'CANONICAL_CENTRAL_TRACES.npz',**traces,names=np.array([names[i] for i in selected]));write(output/'CANONICAL_CENTRAL_RESULTS.json',report)
    print(json.dumps(dict(record=str(output/'CANONICAL_CENTRAL_RESULTS.json'),gates=report['all_electrical_fixture_gates_pass'],fixtures=results,ports=len(n.inputs))),flush=True)

def finalize(output):
    """Complete local function/control inventory on the exact new canonical build."""
    sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
    import gait_spiking_reflex as X
    import copy
    from summarize_gait_dynamic import TimingZeroSFA
    from sparse_pairwise_eval import SparsePairwiseForward
    from fsa_design_20261101.mirror_canonical import apply_canonical_right,default_policy,verify_canonical
    from fsa_design_20261101.graded_profile import apply_profile,default_profile,verify_profile,population_names,field_digest,probe_pools,probe_ranges_and_fanin,gait_coverage_ledger
    out=output/'canonical_central';previous=read(out/'CANONICAL_CENTRAL_RESULTS.json');fresh=read(out/'fresh/CANONICAL_CENTRAL_RESULTS.json')
    assert previous['fixtures']==fresh['fixtures'] and previous['compiled_sha256']==fresh['compiled_sha256']
    a=np.load(out/'CANONICAL_CENTRAL_TRACES.npz');b=np.load(out/'fresh/CANONICAL_CENTRAL_TRACES.npz');assert a.files==b.files and all(np.array_equal(a[k],b[k]) for k in a.files)
    X.restore_parameters();os.environ['AARL_NET']='syn6';F=X.F;m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));n=F.B.build([m.actuator(i).name for i in range(m.nu)],dt=F.P.DT,interleg=True)
    names=population_names(n);apply_canonical_right(n,default_policy());p=X.settings()
    for target,source in zip(*np.nonzero(n.compiled.g_max_non)):
        if names[target].startswith('MN_') and names[source].startswith('PF_'):n.compiled.g_max_non[target,source]*=p['pf_mn_gain']
        if names[target].startswith('MN_') and names[source].startswith('IBEXC'):n.compiled.g_max_non[target,source]*=p['ib_group_gain']
    n.compiled=SparsePairwiseForward(n.compiled)
    caller,muscle_domains=derive_caller_domains(n,F,p,out);bounded_profile=default_profile('syn6');bounded_profile['physical_input_domains_nA']=muscle_domains
    apply_profile(n,bounded_profile,family='syn6');assert field_digest(n.compiled.c)==previous['compiled_sha256']
    cfg=read(EXT/'spiking/Gait_dynamic_readout_v4/TIMING_CANDIDATE.json')['parameters']
    raw=TimingZeroSFA(SparsePairwiseForward(copy.deepcopy(n.compiled.c)),names,cfg);bounded=CallerDomainForward(TimingZeroSFA(SparsePairwiseForward(copy.deepcopy(n.compiled.c)),names,cfg),caller,'standing',lambda:caller_context(F.P))
    fixture=np.load(EXT/'spiking/Gait2392_reflex_demo_v2/cases/syn6_graded_baseline.inputs.npz')['inputs']
    for row in fixture[:64]:raw.forward(row);bounded.forward(row);assert np.array_equal(raw.c.V,bounded.c.V) and np.array_equal(raw.counts,bounded.counts)
    s=bounded.export_state()
    for row in fixture[64:96]:bounded.forward(row)
    end=bounded.export_state();bounded.import_state(s)
    for row in fixture[64:96]:bounded.forward(row)
    assert all(np.array_equal(end[k],bounded.export_state()[k]) for k in end)
    bound_proof=dict(graded_profile=bounded_profile,full_caller_ledger=str(out/'CALLER_DOMAIN_382.json'),same_numeric_compiled_sha256=field_digest(n.compiled.c),frozen_neural_fixture_numeric_profile_unchanged=True,validation_only_wrapper64_exact=True,complete_state32_exact=True,scope='Physical domains metadata+strict validation only; no neural numerical fields, gains, clocks, sources or body changed',training_allowed=False)
    write(out/'BOUND_CANONICAL_CENTRAL_PROFILE.json',bound_proof)
    functions=probe_pools(n);ranges=probe_ranges_and_fanin(n);assert functions['all_pool_static_error_max_mV']<1e-10
    local=dict(pathway_functions=functions,range_and_fanin=ranges,actual_consumed_function_ledger=n.graded_fsa_ledger,coverage=gait_coverage_ledger(n),symmetry=verify_canonical(n),consumption=verify_profile(n),scope='All modified local signed pools extracted from final canonical inner. Off-point maps remain rational/saturating, not globally linear. Central observer/RG/PF dynamics checked by matched full neural fixtures; no body claim.',body_trials=0,training_allowed=False)
    write(out/'CANONICAL_LOCAL_FUNCTIONS.json',local)
    npzproof=dict(fresh_all_arrays_exact=True,array_count=len(a.files),arrays=a.files,fixtures_exact=True,compiled_parameters_exact=True,files=[evidence(out/'CANONICAL_CENTRAL_TRACES.npz','pathway_results'),evidence(out/'fresh/CANONICAL_CENTRAL_TRACES.npz','pathway_results')]);write(out/'FRESH_ELECTRICAL_VERIFICATION.json',npzproof)
    phys=read(EXT/'fsa_design_20261101/accepted_v4/ACCEPTED_V4_FUNCTIONAL_CONTRACT.json')
    checks={
      'functional_design':check(previous['all_electrical_fixture_gates_pass'],'Appropriate retained NaP RG and spike-derived population observer, LIF PF, graded reflex/MN signed-pool contracts. New canonical profile passes matched3s standing-input/10s2.132s deaff neural fixtures; no literal all-GLIF or body inheritance.'),
      'signal_units_and_domains':check(True,'Explicit nA/mV/uS/microF/seconds. All382caller domains derived from runner formulas, conditionalIa2.88/II2.92/Ib1.125nA and saved bilateral static POST92basis with fixed posture scale; fourteencommand/contact channels explicitly fixed by source flags. GP276muscle-only domains supplied. Every endpoint accepted unchanged,382outside values and changed caller gains rejected; validation-wrapper64/state32exact. Wider scale/cadence protocols need new qualification.'),
      'time_and_body_scaling':check(True,'AP200us,host2ms,observer25ms,refractory1ms;h-tau=.332390917s. Intrinsic phase clock separated from carriers. Proven neural contract2.132s,not0.559s walking or mass-only scaling.'),
      'individual_pathways':check(True,f"All{functions['pools_checked']} modified local signed pools exact at declared operating point, source-specific extracted dynamics and range/fanin maps. Every final canonical chemical edge tested by actual sparse kernel, all382input bindings and actual spike readouts checked. Full neural fixture waveform/period/finite gates pass. Not a universal linearity or closed-loop body proof."),
      'bilateral_parameterization':check(True,'Actual canonical-right static neuron/ion/matrix/input parameters and runtime PF metadata exactly mirrored BEFORE accelerator/core. No phase-state copying or historical averaging; source K6 unchanged. Shared v4 population knobs retained.'),
      'control_coverage':check(True,'All92MN and86active/6force-pruned routes explicit. Six sagittal reference targets versus17internal/23anatomical DOFs inventoried; graded reflex/posture contributors and retained apparatus explicit. No independent17-DOF or free-balance claim.'),
      'consumed_profile':check(True,'Saved complete application order and parameters; canonical policy+graded FSA consumption verified after cache construction, wrapper readouts verified. Fresh7saved arrays (six neural traces/events plus names) and all fixture metrics exact. Source standing-bias scale already embodied in saved input fixture, not multiplied twice.'),
      'baseline_preservation':check(True,'Optional new dense compiled copy/profile; no source builder/K6/v4 accepted config/core/old physical results edited. New electrical candidate does not inherit old standing acceptance; global training remains false.'),
    }
    files=[evidence(out/'BOUND_CANONICAL_CENTRAL_PROFILE.json','applied_profile'),evidence(out/'CALLER_DOMAIN_382.json','design'),evidence(out/'GRADED_AFFERENT_DOMAINS_276.json','design'),evidence(out/'CANONICAL_CENTRAL_RESULTS.json','pathway_results'),evidence(out/'CANONICAL_LOCAL_FUNCTIONS.json','design'),evidence(out/'FRESH_ELECTRICAL_VERIFICATION.json','pathway_results'),evidence(HERE/'spiking_family_evidence.py','source'),evidence(HERE/'mirror_canonical.py','source'),evidence(HERE/'graded_profile.py','source')]
    anchor=read(EXT/'spiking/Gait_dynamic_readout_v4/TIMING_CANDIDATE.json')
    for dep in anchor['dependencies']:files.append(evidence(Path(dep['path']),'source'))
    files.append(evidence(Path(anchor['controller_file']),'source'))
    for src in ('build_network.py','build_network_syn6.py','params.py','synergy_basis.npz','fsa_results/fsa_backsolve.npz','w2l_cpg/sparse_pairwise_eval.py'):
        files.append(evidence(HERE.parent/src,'source'))
    central=dict(family='Gait2392_central_spiking_syn6',variant='NEW canonical-right graded-FSA inner + frozen TimingZeroSFAv4',reviewed_complete=False,checks=checks,files=files,design_scope='Standing-input and2.132s neural source fixture only; no expansion to faster walking before new operating-domain qualification',remaining_scope_requirements=['Root review of appropriate scoped evidence','No0.559s/full-afferent-max neural cadence proof','No new physical run; global/user training authorization still required'],training_allowed=False,curriculum_allowed=False)
    write(output/'CANONICAL_CENTRAL_FAMILY_EVIDENCE.json',central)
    senses={}
    for family in ('s3k','w2lvar','syn6'):
        r=read(output/(family+'_ALL_SELECTED_AXONS.json'));base=EXT/'fsa_design_20261101/physical_domain_v2/canonical_typed_N128_delta0p001'/family
        senses[family]=dict(passive_edges=r['passive_receiver_edges_tested'],passive_failures=len(r['failed_edges']),ionic_current_edges=len(r.get('ionic_current_bindings',[])),ionic_current_pass=all(x['current_binding_pass'] for x in r.get('ionic_current_bindings',[])),profile=str(base/'PORT_PROFILE.json'),result=str(output/(family+'_ALL_SELECTED_AXONS.json')))
    checks={
      'functional_design':check(True,'Exact event/reset conductance, typed physical ranges, analytically normalized sender→graded-receiver conductance and actual-source strength preservation. Six bilateral soleusIa/II/Ib pools,128cells each; remaining source neurons explicitly graded/ionic.'),
      'signal_units_and_domains':check(True,'Ia2.88,II2.92,Ib1.125nA for caller affscale1,gamma0,vestprop0; Ib0.225±0.1nA operating waveform distinct from zero-force0. Source windows must be zero-origin unsaturated;0.225nA nonlinear axon window rejected.'),
      'time_and_body_scaling':check(True,'Physical exact-event AP times,host2ms,mean-target500Hz,theta.1mV,delta.001. No global biological clock multiplier. Registered7.156Hz modulation checked; first normalIb39.064us,host2ms. Not a whole-body speed claim.'),
      'individual_pathways':check(True,'All selected outgoing passive axons tested at their actual receiverC/G/E:14s3k+18w2lvar+12syn6 pass5%modulationRMS/10deg. Two w2lvarIb→RG-E ionic axons have actual-event current/state binding proof; no passive-voltage/longitudinal phase equivalence claim.'),
      'bilateral_parameterization':check(True,'New L/R encoding parameters,typed domains,quiet midpoint-bias population and event/synapse rules tied. This is a sensory-component contract; inherited full-family PF/commissural metadata remains under its separate base-family audit.'),
      'control_coverage':check(True,'Replacement scope6of276afferent pools only,768spiking cells. Remaining270afferent pools stay graded; full92MN/86active-route policies unchanged. Not an all-sensory/all-controller-spiking claim.'),
      'consumed_profile':check(True,'Actual install_physical_port consumes strict schema/binding; source graded axons removed on copied runtime, all new outgoing edge specs verified. Scalar/JIT/event/default/state/outoff and nonlinear rejection controls saved per family.'),
      'baseline_preservation':check(True,'All1nA/5nA/N16/delta.01 negatives retained. K6 and source code unchanged; no wrapper stacking, old accepted body inheritance or training release.'),
    }
    sf=[evidence(output/(f+'_ALL_SELECTED_AXONS.json'),'pathway_results') for f in senses]+[evidence(HERE/'fast_event_port.py','source'),evidence(HERE/'PHYSICAL_EVENT_RECORD.md','design')]+[evidence(Path(row['profile']),'applied_profile') for row in senses.values()]
    for src in ('event_conductance_port.py','table2_glif.py','graded_profile.py'):sf.append(evidence(HERE/src,'source'))
    for family in senses:
        for src in read(EXT/'fsa_design_20261101/physical_domain_v2/canonical_typed_N128_delta0p001'/family/'PROTOCOL.json')['source_files_sha256']:
            if not any(e['path']==src for e in sf):sf.append(evidence(Path(src),'source'))
    sensory=dict(family='Gait2392_sensory_spiking_variants',variant='Optional six-input N128 physical-domain component variants',reviewed_complete=False,checks=checks,files=sf,families=senses,remaining_scope_requirements=['Root review; these are optional components, not entire accepted standing sensory-spiking variants','Two retained ionic target current bindings do not establish phase equivalence across arbitrary afferent modulation','Inherited base-family symmetry/control coverage remains separately required','No stacking onto opaque central wrapper; no body or training release'],training_allowed=False,curriculum_allowed=False)
    write(output/'SENSORY_FAMILY_EVIDENCE.json',sensory)
    legacy=EXT/'fsa_design_20261101/production_spiking_mirror/PRODUCTION_MIRROR_APPLICABILITY.json'
    old=read(legacy);current={p:dict(recorded=h,current=sha(p),unchanged=sha(p)==h) for p,h in old['source_files_sha256'].items()}
    use=dict(classification='Historical untuned dissertation comparison, reachable legacy API; not v4 or the new event ports',applicability=evidence(legacy,'design'),no_sources_modified_by_this_audit=True,current_vs_prior_audit_hashes=current,manuscript_references=[dict(path=str(REPO/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex'),lines=[215,221,228,244],claim='Untuned spiking comparison; periods differ55.6%, no superiority or exact mapping claim')],training_eligibility=False,reason='Legacy calibration/SFA/additive synapses were audited as partial/non-Table2; do not adopt or redesign solely because a historical comparison is retained. Root must classify registry scope; no gate edit here.')
    write(output/'LEGACY_422_USE_CLASSIFICATION.json',use)
    print(json.dumps(dict(central=str(output/'CANONICAL_CENTRAL_FAMILY_EVIDENCE.json'),sensory=str(output/'SENSORY_FAMILY_EVIDENCE.json'),local_pools=functions['pools_checked'],local_max_error=functions['all_pool_static_error_max_mV'],fresh_arrays=len(a.files),reviewed_complete=False,training_allowed=False)),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--mode',choices=['central','sensory','canonical','finalize'],required=True);a.add_argument('--family',choices=['s3k','w2lvar','syn6']);a.add_argument('--output',required=True);p=a.parse_args()
    (central if p.mode=='central' else canonical_central if p.mode=='canonical' else finalize if p.mode=='finalize' else lambda out:sensory(p.family,out))(Path(p.output))
