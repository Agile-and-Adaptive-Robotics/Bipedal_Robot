"""Preregistered, bounded sensory-spiking reflex tests in actual Gait standing.

Only soleus-r Ia/Ib input populations change encoding. RG, PF, interneurons,
MNs, native muscles, geometry, and support remain the accepted syn6 condition.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,copy,hashlib,time,argparse,inspect,contextlib
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
FOLDER=REPO/'Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination'
sys.path.insert(0,str(FOLDER))
import gait_standing_family as F
from gait_standing_objective import score
from sparse_pairwise_eval import SparsePairwiseForward
import spiking_sns_inputs as S

OUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2')
PRISTINE_G=copy.deepcopy(F.P.G)
DIAGRAM=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families/syn6_graph_ledger.json'
BEST=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\rustuna_syn6_v1\BEST.json')
CHANNELS=('Ia_soleus_r','Ib_soleus_r');AMP=1.;ONSET=1.6;WIDTH=.1;DURATION=3.
CASES={
 'graded_baseline':('graded',None,None), 'spiking_baseline':('spiking',None,None),
 'graded_Ia_pulse':('graded','Ia',None),'graded_Ib_pulse':('graded','Ib',None),
 'spiking_Ia_pulse':('spiking','Ia',None),'spiking_Ib_pulse':('spiking','Ib',None),
 'spiking_Ia_off_baseline':('spiking',None,'Ia'),'spiking_Ia_off_pulse':('spiking','Ia','Ia'),
 'spiking_Ib_off_baseline':('spiking',None,'Ib'),'spiking_Ib_off_pulse':('spiking','Ib','Ib')}
COLS=['neural_time_s','injected_current_nA','actual_Ia_input_nA','actual_Ib_input_nA',
 'Ia_voltage_mV','Ib_voltage_mV','soleus_MN_mV','tib_ant_MN_mV','Ia_spikes_this_2ms','Ib_spikes_this_2ms',
 'Ia_decoded_fraction','Ib_decoded_fraction','Ia_direct_MN_current_nA','IBIN_MN_current_nA','IBEXC_MN_current_nA','total_soleus_MN_synaptic_current_nA']
OWNER=None;ORIGINAL_KERNEL=S.encode_step

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')

def trace_kernel(*args):
    result=ORIGINAL_KERNEL(*args)
    if OWNER is not None and OWNER.encoder is not None and args[1] is OWNER.encoder.v:
        t=OWNER.steps*OWNER.dt+OWNER.substep*OWNER.encoder.dt
        if ONSET-.2<=t<ONSET+.4:
            OWNER.event_t.append(t);OWNER.event_v.append(args[1].copy());OWNER.event_ref.append(args[2].copy());OWNER.event_spikes.append(args[4].copy())
        OWNER.substep+=1
    return result
S.encode_step=trace_kernel

class TwoSourceEncoder(S.SpikeInputForward):
    def __init__(self,inner,names):
        super().__init__(inner,names,population_size=16)
        wanted=np.array([list(self.indices).index(names.index(name)) for name in CHANNELS],dtype=int)
        for key in ('indices','ports','v','ref','decoded','spikes','counts','total_counts','target','tau','input_gain','span','rest'):
            setattr(self,key,getattr(self,key)[wanted].copy())
        assert [names[i] for i in self.indices]==list(CHANNELS)

def raw_snapshot(base):
    raw=base.c;state={}
    for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last','outputs'):
        if hasattr(raw,key):state['raw_'+key]=np.asarray(getattr(raw,key)).copy()
    for key in ('last_voltage','input_values','i_app','output_values','sums','acc','values'):
        if hasattr(base,key):state['sparse_'+key]=np.asarray(getattr(base,key)).copy()
    return state

def raw_restore(base,state):
    for key,value in state.items():
        if key.startswith('raw_'):setattr(base.c,key[4:],value.copy())
        elif key.startswith('sparse_'):np.copyto(getattr(base,key[7:]),value)
    base.c.V_last=base.last_voltage
    if base.output_onehot:base.c.outputs=base.output_values

class ReflexForward:
    def __init__(self,inner,names,mode='graded',pulse=None,off=None,amplitude=AMP):
        self.names=list(names);self.mode=mode;self.pulse=pulse;self.off=off;self.amplitude=amplitude
        self.steps=0;self.substep=0;self.dt=float(inner.c.dt);self.rows=[];self.event_t=[];self.event_v=[];self.event_ref=[];self.event_spikes=[]
        self.source_indices=np.array([names.index(n) for n in CHANNELS]);self.mn=names.index('MN_soleus_r');self.ant=names.index('MN_tib_ant_r')
        self.ports=np.array([np.flatnonzero(inner.c.input_connectivity[i])[0] for i in self.source_indices])
        self.disabled=[]
        if off:
            col=names.index(off+'_soleus_r');raw=inner.c
            for post in np.flatnonzero(raw.g_max_non[:,col]):
                self.disabled.append(dict(source=names[col],target=names[post],before_gmax_uS=float(raw.g_max_non[post,col]),after_gmax_uS=0.))
            raw.g_max_non[:,col]=0.;inner=SparsePairwiseForward(raw)
            assert not np.any(inner.pre==col)
        self.encoder=TwoSourceEncoder(inner,names) if mode=='spiking' else None
        self.base=self.encoder if self.encoder is not None else inner;self.sparse=inner
        self.config=json.dumps(dict(mode=mode,pulse=pulse,off=off,amplitude=amplitude,onset=ONSET,width=WIDTH,channels=CHANNELS),sort_keys=True)
    def __getattr__(self,name):return getattr(self.base,name)
    def forward(self,x=None):
        global OWNER
        OWNER=self;u=np.array(x,copy=True);t=self.steps*self.dt;self.substep=0
        pulse_start=round(ONSET/self.dt);pulse_steps=round(WIDTH/self.dt)
        injected=self.amplitude if self.pulse and pulse_start<=self.steps<pulse_start+pulse_steps else 0.
        if self.pulse:u[self.ports[0 if self.pulse=='Ia' else 1]]+=injected
        answer=self.base.forward(u);raw=self.sparse.c;sp=self.sparse
        pre=sp.pre;post=sp.post
        frac=np.clip((raw.V_last[pre]-sp.elo)/sp.erange,0.,1.)
        currents=sp.gmax*frac*(sp.reversal-raw.V_last[post]);to_mn=post==self.mn
        direct=float(currents[to_mn&(pre==self.source_indices[0])].sum())
        ibin=float(currents[to_mn&np.array([self.names[i].startswith('IBIN_') for i in pre])].sum())
        ibexc=float(currents[to_mn&np.array([self.names[i].startswith('IBEXC_') for i in pre])].sum())
        counts=self.encoder.counts.copy() if self.encoder else np.zeros(2)
        decoded=self.encoder.decoded.copy() if self.encoder else np.full(2,np.nan)
        self.rows.append([t,injected,*u[self.ports],*raw.V[self.source_indices],raw.V[self.mn],raw.V[self.ant],*counts,*decoded,direct,ibin,ibexc,float(currents[to_mn].sum())])
        self.steps+=1;return answer
    def export_state(self):
        state=self.encoder.export_state() if self.encoder else raw_snapshot(self.sparse)
        state.update(reflex_version=np.array(1),reflex_config=np.array(self.config),reflex_steps=np.array(self.steps),reflex_substep=np.array(self.substep),
           reflex_rows=np.asarray(self.rows).reshape(-1,len(COLS)),reflex_event_t=np.asarray(self.event_t),
           reflex_event_v=np.asarray(self.event_v).reshape(-1,2,16),reflex_event_ref=np.asarray(self.event_ref).reshape(-1,2,16),reflex_event_spikes=np.asarray(self.event_spikes,dtype=np.int8).reshape(-1,2,16))
        return state
    def import_state(self,state):
        assert int(state['reflex_version'])==1 and str(state['reflex_config'])==self.config
        if self.encoder:self.encoder.import_state(state)
        else:raw_restore(self.sparse,state)
        self.steps=int(state['reflex_steps']);self.substep=int(state['reflex_substep']);self.rows=state['reflex_rows'].tolist()
        for key in ('event_t','event_v','event_ref','event_spikes'):setattr(self,key,[x.copy() for x in state['reflex_'+key]])
    def save(self,prefix):
        np.savez_compressed(prefix.with_suffix('.reflex.npz'),rows=np.asarray(self.rows),columns=COLS,
          event_t=np.asarray(self.event_t),event_v=np.asarray(self.event_v).reshape(-1,2,16),event_ref=np.asarray(self.event_ref).reshape(-1,2,16),
          event_spikes=np.asarray(self.event_spikes,dtype=np.int8).reshape(-1,2,16),channels=CHANNELS,metadata_json=json.dumps(self.metadata()))
    def metadata(self):
        return dict(scope='Only two pure-input Ia/Ib source populations spiking; all RG/PF/IN/MN graded.',mode=self.mode,
          spiking=self.encoder.metadata() if self.encoder else None,pulse_current_nA=self.amplitude,pulse_neural_onset_s=ONSET,pulse_nominal_free_physics_onset_s=ONSET-F.R.WARMUP,pulse_width_s=WIDTH,
          disabled_outgoing_edges=self.disabled,voltage_scope='LIF membranes threshold-normalized, threshold1/reset0; decoded source and gradedMN voltages in source relative-rest mV, not absolute membrane mV.',
          physical_scope='Existing0.6s pinned preparation only; thereafter actual free native muscle/contact physics under unchanged inherited support apparatus.')

def settings():return json.loads(BEST.read_text())['trial']['params']
def core_digest():return hashlib.sha256((inspect.getsource(TwoSourceEncoder)+inspect.getsource(ReflexForward)+inspect.getsource(trace_kernel)).encode()).hexdigest()
def restore_parameters():
    ledger=json.loads(DIAGRAM.read_text(encoding='utf-8'));F.P.G.clear();F.P.G.update(copy.deepcopy(PRISTINE_G))
    F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.
    assert F.P.G==ledger['global_G'],'Saved diagram input metadata must match independently restored source defaults'
    assert F.P.G['heel_rge']==F.P.G['toe_rge']==0.
    return ledger
def construct():
    ledger=restore_parameters();os.environ['AARL_NET']='syn6';m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));acts=[m.actuator(i).name for i in range(m.nu)]
    n=F.B.build(acts,dt=F.P.DT,interleg=True);F.tie_weights(n);names=F.canonical_neuron_names(n);p=settings()
    for post,pre in zip(*np.nonzero(n.compiled.g_max_non)):
        if names[post].startswith('MN_') and names[pre].startswith('PF_'):n.compiled.g_max_non[post,pre]*=p['pf_mn_gain']
        if names[post].startswith('MN_') and names[pre].startswith('IBEXC'):n.compiled.g_max_non[post,pre]*=p['ib_group_gain']
    for key,expected in ledger['compiled_matrix_sha256'].items():assert hashlib.sha256(np.asarray(getattr(n.compiled,key)).tobytes()).hexdigest()==expected,key
    return SparsePairwiseForward(n.compiled),names

def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    protocol=dict(machine='easteregg2',source_diagram_ledger=str(DIAGRAM),source_diagram_sha256=digest(DIAGRAM),standing_configuration=str(BEST),configuration_sha256=digest(BEST),
       cases={k:list(v) for k,v in CASES.items()},duration_neural_s=DURATION,pinned_preparation_s=F.R.WARMUP,free_physics_s=DURATION-F.R.WARMUP,
       source_channels=CHANNELS,spiking_cells=32,remaining_graded_cells=792,encoder='Existing deterministic reset/refractory LIF encoder,16cells per population;200µs substeps,5ms membrane,1ms refractory,200Hz max; source-decay-matched spike filters.',
       pulse=dict(current_nA=AMP,width_s=WIDTH,neural_onset_s=ONSET,free_physics_onset_s=ONSET-F.R.WARMUP),
       controls='Matched graded and narrow sensory-spiking no-pulse baselines; Ia/Ib pulses; for each spiking source, a separate outgoing-conductance-zero pulse/no-pulse pair. No direct MN stimulus.',
       predictions='Ia expected direct homonymous excitation with reciprocal inhibition; Ib net sign must be measured because inhibitory IBIN and additive load-sharing IBEXC branches coexist. No multiplicative AND-gate claim.',
       physical_acceptance='Original F.run hyperextension/nonfinite stop and saved standing score; no clamp beyond inherited initialization, no added forces, changed gains, or new tuning.',
       preregistered_metrics=['pre-onset exact physical/neural parity','actual1nA×100ms pulse integral','LIF counts/reset/refractory and decoded source voltage','actual gradedMN response','native muscle force/activation and joint-angle response','measured GRFs/apparatus force','source-output-off difference-of-differences'],
       figures='Separate Ia/Ib columns: stimulus, actual raster/normalized membrane, decoded/graded afferent voltage, MN, force, ankle/knee, and actual contact/support. Pulse effects compared against matched no-pulse cases.',
       no_optimizer=True,script=str(Path(__file__)),script_sha256=digest(__file__),reflex_core_sha256=core_digest(),stock_encoder_sha256=digest(S.__file__))
    path=OUT/'PROTOCOL.json'
    if path.exists():
        old=json.loads(path.read_text());assert old['cases']==protocol['cases'] and old['pulse']==protocol['pulse']
        if not list((OUT/'cases').glob('*.result.json')):
            protocol['prephysical_source_history']=old.get('prephysical_source_history',[])+[old['script_sha256']]
            write(path,protocol)
        else:assert old['reflex_core_sha256']==protocol['reflex_core_sha256']
    else:write(path,protocol)
    record=OUT/'EXPERIMENT_RECORD.md'
    if not record.exists():record.write_text('# Gait2392 Ia/Ib sensory-spiking reflex demonstration\n\nMachine: easteregg2. Protocol registered before physical tests. Frozen accepted syn6 biological standing gains, source92 MNs/86 force-active routes, unchanged support, and0.6s initialization. Only Ia_soleus_r and Ib_soleus_r become16-cell LIF populations; central circuitry and MNs remain graded. No optimizer, new tuning, or whole-controller-spiking claim.\n\nProtocol: `'+str(path)+'`. Raw data remain in this external folder.\n',encoding='utf-8')
    inner,names=construct();u=np.random.default_rng(100101).uniform(0.,.2,(96,inner.c.input_connectivity.shape[1]));reference=SparsePairwiseForward(copy.deepcopy(inner.c));zero=ReflexForward(SparsePairwiseForward(copy.deepcopy(inner.c)),names,'graded',None,None,0.)
    maximum=0.
    for row in u:reference.forward(row);zero.forward(row);maximum=max(maximum,float(np.max(abs(reference.c.V-zero.c.V))))
    assert maximum==0.
    audits={}
    for mode in ('graded','spiking'):
        obj=ReflexForward(SparsePairwiseForward(copy.deepcopy(inner.c)),names,mode,'Ia',None)
        for row in u[:20]:obj.forward(row)
        state=obj.export_state()
        for row in u[20:40]:obj.forward(row)
        expected=obj.export_state();obj.import_state(state)
        for row in u[20:40]:obj.forward(row)
        actual=obj.export_state();equal=set(expected)==set(actual) and all(np.array_equal(expected[k],actual[k],equal_nan=True) if expected[k].dtype.kind in 'fc' else np.array_equal(expected[k],actual[k]) for k in expected)
        assert equal;audits[mode]=dict(complete_state_fields=len(actual),restored_continuation_exact=True)
    write(OUT/'OFFLINE_AUDIT.json',dict(source_matrix_hashes_match_accepted_diagram=True,zero_amplitude_default_max_error_mV=maximum,state_audits=audits,no_physics_run=True))
    print(json.dumps(dict(prepared=True,offline_audit=str(OUT/'OFFLINE_AUDIT.json'))),flush=True)

def run_case(case,fresh=False):
    if not (OUT/'OFFLINE_AUDIT.json').exists():raise ValueError('Offline/default/state gate missing')
    tag=case+('_fresh' if fresh else '');prefix=OUT/'cases'/('syn6_'+tag);comparison=prefix.with_suffix('.comparison.json')
    if comparison.exists():return json.loads(comparison.read_text())
    protocol=json.loads((OUT/'PROTOCOL.json').read_text());assert core_digest()==protocol['reflex_core_sha256']
    assert digest(BEST)==protocol['configuration_sha256'] and digest(DIAGRAM)==protocol['source_diagram_sha256']
    ledger=restore_parameters();p=settings();objects=[];mode,pulse,off=CASES[case]
    def transform(inner,names):
        for key,expected in ledger['compiled_matrix_sha256'].items():assert hashlib.sha256(np.asarray(getattr(inner.c,key)).tobytes()).hexdigest()==expected,key
        obj=ReflexForward(inner,names,mode,pulse,off);objects.append(obj);return obj
    start=time.perf_counter()
    with (OUT/'physics.log').open('a',encoding='utf-8') as log,contextlib.redirect_stdout(log):
        r=F.run('syn6',prefix.parent,duration=DURATION,tag=tag,gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],
          pinned_actuator_settling=True,compiled_transform=transform,capture_neural_inputs=True)
    obj=objects[0];obj.save(prefix);a=score(r,prefix)[1]
    report=dict(case=tag,configuration_case=case,validation_repeat_of=case if fresh else None,result=str(prefix.with_suffix('.result.json')),reflex=str(prefix.with_suffix('.reflex.npz')),assessment=a,metadata=obj.metadata(),wall_s=time.perf_counter()-start,source_protocol_sha256=digest(OUT/'PROTOCOL.json'),source_script_sha256=digest(__file__),reflex_core_sha256=core_digest(),failure=r['failure'])
    write(comparison,report)
    with (OUT/'CASES.jsonl').open('a',encoding='utf-8') as stream:stream.write(json.dumps(report)+'\n')
    print(json.dumps(dict(case=case,accepted=a['accepted'],failure=r['failure'],wall_s=report['wall_s'])),flush=True)
    return report

def load_case(case):
    prefix=OUT/'cases'/('syn6_'+case)
    with np.load(prefix.with_suffix('.reflex.npz')) as a:reflex={k:a[k] for k in a.files}
    with np.load(prefix.with_suffix('.motor_state.npz')) as a:motor={k:a[k] for k in a.files}
    with np.load(prefix.with_suffix('.npz')) as a:native={k:a[k] for k in ('t','q','qfull','ctrl','act','neuro','com','contact','footz','key_joints')}
    with np.load(prefix.with_suffix('.apparatus.npz')) as a:apparatus={k:a[k] for k in ('samples','neuro','columns')}
    report=json.loads(prefix.with_suffix('.comparison.json').read_text())
    result=json.loads(prefix.with_suffix('.result.json').read_text())
    return dict(prefix=prefix,reflex=reflex,motor=motor,native=native,apparatus=apparatus,report=report,result=result)

def prefix_audit():
    saved=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\pinned_preparation\syn6_prepared_bio_nominal_v1.npz')
    current=OUT/'cases/syn6_graded_baseline.npz';rows={}
    with np.load(saved) as old,np.load(current) as new:
        for key in ('t','qfull','q','ctrl','neuro','com'):
            a,b=old[key][:len(new[key])],new[key]
            rows[key]=dict(exact=bool(np.array_equal(a,b)),maximum_absolute_difference=float(np.max(abs(a-b))))
    report=dict(saved=str(saved),current=str(current),arrays=rows,exact=all(x['exact'] for x in rows.values()),scope='Actual default physical/native-neural prefix, including0.6s initialization; no learned or input-gain retuning.')
    write(OUT/'DEFAULT_PHYSICAL_PREFIX_AUDIT.json',report)
    print(json.dumps(dict(default_physical_prefix_exact=report['exact'],arrays=rows)),flush=True)
    if not report['exact']:raise ValueError('Default physical-prefix gate failed; do not run further reflex cases')
    return report

def analyze():
    cases={case:load_case(case) for case in CASES};summary=[];pairs=[]
    for mode in ('graded','spiking'):
        for channel in ('Ia','Ib'):
            pairs.append((mode+'_'+channel+'_pulse',mode+'_baseline',channel))
    for channel in ('Ia','Ib'):pairs.append(('spiking_'+channel+'_off_pulse','spiking_'+channel+'_off_baseline',channel))
    paired={};col={name:i for i,name in enumerate(COLS)}
    for pulse,baseline,channel in pairs:
        a,b=cases[pulse],cases[baseline];ra,rb=a['reflex']['rows'],b['reflex']['rows'];t=ra[:,0]
        pre=t<ONSET;prefix_exact=np.array_equal(ra[pre],rb[pre],equal_nan=True)
        physical_pre=a['motor']['time']<=ONSET-F.R.WARMUP+1e-10
        for key in ('activation','excitation','engine_force_cache_N'):prefix_exact=prefix_exact and np.array_equal(a['motor'][key][physical_pre],b['motor'][key][physical_pre])
        assert prefix_exact,pulse+' pre-onset mismatch'
        nativepre=a['native']['t']<ONSET
        for key in ('qfull','q','ctrl','neuro'):assert np.array_equal(a['native'][key][nativepre],b['native'][key][nativepre]),(pulse,key)
        dose=float(ra[:,col['injected_current_nA']].sum()*F.P.DT);assert abs(dose-AMP*WIDTH)<1e-12
        window=(t>=ONSET)&(t<ONSET+WIDTH+.25);delta=ra[:,col['soleus_MN_mV']]-rb[:,col['soleus_MN_mV']]
        midx=list(a['motor']['actuator_names']).index('soleus_r');fm=-a['motor']['engine_force_cache_N'][:,midx];fb=-b['motor']['engine_force_cache_N'][:,midx]
        assert min(fm)>-1e-8
        tf=a['motor']['time']+F.R.WARMUP-F.P.DT;fw=(tf>=ONSET)&(tf<ONSET+WIDTH+.25);df=fm-fb
        actadr=int(a['motor']['activation_addresses'][midx]);assert actadr>=0
        da=a['motor']['activation'][:,actadr]-b['motor']['activation'][:,actadr]
        de=a['motor']['excitation'][:,midx]-b['motor']['excitation'][:,midx]
        # runner.log_q is already degrees; qfull remains raw generalized units.
        joints=list(a['native']['key_joints']);ankle=a['native']['q'][:,joints.index('ankle_angle_r')]-b['native']['q'][:,joints.index('ankle_angle_r')];knee=a['native']['q'][:,joints.index('knee_angle_r')]-b['native']['q'][:,joints.index('knee_angle_r')]
        sample=dict(case=pulse,matched_baseline=baseline,pre_onset_bit_exact=bool(prefix_exact),actual_test_current_integral_nA_s=dose,
          delta_MN_mV_min=float(delta[window].min()),delta_MN_mV_max=float(delta[window].max()),delta_tension_N_min=float(df[fw].min()),delta_tension_N_max=float(df[fw].max()),
          delta_actual_activation_min=float(da[fw].min()),delta_actual_activation_max=float(da[fw].max()),delta_excitation_min=float(de[fw].min()),delta_excitation_max=float(de[fw].max()),
          maximum_abs_delta_ankle_deg=float(np.max(abs(ankle[window]))),maximum_abs_delta_knee_deg=float(np.max(abs(knee[window]))),
          causal_source_output_off=CASES[pulse][2] is not None)
        if sample['causal_source_output_off']:
            sample['downstream_physics_bit_exact']=all(np.array_equal(a['native'][key],b['native'][key]) for key in ('qfull','q','ctrl','neuro'))
            assert sample['downstream_physics_bit_exact'];assert np.array_equal(fm,fb) and np.max(abs(delta))==0.
        if CASES[pulse][0]=='spiking':
            k=0 if channel=='Ia' else 1;events=a['reflex'];spikes=events['event_spikes'];v=events['event_v'];ref=events['event_ref'];et=events['event_t']
            assert np.all(v[spikes.astype(bool)]==0.) and np.all(ref[spikes.astype(bool)]==.001)
            ew=(et>=ONSET)&(et<ONSET+WIDTH-1e-12);sample['actual_pulsed_population_spikes_during_100ms']=int(spikes[ew,k].sum())
            sample['actual_spike_reset_and_refractory_verified']=True
            decoded_error=np.max(abs(events['rows'][:,col[channel+'_voltage_mV']]-5.*events['rows'][:,col[channel+'_decoded_fraction']]))
            sample['decoded_voltage_identity_max_error_mV']=float(decoded_error);assert decoded_error==0.
        summary.append(sample);paired[pulse]=dict(neural_time_s=t,delta_MN_mV=delta,force_time_neural_s=tf,delta_tension_N=df,delta_ankle_deg=ankle,delta_knee_deg=knee)
    saved=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\pinned_preparation\syn6_prepared_bio_nominal_v1.npz')
    default={}
    with np.load(saved) as old:
        fresh=cases['graded_baseline']['native']
        for key in ('t','qfull','q','ctrl','neuro','com'):
            v=old[key][:len(fresh[key])];different=np.flatnonzero(np.any(v!=fresh[key],axis=tuple(range(1,v.ndim))) if v.ndim>1 else v!=fresh[key]);default[key]=dict(exact=bool(np.array_equal(v,fresh[key])),maximum_absolute_difference=float(np.max(abs(v-fresh[key]))),initial_maximum_difference=float(np.max(abs(v[0]-fresh[key][0]))),first_different_index=int(different[0]) if len(different) else None)
    conditions=[]
    for case,a in cases.items():
        r=a['result'];conditions.append(dict(case=case,native_failure=r['failure'],completed_free_samples=r['completed_samples'],free_physics_s=r['completed_samples']*.002,
          max_knee_angle_deg=r['max_knee_angle_deg'],ground_load_bw=r['mean_ground_load_bw'],upward_support_bw=r['mean_upward_rig_support_bw'],downward_support_bw=r['mean_downward_rig_force_bw'],
          maximum_tilt_deg=r['metrics']['bal_tilt_max'],sway_m=r['metrics']['bal_sway'],wall_s=a['report']['wall_s'],formal_long_standing_gate=a['report']['assessment']['accepted']))
    checkpoint=dict(machine='easteregg2',scope='Frozen syn6 graded-MN biological standing reflex tests: only two input populations are spiking.2.4s free physics each, not the existing long-standing/curriculum gate.',
      conditions=conditions,paired_effects=summary,graded_default_vs_saved_biological_prefix=default,
      formal_duration_note='The unchanged standing scorer requires≥3600 free samples; each short case has1200, so its accepted=False remains preserved. No gate is relaxed or walking/standing promotion asserted.',
      actual_source_topology='Syn6 IBEXC_ankle_pf_r sums nine source Ib inputs. Unlike base/CAD IBEXC, this selected compiled pool has no incoming RG-E edge and is not a multiplicative stance AND gate.',
      data_clock='Reflex rows/event timestamps include0.6s pinned preparation. Motor/apparatus timestamps are free physics after each step; add0.6−0.002s to align with pre-step neural rows.',protocol=str(OUT/'PROTOCOL.json'))
    write(OUT/'REFLEX_CHECKPOINT.json',checkpoint)
    fresh_reports={}
    for channel in ('Ia','Ib'):
        case='spiking_'+channel+'_pulse';fresh=OUT/'cases'/('syn6_'+case+'_fresh')
        if not fresh.with_suffix('.comparison.json').exists():continue
        anchor=OUT/'cases'/('syn6_'+case);checks={}
        for suffix in ('.npz','.apparatus.npz','.apparatus_extended.npz','.motor_state.npz','.inputs.npz','.reflex.npz','.release_state.npz'):
            with np.load(anchor.with_suffix(suffix)) as x,np.load(fresh.with_suffix(suffix)) as y:
                fields={}
                for key in x.files:
                    if key=='cfg':continue
                    fields[key]=bool(np.array_equal(x[key],y[key],equal_nan=True)) if x[key].dtype.kind in 'fc' else bool(np.array_equal(x[key],y[key]))
                checks[suffix]=fields
        exact=all(all(fields.values()) for fields in checks.values());assert exact
        fresh_reports[case]=dict(case=case,all_saved_arrays_exact=exact,checks=checks,scope='Fresh physical process; complete release state, all event/decoder traces, actual MN/muscle/joint histories. cfg is an output-path-bearing log field and is excluded. No byte-identical ZIP-file claim.')
    if fresh_reports:
        write(OUT/'FRESH_REPEATS.json',fresh_reports)
        write(OUT/'FRESH_REPEAT.json',fresh_reports['spiking_Ia_pulse'])
    np.savez_compressed(OUT/'paired_pulse_effects.npz',**{case+'__'+key:value for case,data in paired.items() for key,value in data.items()})
    plot_cases(cases,paired)
    record=OUT/'EXPERIMENT_RECORD.md'
    text='\n## Completed bounded physical comparisons\n\nAll ten cases completed2.4s free physical integration without a native hyperextension/nonfinite stop. The original long-standing scorer remains unpassed because these runs are shorter than its minimum; this is a reflex demonstration, not curriculum promotion. Each pulse/no-pulse pair is bit-exact before injection. The two source-output-off pairs retain bit-exact downstream neural/physical trajectories despite actual encoder spikes, establishing that the measured response travels through source synapses.\n\n'
    text+='| Case | MN delta range (mV) | Tension delta range (N) | Max ankle delta (deg) |\n| --- | ---: | ---: | ---: |\n'
    for row in summary:text+=f"| {row['case']} | {row['delta_MN_mV_min']:.6g} to {row['delta_MN_mV_max']:.6g} | {row['delta_tension_N_min']:.6g} to {row['delta_tension_N_max']:.6g} | {row['maximum_abs_delta_ankle_deg']:.6g} |\n"
    text+='\nActual native tension is minus the signed engine force cache; actual activation and excitation remain separate saved channels. LIF membrane traces are threshold-normalized, not absolute mV; the decoded afferent and MN voltages use the original relative-rest mV frame. Physics, source muscles, all neural gains, and support settings were frozen. Source timing remains2ms neural/physics and0.2ms LIF substeps. See `REFLEX_CHECKPOINT.json` for exact controls, native-condition bounds, performance, and the saved biological-prefix comparison; see `spiking_reflex_pulses.pdf` for aligned actual events and mechanical responses. No source/controller/manuscript was changed.\n'
    old=record.read_text(encoding='utf-8');marker='\n## Completed bounded physical comparisons\n'
    record.write_text(old.split(marker)[0]+text,encoding='utf-8')
    print(json.dumps(dict(conditions=len(conditions),all_pre_onset_exact=all(r['pre_onset_bit_exact'] for r in summary),source_output_off_pairs_exact=all(r.get('downstream_physics_bit_exact',True) for r in summary),default_physical_prefix_exact=all(r['exact'] for r in default.values()),fresh_repeats=list(fresh_reports))),flush=True)

def plot_cases(cases,paired):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(8,2,figsize=(6.5,11.5),sharex=True);colors={'graded':'#0072b2','spiking':'#d55e00'};col={n:i for i,n in enumerate(COLS)}
    for j,channel in enumerate(('Ia','Ib')):
        ax=axes[:,j];ax[0].set_title(channel+' source-current pulse',fontsize=11)
        for mode in ('graded','spiking'):
            key=mode+'_'+channel+'_pulse';a=cases[key];r=a['reflex'];rows=r['rows'];t=rows[:,0]-ONSET;kwargs=dict(color=colors[mode],ls='-' if mode=='graded' else '--',lw=1.1)
            ax[0].plot(t,rows[:,col['injected_current_nA']],label=mode,**kwargs)
            ax[2].plot(t,rows[:,col[channel+'_voltage_mV']],label=mode,**kwargs)
            ax[3].plot(t,paired[key]['delta_MN_mV'],label=mode,**kwargs)
            ax[4].plot(paired[key]['force_time_neural_s']-ONSET,paired[key]['delta_tension_N'],label=mode,**kwargs)
            ax[5].plot(t,paired[key]['delta_ankle_deg'],label=mode,**kwargs);ax[6].plot(t,paired[key]['delta_knee_deg'],label=mode,**kwargs)
            samples=a['apparatus']['samples'];tf=samples[:,0]+F.R.WARMUP-F.P.DT-ONSET;bw=a['result']['body_weight_N']
            ax[7].plot(tf,(samples[:,5]+samples[:,8])/bw,label=mode,**kwargs)
        off=paired['spiking_'+channel+'_off_pulse']
        ax[3].plot(off['neural_time_s']-ONSET,off['delta_MN_mV'],color='#555555',ls=':',lw=1.,label='outputs off')
        ax[4].plot(off['force_time_neural_s']-ONSET,off['delta_tension_N'],color='#555555',ls=':',lw=1.,label='outputs off')
        for row,key in [(5,'delta_ankle_deg'),(6,'delta_knee_deg')]:ax[row].plot(off['neural_time_s']-ONSET,off[key],color='#555555',ls=':',lw=1.)
        r=cases['spiking_'+channel+'_pulse']['reflex'];k=0 if channel=='Ia' else 1;events=r['event_spikes'][:,k,:];et=r['event_t']-ONSET
        for cell in range(16):ax[1].scatter(et[events[:,cell]>0],np.full(np.count_nonzero(events[:,cell]),cell),marker='|',s=12,color='#222222')
        ax[1].set_ylim(-1,16);ax[1].set_yticks([0,8,15])
        for a in ax:
            a.axvspan(0.,WIDTH,color='#999999',alpha=.12);a.grid(alpha=.18);a.set_xlim(-.18,.4)
        ax[-1].set_xlabel('Time relative to pulse onset (s)')
    labels=['Test current\n(nA)','LIF spike raster\n(cell index)','Afferent readout\n(mV, relative rest)','Soleus MN Δ\n(mV)','Soleus tension Δ\n(N)','Right ankle Δ\n(deg)','Right knee Δ\n(deg)','Total foot force\n(BW)']
    for i,label in enumerate(labels):axes[i,0].set_ylabel(label)
    handles,legend_labels=axes[3,1].get_legend_handles_labels()
    fig.legend(handles,legend_labels,loc='lower center',bbox_to_anchor=(.5,.95),ncol=3,frameon=False)
    fig.subplots_adjust(left=.14,right=.98,top=.92,bottom=.105,hspace=.18,wspace=.21)
    fig.text(.04,.02,'Only soleus-r Ia/Ib input encoders are spiking (16 cells each); RG/PF/IN/MN stay graded.\nΔ denotes pulse minus its matched no-pulse trajectory. Source-output-off controls have exact\nzero downstream response. Native muscles/contact and existing support remain live; no clamp\nafter 0.6 s preparation. Short demonstrations do not pass the separate long-standing gate.',fontsize=10)
    for suffix in ('pdf','svg','png'):fig.savefig(OUT/('spiking_reflex_pulses.'+suffix),dpi=200,bbox_inches='tight',facecolor='white')
    for stem,selected in [('spiking_reflex_neural',[0,1,2,3]),('spiking_reflex_mechanics',[4,5,6,7])]:
        small,aa=plt.subplots(4,2,figsize=(6.5,6.8),sharex=True)
        for ri,source_row in enumerate(selected):
            for j in range(2):
                old,new=axes[source_row,j],aa[ri,j]
                for line in old.lines:new.plot(line.get_xdata(),line.get_ydata(),color=line.get_color(),ls=line.get_linestyle(),lw=line.get_linewidth(),label=line.get_label())
                for collection in old.collections:
                    offsets=collection.get_offsets();new.scatter(offsets[:,0],offsets[:,1],marker='|',s=12,color='#222222')
                new.axvspan(0.,WIDTH,color='#999999',alpha=.12);new.set_xlim(old.get_xlim());new.set_ylim(old.get_ylim());new.grid(alpha=.18)
                if j==0:new.set_ylabel(labels[source_row])
                if ri==0:new.set_title(('Ia' if j==0 else 'Ib')+' source-current pulse',fontsize=11)
                if ri==3:new.set_xlabel('Time relative to pulse onset (s)')
        small.legend(handles,legend_labels,loc='lower center',bbox_to_anchor=(.5,.145),ncol=3,frameon=False)
        small.subplots_adjust(left=.15,right=.98,top=.93,bottom=.25,hspace=.18,wspace=.24)
        small.text(.035,.025,'Gait2392 syn6: only Ia/Ib inputs are spiking (16 cells each).\nBlue: graded source; orange: filtered-spike readout/response.\nΔ = pulse minus matched no-pulse; dotted: source outputs off.\nExisting apparatus; short reflex test, not a standing/curriculum gate.',fontsize=10)
        for suffix in ('pdf','svg','png'):small.savefig(OUT/(stem+'.'+suffix),dpi=200,bbox_inches='tight',facecolor='white')
        plt.close(small)
    plt.close(fig)
    fig,aa=plt.subplots(3,2,figsize=(6.5,6.8),sharex=True)
    for j,channel in enumerate(('Ia','Ib')):
        r=cases['spiking_'+channel+'_pulse']['reflex'];k=0 if channel=='Ia' else 1;t=r['event_t']-ONSET
        for cell,color in [(0,'#0072b2'),(8,'#d55e00')]:
            aa[0,j].plot(t,r['event_v'][:,k,cell],color=color,lw=1.,label='cell '+str(cell))
            aa[1,j].plot(t,r['event_ref'][:,k,cell]*1000,color=color,lw=1.)
        aa[0,j].axhline(1.,color='#777777',ls=':',lw=.8);aa[0,j].set_title(channel+' actual LIF states',fontsize=11)
        for cell in range(16):
            fired=r['event_spikes'][:,k,cell]>0;aa[2,j].scatter(t[fired],np.full(fired.sum(),cell),marker='|',color='#222222',s=15)
        aa[2,j].set_xlabel('Time relative to pulse onset (s)')
        for a in aa[:,j]:a.set_xlim(-.005,.065);a.grid(alpha=.18)
    aa[0,0].set_ylabel('Threshold-normalized\nmembrane voltage');aa[1,0].set_ylabel('Refractory time\nremaining (ms)');aa[2,0].set_ylabel('Actual spike events\n(cell index)')
    hh,ll=aa[0,0].get_legend_handles_labels();fig.legend(hh,ll,loc='lower center',bbox_to_anchor=(.5,.145),ncol=2,frameon=False)
    fig.subplots_adjust(left=.17,right=.98,top=.93,bottom=.25,hspace=.22,wspace=.24)
    fig.text(.035,.025,'Actual 200 µs LIF integration: threshold 1, reset 0, membrane 5 ms, refractory 1 ms.\nVoltages are normalized model states, not absolute biological mV. Sixteen\nidentical deterministic cells may synchronize after rest; no RNG or added\nheterogeneity is claimed. Pulse input is +1 nA, not an imposed muscle stretch.',fontsize=10)
    for suffix in ('pdf','svg','png'):fig.savefig(OUT/('spiking_reflex_membranes.'+suffix),dpi=200,bbox_inches='tight',facecolor='white')
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--case',choices=CASES);p.add_argument('--fresh',action='store_true');p.add_argument('--sequence',action='store_true');p.add_argument('--analyze',action='store_true');p.add_argument('--prefix-audit',action='store_true');a=p.parse_args()
    if a.prepare:prepare()
    elif a.case:run_case(a.case,a.fresh)
    elif a.sequence:
        for case in CASES:
            run_case(case)
            if case=='graded_baseline':prefix_audit()
    elif a.analyze:analyze()
    elif a.prefix_audit:prefix_audit()
    else:p.error('Specify --prepare, --case, or --sequence')
