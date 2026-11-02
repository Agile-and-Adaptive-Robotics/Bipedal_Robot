"""Explicit non-Gait design contracts and opt-in passive-path adapters.

No physics, optimization, Hunt training, or automatic model-wide rescaling.
Units at this interface are seconds, mV, nA, nF, and uS. Native AnimatLab
TimeConst and SNS clock-coded capacitances are adapted explicitly below.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace
from types import FunctionType
import importlib.util
import xml.etree.ElementTree as ET
import numpy as np

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent
REPO = SPINAL.parents[2]
ROBOT = SPINAL.parent/'robot_torso'
VIRTUAL = SPINAL/'virtual_walker_20261031'
SIMULINK = REPO/'Code/Matlab/SNS_Simscape'
EXTERNAL = Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')
RAT = EXTERNAL/'subnetworks/Rat132_bounded_diagnostic/Curated134_ReceptorGeometry'
LEDGERS = HERE/'non_gait_ledgers'
RAW_RECORDS = EXTERNAL/'fsa_design/non_gait_electrical'


def write_raw_record(name,payload):
    """External raw evidence; never replace an independent earlier check."""
    RAW_RECORDS.mkdir(parents=True,exist_ok=True)
    path=RAW_RECORDS/name
    if path.exists():raise FileExistsError(f'Preserve existing evidence: {path}. Choose a new --output-suffix.')
    with path.open('x',encoding='utf-8') as stream:json.dump(payload,stream,indent=2)
    return path


def verify_storage(output_suffix=''):
    """Validate the exact preserved move; do not repeat neural experiments."""
    manifest=json.loads((RAW_RECORDS/'MOVE_MANIFEST.json').read_text(encoding='utf-8-sig'))
    checks=[]
    for row in manifest['verified_files']:
        path=Path(row['current_path']);original=Path(row['original_path'])
        checks.append(dict(path=str(path),sha256_matches=sha(path).upper()==row['sha256'].upper(),bytes_match=path.stat().st_size==row['bytes'],original_raw_absent=not original.exists()))
    result=dict(machine='easteregg2',passed=all(all(row[k] for k in ('sha256_matches','bytes_match','original_raw_absent')) for row in checks),files=checks,
                current_helper_sha256=sha(Path(__file__)),historical_proof_source_hashes='Preserved lineage before output-location-only edits; not a new numerical rerun.',physics_steps=0,training_steps=0)
    if not result['passed']:raise AssertionError('Owned raw-record preservation check failed')
    write_raw_record('STORAGE_VERIFICATION'+output_suffix+'.json',result)
    print(json.dumps(dict(passed=result['passed'],verified_files=len(checks),physics_steps=0,training_steps=0)))


def sha(path):
    path = Path(path)
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def evidence(path, needle=None):
    path = Path(path)
    item = dict(path=str(path), exists=path.exists(), sha256=sha(path))
    if needle and path.exists():
        item['line_matches'] = [i for i,line in enumerate(path.read_text(encoding='utf-8',errors='replace').splitlines(),1) if needle in line]
    return item


def positive(name, value):
    value = float(value)
    if not math.isfinite(value) or value <= 0: raise ValueError(name+' must be finite and positive')
    return value


def finite(name,value):
    value=float(value)
    if not math.isfinite(value):raise ValueError(name+' must be finite')
    return value


def passive_design(tau_s, leak_uS, *, operating_synaptic_uS=0., bias_nA=0.):
    """A declared local RC target, not ionic burst or global-bandwidth design."""
    tau = positive('tau_s',tau_s); leak=positive('leak_uS',leak_uS)
    syn=float(operating_synaptic_uS)
    if syn<0 or not math.isfinite(syn):raise ValueError('Nonnegative operating conductance required')
    return dict(tau_s=tau,leak_uS=leak,operating_synaptic_uS=syn,
                capacitance_nF=1000*tau*(leak+syn),bias_nA=float(bias_nA),
                scope='Local passive pole at the declared conductance operating point, not global network bandwidth.')


def endpoint_design(leak_uS, pre_range_mV, post_baseline_mV, delta_post_mV, reversal_absolute_mV):
    """Exact single-path equilibrium endpoint with fixed baseline current.

    This is NOT a proof of linear transmission over the entire presynaptic
    range, simultaneous-input equivalence, or closed-loop stability.
    """
    gm=positive('leak_uS',leak_uS); span=positive('pre_range_mV',pre_range_mV)
    change=float(delta_post_mV); baseline=float(post_baseline_mV); reversal=float(reversal_absolute_mV)
    denominator=reversal-baseline-change
    if not all(math.isfinite(x) for x in (change,baseline,reversal)) or denominator==0:raise ValueError('Invalid endpoint/reversal')
    g=gm*change/denominator
    if g<0:raise ValueError('Desired endpoint requires negative conductance')
    return dict(gmax_uS=g,endpoint_gain=change/span,delta_post_mV=change,
                reversal_absolute_mV=reversal,post_baseline_mV=baseline,
                warning='Endpoint contract only; verify intermediate inputs, competing paths, and transients separately.')


def source_reset_lif_rate_hz(current_nA,*,rest_mV,threshold_mV,reset_mV,leak_uS,capacitance_nF,refractory_s=0.):
    """Exact uncoupled constant-threshold RC rate with the ACTUAL reset.

    Not Table2adaptiveθ, not a rate assertion for recurrent ionic cells.
    """
    gm=positive('leak_uS',leak_uS);cm=positive('capacitance_nF',capacitance_nF)
    if reset_mV>=threshold_mV:raise ValueError('Reset must be below threshold')
    vinf=rest_mV+finite('current_nA',current_nA)/gm
    if vinf<=threshold_mV:return 0.
    period=cm/(1000*gm)*math.log((vinf-reset_mV)/(vinf-threshold_mV))+refractory_s
    return 1/period


def source_reset_lif_current_nA(rate_hz,*,rest_mV,threshold_mV,reset_mV,leak_uS,capacitance_nF,refractory_s=0.):
    """Inverse of the declared source-reset single-cell rate law (>0Hz)."""
    rate=positive('rate_hz',rate_hz);gm=positive('leak_uS',leak_uS);cm=positive('capacitance_nF',capacitance_nF)
    period=1/rate-refractory_s
    if period<=0:raise ValueError('Desired rate violates refractory limit')
    vinf=threshold_mV+(threshold_mV-reset_mV)/math.expm1(period/(cm/(1000*gm)))
    return gm*(vinf-rest_mV)


def _profile(profile):
    if isinstance(profile,(str,Path)):profile=json.loads(Path(profile).read_text())
    return profile.get('opt_in_profile',profile)


def apply_sns_passive_profile(compiled,names,profile=None,*,clock_unit='s',neural_rate_multiplier=1.):
    """Apply selected existing passive cells/edges; update sparse caches too.

    Optional caller hook AFTER compilation and BEFORE any preparation/state
    capture. Existing production runners do not silently call this function.
    """
    if profile is None:return dict(applied=False,reason='Default-off literal no-op')
    profile=_profile(profile);core=getattr(compiled,'c',compiled);names=list(names)
    rate=positive('neural_rate_multiplier',neural_rate_multiplier)
    if clock_unit not in ('s','ms'):raise ValueError('Explicit SNS compile clock required')
    if len(set(names))!=len(names):raise ValueError('Unique runtime neuron names required')
    if len(names)!=len(core.V):raise ValueError('Runtime/source census mismatch')
    if profile.get('adapter')!='sns_passive':raise ValueError('Wrong adapter')
    updates=[]
    for change in profile.get('neurons',[]):
        if set(change)-{'name','tau_s','bias_nA','operating_synaptic_uS','response_policy'}:raise ValueError('Unsupported neuron field; no opaque rescaling')
        i=names.index(change['name'])
        if core.gated and np.any(core.g_ion[:,i]!=0):raise ValueError('Ionic/NaP cell requires a separate contract')
        tau=positive('tau_s',change['tau_s'])
        gop=finite('operating_synaptic_uS',change.get('operating_synaptic_uS',0.))
        if gop<0:raise ValueError('Nonnegative operating conductance required')
        cm=tau*(core.g_m[i]+gop)*rate*(1000 if clock_unit=='ms' else 1)
        policy=change.get('response_policy','source_euler_continuous_RC')
        if policy=='exact_discrete_pole':
            physical_dt=core.dt/(1000 if clock_unit=='ms' else 1)/rate
            cm=core.dt*(core.g_m[i]+gop)/(-math.expm1(-physical_dt/tau))
        elif policy!='source_euler_continuous_RC':raise ValueError('Unknown passive response policy')
        if 'bias_nA' in change:finite('bias_nA',change['bias_nA'])
        updates.append(('neuron',i,change,cm))
    for change in profile.get('synapses',[]):
        if set(change)-{'source','target','gmax_uS'}:raise ValueError('Unsupported synapse field')
        pre=names.index(change['source']);post=names.index(change['target']);g=float(change['gmax_uS'])
        if not math.isfinite(g) or g<0:raise ValueError('Nonnegative finite conductance required')
        if core.g_max_non[post,pre]==0:raise ValueError('Cannot add a missing/zero-source edge through this adapter')
        hits=None
        if hasattr(compiled,'post'):
            hits=np.flatnonzero((compiled.pre==pre)&(compiled.post==post))
            if len(hits)!=1:raise ValueError('Ambiguous sparse source edge')
        updates.append(('synapse',(post,pre),change,hits))
    proof=[]
    for kind,index,change,value in updates:
        if kind=='neuron':
            old=float(core.c_m[index]);core.c_m[index]=value;core.time_factor_membrane[index]=core.dt/value
            if 'bias_nA' in change:core.i_b[index]=float(change['bias_nA'])
            gop=change.get('operating_synaptic_uS',0.);physical_dt=core.dt/(1000 if clock_unit=='ms' else 1)/rate;decay=1-core.dt*(core.g_m[index]+gop)/core.c_m[index]
            proof.append(dict(kind=kind,name=names[index],old_clock_capacitance=old,new_clock_capacitance=float(core.c_m[index]),physical_leak_tau_s=float(core.c_m[index]/core.g_m[index]/rate/(1000 if clock_unit=='ms' else 1)),physical_operating_continuous_RC_tau_s=float(core.c_m[index]/(core.g_m[index]+gop)/rate/(1000 if clock_unit=='ms' else 1)),physical_operating_discrete_tau_s=float(-physical_dt/math.log(decay)) if 0<decay<1 else None,response_policy=change.get('response_policy','source_euler_continuous_RC')))
        else:
            post,pre=index;old=float(core.g_max_non[post,pre]);core.g_max_non[post,pre]=change['gmax_uS']
            if value is not None:compiled.gmax[value]=change['gmax_uS']
            proof.append(dict(kind=kind,source=names[pre],target=names[post],old_gmax_uS=old,new_gmax_uS=float(core.g_max_non[post,pre]),sparse_cache_updated=value is not None))
    return dict(applied=bool(updates),consumed_fields=proof,clock_unit=clock_unit,neural_rate_multiplier=rate,added_neurons=0,added_edges=0,ionic_fields_changed=False)


def apply_native_passive_profile(controller,profile=None,*,neural_rate_multiplier=1.):
    """Source IntegrateFireSim normalized kernel; preserve Size/Ca/gate laws.

    Native edge conductance is the actual normalized kernel coefficient.
    Do not label Size or relative SynAmp as measured physical capacitance/G.
    """
    if profile is None:return dict(applied=False,reason='Default-off literal no-op')
    profile=_profile(profile);names=list(controller.names);rate=positive('neural_rate_multiplier',neural_rate_multiplier)
    if profile.get('adapter')!='native_passive':raise ValueError('Wrong adapter')
    updates=[]
    for change in profile.get('neurons',[]):
        if set(change)-{'name','tau_s','tonic_nA','operating_synaptic_normalized'}:raise ValueError('Unsupported native neuron field')
        if names.count(change['name'])!=1:raise ValueError('Missing/ambiguous native neuron name')
        i=names.index(change['name'])
        if controller.params[i,3]!=0:raise ValueError('Retain ionic source mechanism; pure passive profile rejected')
        if 'tonic_nA' in change:finite('tonic_nA',change['tonic_nA'])
        gop=finite('operating_synaptic_normalized',change.get('operating_synaptic_normalized',0.))
        if gop<0:raise ValueError('Nonnegative operating conductance required')
        updates.append(('neuron',i,change,1000*positive('tau_s',change['tau_s'])*rate*(1+gop)))
    for change in profile.get('synapses',[]):
        if set(change)-{'source','target','gmax_normalized'}:raise ValueError('Unsupported native synapse field')
        if names.count(change['source'])!=1 or names.count(change['target'])!=1:raise ValueError('Ambiguous native edge names')
        pre=names.index(change['source']);post=names.index(change['target'])
        hits=np.flatnonzero((controller.edges[:,0]==pre)&(controller.edges[:,1]==post))
        if len(hits)!=1:raise ValueError('Missing/ambiguous source edge; cannot invent topology')
        if controller.edges[hits[0],2]==0:raise ValueError('Cannot enable sourcezero path without explicit architecture contract')
        g=float(change['gmax_normalized'])
        if not math.isfinite(g) or g<0:raise ValueError('Invalid conductance')
        updates.append(('synapse',int(hits[0]),change,g))
    proof=[]
    for kind,index,change,value in updates:
        if kind=='neuron':
            old=float(controller.params[index,1]);controller.params[index,1]=value
            if 'tonic_nA' in change:controller.tonic[index]=float(change['tonic_nA'])
            proof.append(dict(kind=kind,name=names[index],old_tau_ms=old,new_tau_ms=float(controller.params[index,1]),physical_operating_tau_s=float(controller.params[index,1]/1000/rate/(1+change.get('operating_synaptic_normalized',0.)))))
        else:
            old=float(controller.edges[index,2]);controller.edges[index,2]=value
            proof.append(dict(kind=kind,source=change['source'],target=change['target'],edge_index=index,old_gmax_normalized=old,new_gmax_normalized=float(controller.edges[index,2])))
    return dict(applied=bool(updates),consumed_fields=proof,neural_rate_multiplier=rate,added_neurons=0,added_edges=0,Size_or_ionic_fields_changed=False)


def clock_contract(duration_multiplier=1.,parameter_time_multiplier=1.):
    c=positive('duration_multiplier',duration_multiplier);p=positive('parameter_time_multiplier',parameter_time_multiplier)
    return dict(definition='t_biological=t_physical/c; parameterTimeConstants=p*sourceTimeConstants',duration_multiplier=c,parameter_time_multiplier=p,
                biological_rate_per_physical_rate=1/c,passive_tau_physical_per_source_tau=c*p,
                equations=['tau_physical=C_source/(G*lambda), lambda=1/c','C_source_for_target_physical_tau=tau_target*G*lambda (with clock-unit factor)','NativeTimeConst_ms_for_target=1000*tau_target*lambda*(1+declared_operating_synaptic_ratio)'],
                warning='A biologicalclock change also changes sensory/MN filters, spikecarrier/refractory/readout and startupwindows if inside the same wrapper; it is not a computationalspeedup. Keep physics/actuator clocks explicit.')


def carrier_clock_parameters(*,physical_rate_hz,physical_tau_m_s,physical_refractory_s,physical_filter_s,neural_rate_multiplier=1.):
    """Selected sensory-carrier clock compensation, not global brain scaling."""
    rate=positive('neural_rate_multiplier',neural_rate_multiplier)
    return dict(rate_max=positive('physical_rate_hz',physical_rate_hz)/rate,
                tau_m=positive('physical_tau_m_s',physical_tau_m_s)*rate,
                refractory_s=positive('physical_refractory_s',physical_refractory_s)*rate,
                tau_filter=positive('physical_filter_s',physical_filter_s)*rate,
                neural_rate_multiplier=rate)


def native_input_contract(controller,source,tau_physical_s=.01):
    """Derive admissible current spans from actual pure-input adapter nodes.

    Preserve Size, rests, physical receptor maps, output windows, and graph.
    The10ms physical pole is an explicit engineering candidate, not a claim
    that every biological afferent universally has this time constant.
    """
    root=ET.parse(source).getroot();module=next(x for x in root.iter('NeuralModule') if x.findtext('ModuleName')=='IntegrateFireSim')
    nodes=list(module.iter('Neuron'));ids={n.findtext('ID'):i for i,n in enumerate(nodes)}
    inputs={}
    for adapter in root.iter('Adapter'):
        if adapter.findtext('Enabled')=='False' or adapter.findtext('Type')!='PhysicalToNode':continue
        i=ids.get(adapter.findtext('TargetID'))
        if i is not None:inputs.setdefault(i,[]).append(adapter.findtext('SourceDataType'))
    contracts=[];overrides=[]
    for i,kinds in inputs.items():
        if controller.params[i,3]!=0 or controller.tonic[i]!=0 or np.any(controller.edges[:,1]==i):continue
        outgoing=controller.edges[controller.edges[:,0]==i]
        if not len(outgoing):continue
        rest=controller.params[i,0];lo=outgoing[:,3];hi=outgoing[:,4]
        if np.any(hi<=lo) or np.max(hi)<=rest:continue
        span=float(np.max(hi)-rest);size=float(controller.params[i,2]);name=controller.names[i]
        if controller.names.count(name)!=1:continue
        contracts.append(dict(name=name,actual_physical_kinds=sorted(set(kinds)),input_current_nA=[0.,size*span],output_voltage_mV=[float(rest),float(np.max(hi))],DC_gain_mV_per_source_nA=1/size,physical_leak_tau_s=tau_physical_s,physical_leak_bandwidth_Hz=1/(2*math.pi*tau_physical_s),source_tau_ms=float(controller.params[i,1]),incoming_central_edges=0,physical_adapter_maps_unchanged=True,branch_windows=[dict(target=controller.names[int(e[1])],lower_mV=float(e[3]),upper_mV=float(e[4]),baseline_fraction=float(np.clip((rest-e[3])/(e[4]-e[3]),0,1))) for e in outgoing],heterogeneous_windows_preserved=True))
        overrides.append(dict(name=name,tau_s=tau_physical_s))
    return dict(adapter='native_passive',adopted=False,engineering_candidate=True,neurons=overrides,synapses=[]),contracts


def native_motor_contract(controller,source):
    """Motor Vm domain from the actual source sigmoid and SI adapter.

    Target5..95% is of positive ACTIVE sigmoid plateau before length/series
    effects, never a claim of total tension or measured biological activation.
    Preserve each source MN's physical membrane time while cadence changes.
    """
    root=ET.parse(source).getroot();module=next(x for x in root.iter('NeuralModule') if x.findtext('ModuleName')=='IntegrateFireSim')
    nodes=list(module.iter('Neuron'));ids={n.findtext('ID'):i for i,n in enumerate(nodes)};bodies={b.findtext('ID'):b for b in root.iter('RigidBody')}
    contracts=[];overrides=[]
    for adapter in root.iter('Adapter'):
        if adapter.findtext('Enabled')=='False' or adapter.findtext('Type')!='NodeToPhysical':continue
        i=ids.get(adapter.findtext('SourceID'));body=bodies.get(adapter.findtext('TargetID'))
        if i is None or body is None or controller.params[i,3]!=0:continue
        gain=adapter.find('Gain');coeff=[float(gain.findtext(key,'0')) for key in 'ABCD'];stim=body.find('StimulusTension')
        if stim is None or coeff[0]!=0 or coeff[1]!=0 or coeff[2]==0:continue
        a,b,c,d=[float(stim.findtext(key,'0')) for key in 'ABCD']
        if b<=0 or c==0 or b+d<=0:continue
        name=controller.names[i];points=[]
        for active in (.05,.5,.95):
            sigmoid=((b+d)*active-d)/b
            if not 0<sigmoid<1:break
            muscle_v=a-math.log(1/sigmoid-1)/c;source_v=(muscle_v-coeff[3])/coeff[2]
            limits= gain.findtext('UseLimits')=='True'
            admitted=not limits or float(gain.findtext('LowerLimit'))<=source_v<=float(gain.findtext('UpperLimit'))
            recovered=(((1/(1+math.exp(c*(a-muscle_v))))*b+d)/(b+d))
            points.append(dict(active_plateau_fraction=active,neural_voltage_mV=1000*source_v,muscle_stimulus_volts=muscle_v,adapter_input_volts=source_v,within_source_adapter_limits=admitted,inverse_roundtrip_error=abs(active-recovered)))
        if len(points)!=3:continue
        tau=float(controller.params[i,1]/1000)
        contracts.append(dict(name=name,actuator=body.findtext('Name'),source_sigmoid=dict(A=a,B=b,C=c,D=d),source_adapter_polynomial=coeff,source_voltage_conversion='NeuralmV*0.001 BEFORE polynomial; sourcevoltsandneuronal mVnotinterchangeable.',voltage_active_domain=points,physical_leak_tau_s=tau,scope='Source-preserving motor voltage/currentconditioning; actualforcealsodependson length/velocity/series/passivetissue. Notforcefloor or biologicalactivationproof.'))
        overrides.append(dict(name=name,tau_s=tau))
    return overrides,contracts


class PhysicalCarrierPort:
    """Opt-in source sensory encoder on its own PHYSICAL clock.

    Reuses the existing source encode_step unchanged. Caller must disable
    the corresponding old encoder, preserve source output bindings, and
    hold this decoded signal between physical carrier ticks. No body or
    central-controller hook is installed automatically.
    """
    def __init__(self,source_names,*,population_size=32,physical_dt_s=.0002,
                 physical_rate_hz=200.,physical_tau_m_s=.005,
                 physical_refractory_s=.001,physical_filter_s=.01):
        sys.path.insert(0,str(VIRTUAL))
        from spiking_sensory import encode_step
        self.encoder=encode_step;self.names=list(source_names);self.N=int(population_size)
        if self.N<1 or not self.names:raise ValueError('Nonempty named populations required')
        self.dt=positive('physical_dt_s',physical_dt_s);self.rate=positive('physical_rate_hz',physical_rate_hz)
        self.tau_m=positive('physical_tau_m_s',physical_tau_m_s);self.refractory_s=positive('physical_refractory_s',physical_refractory_s)
        self.tau=np.full(len(self.names),positive('physical_filter_s',physical_filter_s))
        if self.rate*self.refractory_s>=1:raise ValueError('Rate exceeds refractory ceiling')
        self.v=np.tile(np.linspace(.01,.99,self.N),(len(self.names),1));self.ref=np.zeros_like(self.v)
        self.decoded=np.zeros(len(self.names));self.events=np.zeros_like(self.v,dtype=np.int8)
        self.counts=np.zeros(len(self.names),dtype=np.int64);self.total=np.zeros_like(self.counts)
        self.pending=0.;self.steps=0
    def advance(self,physical_interval_s,normalized_target):
        interval=positive('physical_interval_s',physical_interval_s);target=np.asarray(normalized_target,dtype=float)
        if target.shape!=self.decoded.shape or not np.isfinite(target).all() or np.any((target<0)|(target>1)):raise ValueError('Declared normalized input domain is[0,1]; no hidden clipping')
        self.pending+=interval;ticks=int(math.floor((self.pending+1e-14)/self.dt));self.pending-=ticks*self.dt
        if abs(self.pending)<1e-14:self.pending=0.
        self.counts[:]=0
        for _ in range(ticks):self.encoder(target,self.v,self.ref,self.decoded,self.events,self.counts,self.dt,self.rate,self.tau_m,self.tau,self.refractory_s)
        self.steps+=ticks;self.total+=self.counts
        return self.decoded.copy(),self.counts.copy()
    def export_state(self):
        return dict(names=np.array(self.names),configuration=np.r_[self.N,self.dt,self.rate,self.tau_m,self.refractory_s,self.tau],
                    v=self.v.copy(),refractory=self.ref.copy(),decoded=self.decoded.copy(),last_events=self.events.copy(),
                    last_counts=self.counts.copy(),total_counts=self.total.copy(),pending_physical_s=np.array(self.pending),
                    physical_steps=np.array(self.steps,dtype=np.int64),rng_state=np.empty(0,dtype=np.uint8))
    def import_state(self,state):
        expected=self.export_state()
        if set(state)!=set(expected):raise ValueError('Incomplete carrier state')
        for key in ('names','configuration'):
            if not np.array_equal(state[key],expected[key]):raise ValueError('Carrier identity/configuration mismatch')
        if np.asarray(state['rng_state']).size:raise ValueError('Deterministic source carrier has noRNG')
        for key,attribute in [('v','v'),('refractory','ref'),('decoded','decoded'),('last_events','events'),('last_counts','counts'),('total_counts','total')]:
            if state[key].shape!=expected[key].shape or state[key].dtype!=expected[key].dtype:raise ValueError('State shape/dtype mismatch')
            getattr(self,attribute)[:]=state[key]
        self.pending=float(state['pending_physical_s']);self.steps=int(state['physical_steps'])


def deng14_source_module():
    name='aarl_fsa_fixed_deng14_source'
    if name in sys.modules:return sys.modules[name]
    spec=importlib.util.spec_from_file_location(name,SPINAL/'deng_cpg_ode.py');module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module)
    return module


def deng14_profile_function(profile=None):
    """Selective passive-C broadcasting through the ORIGINAL ODE code.

    Original simulate() already divides its voltage vector by global CM;
    an explicit per-cell array permits selected passive capacitances without
    changing that function, the graph, NaP cells, or h-gate law on disk.
    """
    source=deng14_source_module()
    if profile is None:return source.simulate,dict(applied=False,default_function_identity=True)
    profile=_profile(profile)
    if profile.get('adapter')!='deng14_passive' or profile.get('synapses'):raise ValueError('OnlyselectivepassiveCprofileisimplemented')
    cm=np.full(len(source.NAMES),source.CM);changes=[]
    for change in profile.get('neurons',[]):
        if set(change)-{'name','tau_s'}:raise ValueError('UnsupportedDeng14field')
        i=source.NAMES.index(change['name'])
        if i in source.NA:raise ValueError('OriginalNaPcapacitance/gatesmustberetained')
        cm[i]=1000*positive('tau_s',change['tau_s'])*source.GM
        changes.append(dict(name=change['name'],old_C_nF=source.CM,new_C_nF=float(cm[i]),physical_leak_tau_s=float(cm[i]/source.GM/1000)))
    namespace=dict(source.simulate.__globals__);namespace['CM']=cm
    fn=FunctionType(source.simulate.__code__,namespace,'simulate_deng14_selective_passive',source.simulate.__defaults__,source.simulate.__closure__)
    assert np.all(cm[list(source.NA)]==source.CM)
    return fn,dict(applied=bool(changes),consumed_fields=changes,NaP_names=[source.NAMES[i] for i in sorted(source.NA)],NaP_C_nF=source.CM,h_gate_ms=source.TAU_H,source_function_code_unchanged=fn.__code__ is source.simulate.__code__,source_globals_unchanged=isinstance(source.CM,float))


def run_deng14_passive_profile(profile=None,*,duration_s=20.,dt_s=.0001,pulse_time_s=.1,pulse_duration_s=.001,pulse_current_nA=10.):
    fn,proof=deng14_profile_function(profile)
    voltage,h=fn('fixed',t_end=duration_s*1000,dt=dt_s*1000,i_pulse=pulse_current_nA,t_pulse=pulse_time_s*1000,pulse_ms=pulse_duration_s*1000)
    return np.arange(len(voltage))*dt_s,voltage,h,proof


def simulink_mask_recipe(block,*,tau_s,leak_uS):
    """Physical-unit passive block patch only; this does not execute MATLAB."""
    d=passive_design(tau_s,leak_uS)
    return dict(block=str(block),mask_values=dict(Gm=leak_uS,Cm=d['capacitance_nF']),
                source_derivative='1000/Cm in seconds-based Simulink',runtime_consumption_verified=False,
                warning='Apply to a separate passive block/model via set_param; no whole-model or ionic rescaling.')


def simulink_event_mask_recipe(block,*,mean_current_nA,post_mV,reversal_mV,event_rate_hz=100.,tau_s=.01,cap_uS=.12,hybrid=True):
    """Derive the actual additive source primitive in its uncapped regime."""
    drive=float(reversal_mV-post_mV);rate=positive('event_rate_hz',event_rate_hz);tau=positive('tau_s',tau_s);cap=positive('cap_uS',cap_uS)
    if drive==0 or mean_current_nA/drive<=0:raise ValueError('Desired current and reversal drive must have the same sign')
    increment=float(mean_current_nA/drive/(rate*tau));peak=increment/(-math.expm1(-1/(rate*tau)))
    if peak>=cap:raise ValueError('Mean-current linear contract enters cap regime; declare the nonlinear event contract separately')
    values=dict(gmax=cap,ginc=increment,tau_syn=1000*tau,Esyn=reversal_mV)
    if hybrid:values['ThrPre']=-45.
    return dict(block=block,source_primitive='additive/capped rising-event conductance, not reset-to-peak',values=values,mean_current_nA=mean_current_nA,operating_voltage_mV=post_mV,event_rate_hz=rate,steady_peak_uS=peak,physical_voltage_frame='absolute mV source Simulink',adopted=False)


def native_inventory(path):
    root=ET.parse(path).getroot()
    module=next(x for x in root.iter('NeuralModule') if x.findtext('ModuleName')=='IntegrateFireSim')
    nodes=list(module.iter('Neuron'));ids={n.findtext('ID'):i for i,n in enumerate(nodes)}
    def num(node,key,default=0.):return float(node.findtext(key,str(default)))
    fields=('RestingPot','TimeConst','Size','GMaxCa','CaActivation/MidPoint','CaActivation/Slope','CaActivation/TimeConstant','CaDeactivation/MidPoint','CaDeactivation/Slope','CaDeactivation/TimeConstant')
    params=np.array([[num(n,key) for key in fields]+[num(module,'CaEquilPot')] for n in nodes])
    types={t.findtext('ID'):t for t in module.findall('Synapses/NonSpikingSynapses/SynapseType')}
    edges=[]
    for edge in module.iter('Connexion'):
        if edge.findtext('Type')!='1':continue
        syn=types[edge.findtext('SynapseTypeID')]
        edges.append([ids[edge.findtext('SourceID')],ids[edge.findtext('TargetID')]]+[num(syn,key) for key in ('SynAmp','ThreshV','SaturateV','Equil')])
    c=SimpleNamespace(names=[n.findtext('Name') for n in nodes],params=params,edges=np.asarray(edges),tonic=np.array([num(n,'TonicStimulus') for n in nodes]),dt=num(module,'TimeStep'))
    return c


def actual_native_controller(family):
    """Construct designated live neural/adapter objects; never a body model."""
    sys.path.insert(0,str(VIRTUAL))
    if family=='public_rat120_precursor':
        source=REPO/'Neuromechanical_Models/Deng_Biarticular_Rat/Biarticular_net_Standalone.asim'
        audit=source.parent/'build/runs/Biarticular_net/IMPORT_AUDIT.json'
        import native_kernel as K
        c=K.NativeController(source,json.loads(audit.read_text())['muscles'])
        assert len(c.names)==120 and len(c.edges)==238
        return c,source,dict(constructor='NativeController:original public120',source_sha256=sha(source),import_audit=str(audit),physical_input_ports=len(c.inputs),physical_motor_ports=len(c.outputs),source_original_public_not_curated=True,runtime_neural_time_scale=1.,physical_clock_scale=1.)
    if family=='rat134_GA_receptor_variant':
        scripts=RAT/'scripts';sys.path.insert(0,str(scripts))
        import rat134_receptor_runtime as runtime
        c=runtime.variant_controller()
        assert len(c.names)==134 and c.consumed_model_sha256==sha(RAT/'Rat134_ReceptorGeometry.asim')
        return c,RAT/'Rat134_ReceptorGeometry.asim',dict(constructor=str(scripts/'rat134_receptor_runtime.py')+'::variant_controller',source_sha256=c.consumed_model_sha256,physical_input_ports=len(c.inputs),physical_motor_ports=len(c.outputs),actual_curated_adapter=True)
    source=REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'
    cfg=json.loads((VIRTUAL/'GENTLE_HIP_0P6_FLEX_0P5.json').read_text());assert sha(source)==cfg['expected_source_sha256']
    meta=json.loads((EXTERNAL/'spiking/W2L_VirtualWalker_Newer/IMPORT_AUDIT.json').read_text());bank=meta['muscles']
    import native_kernel as K
    if family=='W2L704_hybrid_spiking':
        import spiking_central_runtime as S
        S.SpikingCentralController.configuration=json.loads((EXTERNAL/'spiking/CENTRAL_RESET_REST_CANDIDATE.json').read_text())
        S.SpikingCentralController.sensory_spiking=True;S.SpikingCentralController.clock_scale=1.;S.SpikingCentralController.balance_configuration=None
        c=S.SpikingCentralController(source,bank,contact_si=True,ankle_pf_crossed=cfg['ankle_pf_crossed']);c.apply_output_gain()
    elif family=='W2L108_graded':c=K.NativeController(source,bank,contact_si=True,ankle_pf_crossed=cfg['ankle_pf_crossed'])
    else:raise ValueError('Unknown native anchor family')
    from recruitment import apply_recruitment
    edits=apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale'])
    return c,source,dict(constructor='SpikingCentralController' if family.startswith('W2L704') else 'NativeController',source_sha256=sha(source),configuration=str(VIRTUAL/'GENTLE_HIP_0P6_FLEX_0P5.json'),contact_si=True,authored_clock_already_in_source=2,runtime_neural_time_scale=1.,physical_clock_scale=1.,recruitment=edits,physical_input_ports=len(c.inputs),physical_motor_ports=len(c.outputs))


def native_feedback_budget(c,budget_nA=.5):
    """New explicit first-entry current contract; retained source-side patterns."""
    kinds={int(row[0]):set() for row in c.inputs}
    for row in c.inputs:kinds[int(row[0])].add(('Ia','II','Ib')[int(row[2])] if isinstance(row[2],(int,np.integer)) else str(row[2]))
    groups={}
    for k,e in enumerate(c.edges):
        pre,post=int(e[0]),int(e[1]);op=float(c.params[post,0]+2.5);drive=float(e[5]-op)
        if pre not in kinds or not any(x in c.names[post].upper() for x in ('RG','PF')) or e[2]==0 or drive==0:continue
        groups.setdefault((post,tuple(sorted(kinds[pre])),1 if drive>0 else -1),[]).append(k)
    contracts=[];changes=[]
    for (post,kind,sign),indices in groups.items():
        op=float(c.params[post,0]+2.5);maximum=sum(float(c.edges[k,2])*abs(float(c.edges[k,5])-op)*float(c.params[post,2]) for k in indices);factor=budget_nA/maximum
        contracts.append(dict(target=c.names[post],source_kinds=list(kind),operating_mV=op,old_full_window_current_bound_nA=sign*maximum,new_reference_budget_nA=sign*budget_nA,original_relative_weights_retained=True,scale=factor,indices=indices))
        for k in indices:
            e=c.edges[k];changes.append(dict(source=c.names[int(e[0])],target=c.names[int(e[1])],gmax_normalized=float(e[2]*factor)))
    return dict(adapter='native_passive',neurons=[],synapses=changes),contracts


def native_conditional_input_checks(c,indices,max_currents=None):
    """Retained recurrent inputs have conditional gain, not isolated gain."""
    sys.path.insert(0,str(VIRTUAL));import native_kernel as K
    results=[]
    for post in indices:
        incoming=c.edges[c.edges[:,1]==post].copy()
        if c.params[post,3]!=0:continue
        errors=[];tau_effective=[]
        upper=float((max_currents or {}).get(post,20.))
        for level in (0.,.5,1.):
            for current in (0.,.5*upper,upper):
                p=c.params.copy();p[:,3]=0.;v=p[:,0].copy()
                for e in incoming:v[int(e[0])]=e[3]+level*(e[4]-e[3])
                prior=v.copy();applied=np.zeros(len(v));applied[post]=current;m=np.zeros(len(v));h=m.copy();newv=v.copy();gs=m.copy();gi=m.copy()
                K.neural_steps(v,m,h,newv,gs,gi,p,incoming,c.tonic.copy(),applied,.2,1)
                total=1+gs[post];tau=p[post,1]/1000/total;decay=math.exp(-.0002/tau);drive=(c.tonic[post]+current)/p[post,2]+gi[post]
                expected=p[post,0]+(prior[post]-p[post,0])*decay+drive*(1-decay)/total
                errors.append(abs(float(v[post])-expected));tau_effective.append(float(tau))
        results.append(dict(name=c.names[post],input_domain_nA=[0.,upper],desired_function='Conditional source law: dV/dt=(rest-V+(tonic+I)/Size+sum(g*(E-V)))/tau; input gain1/(Size*(1+sumg)), no background cancellation.',source_membrane_tau_s=float(c.params[post,1]/1000),effective_pole_range_s=[min(tau_effective),max(tau_effective)],incoming_paths=len(incoming),maximum_actual_native_step_error_mV=max(errors),new_parameters=False))
    return results


def native_ionic_mechanism_checks(c):
    """Preserved nonlinear Ca mechanism; no pure-GLIF substitution."""
    sys.path.insert(0,str(VIRTUAL));import native_kernel as K
    rows=[]
    for i,p in enumerate(c.params):
        if p[3]==0:continue
        errors=[];bounded=True;m_inf=[];h_inf=[];tau=[]
        for voltage in (p[0]-20,p[0],p[0]+10,p[0]+30):
            z=math.exp(-p[5]*(voltage-p[4]));mi=1/(1+z);tm=mi*math.sqrt(z)*p[6]
            z=.5*math.exp(-p[8]*(voltage-p[7]));hi=1/(1+z);th=hi*math.sqrt(z)*p[9]
            m_inf.append(mi);h_inf.append(hi);tau.extend((tm/1000,th/1000))
            for level in (0.,.5,1.):
                v=np.array([voltage]);m=np.array([level]);h=np.array([level]);nv=v.copy();gs=np.zeros(1);gi=gs.copy();params=p[None,:].copy()
                vr=v.copy();mr=m.copy();hr=h.copy();nvr=nv.copy();gsr=gs.copy();gir=gi.copy()
                K.neural_steps(v,m,h,nv,gs,gi,params,np.empty((0,6)),np.zeros(1),np.ones(1),.2,1)
                K.neural_steps.py_func(vr,mr,hr,nvr,gsr,gir,params,np.empty((0,6)),np.zeros(1),np.ones(1),.2,1)
                errors.append(float(max(abs(v-vr).max(),abs(m-mr).max(),abs(h-hr).max())));bounded=bounded and bool(0<=m[0]<=1 and 0<=h[0]<=1 and np.isfinite(v).all())
        rows.append(dict(name=c.names[i],mechanism='Source graded voltage-dependent Ca m/h; retained within adaptiveCaLIF if encoded',m_increases_with_voltage=bool(np.all(np.diff(m_inf)>0)),h_decreases_with_voltage=bool(np.all(np.diff(h_inf)<0)),positive_physical_gate_times=bool(min(tau)>0),gate_time_range_s=[min(tau),max(tau)],state_bounds_preserved=bounded,compiled_reference_error=max(errors),Table2_pureGLIF_applicable=False))
    return rows


def native_spiking_port_checks(c):
    """Actual704 carrier/readout modules, not unused graded-shadow cells."""
    assert c.clock_scale==1. and c.population_size==16 and len(c.encoded)==12 and len(c.sensors.indices)==16
    state=c.export_state();carrier=[]
    for level in (0.,.25,.5,1.):
        c.import_state(state);events=np.zeros(16,dtype=np.int64);tail=[]
        for step in range(2000):
            c.applied[c.sensors.indices]=20*c.params[c.sensors.indices,2]*level
            c.sensors.prepare();c.sensors.advance(step*c.dt);events+=c.sensors.counts
            if step>=1500:tail.append(c.v[c.sensors.indices]-c.params[c.sensors.indices,0])
        rates=events/(32*2000*c.dt);requested=200*level
        carrier.append(dict(normalized_input=level,requested_Hz=requested,minimum_actual_rate_Hz=float(rates.min()),maximum_actual_rate_Hz=float(rates.max()),maximum_relative_rate_error=float(np.max(abs(rates-requested))/requested) if requested else float(np.max(abs(rates))),mean_tail_proxy_span_mV=float(np.mean(tail)),actual_events=int(events.sum())))
    c.import_state(state);maximum=0.;events=0;adapt_error=0.;event_identity=True;gate_bounds=True
    for step in range(100):
        previous=c.decoded.copy();before=c.counts.copy();prior_adapt=c.adapt.copy();c.applied[c.encoded]=2.
        c.tick(c.v,c.m,c.h,c.newv,c.gs,c.gi,c.lv,c.lm,c.lh,c.ref,c.adapt,c.decoded,c.counts)
        delta=c.counts-before;events+=int(delta.sum());expected=previous*np.exp(-c.dt/c.settings[3])+delta/(c.population_size*c.settings[4]*c.settings[3]);maximum=max(maximum,float(np.max(abs(c.decoded-expected))))
        new_event=np.isclose(c.ref,.001,atol=1e-16,rtol=0);event_identity=event_identity and np.array_equal(new_event.sum(axis=1),delta)
        taus=c.settings[6]*np.where(c.pf_mask,c.settings[9],1.);jumps=c.settings[5]*np.where(c.pf_mask,c.settings[8],1.)
        expected_adapt=prior_adapt*np.exp(-c.dt/taus[:,None])+new_event*jumps[:,None];adapt_error=max(adapt_error,float(np.max(abs(c.adapt-expected_adapt))))
        gate_bounds=gate_bounds and bool(np.all((c.lm>=0)&(c.lm<=1)) and np.all((c.lh>=0)&(c.lh<=1)))
    c.import_state(state);restored=c.export_state();exact=set(restored)==set(state) and all(np.array_equal(restored[k],state[k]) for k in state)
    return dict(actual_spiking_cells=12*16+16*32,remaining_graded=80,physical_clock_scale=1.,carrier_populations=16,carrier_cells_each=32,carrier_membrane_tau_s=.005,carrier_filter_tau_s=c.sensors.tau.tolist(),carrier_rate_target_hz=[0.,200.],declared_rate_error_tolerance=.1,carrier_checks=carrier,central_populations=12,central_readout_normalizer_Hz=float(c.settings[4]),normalizer_is_not_asserted_Fmax=True,central_filter_tau_s=float(c.settings[3]),actual_central_events=events,central_filter_equation_error=maximum,adaptation_equation_error=adapt_error,actual_event_count_identity=event_identity,local_Ca_gate_bounds=gate_bounds,all_state_restored=exact,passed=all(x['maximum_relative_rate_error']<=.1 for x in carrier) and maximum<1e-12 and adapt_error<1e-12 and event_identity and gate_bounds and events>0 and exact)


def native_adapter_consumption_checks(c):
    """Actual physical-input API with synthetic arrays and zero time advance."""
    count=len(c.bank) if hasattr(c,'bank') else len(c.ibgain);errors=[]
    def value(coeff,x,limits=None):
        if limits is not None and limits[0]>0:
            if x<limits[1]:return limits[3]
            if x>limits[2]:return limits[4]
        return ((coeff[0]*x+coeff[1])*x+coeff[2])*x+coeff[3]
    for level in (0.,.25,1.):
        ia=np.arange(count)*.01+level;ii=np.arange(count)*.02+level;ib=np.arange(count)*.03+level
        expected=np.zeros(len(c.names));contacts={row[4]:level+.1 for row in c.inputs if isinstance(row[2],str) and row[2] in ('ContactCount','Contacts')}
        for row in c.inputs:
            target,body,kind,coeff=row[:4]
            if isinstance(kind,(int,np.integer)):
                x=(ia,ii,ib)[int(kind)][body];expected[target]+=value(coeff,x,row[4])*1e9
            else:
                x=ia[body] if kind=='Ia' else ii[body] if kind=='II' else ib[body] if kind=='Ib' else ib[body]/c.ibgain[body] if kind=='Tension' else contacts[row[4]]
                expected[target]+=value(coeff,x)*1e9
        if hasattr(c,'input_coeff'):actual=c.feedback(ia,ii,ib).copy()
        else:
            c.step(1000.,ia,ii,ib,contacts,0.,True);actual=c.applied.copy()
            if getattr(c,'sensors',None) is not None:actual[c.sensors.indices]=c.sensors.current
        errors.append(float(np.max(abs(actual-expected))))
    return dict(synthetic_physical_array_cases=3,source_actual_API='RatController.feedback' if hasattr(c,'input_coeff') else type(c).__name__+'.step(interval=0)',actual_cached_mapping_vs_source_tuple_error_nA=max(errors),sensory_spiking_dispatch='Actual sensor.current receives its16inputs and native graded applied is cleared; both destinations checked together' if getattr(c,'sensors',None) is not None else None,neural_time_advance_s=0.,physics_steps=0)


def native_role_design(family,output_suffix=''):
    """Electrical prerequisites on actual designated neural/adapter objects."""
    c,source,anchor=actual_native_controller(family)
    sys.path.insert(0,str(VIRTUAL));import native_kernel as K
    original_params=c.params.copy();original_edges=c.edges.copy();original_tonic=c.tonic.copy()
    physical_cache={key:getattr(c,key).copy() for key in ('input_coeff','input_limits','output_coeff','output_limits') if hasattr(c,key)}
    assert not apply_native_passive_profile(c,None)['applied']
    profile,sensors=native_input_contract(c,source);motors,motor_contracts=native_motor_contract(c,source)
    spike_inputs=set(int(i) for i in c.sensors.indices) if getattr(c,'sensors',None) is not None else set()
    profile['neurons'].extend(motors)
    proof=apply_native_passive_profile(c,profile,neural_rate_multiplier=1.)
    cache_updates=[]
    if getattr(c,'sensors',None) is not None:
        c.sensors.tau[:]=c.params[c.sensors.indices,1]*.001
        cache_updates.append(dict(kind='actual704sensory_readout_tau',count=len(c.sensors.tau),cache_matches_parameters=bool(np.array_equal(c.sensors.tau,c.params[c.sensors.indices,1]*.001))))
    output=dict(machine='easteregg2',family=family,anchor=anchor,source=str(source),source_sha256=sha(source),source_neurons=len(c.names),source_edges=len(c.edges),profile=profile,consumption=proof,sensor_contracts=sensors,motor_contracts=motor_contracts,runtime_cache_updates=cache_updates,physics_steps=0,training_steps=0,training_released=False)
    output['all_ionic_tonic_and_non_tau_fields_exact']=np.array_equal(c.params[:,[0,2,3,4,5,6,7,8,9,10]],original_params[:,[0,2,3,4,5,6,7,8,9,10]]) and np.array_equal(c.tonic,original_tonic)
    output['physical_adapter_cache_exact']=all(np.array_equal(value,getattr(c,key)) for key,value in physical_cache.items())
    output['source_edge_parameters_exact']=np.array_equal(c.edges,original_edges)
    output['source_roles']=[];input_ids={int(row[0]) for row in c.inputs};motor_ids={int(row[0]) for row in c.outputs};conditioned={c.names.index(row['name']) for row in sensors}
    for i,name in enumerate(c.names):
        role='retained_CaSFA_spiking_burster' if hasattr(c,'encoded') and i in c.encoded else 'actual32cell_sensory_spiking_carrier' if i in spike_inputs else 'retained_nonlinear_Ca_burster' if c.params[i,3]!=0 else 'physical_motor_voltage_interface' if i in motor_ids else 'independent_fast_encoder' if i in conditioned else 'physical_input_and_recurrent_IN' if i in input_ids else 'retained_passive_conductance_transfer'
        output['source_roles'].append(dict(name=name,role=role,physical_tau_s=float(c.params[i,1]/1000),Size=float(c.params[i,2]),ionic_Table2_applicable=False if c.params[i,3]!=0 else None))
    scalar=[]
    for contract in sensors:
        i=c.names.index(contract['name']);p=c.params[i:i+1].copy()
        if i in spike_inputs:continue
        for fraction in (0.,.25,.5,1.):
            current=contract['input_current_nA'][1]*fraction;v=np.array([p[0,0]]);m=np.zeros(1);h=m.copy();newv=v.copy();gs=m.copy();gi=m.copy()
            K.neural_steps(v,m,h,newv,gs,gi,p,np.empty((0,6)),np.zeros(1),np.array([current]),.2,500)
            target=p[0,0]+current/p[0,2]*(1-math.exp(-.1/contract['physical_leak_tau_s']));scalar.append(abs(float(v[0])-target))
    output['actual_selected_channel_transient_error_mV']=max(scalar,default=0.)
    output['unconditioned_physical_inputs']=[c.names[i] for i in sorted(input_ids-conditioned)]
    output['unconditioned_input_classification']=[]
    for i in sorted(input_ids-conditioned):
        outgoing=c.edges[c.edges[:,0]==i];incoming=c.edges[c.edges[:,1]==i]
        reason='Recurrent input' if len(incoming) else 'Graph-unreachable retained sensor' if not len(outgoing) else 'Source bias/heterogeneous windows: offset function required'
        output['unconditioned_input_classification'].append(dict(name=c.names[i],incoming_edges=len(incoming),outgoing_edges=len(outgoing),tonic_nA=float(c.tonic[i]),ionic=bool(c.params[i,3]!=0),window_lower=sorted(set(outgoing[:,3])) if len(outgoing) else [],window_upper=sorted(set(outgoing[:,4])) if len(outgoing) else [],reason=reason))
    feedback,groups=native_feedback_budget(c);output['feedback_budget_profile']=feedback;output['feedback_budget_contracts']=groups
    feedback_proof=apply_native_passive_profile(c,feedback);output['feedback_consumption']=feedback_proof
    errors=[]
    for group in groups:
        post=c.names.index(group['target']);edges=c.edges[group['indices']].copy()
        for fraction in (0.,.25,.5,1.):
            p=c.params.copy();p[:,3]=0.;v=p[:,0].copy();v[post]=group['operating_mV'];prior=v.copy()
            for e in edges:v[int(e[0])]=e[3]+fraction*(e[4]-e[3])
            prior[:]=v;m=np.zeros(len(v));h=m.copy();newv=v.copy();gs=m.copy();gi=m.copy()
            K.neural_steps(v,m,h,newv,gs,gi,p,edges,np.zeros(len(v)),np.zeros(len(v)),.2,1)
            measured=p[post,2]*(gi[post]-gs[post]*(prior[post]-p[post,0]));errors.append(abs(float(measured)-group['new_reference_budget_nA']*fraction))
    output['actual_fanin_current_error_nA']=max(errors,default=0.)
    motor_currents={c.names.index(row['name']):max((point['neural_voltage_mV']-c.params[c.names.index(row['name']),0])*c.params[c.names.index(row['name']),2]-c.tonic[c.names.index(row['name'])] for point in row['voltage_active_domain']) for row in motor_contracts}
    conditional=native_conditional_input_checks(c,sorted(set(range(len(c.names)))-conditioned-spike_inputs),motor_currents);output['conditional_input_contracts']=conditional
    ionic_checks=native_ionic_mechanism_checks(c);output['retained_ionic_mechanism_checks']=ionic_checks
    output['graph_unreachable_inputs']=[row['name'] for row in output['unconditioned_input_classification'] if row['outgoing_edges']==0]
    if spike_inputs:output['actual_spiking_port_checks']=native_spiking_port_checks(c)
    adapters=native_adapter_consumption_checks(c);output['actual_physical_adapter_consumption']=adapters
    output['new_budget_scope']='Explicit new0.5nA source-equivalent current budget perpost/physical-kind/sign atrest+2.5mV; not an inherited source pass, voltage target, or tied bilateral Hunt fit. Full synaptic saturation bounds the contribution even for unbounded physical polynomials.'
    output['appropriate_scope']='Coordinate accounting, electrical roles, and input/output domains precede scoped training; free3Dwalking/wholebody success is not required. Source nonlinearCa mechanism is retained, not pureGLIF.'
    output['electrical_checks_passed']=output['all_ionic_tonic_and_non_tau_fields_exact'] and output['physical_adapter_cache_exact'] and max(scalar,default=0.)<1e-10 and max(errors,default=0.)<1e-10
    if spike_inputs:output['electrical_checks_passed']=output['electrical_checks_passed'] and output['actual_spiking_port_checks']['passed']
    output['electrical_checks_passed']=output['electrical_checks_passed'] and all(row['maximum_actual_native_step_error_mV']<1e-10 for row in conditional)
    output['electrical_checks_passed']=output['electrical_checks_passed'] and all(row['compiled_reference_error']<1e-10 and row['state_bounds_preserved'] and row['positive_physical_gate_times'] and row['m_increases_with_voltage'] and row['h_decreases_with_voltage'] for row in ionic_checks)
    output['electrical_checks_passed']=output['electrical_checks_passed'] and adapters['actual_cached_mapping_vs_source_tuple_error_nA']<1e-10
    output['helper_sha256']=sha(__file__)
    write_raw_record('NATIVE_ROLE_DESIGN_'+family+output_suffix+'.json',output)
    print(json.dumps(dict(family=family,neurons=len(c.names),sensory_contracts=len(sensors),motor_contracts=len(motor_contracts),actual_adapter=anchor)))
    if not output['electrical_checks_passed']:raise AssertionError('Actual source role/channel check failed; diagnostic retained')


def cad_network():
    sys.path.insert(0,str(ROBOT))
    import run_neural_robot as run
    original_g=run.params.G.copy();original_tau=run.params.TAU.copy()
    try:
        run.params.G.update(c1_gain=.15,v3_gain=0.,heel_rge=.15,toe_rge=.05)
        run.params.TAU['rg_nap_h']=.35
        actuators=[a.get('name') for a in ET.parse(ROBOT/'ah001_neural.xml').getroot().find('actuator')]
        net=run.bn.build(actuators,dt=.002,interleg=True)
    finally:
        run.params.G.clear();run.params.G.update(original_g);run.params.TAU.clear();run.params.TAU.update(original_tau)
    return net


def cad_role_design(output_suffix=''):
    """CAD262 electrical roles; allocator and pressure pair remain distinct."""
    net=cad_network();c=net.compiled;names=[p['name'] for p in net.net.populations]
    import run_neural_robot as R
    initial={key:value.copy() for key,value in c.__dict__.items() if isinstance(value,np.ndarray)}
    mask=np.any(c.g_ion!=0,axis=0);conditional=[]
    for post,name in enumerate(names):
        if mask[post]:continue
        errors=[]
        for fraction in (0.,.5,1.):
            c.V[:]=fraction*5.;c.V[post]=2.5;c.b_gate[:]=initial['b_gate'];c.c_gate[:]=initial['c_gate']
            old=c.V.copy();u=net.make_inputs();u[:]=fraction*5.
            activation=np.clip((old[None,:]-c.e_lo)/(c.e_hi-c.e_lo),0.,1.);g=c.g_max_non*activation;g[~np.isfinite(g)]=0.
            current=float(c.input_connectivity[post]@u);syn=float(np.sum(g[post]*(c.del_e[post]-old[post])))
            expected=old[post]+c.time_factor_membrane[post]*(-c.g_m[post]*(old[post]-c.V_rest[post])+c.i_b[post]+syn+current)
            c.forward(u);errors.append(abs(float(c.V[post])-expected))
        conditional.append(dict(name=name,role='MN engineering/physical output' if name.startswith('MN_') else 'retained graded conditional transfer',calibration_input_nA=[0.,2.5,5.],local_source_tau_s=float(c.c_m[post]/c.g_m[post]),maximum_equation_error_mV=max(errors)))
    for key,value in initial.items():
        target=getattr(c,key)
        if target.shape==value.shape:target[:]=value
    motor_ids=np.array([net.idx['MN_'+act] for act in net.muscles]);u_targets=np.linspace(0,1,len(motor_ids));error=0.
    for tick in range(25):
        u=net.make_inputs();old=c.V[motor_ids].copy();bias=R.motor_target_current(c,motor_ids,u_targets,tau=.005)
        for j,act in enumerate(net.muscles):u[net.inputs.index('POST_'+act)]=bias[j]
        expected=old+c.dt/.005*(5*u_targets-old);c.forward(u);error=max(error,float(np.max(abs(c.V[motor_ids]-expected))))
    result=dict(machine='easteregg2',source_neurons=len(names),source_passive_roles=len(conditional),retained_NaP_cells=[names[i] for i in np.flatnonzero(mask)],NaP_fixed_h_s=[float(c.tau_max_b[0,i]) for i in np.flatnonzero(mask)],NaP_pureGLIF_applicable=False,passive_contracts=conditional,
                external_hybrid_motor_contract=dict(channels=len(motor_ids),normalized_domain=[0.,1.],desired_voltage_mV='5*u',tracking_tau_s=.005,maximum_all54_error_mV=error,background_cancellation_explicit=True,not_neural_coordination=True),
                independent_pressure_pair_contract='Existing separate rational endpoint0.25mV/command0.05 and30ms discrete full-input pole;6literalSNSweights, no pair target clamp. Source pressuregeometry/delay unchanged. See preserved component proof.',
                mechanical_scope='Root/hips/lumbar assisted by measured apparatus/external hybrid control; fourbar constraints reduce six rawhinges to2effective knees. No3Dfreebalance orcleanstanding claim.',
                electrical_checks_passed=max(row['maximum_equation_error_mV'] for row in conditional)<1e-10 and error<1e-10,physics_steps=0,training_steps=0,whole_body_verified=False,helper_sha256=sha(__file__))
    write_raw_record('CAD262_ROLE_DESIGN'+output_suffix+'.json',result)
    print(json.dumps(dict(passive_roles=len(conditional),NaP_cells=int(mask.sum()),motor_channels=len(motor_ids),passive_error=max(row['maximum_equation_error_mV'] for row in conditional),motor_error=error,passed=result['electrical_checks_passed'])))
    if not result['electrical_checks_passed']:raise AssertionError('CAD electrical role check failed')


def base(family,adapter,sources,coverage,blockers):
    return dict(schema='aarl_non_gait_functional_contract_v1',family=family,machine='easteregg2',
                phase='design_only',source_evidence=[evidence(p) for p in sources],
                controlled_subspace=coverage,design_passed=False,scaling_passed=False,pathways_passed=False,
                training_allowed=False,curriculum_allowed=False,
                blockers=blockers,
                limitations=[],design_scope_policy='Account for every source coordinate as neural, passive, constraint, apparatus, or unsupported. Free3Dwalking/wholebody success is not a prerequisite for supported/reduced-scope Hunt training.',
                applicable_design_operations=dict(signal_scale='Declare actualsource and independent desired domains; no inherited globalpass.',passive_membrane='Derive localCfrom physicaltau and leak/backgroundconductance with actualclock.',graded_synapse='Exact single-path endpoint conductance derivation available; multiinput/transient/linearitychecks stillrequired.',ionic_burster='Retain sourceNaP/Ca/gate law; pureGLIFTable2notapplicable.',GLIF_threshold_steps='Notapplicableto gradedcells; separate pureGLIFprofile requiresall activation/transient/rate contracts.',spike_synapse='Declare reset-to-peak vs additive/capped eventlaw beforeusingdelta/Eq28/29.'),
                desired_functional_contract=dict(input_domain=None,output_domain=None,target_gain=None,target_transient_s=None,
                                                target_bandwidth_hz=None,target_rate_hz=None,status='Unassigned targets block application; inherited values are not independent desired functions.'),
                opt_in_profile=dict(adapter=adapter,adopted=False,neurons=[],synapses=[]),
                application_recipe=dict(default_behavior='No profile: literal no-op; no original source or topology changes.',
                                        production_hook_installed=False,
                                        unit_adapter_proof=str(RAW_RECORDS/'NON_GAIT_ADAPTER_ELECTRICAL_CHECK.json'),
                                        candidate_consumption_proof=str(RAW_RECORDS/'NON_GAIT_CANDIDATE_COMPONENT_CHECKS.json'),
                                        fresh_process_proof=str(RAW_RECORDS/'NON_GAIT_FRESH_PROCESS_CHECK.json'),
                                        Simulink_block_consumption_proof=str(RAW_RECORDS/'NON_GAIT_SIMULINK_BLOCK_CHECK.json') if adapter=='simulink_mask_recipe' else None,
                                        shared_GLIF_module=str(HERE/'table2_glif.py'),
                                        ionic_exclusion='NaP/Ca m/h and adaptation currents are not GLIF adaptive-threshold m/tau_theta; retain source mechanism unless a separately verified replacement is declared.'))


def generate():
    LEDGERS.mkdir(exist_ok=True)
    native_source=REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'
    rat_source=RAT/'Rat134_ReceptorGeometry_Standalone.asim'
    sources={}
    for key,path in [('W2L108_graded',native_source),('rat134_GA_receptor_variant',rat_source)]:
        c,path,actual_anchor=actual_native_controller(key)
        coverage=(dict(limb_mechanics='Three source limb hinge coordinates per side: hip/knee/ankle, six bilateral sagittal coordinates; not23DoFhuman.',global_support='WalkingPath/rails plus case-specific apparatus constrain base sway/orientation; match active rung before transfer.',uncontrolled_human_CAD_axes=['hip ab/adduction','hip axial rotation','lumbar three-axis motion','ankle/subtalar frontal rotation','MTP'],free_3D_balance_demonstrated=False)
                  if key.startswith('W2L') else dict(limb_mechanics='Six bilateral source hip/knee/ankle coordinates;16physicalmuscles. Receptor variant adds sensory neurons, not mechanical DOFs.',global_support='Source frozen pelvis/suspended preparation; actual ground/root repair is separate.',uncontrolled_human_CAD_axes=['hip frontal/transverse rotations','trunk3axes','subtalar andMTP'],free_3D_balance_demonstrated=False))
        ledger=base(key,'native_passive',[path,VIRTUAL/'native_kernel.py',VIRTUAL/'SZCZECINSKI2020_SEVEN_STEP_CONTROLLER_AUDIT.md'],coverage,
                    ['Remaining recurrent/passive pathway role and summed-input electrical contracts require verification.',
                     'Actual designated controller physical-adapter and runtime-cache consumption must be checked; raw XML alone is insufficient.'])
        ledger['limitations']=['This supported sagittal source does not establish free3Dhuman/CAD control. That is a transfer-scope limitation, not a prerequisite for its supported-subspace Hunt training.','Physical transduction and coordination can subsequently be calibrated in Hunt B/C within declared domains; old waveform passes are not new design proof.']
        if key.startswith('W2L'):
            imported=json.loads((EXTERNAL/'spiking/W2L_VirtualWalker_Newer/IMPORT_AUDIT.json').read_text())
            entries=[]
            for joint in imported['joints']:
                name=joint['name'];role='neural' if name.startswith(('hip_','knee_','ankle_')) else 'passive' if name.startswith('toe_') else 'apparatus' if name in ('Treadmill Slide','Passive guide X slider') else 'unsupported'
                entries.append(dict(name=name,kind=joint['type'],coverage=role))
            entries.append(dict(name='root6DoF',kind='floating_base',coverage='apparatus',note='Source WalkingPath/rails and selected case harness/guides account for base freedom; no inferred free3Dneural balance.'))
        else:
            entries=[dict(name=j.findtext('Name'),kind=j.findtext('Type'),coverage='neural' if any(k in (j.findtext('Name')or'').lower() for k in ('hip','knee','ankle')) else 'passive' if j.findtext('Name')=='Mid spine' and j.findtext('EnableMotor')=='False' else 'constraint' if j.findtext('Type') in ('Static','Fixed') else 'unsupported',source_lower_limit=j.findtext('LowerLimit/LimitPos'),source_upper_limit=j.findtext('UpperLimit/LimitPos')) for j in ET.parse(path).getroot().iter('Joint')]
            entries.append(dict(name='root6DoF',kind='source_frozen_pelvis',coverage='constraint',note='Frozen/suspended primary reference preparation; no freebalance claim.'))
        ledger['controlled_subspace']['coordinate_accounting']=entries
        ledger['controlled_subspace']['all_source_coordinates_accounted']=not any(x['coverage']=='unsupported' for x in entries)
        ledger['source_census']=dict(neurons=len(c.names),graded_synapses=len(c.edges),ionic_cells=int(np.count_nonzero(c.params[:,3])),neural_step_s=c.dt)
        ledger['actual_runtime_anchor']=actual_anchor
        ledger['units']=dict(voltage='Absolute mV',neural_runtime='seconds converted to ms by dt*1000',TimeConst='ms',current='source nA with explicit Size normalization',conductance='Actual native normalized SynAmp; do not claim measured physicaluS/Cm without normalization audit.')
        ledger['source_operating_ranges']=dict(rest_mV=[float(c.params[:,0].min()),float(c.params[:,0].max())],tau_ms=[float(c.params[:,1].min()),float(c.params[:,1].max())],Size=[float(c.params[:,2].min()),float(c.params[:,2].max())],synaptic_threshold_mV=[float(c.edges[:,3].min()),float(c.edges[:,3].max())],synaptic_saturation_mV=[float(c.edges[:,4].min()),float(c.edges[:,4].max())])
        ledger['source_ionic_names']=[n for n,p in zip(c.names,c.params) if p[3]!=0]
        proposed,contracts=native_input_contract(c,path)
        motor_overrides,motor_contracts=native_motor_contract(c,path)
        proposed['neurons'].extend(motor_overrides)
        ledger['opt_in_profile']=proposed
        ledger['proposed_functional_components']=contracts
        ledger['motor_output_contracts']=motor_contracts
        ledger['motor_clock_goal']='RetainactualsourceMNphysicalleakτ and actualSIvoltage→sigmoidmap while RG/PFcadence is separatelydefined. AdaptercompensatesonlydeclaredMNs byactualneural_rate_multiplier.'
        ledger['engineering_target_rationale']='Preserve eachactual pure-input sensorycurrent→voltage domain and DCgain fromsource Size/outgoingwindows; choose a10msphysicalpole (50samples at.2ms) toseparatefastsensoryconditioning fromsourcecadence. Relative to.292sfastestretainedratAcontext itisT/29.2, notmassscaling or a universalbiologicalbandwidthclaim.'
        ledger['desired_functional_contract'].update(input_domain='Perselectednode actualadmissiblecurrentrange below; out-of-domaincurrents requireexplicitreview, notsilentclamping.',output_domain='Actualsourcepresynapticwindows below.',target_gain='1/Size fordeclaredinput-onlycells',target_transient_s=.01,target_bandwidth_hz=1/(2*math.pi*.01),status='Explicitengineeringcandidate forselectedpassiveinputrelays. Remainingionic/sourcepathwaysretainedandscoped; notwholefamilyvalidated.')
        ledger['clock_contract']=clock_contract()
        ledger['application_recipe']['call']='apply_native_passive_profile(controller, ledger_or_explicit_profile) after constructing actual sourcecontroller and before preparing states; no runner hook installed.'
        if key.startswith('W2L'):
            ledger['existing_reference_configuration']=evidence(VIRTUAL/'GENTLE_HIP_0P6_FLEX_0P5.json')
            ledger['canonical_source_matches_existing_reference_hash']=sha(path)==json.loads((VIRTUAL/'GENTLE_HIP_0P6_FLEX_0P5.json').read_text())['expected_source_sha256']
            ledger['source_evidence'].append(evidence(native_source.parent/'CANONICAL_WORKING_PROFILE.json'))
            ledger['clock_contract'].update(canonical_authored_clock_scale=2,canonical_parameters_already_serialized=True,runtime_neural_time_scale_default=1.,do_not_multiply_authored_clock_twice=True,source_time_scaler='native_kernel.NativeController multipliesTimeConstandCaactivation/deactivationtimeconstants byneural_time_scale, includingpassiveafferents/MNs.',separation_recipe='Retain declaredRG/PFionic clocks; override onlyselectedplain afferent/MN tau_phys throughthisadapter usingactualruntime lambda. Do notscale sourceparamswholesale tohumanmass/cadence.')
            ledger['preserve_explicit_changes']=['Sourceclock and imported sourcehash','anklePF targetswap is routing, not a universal biological correction','sharedPF→MN/contact coefficients and groundedhome remain separate fromFSAderivation']
        else:
            ledger['source_evidence'] += [evidence(RAT/'scripts/rat134_receptor_runtime.py'),evidence(RAT/'scripts/rat134_core.py'),evidence(RAT/'prototype_role_policy_B/refinement12/EXPERIMENT_RECORD.md')]
            ledger['known_scope']='Curated134GA own-receptor engineering variant ofpublic120cellprecursor, notfinal2023replication. AcceptedintrinsicA and rejectedBheldouts remain unchanged.'
            ledger['audit_anchor_primary']=dict(runtime_constructor=evidence(RAT/'scripts/rat134_receptor_runtime.py'),constructor_function='variant_controller()',consumed_source=evidence(RAT/'Rat134_ReceptorGeometry.asim'),checkpoint=evidence(RAT/'FINAL_CHECKPOINT_PROOF.json'),condition_decisions=evidence(RAT/'CONDITION_DECISIONS.json'),secondary_rejected_B=evidence(RAT/'prototype_role_policy_B/refinement12/EXPERIMENT_RECORD.md'))
            actual=native_inventory(RAT/'Rat134_ReceptorGeometry.asim')
            ledger['standalone_vs_consumed_runtime_fields_exact']=dict(names=c.names==actual.names,neuron_parameters=np.array_equal(c.params,actual.params),edges=np.array_equal(c.edges,actual.edges),tonic=np.array_equal(c.tonic,actual.tonic))
            ledger['clock_contract'].update(source_time_scaler='Curated134runtime consumesactualsourceTimeConst/Cagates; replay_steps passesdt_ms andRatController.step passesdt_s*1000. Noadditionalglobalclockmultiplier inthese frozenfunctions.',separation_recipe='Preserve sourceintrinsicfieldfit separately; assign desiredphysical sensory/MN/actuatorresponse independently beforebody transfer.')
        sources[key]=ledger
    net=cad_network();core=net.compiled;names=sorted(net.idx,key=net.idx.get)
    coverage=dict(independent_full_model_DoF=17,coordinates='Floatingpelvis6+hips6+twofourbarcoordinates+lumbar3; physicallooprank4.',original_controller='11internalposecoordinates viaexternalhybridallocator; global6D support uses apparatus.',pressure_pair='Onlybilateralvas_med recruitment/pressure; no independenthipfrontal/transverse/lumbarbalance claim.',welded_foot_channels_without_joint_leverage=20,uncertified_axes=['lumbar3axiscontrolcapacity','hipfrontal/transversecapacity','CADphysicalROM','ankle/subtalar/MTParticulation'],free_3D_balance_demonstrated=False)
    cad=base('CAD262_hybrid','sns_passive',[ROBOT/'run_neural_robot.py',ROBOT/'ah001_neural.xml',SPINAL/'build_network.py',SPINAL/'params.py'],coverage,
             ['Mechanical stance/loop/contact/collision repair remains prerequisite; laptophandoff has fullstatus.',
              'Desired per-axis neural recruitment functions and actualsensor/actuator ranges are unassigned.',
              'Externalallocator/currenttargeting cannot be counted as FSAderived wholebodyneuralcontrol.'])
    cad['source_census']=dict(neurons=len(names),inputs=len(net.inputs),ionic_cells=int(np.count_nonzero(np.any(core.g_ion!=0,axis=0))) if core.gated else 0,active_pairs=int(np.count_nonzero(core.g_max_non)),actuators=54)
    cad['units']=dict(voltage='Restrelative0..5mV signal domain; ionicceiling is separate.',time='Compiled dt=.002seconds and clockcodedCm/Gm yields seconds; physicalnF equivalent is1000*Cmcode, not literalCmcode.',conductance='uS in sourceconductance equation',current='nA',actuator='FixedreferenceF0*filteredstate, not physiological5percentactivation.')
    cad['source_passive_examples']={name:dict(clock_capacitance=float(core.c_m[net.idx[name]]),leak_uS=float(core.g_m[net.idx[name]]),tau_s=float(core.c_m[net.idx[name]]/core.g_m[net.idx[name]]),bias_nA=float(core.i_b[net.idx[name]])) for name in ('MN_vas_med_r','PF_E1_r','Ib_vas_med_r') if name in net.idx}
    cad['application_recipe']['call']='apply_sns_passive_profile(net.compiled, ordered_names, explicit_profile, clock_unit="s") aftercompile/beforewarmup; preserve sparsecache whenusingwrapper.'
    cad['retained_engineering_target_contract']=dict(source_function='run_neural_robot.motor_target_current',desired_normalized_range=[0.,1.],desired_voltage='5*u mV',actual_consumed_tracking_tau_s=.005,actual_discrete_tracking_tau_s=float(-.002/math.log(1-.002/.005)),equation='Vnext=V+(dt/motor_tau)*(5*u−V) by explicitbackgroundcurrentcancellation.',downstream_actuator_state_filter_s=.03,not_pure_neural_coordination=True,note='Sourcecommentmentions30ms, but functiondefault and reviewconsumedCLI motor_tau=.005 define5msEulertracking; passiveMN-Cm30ms andactuator30msareseparate.')
    cad['audit_anchor_primary']=dict(manifest=evidence(EXTERNAL/'review_20261101/cad/REPRODUCTION_MANIFEST.json'),live_protocol=evidence(EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/support085/mixed/mixed/LIVE_PREREGISTRATION.json'),body_summary=evidence(EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/support085/mixed/mixed/cad.json'),scope='Mixed2pressure/52hybrid historical85case; mechanicalrepairfixtures arenotwinners.')
    sources['CAD262_hybrid']=cad
    pair=copy.deepcopy(cad);pair['family']='CAD262_MN_pressure_pair';pair['source_evidence'] += [evidence(ROBOT/'bpa_neural_coordination/cad_pair_unclamped.py'),evidence(EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/RECRUITMENT_PREREGISTRATION.json'),evidence(EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/support085/mixed/ASSUMED_PAIR_CONFIGURATION.json')]
    pair['known_source_function']='Six mirrored existing PF_E1/PF_E2/IBEXC→vas_medMN conductances normalized; pairPOSTzero, no pairoutputfloor. IBEXCdominantRGinput isadditivefeedforward, notmultiplicativeloadgate. Other52hybrid remains.'
    pair['preserved_output_contract']=dict(MN_voltage_range_mV=[0.,5.],normalized_command='u=clip(V_MN/5,0,1)',reference_tension='4530*u N for eachselectedroute',teaching_operating_target_u=.05,teaching_target_not_actual_forcefloor=True,pressure_interface='Sameinverseforce/series-equilibrium pressureconversion,620kPacap,40/60msfillvent,4msdelay. No pressurestate/force projection.',physical_actuator_clock='Remainsphysicaltime; notchangedbyneuralrateprofiles.',calibration_scope='Sixengineeringweightswere fitted tosourceoperatingpoint, not fullfunctionderived range/transientproof.')
    relative={'PF_E1':.2,'PF_E2':.3,'IBEXC_knee_ext':.6};target_v=.05*5.
    actual_reversals=[float(core.del_e[net.idx['MN_vas_med_'+side],net.idx[pre+'_'+side]]) for side in ('r','l') for pre in relative]
    assert len(set(actual_reversals))==1,'Separate reversal-specificfunctioncontract required'
    source_reversal=actual_reversals[0];g_total=endpoint_design(float(core.g_m[net.idx['MN_vas_med_r']]),5.,0.,target_v,source_reversal)['gmax_uS'];pair_edges=[]
    for side in ('r','l'):
        for pre,weight in relative.items():pair_edges.append(dict(source=pre+'_'+side,target='MN_vas_med_'+side,gmax_uS=g_total*weight/sum(relative.values())))
    pair['opt_in_profile']=dict(adapter='sns_passive',adopted=False,engineering_candidate=True,neurons=[dict(name='MN_vas_med_'+side,tau_s=.03,operating_synaptic_uS=g_total,response_policy='exact_discrete_pole') for side in ('r','l')],synapses=pair_edges)
    pair['proposed_functional_components']=dict(role='Isolatedexistingthree-inputfeedforwardrecruitmentbranch only; otherneuralcurrentszeroatthisoperatingpoint.',input='s=(.2*sPF1+.3*sPF2+.6*sIBEXC)/1.1, eachs=clip(Vpre/5,0,1)',desired_endpoint_MN_mV=target_v,desired_endpoint_command=.05,desired_fullinput_pole_s=.03,derived_total_gmax_uS=g_total,actual_compiled_reversal_mV=source_reversal,source_ratios_preserved=relative,exact_function='V_MN(s)=Eactual*Gtotal*s/(Gleak+Gtotal*s); endpointFSA, not exactlinear allrange.',force_reference_at_endpoint_N=4530*.05,not_force_floor=True,not_pressure_state_assignment=True,not_wholebody_gainclaim=True,rationale='Retainprevious5%nominalteachingpointand30msMNresponse asdeclaredlocalengineeringcontract; deriveendpointconductancesexactlyratherthanreusetunednormalization. PreservepressureconversionandliteralSNSrecruitment.')
    pair['blockers'] += ['Only2virtualpressureactuators; frozengeometry fails repairedstance domain. No mount/hardware/dynamic calibration.', 'Signed harness/support andcrossedlegacyfeet disallow cleanstanding/pairautonomous claims.']
    sources['CAD262_MN_pressure_pair']=pair
    hybrid=base('W2L704_hybrid_spiking','retained_CaLIF',[VIRTUAL/'spiking_central_runtime.py',VIRTUAL/'spiking_sensory.py',EXTERNAL/'spiking/COMBINED_SPIKING_CHECKPOINT.json',EXTERNAL/'spiking/CENTRAL_RESET_REST_CANDIDATE.json'],sources['W2L108_graded']['controlled_subspace'],
                ['Actual704-cell controller source-clock/cache and decoded-axon functional contracts require verification.',
                 'Summed sensory/decoded feedback current domains require explicit electrical checks.'])
    hybrid['limitations']=['Ca/SFA nonlinear burst mechanism is retained; pure adaptive-threshold GLIF Table2 is not applicable to it.','Saved55.376Hz is a decoder gain/normalizer, not automatically a biologicalmaximum firingrate.','Supported-subspace scope does not require prior free3Dwalking. Eq27/29 pureGLIF synapse recipes are not claimed for count-filter outputs.']
    hybrid['source_census']=dict(source_populations=108,central_spiking_populations=12,central_cells_each=16,sensory_spiking_populations=16,sensory_cells_each=32,actual_spiking_cells=704,remaining_graded_cells=80)
    hybrid['units']=dict(neural_membrane='Source native tau_ms with dt_s*1000',threshold='AbsolutemVfixed, notGLIFthetaadaptation',adaptation='Inhibitorycurrentjump/decay, not adaptive threshold',output='Decoded spike-countproxyvoltage, notsomaVm',sensor_events='NormalizedLIFreset0/threshold1/refractory1ms; populations encode II/contactinputs')
    hybrid['source_parameters']=json.loads((EXTERNAL/'spiking/CENTRAL_RESET_REST_CANDIDATE.json').read_text())
    hybrid['proposed_functional_components']=dict(role='Existinginput-onlyII/contactcarriers, not centralCa/SFA replacement.',normalized_input=[0.,1.],physical_requested_rate_hz=[0.,200.],physical_LIF_tau_s=.005,physical_refractory_s=.001,physical_readout_tau_s=.01,source_proxy_range_relative_mV=[0.,20.],clock_separation_helper='carrier_clock_parameters; pass derivedarguments intoexistingencode_step on an isolated/reusablecarrierport, notallneuralcells.',Table2_scope='Countfilter readout is not a reset-to-peak conductancesynapse; carrierclock/ratecheckdoesnotcompleteEq27/29.',mechanism_retained='Ca/SFA RG/PF andothergradedsourcecells unchanged.')
    hybrid['proposed_functional_components']['preferred_clock_separation']='PhysicalCarrierPort advances unchangedsourceencode_step on actualphysicaltime. Centralbiologicalclockcannotalteritscarrier/filter/refractory. Existing runtimehooknotinstalled; oldmatchingencoder mustbedisabled, notduplicated.'
    hybrid['clock_contract']=clock_contract()
    hybrid['clock_contract'].update(actual_source=evidence(VIRTUAL/'spiking_central_runtime.py','time/=self.clock_scale'),physical_substep='native_dt*clock_scale',all_affected=['gradedMN/IN/afferentmembranes','CaRG/PFkinetics','CaLIFSFAcurrentdecay','sensoryLIF5msmembraneand1msrefractory','sensoryfilterand200Hzcarrier','central50.886msfilterand55.376Hzreadout','startupstimuluswindows'],unaffected=['MuJoCo/Hillphysicalsteps','physicalafferentsampleintervalunlessexplicitlychanged'],source_neural_only_clock_is_one=True,separation_recipe='NewcontractmustseparateRG/PFtimingfromsensory/MN/actuatorphysicalresponse. Plaincelladaptercancompensateexistingpassivecells; ionic/sensorycarrier/readoutneedtheir ownverifiedpolicies. Do notapplyoneuniversalrat-tohumanclock.')
    hybrid['application_recipe']['status']='Retainionicmechanism; onlyseparatelydeclaredplainports mayusepassive/pureGLIFdesign. No newpureGLIFprofileapplied.'
    sources['W2L704_hybrid_spiking']=hybrid
    for spiking in (False,True):
        family='SNS_Simscape_'+('spiking' if spiking else 'graded')+'_reflexes'
        paths=[SIMULINK/'sns_build_library.m',SIMULINK/'demos/sns_build_demo.m']
        if spiking:paths += [SIMULINK/'sns_build_spiking.m',SIMULINK/'demos/sns_build_knee_spiking.m',SIMULINK/'demos/sns_build_rgmn_spiking.m']
        ledger=base(family,'simulink_mask_recipe',paths,
                    dict(source_control='KneeReflexDemo is1DoFknee antagonisticreflex; RG/MNsubnetwork demo is notfullwalker.',global_support='Analyticaljointorfixture, no free3Dwholebodybalance.',missing_human_CAD_control=['otherjointaxes','pelvis6DoF','lumbar3DoF','frontal/transversebalance'],free_3D_balance_demonstrated=False),
                    ['Perblockfunctionalranges/gains/transienttargets and actualmaskapplication proof remain required.',
                     'Remaining source block/functional-channel verification must be completed within the declared demonstration scope.'])
        ledger['limitations']=['Only the named reflex/mechanism coordinates are modeled; no3Dhuman/CAD coverage is inherited or required before scoped training.']
        ledger['units']=dict(clock='Simulinkseconds',voltage='AbsolutemV',current='nA',conductance='uS',capacitance='PhysicalnF; derivative1000/Cm',synapse_output='[g;g*Esyn], postsynapticcellcomputesIsyn')
        ledger['application_recipe']['example']=simulink_mask_recipe('SEPARATE_MODEL/SN_Ia_ext',tau_s=.03,leak_uS=1.)
        if spiking:
            ledger['known_source_parameters']=dict(rest_mV=-52,threshold_mV=-46,reset_mV=-58,leak_uS=.4,capacitance_nF=2.,membrane_tau_s=.005,synaptic_decay_s=.01)
            lif=dict(rest_mV=-52.,threshold_mV=-46.,reset_mV=-58.,leak_uS=.4,capacitance_nF=2.)
            ledger['source_preserving_rate_contract']=dict(scope='Analyticalinput-onlyLIF atsourcebaseline; noMATLABruntimeclaim orwholefamilyFSApass.',baseline_current_nA=[2.9,3.1,5.2,5.3],analytical_rates_hz=[source_reset_lif_rate_hz(current,**lif) for current in (2.9,3.1,5.2,5.3)],silent_input_upper_nA=2.4,role='Sensorycurrent→spikerate→additivehybridsynapse; lowcurrentdeadbandmustbeaccepted/redesignedexplicitly.',recipe='Retainactualresetandderiveexactcurrent↔rate law forchosenrange, ORexplicitlydeclare newreset-at-restTable2primitive. Do not silentlyreplaceonewiththeother.')
            ledger['desired_functional_contract'].update(input_domain=[2.9,5.3],output_domain='Sourcefixed-resetLIFspikes at analyticallyderivedrate, gradedMNunchanged.',target_rate_hz=[source_reset_lif_rate_hz(2.9,**lif),source_reset_lif_rate_hz(5.3,**lif)],target_transient_s=.005,status='Declaredsource-baselinepreservingengineeringportcontract; no whole-range/Simulinkruntimeproofyet.')
        else:
            ledger['desired_functional_contract'].update(input_domain='ExistingnAafferentport with ownmusclegeometrycalibration; exactphysicalmaxmustbesavedperdemo.',output_domain='Declaredpost-rest0..Rpresynapticsignal withsource[ g;gEsyn ]ports.',target_gain='Endpointderivationavailableforexistingexc/inhibpath only.',target_transient_s=.03,status='Proposed30msphysicalpassiveMN/relayexample; installonlyafteractualblockpath/domainverified.')
            ledger['blockers'] += ['Oldspikingsynapses ADDgincandcapgmax; Table2moduleRESETStopeak. Eq28/29cannotbeinstalledasgincwithoutdeclaringchangedprimitive.',
                                   'Source resetbelowrest andfixedthreshold not automaticallypureGLIFresetrest/equations; adaptiveθblock absent.',
                                   'Existing knee match was firstorder+handtuned; configuredcap/rate notfull2020proof.']
        sources[family]=ledger
    ode=base('reduced_Deng14','deng14_passive',[SPINAL/'deng_cpg_ode.py',SPINAL/'deng_cpg_numpy.py',SIMULINK/'demos/sns_build_deng_cpg.m'],
             dict(mechanics='Reduced neural demonstration only; RG/PF/IN andhipMNs, no completebodyloop.',source_neurons='2RGHC+4PFHC+6IN+2hipMN=14; knee/anklePFsharesonebranch; nofullMNarray.',uncontrolled_axes='All mechanicalaxes requireexplicitbodycoupling; no23DoFclaim.',free_3D_balance_demonstrated=False),
             ['KeepNaP6cells/sourceh-gatemechanism distinctfromGLIFadaptiveθ; fixedtau_h350ms is nottau_theta.',
              'ODE usesms internally(t_end=20000for20s); legacySNSNumPycompileDTsecondswithCM5 requiresclockaudit beforeequivalenceclaim.',
              'PythonODEandSimulinkbuilderPF→MNconductances differ (2.565/3.632 vs2.6563/3.6325); chooseversionandtargetexplicitly, no silentlymergeddefaults.',
              'Selectivepassivewrappernowavailable; wholemechanism/crossplatform/bodyvalidationisnotimplied.'])
    ode['units']=dict(ODE='ms,mV,nA,nF,uS; physicaltau=5nF/1uS=5ms',Simulink='secondswith1000/Cm; gate350ms explicit',legacy_numpy='Do not execute top-level20s script fordesignimport; clock/capacitance mismatch flagged.')
    ode['source_parameters']=dict(capacitance_nF=5,leak_uS=1,rest_HC_mV=-60,rest_MN_mV=-100,NaP_uS=1.5,NaP_reversal_mV=50,h_time_constant_s=.350)
    ode['audit_anchor_primary']=dict(source=evidence(SPINAL/'deng_cpg_ode.py'),function='simulate(fixed,t_end=20000ms,dt=.1ms,i_pulse10nA,t_pulse100ms,pulse_ms1)',reference=evidence(SIMULINK/'results/deng_cpg_ref.mat'),comparators_only=[str(SIMULINK/'demos/sns_build_deng_cpg.m'),str(SPINAL/'deng_cpg_numpy.py')],mechanical_walker_checkpoint=None)
    ds=deng14_source_module();ode['opt_in_profile']=dict(adapter='deng14_passive',adopted=False,engineering_candidate=True,neurons=[dict(name=name,tau_s=.01 if name.startswith('MN_') else .005) for i,name in enumerate(ds.NAMES) if i not in ds.NA],synapses=[])
    dc=[]
    for i,name in enumerate(ds.NAMES):
        if i in ds.NA:continue
        incoming=[row for row in ds.S if int(row[1])==i]
        assert len(incoming)==1,'Declareamultiinputcontractexplicitly'
        row=incoming[0];conductance=float(row[2]);baseline=float(ds.VREST[i]);reversal=float(row[3]);target=(ds.GM*baseline+conductance*reversal)/(ds.GM+conductance)
        derived=endpoint_design(ds.GM,float(row[5]-row[4]),baseline,target-baseline,reversal)
        dc.append(dict(name=name,source=ds.NAMES[int(row[0])],presynaptic_domain_mV=[float(row[4]),float(row[5])],desired_off_mV=baseline,desired_on_mV=target,preserved_DC_gmax_uS=conductance,reconstructed_endpoint_gmax_uS=derived['gmax_uS'],physical_leak_tau_s=.01 if name.startswith('MN_') else .005,rationale='RetainsourceIN5msleakpole; optionalMN10msreadoutfilteris1/35ofthe350mshgate, preservingDCandallupstreamionicclock. Engineeringcandidate, notauthorintent.'))
    ode['proposed_functional_components']=dc
    ode['application_recipe'].update(call='run_deng14_passive_profile(profile,duration_s=...,dt_s=.0001,pulse_time_s=.1,pulse_duration_s=.001,pulse_current_nA=10)',source_anchor='Originalsimulate(fixed,t_end=20000ms,dt=.1ms,i_pulse10nA,t_pulse100ms,pulse_ms1)',default_none_delegates_exact_original=True,selective_C_via_original_vectorized_RHS=True,source_graph_or_ionic_parameters_changed=False)
    sources['reduced_Deng14']=ode
    for key in ('W2L108_graded','W2L704_hybrid_spiking','rat134_GA_receptor_variant'):
        record=RAW_RECORDS/('NATIVE_ROLE_DESIGN_'+key+'_final.json');ledger=sources[key]
        ledger['appropriate_design_evidence']=evidence(record)
        if record.exists():
            report=json.loads(record.read_text());valid=bool(report['electrical_checks_passed']) and ledger['controlled_subspace']['all_source_coordinates_accounted']
            ledger['appropriate_source_role_prerequisites_verified']=valid
            ledger['verified_role_summary']=dict(logical_source_populations=report['source_neurons'],independent_encoder_contracts=len(report['sensor_contracts']),remaining_passive_conditional_contracts=len(report['conditional_input_contracts']),retained_ionic_contracts=len(report['retained_ionic_mechanism_checks']),first_entry_feedback_groups=len(report['feedback_budget_contracts']),actual_physical_adapter_API=report['actual_physical_adapter_consumption'])
            ledger['feedback_budget_profile']=report['feedback_budget_profile']
            ledger['source_roles']={row['name']:row['role'] for row in report['source_roles']}
            ledger['review_required']='Parent must review declared current budgets/source-side patterns and release global training explicitly; no body performance prerequisite.'
            if valid:ledger['blockers']=[]
        else:ledger['appropriate_source_role_prerequisites_verified']=False
    scope_terms=('Mechanical stance/loop','Only2virtualpressure','Signed harness/support','Sagittal supported source','All source-support/3D','KeepNaP6cells','ODE usesms internally','PythonODEandSimulinkbuilder','Selectivepassivewrapper','Oldspikingsynapses','Source resetbelowrest','Existing knee match')
    for key,ledger in sources.items():
        moved=[x for x in ledger['blockers'] if x.startswith(scope_terms)]
        ledger['blockers']=[x for x in ledger['blockers'] if x not in moved];ledger['limitations'].extend(moved)
        if key.startswith('CAD262'):
            ledger['controlled_subspace']['coordinate_accounting']=[dict(name='root6DoF',coverage='apparatus'),dict(name='hip_balls_6axes',coverage='apparatus',note='External hybrid allocator plusSNS signals, not independent3Dneural control'),dict(name='six_fourbar_hinges',coverage='constraint',note='Four closure constraints leave2effective knees; pairSNS and52hybridroutes contribute'),dict(name='lumbar_ball3axes',coverage='apparatus',note='External hybrid posture/rotation guidance')]
            ledger['controlled_subspace']['all_source_coordinates_accounted']=True
            electrical=RAW_RECORDS/'CAD262_ROLE_DESIGN_final.json';ledger['appropriate_design_evidence']=evidence(electrical)
            if electrical.exists() and json.loads(electrical.read_text())['electrical_checks_passed']:
                ledger['blockers']=[];ledger['appropriate_electrical_role_prerequisites_verified']=True
                ledger['review_required']='Parent review of hybrid/2pressure scope; mechanical defects remain separate body-trial gates, not electrical FSA failures.'
        elif key=='reduced_Deng14':
            ledger['controlled_subspace'].update(coordinate_accounting=[],all_source_coordinates_accounted=True)
            ledger['blockers']=[];ledger['review_required']='PrimaryfixedODE selected and passive/sourceclock proof available; no legacy comparator or body outcome required for this mechanism scope.'
            component=RAW_RECORDS/'NON_GAIT_CANDIDATE_COMPONENT_CHECKS.json'
            ledger['appropriate_electrical_role_prerequisites_verified']=bool(component.exists() and json.loads(component.read_text())['components']['reduced_Deng14_passive_profile']['passed'])
        elif key.startswith('SNS_Simscape'):
            ledger['controlled_subspace'].update(coordinate_accounting=[dict(name='knee_fixture',coverage='neural',note='1DoFantagonisticreflexdemo; no other mechanical coordinates instantiated')],all_source_coordinates_accounted=True)
            ledger['desired_functional_contract'].update(input_domain='Declared electrical source input:0..5.3nA; fixedLIF threshold/deadband and exact nonlinear f-I retained.',output_domain='Actual neuron Vm/rate plus chemical [g,gE]; mean signed currents atdeclaredVm fromeventrecipe.',target_gain='Exact source RC and fixed-resetLIF functions; independent +/-1nA event-current target at100Hz.',target_transient_s=.01,status='Actual source block/chemical-port contracts; no pureGLIF adaptive-threshold replacement claimed.')
            ledger['event_mask_recipes']=[simulink_event_mask_recipe('SEPARATE_MODEL/Exc',mean_current_nA=1.,post_mV=-50.,reversal_mV=0.),simulink_event_mask_recipe('SEPARATE_MODEL/Inh',mean_current_nA=-1.,post_mV=-50.,reversal_mV=-80.)]
            synaptic=RAW_RECORDS/'NON_GAIT_SYNAPSE_PORT_CHECK_endpoint_and_cap_v2.json';ledger['actual_synaptic_port_evidence']=evidence(synaptic)
            if synaptic.exists() and json.loads(synaptic.read_text())['passed']:
                ledger['blockers']=[];ledger['appropriate_electrical_block_prerequisites_verified']=True;ledger['review_required']='Parent review ofreflex/mechanism scope; no coupled wholebody or adaptation claim.'
            if key=='SNS_Simscape_graded_reflexes':
                graded=RAW_RECORDS/'NON_GAIT_SYNAPSE_PORT_CHECK_graded_v1.json';ledger['actual_graded_synaptic_port_evidence']=evidence(graded)
                ledger['appropriate_electrical_block_prerequisites_verified']=bool(graded.exists() and json.loads(graded.read_text())['passed'])
                ledger['limitations']=[note for note in ledger['limitations'] if not note.startswith(('Oldspikingsynapses','Source resetbelowrest'))]
            else:
                ledger['limitations'].append('Source fixed-reset LIF and additive/capped event synapses are tested as their exact functions; adaptive-threshold GLIF and reset-to-peak Eq.27/29 synapses are not claimed installed.')
        ledger['appropriate_local_design_checks_verified']=bool(ledger.get('appropriate_source_role_prerequisites_verified') or ledger.get('appropriate_electrical_role_prerequisites_verified') or ledger.get('appropriate_electrical_block_prerequisites_verified'))
        ledger['body_performance_is_not_design_prerequisite']=True
        if key=='rat134_GA_receptor_variant':
            ledger.update(audit_only=True,excluded_from_training_queue=True,appropriate_local_design_checks_verified=False,
                electrical_component_subchecks_passed=True,blockers=['Later unilateral diagnostic prototype, not the submitted public120rat.','Twelve Biar1/HE engineering cells bind left-only signals; MNphysicaloutputs disabled. Ia/II/Ib-named addedsensorports allactually readIa.','Wholecontroller bilateral/pathway adoption is notpassed; preserve as audit-only pending separate redesign.'],
                review_required='Do not adopt ortrain this prototype. Active requested ratdesign is originalpublic120 with explicit NEWbilateral-biarticular122candidate.')
    index=[]
    for key,ledger in sources.items():
        path=LEDGERS/(key+'.json')
        if path.exists():
            prior=json.loads(path.read_text())
            if 'bilateral_design_accounting' in prior:ledger['bilateral_design_accounting']=prior['bilateral_design_accounting']
        path.write_text(json.dumps(ledger,indent=2),encoding='utf-8')
        index.append(dict(family=key,path=str(path),design_passed=False,scaling_passed=False,pathways_passed=False,appropriate_local_design_checks_verified=ledger['appropriate_local_design_checks_verified'],global_review_pending=True,blocker_count=len(ledger['blockers'])))
    gate=json.loads((HERE/'CAMPAIGN_GATE.json').read_text())
    covered=['W2L108_graded','W2L704_hybrid_spiking','CAD262_hybrid_and_neural_pressure_pair','SNS_Simscape_graded_reflexes','SNS_Simscape_spiking_reflexes','reduced_Deng14']
    for family in ('public_rat120_precursor','reduced_Rybak_two_level'):
        path=LEDGERS/(family+'.json')
        if path.exists():
            ledger=json.loads(path.read_text());covered.append(family)
            index.append(dict(family=family,path=str(path),design_passed=False,scaling_passed=False,pathways_passed=False,appropriate_local_design_checks_verified=ledger.get('appropriate_local_design_checks_verified',False),global_review_pending=True,blocker_count=len(ledger.get('blockers',[]))))
    assert set(covered)<=set(gate['required_families'])
    (HERE/'NON_GAIT_DESIGN_INDEX.json').write_text(json.dumps(dict(schema='aarl_non_gait_design_index_v1',covered_global_families=covered,ledgers=index,source_inventory_only=False,optin_contracts_and_component_checks=True,training_allowed=False,all_family_design_verified=False,no_model_copies=True,raw_electrical_evidence_root=str(RAW_RECORDS),storage_move_manifest=str(RAW_RECORDS/'MOVE_MANIFEST.json'),historical_source_hashes_are_lineage=True),indent=2))
    print(json.dumps({key:s.get('source_census','contractrecipe') for key,s in sources.items()},indent=2))


def electrical_probes(output_suffix=''):
    """Isolated electrical fixtures only; no whole-body or training runs."""
    sys.path.insert(0,str(VIRTUAL))
    net=cad_network()
    from exact_fast_eval import ExactFastForward
    import native_kernel as K
    names=sorted(net.idx,key=net.idx.get);wrapped=ExactFastForward(net.compiled)
    before={key:value.copy() for key,value in wrapped.c.__dict__.items() if isinstance(value,np.ndarray)}
    noop=apply_sns_passive_profile(wrapped,names,None)
    zero_exact=all(np.array_equal(value,getattr(wrapped.c,key)) for key,value in before.items())
    configuration=json.loads((EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/RECRUITMENT_PREREGISTRATION.json').read_text())
    weights=configuration['absolute_shared_weights_uS']
    changes=[dict(source=pre+'_'+side,target='MN_vas_med_'+side,gmax_uS=value) for side in ('r','l') for pre,value in weights.items()]
    profile=dict(adapter='sns_passive',neurons=[dict(name='MN_vas_med_r',tau_s=.03)],synapses=changes)
    proof=apply_sns_passive_profile(wrapped,names,profile,clock_unit='s')
    core=wrapped.c;core.V[:]=core.V_rest;pre=net.idx['PF_E1_r'];post=net.idx['MN_vas_med_r'];core.V[pre]=5.;start=float(core.V[post])
    expected=core.dt/core.c_m[post]*(weights['PF_E1']*core.del_e[post,pre]*1.)
    wrapped.forward(np.zeros(len(net.inputs)))
    measured=float(core.V[post]-start)
    ionic=int(np.flatnonzero(np.any(core.g_ion!=0,axis=0))[0]);old=core.c_m.copy();rejected=False
    try:apply_sns_passive_profile(wrapped,names,dict(adapter='sns_passive',neurons=[dict(name=names[ionic],tau_s=.1)],synapses=[]))
    except ValueError:rejected=True
    native=native_inventory(REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim')
    i=int(np.flatnonzero(native.params[:,3]==0)[0]);p=native.params[i].copy();plain=SimpleNamespace(names=[native.names[i]],params=p[None,:].copy(),edges=np.empty((0,6)),tonic=np.zeros(1))
    native_old=plain.params.copy();native_noop=apply_native_passive_profile(plain,None)
    native_zero=np.array_equal(native_old,plain.params)
    tau_s=float(p[1]/2000.)
    diagnostic_rate=5.
    native_proof=apply_native_passive_profile(plain,dict(adapter='native_passive',neurons=[dict(name=plain.names[0],tau_s=tau_s)],synapses=[]),neural_rate_multiplier=diagnostic_rate)
    v=np.array([p[0]]);m=np.zeros(1);h=np.zeros(1);new=np.zeros(1);gs=np.zeros(1);gi=np.zeros(1);applied=np.ones(1);dt_ms=.2
    K.neural_steps(v,m,h,new,gs,gi,plain.params,plain.edges,plain.tonic,applied,dt_ms*diagnostic_rate,1)
    native_expected=float(p[0]+(1-math.exp(-dt_ms/(1000*tau_s)))/p[2])
    native_error=abs(float(v[0])-native_expected)
    lif=dict(rest_mV=-52.,threshold_mV=-46.,reset_mV=-58.,leak_uS=.4,capacitance_nF=2.)
    roundtrips=[]
    for rate in (10.,50.,100.,200.):
        current=source_reset_lif_current_nA(rate,**lif);recovered=source_reset_lif_rate_hz(current,**lif)
        roundtrips.append(dict(desired_Hz=rate,derived_current_nA=current,recovered_Hz=recovered,error_Hz=abs(rate-recovered)))
    # The10Hz inverse lies only9.9e-9nA above the2.4nA threshold: report
    # absolute error and use an explicit relative floating-point criterion.
    lif_relative=max(row['error_Hz']/row['desired_Hz'] for row in roundtrips)
    unit_pass=zero_exact and rejected and np.array_equal(old,core.c_m) and abs(expected-measured)<1e-12 and native_zero and native_error<1e-12 and lif_relative<1e-9
    result=dict(scope='Electrical adapter/sourceconsumption checks only; no wholemodel or training gate promotion.',passed=bool(unit_pass),
                sns_default_noop=noop,sns_default_arrays_exact=zero_exact,sns_consumption=proof,
                source_pair_pathway=dict(source='PF_E1_r',target='MN_vas_med_r',single_step_expected_mV=float(expected),single_step_measured_mV=measured,absolute_error_mV=abs(expected-measured),other_presynaptic_cells_at_rest=True),
                ionic_override_rejected=rejected,ionic_arrays_unchanged=np.array_equal(old,core.c_m),
                native_default_noop=native_noop,native_default_exact=native_zero,native_consumption=native_proof,
                native_singlecell=dict(source_name=plain.names[0],diagnostic_target_tau_s=tau_s,not_adopted_family_target=True,diagnostic_neural_rate_multiplier=diagnostic_rate,physical_dt_s=dt_ms/1000,biological_dt_s=diagnostic_rate*dt_ms/1000,expected_V_mV=native_expected,measured_V_mV=float(v[0]),error_mV=native_error),
                simulink_formula_example=simulink_mask_recipe('UNRUN_UNIT_FIXTURE',tau_s=.03,leak_uS=1.),
                retained_Simulink_reset_LIF_analytical_inverse=roundtrips,
                analytical_inverse_relative_tolerance=1e-9,near_threshold_10Hz_absolute_precision_diagnostic=True,
                MATLAB_runtime_executed=False,physics_steps=0,training_allowed=False,
                sources=[evidence(Path(__file__)),evidence(VIRTUAL/'native_kernel.py')])
    write_raw_record('NON_GAIT_ADAPTER_ELECTRICAL_CHECK'+output_suffix+'.json',result)
    print(json.dumps(dict(passed=bool(unit_pass),sns_error=abs(expected-measured),native_error=native_error,physics_steps=0),indent=2))
    if not unit_pass:raise AssertionError('Adapter electrical check failed; do not promote')


def candidate_component_probes(output_suffix=''):
    """Consume actual candidate profiles on isolated electrical components."""
    sys.path.insert(0,str(VIRTUAL));import native_kernel as K
    results={}
    for family,source in [('W2L108_graded',REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'),('rat134_GA_receptor_variant',RAT/'Rat134_ReceptorGeometry_Standalone.asim')]:
        c=native_inventory(source);old_params=c.params.copy();old_edges=c.edges.copy()
        ledger=json.loads((LEDGERS/(family+'.json')).read_text());profile=ledger['opt_in_profile']
        proof=apply_native_passive_profile(c,profile,neural_rate_multiplier=5.)
        checks=[]
        for contract in ledger['proposed_functional_components']:
            i=c.names.index(contract['name']);p=c.params[i:i+1].copy();duration=.05;dt_physical=.0002;rate=5.;steps=round(duration/dt_physical)
            for fraction in (0.,.5,1.):
                current=contract['input_current_nA'][1]*fraction;v=np.array([p[0,0]]);m=np.zeros(1);h=np.zeros(1);new=np.zeros(1);gs=np.zeros(1);gi=np.zeros(1)
                K.neural_steps(v,m,h,new,gs,gi,p,np.empty((0,6)),np.zeros(1),np.array([current]),dt_physical*1000*rate,steps)
                expected=p[0,0]+current/p[0,2]*(1-math.exp(-duration/contract['physical_leak_tau_s']))
                checks.append(dict(name=contract['name'],normalized_input=fraction,expected_mV=float(expected),measured_mV=float(v[0]),error_mV=abs(float(v[0])-expected)))
        motor_checks=[]
        for contract in ledger['motor_output_contracts']:
            i=c.names.index(contract['name']);p=c.params[i:i+1].copy();tau=contract['physical_leak_tau_s'];physical_dt=.0002;rate=5.;steps=math.ceil(5*tau/physical_dt);duration=steps*physical_dt
            for point in contract['voltage_active_domain']:
                target=point['neural_voltage_mV'];current=p[0,2]*(target-p[0,0])-c.tonic[i]
                v=np.array([p[0,0]]);m=np.zeros(1);h=np.zeros(1);new=np.zeros(1);gs=np.zeros(1);gi=np.zeros(1)
                K.neural_steps(v,m,h,new,gs,gi,p,np.empty((0,6)),np.array([c.tonic[i]]),np.array([current]),physical_dt*1000*rate,steps)
                expected=p[0,0]+(target-p[0,0])*(1-math.exp(-duration/tau))
                motor_checks.append(dict(name=contract['name'],actuator=contract['actuator'],active_plateau_fraction=point['active_plateau_fraction'],source_motor_domain_admissible=point['within_source_adapter_limits'],physical_pole_s=tau,expected_transient_mV=float(expected),measured_transient_mV=float(v[0]),error_mV=abs(float(v[0])-expected),source_sigmoid_inverse_error=point['inverse_roundtrip_error']))
        untouched=np.array_equal(c.params[:,[0,2,3,4,5,6,7,8,9,10]],old_params[:,[0,2,3,4,5,6,7,8,9,10]]) and np.array_equal(c.edges,old_edges)
        passed=bool(checks) and max(x['error_mV'] for x in checks)<1e-10 and untouched and all(row['error_mV']<1e-10 and row['source_sigmoid_inverse_error']<1e-12 for row in motor_checks)
        results[family]=dict(passed=passed,scope='Selectedactualpure-inputpassiverelays plus sourceMNphysicaltiming/active-sigmoiddomain; isolatedcells notwholebody.',contract_count=len(ledger['proposed_functional_components']),motor_contract_count=len(ledger['motor_output_contracts']),source_Size_ionic_gates_and_edges_exact=untouched,consumption=proof,checks=checks,motor_checks=motor_checks,whole_network_pass=False)
    hybrid_net=cad_network();hybrid_core=hybrid_net.compiled
    import run_neural_robot as robot_run
    mn=np.array([hybrid_net.idx['MN_vas_med_'+side] for side in ('r','l')]);ports={name:i for i,name in enumerate(hybrid_net.inputs)};desired=np.array([.05,.1]);expected=np.zeros(2);errors=[];actual=[]
    for tick in range(25):
        hybrid_core.V[hybrid_net.idx['PF_E1_r']]=5 if tick%2==0 else 2.5
        current=np.zeros(len(hybrid_net.inputs));post_current=robot_run.motor_target_current(hybrid_core,mn,desired,tau=.005)
        for j,side in enumerate(('r','l')):current[ports['POST_vas_med_'+side]]=post_current[j]
        expected=expected+hybrid_core.dt/.005*(5*desired-expected)
        hybrid_core.forward(current);errors.append(float(np.max(abs(hybrid_core.V[mn]-expected))));actual.append(hybrid_core.V[mn].tolist())
    results['CAD262_hybrid_target_interface']=dict(passed=max(errors)<1e-12,scope='Originalexplicitengineeringcurrent-target/cancellationinterface, NOTneuralcoordination; arbitrarychangingPFbackground usedonlyinelectricalfixture.',source_function='run_neural_robot.motor_target_current',desired_normalized=desired.tolist(),actual_tracking_argument_s=.005,actual_discrete_pole_s=float(-.002/math.log(1-.002/.005)),maximum_step_equation_error_mV=max(errors),final_MN_mV=actual[-1],physical_actuator_filter_separate_s=.03,source_comment_30ms_is_not_consumed_tracking_tau=True,body_steps=0)
    net=cad_network()
    from exact_fast_eval import ExactFastForward
    names=sorted(net.idx,key=net.idx.get);c=net.compiled;wrapper=ExactFastForward(c)
    pair_ledger=json.loads((LEDGERS/'CAD262_MN_pressure_pair.json').read_text());profile=pair_ledger['opt_in_profile']
    proof=apply_sns_passive_profile(wrapper,names,profile)
    posts=np.array([net.idx['MN_vas_med_'+side] for side in ('r','l')]);held=np.ones(len(names),dtype=bool);held[posts]=False
    pres=np.array([net.idx[change['source']] for change in profile['synapses']]);gt=pair_ledger['proposed_functional_components']['derived_total_gmax_uS'];erev=pair_ledger['proposed_functional_components']['actual_compiled_reversal_mV'];curve=[];first_full=None
    for activation in (0.,.25,.5,.75,1.):
        c.V[:]=c.V_rest
        for tick in range(500):
            c.V[held]=c.V_rest[held];c.V[pres]=5*activation
            wrapper.forward(np.zeros(len(net.inputs)))
            if activation==1. and tick==0:first_full=float(c.V[posts[0]])
        expected=erev*gt*activation/(1+gt*activation)
        curve.append(dict(normalized_weighted_input=activation,expected_MN_mV=expected,measured_MN_mV=c.V[posts].tolist(),maximum_error_mV=float(np.max(abs(c.V[posts]-expected))),command=(c.V[posts]/5).tolist()))
    measured_tau=-c.dt/math.log1p(-first_full/.25)
    s=np.linspace(0,1,1001);max_nonlinearity=float(np.max(abs(erev*gt*s/(1+gt*s)-.25*s))/.25)
    passed=max(row['maximum_error_mV'] for row in curve)<1e-10 and abs(measured_tau-.03)<1e-12
    results['CAD262_MN_pressure_pair']=dict(passed=bool(passed),scope='Actualsource3inputbranch/twoMNvoltageclampisolation; exact rationalendpointfunction+30mspole, NOTlinearallrange. Otherpresynapticcellsheldatrest; notfree-runningwholeCAD.',consumption=proof,curve=curve,actual_source_reversal_mV=erev,physical_fullinput_discrete_pole_s=measured_tau,endpoint_command=.05,endpoint_reference_tension_N=226.5,maximum_graded_linearity_error_fullscale=max_nonlinearity,initial_extra_linearity_diagnostic_limit=.005,extra_linearity_diagnostic_passed=bool(max_nonlinearity<=.005),linearity_scope='Initial20mVassumptionfailedandpreserved; actual8mVrationalfunctionisexplicit. Extra0.5%linearityaspiration remainsmissed; notsilentlylabelledlinear.',pressure_conversion_modified=False,force_floor_applied=False,body_steps=0)
    carrier_rows=[]
    for target in (0.,.25,.5,1.):
        a=PhysicalCarrierPort(['source_input_fixture']);b=PhysicalCarrierPort(['source_input_fixture']);a_bins=[];b_bins=[]
        for outer in range(500):
            a.advance(.002,[target]);before=b.total.copy()
            # Five faster biological caller intervals cover the same2ms physicaltime.
            for _ in range(5):b.advance(.002/5,[target])
            a_bins.append(a.counts.copy());b_bins.append(b.total-before)
        state_a=a.export_state();state_b=b.export_state();causal_keys=['v','refractory','decoded','last_events','total_counts','pending_physical_s','physical_steps','rng_state']
        same={key:np.array_equal(state_a[key],state_b[key]) for key in causal_keys};rate_hz=float(a.total[0]/a.N);requested=200*target
        checkpoint=a.export_state();a.advance(.02,[target]);first=a.export_state();a.import_state(checkpoint);a.advance(.02,[target]);second=a.export_state()
        continuation=all(np.array_equal(value,second[key]) for key,value in first.items())
        error=abs(rate_hz-requested)/requested if requested else float(a.total[0]!=0)
        carrier_rows.append(dict(normalized_input=target,requested_physical_Hz=requested,measured_population_Hz=rate_hz,relative_rate_error=error,relative_rate_limit=.10,cross_clock_causal_states=same,common_physical_bin_counts_exact=np.array_equal(a_bins,b_bins),all_state_continuation_exact=continuation,spikes_observed=int(state_a['total_counts'][0]),passed=all(same.values()) and np.array_equal(a_bins,b_bins) and continuation and error<=.10))
    results['W2L704_physical_carrier_port']=dict(passed=all(row['passed'] for row in carrier_rows),checks=carrier_rows,scope='Reusedsource32cellencoder onindependentphysicalclock; sourcecentralCa/SFAandMNunchanged/notstepped. Oldencoderhooknotreplacedinproduction.',event_semantics='Normalizedspikecountfilter, notTable2conductanceresetprimitive',production_hook_installed=False)
    ds=deng14_source_module();args=dict(duration_s=.2,dt_s=.0001,pulse_time_s=.1,pulse_duration_s=.001,pulse_current_nA=10.)
    _,baseline,bh,no_profile=run_deng14_passive_profile(None,**args)
    _,empty,eh,empty_proof=run_deng14_passive_profile(dict(adapter='deng14_passive',neurons=[],synapses=[]),**args)
    ledger=json.loads((LEDGERS/'reduced_Deng14.json').read_text());_,candidate,ch,proof=run_deng14_passive_profile(ledger['opt_in_profile'],**args)
    dc_error=max(abs(row['preserved_DC_gmax_uS']-row['reconstructed_endpoint_gmax_uS']) for row in ledger['proposed_functional_components'])
    difference=float(np.max(abs(candidate[:,-2:]-baseline[:,-2:])))
    default_exact=np.array_equal(baseline,empty) and np.array_equal(bh,eh);upstream=np.array_equal(baseline[:,:12],candidate[:,:12]) and np.array_equal(bh,ch)
    results['reduced_Deng14_passive_profile']=dict(passed=bool(default_exact and upstream and difference>0 and dc_error<1e-12),scope='OriginalfixedODE code,200mselectricalpulse/readoutcheckonly; no20smechanismorbodytrial.',default_and_empty_profile_exact=default_exact,all12upstream_cells_and_h_state_exact=upstream,MN_readout_difference_mV=difference,DC_endpoint_reconstruction_error_uS=dc_error,consumption=proof,source_pulse=dict(current_nA=10,time_s=.1,duration_s=.001),NaP6_C5nF_and_h350ms_retained=True,physical_clock_seconds_adapter_explicit=True)
    output=dict(schema='aarl_non_gait_candidate_component_checks_v1',passed=all(row['passed'] for row in results.values()),components=results,physics_steps=0,training_steps=0,source=evidence(Path(__file__)),global_gate_promoted=False)
    write_raw_record('NON_GAIT_CANDIDATE_COMPONENT_CHECKS'+output_suffix+'.json',output)
    print(json.dumps({key:dict(passed=row['passed'],contracts=row.get('contract_count'),scope=row['scope']) for key,row in results.items()},indent=2))
    if not output['passed']:raise AssertionError('Candidatecomponentfailed;preservedwithoutgatepromotion')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');parser.add_argument('--electrical-probes',action='store_true');parser.add_argument('--candidate-probes',action='store_true');parser.add_argument('--verify-storage',action='store_true');parser.add_argument('--cad-role-design',action='store_true');parser.add_argument('--native-role-design',choices=('W2L108_graded','W2L704_hybrid_spiking','rat134_GA_receptor_variant','public_rat120_precursor'));parser.add_argument('--output-suffix',default='');args=parser.parse_args()
    if any(c not in '_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789' for c in args.output_suffix):raise ValueError('Simple record suffix only')
    if args.generate:generate()
    if args.electrical_probes:electrical_probes(args.output_suffix)
    if args.candidate_probes:candidate_component_probes(args.output_suffix)
    if args.verify_storage:verify_storage(args.output_suffix)
    if args.native_role_design:native_role_design(args.native_role_design,args.output_suffix)
    if args.cad_role_design:cad_role_design(args.output_suffix)
