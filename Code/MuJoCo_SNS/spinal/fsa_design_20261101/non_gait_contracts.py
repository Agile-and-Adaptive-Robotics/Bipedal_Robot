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
        if not np.allclose(lo,rest,atol=1e-12) or not np.allclose(hi,hi[0],atol=1e-12) or hi[0]<=rest:continue
        span=float(hi[0]-rest);size=float(controller.params[i,2]);name=controller.names[i]
        if controller.names.count(name)!=1:continue
        contracts.append(dict(name=name,actual_physical_kinds=sorted(set(kinds)),input_current_nA=[0.,size*span],output_voltage_mV=[float(rest),float(hi[0])],DC_gain_mV_per_source_nA=1/size,physical_leak_tau_s=tau_physical_s,physical_leak_bandwidth_Hz=1/(2*math.pi*tau_physical_s),source_tau_ms=float(controller.params[i,1]),incoming_central_edges=0,physical_adapter_maps_unchanged=True))
        overrides.append(dict(name=name,tau_s=tau_physical_s))
    return dict(adapter='native_passive',adopted=False,engineering_candidate=True,neurons=overrides,synapses=[]),contracts


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


def simulink_mask_recipe(block,*,tau_s,leak_uS):
    """Physical-unit passive block patch only; this does not execute MATLAB."""
    d=passive_design(tau_s,leak_uS)
    return dict(block=str(block),mask_values=dict(Gm=leak_uS,Cm=d['capacitance_nF']),
                source_derivative='1000/Cm in seconds-based Simulink',runtime_consumption_verified=False,
                warning='Apply to a separate passive block/model via set_param; no whole-model or ionic rescaling.')


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


def base(family,adapter,sources,coverage,blockers):
    return dict(schema='aarl_non_gait_functional_contract_v1',family=family,machine='easteregg2',
                phase='design_only',source_evidence=[evidence(p) for p in sources],
                controlled_subspace=coverage,design_passed=False,scaling_passed=False,pathways_passed=False,
                training_allowed=False,curriculum_allowed=False,
                blockers=blockers,
                applicable_design_operations=dict(signal_scale='Declare actualsource and independent desired domains; no inherited globalpass.',passive_membrane='Derive localCfrom physicaltau and leak/backgroundconductance with actualclock.',graded_synapse='Exact single-path endpoint conductance derivation available; multiinput/transient/linearitychecks stillrequired.',ionic_burster='Retain sourceNaP/Ca/gate law; pureGLIFTable2notapplicable.',GLIF_threshold_steps='Notapplicableto gradedcells; separate pureGLIFprofile requiresall activation/transient/rate contracts.',spike_synapse='Declare reset-to-peak vs additive/capped eventlaw beforeusingdelta/Eq28/29.'),
                desired_functional_contract=dict(input_domain=None,output_domain=None,target_gain=None,target_transient_s=None,
                                                target_bandwidth_hz=None,target_rate_hz=None,status='Unassigned targets block application; inherited values are not independent desired functions.'),
                opt_in_profile=dict(adapter=adapter,adopted=False,neurons=[],synapses=[]),
                application_recipe=dict(default_behavior='No profile: literal no-op; no original source or topology changes.',
                                        production_hook_installed=False,
                                        unit_adapter_proof='NON_GAIT_ADAPTER_ELECTRICAL_CHECK.json',
                                        shared_GLIF_module=str(HERE/'table2_glif.py'),
                                        ionic_exclusion='NaP/Ca m/h and adaptation currents are not GLIF adaptive-threshold m/tau_theta; retain source mechanism unless a separately verified replacement is declared.'))


def generate():
    LEDGERS.mkdir(exist_ok=True)
    native_source=REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'
    rat_source=RAT/'Rat134_ReceptorGeometry_Standalone.asim'
    sources={}
    for key,path in [('W2L108_graded',native_source),('rat134_GA_receptor_variant',rat_source)]:
        c=native_inventory(path)
        coverage=(dict(limb_mechanics='Three source limb hinge coordinates per side: hip/knee/ankle, six bilateral sagittal coordinates; not23DoFhuman.',global_support='WalkingPath/rails plus case-specific apparatus constrain base sway/orientation; match active rung before transfer.',uncontrolled_human_CAD_axes=['hip ab/adduction','hip axial rotation','lumbar three-axis motion','ankle/subtalar frontal rotation','MTP'],free_3D_balance_demonstrated=False)
                  if key.startswith('W2L') else dict(limb_mechanics='Six bilateral source hip/knee/ankle coordinates;16physicalmuscles. Receptor variant adds sensory neurons, not mechanical DOFs.',global_support='Source frozen pelvis/suspended preparation; actual ground/root repair is separate.',uncontrolled_human_CAD_axes=['hip frontal/transverse rotations','trunk3axes','subtalar andMTP'],free_3D_balance_demonstrated=False))
        ledger=base(key,'native_passive',[path,VIRTUAL/'native_kernel.py',VIRTUAL/'SZCZECINSKI2020_SEVEN_STEP_CONTROLLER_AUDIT.md'],coverage,
                    ['Desired component-wise input/output ranges, gains, and transient/bandwidth targets are not yet independently assigned.',
                     'Whole-network pathway signs/net effects and sensory geometric calibration must be checked; retained old waveform passes do not satisfy this gate.',
                     'Sagittal supported source cannot establish frontal/transverse/lumbar coverage in unconstrained human/CAD.'])
        ledger['source_census']=dict(neurons=len(c.names),graded_synapses=len(c.edges),ionic_cells=int(np.count_nonzero(c.params[:,3])),neural_step_s=c.dt)
        ledger['units']=dict(voltage='Absolute mV',neural_runtime='seconds converted to ms by dt*1000',TimeConst='ms',current='source nA with explicit Size normalization',conductance='Actual native normalized SynAmp; do not claim measured physicaluS/Cm without normalization audit.')
        ledger['source_operating_ranges']=dict(rest_mV=[float(c.params[:,0].min()),float(c.params[:,0].max())],tau_ms=[float(c.params[:,1].min()),float(c.params[:,1].max())],Size=[float(c.params[:,2].min()),float(c.params[:,2].max())],synaptic_threshold_mV=[float(c.edges[:,3].min()),float(c.edges[:,3].max())],synaptic_saturation_mV=[float(c.edges[:,4].min()),float(c.edges[:,4].max())])
        ledger['source_ionic_names']=[n for n,p in zip(c.names,c.params) if p[3]!=0]
        proposed,contracts=native_input_contract(c,path)
        ledger['opt_in_profile']=proposed
        ledger['proposed_functional_components']=contracts
        ledger['engineering_target_rationale']='Preserve eachactual pure-input sensorycurrent→voltage domain and DCgain fromsource Size/outgoingwindows; choose a10msphysicalpole (50samples at.2ms) toseparatefastsensoryconditioning fromsourcecadence. Relative to.292sfastestretainedratAcontext itisT/29.2, notmassscaling or a universalbiologicalbandwidthclaim.'
        ledger['desired_functional_contract'].update(input_domain='Perselectednode actualadmissiblecurrentrange below; out-of-domaincurrents requireexplicitreview, notsilentclamping.',output_domain='Actualsourcepresynapticwindows below.',target_gain='1/Size fordeclaredinput-onlycells',target_transient_s=.01,target_bandwidth_hz=1/(2*math.pi*.01),status='Explicitengineeringcandidate forselectedpassiveinputrelays. Remainingionic/sourcepathwaysretainedandscoped; notwholefamilyvalidated.')
        ledger['clock_contract']=clock_contract()
        ledger['application_recipe']['call']='apply_native_passive_profile(controller, ledger_or_explicit_profile) after constructing actual sourcecontroller and before preparing states; no runner hook installed.'
        if key.startswith('W2L'):
            ledger['existing_reference_configuration']=evidence(VIRTUAL/'GENTLE_HIP_0P6_FLEX_0P5.json')
            ledger['source_evidence'].append(evidence(native_source.parent/'CANONICAL_WORKING_PROFILE.json'))
            ledger['clock_contract'].update(canonical_authored_clock_scale=2,canonical_parameters_already_serialized=True,runtime_neural_time_scale_default=1.,do_not_multiply_authored_clock_twice=True,source_time_scaler='native_kernel.NativeController multipliesTimeConstandCaactivation/deactivationtimeconstants byneural_time_scale, includingpassiveafferents/MNs.',separation_recipe='Retain declaredRG/PFionic clocks; override onlyselectedplain afferent/MN tau_phys throughthisadapter usingactualruntime lambda. Do notscale sourceparamswholesale tohumanmass/cadence.')
            ledger['preserve_explicit_changes']=['Sourceclock and imported sourcehash','anklePF targetswap is routing, not a universal biological correction','sharedPF→MN/contact coefficients and groundedhome remain separate fromFSAderivation']
        else:
            ledger['source_evidence'] += [evidence(RAT/'scripts/rat134_receptor_runtime.py'),evidence(RAT/'scripts/rat134_core.py'),evidence(RAT/'prototype_role_policy_B/refinement12/EXPERIMENT_RECORD.md')]
            ledger['known_scope']='Curated134GA own-receptor engineering variant ofpublic120cellprecursor, notfinal2023replication. AcceptedintrinsicA and rejectedBheldouts remain unchanged.'
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
    sources['CAD262_hybrid']=cad
    pair=copy.deepcopy(cad);pair['family']='CAD262_MN_pressure_pair';pair['source_evidence'] += [evidence(ROBOT/'bpa_neural_coordination/cad_pair_unclamped.py'),evidence(EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/RECRUITMENT_PREREGISTRATION.json'),evidence(EXTERNAL/'cad_bpa/pair_neural_recruitment/normalized_unclamped/support085/mixed/ASSUMED_PAIR_CONFIGURATION.json')]
    pair['known_source_function']='Six mirrored existing PF_E1/PF_E2/IBEXC→vas_medMN conductances normalized; pairPOSTzero, no pairoutputfloor. IBEXCdominantRGinput isadditivefeedforward, notmultiplicativeloadgate. Other52hybrid remains.'
    pair['preserved_output_contract']=dict(MN_voltage_range_mV=[0.,5.],normalized_command='u=clip(V_MN/5,0,1)',reference_tension='4530*u N for eachselectedroute',teaching_operating_target_u=.05,teaching_target_not_actual_forcefloor=True,pressure_interface='Sameinverseforce/series-equilibrium pressureconversion,620kPacap,40/60msfillvent,4msdelay. No pressurestate/force projection.',physical_actuator_clock='Remainsphysicaltime; notchangedbyneuralrateprofiles.',calibration_scope='Sixengineeringweightswere fitted tosourceoperatingpoint, not fullfunctionderived range/transientproof.')
    relative={'PF_E1':.2,'PF_E2':.3,'IBEXC_knee_ext':.6};target_v=.05*5.;g_total=endpoint_design(1.,5.,0.,target_v,20.)['gmax_uS'];pair_edges=[]
    for side in ('r','l'):
        for pre,weight in relative.items():pair_edges.append(dict(source=pre+'_'+side,target='MN_vas_med_'+side,gmax_uS=g_total*weight/sum(relative.values())))
    pair['opt_in_profile']=dict(adapter='sns_passive',adopted=False,engineering_candidate=True,neurons=[dict(name='MN_vas_med_'+side,tau_s=.03,operating_synaptic_uS=g_total,response_policy='exact_discrete_pole') for side in ('r','l')],synapses=pair_edges)
    pair['proposed_functional_components']=dict(role='Isolatedexistingthree-inputfeedforwardrecruitmentbranch only; otherneuralcurrentszeroatthisoperatingpoint.',input='s=(.2*sPF1+.3*sPF2+.6*sIBEXC)/1.1, eachs=clip(Vpre/5,0,1)',desired_endpoint_MN_mV=target_v,desired_endpoint_command=.05,desired_fullinput_pole_s=.03,derived_total_gmax_uS=g_total,source_ratios_preserved=relative,force_reference_at_endpoint_N=4530*.05,not_force_floor=True,not_pressure_state_assignment=True,not_wholebody_gainclaim=True,rationale='Retainprevious5%nominalteachingpointand30msMNresponse asdeclaredlocalengineeringcontract; deriveendpointconductancesexactlyratherthanreusetunednormalization. PreservepressureconversionandliteralSNSrecruitment.')
    pair['blockers'] += ['Only2virtualpressureactuators; frozengeometry fails repairedstance domain. No mount/hardware/dynamic calibration.', 'Signed harness/support andcrossedlegacyfeet disallow cleanstanding/pairautonomous claims.']
    sources['CAD262_MN_pressure_pair']=pair
    hybrid=base('W2L704_hybrid_spiking','retained_CaLIF',[VIRTUAL/'spiking_central_runtime.py',VIRTUAL/'spiking_sensory.py',EXTERNAL/'spiking/COMBINED_SPIKING_CHECKPOINT.json',EXTERNAL/'spiking/CENTRAL_RESET_REST_CANDIDATE.json'],sources['W2L108_graded']['controlled_subspace'],
                ['Ca/SFAcurrent burst mechanism is outsidepureGLIFTable2; no automatic m/tau_theta substitution.',
                 'Readoutnormalizer55.376Hz is not verifiedFmax over a declaredinputdomain.',
                 'NoEq27delta/Eq29functionderived synaptic port contract for allreplacedpaths; keep source spike-decoder architecture distinct.',
                 'All source-support/3Dcoverage gaps remain.'])
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
                     'Onlysupportedisolatedjoint/reflex evidence; no3Dhuman/CADcoveragebyinheritance.'])
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
    ode=base('reduced_Deng14','ode14_constants_not_runtime_hooked',[SPINAL/'deng_cpg_ode.py',SPINAL/'deng_cpg_numpy.py',SIMULINK/'demos/sns_build_deng_cpg.m'],
             dict(mechanics='Reduced neural demonstration only; RG/PF/IN andhipMNs, no completebodyloop.',source_neurons='2RGHC+4PFHC+6IN+2hipMN=14; knee/anklePFsharesonebranch; nofullMNarray.',uncontrolled_axes='All mechanicalaxes requireexplicitbodycoupling; no23DoFclaim.',free_3D_balance_demonstrated=False),
             ['KeepNaP6cells/sourceh-gatemechanism distinctfromGLIFadaptiveθ; fixedtau_h350ms is nottau_theta.',
              'ODE usesms internally(t_end=20000for20s); legacySNSNumPycompileDTsecondswithCM5 requiresclockaudit beforeequivalenceclaim.',
              'PythonODEandSimulinkbuilderPF→MNconductances differ (2.565/3.632 vs2.6563/3.6325); chooseversionandtargetexplicitly, no silentlymergeddefaults.',
              'PlainIN/MNendpoint/transienttargetscanbederived, but originalODEglobalCM/GMhasnoselectiveprofilehook.'])
    ode['units']=dict(ODE='ms,mV,nA,nF,uS; physicaltau=5nF/1uS=5ms',Simulink='secondswith1000/Cm; gate350ms explicit',legacy_numpy='Do not execute top-level20s script fordesignimport; clock/capacitance mismatch flagged.')
    ode['source_parameters']=dict(capacitance_nF=5,leak_uS=1,rest_HC_mV=-60,rest_MN_mV=-100,NaP_uS=1.5,NaP_reversal_mV=50,h_time_constant_s=.350)
    sources['reduced_Deng14']=ode
    index=[]
    for key,ledger in sources.items():
        path=LEDGERS/(key+'.json');path.write_text(json.dumps(ledger,indent=2),encoding='utf-8')
        index.append(dict(family=key,path=str(path),design_passed=False,scaling_passed=False,pathways_passed=False,blocker_count=len(ledger['blockers'])))
    gate=json.loads((HERE/'CAMPAIGN_GATE.json').read_text())
    covered=['W2L108_graded','W2L704_hybrid_spiking','rat134_GA_receptor_variant','CAD262_hybrid_and_neural_pressure_pair','SNS_Simscape_graded_reflexes','SNS_Simscape_spiking_reflexes','reduced_Deng14']
    assert set(covered)<=set(gate['required_families'])
    (HERE/'NON_GAIT_DESIGN_INDEX.json').write_text(json.dumps(dict(schema='aarl_non_gait_design_index_v1',covered_global_families=covered,ledgers=index,source_inventory_only=True,training_allowed=False,all_family_design_verified=False,no_model_copies=True),indent=2))
    print(json.dumps({key:s.get('source_census','contractrecipe') for key,s in sources.items()},indent=2))


def electrical_probes():
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
    (HERE/'NON_GAIT_ADAPTER_ELECTRICAL_CHECK.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(passed=bool(unit_pass),sns_error=abs(expected-measured),native_error=native_error,physics_steps=0),indent=2))
    if not unit_pass:raise AssertionError('Adapter electrical check failed; do not promote')


def candidate_component_probes():
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
        untouched=np.array_equal(c.params[:,[0,2,3,4,5,6,7,8,9,10]],old_params[:,[0,2,3,4,5,6,7,8,9,10]]) and np.array_equal(c.edges,old_edges)
        passed=bool(checks) and max(x['error_mV'] for x in checks)<1e-10 and untouched
        results[family]=dict(passed=passed,scope='Selectedactualpure-inputpassiverelays; everyadmissibleendpoint/midpoint at5timeconstants under5xbiologicalrate withphysical10mspole.',contract_count=len(ledger['proposed_functional_components']),source_Size_ionic_gates_and_edges_exact=untouched,consumption=proof,checks=checks,whole_network_pass=False)
    net=cad_network()
    from exact_fast_eval import ExactFastForward
    names=sorted(net.idx,key=net.idx.get);c=net.compiled;wrapper=ExactFastForward(c)
    pair_ledger=json.loads((LEDGERS/'CAD262_MN_pressure_pair.json').read_text());profile=pair_ledger['opt_in_profile']
    proof=apply_sns_passive_profile(wrapper,names,profile)
    posts=np.array([net.idx['MN_vas_med_'+side] for side in ('r','l')]);held=np.ones(len(names),dtype=bool);held[posts]=False
    pres=np.array([net.idx[change['source']] for change in profile['synapses']]);gt=pair_ledger['proposed_functional_components']['derived_total_gmax_uS'];curve=[];first_full=None
    for activation in (0.,.25,.5,.75,1.):
        c.V[:]=c.V_rest
        for tick in range(500):
            c.V[held]=c.V_rest[held];c.V[pres]=5*activation
            wrapper.forward(np.zeros(len(net.inputs)))
            if activation==1. and tick==0:first_full=float(c.V[posts[0]])
        expected=20*gt*activation/(1+gt*activation)
        curve.append(dict(normalized_weighted_input=activation,expected_MN_mV=expected,measured_MN_mV=c.V[posts].tolist(),maximum_error_mV=float(np.max(abs(c.V[posts]-expected))),command=(c.V[posts]/5).tolist()))
    measured_tau=-c.dt/math.log1p(-first_full/.25)
    s=np.linspace(0,1,1001);max_nonlinearity=float(np.max(abs(20*gt*s/(1+gt*s)-.25*s))/.25)
    passed=max(row['maximum_error_mV'] for row in curve)<1e-10 and abs(measured_tau-.03)<1e-12 and max_nonlinearity<=.005
    results['CAD262_MN_pressure_pair']=dict(passed=bool(passed),scope='Actualsource3inputbranch/twoMNvoltageclampisolation; otherpresynapticcellsheldatrest. Thisisnotfree-runningwholeCAD.',consumption=proof,curve=curve,physical_fullinput_discrete_pole_s=measured_tau,endpoint_command=.05,endpoint_reference_tension_N=226.5,maximum_graded_linearity_error_fullscale=max_nonlinearity,linearity_limit=.005,pressure_conversion_modified=False,force_floor_applied=False,body_steps=0)
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
    output=dict(schema='aarl_non_gait_candidate_component_checks_v1',passed=all(row['passed'] for row in results.values()),components=results,physics_steps=0,training_steps=0,source=evidence(Path(__file__)),global_gate_promoted=False)
    (HERE/'NON_GAIT_CANDIDATE_COMPONENT_CHECKS.json').write_text(json.dumps(output,indent=2))
    print(json.dumps({key:dict(passed=row['passed'],contracts=row.get('contract_count'),scope=row['scope']) for key,row in results.items()},indent=2))
    if not output['passed']:raise AssertionError('Candidatecomponentfailed;preservedwithoutgatepromotion')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');parser.add_argument('--electrical-probes',action='store_true');parser.add_argument('--candidate-probes',action='store_true');args=parser.parse_args()
    if args.generate:generate()
    if args.electrical_probes:electrical_probes()
    if args.candidate_probes:candidate_component_probes()
