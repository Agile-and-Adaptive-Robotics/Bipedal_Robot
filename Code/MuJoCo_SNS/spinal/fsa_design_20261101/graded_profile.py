"""Opt-in graded electrical design; no Hunt or body-training entry point.

This repository runs SNS in seconds: capacitance numbers are microfarads,
conductance is microSiemens, voltage is mV, and current is nA. Only true
input encoders receive the engineering fast-relay time constant. NaP RG,
PF shape, laminated timing, MNs, and muscle dynamics are preserved.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import sys
if __name__=='__main__':
    for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np

HERE=Path(__file__).resolve().parent
FAMILIES=('s3k','w2lvar','syn6')
SCHEMA='aarl_graded_fsa_profile_v1'
KINDS=('Ia_homonymous','II_homonymous','Ib_autogenic','Ia_reciprocal_pool')
FIELDS=('c_m','g_m','V_rest','i_b','time_factor_membrane','g_max_non','del_e','e_lo','e_hi','input_connectivity')
ION_FIELDS=('g_ion','e_ion','pow_a','k_a','slope_a','e_a','pow_b','k_b','slope_b','e_b','tau_max_b','pow_c','k_c','slope_c','e_c','tau_max_c')


def exact_keys(value,allowed,where):
    if not isinstance(value,dict) or set(value)!=set(allowed):
        raise ValueError(f'{where}: missing/unused keys {set(allowed)^set(value) if isinstance(value,dict) else type(value)}')


def default_profile(family):
    if family not in FAMILIES:raise ValueError(f'Unsupported graded family {family}')
    return dict(schema=SCHEMA,family=family,units='sns_seconds_microfarad',
        context=dict(fast_stride_s=.559,kinematic_bandwidth_hz=4/.559,phase_budget_deg=10.,
                     interpretation='Engineering per-encoder budget; angle spectrum is not a measured neural bandwidth'),
        probe_current_nA=1.,physical_input_domains_nA=None,operating_post_mV=2.5,
        contracts={k:(.5 if k in KINDS[:2] else -.5) for k in KINDS})


def load_profile(value,family):
    p=json.loads(Path(value).read_text()) if isinstance(value,(str,Path)) else copy.deepcopy(value)
    exact_keys(p,('schema','family','units','context','probe_current_nA','physical_input_domains_nA','operating_post_mV','contracts'),'profile')
    exact_keys(p['context'],('fast_stride_s','kinematic_bandwidth_hz','phase_budget_deg','interpretation'),'context')
    exact_keys(p['contracts'],KINDS,'contracts')
    if p['schema']!=SCHEMA or p['family']!=family or p['units']!='sns_seconds_microfarad':raise ValueError('Profile schema/family/unit mismatch')
    numbers=[p['probe_current_nA'],p['operating_post_mV'],*p['contracts'].values(),p['context']['fast_stride_s'],p['context']['kinematic_bandwidth_hz'],p['context']['phase_budget_deg']]
    if not np.isfinite(numbers).all() or p['probe_current_nA']<=0 or p['context']['fast_stride_s']<=0 or p['context']['kinematic_bandwidth_hz']<=0 or not 0<p['context']['phase_budget_deg']<90:raise ValueError('Invalid profile domain')
    if not p['context']['interpretation'].strip():raise ValueError('Design context must be declared')
    if any(p['contracts'][k]<=0 for k in KINDS[:2]) or any(p['contracts'][k]>=0 for k in KINDS[2:]):raise ValueError('Functional contract sign mismatch')
    return p


def assert_campaign_allowed(action):
    gate=json.loads((HERE/'CAMPAIGN_GATE.json').read_text())
    operation={'Hunt_parameter_training':'Hunt_training','Hunt_training':'Hunt_training','curriculum_simulation':'curriculum','curriculum':'curriculum','whole_body_optimization':'whole_body_optimization'}.get(action)
    if operation:
        from .design_gate import require_calculation_release
        return require_calculation_release(operation)
    if action not in gate['allowed_during_design']:raise ValueError(f'Unknown design action {action}')


def raw_compiled(spinal):
    c=spinal.compiled
    return getattr(c,'c',c)


def field_digest(c):
    h=hashlib.sha256()
    h.update(np.float64(c.dt).tobytes())
    for field in FIELDS+ION_FIELDS:
        if hasattr(c,field):h.update(field.encode());h.update(np.ascontiguousarray(getattr(c,field)).tobytes())
    return h.hexdigest()


def population_names(spinal):
    if any(p['number']!=1 for p in spinal.net.populations):raise ValueError('Scalar Gait populations required')
    names=[p['name'] for p in spinal.net.populations]
    if len(names)!=len(raw_compiled(spinal).V) or len(set(names))!=len(names):raise ValueError('Actual population census mismatch')
    return names


def edge_kind(pre,post,reversal):
    if not post.startswith('MN_'):return None
    act=post[3:]
    if pre=='Ia_'+act and reversal>0:return KINDS[0]
    if pre in ('II_'+act,'IIX_'+act) and reversal>0:return KINDS[1]
    if pre in ('Ib_'+act,'IBIN_'+act) and reversal<0:return KINDS[2]
    if pre.startswith(('Ia_','IaIN_')) and pre.split('_',1)[1]!=act and reversal<0:return KINDS[3]
    return None


def activation_fraction(c,post,pre,voltage):
    span=c.e_hi[post,pre]-c.e_lo[post,pre]
    if span<=0:raise ValueError('Nonpositive synaptic voltage window')
    return float(np.clip((voltage-c.e_lo[post,pre])/span,0.,1.))


def isolated_source_voltage(c,names,pre,current):
    """One afferent channel on, all other sources off; no recurrent-pool claim."""
    name=names[pre]
    if name.startswith(('Ia_','II_','Ib_')):
        if np.any(c.g_max_non[pre]):raise ValueError('Input encoder has incoming conductance; isolated RC assumption invalid')
        return float(c.V_rest[pre]+(current+c.i_b[pre])/c.g_m[pre]),dict(source=name,input_current_nA=current,kind='direct_encoder')
    prefix,act=name.split('_',1)
    aff={'IIX':'II_','IBIN':'Ib_','IaIN':'Ia_'}.get(prefix)
    if aff is None:raise ValueError(f'Unknown sensory terminal stage {name}')
    source=names.index(aff+act)
    u,source_ledger=isolated_source_voltage(c,names,source,current)
    if c.g_max_non[pre,source]<=0:raise ValueError('Missing homologous afferent relay')
    g=c.g_max_non[pre,source]*activation_fraction(c,pre,source,u)
    v=(c.g_m[pre]*c.V_rest[pre]+c.i_b[pre]+g*c.del_e[pre,source])/(c.g_m[pre]+g)
    return float(v),dict(source=name,input_current_nA=current,kind='isolated_afferent_relay',upstream=source_ledger,upstream_voltage_mV=u,relay_active_conductance_uS=float(g),other_inputs='zero; recurrent and PF interactions excluded')


def canonical_edge(pre,post):
    def split(name):return (name[:-2],name[-1]) if name.endswith(('_r','_l')) else (name,None)
    a,sa=split(pre);b,sb=split(post)
    return a,b,('ipsilateral' if sa==sb else 'contralateral'),sb


def apply_profile(spinal,value,*,family):
    """Derive from FINAL compiled values; also usable after an explicit tie pass.

    No edges or populations are created. Missing homologues are rejected,
    not repaired. The original network remains the baseline representation.
    """
    assert_campaign_allowed('optional_profile_implementation')
    p=load_profile(value,family);c=raw_compiled(spinal);names=population_names(spinal)
    domains=p['physical_input_domains_nA']
    if domains is not None:
        expected={n for n in spinal.inputs if n.startswith(('Ia_','II_','Ib_'))};exact_keys(domains,expected,'physical input domains')
        if any(len(v)!=2 or not np.isfinite(v).all() or v[0]>v[1] for v in domains.values()):raise ValueError('Invalid physical input bounds')
    if c.spiking or c.electrical or c.electrical_rectified:raise ValueError('This profile covers graded chemical Gait networks only')
    before=field_digest(c);cap=c.c_m.copy();conductance=c.g_max_non.copy()
    omega=2*math.pi*p['context']['kinematic_bandwidth_hz'];phi=math.radians(p['context']['phase_budget_deg']);theta=omega*c.dt
    continuous_tau=math.tan(phi)/omega
    if phi<=theta:raise ValueError('Phase budget is smaller than one explicit update latency at this neural step')
    tangent=math.tan(phi-theta);pole=tangent/(math.sin(theta)+tangent*math.cos(theta));tau=c.dt/(1-pole)
    encoder=[]
    for i,name in enumerate(names):
        if not name.startswith(('Ia_','II_','Ib_')):continue
        pop=spinal.net.populations[i]['type'].params
        if pop['spiking'] or pop['gated'] or np.any(conductance[i]) or c.V_rest[i]!=0 or c.i_b[i]!=0:raise ValueError(f'Encoder {name} violates passive zero-rest isolated contract')
        if c.dt/tau>1:raise ValueError('Fast encoder is not monotone Euler-resolved at this neural step; do not silently reduce dt')
        cap[i]=tau*c.g_m[i]
        outgoing=np.flatnonzero(conductance[:,i]);windows=sorted({(float(c.e_lo[q,i]),float(c.e_hi[q,i])) for q in outgoing})
        if not windows:raise ValueError('Selected encoder has no functional output')
        interval=[float(c.g_m[i]*(min(a for a,b in windows)-c.V_rest[i])-c.i_b[i]),float(c.g_m[i]*(max(b for a,b in windows)-c.V_rest[i])-c.i_b[i])]
        encoder.append(dict(name=name,index=i,old_tau_s=float(c.c_m[i]/c.g_m[i]),tau_s=tau,capacitance_backend_microF=float(cap[i]),capacitance_nF=float(1000*cap[i]),probe_current_nA=p['probe_current_nA'],physical_input_domain_nA=None if domains is None else domains[name],admissible_engineering_current_domain_nA=interval,outgoing_voltage_windows_mV=windows,domain_evidence='Engineering interval derived from actual outgoing windows and source port leak,not an absolute physical maximum; no current clamp introduced',Euler_leak_factor=float(c.dt/tau)))
    if not encoder:raise ValueError('Profile selected no actual encoders')
    pools={};post,pre=np.nonzero(conductance)
    for dst,src in zip(post,pre):
        kind=edge_kind(names[src],names[dst],c.del_e[dst,src])
        if kind:pools.setdefault((kind,int(dst)),[]).append(int(src))
    counts={k:sum(kind==k for kind,_ in pools) for k in KINDS}
    if any(v==0 for v in counts.values()):raise ValueError(f'Unused functional contracts: {counts}')
    pool_ledger=[];edge_ledger=[]
    for (kind,dst),sources in pools.items():
        base=p['operating_post_mV'];delta=p['contracts'][kind];target=base+delta
        if c.V_rest[dst]!=0 or c.i_b[dst]!=0:raise ValueError('Pool operating-point derivation requires zero intrinsic bias/rest')
        source_values=[];source_evidence=[]
        for src in sources:
            u,evidence=isolated_source_voltage(c,names,src,p['probe_current_nA']);source_values.append(u);source_evidence.append(evidence)
        fractions=np.array([activation_fraction(c,dst,src,u) for src,u in zip(sources,source_values)])
        old=conductance[dst,sources];denominator=float(np.sum(old*fractions*(c.del_e[dst,sources]-target)))
        scale=float(c.g_m[dst]*delta/denominator) if denominator else -1.
        if not np.isfinite(scale) or scale<=0:raise ValueError(f'Unreachable signed operating-point contract: {kind}/{names[dst]}')
        conductance[dst,sources]=old*scale
        for src,u,frac,old_g in zip(sources,source_values,fractions,old):
            edge_ledger.append(dict(kind=kind,pre=names[src],post=names[dst],pre_index=src,post_index=dst,source_voltage_mV=u,fraction_at_unit_input=float(frac),old_g_uS=float(old_g),g_uS=float(conductance[dst,src]),reversal_mV=float(c.del_e[dst,src])))
        pool_ledger.append(dict(kind=kind,post=names[dst],post_index=dst,presynaptic_indices=sources,operating_post_mV=base,target_delta_mV=delta,target_post_mV=target,bias_current_for_probe_nA=float(c.g_m[dst]*base),scale=scale,source_voltage_mV=source_values,source_evidence=source_evidence,fan_in=len(sources),allocation='Preserve relative compiled weights; solve signed shared pool budget including voltage-dependent shunting,actual source range,and reversal.',scope='Isolated terminal pool; simultaneous sensory/PF/recurrent interactions not assumed linear'))
    homologues={}
    for e in edge_ledger:
        a,b,relation,side=canonical_edge(e['pre'],e['post'])
        if side not in ('r','l'):raise ValueError('Modified sensory edge has no explicit bilateral binding')
        homologues.setdefault((a,b,relation),{})[side]=e['g_uS']
    for key,pair in homologues.items():
        if set(pair)!=set(('r','l')) or not math.isclose(pair['r'],pair['l'],rel_tol=1e-12,abs_tol=1e-14):raise ValueError(f'Unmatched/asymmetric derived source edge {key}: {pair}; no silent mirror repair')
    wrapper=spinal.compiled
    if wrapper is not c and not all(hasattr(wrapper,k) for k in ('post','pre','gmax','elo','erange','reversal')):raise ValueError('Unknown accelerated cache; reapply before wrapping')
    # Commit only after every contract and bilateral check passes.
    c.c_m[:]=cap;c.time_factor_membrane[:]=c.dt/cap;c.g_max_non[:]=conductance
    if wrapper is not c:
        wrapper.gmax[:]=c.g_max_non[wrapper.post,wrapper.pre]
        wrapper.elo[:]=c.e_lo[wrapper.post,wrapper.pre];wrapper.erange[:]=(c.e_hi-c.e_lo)[wrapper.post,wrapper.pre];wrapper.reversal[:]=c.del_e[wrapper.post,wrapper.pre]
    ledger=dict(schema=SCHEMA,family=family,profile=p,compiled_before_sha256=before,compiled_after_sha256=field_digest(c),encoder_tau_s=tau,continuous_encoder_tau_bound_s=continuous_tau,discrete_pole_time_constant_s=-c.dt/math.log(pole),phase_convention='Input sampled at t_n,held during update; output state timestamp t_n+dt. One update latency included. An additional stale sensor sample of age0..dt adds0..360*f*dt degrees,not included in per-encoder budget.',encoders=encoder,pools=pool_ledger,modified_edges=edge_ledger,unit_adapter=dict(time='s',voltage='mV',current='nA',conductance='uS',capacitance_backend='microF',physical_capacitance='nF=1000*backend',example='0.03microF/1uS=0.03s=30ms; physical30nF'),applied_to='Final compiled arrays; builder graph numerical labels are baseline values',preserved='All NaP parameters and RG/PF/lamination/MN capacitances,PF shape,topology,and source NMF K6',bilateral_modified_edge_pairs=len(homologues),global_training_gate=False,cost_ledger=dict(changed_capacitances=len(encoder),changed_synapses=len(edge_ledger),new_neurons=0,new_edges=0,extra_forward_operations=0,unit_adaptation='Optional input conversion/validation is outside the neural kernel'))
    spinal.graded_fsa_profile=copy.deepcopy(p);spinal.graded_fsa_ledger=ledger
    ledger.update(physical_input_domains_verified=False,functional_target_scope='Engineering isolated +/-0.1 of the5mV MN activation span around50% activation,not a physiological observation or Hunt fit; local conductance-dependent operating-point contract only.')
    verify_profile(spinal)
    return ledger


def verify_profile(spinal):
    """Must run after final mirror/tie/accelerator installation, before use."""
    if not hasattr(spinal,'graded_fsa_ledger'):raise ValueError('No applied graded profile')
    c=raw_compiled(spinal);ledger=spinal.graded_fsa_ledger
    if field_digest(c)!=ledger['compiled_after_sha256']:raise ValueError('Final consumed fields changed after profile (including mirror/tie pass); explicitly reapply and reverify')
    wrapper=spinal.compiled
    if wrapper is not c:
        for actual,expected in ((wrapper.gmax,c.g_max_non[wrapper.post,wrapper.pre]),(wrapper.elo,c.e_lo[wrapper.post,wrapper.pre]),(wrapper.erange,(c.e_hi-c.e_lo)[wrapper.post,wrapper.pre]),(wrapper.reversal,c.del_e[wrapper.post,wrapper.pre])):
            if not np.array_equal(actual,expected):raise ValueError('Accelerated edge cache differs from consumed profile')
    return dict(final_consumed_fields_match=True,compiled_sha256=field_digest(c))


def adapt_inputs(spinal,values,*,unit='nA',enforce_domain=True):
    verify_profile(spinal)
    scales={'nA':1.,'pA':.001,'uA':1000.}
    if unit not in scales:raise ValueError('Input adapter requires a declared current unit')
    x=np.asarray(values,float)*scales[unit]
    if x.shape!=(len(spinal.inputs),) or not np.isfinite(x).all():raise ValueError('Input shape/nonfinite error')
    selected=np.array([i for i,n in enumerate(spinal.inputs) if n.startswith(('Ia_','II_','Ib_'))]);domains=spinal.graded_fsa_profile['physical_input_domains_nA']
    if enforce_domain=='engineering':domains={v['name']:v['admissible_engineering_current_domain_nA'] for v in spinal.graded_fsa_ledger['encoders']}
    elif enforce_domain not in (True,False):raise ValueError('Domain selector must be physical(True),engineering,or measurement-only(False)')
    if enforce_domain and domains is None:raise ValueError('Physical input domain unresolved; a1nA calibration pulse does not establish it')
    if enforce_domain and any(not domains[spinal.inputs[i]][0]<=x[i]<=domains[spinal.inputs[i]][1] for i in selected):raise ValueError('Measured afferent current outside declared domain; no silent clipping')
    return x,dict(unit_in=unit,unit_out='nA',afferent_min_nA=float(min(x[selected])),afferent_max_nA=float(max(x[selected])),range_evidence='actual supplied input vector',body_range_validated=False,physical_input_domains_resolved=domains is not None)


def gait_coverage_ledger(spinal):
    names=population_names(spinal)
    return dict(anatomical_floating_base=6,anatomical_internal=17,anatomical_total=23,
        primary_sagittal_reference=['hip_flexion_'+s for s in ('r','l')]+['knee_angle_'+s for s in ('r','l')]+['ankle_angle_'+s for s in ('r','l')],
        other_internal=['hip_adduction_'+s for s in ('r','l')]+['hip_rotation_'+s for s in ('r','l')]+['subtalar_angle_'+s for s in ('r','l')]+['mtp_angle_'+s for s in ('r','l')]+['lumbar_extension','lumbar_bending','lumbar_rotation'],
        primary_reference_scope='Six sagittal coordinates; not seventeen independently validated joint targets',
        actual_populations=dict(graded_afferent_encoders=sum(n.startswith(('Ia_','II_','Ib_')) for n in names),motoneurons=sum(n.startswith('MN_') for n in names),pattern_nodes=[n for n in names if n.startswith('PF_')],balance_nodes=[n for n in names if n.startswith(('BAL_','VEST_'))]),
        route_scope='Muscle paths can create off-sagittal torques; route coverage or NMF43-muscle representation does not establish closed-loop rank or independent control',
        passive_scope='Joint/contact/passive-muscle actions are mechanical contributions,not neural PF controllers',
        apparatus_scope='Existing runs retain horizontal/rotational guides; zero prescribed vertical support does not imply free balance. Soft guides do not remove anatomical DoFs.',
        source_conversion_scope='Observed native105-coordinate model has82 auxiliary follower coordinates in addition to23 anatomical coordinates; native count is not anatomical DoF count',
        static_span_evidence=str(HERE/'syn6_torque_span.py'),static_span_scope='Parent audit: source92/represented86 local activation torque rank17; frozen bilateral K6 activation basis rank12 at4walking poses. This excludes extra reflex/posture channels and proves neither cone feasibility nor closed-loop control.',
        dimensional_extension_verified=False,closed_loop_stability_proved=False)


def probe_pools(spinal):
    """All consumed pools get equilibrium checks; one per kind gets Toolbox dynamics."""
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    from sns_toolbox.connections import NonSpikingSynapse
    verify_profile(spinal);c=raw_compiled(spinal);ledger=spinal.graded_fsa_ledger;rows=[];representatives={}
    for pool in ledger['pools']:
        dst=pool['post_index'];src=np.array(pool['presynaptic_indices']);v=np.array(pool['source_voltage_mV']);frac=np.array([activation_fraction(c,dst,int(q),u) for q,u in zip(src,v)]);g=c.g_max_non[dst,src]*frac
        steady=(pool['bias_current_for_probe_nA']+sum(g*c.del_e[dst,src]))/(c.g_m[dst]+sum(g));error=abs(steady-pool['target_post_mV'])
        rows.append(dict(kind=pool['kind'],post=pool['post'],predicted_post_mV=float(steady),contract_error_mV=float(error),domain='Unit afferent current through isolated homonymous input/relay,then frozen presynaptic terminal voltages'))
        if pool['kind'] in representatives:continue
        n=Network(name='isolated consumed signed pool');count=len(src)
        for k in range(count):n.add_neuron(NonSpikingNeuron(membrane_capacitance=.02),name='pre'+str(k))
        n.add_neuron(NonSpikingNeuron(membrane_capacitance=float(c.c_m[dst]),membrane_conductance=float(c.g_m[dst])),name='post');n.add_input('post')
        for k,q in enumerate(src):n.add_connection(NonSpikingSynapse(max_conductance=float(c.g_max_non[dst,q]),reversal_potential=float(c.del_e[dst,q]),e_lo=float(c.e_lo[dst,q]),e_hi=float(c.e_hi[dst,q])),'pre'+str(k),'post')
        tiny=n.compile(dt=c.dt,backend='numpy');trace=[];maximum_error=0.;G=float(c.g_m[dst]+sum(g));a=1-c.dt*G/c.c_m[dst]
        if not 0<=a<1:raise ValueError('Consumed pool is not monotone Euler-resolved in its declared isolated domain')
        for step in range(500):
            tiny.V[:count]=v;tiny.forward([pool['bias_current_for_probe_nA']]);expected=steady*(1-a**(step+1));maximum_error=max(maximum_error,abs(tiny.V[-1]-expected));trace.append(float(tiny.V[-1]))
        tiny.V[:count]=c.e_hi[dst,src]*2;tiny.forward([pool['bias_current_for_probe_nA']]);saturated=c.g_max_non[dst,src];saturation_target=(pool['bias_current_for_probe_nA']+sum(saturated*c.del_e[dst,src]))/(c.g_m[dst]+sum(saturated))
        representatives[pool['kind']]=dict(fan_in=count,measured_post_min_mV=min(trace),measured_post_max_mV=max(trace),final_contract_error_mV=float(abs(trace[-1]-steady)),transient_error_mV=float(maximum_error),effective_tau_s=float(c.c_m[dst]/G),saturated_presynaptic_fraction_max=float(max(activation_fraction(c,dst,int(q),2*c.e_hi[dst,q]) for q in src)),saturated_equilibrium_mV=float(saturation_target),scope='Actual SNS-Toolbox isolated extracted pathway; presynaptic voltages clamped. No whole-network or body outcome.')
    return dict(all_pool_static_error_max_mV=max(v['contract_error_mV'] for v in rows),pools_checked=len(rows),representatives=representatives,static_rows=rows)


def probe_encoder(spinal):
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    c=raw_compiled(spinal);ledger=spinal.graded_fsa_ledger;i=ledger['encoders'][0]['index'];results=[]
    for frequency in (1/ledger['profile']['context']['fast_stride_s'],ledger['profile']['context']['kinematic_bandwidth_hz']):
        n=Network();n.add_neuron(NonSpikingNeuron(membrane_capacitance=float(c.c_m[i]),membrane_conductance=float(c.g_m[i])),name='encoder');n.add_input('encoder');tiny=n.compile(dt=c.dt,backend='numpy');observed=[];inputs=[]
        for step in range(1500):
            current=.5+.5*math.sin(2*math.pi*frequency*step*c.dt);tiny.forward([current]);observed.append(float(tiny.V[0]));inputs.append(current)
        time=(np.arange(1500)+1)*c.dt;fit=np.linalg.lstsq(np.column_stack((np.ones(1000),np.sin(2*math.pi*frequency*time[500:]),np.cos(2*math.pi*frequency*time[500:]))),observed[500:],rcond=None)[0]
        measured_phase=-math.degrees(math.atan2(fit[2],fit[1]));a=1-c.dt*c.g_m[i]/c.c_m[i];theta=2*math.pi*frequency*c.dt
        H=(1-a)*np.exp(-1j*theta)/(1-a*np.exp(-1j*theta));predicted_phase=-float(np.angle(H))*180/math.pi
        if abs(measured_phase-predicted_phase)>1e-9:raise AssertionError('Discrete encoder phase mismatch')
        results.append(dict(frequency_Hz=frequency,measured_phase_deg=measured_phase,predicted_phase_deg=predicted_phase,measured_gain_mV_per_nA=float(np.hypot(fit[1],fit[2])/.5),input_range_nA=[min(inputs),max(inputs)],output_range_mV=[min(observed),max(observed)],additional_input_sample_age_0_to_dt_phase_deg=[0.,360*frequency*c.dt]))
    pulse=[];tiny.V[:]=0.
    for step in range(75):
        current=20. if 10<=step<25 else 0.;tiny.forward([current]);pulse.append([float((step+1)*c.dt),current,float(tiny.V[0]),float(np.clip(tiny.V[0]/5,0,1))])
    release=[v for v in pulse if v[0]>.05 and v[2]<=5]
    return dict(source='Actual SNS-Toolbox isolated encoder,1500 steps/frequency',results=results,
        out_of_range_pulse=dict(current_nA=20.,current_clipped=False,voltage_peak_mV=max(v[2] for v in pulse),outgoing_fraction_peak=max(v[3] for v in pulse),first_post_release_voltage_below5_mV_s=release[0][0] if release else None,release_time_s=.05,trace=pulse),
        engineering_domain_verified=True,absolute_physical_maximum_inferred=False)


def probe_ranges_and_fanin(spinal):
    c=raw_compiled(spinal);ledger=spinal.graded_fsa_ledger;names=population_names(spinal);representatives={};post_pools={}
    for pool in ledger['pools']:
        post_pools.setdefault(pool['post_index'],[]).append(pool)
        if pool['kind'] in representatives:continue
        dst=pool['post_index'];src=pool['presynaptic_indices'];cases=[]
        for operating in (0.,1.25,2.5,3.75):
            for amplitude in (0.,.25,.5,1.,2.,5.,20.):
                voltage=np.array([isolated_source_voltage(c,names,q,amplitude*ledger['profile']['probe_current_nA'])[0] for q in src]);g=np.array([c.g_max_non[dst,q]*activation_fraction(c,dst,q,u) for q,u in zip(src,voltage)])
                value=(c.g_m[dst]*operating+sum(g*c.del_e[dst,src]))/(c.g_m[dst]+sum(g));linear=operating+amplitude*pool['target_delta_mV']
                cases.append(dict(operating_post_mV=operating,probe_amplitude_ratio=amplitude,post_mV=float(value),delta_mV=float(value-operating),linear_target_error_mV=float(value-linear),within_MN_activation_span=bool(0<=value<=5)))
        representatives[pool['kind']]=dict(cases=cases,maximum_local_linear_target_error_mV=max(abs(v['linear_target_error_mV']) for v in cases),interpretation='Local unit-amplitude target is not a global linear map. Off-point and simultaneous shunting errors remain design evidence,not a family/domain pass.')
    simultaneous=[]
    for dst,pools in post_pools.items():
        for label,subset in (('all',pools),('positive',[p for p in pools if p['target_delta_mV']>0]),('negative',[p for p in pools if p['target_delta_mV']<0])):
            if len(subset)<2:continue
            base=ledger['profile']['operating_post_mV'];G=E=0.
            for pool in subset:
                for q,u in zip(pool['presynaptic_indices'],pool['source_voltage_mV']):
                    g=c.g_max_non[dst,q]*activation_fraction(c,dst,q,u);G+=g;E+=g*c.del_e[dst,q]
            actual=(c.g_m[dst]*base+E)/(c.g_m[dst]+G);additive=base+sum(p['target_delta_mV'] for p in subset)
            simultaneous.append(dict(post=names[dst],subset=label,terminal_pools=len(subset),post_mV=float(actual),sum_of_isolated_targets_mV=float(additive),nonadditivity_mV=float(actual-additive)))
    return dict(amplitude_and_operating_points=representatives,simultaneous_signed_pools=simultaneous,maximum_simultaneous_nonadditivity_mV=max(abs(v['nonadditivity_mV']) for v in simultaneous),engineering_current_domain_nA=[0,5],out_of_range_probe_nA=20,full_physical_domain_verified=False,contract_type='Declared rational/piecewise-saturating terminal map; not a globally linear gain')


def syn6_component_homology(spinal,before_projection):
    """Read-only source-factor timing/support audit; no permutation or refit."""
    import build_network_syn6 as S
    z=np.load(S.SYNERGY_NPZ);f=np.load(S.FSA_PHASE_NPZ);tables=S._load_synergy_tables();names=population_names(spinal);c=raw_compiled(spinal);rows=[]
    right_names=[str(n)[:-2] for n in z['muscle_names_r']];left_names=[str(n)[:-2] for n in z['muscle_names_l']]
    Wr=z['W_r'];Wl=z['W_l'][[left_names.index(n) for n in right_names]]
    cosine=(Wr.T@Wl)/(np.linalg.norm(Wr,axis=0)[:,None]*np.linalg.norm(Wl,axis=0)[None,:])
    phase_corr=np.corrcoef(f['r_pf_phase_mean'].T,f['l_pf_phase_mean'].T)[:6,6:]
    for side in ('r','l'):
        grid=f[side+'_phase_grid'];pm=f[side+'_pf_phase_mean'];stance=pm[grid<50].mean(axis=0);swing=pm[grid>=50].mean(axis=0)
        for k in range(6):
            name=f'PF_S{k+1}_{side}';source=names.index(name);old=np.flatnonzero(before_projection[:,source]);new=np.flatnonzero(c.g_max_non[:,source]);added=[names[q] for q in new if before_projection[q,source]==0]
            rows.append(dict(side=side,component=k+1,raw_H_peak_time_s=float(z['time'][np.argmax(z['H_'+side][k])]),per_leg_FSA_peak_phase_percent=float(tables[side]['peak_phase'][k]),per_leg_FSA_stance_fraction=float(stance[k]/max(stance[k]+swing[k],1e-12)),builder_shared_family=tables[side]['family'][k],builder_tau_multiplier=float(tables[side]['tau'][k]),actual_compiled_tau_s=float(c.c_m[source]/c.g_m[source]),source_outgoing_edges=int(len(old)),projected_outgoing_edges=int(len(new)),added_targets=added))
    return dict(source_basis=str(S.SYNERGY_NPZ),source_basis_sha256=hashlib.sha256(S.SYNERGY_NPZ.read_bytes()).hexdigest(),source_phase_profiles=str(S.FSA_PHASE_NPZ),source_phase_profiles_sha256=hashlib.sha256(S.FSA_PHASE_NPZ.read_bytes()).hexdigest(),source_components=rows,canonical_muscle_aligned_W_cosine=cosine.tolist(),per_leg_FSA_phase_correlation=phase_corr.tolist(),same_index_W_cosine=np.diag(cosine).tolist(),same_index_phase_correlation=np.diag(phase_corr).tolist(),functional_homology_verified=False,permutation_applied=False,release_blocked=True,reason='Same numerical PF index is not verified functional homology across independently fitted NMF factors. Historical mean-family assignment and index-based conductance averaging are engineered choices. Source factors remain unchanged; no canonical-side mirror/permutation candidate is introduced here.')


def electrical_probe(family,model,out,baseline=False):
    # Fresh process per family prevents variant params overlays leaking across families.
    os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['AARL_NET']='' if family=='s3k' else family
    sys.path[:0]=[str(HERE.parent),str(HERE.parent/'w2l_cpg')]
    import mujoco
    import build_network as BN
    nmf_path=HERE.parent/'synergy_basis.npz';nmf_before=hashlib.sha256(nmf_path.read_bytes()).hexdigest()
    m=mujoco.MjModel.from_binary_path(str(model));actuators=[m.actuator(i).name for i in range(m.nu)]
    net=BN.build(actuators);c=raw_compiled(net);rng=np.random.default_rng(1101);h=hashlib.sha256()
    for step in range(128):
        x=rng.uniform(0.,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;c.forward(x);h.update(c.V.tobytes())
    reference=dict(family=family,census=[len(c.V),len(net.inputs),len(net.net.connections)],neural_state_sha256=h.hexdigest(),compiled_sha256=field_digest(c),scope='128 seeded neural-only steps; not a body anchor or rhythm/curriculum pass')
    out.mkdir(parents=True,exist_ok=True)
    if baseline:
        (out/(family+'_baseline.json')).write_text(json.dumps(reference,indent=2));print(json.dumps(reference));return
    saved=json.loads((out/(family+'_baseline.json')).read_text());reference['default_off_exact']=reference==saved
    if not reference['default_off_exact']:raise AssertionError('Default-off neural anchor changed')
    net=BN.build(actuators,graded_profile=default_profile(family));ledger=net.graded_fsa_ledger;initial_ledger=copy.deepcopy(ledger);checks=probe_pools(net);encoder_checks=probe_encoder(net);range_checks=probe_ranges_and_fanin(net)
    if initial_ledger['compiled_before_sha256']!=reference['compiled_sha256']:raise AssertionError('Profile did not start from the canonical compiled baseline')
    # Final mirror/tie pass must invalidate the old ledger,then explicit reapplication must recover it.
    edge=ledger['modified_edges'][0];dst,src=edge['post_index'],edge['pre_index'];net.compiled.g_max_non[dst,src]*=.99
    try:verify_profile(net);rejected=False
    except ValueError:rejected=True
    if not rejected:raise AssertionError('Post-profile weight mutation was not detected')
    ledger=apply_profile(net,default_profile(family),family=family)
    # Exercise the actual standing-family projection without calling its runner.
    sys.path.insert(0,str(HERE.parent.parent/'robot_torso'/'bpa_neural_coordination'))
    import gait_standing_family as F
    before_projection=net.compiled.g_max_non.copy();mask_before=before_projection!=0;tie_before=field_digest(net.compiled);changes=F.tie_weights(net);changed=field_digest(net.compiled)!=tie_before
    try:verify_profile(net);tie_invalidated=False
    except ValueError:tie_invalidated=True
    if changed!=tie_invalidated:raise AssertionError('Actual mirror projection consumption check failed')
    tie_ledger=dict(source=str(Path(F.__file__)),source_sha256=hashlib.sha256(Path(F.__file__).read_bytes()).hexdigest(),changed_pairs=len(changes),added_nonzero_edges=int(np.sum((net.compiled.g_max_non!=0)&~mask_before)),profile_invalidated=tie_invalidated,scope='Existing F.tie_weights behavior audited explicitly; this profile creates no missing homologue')
    homology=syn6_component_homology(net,before_projection) if family=='syn6' else None
    if homology is not None:tie_ledger.update(functional_PF_homology_verified=False,strict_symmetric_release_blocked=True)
    ledger=apply_profile(net,default_profile(family),family=family)
    x,adapter=adapt_inputs(net,np.zeros(len(net.inputs)),unit='pA',enforce_domain=False)
    bad=default_profile(family);bad['unused_typo']=1
    try:load_profile(bad,family);unused_rejected=False
    except ValueError:unused_rejected=True
    from sparse_pairwise_eval import SparsePairwiseForward
    dense=copy.deepcopy(net.compiled);net.compiled=SparsePairwiseForward(net.compiled);verify_profile(net);rng=np.random.default_rng(1102)
    for step in range(128):
        x=rng.uniform(0.,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;dense.forward(x);net.compiled.forward(x)
        if not np.array_equal(dense.V,net.compiled.V):raise AssertionError('Profile dense/accelerated neural trajectory differs')
    if hashlib.sha256(nmf_path.read_bytes()).hexdigest()!=nmf_before:raise AssertionError('Read-only K6 basis changed')
    record=dict(source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reference=reference,initial_ledger=initial_ledger,ledger=ledger,pathway_checks=checks,encoder_checks=encoder_checks,range_and_fanin_checks=range_checks,actual_mirror_projection=tie_ledger,unit_adapter=adapter,unused_key_rejected=unused_rejected,post_tie_mutation_rejected=rejected,final_consumption=verify_profile(net),accelerated_neural128_exact=True,coverage=gait_coverage_ledger(net),scope='Design/unit/isolated-pathway checks only; no Hunt/body/curriculum simulation')
    record['source_NMF_read_only']=dict(path=str(nmf_path),sha256=nmf_before,unchanged=True,K=6,ZCode_K8_untouched=True)
    record['syn6_component_homology']=homology
    (out/(family+'_profile.json')).write_text(json.dumps(default_profile(family),indent=2))
    (out/(family+'_graded_design.json')).write_text(json.dumps(record,indent=2));print(json.dumps(dict(family=family,default_off_exact=True,census=reference['census'],tau_s=ledger['encoder_tau_s'],pools=checks['pools_checked'],maximum_pool_error_mV=checks['all_pool_static_error_max_mV'],final_consumed=True)))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--family',choices=FAMILIES,required=True);parser.add_argument('--model',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--capture-baseline',action='store_true');a=parser.parse_args()
    electrical_probe(a.family,a.model,a.out,a.capture_baseline)
