"""Opt-in Hunt timing encoding: range and feedback strength are independent.

Source controller.py and historical conditions remain unchanged. A bound
timing_features wrapper is consumed by the ORIGINAL live HuntNetwork.step;
only12input-relay capacitances and12existing timing→RG weights change.
No body binding, Hunt fitting, or curriculum entry point is implemented.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from types import MethodType,SimpleNamespace
import numpy as np

HERE=Path(__file__).resolve().parent
DEFAULT_CONDITION=HERE.parent/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/aligned_support_ladder_v1/packaged_configuration.json'
DEFAULT_OUTPUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing\source_equivalent_v2')
SCHEMA='aarl_hunt_normalized_timing_candidate_v2'
LEGACY_SCHEMA='aarl_hunt_normalized_timing_candidate_v1'


def default_profile(mode='function_equivalent',relay_timing='source'):
    if mode not in ('function_equivalent','saturated_source_equivalent','canonical_right_role_budget','balanced_0p5','custom_budgets') or relay_timing not in ('source','fast_role_budget'):raise ValueError('Unknown explicit timing design mode')
    return dict(schema=SCHEMA,feature_clip_max=1.5,reference_side='r',
                mode=mode,relay_timing=relay_timing,operating_RG_voltage_mV=2.5,
                per_RG_budget_nA={'F':.5,'E':.5} if mode in ('balanced_0p5','canonical_right_role_budget') else None,
                context=dict(fast_stride_s=.559,kinematic_bandwidth_hz=4/.559,phase_budget_deg=10.),
                input_policy='feature/1.5 -> normalized 0..1; actual source outgoing window sets R; I=Gm*R*x',
                strength_policy='Separate existing HF→RG conductance budgets; preserve right-half relative weights and mirror left.',
                bandwidth_scope='Engineering fast-relay budget from fourth harmonic of a .559s stride; not measured afferent bandwidth or an intrinsic RG/PF clock change.',
                no_body_or_training=True)


def digest_arrays(core):
    h=hashlib.sha256()
    for key in sorted(core.__dict__):
        value=getattr(core,key)
        if isinstance(value,np.ndarray):h.update(key.encode());h.update(str(value.dtype).encode());h.update(str(value.shape).encode());h.update(value.tobytes())
    return h.hexdigest()


def _clone_forward(source):
    # Wrapper.__getattr__ delegates to .c and is unsafe while deepcopy is
    # constructing an empty instance. Clone its complete dictionary together,
    # retaining shared-buffer aliases through deepcopy's single memo instead.
    target=object.__new__(type(source));target.__dict__.update(copy.deepcopy(source.__dict__))
    return target


def _normalized_features(self,data,side):
    design=self.hunt_timing_normalized_design
    if self.options.get('timing_scale')!=1.:
        raise ValueError('Normalized timing profile owns encoding scale=1; feedback strength belongs in declared synaptic budgets.')
    raw=np.asarray(self._hunt_timing_source_features(data,side),dtype=float)
    if raw.shape!=(6,) or not np.isfinite(raw).all() or np.any(raw<0) or np.any(raw>design['feature_clip_max']+1e-12):
        raise ValueError('Source timing feature contract changed; refuse a second/unknown mapping')
    return raw*np.asarray(design['returned_feature_coefficients'][side])


def install(net,profile=None):
    if profile is None:return dict(applied=False,default_noop=True)
    if hasattr(net,'hunt_timing_normalized_design'):raise ValueError('Timing profile already installed; no double normalization')
    p=json.loads(Path(profile).read_text()) if isinstance(profile,(str,Path)) else copy.deepcopy(profile)
    if p.get('schema')==LEGACY_SCHEMA:
        p.update(schema=SCHEMA,mode='balanced_0p5',relay_timing='fast_role_budget',per_RG_budget_nA={phase:p['coactive_reference_current_per_RG_pool_nA'] for phase in ('F','E')},legacy_profile_schema=LEGACY_SCHEMA)
    if p.get('schema')!=SCHEMA or p.get('feature_clip_max')!=1.5 or p.get('reference_side')!='r' or p.get('mode') not in ('function_equivalent','saturated_source_equivalent','canonical_right_role_budget','balanced_0p5','custom_budgets') or p.get('relay_timing') not in ('source','fast_role_budget'):raise ValueError('Explicit supported timing-profile schema/domain required')
    numbers=[p['operating_RG_voltage_mV'],*p['context'].values()]
    budgets=p.get('per_RG_budget_nA')
    if not np.isfinite(numbers).all() or (budgets is not None and (set(budgets)!=set(('F','E')) or not np.isfinite(list(budgets.values())).all() or min(budgets.values())<=0)):raise ValueError('Invalid timing design')
    if p['mode']=='custom_budgets' and budgets is None:raise ValueError('Custom-budget design requires both explicit RG budgets')
    if p['mode'] in ('balanced_0p5','canonical_right_role_budget') and budgets!={'F':.5,'E':.5}:raise ValueError('Canonical/balanced0.5 is explicitly0.5 per pool; use custom_budgets for other changes')
    if p['mode']=='saturated_source_equivalent' and budgets is not None:raise ValueError('Saturated equivalence retains each original side; no changed current budgets')
    c=net.compiled;core=c.c
    if not all(hasattr(net,key) for key in ('original_g','sync','target_meta','timing_features')):raise ValueError('Actual installed Hunt candidate required')
    method=net.timing_features
    if getattr(method,'__name__','')!='timing_features' or getattr(method,'__module__','')!='controller':raise ValueError('Unknown pre-existing timing wrapper; refuse double mapping')
    omega=2*math.pi*p['context']['kinematic_bandwidth_hz'];phi=math.radians(p['context']['phase_budget_deg']);theta=omega*core.dt
    if phi<=theta:raise ValueError('One source update already exceeds the declared phase budget')
    tangent=math.tan(phi-theta);pole=tangent/(math.sin(theta)+tangent*math.cos(theta));tau=core.dt/(1-pole)
    old_scale=float(net.options.get('timing_scale',1.))
    if not np.isfinite(old_scale) or old_scale<=0:raise ValueError('Positive original input scale required for function-preserving state conversion')
    relays=[];edges=[];coefficients={};groups={};selected=[]
    for side in ('r','l'):
        coefficients[side]=[]
        for feature in range(6):
            source=f'HUNT_AF_{feature}_{side}';pre=net.idx[source];port=net.ports[source]
            if np.any(c.post==pre) or core.g_ion[:,pre].any() or core.V_rest[pre]!=0 or core.i_b[pre]!=0:raise ValueError('Timing relay is not an isolated passive zero-rest input')
            out=np.flatnonzero(c.pre==pre)
            if len(out)!=1:raise ValueError('Expected one existing timing→RG edge per feature')
            edge=int(out[0]);post=int(c.post[edge]);expected=f'RG_'+('F' if feature in (0,2) else 'E')+'_'+side
            if net.names[post]!=expected:raise ValueError('Timing target graph differs from actual source derivation')
            lo=float(core.e_lo[post,pre]);hi=float(core.e_hi[post,pre]);gain=float(core.input_connectivity[pre,port]);gm=float(core.g_m[pre])
            if lo!=0 or hi<=0 or gain<=0 or gm<=0 or np.count_nonzero(core.input_connectivity[pre])!=1:raise ValueError('Unsupported timing input/window convention')
            coefficient=gm*(hi-lo)/(5*gain*p['feature_clip_max']);coefficients[side].append(coefficient)
            old_max_voltage=5*p['feature_clip_max']*old_scale*gain/gm
            if old_max_voltage>hi+1e-12 and p['mode'] not in ('saturated_source_equivalent','canonical_right_role_budget'):raise ValueError('Original feature domain exceeds the linear source window; simple equivalent re-encoding unsupported')
            if p['mode']=='canonical_right_role_budget' and hi!=5.:raise ValueError('Canonical new relay output range requires the declared5mV source interface')
            state_scale=hi/old_max_voltage
            # last_voltage is allocated np.empty by the fast backend and is
            # overwritten before every read. Until it aliases V_last, it is
            # scratch, not a valid prepared neuron state to range-check.
            for array in (core.V,core.V_last):
                if not -1e-12<=array[pre]<=old_max_voltage+1e-12:raise ValueError('Existing timing state lies outside the declared source feature domain; do not project it')
            new_cap=tau*gm if p['relay_timing']=='fast_role_budget' else float(core.c_m[pre])
            relays.append(dict(name=source,index=int(pre),input_port=int(port),input_gain=gain,leak_uS=gm,output_range_mV=[lo,hi],old_max_feature_voltage_mV=old_max_voltage,state_scale=state_scale,old_capacitance_backend=float(core.c_m[pre]),new_capacitance_backend=new_cap,capacitance_backend_unit='uF with second-based dt and uS leak',new_capacitance_nF=1000*new_cap,old_leak_tau_s=float(core.c_m[pre]/gm),new_backend_leak_tau_s=new_cap/gm))
            edges.append(dict(feature=feature,side=side,index=edge,source=source,target=expected,pre=int(pre),post=post,reversal_mV=float(c.reversal[edge]),source_relative_conductance=float(net.original_g[edge]),old_active_conductance_uS=float(c.gmax[edge]),old_threshold_mV=lo,old_saturation_mV=hi,new_threshold_mV=lo*state_scale if p['mode']=='saturated_source_equivalent' else lo,new_saturation_mV=hi*state_scale if p['mode']=='saturated_source_equivalent' else hi))
            selected.append(edge)
    for phase in ('F','E'):
        right=[row for row in edges if row['side']=='r' and row['target']==f'RG_{phase}_r']
        raw=np.array([row['source_relative_conductance'] for row in right])
        if np.any(raw<0) or raw.sum()<=0:raise ValueError('Cannot infer source-relative timing budget from absent/negative right-half weights')
        relative=raw/raw.sum();driving=np.array([row['reversal_mV']-p['operating_RG_voltage_mV'] for row in right])
        if np.any(driving<=0):raise ValueError('Selected source timing edge is not excitatory at the declared operating voltage')
        old_fractions=np.clip(np.array([next(r for r in relays if r['name']==row['source'])['old_max_feature_voltage_mV']/next(r for r in relays if r['name']==row['source'])['output_range_mV'][1] for row in right]),0,1)
        old_budget=float(np.sum(raw*old_fractions*driving))
        requested=old_budget if budgets is None else float(budgets[phase])
        if p['mode']=='function_equivalent' and not np.isclose(requested,old_budget,rtol=1e-12,atol=1e-14):raise ValueError('Function-equivalent mode cannot silently change an old current budget')
        total=requested/float(relative@driving)
        for row,fraction in zip(right,relative):
            for side in ('r','l'):
                target=next(x for x in edges if x['side']==side and x['feature']==row['feature'])
                if target['reversal_mV']!=row['reversal_mV']:raise ValueError('Mirrored reversal mismatch: requires a separate contract')
                if p['mode']=='canonical_right_role_budget':
                    target['new_threshold_mV']=row['old_threshold_mV'];target['new_saturation_mV']=row['old_saturation_mV']
                relay=next(r for r in relays if r['name']==row['source'])
                if p['mode']=='saturated_source_equivalent':value=target['old_active_conductance_uS']
                elif p['mode']=='function_equivalent':
                    mirror_relay=next(r for r in relays if r['name']==target['source'])
                    if target['old_active_conductance_uS']!=row['old_active_conductance_uS'] or mirror_relay['state_scale']!=relay['state_scale']:raise ValueError('Cannot claim both mirrored weights and exact old function when the source is asymmetric')
                    value=row['old_active_conductance_uS']/relay['state_scale']
                else:value=float(total*fraction)
                target['new_conductance_uS']=value
                target['reference_current_nA']=float(value*(target['reversal_mV']-p['operating_RG_voltage_mV']))
        groups[phase]=dict(features=[row['feature'] for row in right],relative_conductances=relative.tolist(),total_conductance_uS=total,old_reference_current_nA=old_budget,requested_reference_current_nA=requested,reference_operating_current_nA=float(total*(relative@driving)),reference_voltage_mV=p['operating_RG_voltage_mV'],not_a_global_current_clamp=True,current_at_RG0mV_nA=float(total*sum(f*r['reversal_mV'] for f,r in zip(relative,right))),current_at_RG5mV_nA=float(total*sum(f*(r['reversal_mV']-5) for f,r in zip(relative,right))))
    pools_by_side={}
    for side in ('r','l'):
        pools_by_side[side]={}
        for phase in ('F','E'):
            rows=[r for r in edges if r['target']==f'RG_{phase}_{side}'];old_current=0.;new_current=0.
            for row in rows:
                relay=next(r for r in relays if r['name']==row['source']);drive=row['reversal_mV']-p['operating_RG_voltage_mV']
                old_fraction=float(np.clip((relay['old_max_feature_voltage_mV']-row['old_threshold_mV'])/(row['old_saturation_mV']-row['old_threshold_mV']),0,1))
                new_fraction=float(np.clip((relay['output_range_mV'][1]-row['new_threshold_mV'])/(row['new_saturation_mV']-row['new_threshold_mV']),0,1))
                old_current+=row['old_active_conductance_uS']*old_fraction*drive;new_current+=row['new_conductance_uS']*new_fraction*drive
                row['reference_current_nA']=row['new_conductance_uS']*new_fraction*drive
            pools_by_side[side][phase]=dict(old_reference_current_nA=old_current,new_reference_current_nA=new_current)
    if p['mode']=='saturated_source_equivalent':
        for phase in ('F','E'):
            rows=[r for r in edges if r['target']==f'RG_{phase}_r'];groups[phase]['total_conductance_uS']=sum(r['new_conductance_uS'] for r in rows)
            groups[phase]['reference_operating_current_nA']=pools_by_side['r'][phase]['new_reference_current_nA']
            for voltage in (0,5):groups[phase][f'current_at_RG{voltage}mV_nA']=sum(r['new_conductance_uS']*np.clip((next(v for v in relays if v['name']==r['source'])['output_range_mV'][1]-r['new_threshold_mV'])/(r['new_saturation_mV']-r['new_threshold_mV']),0,1)*(r['reversal_mV']-voltage) for r in rows)
    design=dict(schema=SCHEMA,profile=p,feature_clip_max=p['feature_clip_max'],returned_feature_coefficients=coefficients,
                old_timing_scale=float(net.options.get('timing_scale',1.)),new_timing_scale=1.,relays=relays,edges=edges,RG_pools=groups,RG_pools_by_side=pools_by_side,
                symmetry_scope='Saturated mode retains original per-side weights, including historical asymmetry. Canonical-right budget copies only the right pattern/windows toleft, neveraverages. Unsaturated equivalent mode requires tied source weights.',
                last_voltage_cache_policy='When not aliased to authoritative V_last, preallocated np.empty scratch is initialized from current V at selected relays; no neuronal-state projection.',
                role_budget_backend_tau_s=tau,encoder_backend_tau_s=relays[0]['new_backend_leak_tau_s'],encoder_discrete_pole=1-core.dt/relays[0]['new_backend_leak_tau_s'],encoder_effective_pole_tau_s=float(-core.dt/math.log(1-core.dt/relays[0]['new_backend_leak_tau_s'])),
                phase_convention='Input at t_n, output state at t_n+dt; one update latency included. Physical sensor sample age is additional.',
                feature_source='Original timing_features: positivegroupmean oflength(mid/half),velocity/.6,or tensileforce/Fmax;clip0..1.5 unchanged.',
                strength='Unsaturated mode uses g_new=g_old/state_scale with fixed windows. Saturated mode retains Gmax and scales presynaptic windows by state_scale, preserving clipping exactly. Budget variants are separate changes.',
                preserved='NaP/PF/banks/phase/physicalsensor maps and all other classes; no bodybinding, no pressure/currentprojection.',
                runtime_hook='Original HuntNetwork.step calls the installed stateless timing_features method; no replacement step or compiled forward.',
                training_allowed=False,curriculum_allowed=False)
    # All source/range/graph checks happen before any mutation.
    backup=dict(cap=core.c_m.copy(),time=core.time_factor_membrane.copy(),g=c.gmax.copy(),original=net.original_g.copy(),options=copy.deepcopy(net.options),voltage=core.V.copy(),voltage_last=core.V_last.copy(),wrapper_last=c.last_voltage.copy(),elo=core.e_lo.copy(),ehi=core.e_hi.copy(),cached_elo=c.elo.copy(),cached_erange=c.erange.copy())
    try:
        if not np.shares_memory(c.last_voltage,core.V_last):
            for relay in relays:c.last_voltage[relay['index']]=core.V[relay['index']]
        for relay in relays:
            i=relay['index'];core.c_m[i]=relay['new_capacitance_backend'];core.time_factor_membrane[i]=core.dt/core.c_m[i]
        seen=set()
        for array in (core.V,core.V_last,c.last_voltage):
            if id(array) in seen:continue
            seen.add(id(array))
            for relay in relays:array[relay['index']]*=relay['state_scale']
        for row in edges:
            c.gmax[row['index']]=row['new_conductance_uS'];net.original_g[row['index']]=row['new_conductance_uS']
            core.e_lo[row['post'],row['pre']]=row['new_threshold_mV'];core.e_hi[row['post'],row['pre']]=row['new_saturation_mV']
            c.elo[row['index']]=row['new_threshold_mV'];c.erange[row['index']]=row['new_saturation_mV']-row['new_threshold_mV']
            net.options['sensory'][row['side']][row['feature']]=row['new_conductance_uS']
        net.options['timing_scale']=1.
        net._hunt_timing_source_features=method
        net.hunt_timing_normalized_design=design
        net.timing_features=MethodType(_normalized_features,net)
        net.sync();verify_runtime(net)
    except Exception:
        core.c_m[:]=backup['cap'];core.time_factor_membrane[:]=backup['time'];c.gmax[:]=backup['g'];net.original_g[:]=backup['original'];net.options=backup['options'];net.timing_features=method
        core.V[:]=backup['voltage'];core.V_last[:]=backup['voltage_last'];c.last_voltage[:]=backup['wrapper_last']
        core.e_lo[:]=backup['elo'];core.e_hi[:]=backup['ehi'];c.elo[:]=backup['cached_elo'];c.erange[:]=backup['cached_erange']
        for field in ('_hunt_timing_source_features','hunt_timing_normalized_design'):net.__dict__.pop(field,None)
        net.sync();raise
    return design


def verify_runtime(net):
    design=net.hunt_timing_normalized_design;c=net.compiled
    if net.options['timing_scale']!=1.:raise AssertionError('Normalized input strength was recompressed')
    for row in design['edges']:
        value=row['new_conductance_uS'];i=row['index']
        if c.gmax[i]!=value or c.c.g_max_non[row['post'],row['pre']]!=value or net.original_g[i]!=value:raise AssertionError('Timing budget overwritten by runtime policy')
        if c.c.e_lo[row['post'],row['pre']]!=row['new_threshold_mV'] or c.c.e_hi[row['post'],row['pre']]!=row['new_saturation_mV'] or c.elo[i]!=row['new_threshold_mV'] or c.erange[i]!=row['new_saturation_mV']-row['new_threshold_mV']:raise AssertionError('Timing presynaptic window/cache overwritten')
        mirror=next(x for x in design['edges'] if x['feature']==row['feature'] and x['side']!=row['side'])
        if design['profile']['mode']!='saturated_source_equivalent' and value!=mirror['new_conductance_uS']:raise AssertionError('Timing budget is not bilaterally tied')
    tied=all(row['new_conductance_uS']==next(x for x in design['edges'] if x['feature']==row['feature'] and x['side']!=row['side'])['new_conductance_uS'] for row in design['edges'])
    return dict(selected_edges=12,selected_encoders=12,bilateral_tied=tied,source_asymmetry_preserved=design['profile']['mode']=='saturated_source_equivalent',authoritative_source_and_sparse_arrays_match=True,presynaptic_window_caches_match=True)


def run_checks(condition,output,profile=None):
    from hunt_pathway_audit import C,MODEL,options_from_physical_condition,mujoco
    cfg=json.loads(condition.read_text());model=mujoco.MjModel.from_xml_path(str(MODEL));acts=[model.actuator(i).name for i in range(model.nu)]
    net=C.build(acts,dt=.002,options=options_from_physical_condition(cfg))
    # PT mechanical calibration is neither invented nor executed in this
    # timing-only electrical fixture. The real candidate keeps its PT object.
    net.__dict__.pop('balance',None);net.apply_weights(1.)
    digest=digest_arrays(net.compiled.c);method=net.timing_features;old_options=copy.deepcopy(net.options)
    noop=install(net,None)
    assert noop['default_noop'] and digest_arrays(net.compiled.c)==digest and net.timing_features==method and net.options==old_options
    physical_maps=(net.sensor_mid.copy(),net.sensor_half.copy(),copy.deepcopy(net.target_meta['feature_members']))
    caps=net.compiled.c.c_m.copy();g=net.compiled.gmax.copy();ionic={key:value.copy() for key,value in net.compiled.c.__dict__.items() if isinstance(value,np.ndarray) and key in ('g_ion','e_ion','tau_max_b','tau_max_c','k_a','k_b','slope_a','slope_b','pow_a','pow_b')}
    source_forward=_clone_forward(net.compiled)
    design=install(net,default_profile() if profile is None else profile);policy=[]
    candidate_forward=_clone_forward(net.compiled)
    changed_edges=np.array([row['index'] for row in design['edges']])
    unchanged_other_weights=np.delete(g,changed_edges).tobytes()==np.delete(net.compiled.gmax,changed_edges).tobytes()
    for gate in (0.,.5,1.,0.):net.apply_weights(gate);policy.append(dict(command_gate=gate,**verify_runtime(net)))
    selected=np.array([row['index'] for row in design['edges']]);relay_ids=np.array([row['index'] for row in design['relays']])
    unchanged_caps=np.delete(caps,relay_ids).tobytes()==np.delete(net.compiled.c.c_m,relay_ids).tobytes()
    # Weight-policy changes may affect other pathways; installation itself was
    # already limited to timing indices, and source physical maps stay exact.
    unchanged_maps=np.array_equal(physical_maps[0],net.sensor_mid) and np.array_equal(physical_maps[1],net.sensor_half) and physical_maps[2]==net.target_meta['feature_members']
    unchanged_ionic=all(np.array_equal(value,getattr(net.compiled.c,key)) for key,value in ionic.items())
    # Actual C.step consumption with synthetic sensor arrays and a fixed zero
    # geometry witness. No bind_mujoco, body state, or physical sensor claim.
    net.Fmax=np.maximum(model.actuator_gainprm[:,2],5.)
    count=len(acts);fixture=SimpleNamespace(actuator_length=net.sensor_mid.copy(),actuator_velocity=np.zeros(count),actuator_force=np.zeros(count),qpos=np.zeros(17),qvel=np.zeros(17))
    net.arm_model=SimpleNamespace(addresses=np.arange(17),rows=np.arange(17),update=lambda data:np.zeros((17,count)))
    net.goal=np.zeros(17);net.live=(None,fixture);net.options['axis_feedback']=0.
    ranges=[]
    for raw in (-1.,0.,.375,.75,1.125,1.5,2.):
        fixture.actuator_length[:]=net.sensor_mid+raw*net.sensor_half;fixture.actuator_velocity[:]=.6*raw;fixture.actuator_force[:]=-net.Fmax*raw
        u=np.zeros(len(net.inputs));u[net.ports['DRIVE']]=2.5
        raw_source={side:net._hunt_timing_source_features(fixture,side).copy() for side in net.sides}
        net.step(u)
        observed={side:[float(u[net.ports[f'HUNT_AF_{i}_{side}']]) for i in range(6)] for side in net.sides}
        expected={side:[float(row['leak_uS']*row['output_range_mV'][1]*raw_source[side][i]/1.5/row['input_gain']) for i,row in enumerate([r for r in design['relays'] if r['name'].endswith('_'+side)])] for side in net.sides}
        error=max(np.max(abs(np.array(observed[side])-expected[side])) for side in net.sides)
        ranges.append(dict(synthetic_source_level=raw,raw_source_features={s:x.tolist() for s,x in raw_source.items()},actual_step_port_currents_nA=observed,maximum_error_nA=float(error),source_clip_preserved=True))
    # Endpoint settling of the actual12sourcecompiledrelaycells. This is a
    # short electrical input-encoder check, not a Hunt rhythm/body test.
    fixture.actuator_length[:]=net.sensor_mid+1.5*net.sensor_half;fixture.actuator_velocity[:]=.9;fixture.actuator_force[:]=-1.5*net.Fmax
    for _ in range(250):
        u=np.zeros(len(net.inputs));u[net.ports['DRIVE']]=2.5;net.step(u)
    endpoint=net.compiled.V[relay_ids].copy();endpoint_error=max(abs(endpoint[i]-row['output_range_mV'][1]) for i,row in enumerate(design['relays']))
    frequency=design['profile']['context']['kinematic_bandwidth_hz'];omega=2*np.pi*frequency;dt=net.compiled.c.dt
    relay_trace=[]
    for step in range(600):
        raw=.75+.6*np.sin(omega*step*dt)
        fixture.actuator_length[:]=net.sensor_mid+raw*net.sensor_half;fixture.actuator_velocity[:]=.6*raw;fixture.actuator_force[:]=-net.Fmax*raw
        u=np.zeros(len(net.inputs));u[net.ports['DRIVE']]=2.5;net.step(u);relay_trace.append(net.compiled.V[relay_ids].copy())
    # Output states belong at t_n+dt; this prevents dropping the actual
    # source's one-update latency by plotting input and output at one timestamp.
    discard=max(100,int(np.ceil(40*max(r['new_backend_leak_tau_s'] for r in design['relays'])/dt)))
    if discard>=550:raise ValueError('Declared relay requires a separately planned longer frequency fixture')
    times=(np.arange(600)+1)*dt;fit=np.column_stack((np.sin(omega*times),np.cos(omega*times),np.ones(600)))[discard:]
    sine_fit=np.linalg.lstsq(fit,np.asarray(relay_trace)[discard:],rcond=None)[0]
    phase_deg=-np.degrees(np.arctan2(sine_fit[1],sine_fit[0]))
    expected_phase=np.array([np.degrees(np.arctan2(np.sin(omega*dt),np.cos(omega*dt)-(1-dt/r['new_backend_leak_tau_s']))) for r in design['relays']])
    phase_error=float(max(abs(phase_deg-expected_phase)))
    # On/off matched single-step states isolate the EXISTING timing currents
    # at RG2.5mV, retaining all ionic fields and the actual sparse kernel.
    base=_clone_forward(net.compiled);base.c.V[:]=base.c.V_rest
    for relay in design['relays']:base.c.V[relay['index']]=relay['output_range_mV'][1]
    rg_ids=sorted(set(row['post'] for row in design['edges']))
    base.c.V[rg_ids]=2.5;on=_clone_forward(base);off=_clone_forward(base);off.gmax[selected]=0.
    zeros=np.zeros(len(net.inputs));on.forward(zeros);off.forward(zeros)
    currents=(on.V[rg_ids]-off.V[rg_ids])/on.c.time_factor_membrane[rg_ids]
    budget=[]
    for i,value in zip(rg_ids,currents):
        wanted=design['RG_pools_by_side'][net.names[i].split('_')[-1]][net.names[i].split('_')[1]]['new_reference_current_nA']
        budget.append(dict(target=net.names[i],measured_incremental_current_nA=float(value),desired_nA=wanted,error_nA=abs(float(value)-wanted)))
    onehot=[]
    for feature in range(6):
        reference=_clone_forward(net.compiled);reference.c.V[relay_ids]=0.;reference.c.V[rg_ids]=2.5
        rows=[row for row in design['edges'] if row['feature']==feature]
        for row in rows:reference.c.V[row['pre']]=5.
        on=_clone_forward(reference);off=_clone_forward(reference);off.gmax[selected]=0.;on.forward(zeros);off.forward(zeros)
        measured={row['side']:float((on.V[row['post']]-off.V[row['post']])/on.c.time_factor_membrane[row['post']]) for row in rows}
        expected={row['side']:row['new_conductance_uS']*np.clip((5.-row['new_threshold_mV'])/(row['new_saturation_mV']-row['new_threshold_mV']),0,1)*(row['reversal_mV']-2.5) for row in rows}
        onehot.append(dict(feature=feature,actual_increment_nA=measured,expected_increment_nA=expected,maximum_error_nA=max(abs(measured[s]-expected[s]) for s in ('r','l')),left_right_difference_nA=abs(measured['r']-measured['l'])))
    equivalence=[]
    if design['profile']['mode'] in ('function_equivalent','saturated_source_equivalent'):
        for feature in (0.,.25,.5,.75,1.,1.25,1.5):
            for voltage in (-2.5,0.,2.5,5.,8.):
                old=_clone_forward(source_forward);new=_clone_forward(candidate_forward)
                for relay in design['relays']:
                    i=relay['index'];old.c.V[i]=5*feature*design['old_timing_scale']*relay['input_gain']/relay['leak_uS'];new.c.V[i]=old.c.V[i]*relay['state_scale']
                old.c.V[rg_ids]=voltage;new.c.V[rg_ids]=voltage
                old.forward(zeros);new.forward(zeros)
                delta=(new.V[rg_ids]-old.V[rg_ids])/old.c.time_factor_membrane[rg_ids]
                equivalence.append(dict(raw_feature=feature,RG_voltage_mV=voltage,maximum_existing_kernel_current_difference_nA=float(max(abs(delta)))))
        # Dynamic encoder/state identity, with the original source pole retained
        # or deliberately changed only in the separate fast-relay variant.
        old=_clone_forward(source_forward);new=_clone_forward(candidate_forward);wave_error=0.;g_error=0.
        for step in range(120):
            feature=(0.,.25,.5,1.5)[(step//15)%4];uo=zeros.copy();un=zeros.copy()
            for relay in design['relays']:
                port=relay['input_port'];uo[port]=5*feature*design['old_timing_scale'];un[port]=relay['leak_uS']*relay['output_range_mV'][1]*feature/1.5/relay['input_gain']
            old.forward(uo);new.forward(un)
            for relay in design['relays']:
                i=relay['index'];wave_error=max(wave_error,abs(float(new.V[i]/relay['state_scale']-old.V[i])))
            for row in design['edges']:
                old_g=row['old_active_conductance_uS']*np.clip((old.V[row['pre']]-row['old_threshold_mV'])/(row['old_saturation_mV']-row['old_threshold_mV']),0,1);new_g=row['new_conductance_uS']*np.clip((new.V[row['pre']]-row['new_threshold_mV'])/(row['new_saturation_mV']-row['new_threshold_mV']),0,1)
                g_error=max(g_error,abs(float(new_g-old_g)))
        dynamic=dict(source_tau_retained=design['profile']['relay_timing']=='source',maximum_rescaled_relay_voltage_error_mV=wave_error,maximum_delivered_conductance_error_uS=g_error,equivalence_required=design['profile']['relay_timing']=='source')
    else:dynamic=dict(equivalence_required=False,reason='Explicit new current-budget architecture, not function-equivalent')
    # Separate source-built fixture checks nonzero prepared-state conversion
    # and rejection before mutation. It never binds a physical body.
    prepared=C.build(acts,dt=.002,options=options_from_physical_condition(cfg));prepared.__dict__.pop('balance',None);prepared.apply_weights(1.)
    first=design['relays'][0];prepared.compiled.c.V[first['index']]=first['old_max_feature_voltage_mV']+.001
    invalid_digest=digest_arrays(prepared.compiled.c);invalid_options=copy.deepcopy(prepared.options);invalid_rejected=False
    try:install(prepared,design['profile'])
    except ValueError:invalid_rejected=True
    invalid_unchanged=digest_arrays(prepared.compiled.c)==invalid_digest and prepared.options==invalid_options
    for array in (prepared.compiled.c.V,prepared.compiled.c.V_last,prepared.compiled.last_voltage):
        for relay in design['relays']:array[relay['index']]=.5*relay['old_max_feature_voltage_mV']
    untouched_voltage=np.delete(prepared.compiled.c.V,relay_ids).copy();explicit_profile=copy.deepcopy(design['profile'])
    if explicit_profile['mode']=='function_equivalent':explicit_profile['per_RG_budget_nA']={phase:design['RG_pools'][phase]['old_reference_current_nA'] for phase in ('F','E')}
    install(prepared,explicit_profile)
    prepared_error=max(max(abs(array[relay_ids]-2.5)) for array in (prepared.compiled.c.V,prepared.compiled.c.V_last,prepared.compiled.last_voltage))
    prepared_audit=dict(nonzero_source_state_mapped_to_half_range=True,maximum_error_mV=float(prepared_error),all_other_voltage_states_exact=np.array_equal(untouched_voltage,np.delete(prepared.compiled.c.V,relay_ids)),out_of_domain_state_rejected=invalid_rejected,rejected_state_unchanged=invalid_unchanged,explicit_per_RG_budgets_consumed=explicit_profile.get('per_RG_budget_nA'))
    duplicate_rejected=False
    try:install(net,default_profile())
    except ValueError:duplicate_rejected=True
    equivalent_grid_pass=not equivalence or max(row['maximum_existing_kernel_current_difference_nA'] for row in equivalence)<1e-10
    equivalent_dynamic_pass=not dynamic['equivalence_required'] or (dynamic['maximum_rescaled_relay_voltage_error_mV']<1e-12 and dynamic['maximum_delivered_conductance_error_uS']<1e-12)
    onehot_pass=max(row['maximum_error_nA'] for row in onehot)<1e-10 and (design['profile']['mode']=='saturated_source_equivalent' or max(row['left_right_difference_nA'] for row in onehot)<1e-10)
    passed=unchanged_caps and unchanged_other_weights and unchanged_maps and unchanged_ionic and max(r['maximum_error_nA'] for r in ranges)<1e-12 and endpoint_error<1e-10 and phase_error<1e-8 and max(row['error_nA'] for row in budget)<1e-10 and duplicate_rejected and equivalent_grid_pass and equivalent_dynamic_pass and prepared_error<1e-12 and prepared_audit['all_other_voltage_states_exact'] and invalid_rejected and invalid_unchanged and onehot_pass
    output.mkdir(parents=True,exist_ok=True)
    result=dict(schema=SCHEMA,machine='easteregg2',passed=bool(passed),condition=str(condition),condition_sha256=hashlib.sha256(condition.read_bytes()).hexdigest(),source_controller=str(Path(C.__file__).resolve()),source_controller_sha256=hashlib.sha256(Path(C.__file__).read_bytes()).hexdigest(),helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),neurons=len(net.compiled.V),default_None_exact=True,command_policy=policy,live_original_C_step_consumed=True,physical_sensor_fixture_only=True,body_binding_calls=0,physics_steps=0,Hunt_training_steps=0,scope='Timing encoding and12existingHF→RG budgets only; no locomotion/wholefamily pass.',ranges=ranges,relay_endpoint_error_mV=float(endpoint_error),RG_reference_current_budgets=budget,all_non_timing_capacitances_preserved=bool(unchanged_caps),ionic_fields_preserved=bool(unchanged_ionic),physical_maps_preserved=bool(unchanged_maps),double_install_rejected=duplicate_rejected,training_allowed=False,curriculum_allowed=False)
    result['all_other_sparse_weights_preserved_at_install']=bool(unchanged_other_weights)
    result['actual_encoder_frequency_check']=dict(frequency_hz=float(frequency),lag_deg_each_relay=phase_deg.tolist(),expected_from_selected_pole_deg=expected_phase.tolist(),maximum_error_vs_selected_pole=phase_error,output_timestamp='t_n+dt; physical sensor sample age not included',synthetic_electrical_samples=600,discarded_startup_samples=discard)
    result['mode']=design['profile']['mode'];result['relay_timing']=design['profile']['relay_timing'];result['old_transfer_current_grid']=equivalence;result['dynamic_encoding_equivalence']=dynamic
    result['prepared_state_conversion']=prepared_audit
    result['selected_presynaptic_windows_verified']=True
    result['bilateral_conductances_tied']=verify_runtime(net)['bilateral_tied']
    result['onehot_feature_currents']=onehot
    (output/'RESULTS.json').write_text(json.dumps(result,indent=2))
    (output/'APPLIED_TIMING_PROFILE.json').write_text(json.dumps(design,indent=2))
    (output/'TIMING_PROFILE.json').write_text(json.dumps(design['profile'],indent=2))
    print(json.dumps({k:result[k] for k in ('passed','neurons','relay_endpoint_error_mV','RG_reference_current_budgets','live_original_C_step_consumed','body_binding_calls')}))
    if not passed:raise AssertionError('Timing-only electrical contract failed; preserve diagnostic, no promotion')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--condition',type=Path,default=DEFAULT_CONDITION);parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT);parser.add_argument('--profile',type=Path);parser.add_argument('--mode',default='function_equivalent',choices=('function_equivalent','saturated_source_equivalent','canonical_right_role_budget','balanced_0p5'));parser.add_argument('--relay-timing',default='source',choices=('source','fast_role_budget'));args=parser.parse_args()
    run_checks(args.condition.resolve(),args.output.resolve(),args.profile.resolve() if args.profile else default_profile(args.mode,args.relay_timing))
