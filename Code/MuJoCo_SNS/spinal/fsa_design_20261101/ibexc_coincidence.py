"""Optional thresholded Ib/phase coincidence in existing IBEXC neurons.

Binary engineering AND, not multiplication or physiological validation.
Explicit per-source current contracts bind F/F0 to the route-specific
Ib->IBEXC window. Autogenic Ib routes and all dynamic states are preserved.
No Hunt runtime hook or body-training entry point is installed.
"""
from __future__ import annotations
import argparse,ast,copy,hashlib,inspect,json,os,sys,textwrap
from pathlib import Path
if __name__=='__main__':
    for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
HERE=Path(__file__).resolve().parent


def utilities():
    try:from . import graded_profile as G
    except ImportError:import graded_profile as G
    return G


def make_policy(family,load_high_current_nA,origin):
    return dict(schema='aarl_ibexc_binary_coincidence_v1',family=family,phase_budget_over_leak=.5,load_budget_over_leak=.5,lower_gap_fraction=.5,upper_gap_fraction=.8,load_high_current_nA=load_high_current_nA,load_contract_origin=origin)


def hunt_identity(spinal):
    cls=type(spinal);module=sys.modules.get(cls.__module__)
    expected=HERE.parent/'hunt_training/controller.py'
    if cls.__name__!='HuntNetwork' or module is None or getattr(module,'HuntNetwork',None) is not cls or Path(module.__file__).resolve()!=expected.resolve():raise ValueError('Actual controller.HuntNetwork identity required')
    lineage={c.__name__:c for c in cls.__mro__}
    for name,path in (('JointIntegrationNetwork',HERE.parent/'functional_gait2392/deng2023_joint_pf.py'),('FunctionalNetwork',HERE.parent/'functional_gait2392/network.py'),('W2LVarnet',HERE.parent/'build_network_w2lvar.py')):
        base=lineage.get(name);owner=sys.modules.get(base.__module__) if base is not None else None
        if owner is None or getattr(owner,name,None) is not base or Path(owner.__file__).resolve()!=path.resolve():raise ValueError('Unverified W2L-derived Hunt MRO')
    G=utilities();names=G.population_names(spinal)
    if len(names) not in (1033,1027) or len([n for n in names if n.startswith('MN_')])!=92:raise ValueError('Unsupported actual Hunt census')
    required={f'RG_{phase}_{side}' for phase in ('E','F') for side in ('r','l')}|{f'HUNT_AF_{i}_{side}' for i in range(6) for side in ('r','l')}
    if not required.issubset(names) or not all(hasattr(spinal,k) for k in ('apply_weights','sync','reflex_gain','mn_edges','original_g')):raise ValueError('Missing actual Hunt source cells/runtime buffers')
    return {'class':str(expected),'neurons':len(names),'W2L_derived_MRO_verified':True}


def eligible_load_sources(spinal,family):
    """Exact policy keys; Hunt includes read-only inactive homologous input contracts."""
    G=utilities();c=G.raw_compiled(spinal);names=G.population_names(spinal)
    if family=='Hunt_PT':hunt_identity(spinal)
    ids=[i for i,n in enumerate(names) if n.startswith('IBEXC_')]
    return sorted({names[j] for i in ids for j in np.flatnonzero(c.g_max_non[i]) if names[j].startswith('Ib_')})


def hunt_weight_write_contract(spinal,gate_ids):
    """Audit actual apply_weights index sets,not merely their current zero values."""
    method=type(spinal).apply_weights
    if spinal.apply_weights.__func__ is not method or spinal.sync.__func__ is not type(spinal).sync:raise ValueError('Unaudited Hunt runtime method override')
    source=textwrap.dedent(inspect.getsource(method));tree=ast.parse(source);fields=set();masks=set()
    for node in ast.walk(tree):
        targets=node.targets if isinstance(node,ast.Assign) else [node.target] if isinstance(node,ast.AugAssign) else []
        for target in targets:
            if isinstance(target,ast.Subscript) and isinstance(target.value,ast.Attribute) and isinstance(target.value.value,ast.Name) and target.value.value.id=='c':
                fields.add(target.value.attr)
                if target.value.attr!='gmax' or not isinstance(target.slice,ast.Attribute) or not isinstance(target.slice.value,ast.Name) or target.slice.value.id!='self':raise ValueError('Unrecognized Hunt array write; explicit re-audit required')
                masks.add(target.slice.attr)
    c=spinal.compiled;written=set()
    for name in masks:written.update(np.asarray(getattr(spinal,name),dtype=int).tolist())
    incoming=set(np.flatnonzero(np.isin(c.post,gate_ids)).tolist());intersection=sorted(written&incoming)
    if intersection:raise ValueError('Hunt weight updates write incoming coincidence conductances')
    outputs=np.flatnonzero(np.isin(c.pre,gate_ids));rewritten=[i for i in outputs if int(i) in written]
    if any(i not in spinal.reflex_edges or c.gmax[i]!=spinal.reflex_gain[i] for i in rewritten):raise ValueError('Gate motor weights lack stable reflex-gain authority')
    return {'apply_weights_source_sha256':hashlib.sha256(source.encode()).hexdigest(),'array_fields_written':sorted(fields),'runtime_index_masks':sorted(masks),'incoming_gate_write_intersection':intersection,'gate_window_writes':False,'rewritten_gate_motor_slots':len(rewritten),'rewritten_gate_motor_source':'Unchanged command-independent reflex_gain','PT_balance_external_writes_covered':False}


def gate_hash(c,rows):
    h=hashlib.sha256();h.update(np.float64(c.dt).tobytes())
    for row in rows:
        i=row['index']
        for key in ('c_m','g_m','V_rest','i_b','time_factor_membrane'):h.update(np.float64(getattr(c,key)[i]).tobytes())
        for key in ('g_max_non','del_e','e_lo','e_hi'):
            a=getattr(c,key);h.update(a[i].tobytes());h.update(a[:,i].tobytes())
        for source in row['candidate']['loads']:
            j=source['index']
            for key in ('c_m','g_m','V_rest','i_b','time_factor_membrane'):h.update(np.float64(getattr(c,key)[j]).tobytes())
            h.update(c.input_connectivity[j].tobytes())
    return h.hexdigest()


def snapshot_gate(c,names,i,load_contract):
    sources=np.flatnonzero(c.g_max_non[i]);loads=[j for j in sources if names[j].startswith('Ib_')];phase=[j for j in sources if names[j].startswith('RG_E_')]
    outputs=np.flatnonzero(c.g_max_non[:,i]);out=[]
    for j in outputs:out.append(dict(name=names[j],index=int(j),g=float(c.g_max_non[j,i]),E=float(c.del_e[j,i]),lo=float(c.e_lo[j,i]),hi=float(c.e_hi[j,i]),C=float(c.c_m[j]),G=float(c.g_m[j])))
    incoming=[]
    for j in loads:
        incoming.append(dict(name=names[j],index=int(j),g=float(c.g_max_non[i,j]),E=float(c.del_e[i,j]),lo=float(c.e_lo[i,j]),hi=float(c.e_hi[i,j]),C=float(c.c_m[j]),G=float(c.g_m[j]),rest=float(c.V_rest[j]),bias=float(c.i_b[j]),**load_contract[names[j]]))
    p=None if not phase else dict(name=names[phase[0]],index=int(phase[0]),g=float(c.g_max_non[i,phase[0]]),E=float(c.del_e[i,phase[0]]),lo=float(c.e_lo[i,phase[0]]),hi=float(c.e_hi[i,phase[0]]))
    return dict(name=names[i],index=int(i),dt_s=float(c.dt),C=float(c.c_m[i]),G=float(c.g_m[i]),rest=float(c.V_rest[i]),bias=float(c.i_b[i]),loads=incoming,phase=p,outputs=out)


def apply_coincidence(spinal,policy=None):
    if policy is None:return dict(applied=False)
    G=utilities();G.assert_campaign_allowed('optional_profile_implementation')
    p=json.loads(Path(policy).read_text()) if isinstance(policy,(str,Path)) else copy.deepcopy(policy)
    G.exact_keys(p,('schema','family','phase_budget_over_leak','load_budget_over_leak','lower_gap_fraction','upper_gap_fraction','load_high_current_nA','load_contract_origin'),'coincidence policy')
    if p['schema']!='aarl_ibexc_binary_coincidence_v1' or p['family'] not in ('s3k','w2lvar','syn6','Hunt_PT'):raise ValueError('Unsupported coincidence policy')
    actual_module=type(spinal).__module__.split('.')[-1]
    hunt=p['family']=='Hunt_PT';identity=hunt_identity(spinal) if hunt else None
    if not hunt and {'build_network':'s3k','build_network_w2lvar':'w2lvar','build_network_syn6':'syn6'}.get(actual_module)!=p['family']:raise ValueError('Coincidence policy family does not match actual builder')
    values=[p[k] for k in ('phase_budget_over_leak','load_budget_over_leak','lower_gap_fraction','upper_gap_fraction')]
    if not np.isfinite(values).all() or min(values[:2])<=0 or not 0<values[2]<values[3]<1 or not p['load_contract_origin'].strip():raise ValueError('Invalid engineering budget/window contract')
    wrapper=spinal.compiled;c=G.raw_compiled(spinal)
    if wrapper is not c and not hunt:raise ValueError('Change coincidence topology before constructing accelerator')
    if c.spiking or c.electrical or c.electrical_rectified:raise ValueError('Graded chemical gate only')
    names=G.population_names(spinal);all_ids=[i for i,n in enumerate(names) if n.startswith('IBEXC_')];ids=[i for i in all_ids if not hunt or np.any(c.g_max_non[:,i])];inactive=[names[i] for i in all_ids if i not in ids]
    if not ids:raise ValueError('No actual IBEXC neurons')
    asymmetric=[]
    if hunt:
        if len(all_ids)!=10:raise ValueError('Expected ten allocated Hunt IBEXC source cells')
        for i in ids:
            other=names[i][:-2]+('_l' if names[i].endswith('_r') else '_r')
            if names.index(other) not in ids:asymmetric.append({'active_gate':names[i],'masked_homologue':other,'action':'Preserve existing mask; no bilateral completeness claim'})
    source_names={names[j] for i in (all_ids if hunt else ids) for j in np.flatnonzero(c.g_max_non[i]) if names[j].startswith('Ib_')}
    G.exact_keys(p['load_high_current_nA'],source_names,'per-source load-high current contracts')
    if any(not np.isfinite(v) or v<=0 for v in p['load_high_current_nA'].values()):raise ValueError('Load-high current must be explicit,finite,andpositive')
    for name,current in p['load_high_current_nA'].items():
        other=name[:-2]+('_l' if name.endswith('_r') else '_r')
        if other not in p['load_high_current_nA'] or p['load_high_current_nA'][other]!=current:raise ValueError('Shared bilateral load-current contracts required')
    write_contract=hunt_weight_write_contract(spinal,ids) if hunt else None
    candidate={k:getattr(c,k).copy() for k in ('g_max_non','del_e','e_lo','e_hi')};rows=[];new_edges=[];source_contract={}
    for source in sorted(source_names):
        j=names.index(source)
        if np.any(c.g_max_non[j]) or c.V_rest[j]!=0 or c.i_b[j]!=0:raise ValueError('Load-current adapter requires isolated zero-rest encoder')
        ports=np.flatnonzero(c.input_connectivity[j])
        if len(ports)!=1 or spinal.inputs[ports[0]]!=source:raise ValueError('Ambiguous physical-current input port')
        weight=float(c.input_connectivity[j,ports[0]]);current=float(p['load_high_current_nA'][source]);voltage=current*weight/c.g_m[j]
        if not np.isfinite(voltage) or voltage<=0:raise ValueError('Invalid load current-to-voltage contract')
        source_contract[source]=dict(high_current_nA=current,input_port_weight=weight,high_voltage_mV=float(voltage),physical_meaning=p['load_contract_origin'],applied_to_active_gate=bool(np.any(c.g_max_non[ids,j])))
    for i in ids:
        name=names[i];side=name[-1];phase=names.index('RG_E_'+side);incoming=np.flatnonzero(c.g_max_non[i]);loads=[j for j in incoming if names[j].startswith('Ib_')]
        unsupported=[names[j] for j in incoming if j not in loads and j!=phase]
        if unsupported or not loads or c.i_b[i]!=0 or c.V_rest[i]!=0 or np.any(c.input_connectivity[i]) or spinal.net.populations[i]['type'].params['gated'] or (c.gated and np.any(c.g_ion[:,i])) or c.g_m[i]<=0 or c.c_m[i]<=0:raise ValueError(f'Unsupported IBEXC inputs/bias: {name}/{unsupported}')
        group=name[len('IBEXC_'):-2]
        for j in loads:
            act=names[j][3:]
            if not act.endswith('_'+side) or spinal.muscles[act].groups[0]!=group:raise ValueError('Unverified group/side load binding')
        if c.g_max_non[i,phase]:
            E=float(c.del_e[i,phase]);lo=float(c.e_lo[i,phase]);hi=float(c.e_hi[i,phase])
        else:
            if p['family']!='syn6' or not hasattr(spinal,'canonical_right_ledger'):raise ValueError('Only explicit canonical-right syn6 may add a phase edge')
            destinations=np.flatnonzero((c.g_max_non[:,phase]>0)&(c.del_e[:,phase]>0));domains={(float(c.del_e[q,phase]),float(c.e_lo[q,phase]),float(c.e_hi[q,phase])) for q in destinations}
            if len(domains)!=1:raise ValueError('Ambiguous actual RG_E excitatory reversal/window')
            E,lo,hi=domains.pop();new_edges.append(dict(source=names[phase],target=name,reason='Explicit engineering phase predicate absent from original syn6'))
        if E<=0 or hi<=lo or np.any(c.del_e[i,loads]<=0):raise ValueError('Excitatory phase and load routes required')
        gm=float(c.g_m[i]);gp=p['phase_budget_over_leak']*gm;gl=p['load_budget_over_leak']*gm;weights=c.g_max_non[i,loads];weights=weights/weights.sum()
        if not np.isfinite([gm,gp,gl,c.c_m[i],E,lo,hi]).all() or not np.isfinite(weights).all() or c.dt*(gm+gp+gl)/c.c_m[i]>1:raise ValueError('Invalid or nonmonotone discrete coincidence gate')
        candidate['g_max_non'][i,loads]=gl*weights;candidate['g_max_non'][i,phase]=gp;candidate['del_e'][i,phase]=E;candidate['e_lo'][i,phase]=lo;candidate['e_hi'][i,phase]=hi
        for j in loads:candidate['e_lo'][i,j]=0.;candidate['e_hi'][i,j]=source_contract[names[j]]['high_voltage_mV']
        alone_phase=gp*E/(gm+gp);alone_load_bound=gl*float(max(c.del_e[i,loads]))/(gm+gl);both=(gp*E+gl*float(weights@c.del_e[i,loads]))/(gm+gp+gl);alone=max(alone_phase,alone_load_bound);gap=both-alone
        if gap<=1e-6:raise ValueError('No feasible strict coincidence window for actual reversals')
        lower=alone+p['lower_gap_fraction']*gap;upper=alone+p['upper_gap_fraction']*gap;outputs=np.flatnonzero(c.g_max_non[:,i])
        if not len(outputs):raise ValueError('Unused IBEXC output')
        candidate['e_lo'][outputs,i]=lower;candidate['e_hi'][outputs,i]=upper
        rows.append(dict(name=name,index=i,phase_source=names[phase],phase_g_uS=gp,maximum_load_pool_g_uS=gl,allocation='Functional total conductance budget with preserved relative old Ib weights; not blind division by fan-in',alone_phase_mV=alone_phase,alone_load_upper_bound_mV=alone_load_bound,both_on_mV=both,output_lo_mV=lower,output_hi_mV=upper,strict_off_margin_mV=lower-alone,strict_on_margin_mV=both-upper,original=snapshot_gate(c,names,i,source_contract)))
    other_routes=[]
    for source in sorted(source_names):
        src=names.index(source)
        for dst in np.flatnonzero(c.g_max_non[:,src]):
            if dst in ids:continue
            preserved=all(candidate[k][dst,src]==getattr(c,k)[dst,src] for k in candidate)
            if not preserved:raise AssertionError('Non-IBEXC Ib route changed')
            other_routes.append(dict(source=source,target=names[dst],g_uS=float(c.g_max_non[dst,src]),reversal_mV=float(c.del_e[dst,src]),voltage_window_mV=[float(c.e_lo[dst,src]),float(c.e_hi[dst,src])],preserved=True,central_or_load_bypass=names[dst].startswith(('RG_','PF_','InE','InF','LBIN'))))
    states={k:getattr(c,k).copy() for k in ('V','V_last','V_0','b_gate','c_gate','b_gate_last','c_gate_last') if hasattr(c,k)}
    if hunt and not np.array_equal(candidate['g_max_non']!=0,c.g_max_non!=0):raise ValueError('Hunt coincidence may not reactivate masked paths or add topology')
    for key,value in candidate.items():getattr(c,key)[:]=value
    if wrapper is not c:
        # Existing sparse slots only; never insert a masked/template edge.
        wrapper.gmax[:]=c.g_max_non[wrapper.post,wrapper.pre];wrapper.reversal[:]=c.del_e[wrapper.post,wrapper.pre];wrapper.elo[:]=c.e_lo[wrapper.post,wrapper.pre];wrapper.erange[:]=(c.e_hi-c.e_lo)[wrapper.post,wrapper.pre]
    for row in rows:row['candidate']=snapshot_gate(c,names,row['index'],source_contract)
    if any(not np.array_equal(getattr(c,k),v) for k,v in states.items()):raise AssertionError('Gate policy changed dynamic state')
    ledger=dict(schema=p['schema'],policy=p,rows=rows,added_phase_edges=new_edges,new_neurons=0,Hunt_identity=identity,Hunt_weight_write_contract=write_contract,inactive_gates_preserved=inactive,asymmetric_active_gate_masks_preserved=asymmetric,source_current_contracts=source_contract,unchanged_other_Ib_routes=other_routes,central_or_load_bypass_count=sum(v['central_or_load_bypass'] for v in other_routes),protected_gate_sha256=gate_hash(c,rows),state_unchanged=True,scope='Binary thresholded coincidence,not multiplication. Load-high is explicit F/F0-to-current encoding,not a GRF sensor or absolute physical maximum.',Hunt_runtime_hook_installed=False,persistence_after_Hunt_dynamic_apply_verified=False,training_allowed=False)
    spinal.ibexc_coincidence_ledger=ledger;verify_coincidence(spinal)
    return ledger


def verify_coincidence(spinal):
    G=utilities();c=G.raw_compiled(spinal);ledger=spinal.ibexc_coincidence_ledger
    if gate_hash(c,ledger['rows'])!=ledger['protected_gate_sha256']:raise ValueError('Coincidence fields overwritten; no Hunt/dynamic-update persistence claim')
    wrapper=spinal.compiled
    if wrapper is not c:
        if not np.array_equal(wrapper.gmax,c.g_max_non[wrapper.post,wrapper.pre]) or not np.array_equal(wrapper.reversal,c.del_e[wrapper.post,wrapper.pre]) or not np.array_equal(wrapper.erange,c.e_hi[wrapper.post,wrapper.pre]-c.e_lo[wrapper.post,wrapper.pre]) or not np.array_equal(wrapper.elo,c.e_lo[wrapper.post,wrapper.pre]):raise ValueError('Coincidence accelerator cache mismatch')
    return dict(actual_gate_fields_match=True,protected_gate_sha256=ledger['protected_gate_sha256'])


def isolated(circuit,phase_on,load_ratio,release=False):
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    from sns_toolbox.connections import NonSpikingSynapse
    dt=circuit['dt_s'];n=Network(name='actual copied IBEXC circuit');loads=circuit['loads'];count=len(loads);phase=count;gate=count+1;post=count+2
    for j,src in enumerate(loads):n.add_neuron(NonSpikingNeuron(membrane_capacitance=src['C'],membrane_conductance=src['G']),name='load'+str(j));n.add_input('load'+str(j))
    n.add_neuron(NonSpikingNeuron(),name='phase');n.add_neuron(NonSpikingNeuron(membrane_capacitance=circuit['C'],membrane_conductance=circuit['G']),name='gate');output=circuit['outputs'][0];n.add_neuron(NonSpikingNeuron(membrane_capacitance=output['C'],membrane_conductance=output['G']),name='post')
    for j,src in enumerate(loads):n.add_connection(NonSpikingSynapse(max_conductance=src['g'],reversal_potential=src['E'],e_lo=src['lo'],e_hi=src['hi']),'load'+str(j),'gate')
    source=circuit['phase']
    if source is not None:n.add_connection(NonSpikingSynapse(max_conductance=source['g'],reversal_potential=source['E'],e_lo=source['lo'],e_hi=source['hi']),'phase','gate')
    n.add_connection(NonSpikingSynapse(max_conductance=output['g'],reversal_potential=output['E'],e_lo=output['lo'],e_hi=output['hi']),'gate','post');c=n.compile(dt=dt,backend='numpy');trace=[]
    warm_steps=round(.5/dt)
    for step in range(round((1. if release else .8)/dt)):
        phase_value=phase_on and (not release or step<warm_steps);off=release and step>=warm_steps
        # Held force with larger phase-off encoding gain is an explicit stress test.
        x=np.array([s['high_current_nA']*s['input_port_weight']*load_ratio*(10/3 if off else 1) for s in loads]);c.V[phase]=0. if not phase_value else (source['hi'] if source is not None else 5.)
        c.forward(x);fraction=float(np.clip((c.V[gate]-output['lo'])/(output['hi']-output['lo']),0,1));old_fraction=float(np.clip((c.V_last[gate]-output['lo'])/(output['hi']-output['lo']),0,1));current=output['g']*old_fraction*(output['E']-c.V_last[post]);trace.append([float((step+1)*dt),float(c.V[gate]),fraction,float(c.V[post]),float(current)])
    a=np.array(trace);offmask=a[:,0]>.5
    return dict(phase_on=phase_on,force_over_F0=load_ratio,final_gate_mV=float(a[-1,1]),final_fraction=float(a[-1,2]),final_post_mV=float(a[-1,3]),gate_release_s=float(a[np.flatnonzero(offmask&(a[:,2]==0))[0],0]-.5) if release and (offmask&(a[:,2]==0)).any() else None,synaptic_current_release_s=float(a[np.flatnonzero(offmask&(a[:,4]==0))[0],0]-.5) if release and (offmask&(a[:,4]==0)).any() else None,sparse_trace=a[::10].tolist())


def probe(family,model,out):
    sys.path[:0]=[str(HERE.parent),str(HERE.parent/'w2l_cpg'),str(HERE.parent.parent/'robot_torso'/'bpa_neural_coordination')];os.environ['AARL_NET']='' if family=='s3k' else family
    import mujoco,build_network as B,params as P
    from graded_profile import default_profile,apply_profile,verify_profile,field_digest,population_names
    from mirror_canonical import apply_canonical_right,default_policy,verify_canonical
    from sparse_pairwise_eval import SparsePairwiseForward
    import gait_standing_family as F
    m=mujoco.MjModel.from_binary_path(str(model));net=B.build([m.actuator(j).name for j in range(m.nu)])
    if family=='syn6':apply_canonical_right(net,default_policy())
    apply_profile(net,default_profile(family),family=family)
    names=population_names(net);c=net.compiled;ids=[i for i,n in enumerate(names) if n.startswith('IBEXC_')];sources={names[j] for i in ids for j in np.flatnonzero(c.g_max_non[i]) if names[j].startswith('Ib_')}
    # Actual canonical runner settings,declared stance=1/F=F0. No saved-winner binding.
    current=.5*float(P.AFF['ib_gain'])*.3;policy=make_policy(family,{n:current for n in sorted(sources)},f'Canonical runner assumption: uIb=.5*g_ib_eff*F/F0; AFF.ib_gain={P.AFF["ib_gain"]},stance=1,g_ib_eff=.3*AFF.ib_gain; F/F0=1. Not a saved winner/Hunt binding.')
    negative_fixtures={};before=field_digest(c);i=ids[-1]
    for label in ('gate_bias','unsupported_drive_input'):
        if label=='gate_bias':old=c.i_b[i];c.i_b[i]=old+.01
        else:j=names.index('DRIVE');old=c.g_max_non[i,j];c.g_max_non[i,j]=.1
        try:apply_coincidence(net,policy);rejected=False
        except ValueError:rejected=True
        if label=='gate_bias':c.i_b[i]=old
        else:c.g_max_non[i,j]=old
        if not rejected or field_digest(c)!=before:raise AssertionError('Negative fixture accepted or partially mutated network')
        negative_fixtures[label]=True
    ledger=apply_coincidence(net,policy);truth=[]
    representative=max(ledger['rows'],key=lambda r:len(r['candidate']['loads']))
    for variant in ('original','candidate'):
        for phase,ratio in ((False,0),(True,0),(False,1),(True,1)):
            value=isolated(representative[variant],phase,ratio);value.update(variant=variant)
            if variant=='candidate' and value['final_fraction']!=float(phase and ratio):raise AssertionError('Normal-F0 coincidence truth table failed')
            truth.append(value)
    releases=[isolated(representative['candidate'],True,r,True) for r in (1.,1.5,5.)]
    load_sweep=[isolated(representative['candidate'],True,r) for r in (0.,.1,.25,.5,.75,1.,1.5,5.)]
    if any(r['gate_release_s'] is None or r['synaptic_current_release_s'] is None for r in releases):raise AssertionError('Held-load phase-off release failed')
    # Test every actual gate analytically acrossnormal andsupra-normal force encodings.
    for row in ledger['rows']:
        if not row['alone_load_upper_bound_mV']<row['output_lo_mV']<row['output_hi_mV']<row['both_on_mV']:raise AssertionError('Gate margins invalid')
    i=ledger['rows'][0]['index'];j=names.index(ledger['rows'][0]['phase_source']);value=c.g_max_non[i,j];c.g_max_non[i,j]=value*1.01
    try:verify_coincidence(net);rejected=False
    except ValueError:rejected=True
    c.g_max_non[i,j]=value
    if not rejected:raise AssertionError('Dynamic overwrite not detected')
    net.compiled=SparsePairwiseForward(c);apply_profile(net,default_profile(family),family=family);verify_profile(net);verify_coincidence(net)
    if F.tie_weights(net):raise AssertionError('Final historical mirror projection changes candidate')
    verify_profile(net);verify_coincidence(net)
    if family=='syn6':verify_canonical(net)
    a=copy.deepcopy(c);fast=SparsePairwiseForward(copy.deepcopy(c));rng=np.random.default_rng(1105);h=hashlib.sha256()
    for step in range(128):
        x=rng.uniform(0,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;a.forward(x);fast.forward(x)
        if not np.array_equal(a.V,fast.V):raise AssertionError('Coincidence dense/sparse mismatch')
        h.update(a.V.tobytes())
    result=dict(source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),family=family,policy=policy,ledger=ledger,representative=representative['name'],truth_table=truth,held_load_release=releases,phase_on_analog_load_sweep=load_sweep,negative_fixtures_rejected_without_partial_mutation=negative_fixtures,final_consumption=verify_coincidence(net),FSA_final_consumption=verify_profile(net),final_historical_mirror_projection_changes=0,dense_accelerated128_exact=True,neural_voltage_sha256=h.hexdigest(),overwrite_rejected=rejected,normal_F0_current_nA=current,all_gate_static_margins_verified=True,Hunt_dynamic_hook=False,body_trials=0,training_allowed=False)
    out.mkdir(parents=True,exist_ok=True);(out/(family+'_COINCIDENCE.json')).write_text(json.dumps(result,indent=2));(out/(family+'_POLICY.json')).write_text(json.dumps(policy,indent=2));print(json.dumps(dict(family=family,gates=len(ids),new_phase_edges=len(ledger['added_phase_edges']),load_current_nA=current,gate_release_s=[r['gate_release_s'] for r in releases],fast_exact=True)))


def audit_hunt_gain_pairing(n):
    """Read-only lexical/declared pairing audit; never infer new axis reflections."""
    c=n.compiled;names=utilities().population_names(n);lookup={(names[p],names[q]):i for i,(p,q) in enumerate(zip(c.pre,c.post))}
    if len(lookup)!=len(c.pre):raise ValueError('Duplicate sparse endpoints require separate pairing audit')
    def mirror(name):return __import__('re').sub(r'_([rl])(?=$|[^a-zA-Z0-9])',lambda m:'_l' if m[1]=='r' else '_r',name)
    geometry=set(n.geom_edges.tolist());explicit={}
    for a,b in n.midline_drive_pairs:explicit[a]=b;explicit[b]=a
    masks={'pf_gain':np.union1d(n.pfbank_edges,n.pf_internal_edges),'ff_gain':np.union1d(n.pf_edges,n.geom_edges),
           'reflex_gain':np.union1d(n.reflex_edges,n.reciprocal_edges) if n.options.get('learn_reciprocal',False) else n.reflex_edges,'bi_gain':n.bi_edges}
    result={}
    for buffer,indices in masks.items():
        values=getattr(n,buffer);selected=set(np.asarray(indices,dtype=int).tolist());seen=set();pairs=[];unpaired=[];shared=[]
        for i in sorted(selected):
            if i in seen:continue
            pre,post=names[c.pre[i]],names[c.post[i]];endpoint=(pre,post);mapped=(mirror(pre),mirror(post));j=explicit.get(i,lookup.get(mapped));rule='source-declared midline channel reflection' if i in explicit else 'lexical R/L label reflection only'
            if mapped==endpoint and i not in explicit:shared.append({'source':pre,'target':post,'value':float(values[i])});seen.add(i);continue
            if j is None or j not in selected:
                unpaired.append({'source':pre,'target':post,'value':float(values[i]),'nonzero':bool(values[i]!=0),'above_pairing_tolerance':bool(abs(values[i])>1e-12),'scope':'shared/midline' if '_mid' in pre+post else 'missing counterpart or counterpart outside authoritative mask'});seen.add(i);continue
            seen.update((i,j));delta=float(abs(values[i]-values[j]));scope='geometry-recruitment; lexical only,anatomical direction not verified' if i in geometry or j in geometry else 'declared midline' if i in explicit else 'other labeled bilateral paths'
            pairs.append({'source':pre,'target':post,'mirror_source':names[c.pre[j]],'mirror_target':names[c.post[j]],'value':float(values[i]),'mirror_value':float(values[j]),'absolute_mismatch':delta,'rule':rule,'scope':scope})
        active_pair=[p for p in pairs if p['value']!=0 or p['mirror_value']!=0];bad=[p for p in pairs if p['absolute_mismatch']>1e-12]
        result[buffer]={'authoritative_slots':len(selected),'authoritative_nonzero_slots':int(np.count_nonzero(values[list(selected)])),'all_buffer_nonzero_but_not_authoritative':int(np.count_nonzero(values))-int(np.count_nonzero(values[list(selected)])),
                        'paired_parameters':len(pairs),'nonzero_paired_parameters':len(active_pair),'comparison_tolerance_uS':1e-12,'mismatched_pairs':len(bad),'maximum_absolute_mismatch_uS':max((p['absolute_mismatch'] for p in pairs),default=0.),
                        'full_buffer_sha256':hashlib.sha256(values.tobytes()).hexdigest(),'authoritative_mask_sha256':hashlib.sha256(np.array(sorted(selected),dtype=np.int64).tobytes()).hexdigest(),
                        'largest_mismatches':sorted(bad,key=lambda p:p['absolute_mismatch'],reverse=True)[:12],
                        'pair_scope_counts':{scope:sum(p['scope']==scope for p in pairs) for scope in sorted({p['scope'] for p in pairs})},
                        'unpaired_count':len(unpaired),'unpaired_literal_nonzero_count':sum(p['nonzero'] for p in unpaired),'unpaired_above_tolerance_count':sum(p['above_pairing_tolerance'] for p in unpaired),'maximum_unpaired_absolute_gain_uS':max((abs(p['value']) for p in unpaired),default=0.),'unpaired':unpaired,'shared_self_edges_count':len(shared),'shared_self_edges':shared,
                        'interpretation':'Exact named-buffer comparison. Geometry-dependent physical signs/authority and PT calibration are outside this pairing proof; no values changed.'}
    return {'buffers':result,'require_bilateral_symmetry_option':bool(n.options.get('require_bilateral_symmetry',False)),'declared_midline_drive':bool(n.options.get('bilateral_midline_drive',False)),
            'declared_midline_edge_pairs':len(n.midline_drive_pairs),'geometry_slots':len(n.geom_edges),'geometry_factor_range_unbound':[float(n.geometry_factor.min()),float(n.geometry_factor.max())],
            'geometry_scope':'Only static authoritative gains audited. Live geometry_factor is unbound; no physical mirroring or moment-arm/coordinate sign repair inferred.',
            'PT_scope':'Unbound mechanical calibration explicitly excluded. PT gain buffers are not one of these four fitted parameter sets.',
            'RF_HAM_fourcell_labels':{side:[name for name in (f'PF_BI_JOINT_HAM_{side}',f'PF_BI_JOINT_RECT_{side}',f'PF_IN_BI_PAIR_HAM_{side}',f'PF_IN_BI_PAIR_RECT_{side}') if name in names] for side in ('r','l')},'no_parameter_repair':True}


def probe_hunt(condition,out,timing_mode='canonical_right_role_budget'):
    sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
    from hunt_pathway_audit import C,MODEL,options_from_physical_condition,mujoco
    import hunt_graded_profile as local,hunt_timing_profile as timing
    from graded_profile import default_profile,population_names
    cfg=json.loads(condition.read_text());weight_dir=Path(cfg.get('weights_directory',str(HERE.parent/'hunt_training/runs/fitted')))
    dependencies=(Path(__file__),Path(C.__file__),Path(local.__file__),Path(timing.__file__),condition)+tuple(weight_dir/name for name in ('pf_weights.json','mn_weights.json','reflex_weights.json','bi_weights.json','timing_weights.json','circuit_options.json') if (weight_dir/name).exists())
    frozen={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}
    m=mujoco.MjModel.from_xml_path(str(MODEL));n=C.build([m.actuator(j).name for j in range(m.nu)],dt=.002,options=options_from_physical_condition(cfg));had_PT=n.__dict__.pop('balance',None) is not None
    n.apply_weights(1.);source_pairing=audit_hunt_gain_pairing(n);local.install(n,default_profile('w2lvar'));timing_profile=timing.default_profile(mode=timing_mode,relay_timing='fast_role_budget');timing.install(n,timing_profile)
    sources=eligible_load_sources(n,'Hunt_PT');policy=make_policy('Hunt_PT',{name:.225 for name in sources},'Explicit canonical runner engineering contract: AFF.ib_gain=1.5,full-stance attenuation0.3,multiplier0.5,F/Fmax=1. Actual saved-condition force/gain encoding is not body-validated.')
    c=n.compiled.c;before=c.g_max_non.copy();ledger=apply_coincidence(n,policy)
    if not np.array_equal(before!=0,c.g_max_non!=0):raise AssertionError('Hunt masked topology changed')
    checks=[]
    for command in (0.,.5,1.,0.):
        n.apply_weights(command);n.sync();checks.append({'command':command,'gate':verify_coincidence(n),'local':local.verify_runtime(n),'timing':timing.verify_runtime(n)})
    u=np.zeros(len(n.inputs))
    for step in range(32):n.step(u);verify_coincidence(n);local.verify_runtime(n);timing.verify_runtime(n)
    ledger['persistence_after_Hunt_dynamic_apply_verified']=True
    if frozen!={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}:raise RuntimeError('Composition source changed during audit')
    result={'schema':'aarl_Hunt_PT_IBEXC_composition_v1','source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha256':frozen,'condition':str(condition),'condition_sha256':hashlib.sha256(condition.read_bytes()).hexdigest(),'neurons':len(n.compiled.V),'ledger':ledger,'source_gain_pairing_before_profiles':source_pairing,'composed_gain_pairing':audit_hunt_gain_pairing(n),'timing_profile':timing_profile,'command_checks':checks,'neural_steps':32,'finite':bool(np.isfinite(n.compiled.V).all()),'unbound_PT_mechanical_pathway_excluded':had_PT,'masked_path_reactivation':False,'actual_body_force_encoding_bound':False,'source_current_contract_nA':.225,'Hunt_training':False,'body_trials':0,'global_release':False}
    if not result['finite']:raise AssertionError('Hunt neural composition became nonfinite')
    out.mkdir(parents=True,exist_ok=True);(out/'HUNT_COINCIDENCE.json').write_text(json.dumps(result,indent=2));print(json.dumps({'neurons':len(n.compiled.V),'active_gates':len(ledger['rows']),'inactive_gates_preserved':ledger['inactive_gates_preserved'],'command_persistence':True,'steps32_finite':True,'output':str(out)}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--family',choices=('s3k','w2lvar','syn6','Hunt_PT'),required=True);p.add_argument('--model',type=Path);p.add_argument('--condition',type=Path);p.add_argument('--hunt-timing-mode',choices=('function_equivalent','saturated_source_equivalent','canonical_right_role_budget'),default='canonical_right_role_budget');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.family=='Hunt_PT':
        if a.condition is None:p.error('Hunt_PT requires --condition')
        probe_hunt(a.condition.resolve(),a.out,a.hunt_timing_mode)
    else:
        if a.model is None:p.error('Canonical families require --model')
        probe(a.family,a.model,a.out)
