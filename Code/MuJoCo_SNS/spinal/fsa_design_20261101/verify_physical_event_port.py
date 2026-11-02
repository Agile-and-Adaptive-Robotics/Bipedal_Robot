"""Physical-domain-v2 electrical checks/benchmarks, never body integration."""
import os,sys,copy,json,math,time,argparse,hashlib
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).parent;sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import numpy as np
from fsa_design_20261101.fast_event_port import FastEventConductanceInputPort as Fast,physical_domain_profile,physical_population_profile,canonical_typed_profile,canonical_tight_tail_profile,install_physical_port,event_conductance_tick
from fsa_design_20261101.event_conductance_port import EventConductanceInputPort,steady_carrier_fraction
from fsa_design_20261101.graded_profile import default_profile,apply_profile,population_names
from fsa_design_20261101.table2_glif import constant_threshold_rate_hz

def write(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
def array_error(a,b):return max([float(np.max(abs(a[k]-b[k]))) for k in a if a[k].dtype.kind in 'fc' and a[k].size]+[0.])
def main(family,output,population=False,typed=False,tight=False):
    os.environ['AARL_NET']='' if family=='s3k' else family;out=Path(output);out.mkdir(parents=True,exist_ok=True)
    import mujoco,build_network as B
    files=[Path(__file__),HERE/'fast_event_port.py',HERE/'event_conductance_port.py',HERE/'graded_profile.py']+[HERE.parent/q for q in ('build_network.py','build_network_syn6.py','build_network_w2lvar.py','params.py','synergy_basis.npz','fsa_results/fsa_backsolve.npz')];hashes={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
    pin=json.loads(Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\gait_table2\source_profiles\syn6_DESIGN_PROFILE.json').read_text())['source_NMF_basis']['sha256'];assert hashes[str(HERE.parent/'synergy_basis.npz')]==pin
    m=mujoco.MjModel.from_xml_path(str(HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'));n=B.build([m.actuator(i).name for i in range(m.nu)]);grade=default_profile(family);apply_profile(n,grade,family=family);raw=n.compiled;names=population_names(n);p=(canonical_tight_tail_profile if tight else (canonical_typed_profile if typed else (physical_population_profile if population else physical_domain_profile)))(family);args={k:v for k,v in p.items() if k not in ('schema','family','binding','capture_events')}
    protocol=dict(profile=p,graded_profile=grade,scope='Actual compiled neural family; static MuJoCo compile supplies actuator names only. No body, optimization, Hunt or curriculum.',gates=dict(scalar_JIT_state_max=1e-9,Python_reference_state_max=1e-9,event_time_max_s=1e-11,complete_state32_exact=True,nonlinear_short_window_rejected=True,physical_rate_error_max_Hz=1.,steady_fullscale_G_error_max=.013,normal_Ib_receiver_modulation_relative_RMS_max=.05,normal_Ib_receiver_extra_phase_deg_max=10.),criterion_provenance='The modulation-relative receiver criterion was declared to root before the typed/N128 extension, after the5nA/N16 raw-shot-noise issue was observed. It does not rewrite the older fullscale-conductance gate.',latency_policy='Measure from exact steady-quiet membrane state, not a precharged phase initializer. Report normalIb delay even if it exceeds10deg at7.156Hz; no threshold relaxation or fabricated spikes.',benchmark='Neural-only. JIT warmup excluded. Same bounded-logging scalar/JIT algorithms and same network; no full-body speed claim.',source_files_sha256=hashes,training_allowed=False)
    write(out/'PENDING_PROTOCOL.json',protocol) # Failed repeats never replace the last completed protocol.
    def make(**kw):return Fast(raw,names,**dict(args,**kw))
    z=np.zeros(len(n.inputs));warm=make();t=time.perf_counter();warm.forward(z);warmup=time.perf_counter()-t
    a=install_physical_port(n,p,family=family);b=make(backend='reference');ref=make(backend='reference',capture_events=True) if population else EventConductanceInputPort(raw,names,**args);max_state=max_python=max_events=0.;data=[];reference_event_cursor=0
    for k in range(600):
        x=z.copy();fraction=(0.,.2 if typed else .045,.002,1.,.584,.1)[k//100]
        for q,(port,w) in enumerate(a.ports):x[int(port)]=fraction*a.input_limits[q]/w
        n.step(x);b.forward(x);ref.forward(x)
        assert np.isfinite(a.V).all() and np.isfinite(b.V).all() and np.isfinite(ref.V).all()
        max_state=max(max_state,float(abs(a.V-b.V).max()),float(abs(a.U-b.U).max()),float(abs(a.gstate-b.gstate).max()));max_python=max(max_python,float(abs(a.V-ref.V).max()),float(abs(a.U-ref.U).max()),float(abs(a.gstate-ref.gstate).max()))
        block=np.array(ref.events[reference_event_cursor:]).reshape(-1,3);reference_event_cursor=len(ref.events);assert block.shape==a.last_events().shape
        if len(block):max_events=max(max_events,float(abs(block-a.last_events()).max()))
        assert np.array_equal(a.counts,b.counts) and np.array_equal(a.counts,ref.counts)
        data.append([a.time,*[x[int(port)]*w for port,w in a.ports],*a.counts,float(a.V[names.index('MN_soleus_r')])])
    assert max_state<1e-9 and max_python<1e-9 and max_events<1e-11 and not a.events
    state=a.export_state();u=z.copy();u[int(a.ports[0,0])]=.225
    for _ in range(32):a.forward(u)
    expected=a.export_state();a.import_state(state)
    for _ in range(32):a.forward(u)
    restored=a.export_state();assert all(np.array_equal(expected[k],restored[k]) for k in expected)
    # Route-specific0..0.225nA window is a different nonlinear axon function.
    bad=copy.deepcopy(raw);ii=names.index('Ib_soleus_r');post=int(np.flatnonzero(bad.g_max_non[:,ii])[0]);bad.e_hi[post,ii]=.225;rejected=False
    try:Fast(bad,names,**args)
    except ValueError:rejected=True
    assert rejected
    # Out-of-domain checks precede any encoder mutation.
    c=make();before=c.export_state();bad_u=z.copy();bad_u[int(c.ports[-1,0])]=c.input_limits[-1]*1.001
    try:c.forward(bad_u)
    except ValueError:pass
    else:raise AssertionError('Out-of-domain current accepted')
    after_invalid=c.export_state();mutated=[k for k in before if before[k].tobytes()!=after_invalid[k].tobytes()];assert not mutated,mutated # Scratch NaN bitpatterns must remain unchanged too.
    # Disabled path and actual-spike/output-path ablation remain distinct.
    disabled=make(enabled=False);plain=type(disabled.inner)(copy.deepcopy(disabled.inner.c));off_base=make(outgoing_off=True);off_pulse=make(outgoing_off=True);non=np.setdiff1d(np.arange(len(names)),a.indices);pulse_spikes=0
    for k in range(160):
        u=z.copy();u[int(a.ports[4,0])]=.225 if 50<=k<100 else 0.;plain.forward(u);disabled.forward(u);assert np.array_equal(plain.V,disabled.V)
        off_base.forward(z);off_pulse.forward(u);pulse_spikes+=off_pulse.counts[4];assert np.array_equal(off_base.V[non],off_pulse.V[non])
    assert pulse_spikes>0
    def component(c,u):
        c.event_block_count=event_conductance_tick(u,c.ports,c.input_limits,c.activation.range_mV,c.activation.theta0_mV,c.cell_bias,c.design.tau_mem_s,c.dt,c.time,c.U,c.counts,c.total_counts,c.cell_counts,c.event_times,c.event_block,c.edge_specs,c.peaks,c.gstate,c.gs,c.ge);c.time+=c.dt
    points=[];source=4;limit=a.input_limits[source]
    for current in sorted(set((0.,limit*.002,limit*.0025,limit*.004,.225,limit*.5,limit))):
        c=make();u=z.copy();u[int(c.ports[source,0])]=current;receiver=int(next(e[1] for e in c.edge_specs if int(e[0])==source and e[4]>0));mean=sum(e[4] for e in c.edge_specs if int(e[0])==source and int(e[1])==receiver);count=0;g=0.
        for k in range(1500):
            component(c,u)
            if k>=500:count+=c.counts[source];g+=c.gs[receiver]
        rate=count/(c.N*2.);inf=current/limit*20+c.cell_bias[source];rates=np.zeros(c.N);active=inf>c.activation.theta0_mV;rates[active]=1/(c.design.tau_mem_s*np.log(inf[active]/(inf[active]-c.activation.theta0_mV)));theory=float(rates.mean());normalized=g/(1000*mean);ek=next(k for k,e in enumerate(c.edge_specs) if int(e[0])==source and int(e[1])==receiver);tau=c.edge_specs[ek,3];factors=np.zeros(c.N);factors[active]=tau*rates[active]*(-np.expm1(-1/(tau*rates[active])));target=float(factors.mean()*c.peaks[ek]/mean);assert abs(rate-theory)<=1 and abs(normalized-current/limit)<=.013
        points.append(dict(current_nA=current,declared_type_max_nA=float(limit),actual_Hz=rate,exact_Hz=theory,normalized_mean_G=normalized,exact_steady_G=target,encoding_domain_G_error=normalized-current/limit))
    # Exact population DC curve and an analytic full-domain error envelope.
    c=make();ek=next(k for k,e in enumerate(c.edge_specs) if int(e[0])==source and e[4]>0);tau=c.edge_specs[ek,3];D=c.edge_specs[ek,4]/c.peaks[ek];theta=c.activation.theta0_mV;R=c.activation.range_mV;F=R/(c.design.tau_mem_s*theta);beta=c.cell_bias[source]/theta;eps=theta/(2*R);A=tau*(1-p['delta'])*F/D;B=tau*F/D;extra=(theta/R)*float(np.maximum(beta-.5,0).mean());dc_bounds=[min(-eps,A*(1-eps)-1),max(0.,B-1)+B*extra];curve=[]
    for current in np.linspace(0,limit,1001):
        inf=current/limit*R+c.cell_bias[source];active=inf>theta;f=np.zeros(c.N);f[active]=1/(c.design.tau_mem_s*np.log(inf[active]/(inf[active]-theta)));reset=np.zeros(c.N);reset[active]=tau*f[active]*(-np.expm1(-1/(tau*f[active])));curve.append(float(reset.mean()/D-current/limit))
    assert min(curve)>=dc_bounds[0]-1e-12 and max(curve)<=dc_bounds[1]+1e-12
    delays=[]
    for current in sorted(set((.225,limit*.5,limit))):
        c=make();c.U[:]=c.cell_bias;u=z.copy();u[int(c.ports[source,0])]=current;first=None
        for _ in range(100):
            component(c,u)
            if c.event_block_count and first is None:first=float(c.last_events()[:,0].min()) # Block is cell-ordered, not chronological.
        applied=current/c.input_limits[source]*c.activation.range_mV;eligible=applied+c.cell_bias[source]>c.activation.theta0_mV;expected=float(np.min(c.design.tau_mem_s*np.log(applied/(applied+c.cell_bias[source,eligible]-c.activation.theta0_mV))));assert first is not None and abs(first-expected)<1e-11;delays.append(dict(current_nA=current,first_actual_spike_s=first,independent_exact_first_spike_s=expected,first_host_output_endpoint_s=math.ceil(first/c.dt)*c.dt,equivalent_7p156Hz_phase_deg=first*(4/.559)*360,initialization='Exact quiet membrane equilibrium; no precharge'))
    bandwidth={};waveforms={};frequency=4/.559;omega=2*math.pi*frequency
    role_cases=[(4,'high_domain',float(a.input_limits[4])*.5,float(a.input_limits[4])*.4),(4,'normal_Ib',.225,.1),(0,'normal_Ia',1.,.5),(2,'normal_II',1.,.5)]
    for source,label,center,amplitude in role_cases:
        limit=a.input_limits[source]
        c=make();c.U[:]=c.cell_bias;candidates=[int(e[1]) for e in c.edge_specs if int(e[0])==source and e[4]>0 and not np.any(raw.g_ion[:,int(e[1])])];assert candidates,'No passive receiver; a separate ionic response contract is required';receiver=next((i for i in candidates if names[i].startswith('MN_')),candidates[0]);mean=sum(e[4] for e in c.edge_specs if int(e[0])==source and int(e[1])==receiver);wave=[];vs=vg=2.5;gm=raw.g_m[receiver];cm=raw.c_m[receiver];rest=raw.V_rest[receiver];baseline=gm*(2.5-rest);rev=float(next(e[2] for e in c.edge_specs if int(e[0])==source and int(e[1])==receiver))
        for k in range(1700):
            t=k*c.dt;u=z.copy();u[int(c.ports[source,0])]=center+amplitude*math.sin(omega*t);component(c,u);gf=mean*u[int(c.ports[source,0])]/limit;vs+=c.dt/cm*(-gm*(vs-rest)+baseline+c.gs[receiver]*(rev-vs));vg+=c.dt/cm*(-gm*(vg-rest)+baseline+gf*(rev-vg));wave.append([c.time,u[int(c.ports[source,0])]/limit,c.gs[receiver]/mean,vs,vg])
        wave=np.array(wave);mask=wave[:,0]>=1.;basis=np.c_[np.sin(omega*wave[mask,0]),np.cos(omega*wave[mask,0]),np.ones(mask.sum())];co=np.linalg.lstsq(basis,wave[mask,2],rcond=None)[0];phase=float(math.degrees(math.atan2(co[1],co[0])));rmse=float(np.sqrt(np.mean((wave[mask,2]-wave[mask,1])**2)));cs=np.linalg.lstsq(basis,wave[mask,3],rcond=None)[0];cg=np.linalg.lstsq(basis,wave[mask,4],rcond=None)[0];receiver_phase=float(np.degrees(np.angle(np.exp(1j*(math.atan2(cs[1],cs[0])-math.atan2(cg[1],cg[0]))))));err=float(np.sqrt(np.mean((wave[mask,3]-wave[mask,4])**2)));mod=float(np.hypot(cg[0],cg[1])/math.sqrt(2));relative=err/mod;bandwidth[label]=dict(physical_frequency_Hz=frequency,current_center_nA=center,current_amplitude_nA=amplitude,declared_type_max_nA=float(limit),phase_at_end_of_host_interval_deg=phase,phase_relative_held_sample_deg=phase+math.degrees(omega*c.dt),normalized_encoding_domain_G_RMSE=rmse,source5mV_window_G_RMSE=rmse*limit/5,raw_G_modulation_relative_RMSE=rmse/(amplitude/limit),ten_degree_phase_gate=abs(phase)<=10,actual_first_receiver=names[receiver],receiver_probe=dict(source_C_microF=float(cm),source_G_microS=float(gm),source_reversal_mV=rev,baseline_mV=2.5,spike_vs_matched_graded_voltage_RMSE_mV=err,graded_modulation_RMS_mV=mod,modulation_relative_RMSE=relative,spike_vs_matched_graded_extra_phase_deg=receiver_phase,design_dynamic_gate=relative<=.05 and abs(receiver_phase)<=10,reference='Same passive receiver driven by continuousG proportional to held physicalI; excludes presynaptic graded encoder pole and all other network currents.',scope='Isolated actual passive membrane parameters at declared2.5mV operating bias; no body or voltageclamp'),scope='Actual event conductance component, from exact quiet initialization; no neural/body standing claim');waveforms[label]=wave
    bench={}
    for mode in ('reference','numba'):
        c=make(backend=mode);u=z.copy();u[int(c.ports[0,0])]=c.input_limits[0]*.5;c.forward(u);start=time.perf_counter()
        for _ in range(1000):c.forward(u)
        bench[mode]=time.perf_counter()-start
    if not population:
        c=EventConductanceInputPort(raw,names,**args);u=z.copy();u[int(c.ports[0][0])]=2.5;start=time.perf_counter()
        for _ in range(1000):c.forward(u);c.events.clear()
        bench['old_equation_reference_bounded_history']=time.perf_counter()-start
    assert hashes=={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
    result=dict(family=family,profile=p,scalar_JIT_max_state_error=max_state,reference_kind='Scalar same-equation midpoint-bias extension' if population else 'Preserved v1 Python equation implementation',Python_reference_max_state_error=max_python,event_time_max_error_s=max_events,complete_state32_exact=True,disabled160_exact=True,outgoing_off_nonencoded160_exact=True,outgoing_off_actual_spikes=int(pulse_spikes),short_nonlinear_window_rejected=rejected,out_of_domain_rejected_before_mutation=True,default_unbounded_event_history_entries=len(a.events),bounded_last_block_capacity=a.event_block.shape[0],DC_points=points,quiet_onset_delays=delays,bandwidth=bandwidth,steady_DC_analytical_fullscale_bounds=dc_bounds,exact_population_DC_1001_error_range=[min(curve),max(curve)],population_DC_bound_derivation='f_j/F<=max(0,x+(b_j/theta-.5)*theta/R),lowermean>=max(0,x-theta/(2R));summean bias=theta/2. Boundresetfactor by(1-delta)..1 then divide exact full-population normalizationD.',population_policy=a.policy,bounded_neural_benchmark=dict(steps=1000,physical_neural_duration_s=2.,warmup_excluded_s=warmup,wall_s=bench,speedup=bench['reference']/bench['numba'],vs_old_equation_reference=bench.get('old_equation_reference_bounded_history',None)/bench['numba'] if not population else None,full_body_claim=False),source_files_unchanged=True,K6_basis_unchanged=True,consumption=n.fsa_spiking_port_ledger,training_allowed=False)
    result['all_tested_role_components_pass']=all(v['receiver_probe']['design_dynamic_gate'] and v['ten_degree_phase_gate'] for v in bandwidth.values())
    write(out/'PROTOCOL.json',protocol);write(out/'RESULTS.json',result);write(out/'PORT_PROFILE.json',p);np.savez_compressed(out/'SELECTED_TRACES.npz',trace=np.array(data),**{label+'_bandwidth':v for label,v in waveforms.items()},profile_json=json.dumps(p,sort_keys=True),result_json=json.dumps(result,sort_keys=True),protocol_json=json.dumps(protocol,sort_keys=True),columns=['time_s',*['current_nA_'+s for s in p['input_names']],*['APcount_'+s for s in p['input_names']],'right_soleus_MN_mV'])
    write(out/'COMPLETED.json',dict(artifacts={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ('PROTOCOL.json','RESULTS.json','PORT_PROFILE.json','SELECTED_TRACES.npz')},electrical_run_completed=True,all_tested_role_components_pass=result['all_tested_role_components_pass'],training_allowed=False));(out/'PENDING_PROTOCOL.json').unlink()
    print(json.dumps(dict(family=family,parity_max=max_state,event_error=max_events,normal_Ib_onset_s=delays[0]['first_actual_spike_s'],role_receiver_relative_error={k:v['receiver_probe']['modulation_relative_RMSE'] for k,v in bandwidth.items()},all_roles=result['all_tested_role_components_pass'],neural_speedup=result['bounded_neural_benchmark']['speedup'],short_window_rejected=True)),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--family',required=True,choices=['s3k','w2lvar','syn6']);a.add_argument('--output',required=True);a.add_argument('--population-bias',action='store_true');a.add_argument('--typed',action='store_true');a.add_argument('--tight-tail',action='store_true');v=a.parse_args();main(v.family,v.output,v.population_bias or v.typed or v.tight_tail,v.typed or v.tight_tail,v.tight_tail)
