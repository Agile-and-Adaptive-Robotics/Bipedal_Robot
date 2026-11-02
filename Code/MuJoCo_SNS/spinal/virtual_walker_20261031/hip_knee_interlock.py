"""Opt-in heteronymous II timing hypothesis, not a literature replication.

Hip-extensor stretch excites a dedicated graded inhibitory interneuron, which
delays homonymous knee-extensor recruitment until hip retraction. Existing
Ia cells are not reused because their source RG/PF collaterals change phase.
"""
import inspect,json,hashlib
from pathlib import Path
import numpy as np
from import_virtual_walker import mujoco
import native_kernel as K

OPTIONAL_KEYS=('hip_II_knee_interlock_gmax_uS','hip_II_knee_interlock_threshold_mV','hip_II_knee_interlock_saturation_mV','hip_II_knee_flexor_gmax_uS','hip_II_flexor_primitive_input_gmax_uS','hip_II_primitive_load_inhibition_gmax_uS')


def build(c):
    """Literal chemical edges; gain0 creates no extra cells or pathways."""
    g=c.configuration.get(OPTIONAL_KEYS[0],0.);flex=c.configuration.get(OPTIONAL_KEYS[3],0.);primitive=c.configuration.get(OPTIONAL_KEYS[4],0.);load=c.configuration.get(OPTIONAL_KEYS[5],0.);assert 0<=g<=.5 and 0<=flex<=1.98 and 0<=primitive<=.8 and (primitive==0 or flex>0) and 0<=load<=2.05 and (load==0 or primitive>0)
    c.heteronymous_nodes=0;c.heteronymous_edges=0;c.heteronymous_edits=[]
    c.heteronymous_primitive_indices=[]
    c.heteronymous_load_indices=[];c.heteronymous_load_current_nA=np.zeros(0);c.heteronymous_load_logs=[]
    if g==0 and flex==0:return
    assert c.sensors is not None,'Requires the actual spike-driven source II encoders'
    threshold=c.configuration[OPTIONAL_KEYS[1]];sat=c.configuration[OPTIONAL_KEYS[2]];assert -60<=threshold<sat<=-39
    templates=[c.names.index(s+' Knee flx Ia') for s in 'LR'];assert np.all(c.params[templates,3]==0) and np.all(c.tonic[templates]==0)
    assert np.array_equal(c.params[templates[0],:4],c.params[templates[1],:4]),'New paired IN membrane dynamics must match'
    indices=np.arange(len(c.v),len(c.v)+2);p=c.params[templates].copy();rows=[]
    def append_cells(parameters,names):
        start=len(c.v);c.params=np.vstack([c.params,parameters]);rest=parameters[:,0]
        for name in ('v','newv','m','h','gs','gi','applied','tonic','last_test_impulse_nA_s'):
            setattr(c,name,np.concatenate([getattr(c,name),rest.copy() if name in ('v','newv') else np.zeros(2)]))
        c.names.extend(names);c.heteronymous_nodes+=2;return np.arange(start,start+2)
    if g:indices=append_cells(p,[s+' Hip-II knee-extension IN' for s in 'LR'])
    if primitive:
        parameters=c.params[[c.names.index(s+' Knee PF flx') for s in 'LR']].copy();parameters[:,3]=0.;assert np.all(parameters[:,0]==-60) and np.all(parameters[:,2]==1)
        pf=append_cells(parameters,[s+' Hip-II flexor primitive PF' for s in 'LR']);c.heteronymous_primitive_indices=pf.tolist();c.encoded=np.concatenate([c.encoded,pf]);c.pf_mask=np.concatenate([c.pf_mask,np.ones(2,dtype=bool)])
        for name in ('lv','lm','lh','ref','adapt'):
            extra=np.zeros((2,c.population_size))
            if name=='lv':extra[:]=parameters[:,0,None]+.2*(np.arange(c.population_size)/(c.population_size-1)-.5)
            setattr(c,name,np.vstack([getattr(c,name),extra]))
        c.decoded=np.concatenate([c.decoded,np.zeros(2)]);c.counts=np.concatenate([c.counts,np.zeros(2,dtype=np.int64)]);c.total_counts=np.concatenate([c.total_counts,np.zeros(2,dtype=np.int64)])
    if load:
        loads=append_cells(p,[s+' Foot-load primitive IN' for s in 'LR']);c.heteronymous_load_indices=loads.tolist();c.heteronymous_load_current_nA=np.zeros(2)
    for side,i,template in zip('LR',indices,templates):
        pre=c.names.index(side+' Hip extensor II relay');post=c.names.index(side+' Knee MN ext');down=c.edges[(c.edges[:,0]==template)&(c.edges[:,1]==post)];assert len(down)==1 and down[0,5]<-60
        edge_rows=([pre,i,g,threshold,sat,-10.],[i,post,*down[0,2:]]) if g else []
        if primitive:
            cell=pf[0 if side=='L' else 1];edge_rows=list(edge_rows)+[[pre,cell,primitive,threshold,sat,-10.],[cell,c.names.index(side+' Knee MN flx'),flex,-60.,-40.,-10.]]
            if load:edge_rows+=[[loads[0 if side=='L' else 1],cell,load,-60.,-40.,-100.]]
        elif flex:edge_rows=list(edge_rows)+[[pre,c.names.index(side+' Knee MN flx'),flex,threshold,sat,-10.]]
        for edge in edge_rows:
            rows.append(edge);c.heteronymous_edits.append(dict(source=c.names[int(edge[0])],target=c.names[int(edge[1])],gmax_uS=float(edge[2]),threshold_mV=float(edge[3]),saturation_mV=float(edge[4]),reversal_mV=float(edge[5]),template=c.names[template]))
    c.edges=np.vstack([c.edges,rows]);c.heteronymous_edges=len(rows)


def audit(root):
    root=Path(root);fit=json.loads((root/'OFFLINE_NEURAL_FIT.json').read_text());base=Path(fit['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());cfg=rec['source_configuration'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];params=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];allowed=set(inspect.signature(K.NativeController).parameters)-{'source','bank'};c=K.NativeController(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);q=np.load(root/'LANDING_PHASE_PAIR_TARGETS.npz')['qpos'];meta=json.loads((root/'LANDING_PHASE_PAIR_TARGETS.json').read_text());tids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+b['name']) for b in bank]);bindings={c.names[pre]:i for pre,i,coeff in c.outputs};rows=[];stance=[];swing=[]
    for phase,pose in zip(meta['targets'],q):
        d.qpos[:]=pose;d.qvel[:]=0;mujoco.mj_forward(m,d);L=d.ten_length[tids].copy();internal=np.zeros(len(bank));force=internal.copy();ia=internal.copy();ii=internal.copy();ib=internal.copy()
        for k,p in enumerate(params):internal[k]=(p[3]*(L[k]-p[0])+max(0.,p[7]/(1+np.exp(np.clip(p[8]*(p[6]+.1),-700,700)))+p[9])*max(0.,1-((L[k]-p[0])/p[1])**2))/(1+p[3]/p[2])
        K.hill_step(L,np.zeros(len(bank)),np.full(len(bank),-.1),params,internal,force,ia,ii,ib,0.)
        for side in 'LR':
            sensor=c.names.index(side+' Hip extensor II relay');inputs=[v for v in c.inputs if v[0]==sensor];assert len(inputs)==1 and inputs[0][2]=='II';_,i,kind,coeff,body=inputs[0];act=bindings[side+' Hip MN ext'];dof=m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'hip_'+side)];J=d.ten_J.reshape(m.ntendon,m.nv);assert J[tids[i],dof]>0 and J[tids[act],dof]>0,'Both routes must stretch under hip flexion';current=K.polynomial(coeff,ii[i])*1e9;proxy=c.params[sensor,0]+np.clip(current/c.params[sensor,2],0,20);is_stance=side==phase['side'];(stance if is_stance else swing).append(proxy);rows.append(dict(stance_side=phase['side'],side=side,sensor=c.names[sensor],receptor=bank[i]['name'],actuator=bank[act]['name'],receptor_dLength_dHip_m_rad=float(J[tids[i],dof]),actuator_dLength_dHip_m_rad=float(J[tids[act],dof]),raw_II=float(ii[i]),input_nA=float(current),expected_steady_proxy_mV=float(proxy),condition='stance' if is_stance else 'swing',coefficient=coeff.tolist()))
    low=float(np.mean(stance));high=float(np.mean(swing));assert high-low>1.;threshold=low+.25*(high-low);sat=low+.75*(high-low);collaterals=[]
    for side in 'LR':
        pre=c.names.index(side+' Knee flx Ia')
        collaterals.extend(dict(source=c.names[pre],target=c.names[int(e[1])],gmax_uS=float(e[2]),reversal_mV=float(e[5])) for e in c.edges[c.edges[:,0]==pre])
    sensory=np.load(base/'sensory_spikes.npz');z=np.load(base/'trajectory.npz');actual=[]
    for side in 'LR':
        i=sensory['names'].tolist().index(side+' Hip extensor II relay');mask=sensory['t']>=1;current=sensory['physical_input_nA'][mask,i];proxy=-60+20*sensory['decoded_activation'][mask,i];angles=np.interp(sensory['t'][mask],z['t'],z['angles'][:,0 if side=='L' else 3]);actual.append(dict(side=side,actual_input_nA_range=[float(current.min()),float(current.max())],actual_filtered_proxy_mV_range=[float(proxy.min()),float(proxy.max())],hip_angle_deg_range=[float(angles.min()),float(angles.max())],hip_angle_raw_current_correlation=float(np.corrcoef(angles,current)[0,1])))
    report=dict(machine='easteregg2',source_sha256=rec['source_sha256'],baseline=str(base),phase_calibration=rows,baseline_observability=actual,source_Ia_collaterals_not_reused=collaterals,configuration={OPTIONAL_KEYS[0]:.1,OPTIONAL_KEYS[1]:threshold,OPTIONAL_KEYS[2]:sat},scope='Engineering pre-contact timing hypothesis: actual hip-extensor receptor route length increases with hip flexion, and differs between two mechanically feasible landing anchors. Expected steady proxy is an encoder-limit approximation, not a measured target-phase state. Two dedicated graded INs copy only the source inhibitory knee-MN output; source Ia RG/PF collaterals remain untouched. No force/torque allocator, pose clamp, source-model mutation, or reduced-support claim.');(root/'HIP_II_KNEE_INTERLOCK_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def probe(root):
    from spiking_central_runtime import SpikingCentralController as C
    from recruitment import apply_recruitment
    root=Path(root);audit=json.loads((root/'HIP_II_KNEE_INTERLOCK_AUDIT.json').read_text());base=Path(audit['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());cfg=rec['source_configuration'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'};controls={};peak=0.;windows={};state_checks=None
    for label,gain in [('absent',None),('zero',0.),('enabled',.1)]:
        configuration=rec['configuration'].copy()
        if gain is not None:configuration.update(audit['configuration']);configuration[OPTIONAL_KEYS[0]]=gain
        C.configuration=configuration;C.sensory_spiking=True;C.balance_configuration=None;C.clock_scale=1.;C.audit_state=False;C.reference_calculations=False;c=C(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});c.configuration=configuration.copy();apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();rows=[];zero=np.zeros(len(bank))
        def advance(k):
            ii=zero.copy();condition='swing' if k<500 else 'stance'
            for side in 'LR':
                target=c.names.index(side+' Hip extensor II relay');idx=next(v[1] for v in c.inputs if v[0]==target);ii[idx]=next(v['raw_II'] for v in audit['phase_calibration'] if v['side']==side and v['condition']==condition)
            c.step(k*.002,zero,ii,zero,{},.002,True)
        for k in range(1250):
            advance(k)
            if gain:
                edges=c.edges[-2:] # Select outputs by actual target below, not row order.
                values=[]
                for side in 'LR':
                    pre=c.names.index(side+' Hip-II knee-extension IN');post=c.names.index(side+' Knee MN ext');edge=c.edges[(c.edges[:,0]==pre)&(c.edges[:,1]==post)][0];I=edge[2]*np.clip((c.v[pre]-edge[3])/(edge[4]-edge[3]),0,1)*(edge[5]-c.v[post]);values.extend([c.v[pre],c.v[post],I]);peak=max(peak,abs(float(I)))
                rows.append([k*.002,*values])
                if k==1000:
                    initial=c.export_state();states=[]
                    for repeat in range(2):
                        c.import_state(initial)
                        for j in range(1,33):advance(k+j)
                        states.append(c.export_state())
                    state_checks={key:bool(value.shape==states[1][key].shape and value.dtype==states[1][key].dtype and value.tobytes()==states[1][key].tobytes()) for key,value in states[0].items()};assert all(state_checks.values());c.import_state(initial)
        controls[label]=c.export_state()
        if gain:
            arr=np.array(rows)
            for name,lo,hi in [('hip_flexed',.7,1.),('hip_retracted',2.2,2.5)]:windows[name]=dict(zip(['L_IN_mV','L_KneeMNExt_mV','L_inhibitory_current_nA','R_IN_mV','R_KneeMNExt_mV','R_inhibitory_current_nA'],arr[(arr[:,0]>=lo)&(arr[:,0]<hi),1:].mean(axis=0).tolist()))
            actual=c.metadata()
    keys=set(controls['absent'])-{'configuration'};eq={key:bool(controls['absent'][key].dtype==controls['zero'][key].dtype and controls['absent'][key].shape==controls['zero'][key].shape and controls['absent'][key].tobytes()==controls['zero'][key].tobytes()) for key in keys};assert all(eq.values());result=dict(zero_gain_default_byte_equivalence=eq,state_continuation_byte_equivalence=state_checks,sampled_2ms_added_current_peak_nA=peak,steady_input_windows=windows,architecture=actual,scope='2.5s isolated controller; only actual hip-extensor-II physical inputs switch from flexed-anchor to retracted-anchor values at1s. Other afferents/contact zero; no mechanical simulation or stepping claim. Existing source neural startup stimulus retained. Currents are sampled at2ms controller boundaries, not exact neural-substep maxima.');(root/'INTERLOCK_ISOLATED_CHECK.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result


def summarize(root):
    """Apply identical recorded-force touchdown and timing definitions to all cases."""
    root=Path(root);audit=json.loads((root/'HIP_II_KNEE_INTERLOCK_AUDIT.json').read_text());folders=[Path(audit['baseline'])]+sorted((root/'W2L_VirtualWalker_Newer').glob('hipII_*_v1'));rows=[]
    for folder in folders:
        rec=json.loads((folder/'SPIKING_COMPARISON.json').read_text());r=rec['result'];a=rec['assessment'];z=np.load(folder/'trajectory.npz');bw=r['robot_mass_kg']*9.81;t=z['t'];forces=np.stack([z[k][:,1] for k in ('rail_force_world','guide_force_world','harness_force')],axis=1);row=dict(path=str(folder),mode=rec['mode'],configuration=rec['configuration'],air_gate=rec.get('air_gate'),ground_gate=a['alternating_supported_stepping_gate'],swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],net_support_BW=r['mean_net_vertical_assistance_BW'],upward_BW=float(np.maximum(forces,0).sum(1).mean()/bw),downward_BW=float(np.maximum(-forces,0).sum(1).mean()/bw),progress_m=r['root_forward_relative_to_belt_m'],failure=r['failure'],knee_min_deg=r['anatomical_knee_min_deg'],warmed_wall_s=r['warmed_wall_s'],spiking=rec['spiking']['actual_spiking_cells'],graded=rec['spiking']['remaining_graded_cells'],landings=[]);names=z['neuron_names'].tolist();peaks=[]
        for edge in rec['spiking'].get('heteronymous_edges',[]):
            pre=z['neuro'][:,names.index(edge['source'])];post=z['neuro'][:,names.index(edge['target'])];I=edge['gmax_uS']*np.clip((pre-edge['threshold_mV'])/(edge['saturation_mV']-edge['threshold_mV']),0,1)*(edge['reversal_mV']-post);peaks.append(dict(**edge,sampled_2ms_peak_abs_current_nA=float(np.abs(I).max())))
        row['added_currents']=peaks
        if rec['mode']!='air':
            for s,side in enumerate('LR'):
                knee=z['angles'][:,3*s+1];hip=z['angles'][:,3*s];vel=np.gradient(knee,t)
                for swing in a['legs'][side]['swing_intervals']:
                    start=int(np.searchsorted(t,swing['start_s']));end=int(np.searchsorted(t,swing['end_s']));dwell=max(1,round(.020/r['physics_dt']));candidates=[k for k in range(end,min(len(t)-dwell,int(np.searchsorted(t,swing['end_s']+.15)))) if np.all(z['foot_vertical'][k:k+dwell,s]>.05*bw)]
                    if not candidates:continue
                    strike=candidates[0];peak=start+int(np.argmax(knee[start:end+1]));extension=[k for k in range(peak,max(peak,strike-dwell)) if np.all(vel[k:k+dwell]<-20.)];onset=extension[0] if extension else None
                    row['landings'].append(dict(side=side,touchdown_s=float(t[strike]),hip_at_touchdown_deg=float(hip[strike]),knee_at_touchdown_deg=float(knee[strike]),foot_origin_anterior_to_COM_m=float(z['foot_body_position'][strike,s,0]-z['com_position'][strike,0]),CoP_anterior_to_COM_m=float(z['foot_cop_world'][strike,s,0]-z['com_position'][strike,0]),extension_onset_to_touchdown_s=None if onset is None else float(t[strike]-t[onset]),hip_at_extension_onset_deg=None if onset is None else float(hip[onset])))
            if row['landings']:
                row['median_landing']={key:float(np.median([v[key] for v in row['landings'] if v[key] is not None])) for key in ('hip_at_touchdown_deg','knee_at_touchdown_deg','foot_origin_anterior_to_COM_m','CoP_anterior_to_COM_m','extension_onset_to_touchdown_s','hip_at_extension_onset_deg') if any(v[key] is not None for v in row['landings'])}
        rows.append(row)
    report=dict(cases=rows,promoted=False,timing_rule='Qualified unloaded-clear swing intervals from unchanged gate; first following foot load>.05BW sustained20ms within150ms. Knee extension onset: after maximum swing flexion, first velocity<âˆ’20deg/s sustained20ms. Same definitions applied to baseline and all variants; missing events remain missing. Currents reconstructed at saved2msstates are sampled diagnostics, not guaranteed substep maxima.');(root/'INTERLOCK_GATES.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def fit_spiking_primitive(root):
    """Bounded expected-input readout calibration, before coupled validation."""
    from spiking_central_runtime import SpikingCentralController as C
    from recruitment import apply_recruitment
    root=Path(root);audit=json.loads((root/'HIP_II_KNEE_INTERLOCK_AUDIT.json').read_text());base=Path(audit['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());cfg=rec['source_configuration'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'};target=json.loads((root/'NEURAL_PHASE_CURRENT_MAP.json').read_text());output=float(np.mean([v['exc_gmax_uS'] for v in target['rows'] if v['role']=='knee_flexor' and v['phase']=='stance']));attempts=[];selected=None
    for gain in (.2,.4,.8):
        config=dict(rec['configuration'],**audit['configuration']);config.update(hip_II_knee_interlock_gmax_uS=.3,hip_II_knee_flexor_gmax_uS=output,hip_II_flexor_primitive_input_gmax_uS=gain);C.configuration=config;C.sensory_spiking=True;C.balance_configuration=None;C.clock_scale=1.;C.audit_state=False;C.reference_calculations=False;c=C(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});c.configuration=config.copy();apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();zero=np.zeros(len(bank));trace=[];continuation=None
        def advance(k):
            ii=zero.copy();condition='swing' if k<500 else 'stance'
            for side in 'LR':
                idx=next(v[1] for v in c.inputs if v[0]==c.names.index(side+' Hip extensor II relay'));ii[idx]=next(v['raw_II'] for v in audit['phase_calibration'] if v['side']==side and v['condition']==condition)
            c.step(k*.002,zero,ii,zero,{},.002,True)
        for k in range(1250):
            advance(k);trace.append([k*.002,*np.clip((c.v[c.heteronymous_primitive_indices]+60)/20,0,1)])
            if k==1000:
                initial=c.export_state();states=[]
                for repeat in range(2):
                    c.import_state(initial)
                    for j in range(1,33):advance(k+j)
                    states.append(c.export_state())
                continuation={key:bool(value.dtype==states[1][key].dtype and value.shape==states[1][key].shape and value.tobytes()==states[1][key].tobytes()) for key,value in states[0].items()};assert all(continuation.values());c.import_state(initial)
        a=np.array(trace);on=a[(a[:,0]>=.7)&(a[:,0]<1),1:].mean(axis=0);off=a[a[:,0]>=2.2,1:].mean(axis=0);passed=bool(np.all(on>=.85) and np.all(off<=.1));row=dict(input_gmax_uS=gain,output_gmax_uS=output,mean_flexed_chemical_activation=on.tolist(),mean_retracted_chemical_activation=off.tolist(),passed=passed,state_continuation_all_byte_identical=all(continuation.values()),state_field_count=len(continuation),spiking_cells=c.metadata()['actual_spiking_cells'],graded_cells=c.metadata()['remaining_graded_cells'],PF_spike_counts=c.total_counts[-2:].tolist());attempts.append(row)
        if passed:selected=config;break
    report=dict(attempts=attempts,selected_configuration=selected,source_sha256=rec['source_sha256'],expected_source_current_envelope_nA=178.17201581560943,new_projection_max_current_at_minus100mV_nA=90*output,scope='Two extra16-cell adapting LIF PF populations, without intrinsic calcium, integrate actual spike-driven HipExtII input and recruit knee-flexor MNs. All original central/sensory settings and physical parameters frozen. Output projection is the mechanically feasible stance-anchor source-current fit, not multiplied by original-PF attenuation. Input gain chosen from three predeclared expected-input tests: mean chemical activation>=.85 flexed,<=.1 retracted. This is a limited neural input/readout calibration, not complete-cycle Hunt training or a coupled walking pass.');(root/'HIP_II_SPIKING_PRIMITIVE_FIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    if selected:(root/'HIP_II_SPIKING_PRIMITIVE_CANDIDATE.json').write_text(json.dumps(selected,indent=2),encoding='utf-8')
    return report


def probe_load_gate(root):
    from spiking_central_runtime import SpikingCentralController as C
    from recruitment import apply_recruitment
    root=Path(root);audit=json.loads((root/'HIP_II_KNEE_INTERLOCK_AUDIT.json').read_text());fit=json.loads((root/'HIP_II_SPIKING_PRIMITIVE_FIT.json').read_text());config=dict(fit['selected_configuration'],hip_II_primitive_load_inhibition_gmax_uS=2.);base=Path(audit['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());cfg=rec['source_configuration'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'};C.configuration=config;C.sensory_spiking=True;C.balance_configuration=None;C.clock_scale=1.;C.audit_state=False;C.reference_calculations=False;c=C(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});c.configuration=config.copy();apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));bw=rec['result']['robot_mass_kg']*9.81;zero=np.zeros(len(bank));trace=[];checks=None
    def advance(k):
        ii=zero.copy()
        for side in 'LR':
            idx=next(v[1] for v in c.inputs if v[0]==c.names.index(side+' Hip extensor II relay'));ii[idx]=next(v['raw_II'] for v in audit['phase_calibration'] if v['side']==side and v['condition']=='swing')
        load=.2*bw if 1<=k*.002<1.8 else 0.;c.observe_load_feedback(m,dict(foot_force_world=np.array([[0.,load,0.],[0.,load,0.]])));c.step(k*.002,zero,ii,zero,{},.002,True)
    for k in range(1250):
        advance(k);trace.append([k*.002,*np.clip((c.v[c.heteronymous_primitive_indices]+60)/20,0,1),*c.applied[c.heteronymous_load_indices]])
        if k==1000:
            initial=c.export_state();states=[]
            for repeat in range(2):
                c.import_state(initial)
                for j in range(1,33):advance(k+j)
                states.append(c.export_state())
            checks={key:bool(value.dtype==states[1][key].dtype and value.shape==states[1][key].shape and value.tobytes()==states[1][key].tobytes()) for key,value in states[0].items()};assert all(checks.values());c.import_state(initial)
    a=np.array(trace);windows={}
    for name,lo,hi in [('unloaded',.7,1.),('loaded',1.4,1.8),('released',2.2,2.5)]:windows[name]=a[(a[:,0]>=lo)&(a[:,0]<hi),1:3].mean(axis=0).tolist()
    passed=bool(min(windows['unloaded'])>=.85 and max(windows['loaded'])<=.1 and min(windows['released'])>=.85);report=dict(passed=passed,windows=windows,actual_load_IN_current_peak_nA=a[:,3:].max(axis=0).tolist(),complete_state_byte_checks=checks,architecture=c.metadata(),scope='Isolated expected-input test with HipExtII held flexed throughout. Actual normal-load observation mapping0â†’.2BWâ†’0 independently suppresses/re-enables only the extra spiking PF; two graded load-INs, input ceiling8nA, output2ÂµS. No physical standing/stepping claim. Original heel/toe feedback unchanged.');(root/'HIP_II_LOAD_GATE_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    if passed:(root/'HIP_II_LOAD_GATED_PRIMITIVE_CANDIDATE.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    return report
