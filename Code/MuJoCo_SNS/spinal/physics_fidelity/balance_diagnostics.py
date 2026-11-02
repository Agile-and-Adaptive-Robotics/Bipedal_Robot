"""Isolated source-consistent task calibration and signed standing diagnostics."""
import json
from pathlib import Path
import numpy as np,mujoco
from thelen_dynamic import DynamicThelenBank,RigidThelenBank,parameters_from_osim

def install_posture_scales(scales):
    """Scale only declared POST input ports, tied by functional group."""
    import controller as C
    allowed={'ankle_pf','ankle_df','hip_flex','hip_ext','knee_flex','knee_ext','trunk_ext','trunk_flex','all'}
    if set(scales)-allowed or any(not np.isfinite(v) or not 0.<=v<=2. for v in scales.values()):raise ValueError('Invalid isolated posture scale')
    original=C.build
    def build(*args,**kwargs):
        n=original(*args,**kwargs);ports=[];gains=[];bindings=[]
        for act,info in n.muscles.items():
            scale=float(scales.get('all',1.))*float(scales.get(info.groups[0],1.))
            if scale!=1.:
                ports.append(n.ports['POST_'+act]);gains.append(scale);bindings.append(dict(actuator=act,group=info.groups[0],scale=scale))
        ports=np.array(ports,dtype=int);gains=np.array(gains);native=n.step
        def step(u):u[ports]*=gains;return native(u)
        n.step=step;n.isolated_posture_bindings=bindings
        return n
    C.build=build

def install_source_task_basis(osim,mode='compliant'):
    """Patch only this subprocess's TaskBalance.bind; restore physics afterward.

    Finite activation probes use source static tendon/fiber equilibrium at the
    actual pose. This changes a task-recruitment calibration, not source muscle
    dimensions or runtime controller targets. Native initialization remains.
    """
    import balance as B
    original_bind=B.TaskBalance.bind
    def bind(self,m,d):
        original_forward=mujoco.mj_forward;fields=('actuator_gaintype','actuator_gainprm','actuator_biastype','actuator_biasprm');saved={k:getattr(m,k).copy() for k in fields}
        p=parameters_from_osim(osim,[m.actuator(i).name for i in range(m.nu)],0.)
        bank=(DynamicThelenBank if mode=='compliant' else RigidThelenBank)(p,float(m.opt.timestep));calls=0;error=0.
        def forward(mm,dd):
            nonlocal calls,error
            original_forward(mm,dd)
            if mm is m and dd is not d:
                bank.initialize(dd.actuator_length,dd.act)
                if np.any(bank.bounds):raise RuntimeError('Task calibration source-fiber bound')
                error=max(error,float(np.max(abs(bank.error))))
                mm.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;mm.actuator_gainprm[:]=0.
                mm.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE;mm.actuator_biasprm[:]=0.;mm.actuator_biasprm[:,0]=-bank.force
                original_forward(mm,dd);calls+=1
        try:
            mujoco.mj_forward=forward;original_bind(self,m,d)
            self.calibration.update(force_probe_backend='source_'+mode,static_equilibrium_probe_calls=calls,max_probe_equilibrium_error_N=error,probe_scope='Source steady fiber/tendon equilibrium at fixed pose and activation; native source dimensions, no runtime inverse allocator.')
        finally:
            mujoco.mj_forward=original_forward
            for key,value in saved.items():getattr(m,key)[:]=value
            original_forward(m,d)
    B.TaskBalance.bind=bind

def summarize(root):
    from tendon_replay_campaign import build
    result={}
    for folder in root.iterdir():
        if not all((folder/name).exists() for name in ('controller.result.json','backend_trace.npz','BACKEND_RESULT.json')):continue
        r=json.loads((folder/'controller.result.json').read_text());z=np.load(folder/'controller.npz',allow_pickle=True);a=np.load(folder/'controller.assistance.npz')['samples'];cfg=json.loads((folder/'controller.json').read_text());m,d=build(z['qfull'][0],cfg,1.)
        ref=np.array([np.mean([d.xpos[m.body('talus_'+s).id,0] for s in ('r','l')]),z['com'][0,1]])
        q=z['com'][:,:2];rel=q-ref;delta=q-q[0];k=int(np.argmax(np.linalg.norm(rel,axis=1)));late=a[:,0]>=4.;bw=r['measured']['body_weight_N'];fz=(a[late,1]+a[late,2])/bw
        result[folder.name]=dict(assessment=r['assessment'],initial_com_minus_ankle_reference_xy_m=rel[0].tolist(),maximum_com_excursion_from_initial_m=float(np.linalg.norm(delta,axis=1).max()),last_com_minus_reference_xy_m=rel[-1].tolist(),time_of_maximum_reference_radius_s=float(z['t'][k]),signed_apparatus_vertical_mean_bw=float(fz.mean()),upward_apparatus_mean_bw=float(np.maximum(fz,0).mean()),downward_apparatus_mean_bw=float(np.maximum(-fz,0).mean()),measured_horizontal_rms_bw=r['measured']['horizontal_harness_rms_bw'])
        trace=np.load(folder/'backend_trace.npz');dt=float(np.median(np.diff(trace['t'])));pulse=trace['applied_pulse'];active=np.any(pulse!=0,axis=1)
        result[folder.name].update(apparatus_analysis_window_s=[4.,float(a[-1,0])],actual_pulse_nonzero_steps=int(active.sum()),actual_pulse_starts=int(np.count_nonzero(np.diff(np.r_[False,active].astype(int))==1)),actual_pulse_impulse_Ns=(pulse.sum(axis=0)*dt).tolist())
    (root/'BALANCE_DIAGNOSIS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

def neural_currents(folder):
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
    import params as P
    z=np.load(folder/'controller.neurons.npz');a=np.load(folder/'controller.assistance.npz')['samples'];r=json.loads((folder/'controller.result.json').read_text());g=json.loads((folder/'controller.graph.json').read_text())
    names={int(n['id'][1:]):n['label'] for n in g['nodes']};cells=z['cells'].tolist();columns={n:i for i,n in enumerate(cells)};v=z['voltage'];old=np.vstack([np.zeros((1,v.shape[1])),v[:-1]]);keep=(z['t']>=5)&(z['t']<13);bw=r['measured']['body_weight_N'];out={}
    for side,si in (('r',0),('l',1)):
        stance=keep&(a[:,7+si]>.05*bw)&(old[:,columns['RG_E_'+side]]>old[:,columns['RG_F_'+side]])
        for muscle in ('vas_lat','soleus'):
            target='MN_'+muscle+'_'+side;currents={k:np.zeros(len(v)) for k in ('Ia_homonymous','Ib_inhibitory','Ib_excitatory','PT_excitatory_recorded')};edges=[]
            for e in g['edges']:
                pre,post=names[e['source']],names[e['target']]
                if post!=target or pre not in columns:continue
                cat='Ia_homonymous' if pre.startswith('Ia_') else 'Ib_inhibitory' if pre.startswith('IBIN_') else 'Ib_excitatory' if pre.startswith('IBEXC_') else 'PT_excitatory_recorded' if pre.startswith('PT_BAL_') else None
                if cat is None:continue
                current=e['walking']*np.clip(old[:,columns[pre]]/P.E_HI,0.,1.)*(e['reversal']-old[:,columns[post]])
                currents[cat]+=current;edges.append(dict(source=pre,gmax=e['walking'],reversal_mV=e['reversal']))
            out[target]=dict(loaded_extensor_phase_samples=int(stance.sum()),selected_currents_nA={k:dict(walking_mean=float(x[keep].mean()),loaded_extensor_phase_mean=float(x[stance].mean()) if np.any(stance) else None,walking_peak_abs=float(abs(x[keep]).max())) for k,x in currents.items()},edges=edges)
    result=dict(E_HI_mV=P.E_HI,scope='Reconstructed selected synaptic currents using previous recorded voltage, BN._syn thresholds0/E_HI, actual compiled reversals, and fixed fitted reflex/PT conductances. Not total MN current; PT inhibitory relays and reciprocal IaIN were not recorded.',neurons_max_walking_mV={n:float(v[keep,columns[n]].max()) for n in cells if n.startswith(('IBEXC_','PT_BAL_'))},motor_neurons=out)
    (folder/'SELECTED_SYNAPTIC_CURRENTS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

def pool_bias_probe(folder):
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent));import params as P
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    from sns_toolbox.connections import NonSpikingSynapse
    n=Network(name='isolated_phase_bias');n.add_neuron(NonSpikingNeuron(membrane_capacitance=P.TAU['ib_exc'],membrane_conductance=1.,bias=P.E_HI),name='phase');n.add_neuron(NonSpikingNeuron(membrane_capacitance=P.TAU['ib_exc'],membrane_conductance=1.),name='pool');n.add_connection(NonSpikingSynapse(max_conductance=1.,reversal_potential=8.,e_lo=0.,e_hi=P.E_HI),'phase','pool');n.add_output('pool');model=n.compile(dt=.002,backend='numpy')
    for _ in range(1000):model.forward(np.empty(0))
    isolated=float(model.V[1]);assert abs(isolated-4.)<1e-10
    z=np.load(folder/'controller.neurons.npz');cells=z['cells'].tolist();v=z['voltage'];columns={s:i for i,s in enumerate(cells)};dt=float(np.median(np.diff(z['t'])));keep=(z['t']>=5)&(z['t']<13);out={}
    for side in ('r','l'):
        rg=v[:,columns['RG_E_'+side]]
        for group in ('knee_ext','ankle_pf'):
            name='IBEXC_'+group+'_'+side;actual=v[:,columns[name]];cf=np.empty(len(v));cf[0]=actual[0]
            for k in range(1,len(v)):
                conductance=np.clip(rg[k-1]/P.E_HI,0.,1.);cf[k]=cf[k-1]+dt/P.TAU['ib_exc']*(-cf[k-1]+conductance*(8.-cf[k-1]))
            out[name]=dict(recorded_peak_mV=float(actual[keep].max()),RG_only_peak_mV=float(cf[keep].max()),recorded_mean_mV=float(actual[keep].mean()),RG_only_mean_mV=float(cf[keep].mean()),direct_Ib_increment_mean_mV=float((actual[keep]-cf[keep]).mean()),RG_only_fraction_of_pool_activity=float(cf[keep].mean()/max(actual[keep].mean(),1e-12)))
    result=dict(isolated_SNS_no_Ib_RG_at_5mV_pool_steady_mV=isolated,scope='Direct-input causal probe: freeze recorded RG_E waveform, set all Ib afferent inputs to the pool to zero, replay exact RC Euler update with leak1,tau0.05s,g_RG1,rev8mV. Not a whole-body closed-loop ablation.',pools=out)
    (folder/'IB_POOL_PHASE_BIAS_PROBE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

def ib_runtime_currents(folder):
    """Summarize recorded, time-varying Ib currents rather than final graph gains."""
    z=np.load(folder/'source_ib_currents.npz');s=np.load(folder/'source_ib_signals.npz')['samples'];t=z['neural_t_s'];keep=(t>=5)&(t<13);edges=z['edges'].tolist();rows={}
    for side,offset in [('r',0),('l',1)]:
        masks=dict(walking=keep,loaded_extensor=keep&(s[:,2+offset]>0)&(s[:,4+offset]>=1.),unloaded_flexor=keep&(s[:,2+offset]==0)&(s[:,4+offset]==0))
        for family in ('vas_lat','vas_med','vas_int','soleus','med_gas','lat_gas','per_tert','rect_fem','bifemsh','tib_ant'):
            act=family+'_'+side;groups={cat:[i for i,e in enumerate(edges) if e.endswith(' -> MN_'+act) and e.startswith(prefix)] for cat,prefix in [('excitatory_pool','IBEXC_'),('autogenic_inhibitory_relay','IBIN_')]};states={}
            for state,mask in masks.items():
                if not np.any(mask):continue
                record={'samples':int(mask.sum())};net=0.
                for cat,indices in groups.items():
                    current=z['current_nA'][mask][:,indices].sum(axis=1);mean=float(current.mean());net+=mean;record[cat]=dict(mean_current_nA=mean,range_current_nA=[float(current.min()),float(current.max())],maximum_gmax=float(z['gmax'][mask][:,indices].max()) if indices else 0.)
                record['net_selected_Ib_current_nA']=net;states[state]=record
            rows[act]=states
    result=dict(scope='Recorded pre-step synaptic currents with actual time-varying conductances and compiled voltage thresholds/reversals. Selected Ib sources only, not total MN current. Loaded-extensor means actual regional contact>0.05BW and endogenous RG_E>RG_F.',motor_neurons=rows)
    (folder/'IB_RUNTIME_CURRENTS.json').write_text(json.dumps(result,indent=2));print(json.dumps({name:row.get('loaded_extensor',{}) for name,row in rows.items() if name.startswith(('soleus','vas_lat'))},indent=2));return result

def backend_parity(left,right,target):
    """Compare actual saved arrays, preserving both source simulations."""
    x=np.load(left/'backend_trace.npz');y=np.load(right/'backend_trace.npz');rows={}
    for key in ('t','qpos','force','activation','ctrl','fiber_state','fiber_force_components'):
        equal=x[key].shape==y[key].shape;rows[key]=dict(same_shape=equal,bit_exact=bool(equal and np.array_equal(x[key],y[key])),maximum_absolute_difference=float(np.max(abs(x[key]-y[key]))) if equal else None)
    result=dict(left=str(left),right=str(right),all_bit_exact=all(row['bit_exact'] for row in rows.values()),arrays=rows);target.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result

def ib_checkpoint(root):
    """Compact completion record for the bounded force-feedback diagnosis."""
    import hashlib
    here=Path(__file__).resolve().parent;study=json.loads((root/'load_recruitment_v1/SUMMARY.json').read_text());trials=[json.loads(line) for line in (root/'load_recruitment_v1/TRIALS.jsonl').read_text().splitlines()];best=Path(study['best']['diagnostics']['case']);rows=[]
    for trial in trials:
        d=trial['diagnostics'];case=Path(d['case']);r=json.loads((case/'controller.result.json').read_text());b=json.loads((case/'BACKEND_RESULT.json').read_text());a=np.load(case/'controller.assistance.npz');samples=a['samples'];bw=r['measured']['body_weight_N'];window=(samples[:,0]>=5)&(samples[:,0]<13);vertical=(samples[window,1]+samples[window,2])/bw
        rows.append(dict(number=trial['number'],parameters=trial['params'],objective=trial['value'],engineering_accepted=r['assessment']['accepted'],reference_envelope_accepted=d['biomechanical_envelope']['accepted_reference_envelope'],envelope_cost=d['biomechanical_envelope']['cost'],sustained_swings=[foot['sustained_swings'] for foot in r['foot_clearance']['feet']],signed_apparatus_vertical_mean_bw=float(vertical.mean()),upward_mean_bw=float(np.maximum(vertical,0.).mean()),downward_mean_bw=float(np.maximum(-vertical,0.).mean()),maximum_knee_extension_deg=b['backend_diagnostics'][0]['maximum_knee_extension_deg'],path=str(case)))
    fast=json.loads((root/'best_fast_fresh/BACKEND_RESULT.json').read_text());regular=json.loads((root/'best_native_fresh/BACKEND_RESULT.json').read_text());fr=json.loads((root/'best_fast_fresh/controller.result.json').read_text());rr=json.loads((root/'best_native_fresh/controller.result.json').read_text());pt=json.loads((root/'pt_current_chain_v2/PT_CURRENT_CHAIN.json').read_text())
    profiles={'shared_force_law_functional_v1.json':'67f0cad8e69d81f1042fac30e1ebe023e93cd6d0a924897ef8bf1a2e19d9bb5a','shared_force_law_functional_v2.json':'a288dd88ceb3888e69b13f40bcbd140189f6d7aad1e774d45499413bf65271a7','shared_force_law_functional_v3.json':'9dc7d7ba444a961e06fd6efd15210ea324a84ee70992b251994cb7d7ee7165cf'}
    frozen={name:hashlib.sha256((here/name).read_bytes()).hexdigest()==expected for name,expected in profiles.items()};frozen['runner.py']=hashlib.sha256((here.parent/'runner.py').read_bytes()).hexdigest()=='585f4e6274df8d6e9aebf1160400e218420741bde1c34c11c6294c7755fad4a5'
    result=dict(model='Gait2392 human model, source-compliant Thelen2003, not a BPA robot',completed_trials=len(trials),fully_accepted_trials=sum(row['engineering_accepted'] and row['reference_envelope_accepted'] for row in rows),advanced_support_ladder=False,best_diagnostic_case=str(best),scope='All anatomical source dimensions and existing FF/POST/PT/timing fixed. Only two common bilateral load-output conductances were searched. Reference envelope is the saved retargeted subject01 IK, not universal human ROM. Original engineering pass and anatomical failure are separate.',frozen_spec_sha256=hashlib.sha256((root/'IK_ENVELOPE_SPEC_V1.json').read_bytes()).hexdigest(),trials=rows,fresh_parity=json.loads((root/'BEST_FAST_FRESH_PARITY.json').read_text()),regular_backend_parity=json.loads((root/'BEST_NATIVE_PARITY.json').read_text()),runtime=dict(simulated_s=16.,neurons=fr['compiled']['neurons'],accelerated_complete_wall_s=fast['total_wall_s'],regular_complete_wall_s=regular['total_wall_s'],complete_speedup=regular['total_wall_s']/fast['total_wall_s'],accelerated_neural_controller_s=fr['compiled']['neural_step_wall_s'],regular_neural_controller_s=rr['compiled']['neural_step_wall_s'],neural_controller_speedup=rr['compiled']['neural_step_wall_s']/fr['compiled']['neural_step_wall_s'],scope='Matched full16s source-compliant simulations with source-Ib instrumentation. Neural/controller timing includes sensory adapters, not a neural-only kernel benchmark.'),pt_shadow_audit=dict(source=str(root/'pt_current_chain_v2/PT_CURRENT_CHAIN.json'),current_scaling_rows=[row for row in pt['rows'] if row['excitatory_gain_multiplier']==1. and row['command']==.05],gain_sweep_errors=pt['common_excitation_multiplier_mean_squared_errors'],live_trace_parity=json.loads((root/'PT_AUDIT_LIVE_TRAJECTORY_PARITY.json').read_text())['all_bit_exact']),frozen_source_files_unchanged=frozen,next_checkpoint='Re-fit source-consistent motor recruitment/afferent length normalization against the saved reference, particularly persistent unloaded dorsiflexion. Diagnose backward PT vertical cross-coupling separately; a global PT gain increase is unsupported.')
    (root/'CHECKPOINT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(completed_trials=result['completed_trials'],fully_accepted_trials=result['fully_accepted_trials'],runtime=result['runtime'],frozen_source_files_unchanged=frozen),indent=2));return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path);p.add_argument('--neural-folder',type=Path);p.add_argument('--pool-probe-folder',type=Path);p.add_argument('--ib-folder',type=Path);p.add_argument('--parity-left',type=Path);p.add_argument('--parity-right',type=Path);p.add_argument('--parity-out',type=Path);p.add_argument('--ib-checkpoint',type=Path);a=p.parse_args()
    if a.ib_checkpoint:ib_checkpoint(a.ib_checkpoint)
    elif a.parity_left and a.parity_right and a.parity_out:backend_parity(a.parity_left,a.parity_right,a.parity_out)
    elif a.ib_folder:ib_runtime_currents(a.ib_folder)
    elif a.pool_probe_folder:pool_bias_probe(a.pool_probe_folder)
    elif a.neural_folder:neural_currents(a.neural_folder)
    elif a.root:summarize(a.root)
    else:p.error('Provide --root or --neural-folder')
