"""Fresh-process air curriculum probe of the frozen four-command Hunt fit.

Uses the existing physical evaluator and its mechanical repair stack. The
neural factory is composed explicitly; archived controller.py is unchanged.
Expected-replay and live-body afferents use identical declared normalization.
Stages 3/4 are independent standing and supported-contact diagnostics.
"""
import argparse,hashlib,json,os,runpy,sys,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;S=HERE.parent;H=S/'hunt_training'
sys.path[:0]=[str(H),str(S),str(S/'functional_gait2392'),str(S/'w2l_cpg')]
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--reflex-fit',type=Path,required=True);p.add_argument('--stage',type=int,choices=(1,2,3,4),required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--expression',type=float,default=1.);p.add_argument('--feedback-scale',type=float,default=1.);p.add_argument('--local-ii-scale',type=float,default=1.);p.add_argument('--stim-cell');p.add_argument('--stim-nA',type=float,default=.5);p.add_argument('--stim-start',type=float,default=7.1);p.add_argument('--stim-duration',type=float,default=.1);p.add_argument('--noise-rg-sigma',type=float,default=0.);p.add_argument('--noise-seed',type=int,default=0);p.add_argument('--ankle-ia-scale',type=float,default=1.);p.add_argument('--ankle-ii-scale',type=float,default=1.);p.add_argument('--hip-expression',type=float,default=1.);p.add_argument('--audit-electrical',action='store_true');p.add_argument('--preserve-knee-crossing',action='store_true')
    p.add_argument('--lumbar-fit',type=Path)
    p.add_argument('--lumbar-kp',type=float,default=180.)
    p.add_argument('--lumbar-kd',type=float,default=18.)
    p.add_argument('--lumbar-feedback-scale',type=float,default=1.)
    p.add_argument('--lumbar-load-compensation',action='store_true')
    p.add_argument('--lumbar-pose-recruitment',action='store_true')
    p.add_argument('--lumbar-precharge',action='store_true')
    p.add_argument('--lumbar-phase-teacher',action='store_true')
    p.add_argument('--limb-onset-ramp',type=float,default=0.)
    p.add_argument('--torso-guide-scale',type=float,default=1.)
    p.add_argument('--native-force-law',action='store_true',
                   help='Keep the selected body model’s native actuator parameters')
    a=p.parse_args()
    if a.lumbar_fit and (a.stage not in (1,2) or not 0<=a.lumbar_kp<=800 or not 0<=a.lumbar_kd<=100):
        raise ValueError('Lumbar PF probes require an early stage and bounded gains')
    if a.lumbar_phase_teacher and (a.stage!=1 or not a.lumbar_fit):
        raise ValueError('The phase teacher is a separate Stage 1 diagnostic')
    if not 0. <= a.lumbar_feedback_scale <= 1. or (a.lumbar_feedback_scale!=1. and
            (not a.lumbar_fit or (a.stage!=2 and not a.lumbar_phase_teacher))):
        raise ValueError('Scale only a lumbar feedback controller, between zero and one')
    if a.lumbar_pose_recruitment and (not a.lumbar_fit or not a.lumbar_load_compensation
                                    or (a.stage!=2 and not a.lumbar_phase_teacher)):
        raise ValueError('Pose-dependent recruitment requires a lumbar feedback fit with load compensation')
    if not 0. <= a.limb_onset_ramp <= 3.:
        raise ValueError('Use a limb-onset ramp between zero and three seconds')
    native_spec = os.environ.get('AARL_SYN6_NATIVE_PF')
    start_knee = float(os.environ.get('AARL_AIR_START_KNEE_DEG', '-10'))
    if not -45. <= start_knee <= -5.:
        raise ValueError('Air probe initial knee must be flexed between 5 and 45 degrees')
    if native_spec and (a.expression != 1. or a.hip_expression != 1. or
                        a.ankle_ia_scale != 1. or a.ankle_ii_scale != 1. or
                        a.feedback_scale != 1. or a.local_ii_scale != 1.):
        raise ValueError('Native PF validation requires the frozen readout and sensory gains')
    if not .05<=a.expression<=1.:raise ValueError('Expression must be in [.05,1]')
    if not .05<=a.hip_expression<=1.:raise ValueError('Hip expression must be in [.05,1]')
    if not .01<=a.feedback_scale<=1.:raise ValueError('Feedback scale must be in [.01,1]')
    if not .01<=a.local_ii_scale<=1.:raise ValueError('Local II scale must be in [.01,1]')
    if not 0<=a.noise_rg_sigma<=.05:raise ValueError('RG current noise sigma must be in [0,.05] nA')
    if not all(.01<=g<=1. for g in(a.ankle_ia_scale,a.ankle_ii_scale)):raise ValueError('Ankle-crossing sensory scales must be in [.01,1]')
    if a.stim_cell and(not 0<=a.stim_start<=14 or not 0<a.stim_duration<=.5 or abs(a.stim_nA)>2):raise ValueError('Use a bounded declared stimulus')
    d_final=json.loads((a.reflex_fit/'FINAL.json').read_text());d_manifest=json.loads((a.reflex_fit/'MANIFEST.json').read_text())
    a.local_ii_scale*=d_final['parameters'].get('local_II',1.)
    if not(d_final['accepted']and d_final['fresh_process_matched']):raise ValueError('Fresh-process Hunt D pass required')
    peaks=Path(d_manifest['fit']);c_fit=json.loads((peaks/'FINAL.json').read_text());c_manifest=json.loads((peaks/'MANIFEST.json').read_text());neural_cfg=c_manifest['configuration'];timing=Path(c_manifest['source_timing_block']);fixed=json.loads((timing/'FINAL.json').read_text())['parameters'];family=neural_cfg['family'];hunt_options=neural_cfg.get('hunt_options')or{}
    if native_spec:
        neural_cfg['controller_profile']['native_pf'] = json.loads(Path(native_spec).read_text())
    if a.lumbar_fit:
        neural_cfg['controller_profile']['lumbar_pf'] = json.loads(a.lumbar_fit.read_text())
    if a.stage >= 3 and (family != 'syn6' or not native_spec):
        raise ValueError('Ground diagnostics require an explicit native syn6 specification')
    a.output.mkdir(parents=True,exist_ok=False)
    import controller as C,runner as R
    from controller_profiles_20261104 import profile as P
    import params
    original=C.build;busy=False
    lumbar_instance=None
    electrical={}
    cadence = None
    target=np.load(neural_cfg['targets'],allow_pickle=False);mid=np.median(target['length'],axis=0);half=np.maximum(np.ptp(target['length'],axis=0)/2.,.02)
    Fmax=np.maximum(target['Fmax'],1.)
    def factory(acts,dt=.002,interleg=False,options=None):
        nonlocal busy, cadence, lumbar_instance
        if busy:return original(acts,dt=dt,interleg=interleg,options=options)
        if list(acts)!=target['act_names'].tolist():raise ValueError('Target/live actuator order differs')
        busy=True
        try:n=P.build(family,acts,dt=dt,options=neural_cfg.get('hunt_options'),config=neural_cfg['controller_profile'])
        finally:busy=False
        cadence = getattr(n, 'native_pf_cadence', None)
        # The historical evaluator selects its PF logger using this selector.
        # Bind it to the network actually returned, rather than aliasing cells.
        os.environ['AARL_NET']='w2lvar'if family=='hunt'else family
        c=n.compiled;names=P.nodes(n);rg=np.array([n.idx['RG_'+phase+'_'+side]for side in ('l','r')for phase in ('E','F')])
        if family!='hunt':
            n.names=names;n.ports=dict(n.command_runtime.ports);n.fast_io=True;n.live=None
            # Compatibility names used by the historical evaluator. These are
            # aliases of existing command ports, not extra cells or currents;
            # Runtime.prepare_inputs writes the declared commands afterward.
            for side in ('r','l'):
                for phase in ('E','F'):n.ports['HUNT_INIT_'+phase+'_'+side]=n.ports['RG_'+phase+'_drive_'+side]
            n.original_g=c.gmax.copy();n.ff_gain=c.gmax.copy();n.reflex_gain=c.gmax.copy()
            n.mn_edges=np.array([i for i,p in enumerate(c.post)if names[p].startswith('MN_')],dtype=int)
            n.pf_edges=np.array([i for i in n.mn_edges if names[c.pre[i]].startswith('PF_')and not names[c.pre[i]].startswith('PF_IN')],dtype=int)
            n.reflex_edges=np.array([i for i in n.mn_edges if names[c.pre[i]].startswith(('Ia_','II_','Ib_','IaIN_','IBIN_','IBEXC_','IIX_','IIIN_'))],dtype=int)
            n.stand_edges=np.array([],dtype=int)
            def sync():c.c.g_max_non[c.post,c.pre]=c.gmax
            def apply_weights(gate):
                c.gmax[n.pf_edges]=n.ff_gain[n.pf_edges]*gate
                c.gmax[n.reflex_edges]=n.reflex_gain[n.reflex_edges]
                sync()
            def bind(m,d):n.live=(m,d)
            n.sync=sync;n.apply_weights=apply_weights;n.bind_mujoco=bind
            def static_advance(u):
                n.apply_weights(1.);c.forward(u);return c.V
            n.command_runtime.original_step=static_advance
        c.c.tau_max_b[:,rg]*=fixed.get('tau_scale',1.);c.c.g_ion[:,rg]*=fixed.get('g_ion_scale',1.);c.refresh_parameters()
        from muscle_map import classify
        from controller_profiles_20261104.readout_keys import multiplier
        for i,(pre,post)in enumerate(zip(c.pre,c.post)):
            source,target_name=names[pre],names[post]
            if source.startswith('V3E_')or(source.startswith('RG_E_')and target_name.startswith('V3E_')):c.gmax[i]*=fixed.get('commissural_E',1.)
            elif source.startswith(('CIN_F_','V2a_','V0V_','V0D_','INI_'))or(source.startswith('RG_F_')and target_name.startswith(('CIN_F_','V2a_','V0D_'))):c.gmax[i]*=fixed.get('commissural_F',1.)
            if target_name.startswith('MN_')and source.startswith('PF_')and not source.startswith('PF_IN'):
                n.ff_gain[i]*=multiplier(source,target_name,c_fit)
                if classify(target_name[3:]).groups[0]in('hip_ext','hip_flex'):n.ff_gain[i]*=a.hip_expression
            if not P.central(target_name):
                for aff,scale in(('Ia',a.ankle_ia_scale),('II',a.ankle_ii_scale)):
                    if scale!=1. and source.startswith(aff+'_'):
                        info=classify(source[len(aff)+1:])
                        if info and any(group in('ankle_pf','ankle_df')for group in info.groups)and not(a.preserve_knee_crossing and any(group.startswith('knee_')for group in info.groups)):
                            c.gmax[i]*=scale;n.original_g[i]*=scale
                            if i in n.reflex_edges:n.reflex_gain[i]*=scale
                if source.startswith('II_'):
                    c.gmax[i]*=a.local_ii_scale;n.original_g[i]*=a.local_ii_scale
                    if i in n.reflex_edges:n.reflex_gain[i]*=a.local_ii_scale
                for aff in ('Ia','Ib'):
                    if source.startswith(aff+'_'):
                        gain=d_final['parameters'].get('local_'+aff,1.)
                        c.gmax[i]*=gain;n.original_g[i]*=gain
                        if i in n.reflex_edges:n.reflex_gain[i]*=gain
            if target_name.startswith('MN_'):
                key='relay_Ib_output'if source.startswith(('IBIN_','IBEXC_'))else'relay_II_output'if source.startswith(('IIX_','IIIN_'))else None
                if key:
                    gain=d_final['parameters'].get(key,1.)
                    c.gmax[i]*=gain;n.original_g[i]*=gain
                    if i in n.reflex_edges:n.reflex_gain[i]*=gain
        lumbar = None
        if a.lumbar_fit:
            from controller_profiles_20261104.lumbar_pf import Runtime
            lumbar = Runtime(n,a.stage,a.lumbar_kp,a.lumbar_kd,a.lumbar_load_compensation,
                             a.lumbar_precharge,a.lumbar_phase_teacher,
                             pose_recruitment=a.lumbar_pose_recruitment,
                             feedback_scale=a.lumbar_feedback_scale)
            n.lumbar_runtime=lumbar
            lumbar_instance=lumbar
            native_weights = n.apply_weights
            def apply_with_lumbar(gate):
                native_weights(gate)
                lumbar.apply_weights()
            n.apply_weights = apply_with_lumbar
        n.apply_weights(1.)
        if a.limb_onset_ramp:
            limb_edges = np.array([i for i in n.pf_edges if not names[c.post[i]][3:].startswith(('ercspn','intobl','extobl'))],dtype=int)
            limb_gains = n.ff_gain[limb_edges].copy()
        if family == 'syn6':
            from controller_profiles_20261104.native_pf import motor_digest
            (a.output / 'MOTOR_SHA256.txt').write_text(motor_digest(c, names) + '\n')
        n.session_stimulus=np.zeros(len(c.V))
        stimulus_index=n.idx[a.stim_cell]if a.stim_cell else None
        stimulus_bias=float(c.c.i_b[stimulus_index])if stimulus_index is not None else None
        noise_rng=np.random.default_rng(a.noise_seed)
        noise_ids=rg.copy();noise_bias=c.c.i_b[noise_ids].copy()
        noise_state=np.zeros(len(noise_ids));noise_draw=np.zeros(len(noise_ids))
        noise_decay=np.exp(-dt/.02)
        noise_scale=a.noise_rg_sigma*np.sqrt(1.-noise_decay**2)
        runtime=n.command_runtime
        for key in runtime.shared_gains:runtime.shared_gains[key]=fixed.get(key,1.)*(a.feedback_scale if key.endswith('_RG')else 1.)
        runtime.sensory_enabled=a.stage>=2;runtime.contact_enabled=a.stage>=3
        disabled=np.array([i for name,i in runtime.ports.items()if name.startswith(('POST','BAL_','FUNC_'))],dtype=int)
        ia=np.array([runtime.ports['Ia_'+act]for act in acts]);ii=np.array([runtime.ports['II_'+act]for act in acts]);ib=np.array([runtime.ports['Ib_'+act]for act in acts])
        load_ids={side:np.array([i for i,act in enumerate(acts)if act.endswith('_'+side)and classify(act).groups[0]in('hip_ext','knee_ext','ankle_pf')])for side in ('l','r')}
        # Use the same reference load scale as the replay adapter.
        load_scale={side:max(float(np.percentile(np.maximum(target['force'][:,ids],0.).sum(axis=1),95)),1.)for side,ids in load_ids.items()}
        native=n.step;steps=0
        protection = runtime.config.get('native_pf', {}).get('knee_premotor')
        prediction_s = float(protection.get('prediction_s', .065)) if protection else 0.
        if protection and not .02 <= prediction_s <= .15:
            raise ValueError('Knee premotor prediction must be between 20 and 150 ms')
        knee_sensors = None
        ratio_max=np.zeros(len(c.V));voltage_min=np.full(len(c.V),np.inf);voltage_max=-voltage_min.copy();ratio_over2=np.zeros(len(c.V),dtype=int);electrical_samples=0
        if a.audit_electrical:electrical.update(neuron_names=names,maximum_Euler_load_ratio=ratio_max,steps_above2=ratio_over2,minimum_raw_voltage_mV=voltage_min,maximum_raw_voltage_mV=voltage_max,samples=0,scope='Actual live-body synaptic-load Euler ratio, excluding ionic conductance; not a full recurrent stability proof')
        def advance(u):
            nonlocal steps,electrical_samples,knee_sensors
            m,d=n.live
            if stimulus_index is not None:
                amplitude=a.stim_nA if a.stim_start<=steps*dt<a.stim_start+a.stim_duration else 0.
                n.session_stimulus[stimulus_index]=amplitude;c.c.i_b[stimulus_index]=stimulus_bias+amplitude
            if a.noise_rg_sigma:
                noise_rng.standard_normal(out=noise_draw)
                noise_state[:]=noise_decay*noise_state+noise_scale*noise_draw
                c.c.i_b[noise_ids]=noise_bias+noise_state
                n.session_stimulus[noise_ids]=noise_state
                if stimulus_index is not None and stimulus_index in noise_ids:
                    c.c.i_b[stimulus_index]+=amplitude;n.session_stimulus[stimulus_index]+=amplitude
            u[disabled]=0.
            u[ia]=params.AFF['ia_gain']*np.clip(d.actuator_velocity/.6,-1.,1.)
            u[ii]=params.AFF['ii_gain']*np.clip((d.actuator_length-mid)/half,-1.,1.)
            force=np.maximum(-d.actuator_force,0.)
            u[ib]=params.AFF['ib_gain']*force/Fmax
            if a.stage < 3:
                for side,ids in load_ids.items():u[runtime.ports['LOAD_c_'+side]]=np.clip(force[ids].sum()/load_scale[side],0.,1.)
            if protection:
                if knee_sensors is None:
                    knee_sensors = [(runtime.ports['KNEE_PREMOTOR_' + side],
                                     int(m.joint('knee_angle_' + side).qposadr[0]),
                                     int(m.joint('knee_angle_' + side).dofadr[0]))
                                    for side in ('l', 'r')]
                for port, q, v in knee_sensors:
                    velocity = max(float(d.qvel[v]), 0.)
                    projected = float(d.qpos[q]) + prediction_s * velocity
                    proximity = min(1., max(0., (projected + .40) / .30))
                    u[port] = 5. * proximity * min(1., velocity / 3.)
            for side in ('l','r'):
                runtime.commands['F_'+side]=0. if a.stage in (1,3) else neural_cfg['command_F']
                runtime.commands['E_'+side]=0. if a.stage in (1,3) else neural_cfg['command_E']
            if a.limb_onset_ramp:
                progress=min(1.,d.time/a.limb_onset_ramp)
                expression=.5-.5*np.cos(np.pi*progress)
                n.ff_gain[limb_edges]=expression*limb_gains
                for key in runtime.commands:
                    runtime.commands[key]*=expression
            if steps<50 and a.stage!=3:runtime.commands['E_r']+=.4;runtime.commands['F_l']+=.4
            if lumbar is not None:
                lumbar.update(u,m,d)
            result=native(u)
            if a.audit_electrical and 5.<=steps*dt<14.:
                conductance=c.gmax*np.clip((c.last_voltage[c.pre]-c.elo)/c.erange,0.,1.);ratio=c.c.time_factor_membrane*(c.c.g_m+np.bincount(c.post,weights=conductance,minlength=len(c.V)));np.maximum(ratio_max,ratio,out=ratio_max);ratio_over2[:]+=ratio>2.;np.minimum(voltage_min,c.V,out=voltage_min);np.maximum(voltage_max,c.V,out=voltage_max);electrical_samples+=1
                electrical['samples']=electrical_samples
            steps+=1;return result
        n.step=advance;return n
    C.build=factory
    base=json.loads(Path(r'D:\Users\Ben Bolen\NeuroAcceleration_20261104\body\compact_final.json').read_text())
    cfg={key:base[key]for key in ('force_law_parameters_file','knee_training_limits_deg','standing_solver')if key in base}
    if a.native_force_law:
        cfg.pop('force_law_parameters_file',None)
    trace=['RG_F_r','RG_F_l','MN_iliacus_r','MN_iliacus_l','MN_bifemsh_r','MN_bifemsh_l']
    trace+=['PF_S1_r','PF_S2_r','PF_S1_l','PF_S2_l']if family=='syn6'else['PF_HIP-F_r','PF_KNEE-F_r']if family in('hunt','w2lvar')else['PF_F1_r','PF_E1_r']
    cfg.update(stage=a.stage,backend='compact',controller_family=family,weights_directory=hunt_options.get('weights_directory',str(H/'runs/fitted')),balance=None,biarticular_pair_gain=hunt_options.get('biarticular_pair_gain',0.),learn_reciprocal=hunt_options.get('learn_reciprocal',False),rg_tonic_current=0.,axis_feedback=0.,timing_afferents=False,restriction=1.,support_fraction=0.,start_hip_deg=-7.,start_knee_deg=-10.,start_ankle_deg=0.,trace_cells=trace,verification_capture=str(a.output/'capture'),verification_sparse_currents=False,diagnostic_scope='Frozen four-command '+family+' fit, matched live afferent encoder, no direct posture/MN command; air curriculum probe')
    cfg['motor_scale']=a.expression
    if a.lumbar_fit:
        cfg['lumbar_pf_fit']=str(a.lumbar_fit.resolve())
        cfg['lumbar_pf_fit_sha256']=hashlib.sha256(a.lumbar_fit.read_bytes()).hexdigest()
        feedback=a.stage==2 or a.lumbar_phase_teacher
        cfg['lumbar_feedback_gains']=[a.lumbar_kp,a.lumbar_kd] if feedback else [0.,0.]
        cfg['lumbar_feedback_scale']=a.lumbar_feedback_scale if feedback else 0.
        cfg['lumbar_load_compensation']=bool(feedback and a.lumbar_load_compensation)
        cfg['lumbar_pose_recruitment']=a.lumbar_pose_recruitment
        cfg['lumbar_precharge']=a.lumbar_precharge
        cfg['lumbar_phase_teacher']=a.lumbar_phase_teacher
    cfg['torso_guide_scale']=a.torso_guide_scale
    if a.limb_onset_ramp:
        cfg['limb_onset_ramp_s']=a.limb_onset_ramp
    cfg['hip_PF_to_MN_expression']=a.hip_expression
    cfg['central_RG_feedback_scale']=a.feedback_scale
    cfg['local_II_feedback_scale']=a.local_ii_scale
    cfg['ankle_crossing_local_Ia_scale']=a.ankle_ia_scale
    cfg['ankle_crossing_local_II_scale']=a.ankle_ii_scale
    cfg['preserve_knee_crossing_afferents']=a.preserve_knee_crossing
    cfg['local_relay_output_gains']={key:d_final['parameters'].get(key,1.)for key in('relay_Ib_output','relay_II_output')}
    if a.stim_cell:cfg['neural_stimulus']=dict(cell=a.stim_cell,current_nA=a.stim_nA,start_s=a.stim_start,duration_s=a.stim_duration,implementation='Temporary external bias-current addition; no voltage clamp')
    if a.noise_rg_sigma:cfg['neural_noise']=dict(targets=['RG_E_r','RG_F_r','RG_E_l','RG_F_l'],stationary_sigma_nA=a.noise_rg_sigma,tau_s=.02,seed=a.noise_seed,distribution='Independent zero-mean Ornstein-Uhlenbeck current per RG cell; starts at zero',scope='Fixed parameters; continuous current disturbance, not parameter perturbation')
    if family=='syn6':
        from verification_capture import requested_cells
        cfg['verification_cells']=requested_cells()+['PF_S'+str(i)+'_'+side for side in('r','l')for i in range(1,7)]
        if neural_cfg['controller_profile'].get('native_pf', {}).get('knee_premotor'):
            cfg['verification_cells'] += ['TimingKneePremotor_l', 'TimingKneePremotor_r']
    elif family=='hunt':
        from verification_capture import requested_cells
        cfg['verification_cells']=requested_cells()+['PF_'+joint+'-'+phase+'_'+side for joint in('HIP','KNEE')for phase in('E','F')for side in('r','l')]
    if a.stim_cell:
        from verification_capture import requested_cells
        cfg['verification_cells']=list(dict.fromkeys(cfg.get('verification_cells',requested_cells())+[a.stim_cell]))
    if a.lumbar_fit:
        from verification_capture import requested_cells
        from controller_profiles_20261104.lumbar_pf import LABELS, ROUTES
        prefix='PF_FUNC_' if family=='hunt' else 'PF_'
        lumbar_cells=[f'{prefix}{label}_{direction}_mid' for label in LABELS for direction in ('P','N')]
        cfg['verification_cells']=list(dict.fromkeys(cfg.get('verification_cells',requested_cells())+
            lumbar_cells+['MN_'+route for route in ROUTES]))
    if native_spec:
        cfg['native_pf_spec'] = native_spec
    cfg['start_knee_deg'] = start_knee
    if a.stage >= 3:
        cfg.update(standing_solver='contact_symmetric', regional_contacts=True,
                   support_fraction=.65 if a.stage == 4 else 0.,
                   diagnostic_scope='Frozen syn6 ground diagnostic with measured hindfoot and load feedback')
    cfg_path=a.output/'body.json';cfg_path.write_text(json.dumps(cfg,indent=2))
    provenance=dict(reflex_fit=str(a.reflex_fit),peak_fit=str(peaks),timing_fit=str(timing),
        model_xml=str(Path(R.MODEL).resolve()),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        frame='Existing converted Gait2392 world x anterior, y left, z superior; explicit legacy frame retained',
        command_scope='Normalized F/E commands, no sustained direct RG current; inherited POST/BAL/FUNC inputs clamped off',
        air_expression_scale=a.expression,
        expression_scope='Shared PF->MN conductance scaling for suspended-air calibration; separate from the original ground-reference peak fit')
    if a.lumbar_fit:
        provenance['lumbar_pf_implementation_sha256']=hashlib.sha256((HERE/'lumbar_pf.py').read_bytes()).hexdigest()
    (a.output/'SOURCE.json').write_text(json.dumps(provenance,indent=2))
    sys.argv=[str(H/'physical_eval.py'),str(cfg_path)]
    def compact_print(*args,**kwargs):
        if args and isinstance(args[0],str)and len(args[0])>10000:print('Detailed structured body result saved beside configuration.')
        else:print(*args,**kwargs)
    import mujoco
    import verification_capture as VC
    original_sample=VC.AlignedCapture.sample
    if a.stim_cell or a.noise_rg_sigma:
        def sample_with_stimulus(recorder,data,*args,**kwargs):
            kwargs['stimulus_nA']=recorder.network.session_stimulus
            return original_sample(recorder,data,*args,**kwargs)
        VC.AlignedCapture.sample=sample_with_stimulus
    native_physics=mujoco.mj_step;guard_info={}
    def protected_physics(m,d,*args,**kwargs):
        for side in ('r','l'):
            j=m.joint('knee_angle_'+side);q=float(d.qpos[j.qposadr[0]]);v=float(d.qvel[j.dofadr[0]])
            if v>0. and q+m.opt.timestep*v>=np.radians(-.1):
                guard_info.update(side=side,time_s=float(d.time),angle_deg=float(np.degrees(q)),velocity_deg_s=float(np.degrees(v)),rule='Stop before predicted extension reaches -0.1 degree; no state clipping')
                raise RuntimeError('Predictive knee guard')
        return native_physics(m,d,*args,**kwargs)
    mujoco.mj_step=protected_physics
    try:
        runpy.run_path(str(H/'physical_eval.py'),run_name='__main__',init_globals={'print':compact_print})
        from kinematic_response import score_recording
        response,_,_=score_recording(a.output)
        (a.output/'body.kinematics.json').write_text(json.dumps(response,indent=2,allow_nan=False))
    except Exception as exc:
        (a.output/'STOP.json').write_text(json.dumps(dict(error=repr(exc),knee_guard=guard_info,scope='Failed body probe; inspect capture/log'),indent=2));raise
    finally:
        mujoco.mj_step=native_physics;C.build=original;VC.AlignedCapture.sample=original_sample
        if lumbar_instance is not None:
            np.savez_compressed(a.output/'lumbar_pf_trace.npz',samples=np.asarray(lumbar_instance.history),
                                columns=['time_s','extension_command','bending_command','rotation_command',
                                         'RG_E_r','RG_F_r','RG_E_l','RG_F_l','extension_rad',
                                         'bending_rad','rotation_rad'])
        if cadence is not None:
            (a.output/'CADENCE.json').write_text(json.dumps(dict(
                period_s=cadence.period.tolist(), measured=cadence.measured,
                native_fraction=cadence.mix, full_native_since_neural_s=cadence.full_native_since_s,
                neural_time_s=cadence.time), indent=2))
        if a.audit_electrical:
            count=electrical.get('samples',0);saved={key:(value.tolist()if count else None)if isinstance(value,np.ndarray)else value for key,value in electrical.items()}
            if count:saved['fraction_above2']=(electrical['steps_above2']/count).tolist()
            (a.output/'ELECTRICAL_AUDIT.json').write_text(json.dumps(saved,indent=2))
if __name__=='__main__':main()
