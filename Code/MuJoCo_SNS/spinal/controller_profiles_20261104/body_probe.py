"""Fresh-process air curriculum probe of the frozen four-command Hunt fit.

Uses the existing physical evaluator and its mechanical repair stack. The
neural factory is composed explicitly; archived controller.py is unchanged.
Expected-replay and live-body afferents use identical declared normalization.
Only stages 1/2 are allowed; posture integration is a separate experiment.
"""
import argparse,hashlib,json,os,runpy,sys,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;S=HERE.parent;H=S/'hunt_training'
sys.path[:0]=[str(H),str(S),str(S/'functional_gait2392'),str(S/'w2l_cpg')]
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--reflex-fit',type=Path,required=True);p.add_argument('--stage',type=int,choices=(1,2),required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--expression',type=float,default=1.);p.add_argument('--feedback-scale',type=float,default=1.);p.add_argument('--local-ii-scale',type=float,default=1.);p.add_argument('--stim-cell');p.add_argument('--stim-nA',type=float,default=.5);p.add_argument('--stim-start',type=float,default=7.1);p.add_argument('--stim-duration',type=float,default=.1);p.add_argument('--noise-rg-sigma',type=float,default=0.);p.add_argument('--noise-seed',type=int,default=0);a=p.parse_args()
    if not .05<=a.expression<=1.:raise ValueError('Expression must be in [.05,1]')
    if not .01<=a.feedback_scale<=1.:raise ValueError('Feedback scale must be in [.01,1]')
    if not .01<=a.local_ii_scale<=1.:raise ValueError('Local II scale must be in [.01,1]')
    if not 0<=a.noise_rg_sigma<=.05:raise ValueError('RG current noise sigma must be in [0,.05] nA')
    if a.stim_cell and(not 0<=a.stim_start<=14 or not 0<a.stim_duration<=.5 or abs(a.stim_nA)>2):raise ValueError('Use a bounded declared stimulus')
    d_final=json.loads((a.reflex_fit/'FINAL.json').read_text());d_manifest=json.loads((a.reflex_fit/'MANIFEST.json').read_text())
    if not(d_final['accepted']and d_final['fresh_process_matched']):raise ValueError('Fresh-process Hunt D pass required')
    peaks=Path(d_manifest['fit']);c_fit=json.loads((peaks/'FINAL.json').read_text());c_manifest=json.loads((peaks/'MANIFEST.json').read_text());neural_cfg=c_manifest['configuration'];timing=Path(c_manifest['source_timing_block']);fixed=json.loads((timing/'FINAL.json').read_text())['parameters'];family=neural_cfg['family'];hunt_options=neural_cfg.get('hunt_options')or{}
    a.output.mkdir(parents=True,exist_ok=False)
    import controller as C,runner as R
    from controller_profiles_20261104 import profile as P
    import params
    original=C.build;busy=False
    target=np.load(neural_cfg['targets'],allow_pickle=False);mid=np.median(target['length'],axis=0);half=np.maximum(np.ptp(target['length'],axis=0)/2.,.02)
    Fmax=np.maximum(target['Fmax'],1.)
    def factory(acts,dt=.002,interleg=False,options=None):
        nonlocal busy
        if busy:return original(acts,dt=dt,interleg=interleg,options=options)
        if list(acts)!=target['act_names'].tolist():raise ValueError('Target/live actuator order differs')
        busy=True
        try:n=P.build(family,acts,dt=dt,options=neural_cfg.get('hunt_options'),config=neural_cfg['controller_profile'])
        finally:busy=False
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
        for i,(pre,post)in enumerate(zip(c.pre,c.post)):
            source,target_name=names[pre],names[post]
            if source.startswith('V3E_')or(source.startswith('RG_E_')and target_name.startswith('V3E_')):c.gmax[i]*=fixed.get('commissural_E',1.)
            elif source.startswith(('CIN_F_','V2a_','V0V_','V0D_','INI_'))or(source.startswith('RG_F_')and target_name.startswith(('CIN_F_','V2a_','V0D_'))):c.gmax[i]*=fixed.get('commissural_F',1.)
            if target_name.startswith('MN_')and source.startswith('PF_')and not source.startswith('PF_IN'):
                n.ff_gain[i]*=c_fit['best']['gains'].get(classify(target_name[3:]).base,1.)
            if not P.central(target_name):
                if source.startswith('II_'):
                    c.gmax[i]*=a.local_ii_scale;n.original_g[i]*=a.local_ii_scale
                    if i in n.reflex_edges:n.reflex_gain[i]*=a.local_ii_scale
                for aff in ('Ia','Ib'):
                    if source.startswith(aff+'_'):
                        gain=d_final['parameters'].get('local_'+aff,1.)
                        c.gmax[i]*=gain;n.original_g[i]*=gain
                        if i in n.reflex_edges:n.reflex_gain[i]*=gain
        n.apply_weights(1.)
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
        runtime.sensory_enabled=a.stage==2;runtime.contact_enabled=False
        disabled=np.array([i for name,i in runtime.ports.items()if name.startswith(('POST','BAL_','FUNC_'))],dtype=int)
        ia=np.array([runtime.ports['Ia_'+act]for act in acts]);ii=np.array([runtime.ports['II_'+act]for act in acts]);ib=np.array([runtime.ports['Ib_'+act]for act in acts])
        load_ids={side:np.array([i for i,act in enumerate(acts)if act.endswith('_'+side)and classify(act).groups[0]in('hip_ext','knee_ext','ankle_pf')])for side in ('l','r')}
        # Use the same reference load scale as the replay adapter.
        load_scale={side:max(float(np.percentile(np.maximum(target['force'][:,ids],0.).sum(axis=1),95)),1.)for side,ids in load_ids.items()}
        native=n.step;steps=0
        def advance(u):
            nonlocal steps
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
            for side,ids in load_ids.items():u[runtime.ports['LOAD_c_'+side]]=np.clip(force[ids].sum()/load_scale[side],0.,1.)
            for side in ('l','r'):
                runtime.commands['F_'+side]=0. if a.stage==1 else neural_cfg['command_F']
                runtime.commands['E_'+side]=0. if a.stage==1 else neural_cfg['command_E']
            if steps<50:runtime.commands['E_r']+=.4;runtime.commands['F_l']+=.4
            result=native(u);steps+=1;return result
        n.step=advance;return n
    C.build=factory
    base=json.loads(Path(r'D:\Users\Ben Bolen\NeuroAcceleration_20261104\body\compact_final.json').read_text())
    cfg={key:base[key]for key in ('force_law_parameters_file','knee_training_limits_deg','standing_solver')if key in base}
    trace=['RG_F_r','RG_F_l','MN_iliacus_r','MN_iliacus_l','MN_bifemsh_r','MN_bifemsh_l']
    trace+=['PF_S1_r','PF_S2_r','PF_S1_l','PF_S2_l']if family=='syn6'else['PF_HIP-F_r','PF_KNEE-F_r']if family in('hunt','w2lvar')else['PF_F1_r','PF_E1_r']
    cfg.update(stage=a.stage,backend='compact',controller_family=family,weights_directory=hunt_options.get('weights_directory',str(H/'runs/fitted')),balance=None,biarticular_pair_gain=hunt_options.get('biarticular_pair_gain',0.),learn_reciprocal=hunt_options.get('learn_reciprocal',False),rg_tonic_current=0.,axis_feedback=0.,timing_afferents=False,restriction=1.,support_fraction=0.,start_hip_deg=-7.,start_knee_deg=-10.,start_ankle_deg=0.,trace_cells=trace,verification_capture=str(a.output/'capture'),verification_sparse_currents=False,diagnostic_scope='Frozen four-command '+family+' fit, matched live afferent encoder, no direct posture/MN command; air curriculum probe')
    cfg['motor_scale']=a.expression
    cfg['central_RG_feedback_scale']=a.feedback_scale
    cfg['local_II_feedback_scale']=a.local_ii_scale
    if a.stim_cell:cfg['neural_stimulus']=dict(cell=a.stim_cell,current_nA=a.stim_nA,start_s=a.stim_start,duration_s=a.stim_duration,implementation='Temporary external bias-current addition; no voltage clamp')
    if a.noise_rg_sigma:cfg['neural_noise']=dict(targets=['RG_E_r','RG_F_r','RG_E_l','RG_F_l'],stationary_sigma_nA=a.noise_rg_sigma,tau_s=.02,seed=a.noise_seed,distribution='Independent zero-mean Ornstein-Uhlenbeck current per RG cell; starts at zero',scope='Fixed parameters; continuous current disturbance, not parameter perturbation')
    if family=='syn6':
        from verification_capture import requested_cells
        cfg['verification_cells']=requested_cells()+['PF_S'+str(i)+'_'+side for side in('r','l')for i in range(1,7)]
    elif family=='hunt':
        from verification_capture import requested_cells
        cfg['verification_cells']=requested_cells()+['PF_'+joint+'-'+phase+'_'+side for joint in('HIP','KNEE')for phase in('E','F')for side in('r','l')]
    if a.stim_cell:
        from verification_capture import requested_cells
        cfg['verification_cells']=list(dict.fromkeys(cfg.get('verification_cells',requested_cells())+[a.stim_cell]))
    cfg_path=a.output/'body.json';cfg_path.write_text(json.dumps(cfg,indent=2));(a.output/'SOURCE.json').write_text(json.dumps(dict(reflex_fit=str(a.reflex_fit),peak_fit=str(peaks),timing_fit=str(timing),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),frame='Existing converted Gait2392 world x anterior, y left, z superior; explicit legacy frame retained',command_scope='Normalized F/E commands, no sustained direct RG current; all POST/BAL/FUNC channels clamped off',air_expression_scale=a.expression,expression_scope='Shared PF->MN conductance scaling for suspended-air calibration; separate from the original ground-reference peak fit'),indent=2))
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
    try:runpy.run_path(str(H/'physical_eval.py'),run_name='__main__',init_globals={'print':compact_print})
    except Exception as exc:
        (a.output/'STOP.json').write_text(json.dumps(dict(error=repr(exc),knee_guard=guard_info,scope='Failed body probe; inspect capture/log'),indent=2));raise
    finally:mujoco.mj_step=native_physics;C.build=original;VC.AlignedCapture.sample=original_sample
if __name__=='__main__':main()
