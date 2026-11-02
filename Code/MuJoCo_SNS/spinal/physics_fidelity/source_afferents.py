"""Isolated source-fiber spindle and complete ankle MN-current diagnostics."""
import os,json,hashlib
from pathlib import Path
import numpy as np
from thelen_dynamic import parameters_from_osim

CATEGORIES=('PF_feedforward','Ia_homonymous','Ia_reciprocal','II_direct','Ib_inhibitory','Ib_pool','PT_excitatory','PT_inhibitory','inherited_balance','posture_synaptic','other_synaptic','direct_POST','other_direct','leak','bias')

def category(name):
    for prefix,label in [('PF_','PF_feedforward'),('Ia_','Ia_homonymous'),('IaIN_','Ia_reciprocal'),('II_','II_direct'),('IBIN_','Ib_inhibitory'),('IBEXC_','Ib_pool'),('PT_BAL_','PT_excitatory'),('PT_IN_','PT_inhibitory'),('BAL_','inherited_balance'),('POST','posture_synaptic')]:
        if name.startswith(prefix):return CATEGORIES.index(label)
    return CATEGORIES.index('other_synaptic')

class Session:
    def __init__(self):self.networks=[];self.backends={}
    def bind_backend(self,m,d,backend):
        if backend is None:raise ValueError('Source-fiber diagnostics require source physics')
        self.backends[(id(m),id(d))]=backend
    def save(self,out):
        records=[]
        for i,n in enumerate(self.networks):
            suffix='' if i==0 else str(i);k=n.source_audit_count
            np.savez_compressed(out/('source_afferents'+suffix+'.npz'),neural_t_s=n.source_audit_t[:k],kinematics=n.source_audit_kin[:k],kinematic_fields=['MTU_length_m','MTU_velocity_mps','fiber_length_m','fiber_velocity_mps'],source_fiber_ready=n.source_audit_ready[:k],actuators=list(n.muscles),selected_actuators=n.source_audit_acts,raw_local_inputs=n.source_audit_inputs[:k],input_fields=['Ia','II'],mn_currents=n.source_audit_currents[:k],current_categories=CATEGORIES,voltage=n.source_audit_voltage[:k],command=n.source_audit_command[:k],delivered_activation=n.source_audit_activation[:k],PF_voltage=n.source_audit_pf[:k],PF_cells=n.source_audit_pf_names)
            records.append(n.source_audit_metadata)
        (out/'SOURCE_AFFERENT_AUDIT.json').write_text(json.dumps(records,indent=2));return records

def install(options,physical_source,stage):
    import controller as C
    allowed={'format','mode','functional_roles_file','motor_functional_scaling','silence_dorsiflexor_Ia','local_encoder'}
    if options.get('format')!='aarl_source_afferent_v1' or options.get('mode') not in ('audit','calibrate') or set(options)-allowed:raise ValueError('Explicit source-afferent options required')
    scaling=options.get('motor_functional_scaling',{'ankle_df':1.,'ankle_pf':1.})
    if set(scaling)!={'ankle_df','ankle_pf'} or any(not np.isfinite(v) or not .1<=v<=3. for v in scaling.values()):raise ValueError('Declared shared ankle motor recruitment multipliers in[.1,3] only')
    encoder=options.get('local_encoder','legacy')
    if encoder not in ('legacy','source_fiber'):raise ValueError('Unknown local spindle encoder')
    if options['mode']=='audit' and (any(v!=1. for v in scaling.values()) or options.get('silence_dorsiflexor_Ia') or encoder!='legacy'):raise ValueError('Observational audit must not alter the controller')
    roles=json.loads(Path(options['functional_roles_file']).read_text())['roles']['ankle_pf'] if options.get('functional_roles_file') else None
    if roles is None and (any(v!=1. for v in scaling.values()) or options.get('silence_dorsiflexor_Ia') or encoder!='legacy'):raise ValueError('Frozen equality-aware functional roles required')
    session=Session();original=C.build
    def build(*args,**kwargs):
        n=original(*args,**kwargs);c=n.compiled;names=n.names;acts=list(n.muscles);p=parameters_from_osim(physical_source,acts,0.)
        if encoder=='source_fiber' and C.P.G.get('gamma_coact',0.)!=0.:raise ValueError('Source-fiber conversion currently validated for fixed fusimotor gain only')
        function={}
        if roles is not None:
            for a in acts:
                mirror=a[:-2]+('_l' if a.endswith('_r') else '_r');both=[roles[a],roles[mirror]]
                if all(row['min_moment_arm_m']>.001 for row in both):function[a]='ankle_df'
                elif all(row['max_moment_arm_m']<-.001 for row in both):function[a]='ankle_pf'
            for edge in n.pf_edges:
                a=names[c.post[edge]][3:]
                if a in function:n.ff_gain[edge]*=float(scaling[function[a]])
        families=('tib_ant','ext_dig','ext_hal','per_tert','soleus','med_gas','lat_gas','per_long','per_brev','tib_post','flex_dig','flex_hal','vas_lat','bifemsh','rect_fem','semimem')
        selected=[i for i,a in enumerate(acts) if a[:-2] in families];selected_acts=[acts[i] for i in selected];mn=np.array([n.idx[n.mn_names[a]] for a in selected_acts]);rows={int(node):i for i,node in enumerate(mn)};postports=np.array([n.ports['POST_'+a] for a in selected_acts]);inputmatrix=c.c.input_connectivity[mn].copy();postweights=inputmatrix[np.arange(len(mn)),postports]
        edges=np.array([i for i,post in enumerate(c.post) if int(post) in rows],dtype=int);slots=np.array([rows[int(c.post[i])]*len(CATEGORIES)+category(names[c.pre[i]]) for i in edges]);ia=np.array([n.ports['Ia_'+a] for a in selected_acts]);ii=np.array([n.ports['II_'+a] for a in selected_acts]);pf_names=[name for name in names if name.startswith('PF_FUNC_ANKLE_Z_') or name.startswith('RG_') or name=='PF_HUNT_TONIC'];pf=np.array([n.idx[name] for name in pf_names])
        size=8500;n.source_audit_count=0;n.source_audit_t=np.zeros(size);n.source_audit_kin=np.full((size,4,len(acts)),np.nan);n.source_audit_ready=np.zeros(size,dtype=bool);n.source_audit_currents=np.zeros((size,len(mn),len(CATEGORIES)));n.source_audit_inputs=np.zeros((size,2,len(mn)));n.source_audit_voltage=np.zeros((size,len(mn)));n.source_audit_command=np.zeros_like(n.source_audit_voltage);n.source_audit_activation=np.zeros_like(n.source_audit_voltage);n.source_audit_pf=np.zeros((size,len(pf)));n.source_audit_acts=selected_acts;n.source_audit_pf_names=pf_names
        if np.any(c.c.g_ion[:,mn]!=0.):raise RuntimeError('Current accounting expects non-gated MNs')
        n.source_audit_metadata=dict(format='aarl_source_afferent_audit_v1',mode='audit',stage=stage,module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(physical_source).read_bytes()).hexdigest(),source_dimensions=p.tolist(),selected_actuators=selected_acts,direct_POST_input_weights=postweights.tolist(),maximum_MN_update_reconstruction_error_mV=0.,source_fiber_scope='Actual committed source fiber state and last-step velocity after source backend starts. Native pose-held initialization is retained; unavailable fiber rows are explicitly marked.',encoding_scope='Observational only: existing Ia/II and Hunt timing definitions unchanged.',sparse_forward_required=True)
        n.source_audit_metadata.update(mode=options['mode'],functional_motor_scaling=scaling,verified_functional_actuators=function,silence_dorsiflexor_Ia=bool(options.get('silence_dorsiflexor_Ia')),local_encoder=encoder,encoding_scope='Hunt B timing remains frozen. Optional local spindle conversion uses strain(Lf/l0-1)/0.2 and velocity Lfdot/(Vmax*l0), with original source dimensions;0.2 is a declared engineering transduction reference,not a measured spindle constant. Existing phase,speed,andtonic current gains retained.')
        allports={prefix:np.array([n.ports[prefix+'_'+a] for a in acts]) for prefix in ('Ia','II')};df=np.array([i for i,a in enumerate(acts) if function.get(a)=='ankle_df'],dtype=int);ankles=np.array([i for i,a in enumerate(acts) if a in function],dtype=int)
        disabled=set(os.environ.get('AARL_DEAFF','').split(','))-{''};off=np.array([i for i,a in enumerate(acts) if n.muscles[a].side in disabled],dtype=int)
        native_step=n.step
        def step(u):
            k=n.source_audit_count;m,d=n.live
            if k>=size:raise RuntimeError('Audit preallocation exceeded')
            n.source_audit_t[k]=k*C.P.DT;n.source_audit_kin[k,0]=d.actuator_length;n.source_audit_kin[k,1]=d.actuator_velocity;n.source_audit_inputs[k]=[u[ia],u[ii]];n.source_audit_activation[k]=d.act[selected]
            backend=session.backends.get((id(m),id(d)))
            if backend is not None:
                s=backend.bank.state;lf=np.hypot(s,p[:,3]);n.source_audit_kin[k,2]=lf;n.source_audit_kin[k,3]=(s-backend.old_state)*s/lf/backend.dt;n.source_audit_ready[k]=True
                if encoder=='source_fiber' and stage!=1:
                    from muscle_map import EXTENSOR_STANCE_GROUPS
                    stance=np.array([np.clip(c.V[n.idx['RG_E_'+n.muscles[a].side]]/C.P.E_HI,0.,1.) for a in acts]);ext=np.array([n.muscles[a].groups[0] in EXTENSOR_STANCE_GROUPS for a in acts]);drive=np.clip(c.V[n.idx['DRIVE']]/C.P.E_HI,0.,1.)
                    ia_gain=C.P.AFF['ia_gain']*(1+C.P.MOD['ia']*drive)*np.where(ext,.5+.5*stance,1.);ii_gain=C.P.AFF['ii_gain']*(1+C.P.MOD['ii']*drive)*np.where(ext,.3+.7*stance,.5)
                    velocity=np.clip(n.source_audit_kin[k,3]/(p[:,1]*p[:,8]),-1.,1.);length=np.clip((lf/p[:,1]-1.)/.2,-1.,1.)
                    u[allports['Ia'][ankles]]=np.maximum(ia_gain[ankles]*velocity[ankles],0.);u[allports['II'][ankles]]=np.maximum(C.P.AFF['i0_ii']*(.3+.7*stance[ankles])+ii_gain[ankles]*length[ankles],0.)
            if options.get('silence_dorsiflexor_Ia') and stage!=1:u[allports['Ia'][df]]=0.
            if len(off):u[allports['Ia'][off]]=0.;u[allports['II'][off]]=0.
            n.source_audit_inputs[k]=[u[ia],u[ii]]
            return native_step(u)
        n.step=step;native_forward=c.forward
        def forward(u=None):
            if u is None:raise RuntimeError('Explicit neural input vector required')
            k=n.source_audit_count;old=c.V[mn].copy();current=c.gmax[edges]*np.clip((c.V[c.pre[edges]]-c.elo[edges])/c.erange[edges],0.,1.)*(c.reversal[edges]-c.V[c.post[edges]])
            n.source_audit_currents[k]=np.bincount(slots,weights=current,minlength=len(mn)*len(CATEGORIES)).reshape(len(mn),len(CATEGORIES));post=postweights*u[postports];direct=inputmatrix@u
            n.source_audit_currents[k,:,CATEGORIES.index('direct_POST')]=post;n.source_audit_currents[k,:,CATEGORIES.index('other_direct')]=direct-post;n.source_audit_currents[k,:,CATEGORIES.index('leak')]=-c.c.g_m[mn]*(old-c.c.V_rest[mn]);n.source_audit_currents[k,:,CATEGORIES.index('bias')]=c.c.i_b[mn]
            predicted=old+c.c.time_factor_membrane[mn]*n.source_audit_currents[k].sum(axis=1);result=native_forward(u);error=float(np.max(abs(predicted-c.V[mn])));n.source_audit_metadata['maximum_MN_update_reconstruction_error_mV']=max(n.source_audit_metadata['maximum_MN_update_reconstruction_error_mV'],error)
            if error>1e-10:raise RuntimeError('Incomplete MN-current accounting')
            n.source_audit_voltage[k]=c.V[mn];n.source_audit_command[k]=np.clip(c.V[mn]/C.P.E_HI,0.,1.);n.source_audit_pf[k]=c.V[pf];n.source_audit_count+=1;return result
        c.forward=forward;session.networks.append(n);return n
    C.build=build;return session

def summarize(folder):
    archive=np.load(folder/'source_afferents.npz');z={key:archive[key] for key in archive.files};meta=json.loads((folder/'SOURCE_AFFERENT_AUDIT.json').read_text())[0];p=np.array(meta['source_dimensions']);acts=z['actuators'].tolist();selected=z['selected_actuators'].tolist();t=z['neural_t_s'];keep=(t>=5)&(t<13)&z['source_fiber_ready'];kin=z['kinematics'];lf=kin[:,2]/p[:,1];vf=kin[:,3]/(p[:,1]*p[:,8]);rows={};categories=z['current_categories'].tolist()
    for family in ('tib_ant','soleus','med_gas','per_tert','vas_lat','bifemsh'):
        for side in ('r','l'):
            name=family+'_'+side;i=acts.index(name);j=selected.index(name)
            rows[name]=dict(MTU_length_range_m=[float(kin[keep,0,i].min()),float(kin[keep,0,i].max())],fiber_length_over_optimal_range=[float(lf[keep,i].min()),float(lf[keep,i].max())],MTU_velocity_absolute_p95_mps=float(np.percentile(abs(kin[keep,1,i]),95)),fiber_velocity_absolute_p95_mps=float(np.percentile(abs(kin[keep,3,i]),95)),fiber_velocity_over_source_Vmax_p95=float(np.percentile(abs(vf[keep,i]),95)),local_Ia_II_mean_current_nA=z['raw_local_inputs'][keep,:,j].mean(axis=0).tolist(),mean_MN_voltage_mV=float(z['voltage'][keep,j].mean()),mean_neural_command=float(z['command'][keep,j].mean()),mean_delivered_activation=float(z['delivered_activation'][keep,j].mean()),mean_current_nA={label:float(z['mn_currents'][keep,j,c].mean()) for c,label in enumerate(categories)},peak_abs_current_nA={label:float(np.max(abs(z['mn_currents'][keep,j,c]))) for c,label in enumerate(categories)})
    pf={name:dict(min_mV=float(z['PF_voltage'][keep,i].min()),max_mV=float(z['PF_voltage'][keep,i].max())) for i,name in enumerate(z['PF_cells'].tolist())}
    S=Path(__file__).resolve().parent.parent;target=np.load(S/'hunt_training/runs/targets/targets.npz');tm=json.loads((S/'hunt_training/runs/targets/metadata.json').read_text());feature={}
    for si,side in enumerate(('r','l')):
        indices=tm['feature_members'][side]
        feature[side]={}
        for k in range(4):
            members=indices[k];legacy=(kin[:,0]-np.array(tm['sensor_mid']))/np.array(tm['sensor_half']) if k<2 else kin[:,1]/.6
            proposed=(lf-1.)/.2 if k<2 else vf
            old=np.maximum(legacy[:,members],0.).mean(axis=1);new=np.maximum(proposed[:,members],0.).mean(axis=1)
            feature[side][tm['feature_names'][k]]=dict(original_expected_feature_range=[float(target['features'][:,si,k].min()),float(target['features'][:,si,k].max())],current_legacy_MTUbased_range=[float(old[keep].min()),float(old[keep].max())],candidate_source_fiber_range=[float(new[keep].min()),float(new[keep].max())],candidate_scope='Diagnostic dimensionless strain relative to original l0,scaled by0.2; speed relative to original Vmax*l0. Transduction scaling needs validation; not installed.')
    result=dict(scope='Observed saved-controller currents and actual committed fiber kinematics,logical5–13s. No encoder change or fit. B-feature comparisons are range diagnostics on different trajectories, not phase-accuracy validation.',maximum_MN_update_reconstruction_error_mV=meta['maximum_MN_update_reconstruction_error_mV'],motor_neurons=rows,phase_population_voltage=pf,Hunt_B_feature_comparison=feature,source_fiber_ready_samples=int(keep.sum()))
    (folder/'SOURCE_AFFERENT_DIAGNOSIS.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(error=meta['maximum_MN_update_reconstruction_error_mV'],ankle_PF=pf,tib_ant_r=rows['tib_ant_r']),indent=2));return result

def handle_request(request):
    """Bounded source-consistent ankle motor recruitment with timing frozen."""
    import subprocess
    from biomechanical_envelope import assess
    here=Path(__file__).resolve().parent;manifest=request['manifest'];study=Path(manifest['output_dir']);root=study.parent;spec=Path(manifest['ik_envelope_file']);roles=Path(manifest['functional_roles_file']);bases={1:'df0p3_air1',2:'df0p3_air2',3:'df0p3_stand',4:'df0p3_ground'}
    if hashlib.sha256(spec.read_bytes()).hexdigest()!=manifest['ik_envelope_sha256'] or hashlib.sha256(roles.read_bytes()).hexdigest()!=manifest['functional_roles_sha256']:raise RuntimeError('Frozen biomechanical specifications changed')
    if request['op']=='preflight':
        parity=json.loads((root/'CURRENT_AUDIT_PARITY.json').read_text());audit=json.loads((root/'current_audit_ground/SOURCE_AFFERENT_DIAGNOSIS.json').read_text())
        if not parity['all_bit_exact'] or audit['maximum_MN_update_reconstruction_error_mV']>1e-10:raise RuntimeError('Observational current audit not validated')
        for stage in (1,2,3):
            b=json.loads((root/bases[stage]/'BACKEND_RESULT.json').read_text())
            if b['guard'] is not None or not b['protocol_assessment']['accepted']:raise RuntimeError('Reduced-DF prerequisite failed')
        return dict(ok=True,ready_for_optimization=True,scope='Two common bilateral,route-verified ankle PF-to-MN conductance multipliers. Hunt A/B,timing,reflex,source dimensions,and IK envelope fixed. Every new candidate must pass air1/2 and standing before ground.')
    if request['op']!='evaluate':raise ValueError('Unknown motor calibration operation')
    p=request['params']
    if set(p)!={'ankle_df','ankle_pf'} or not .1<=float(p['ankle_df'])<=1. or not .75<=float(p['ankle_pf'])<=3.:raise ValueError('Declared functional ankle recruitment range only')
    rows=[];accepted=False;value=1000000.;identity=float(p['ankle_df'])==.3 and float(p['ankle_pf'])==1.
    for stage in (1,2,3,4):
        case=root/bases[stage] if identity else study/request['evaluation_id']/('stage'+str(stage))
        if not (case/'BACKEND_RESULT.json').exists():
            case.mkdir(parents=True,exist_ok=True);cfg=json.loads((root/('df0p3_stage'+str(stage)+'.json')).read_text());cfg['isolated_afferents']['motor_functional_scaling']={key:float(v) for key,v in p.items()};target=case/'INPUT.json';target.write_text(json.dumps(cfg,indent=2));protocol='hunt-air'+str(stage) if stage<3 else 'hunt-stand' if stage==3 else 'hunt-ground'
            with (case.parent/'simulation_output.log').open('a',encoding='utf-8') as log:subprocess.run([r'D:\Anaconda\envs\myo\python.exe',str(here/'live_tendon_trial.py'),'--protocol',protocol,'--physics','compliant','--source-task-basis','--config',str(target),'--out',str(case)],stdout=log,stderr=subprocess.STDOUT,timeout=180,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        b=json.loads((case/'BACKEND_RESULT.json').read_text());r=json.loads((case/'controller.result.json').read_text()) if (case/'controller.result.json').exists() else None;ok=b['guard'] is None and r is not None and r['assessment']['accepted'];row=dict(stage=stage,path=str(case),accepted=ok,guard=b['guard'],assessment=None if r is None else r['assessment']);rows.append(row)
        if b['guard'] is not None or r is None:
            trace=np.load(case/'backend_trace.npz');value=1000000.-float(trace['t'][-1]);break
        if stage<4 and not ok:value=5000.-float(r['assessment']['score']);break
        if stage==4:
            envelope=assess(case,spec);row['biomechanical_envelope']=envelope;swings=[foot['sustained_swings'] for foot in r['foot_clearance']['feet']];phase=float(r['ground_contacts']['phase_error']);accepted=ok and envelope['accepted_reference_envelope'];row.update(sustained_swings=swings,actual_assistance_bw=r['ground_contacts']['mean_harness_vertical_bw']);value=-100.-sum(swings)+phase if accepted else 200.+float(envelope['cost'])+50.*sum(max(0,3-n) for n in swings)+200.*phase+100.*(not ok)
    return dict(ok=True,value=float(value),diagnostics=dict(fully_accepted=bool(accepted),stages=rows,bilateral_parameters_tied=True,source_dimensions_fixed=True,Hunt_B_frozen=True))

def prepare_study(root):
    import shutil
    here=Path(__file__).resolve().parent;ib=root.parent/'hunt_ib_load_v1';folder=root/'motor_recruitment_v1';folder.mkdir(exist_ok=True);target=folder/'INPUT_MANIFEST.json'
    if target.exists():raise FileExistsError('Preserve original study')
    spec=ib/'IK_ENVELOPE_SPEC_V1.json';roles=ib/'FUNCTIONAL_ROLES_V1.json';manifest=dict(adapter_path=str(Path(__file__).resolve()),output_dir=str(folder),max_trials=6,wall_budget_s=1800,trial_timeout_s=540,seed=3101035,n_startup_trials=3,fixed_parameters={},search_space=dict(ankle_df=dict(low=.1,high=1.,log=True),ankle_pf=dict(low=.75,high=3.,log=True)),initial_candidates=[dict(ankle_df=d,ankle_pf=p) for d,p in [(.3,1.),(.3,1.5),(.6,1.5)]],ik_envelope_file=str(spec),ik_envelope_sha256=hashlib.sha256(spec.read_bytes()).hexdigest(),functional_roles_file=str(roles),functional_roles_sha256=hashlib.sha256(roles.read_bytes()).hexdigest(),acceptance='Same controller must pass actual air1/2,standing,ground engineering and frozen six-coordinate IK envelope. All physical knee states<=0. No source dimension or gate changes.',parameter_rationale='Complete current audit: TA PF current0.894nA versus Ia0.0129nA,delivered activation0.199; source-compliant ankles exceed reference through unloaded swing. Calibrate motor recruitment rather than alter B timing based on a weak local spindle contribution.')
    snapshot=folder/'source_snapshot';snapshot.mkdir()
    for name in ('source_afferents.py','stance_ib.py','live_tendon_trial.py','biomechanical_envelope.py','thelen_dynamic.py','balance_diagnostics.py'):shutil.copy2(here/name,snapshot/name)
    target.write_text(json.dumps(manifest,indent=2));print(target)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--summarize',type=Path);group.add_argument('--prepare-study',type=Path);args=parser.parse_args()
    if args.summarize:summarize(args.summarize)
    else:prepare_study(args.prepare_study)
