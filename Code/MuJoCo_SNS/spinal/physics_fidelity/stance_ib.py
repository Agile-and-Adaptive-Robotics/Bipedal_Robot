"""Opt-in source-normalized, phase/contact-permitted Ib load recruitment.

Existing excitatory Ib pools are reused. In force-pool modes their additive
RG_E drive is removed; RG/contact/descending permission instead modulates
their outputs. This is an explicit engineering presynaptic-gating adapter,
not a literal replication of a published interneuron circuit.
"""
import os,json,hashlib
from pathlib import Path
import numpy as np
from thelen_dynamic import parameters_from_osim

class Session:
    def __init__(self):self.networks=[]
    def save(self,out):
        records=[]
        for i,n in enumerate(self.networks):
            records.append(n.source_ib_metadata)
            np.savez_compressed(out/('source_ib_signals'+('' if i==0 else str(i))+'.npz'),samples=n.source_ib_log[:n.source_ib_count],columns=['neural_t_s','descending_gate','right_phase_gate','left_phase_gate','right_contact_gate','left_contact_gate','right_knee_group_force_F0','left_knee_group_force_F0','right_ankle_group_force_F0','left_ankle_group_force_F0'])
            np.savez_compressed(out/('source_ib_currents'+('' if i==0 else str(i))+'.npz'),neural_t_s=n.source_ib_log[:n.source_ib_count,0],current_nA=n.source_ib_current[:n.source_ib_count],gmax=n.source_ib_gmax[:n.source_ib_count],edges=n.source_ib_edge_names)
        (out/'SOURCE_IB_PATHWAY.json').write_text(json.dumps(records,indent=2));return records

def install(options,physical_source,stage,global_reflex_scale=1.):
    import controller as C
    import params as P
    from muscle_map import EXTENSOR_STANCE_GROUPS
    if options.get('format')!='aarl_stance_ib_v1':raise ValueError('Explicit source-Ib format required')
    mode=options['mode']
    if mode not in ('source_normalized','force_pool','stance_load'):raise ValueError('Unknown source-Ib mode')
    rolepath=Path(options['functional_roles_file']);roles=json.loads(rolepath.read_text())['roles']
    output=options.get('output_gmax',{})
    if mode=='stance_load' and (set(output)!={'knee_ext','ankle_pf'} or any(not np.isfinite(v) or not 0.<v<=4. for v in output.values())):raise ValueError('Declare both load conductances in (0,4]')
    if options.get('contact_permission_bw')!=[.025,.05]:raise ValueError('Frozen measured contact-permission thresholds required')
    session=Session();original=C.build
    def build(*args,**kwargs):
        n=original(*args,**kwargs);c=n.compiled;names=n.names;acts=list(n.muscles);n.source_ib_log=np.zeros((8500,10));n.source_ib_count=0
        ib_edges=np.array([i for i in n.reflex_edges if names[c.pre[i]].startswith(('IBEXC_','IBIN_'))],dtype=int)
        n.source_ib_edge_names=[names[c.pre[i]]+' -> '+names[c.post[i]] for i in ib_edges];n.source_ib_current=np.zeros((8500,len(ib_edges)));n.source_ib_gmax=np.zeros_like(n.source_ib_current)
        groups={};removed=[];excluded=[]
        def valid(group,act):
            mirror=act[:-2]+('_l' if act.endswith('_r') else '_r')
            return bool(roles[group].get(act,{}).get('valid_extension_route',False) and roles[group].get(mirror,{}).get('valid_extension_route',False))
        for group in ('knee_ext','ankle_pf'):
            edges=np.array([i for i in n.reflex_edges if names[c.pre[i]].startswith('IBEXC_'+group+'_')],dtype=int)
            weights=np.array([n.reflex_gain[i] if valid(group,names[c.post[i]][3:]) else 0. for i in edges]);maximum=float(weights.max()) if len(weights) else 0.
            if maximum<=0.:raise RuntimeError('No route-valid fitted load outputs: '+group)
            groups[group]=dict(edges=edges,original=weights.copy(),relative=weights/maximum,side=np.array([0 if names[c.post[i]].endswith('_r') else 1 for i in edges]),members={side:np.array([j for j,a in enumerate(acts) if a.endswith('_'+side) and valid(group,a) and n.muscles[a].groups[0]==group],dtype=int) for side in ('r','l')})
        if mode!='source_normalized':
            for i,(pre,post) in enumerate(zip(c.pre,c.post)):
                source,target=names[pre],names[post]
                for group in groups:
                    if not target.startswith('IBEXC_'+group+'_'):continue
                    if source=='RG_E_'+target[-1]:c.gmax[i]=n.original_g[i]=0.;removed.append(source+' -> '+target)
                    elif source.startswith('Ib_') and not valid(group,source[3:]):c.gmax[i]=n.original_g[i]=0.;excluded.append(source+' -> '+target)
            for group,g in groups.items():
                n.reflex_gain[g['edges']]=g['original']
                for i,w in zip(g['edges'],g['original']):
                    if w==0. and not valid(group,names[c.post[i]][3:]):excluded.append(names[c.pre[i]]+' -> '+names[c.post[i]])
            n.sync()
        original_bind=n.bind_mujoco
        def bind(m,d):
            original_bind(m,d);p=parameters_from_osim(physical_source,[m.actuator(i).name for i in range(m.nu)],0.);native=n.Fmax.copy();n.Fmax=p[:,0].copy();n.source_ib_F0=n.Fmax.copy()
            n.source_ib_metadata.update(native_vs_source_F0_ratio_range=[float((native/n.Fmax).min()),float((native/n.Fmax).max())],source_force_normalization='Original bilateral source max_isometric_force; native initializer retained; actual force field is source law after physical hook starts.')
        n.bind_mujoco=bind
        n.source_ib_metadata=dict(mode=mode,stage=stage,removed_phase_bias_edges=removed,excluded_wrong_function_routes=excluded,output_gmax=output,source_force_normalized=True,functional_roles_sha256=hashlib.sha256(rolepath.read_bytes()).hexdigest(),permission='clip((RG_E-RG_F)/E_HI,0,1) * clip((whole_foot_BW-.025)/.025,0,1) * clip(DRIVE/2.5,0,1)',encoding='Existing0.5*AFF.ib_gain,1-.7*RG_E/E_HI extensor attenuation,force/F0 clip0..1.5 retained; deafferented air stage never overwritten.',output_weighting='Route-valid, bilaterally common relative nonzero fitted Ib output pattern; zero fitted outputs remain zero.',runtime_current_scope='Existing positive Ib-pool→MN synapses; output conductance is gated, no target-voltage/torque cancellation.')
        n.source_ib_metadata.update(module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),actual_current_record='Pre-step voltage, actual fitted/gated gmax, compiled elo/erange/reversal; excitatory pools and autogenic inhibitory relays separately recorded.',global_reflex_scale=float(global_reflex_scale))
        phase=np.zeros(2);contact=np.zeros(2);force=np.zeros(len(acts));ports=np.array([n.ports['Ib_'+a] for a in acts]);isext=np.array([n.muscles[a].groups[0] in EXTENSOR_STANCE_GROUPS for a in acts]);sides=np.array([0 if n.muscles[a].side=='r' else 1 for a in acts]);rgE=np.array([n.idx['RG_E_'+s] for s in ('r','l')]);rgF=np.array([n.idx['RG_F_'+s] for s in ('r','l')]);native_step=n.step;deaff=set(os.environ.get('AARL_DEAFF','').split(','))-{''}
        def step(u):
            m,d=n.live;stance=np.clip(c.V[rgE]/P.E_HI,0.,1.);phase[:]=np.clip((c.V[rgE]-c.V[rgF])/P.E_HI,0.,1.);drive=float(np.clip(u[n.ports['DRIVE']]/2.5,0.,1.))
            np.negative(d.actuator_force,out=force);np.maximum(force,0.,out=force);np.divide(force,n.source_ib_F0,out=force);np.clip(force,0.,1.5,out=force)
            if stage!=1:
                gain=.5*P.AFF['ib_gain']*np.where(isext,1.-.7*stance[sides],1.)
                u[ports]=gain*force
                for j,a in enumerate(acts):
                    if n.muscles[a].side in deaff:u[ports[j]]=0.
            # The legacy LOAD port can have zero gain. Use the regional
            # sampler's once-counted force cache, populated immediately above
            # this wrapper by physical_eval; never infer force from a zero port.
            contact[:]=0.
            if hasattr(n,'regional_sensor'):
                load=n.regional_sensor.loads.sum(axis=1)/(float(m.body_mass.sum())*9.81)
                contact[:]=np.clip((load-.025)/.025,0.,1.)
            elif mode=='stance_load' and stage==4:raise RuntimeError('Actual regional force sampler required for ground load permission')
            if mode=='stance_load':
                permission=drive*phase*contact if stage!=1 else np.zeros(2)
                for group,g in groups.items():n.reflex_gain[g['edges']]=float(output[group])*g['relative']*permission[g['side']]*float(global_reflex_scale)
            k=n.source_ib_count
            if k<len(n.source_ib_log):
                means=[float(force[groups[group]['members'][side]].mean()) for group in ('knee_ext','ankle_pf') for side in ('r','l')]
                n.source_ib_log[k]=[k*P.DT,drive,*phase,*contact,*means]
                n.source_ib_gmax[k]=n.reflex_gain[ib_edges]
                n.source_ib_current[k]=n.reflex_gain[ib_edges]*np.clip((c.V[c.pre[ib_edges]]-c.elo[ib_edges])/c.erange[ib_edges],0.,1.)*(c.reversal[ib_edges]-c.V[c.post[ib_edges]])
                n.source_ib_count+=1
            return native_step(u)
        n.step=step;session.networks.append(n);return n
    C.build=build;return session

def handle_request(request):
    """Bounded Rustuna adapter; no physical dimensions or gates are fitted."""
    from biomechanical_envelope import assess
    import subprocess
    here=Path(__file__).resolve().parent;manifest=request['manifest'];study=Path(manifest['output_dir']);root=study.parent
    prerequisites=[]
    for stage,case in [(1,'stance_load_air1'),(2,'stance_load_air2'),(3,'stance_load_stand')]:
        folder=root/case;b=json.loads((folder/'BACKEND_RESULT.json').read_text());r=json.loads((folder/'controller.result.json').read_text());s=np.load(folder/'source_ib_signals.npz')['samples']
        permission=s[:,1,None]*s[:,2:4]*s[:,4:6]
        if b['guard'] is not None or not r['assessment']['accepted'] or np.any(permission!=0.):raise RuntimeError('Unchanged prerequisite must pass with identically inactive new load conductances: '+case)
        prerequisites.append(dict(stage=stage,path=str(folder),accepted=True,load_permission_maximum=float(permission.max()),reuse_proof='Both tuned conductances multiply exactly zero in every recorded prerequisite neural step. All other controller, sensory, and physical parameters are fixed. Winner still requires fresh explicit stage configs.'))
    spec=root/'IK_ENVELOPE_SPEC_V1.json';roles=root/'FUNCTIONAL_ROLES_V1.json'
    if hashlib.sha256(spec.read_bytes()).hexdigest()!=manifest['ik_envelope_sha256'] or hashlib.sha256(roles.read_bytes()).hexdigest()!=manifest['functional_roles_sha256']:raise RuntimeError('Frozen reference/role specification changed')
    if request['op']=='preflight':return dict(ok=True,ready_for_optimization=True,prerequisites=prerequisites,scope='Two shared stance-specific Ib output conductances only; source dimensions, FF/POST/PT, timing, physical support, and strict gates fixed.')
    if request['op']!='evaluate':raise ValueError('Unsupported Ib-study operation')
    p=request['params']
    if set(p)!={'g_knee','g_ankle'} or any(not .05<=float(v)<=4. for v in p.values()):raise ValueError('Predeclared common bilateral Ib conductances only')
    case=study/request['evaluation_id']/'stage4';case.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((root/'stance_load_stage4.json').read_text());cfg['isolated_ib_load']['output_gmax']={'knee_ext':float(p['g_knee']),'ankle_pf':float(p['g_ankle'])};target=case/'INPUT.json'
    if not (case/'BACKEND_RESULT.json').exists():
        target.write_text(json.dumps(cfg,indent=2))
        with (case.parent/'simulation_output.log').open('a',encoding='utf-8') as log:
            subprocess.run([r'D:\Anaconda\envs\myo\python.exe',str(here/'live_tendon_trial.py'),'--protocol','hunt-ground','--physics','compliant','--source-task-basis','--config',str(target),'--out',str(case)],stdout=log,stderr=subprocess.STDOUT,timeout=180,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    elif json.loads(target.read_text())!=cfg:raise RuntimeError('Existing trial has different parameters')
    b=json.loads((case/'BACKEND_RESULT.json').read_text());r=json.loads((case/'controller.result.json').read_text()) if (case/'controller.result.json').exists() else None
    envelope=assess(case,spec) if (case/'controller.npz').exists() else None;accepted=False;value=1000000.
    if b['guard'] is None and r is not None and b['complete']:
        engineering=bool(r['assessment']['accepted']);biomechanical=bool(envelope['accepted_reference_envelope']);swings=[f['sustained_swings'] for f in r['foot_clearance']['feet']];phase=float(r['ground_contacts']['phase_error']);accepted=engineering and biomechanical
        value=(-100.-sum(swings)+phase if accepted else 200.+float(envelope['cost'])+50.*sum(max(0,3-n) for n in swings)+200.*phase+100.*(not engineering))
    else:
        z=np.load(case/'backend_trace.npz');value-=min(float(z['t'][-1]),17.)
    return dict(ok=True,value=float(value),diagnostics=dict(accepted=accepted,prerequisites=prerequisites,guard=b['guard'],engineering_assessment=None if r is None else r['assessment'],biomechanical_envelope=envelope,source_dimensions_fixed=True,bilateral_parameters_tied=True,case=str(case)))

def prepare_study(root):
    """Freeze a bounded two-conductance search before evaluating any points."""
    import shutil
    here=Path(__file__).resolve().parent;study=root/'load_recruitment_v1';study.mkdir(exist_ok=True);target=study/'INPUT_MANIFEST.json'
    if target.exists():raise FileExistsError('Preserve existing study')
    cfg=dict(adapter_path=str(Path(__file__).resolve()),output_dir=str(study),max_trials=8,wall_budget_s=1800,trial_timeout_s=240,seed=3101034,n_startup_trials=6,fixed_parameters={},search_space={key:dict(low=.05,high=4.,log=True) for key in ('g_knee','g_ankle')},initial_candidates=[dict(g_knee=k,g_ankle=a) for k,a in [(.05,.05),(.3,.3),(1.,1.),(3.,3.),(.3,1.),(1.,.3)]],ik_envelope_sha256=hashlib.sha256((root/'IK_ENVELOPE_SPEC_V1.json').read_bytes()).hexdigest(),functional_roles_sha256=hashlib.sha256((root/'FUNCTIONAL_ROLES_V1.json').read_bytes()).hexdigest(),acceptance='Original contact/clearance/frontal/phase/safety gates AND frozen per-side saved IK envelope. No ROM margin, no post-hoc gate change. These do not establish normal gait.',objective='Accepted: -100-total sustained swings+contact phase error. Otherwise:200+envelope cost+50*missing swings+200*phase error+100*engineering failure. Guard/incomplete:1e6-last physical time.',fixed_source_dimensions='Original bilateral pooled Thelen2003 F0, optimal fiber length, slack tendon length, pennation, source FV and tendon strain; beta0. Existing compliant backend, 2ms neural/physics step.')
    snapshot=study/'source_snapshot';snapshot.mkdir()
    for name in ('stance_ib.py','biomechanical_envelope.py','live_tendon_trial.py','balance_diagnostics.py','thelen_dynamic.py'):shutil.copy2(here/name,snapshot/name)
    target.write_text(json.dumps(cfg,indent=2));print(str(target))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--prepare-study',type=Path);args=parser.parse_args()
    if args.prepare_study:prepare_study(args.prepare_study)
    else:parser.error('Use --prepare-study, or import the isolated adapter')
