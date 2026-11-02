"""Read-only controlled-II audit of saved Hunt-D trial 0; no body steps."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import mujoco
sys.path.insert(0,str(Path(__file__).resolve().parent))
import integrated_hunt_design as I


def differences(a,b,prefix=''):
    if isinstance(a,dict) and isinstance(b,dict):
        return sum((differences(a.get(k),b.get(k),prefix+'/'+k) for k in sorted(set(a)|set(b))),[])
    return [] if a==b else [{'path':prefix,'control':a,'trial':b}]


def load_network(cfg):
    acts=[a.attrib['name'] for a in ET.parse(I.MODEL).getroot().find('actuator')]
    return I.build(acts,options=I.condition_options(cfg),profile=cfg['integrated_fsa_profile'],protocol_stage=2)


def currents(net,path):
    with I.unbound_electrical(net):net.apply_weights(1.)
    z=np.load(path,allow_pickle=True);names=z['names'].tolist();order=z['acts'].tolist();c=net.compiled;v=z['voltage']
    if names!=net.names[:len(c.V)]:raise ValueError('Recorded neuronal order mismatch')
    groups={k:np.zeros((len(v),len(order))) for k in (*I.GP.KINDS,'Ib_group')}
    for edge in net.mn_edges:
        pre=net.names[c.pre[edge]];post=net.names[c.post[edge]];kind=I.GP.edge_kind(pre,post,c.reversal[edge])
        if pre.startswith('IBEXC_'):kind='Ib_group'
        if kind is None:continue
        groups[kind][:,order.index(post[3:])]+=c.gmax[edge]*np.clip((v[:,c.pre[edge]]-c.elo[edge])/c.erange[edge],0.,1.)*(c.reversal[edge]-v[:,c.post[edge]])
    errors={'Ia_Ib_recorded_group_max_error_nA':float(np.max(abs(groups['Ia_homonymous']+groups['Ib_autogenic']+groups['Ib_group']-z['currents'][:,1,:]))),'reciprocal_recorded_group_max_error_nA':float(np.max(abs(groups['Ia_reciprocal_pool']-z['currents'][:,2,:])))}
    if max(errors.values())>1e-10:raise AssertionError('Recorded delivered current mismatch')
    selected={}
    for side in ('r','l'):
        f=v[:,names.index('RG_F_'+side)]>v[:,names.index('RG_E_'+side)]
        for label,phase in (('F_dominant',f),('E_dominant',~f)):
            keep=(z['t']>=5.)&(z['t']<14.)&phase;rows={}
            for muscle in ('psoas','iliacus','bifemlh','soleus','tib_ant'):
                act=muscle+'_'+side;q=order.index(act)
                rows[act]={'samples':int(keep.sum()),'II_current_mean_nA':float(groups['II_homonymous'][keep,q].mean()),'II_current_rms_nA':float(np.sqrt(np.mean(groups['II_homonymous'][keep,q]**2))),'all_local_mean_nA':{k:float(val[keep,q].mean()) for k,val in groups.items()},'recorded_total_synaptic_current_mean_nA':float(z['currents'][keep,:,q].sum(axis=1).mean()),'MN_voltage_mean_mV':float(v[keep,names.index('MN_'+act)].mean()),'activation_mean':float(z['act'][keep,q].mean())}
            selected[side+'_'+label]=rows
    return {'source':str(path),'sha256':I.sha(path),'aggregate_replay_checks':errors,'RG_phase_summaries':selected},z,groups


def run(root,out):
    row=json.loads((root/'training/hunt_D_local_v1/TRIALS.jsonl').read_text().splitlines()[0]);trial_cfg_path=Path(row['diagnostics']['configuration']);control_cfg_path=root/'configurations/fsa_v2_equivalent_s2.json'
    cfgs=[json.loads(p.read_text()) for p in (control_cfg_path,trial_cfg_path)];config_diff=differences(*cfgs)
    if {q['path'] for q in config_diff}!={'/integrated_fsa_profile','/verification_capture'}:raise AssertionError('Other body/runtime configuration differs')
    profiles=[I.load_profile(c['integrated_fsa_profile']) for c in cfgs];profile_diff=differences(*profiles)
    if {q['path'] for q in profile_diff}!={'/parameter_checkpoint','/local_profile/contracts/II_homonymous'}:raise AssertionError('Not an isolated II target change')
    nets=[load_network(c) for c in cfgs];a,b=nets;ca,cb=a.compiled,b.compiled
    if a.names!=b.names or a.inputs!=b.inputs:raise AssertionError('Neuronal/input order differs')
    exact_fields=[];changed={}
    fields=I.GP.FIELDS+I.GP.ION_FIELDS+('V','V_last','b_gate','b_gate_last','c_gate','c_gate_last')
    for field in fields:
        av,bv=np.asarray(getattr(ca.c,field)),np.asarray(getattr(cb.c,field))
        if np.array_equal(av,bv):exact_fields.append(field)
        else:changed[field]=int(np.count_nonzero(av!=bv))
    if changed!={'g_max_non':92}:raise AssertionError(('Unexpected parameter/initial-state changes',changed))
    for name in ('pf_gain','ff_gain','bi_gain'):
        if not np.array_equal(getattr(a,name),getattr(b,name)):raise AssertionError('Other authoritative gain buffer changed')
    if not np.array_equal(ca.post,cb.post) or not np.array_equal(ca.pre,cb.pre):raise AssertionError('Sparse topology changed')
    ids=np.flatnonzero(ca.gmax!=cb.gmax)
    if len(ids)!=92 or any(I.GP.edge_kind(a.names[ca.pre[i]],a.names[ca.post[i]],ca.reversal[i])!='II_homonymous' for i in ids):raise AssertionError('Sparse changed edges not exactly homonymous II92')
    ratios=cb.gmax[ids]/ca.gmax[ids]
    command_checks=[]
    for command in (0.,.5,1.,0.):
        for n in nets:
            with I.unbound_electrical(n):n.apply_weights(command)
            I.verify(n,full=True)
        if not np.array_equal(np.flatnonzero(ca.gmax!=cb.gmax),ids):raise AssertionError('Command changed isolated scope')
        command_checks.append(command)
    physics={}
    caps=[Path(c['verification_capture']) for c in cfgs]
    modelpaths=[p/'model.mjb' for p in caps]
    models=[mujoco.MjModel.from_binary_path(str(p)) for p in modelpaths];model_diff=[];model_fields=0
    for name in dir(models[0]):
        val=getattr(models[0],name)
        if isinstance(val,np.ndarray):
            model_fields+=1
            if not np.array_equal(val,getattr(models[1],name)):model_diff.append(name)
    if model_diff:raise AssertionError(('Compiled physical parameter mismatch',model_diff))
    if models[0].opt.timestep!=models[1].opt.timestep:raise AssertionError('Physics timestep changed')
    if I.sha(modelpaths[0])!=I.sha(modelpaths[1]):raise AssertionError('Complete compiled model binaries differ')
    initial={}
    zs=[np.load(p/'aligned_trace.npz',allow_pickle=True) for p in caps]
    for key in ('qpos','qvel','activation','control','neural_voltage'):
        if key in zs[0] and key in zs[1]:initial[key]={'first_record_exact':bool(np.array_equal(zs[0][key][0],zs[1][key][0])),'max_abs_difference':float(np.max(abs(zs[0][key][0]-zs[1][key][0])))}
    saved_results=[json.loads(p.with_suffix('.result.json').read_text()) for p in (control_cfg_path,trial_cfg_path)]
    runtime_source_differences=differences(*(r['source_sha256'] for r in saved_results))
    if runtime_source_differences:raise AssertionError(('Saved runtime source hashes differ',runtime_source_differences))
    physics={'model_files':[str(p) for p in modelpaths],'model_sha256':[I.sha(p) for p in modelpaths],'complete_compiled_model_binary_exact':True,'saved_runtime_source_hashes_exact':True,'all_compiled_numeric_model_arrays_exact':True,'arrays_checked':model_fields,'timestep_s':models[0].opt.timestep,'first_aligned_saved_state':initial,'capture_keys':zs[0].files,'first_state_time_s':[float(z['t_state'][0]) for z in zs]}
    paths=[control_cfg_path.with_suffix('.mechanics.npz'),trial_cfg_path.with_suffix('.mechanics.npz')]
    replay=[currents(n,p) for n,p in zip(nets,paths)];counterfactual={}
    saved_initial={key:bool(np.array_equal(replay[0][1][key][0],replay[1][1][key][0]))for key in ('voltage','qpos','qvel','act','ctrl','force')}
    if not all(saved_initial.values()):raise AssertionError(('Pre-warmup initial saved state differs',saved_initial))
    for label,(record,z,groups) in zip(('control_trajectory','trial_trajectory'),replay):
        order=z['acts'].tolist();q=order.index('psoas_r');v=z['voltage'];pre=ca.pre[ids];post=ca.post[ids]
        errors=[]
        for i in ids:
            alpha=np.clip((v[:,ca.pre[i]]-ca.elo[i])/ca.erange[i],0,1);drive=ca.reversal[i]-v[:,ca.post[i]]
            errors.append(np.max(abs(cb.gmax[i]*alpha*drive-ratios[0]*ca.gmax[i]*alpha*drive)))
        counterfactual[label]={'same_state_only_II_current_changes':True,'II_current_ratio':float(ratios[0]),'maximum_ratio_identity_error_nA':float(max(errors)),'scope':'Identical recorded presynaptic and MN voltages; actual trajectory-dependent mean currents need not share this ratio.'}
    result={'schema':'controlled_homonymous_II_refinement_v1','machine':'easteregg2','body_steps_run_by_audit':0,'training_run_by_audit':False,'trial_number':row['number'],'configuration_differences':config_diff,'profile_differences':profile_diff,'all_other_compiled_parameter_and_initial_neural_arrays_exact':True,'exact_compiled_fields':exact_fields,'changed_compiled_fields':changed,'changed_sparse_edges':len(ids),'changed_edge_kind':'II_homonymous','gain_ratio_min_max':list(map(float,(ratios.min(),ratios.max()))),'other_authoritative_gain_buffers_exact':['pf_gain','ff_gain','bi_gain'],'command_checks':command_checks,'physics':physics,'delivered_current_replay':{'window_s':[5.,14.],'phase_scope':'RG dominance in air, not physical stance/swing','conditions':dict(zip(('control','trial0'),(r[0] for r in replay)))},'same_state_counterfactual':counterfactual,'saved_assessments':{'control':json.loads(control_cfg_path.with_suffix('.result.json').read_text())['assessment'],'trial0':row['diagnostics']['assessment']},'source_hashes':{'audit':I.sha(__file__),'control_configuration':I.sha(control_cfg_path),'trial_configuration':I.sha(trial_cfg_path)},'interpretation':'The controlled input change is the92 homonymous-II conductances derived from +0.5 to +0.02mV targets; removing the old parameter checkpoint has separately proven neutral effect at unchanged targets. Saved stage2 acceptance changes from fail to pass. Current/activation changes measured on different trajectories include the closed-loop kinematic response; same-state counterfactuals isolate the immediate parameter effect. This is one supported air-stage test, not ground stepping or basin/fresh robustness.'}
    result['recorded_initial_state_before_warmup']={'time_s':float(replay[0][1]['t'][0]),'fields_exact':saved_initial}
    result['held_preparation_actual_activation_checks']={label:{'last_preparation_mechanics_sample_s':float(z['t'][z['t']<.6][-1]),'maximum_actual_activation_change_from_t0':float(np.max(abs(z['act'][z['t']<.6]-z['act'][0])))}for label,(_,z,_)in zip(('control','trial0'),replay)}
    result['initialization_scope']='Recorded neural voltage, qpos, qvel, actual activation, control, and force agree exactly before preparation at logical 0 s. During logical t < 0.6 s, neural updates and changing ctrl continue, but runner calls mj_forward rather than mj_step. Do not infer activation integration from the stale warmup comment or legacy NPZ act (which logs ctrl). At the first actual physics interval ending at logical 0.602 s, activation, control, and neural voltage differ, while qpos and qvel remain exact; activation can first diverge as ctrl-driven dynamics integrate. These are consequences of the controlled II parameter change, not independently changed initialization. This test does not separate ongoing feedback from the neural preparation-state effect.'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({'output':str(out),'only_II92_changed':True,'control_accepted':result['saved_assessments']['control']['accepted'],'trial0_accepted':result['saved_assessments']['trial0']['accepted']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args();run(args.root,args.out)
