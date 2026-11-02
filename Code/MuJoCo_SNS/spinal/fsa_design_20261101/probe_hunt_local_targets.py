"""Static proof of shared Hunt-D functional targets and saved-current replay."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,copy,hashlib,json,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import integrated_hunt_design as I


def digest(array):return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def run(profile,out):
    p=I.load_profile(profile);cfg=json.loads(Path(p['condition_file']).read_text());acts=[q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')]
    base=I.build(acts,options=I.condition_options(cfg),profile=p);base_hash=I.GP.field_digest(base.compiled.c);arrays={};results={}
    controls={'baseline_without_checkpoint':{k:(.5 if k in I.GP.KINDS[:2] else -.5) for k in I.GP.KINDS},'II_0p02':{k:(.02 if k=='II_homonymous' else v) for k,v in p['local_profile']['contracts'].items()},'weak_all':dict(zip(I.GP.KINDS,(.1,.02,-.1,-.1)))}
    for label,contracts in controls.items():
        candidate=copy.deepcopy(p);candidate['parameter_checkpoint']=None;candidate['local_profile']['contracts']=contracts
        net=I.build(acts,options=I.condition_options(cfg),profile=candidate);c=net.compiled;design=net.hunt_graded_fsa['design_ledger'];rows=[]
        for pool in design['pools']:
            q=pool['post_index'];sources=np.array(pool['presynaptic_indices']);voltage=pool['source_voltage_mV'];effective=c.gmax[[np.flatnonzero((c.post==q)&(c.pre==s))[0] for s in sources]]*np.array([I.GP.activation_fraction(c.c,q,int(s),v) for s,v in zip(sources,voltage)])
            steady=(pool['bias_current_for_probe_nA']+sum(effective*c.c.del_e[q,sources]))/(c.c.g_m[q]+sum(effective));error=abs(steady-(candidate['local_profile']['operating_post_mV']+contracts[pool['kind']]))
            if error>1e-12:raise AssertionError('New target not consumed at actual final pool')
            rows.append({'kind':pool['kind'],'post':pool['post'],'error_mV':float(error)})
        for command in (0.,.5,1.,0.):
            with I.unbound_electrical(net):net.apply_weights(command)
            I.verify(net,full=True)
        with I.unbound_electrical(net):net.apply_weights(1.)
        exact=all(np.array_equal(getattr(base,field),getattr(net,field)) for field in ('pf_gain','ff_gain','bi_gain'))
        if not exact:raise AssertionError('Non-reflex learned gains changed')
        for key in ('pre','post','elo','erange','reversal'):
            if not np.array_equal(getattr(base.compiled,key),getattr(c,key)):raise AssertionError('Topology or windows changed with local targets')
        timing_exact=all(a['new_conductance_uS']==b['new_conductance_uS'] for a,b in zip(base.hunt_timing_normalized_design['edges'],net.hunt_timing_normalized_design['edges']))
        if not timing_exact:raise AssertionError('Timing was not frozen')
        baseline_exact=I.GP.field_digest(c.c)==base_hash
        if label=='baseline_without_checkpoint' and (not baseline_exact or not np.array_equal(net.reflex_gain,base.reflex_gain)):raise AssertionError('Clearing checkpoint changed baseline design')
        ids=net.hunt_graded_fsa['protected_indices'];arrays[label]=c.gmax[ids].copy()
        if 'source_names' not in arrays:
            arrays['source_names']=np.array([net.names[c.pre[i]] for i in ids]);arrays['target_names']=np.array([net.names[c.post[i]] for i in ids]);arrays['kinds']=np.array([I.GP.edge_kind(net.names[c.pre[i]],net.names[c.post[i]],c.reversal[i]) for i in ids])
        summary={kind:{'edges':int(sum(arrays['kinds']==kind)),'minimum_g_uS':float(arrays[label][arrays['kinds']==kind].min()),'maximum_g_uS':float(arrays[label][arrays['kinds']==kind].max())} for kind in I.GP.KINDS}
        results[label]={'contracts_mV':contracts,'parameter_checkpoint':None,'static_pools_checked':len(rows),'maximum_error_mV':max(r['error_mV'] for r in rows),'other_gain_buffers_exact':exact,'timing_frozen':timing_exact,'baseline_parameter_graph_exact':baseline_exact,'actual_authoritative_gain_summary':summary,'command_persistence_verified':[0.,.5,1.,0.],'neural_steps':0,'body_steps':0}
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);gainfile=out.with_name('LOCAL_TARGET_GAINS.npz');np.savez_compressed(gainfile,**arrays)
    record={'local_targets_consumed':True,'other_gain_buffers_exact':True,'timing_frozen':True,'baseline_without_checkpoint_exact':True,'cases':results,'exact_named_gain_ledger':str(gainfile),'gain_ledger_sha256':I.sha(gainfile),'profile_source':str(profile),'profile_sha256':I.sha(profile),'test_source_sha256':I.sha(__file__),'suggested_initial_engineering_bounds_mV':{'Ia_homonymous':[.05,.5],'II_homonymous':[.005,.15],'Ib_autogenic':[-.5,-.05],'Ia_reciprocal_pool':[-.5,-.05]},'bounds_scope':'Small source-centered engineering search ranges,not physiological identified bounds; zero excluded by signed-target definition','training_launched':False,'body_steps':0}
    out.write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps({'local_targets_consumed':True,'other_gain_buffers_exact':True,'timing_frozen':True,'baseline_without_checkpoint_exact':True,'output':str(out)}))


def replay_saved_currents(profile,root,out):
    """Recompute fixed local gains against actual recorded pre-update voltages.

    Geometry-dependent FF gains are not inferred. Two recorded aggregate groups
    independently verify the fixed-gain replay before the finer II split is used.
    """
    p=I.load_profile(profile);cfg=json.loads(Path(p['condition_file']).read_text());acts=[q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')];results=[];root=Path(root)
    for label in ('source_baseline_s2','fsa_equivalent_s2'):
        net=I.C.build(acts,options=I.condition_options(cfg)) if label.startswith('source') else I.build(acts,options=I.condition_options(cfg),profile=p)
        with I.unbound_electrical(net):net.apply_weights(1.)
        path=root/'configurations'/(label+'.mechanics.npz');z=np.load(path,allow_pickle=True);names=z['names'].tolist();order=z['acts'].tolist()
        if names!=net.names[:len(net.compiled.V)]:raise ValueError('Recorded true-node order does not match actual source constructor')
        c=net.compiled;v=z['voltage'];groups={kind:np.zeros((len(v),len(order))) for kind in (*I.GP.KINDS,'Ib_group')}
        for edge in net.mn_edges:
            pre=net.names[c.pre[edge]];post=net.names[c.post[edge]];kind=I.GP.edge_kind(pre,post,c.reversal[edge])
            if pre.startswith('IBEXC_'):kind='Ib_group'
            if kind is None:continue
            value=c.gmax[edge]*np.clip((v[:,c.pre[edge]]-c.elo[edge])/c.erange[edge],0.,1.)*(c.reversal[edge]-v[:,c.post[edge]])
            groups[kind][:,order.index(post[3:])]+=value
        check1=groups['Ia_homonymous']+groups['Ib_autogenic']+groups['Ib_group'];check2=groups['Ia_reciprocal_pool'];errors={'recorded_Ia_Ib_group_max_error_nA':float(np.max(abs(check1-z['currents'][:,1,:]))),'recorded_IaIN_group_max_error_nA':float(np.max(abs(check2-z['currents'][:,2,:])))}
        if max(errors.values())>1e-10:raise AssertionError('Fixed-gain replay does not match recorded delivered currents')
        summary={}
        for side in ('r','l'):
            dominant=v[:,names.index('RG_F_'+side)]>v[:,names.index('RG_E_'+side)];summary[side]={}
            for phase,condition in (('F_dominant',dominant),('E_dominant',~dominant)):
                keep=(z['t']>=5.)&(z['t']<14.)&condition;phase_rows={}
                for muscle in ('iliacus','psoas','bifemlh','semimem','bifemsh','rect_fem','vas_lat','vas_int','soleus','tib_ant'):
                    act=muscle+'_'+side;index=order.index(act);phase_rows[act]={'local_current_mean_nA':{kind:float(value[keep,index].mean()) for kind,value in groups.items()},'recorded_total_synaptic_current_mean_nA':float(z['currents'][keep,:,index].sum(axis=1).mean()),'recorded_MN_mean_mV':float(v[keep,names.index('MN_'+act)].mean()),'activation_mean':float(z['act'][keep,index].mean())}
                summary[side][phase]=phase_rows
        results.append({'condition':label,'source':str(path),'source_sha256':I.sha(path),'recorded_aggregate_replay_checks':errors,'selected_motor_neurons_by_RG_phase':summary,'scope':'Actual pre-update recorded voltages plus fixed actual local gains/windows; aggregate checks independent. II is split from recorded other group; geometry-dependent FF current is not reconstructed.'})
    record={'window_s':[5.,14.],'phase_scope':'RG dominance in suspended-air trials,not measured physical stance/swing','rows':results,'new_body_steps':0,'training':False,'profile_source':str(profile),'profile_sha256':I.sha(profile),'test_source_sha256':I.sha(__file__)};Path(out).write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps({'actual_recorded_current_replay_verified':True,'output':str(out)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--profile',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--trace-root',type=Path);a=parser.parse_args();replay_saved_currents(a.profile,a.trace_root,a.out) if a.trace_root else run(a.profile,a.out)
