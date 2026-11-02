"""Electrical/source-scope proof of actual default-off caller posture policies."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,ast,copy,json,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,mujoco
sys.path.insert(0,str(Path(__file__).resolve().parent))
import integrated_hunt_design as I


def policy_node(tree,key):
    nodes=[n for n in ast.walk(tree)if isinstance(n,ast.If)and any(isinstance(x,ast.Constant)and x.value==key for x in ast.walk(n.test))]
    # The mutual-exclusion guard contains both keys; choose the single-key block.
    return next(n for n in nodes if not isinstance(n.test,ast.BoolOp))


def execute(node,cfg,n,u=None):
    env={'cfg':cfg,'n':n,'np':np,'capture':{}}
    if u is not None:env['u']=u
    exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(node)],type_ignores=[])),'actual_physical_eval_policy','exec'),env)
    return env['capture']


def run(root,out):
    runtime=I.H/'physical_eval.py';tree=ast.parse(runtime.read_text());full=policy_node(tree,'suppress_inherited_posture');mn_only=policy_node(tree,'suppress_inherited_posture_mn');guard=next(n for n in tree.body if isinstance(n,ast.If)and isinstance(n.test,ast.BoolOp)and 'suppress_inherited_posture_mn' in ast.unparse(n.test))
    archive=root/'integration/input_policy_v1/pre_mask_physical_eval.py';old=ast.parse(archive.read_text());adv=lambda t:next(n for n in ast.walk(t)if isinstance(n,ast.FunctionDef)and n.name=='advance');newadv=copy.deepcopy(adv(tree));newadv.body=[n for n in newadv.body if ast.dump(n)!=ast.dump(full)]
    if ast.dump(newadv)!=ast.dump(adv(old)):raise AssertionError('Another numerical advance change occurred')
    rejected=False
    try:execute(guard,{'suppress_inherited_posture':True,'suppress_inherited_posture_mn':True},None)
    except ValueError:rejected=True
    if not rejected:raise AssertionError('Conflicting policies not rejected')
    cfg=json.loads((root/'configurations/candidate_ii02_s3.json').read_text());p=I.load_profile(cfg['integrated_fsa_profile']);acts=[e.attrib['name']for e in ET.parse(I.MODEL).getroot().find('actuator')];nets={};metadata={}
    policies={'default':{},'false':{'suppress_inherited_posture':False,'suppress_inherited_posture_mn':False},'MN_only':{'suppress_inherited_posture_mn':True},'full_source_diagnostic':{'suppress_inherited_posture':True}}
    for name,policy in policies.items():
        n=I.C.build(acts,options=I.condition_options(cfg));metadata[name]=execute(mn_only,policy,n);I.apply_to_network(n,p,protocol_stage=3);nets[name]=n
    base=nets['default'];b=base.compiled;posture=base.idx['POSTURE'];ids=np.flatnonzero((b.pre==posture)&np.array([base.names[k].startswith('MN_')for k in b.post]));other=np.setdiff1d(np.arange(len(b.gmax)),ids);nonmn=np.flatnonzero((b.pre==posture)&~np.array([base.names[k].startswith('MN_')for k in b.post]));paired=nets['MN_only'];c=paired.compiled
    if len(ids)!=88:raise AssertionError(('Unexpected global MN posture count',len(ids)))
    if not np.array_equal(np.flatnonzero(base.original_g!=paired.original_g),ids):raise AssertionError('original_g changes exceed88MN slots')
    for field in ('pf_gain','ff_gain','reflex_gain','bi_gain'):
        if not np.array_equal(getattr(base,field),getattr(paired,field)):raise AssertionError('Learned gain buffer changed')
    for field in I.GP.FIELDS+I.GP.ION_FIELDS:
        if field!='g_max_non' and not np.array_equal(getattr(b.c,field),getattr(c.c,field)):raise AssertionError('Other neural parameters changed')
    commands=[]
    for command in (0.,.5,1.,0.):
        for n in nets.values():
            with I.unbound_electrical(n):n.apply_weights(command)
            I.verify(n,full=True)
        if np.any(c.gmax[ids])or np.any(c.c.g_max_non[c.post[ids],c.pre[ids]])or np.any(paired.original_g[ids]):raise AssertionError('MN-only zeros reactivated or caches stale')
        if not np.array_equal(b.gmax[other],c.gmax[other]):raise AssertionError('Non-MN-posture weights changed')
        if not np.array_equal(b.gmax,nets['false'].compiled.gmax):raise AssertionError('Explicitfalse differs from default')
        commands.append(command)
    direct=json.loads((root/'integration/v2/STANDING_PATHWAY_AUDIT.json').read_text())['posture_inputs']['all_direct_POST_currents_nA'];u=np.zeros(len(base.inputs));u[base.ports['POSTURE']]=1.
    for act,value in zip(acts,direct):u[base.ports['POST_'+act]]=value
    postports=np.array([base.ports['POST_'+a]for a in acts]);maxfalse=0.;MNcurrent=0.;fullcurrent=0.;sourcevolt={key:[]for key in nets};rgerror=0.
    for step in range(128):
        results={}
        for name,n in nets.items():
            inp=u.copy();execute(full,policies[name],n,inp)
            if not np.array_equal(inp[postports],u[postports]):raise AssertionError('DirectPOST inputs changed')
            changed=np.flatnonzero(inp!=u)
            if name=='full_source_diagnostic':
                if not np.array_equal(changed,np.array([n.ports['POSTURE']])):raise AssertionError('Full-source mask changes another input')
            elif len(changed):raise AssertionError('MN-only/default policy mutates input')
            with I.unbound_electrical(n):n.step(inp)
            if step%32==0:I.verify(n,full=True)
            results[name]=n.compiled.c.V.copy();cc=n.compiled
            active=cc.gmax[ids]*np.clip((cc.c.V[cc.pre[ids]]-cc.elo[ids])/cc.erange[ids],0,1)*(cc.reversal[ids]-cc.c.V[cc.post[ids]])
            if name=='MN_only':MNcurrent=max(MNcurrent,float(np.max(abs(active))))
            if name=='full_source_diagnostic':fullcurrent=max(fullcurrent,float(np.max(abs(active))))
            sourcevolt[name].append(float(cc.c.V[n.idx['POSTURE']]))
        maxfalse=max(maxfalse,float(np.max(abs(results['default']-results['false']))))
        rgids=[base.idx['RG_E_'+s]for s in ('r','l')];rgerror=max(rgerror,float(np.max(abs(results['default'][rgids]-results['MN_only'][rgids]))))
    if maxfalse!=0 or MNcurrent!=0 or fullcurrent!=0 or max(sourcevolt['full_source_diagnostic'])!=0:raise AssertionError('Default parity or zero-source/current proof failed')
    model=root/'figures/candidate_ii02_s3/model.mjb';before=I.sha(model);m=mujoco.MjModel.from_binary_path(str(model));arrays={k:getattr(m,k).copy()for k in dir(m)if isinstance(getattr(m,k),np.ndarray)}
    execute(mn_only,policies['MN_only'],paired);execute(full,policies['full_source_diagnostic'],nets['full_source_diagnostic'],u.copy())
    if I.sha(model)!=before or any(not np.array_equal(v,getattr(m,k))for k,v in arrays.items()):raise AssertionError('Physics changed')
    record={'schema':'caller_posture_policies_v1','physical_eval_source':str(runtime),'physical_eval_sha256':I.sha(runtime),'audit_source_sha256':I.sha(__file__),'pre_mask_archive':str(archive),'pre_mask_sha256':I.sha(archive),'advance_excluding_input_mask_AST_exact':True,'conflicting_policies_rejected':rejected,'profiles_rebound':False,'default_and_explicit_false_128steps_max_error_mV':maxfalse,'MN_only':{'masked_edges':len(ids),'original_g_owns_zeros':True,'command_checks':commands,'actual_MN_current_max_nA':MNcurrent,'POSTURE_source_final_mV':sourcevolt['MN_only'][-1],'POSTURE_nonMN_edges_preserved':[{'source':'POSTURE','target':base.names[b.post[i]],'g_uS':float(b.gmax[i])}for i in nonmn],'RG_E_voltage_difference_in_electrical_fixture_max_mV':rgerror,'other_parameters_and_learned_gain_buffers_exact':True,'metadata':metadata['MN_only']},'full_source_diagnostic':{'POSTURE_source_max_mV':max(sourcevolt['full_source_diagnostic']),'actual_POSTURE_to_MN_current_max_nA':fullcurrent,'scope':'Separate diagnostic; removes source drive to RG as well as MNs,not the selected duplicate-MN repair. Switchingduring a live run would permit existing source voltage to decay,not hard-reset it.'},'direct_POST_inputs_preserved':len(postports),'actual_physics_model_arrays_unchanged':len(arrays),'archived_model_sha256':before,'body_steps':0,'training':False,'scope':'Synthetic electrical fixture uses actual saved directPOST currents,POSTURE1nA,andotherinputs0. This verifies source-policy consumption/cache ownership,not standing or walking.'}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps({'output':str(out),'default_false_error':maxfalse,'MN_only_current':MNcurrent,'MN_edges':len(ids),'RG_edges_preserved':len(nonmn)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.root,a.out)
