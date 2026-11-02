"""Current-code static reconstruction of dissertation prune/full/WALK recipes.

Extracts only loader functions from source AST: never imports/calls runner or
curriculum main, never advances SNS, and never loads a MuJoCo model.
Historical JSONL metrics remain evidence of their recorded run, not a new run.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,ast,copy,hashlib,json,sys,types,math
from pathlib import Path
import xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent;S=HERE.parent;sys.path[:0]=[str(HERE),str(S),str(S/'w2l_cpg')]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def extract(path,names,environment):
    tree=ast.parse(path.read_text(encoding='utf-8-sig'));selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    if {n.name for n in selected}!=set(names):raise ValueError('Missing authoritative loader')
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),environment)
    return types.SimpleNamespace(**{n:environment[n] for n in names})


def static_profile_check(network,family,P):
    """Optional NEW FSA candidate,with only algebra/cache checks; no steps."""
    import numpy as np,graded_profile as GP
    baseline=GP.field_digest(network.compiled);n=copy.deepcopy(network);c=n.compiled;caps=c.c_m.copy();ions={k:getattr(c,k).copy() for k in GP.ION_FIELDS if hasattr(c,k)}
    try:
        p=GP.default_profile(family)
        if P.G.get('gamma_coact',0) or P.G.get('vest_prop',0):raise ValueError('Caller gamma/vest modulation needs its own bound before profile check')
        maxima={'Ia':P.AFF['ia_gain']*(1+P.MOD['ia']),'II':P.AFF['i0_ii']+P.AFF['ii_gain']*(1+P.MOD['ii']),'Ib':.5*P.AFF['ib_gain']*1.5}
        p['physical_input_domains_nA']={name:[0.,maxima[name.split('_')[0]]] for name in n.inputs if name.startswith(('Ia_','II_','Ib_'))}
        ledger=GP.apply_profile(n,p,family=family);rows=[]
        for pool in ledger['pools']:
            dst=pool['post_index'];src=np.array(pool['presynaptic_indices']);voltage=pool['source_voltage_mV'];active=c.g_max_non[dst,src]*[GP.activation_fraction(c,dst,int(q),u) for q,u in zip(src,voltage)]
            steady=(pool['bias_current_for_probe_nA']+sum(active*c.del_e[dst,src]))/(c.g_m[dst]+sum(active));error=abs(steady-pool['target_post_mV'])
            if error>1e-12:raise AssertionError('Actual compiled isolated-pool contract failed')
            rows.append({'kind':pool['kind'],'post':pool['post'],'fan_in':pool['fan_in'],'static_error_mV':float(error)})
        names=GP.population_names(n);nonencoder=[i for i,name in enumerate(names) if not name.startswith(('Ia_','II_','Ib_'))]
        if not np.array_equal(c.c_m[nonencoder],caps[nonencoder]) or any(not np.array_equal(getattr(c,k),v) for k,v in ions.items()):raise AssertionError('Intrinsic RG/PF/MN or NaP dynamics changed')
        theta=2*math.pi*p['context']['kinematic_bandwidth_hz']*c.dt;tau=ledger['encoder_tau_s'];pole=1-c.dt/tau;phase=math.degrees(math.atan2(math.sin(theta),math.cos(theta)-pole))
        upper=np.zeros(len(n.inputs));indices={name:i for i,name in enumerate(n.inputs)}
        for name,interval in p['physical_input_domains_nA'].items():upper[indices[name]]=interval[1]
        if not np.array_equal(GP.adapt_inputs(n,upper)[0],upper) or not np.array_equal(GP.adapt_inputs(n,upper*1000,unit='pA')[0],upper):raise AssertionError('Physical/unit endpoint adapter changed current')
        for kind in ('Ia','II','Ib'):
            bad=upper.copy();name=next(q for q in n.inputs if q.startswith(kind+'_'));bad[indices[name]]+=1
            try:GP.adapt_inputs(n,bad)
            except ValueError:pass
            else:raise AssertionError('Out-of-domain afferent current accepted')
        import central_afferent_budget as CB
        central_proposal=CB.propose(n,CB.default_profile())
        central_inventory={'proposed_only_not_applied':True,'qualifying_groups':len(central_proposal['groups']),'qualifying_existing_edges':sum(len(q['edges']) for q in central_proposal['groups']),'unmatched_groups':central_proposal['unmatched_groups'],'unqualified_direct_edges':central_proposal['unqualified_direct_edges'],'composite_inventory':central_proposal['composite_inventory'],'groups':[{'target':q['target'],'kind':q['kind'],'sign':q['sign'],'fan_in':len(q['edges']),'old_reference_maximum_current_nA':q['old_maximum_reference_current_nA'],'optional_new_reference_budget_nA':q['new_maximum_reference_current_nA']} for q in central_proposal['groups']],'scope':'Existing direct routes only. Composite/missing paths remain unresolved; no central-current budget or new source route applied.'}
        from sparse_pairwise_eval import SparsePairwiseForward
        n.compiled=SparsePairwiseForward(c);verification=GP.verify_profile(n)
        if GP.field_digest(network.compiled)!=baseline:raise AssertionError('Historical baseline modified')
        return {'compatible':True,'candidate_scope':'NEW opt-in graded FSA candidate cloned from actual dissertation full/pruned recipe; not a changed historical record','profile':p,'ledger':ledger,'static_pools_checked':len(rows),'static_pool_error_max_mV':max(q['static_error_mV'] for q in rows),'pool_rows':rows,'source_physical_muscle_bounds_nA':maxima,'bound_scope':'Muscle ports only; clipped length/velocity/force source formulas,default gamma0/vest_prop0. Not a full command/contact domain certificate.','AFF':dict(P.AFF),'MOD':dict(P.MOD),'all_muscle_upper_endpoint_identity':True,'pA_to_nA_identity':True,'outside_each_afferent_type_rejected':True,'discrete_encoder_phase_deg_at_design_band':phase,'nonencoder_clocks_and_NaP_parameters_unchanged':True,'final_sparse_parameter_cache':verification,'sparse_dynamic_parity_claimed':False,'baseline_parameter_sha256':baseline,'baseline_unchanged':True,'central_current_inventory_and_unapplied_proposal':central_inventory,'new_mirror_CB_coincidence_or_heel_policy_applied':False,'global_training_release':False,'neural_steps':0,'physics_steps':0}
    except (ValueError,AssertionError) as e:
        return {'compatible':False,'failure_type':type(e).__name__,'failure':str(e),'stopped_without_adapting':True,'baseline_unchanged':GP.field_digest(network.compiled)==baseline,'neural_steps':0,'physics_steps':0}


def audit(family,condition='full',profile_checks=False):
    if family=='s3k':os.environ.pop('AARL_NET',None)
    else:os.environ['AARL_NET']=family
    import params as P,build_network as B
    from gait_topology_audit import raw
    recipe_file=S/({'s3k':'reports_20260923/s3k_trial34_full_params.json','syn6':'curriculum_syn6_stage5.json','w2lvar':'curriculum_w2lvar_stage5.json'}[family]);t34_file=S/'reports_20260923/s3k_trial34_full_params.json'
    stagefile=S/({'s3k':'_curriculum.py','syn6':'_curriculum_syn6.py','w2lvar':'_curriculum_w2lvar.py'}[family]);basefile=S/'best_walk_params_v10.json'
    baseline=json.loads(basefile.read_text())['multipliers'];baseline=copy.deepcopy(baseline);baseline['renshaw']=0. if family=='syn6' else .5
    if family=='syn6':baseline.update(syn6=1.,syn6_brainstem=0.)
    recipe={**baseline,**json.loads(recipe_file.read_text())['params']}
    owfile=S/'optuna_walk.py';env={'params':P,'HERE':S,'json':json,'FIT':None,'BASE':{}}
    OW=extract(owfile,('load_fitted_baseline','set_params'),env)
    stage_env={'OW':OW,'P':P,'os':os};loader=extract(stagefile,('set_stage',),stage_env)
    cuts={};interleg=True;optimization_record=None;executed=['load_fitted_baseline','set_params','set_stage']
    if condition.startswith('retuned_walk'):
        if family!='s3k':raise ValueError('Only actual S3K pruned-retune source defined')
        number=condition[len('retuned_walk')];recipe_file=S/('curriculum_pruned_walk_'+number+'_fix_20261003.json');saved=json.loads(recipe_file.read_text());t34=json.loads(t34_file.read_text())['params'];prunefile=S/'_curriculum_pruned.py';tree=ast.parse(prunefile.read_text(encoding='utf-8-sig'))
        penv={'T34':t34}
        for key in ('PIN','PRUNE_KEYS','PRUNE_ARGV'):
            selected=[q for q in tree.body if isinstance(q,ast.Assign) and any(isinstance(t,ast.Name) and t.id==key for t in q.targets)]
            if len(selected)!=1:raise ValueError('Unrecognized pruned configuration constant '+key)
            exec(compile(ast.Module(body=selected,type_ignores=[]),str(prunefile),'exec'),penv)
        if saved['pin']!=penv['PIN'] or saved['prune_keys']!=penv['PRUNE_KEYS'] or saved['prune_argv']!=penv['PRUNE_ARGV']:raise ValueError('Saved pruned recipe no longer matches actual source overlay')
        config_env={**penv,'BASE':baseline,'CS3':loader,'P':P};config=extract(prunefile,('_config',),config_env);argv=config._config(saved['params'])
        if argv!=['--no-interleg']:raise ValueError('Unrecognized pruned interleg command')
        recipe={**baseline,**{k:float(v) for k,v in t34.items()},**penv['PIN'],**{k:float(v) for k,v in saved['params'].items()}};cuts=penv['PRUNE_KEYS'];interleg=False;executed.append('_config')
        optimization_record={'file':str(recipe_file),'score':saved['score'],'trial':saved.get('trial'),'study':saved.get('study'),'metric_scope':'Saved optimization objective,not asserted to be a raw full/WALK kine_score or independently reproduced run'}
    else:loader.set_stage(5,dict(recipe))
    if condition=='combo':
        if family!='s3k':raise ValueError('Only the actual S3K combo recipe is defined')
        tree=ast.parse((S/'prune_combo.py').read_text(encoding='utf-8-sig'));matches=[q.value for q in tree.body if isinstance(q,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='COMBO_KEYS' for t in q.targets)]
        if len(matches)!=1:raise ValueError('Unrecognized literal combo recipe')
        cuts=ast.literal_eval(matches[0]);P.G.update(cuts);interleg=False
    before_G=copy.deepcopy(P.G)
    # Evaluate only runner's MODEL assignment with a deliberately empty override
    # environment: the historical default,not a possibly unrelated live override.
    rfile=S/'runner.py';tree=ast.parse(rfile.read_text(encoding='utf-8-sig'));assignment=[n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='MODEL' for t in n.targets)]
    if len(assignment)!=1:raise ValueError('Unrecognized source MODEL declaration')
    model_env={'Path':Path,'_os':types.SimpleNamespace(environ={})};exec(compile(ast.Module(body=assignment,type_ignores=[]),str(rfile),'exec'),model_env);model=model_env['MODEL']
    acts=[a.attrib['name'] for a in ET.parse(model).getroot().find('actuator')]
    n=B.build(acts,dt=P.DT,interleg=interleg);summary=raw(n);profile=static_profile_check(n,family,P) if profile_checks else {'requested':False}
    trace_file=S/('prune_results_'+family+'.jsonl');records=[json.loads(line) for line in trace_file.read_text().splitlines() if line.strip()];records=[r for r in records if r.get('config')==condition and r.get('mode')=='WALK']
    if not records and optimization_record is None:raise ValueError('Missing original full/WALK record')
    relevant=('full_rules','ia_in','renshaw','heel_rge','toe_rge','ib_rge','ia_f_central','ii_f_central','ii_e_central','ib_e_central','joint_pf','heel_in_f_exc','ia_pf_f','ii_pf_f','contact_onset','contra_swing','c1_gain','v3_gain')
    if 'mujoco' in sys.modules or 'runner' in sys.modules:raise AssertionError('Source-only audit unexpectedly imported physics/runner')
    sources=[Path(__file__),S/'prune_matrix.py',S/'prune_combo.py',S/'_curriculum_pruned.py',stagefile,owfile,S/'fitted_walk_params.json',basefile,t34_file,recipe_file,trace_file,rfile,S/'params.py',S/'build_network.py',S/'build_network_w2lvar.py',S/'build_network_syn6.py',S/'muscle_map.py',S/'synergy_basis.npz',HERE/'graded_profile.py',HERE/'central_afferent_budget.py',S/'w2l_cpg/sparse_pairwise_eval.py',model]
    return {'role':'Gait neural pathway composition reviewer','family':family,'condition':'Dissertation prune/'+condition+'/WALK recipe,current-code static reconstruction',
            'historical_metrics_records':records,'saved_optimization_record':optimization_record,'historical_metrics_scope':'Original saved result,not reproduced or re-scored here. JSONL lacks historical source/model/config hashes; current source/model hashes do not prove historical byte identity.',
            'recipe':{'merge':'best_walk_params_v10.multipliers plus family source recipe and actual extracted loader; retunes additionally apply actual _config/PIN/prune overlay','recipe_file':str(recipe_file),'merged_parameters':recipe,'combo_cuts':cuts,'stage_loader':str(stagefile),'source_function_execution_only':executed,'model_policy':'Current runner default declaration with empty AARL_MODEL override environment; historical override not recorded','model':str(model),'actuators':len(acts),'interleg_constructor':interleg},
            'new_optional_FSA_static_compatibility':profile,'current_runner_connectome_gains_overlay_present':(S/'connectome_gains.json').exists(),'historical_external_overlays_or_environment_byte_proven':False,
            'actual_neural_topology':summary,'G_after_stage_before_constructor':{k:before_G[k] for k in relevant},'G_after_constructor':{k:P.G[k] for k in relevant},'constructor_G_changes':{k:{'before':before_G[k],'after':P.G[k]} for k in before_G if before_G[k]!=P.G[k]},
            'neural_steps':0,'physics_steps':0,'runner_imported':False,'mujoco_imported':False,'training':False,'historical_artifacts_edited':False,'source_and_artifact_sha256':{str(p):sha(p) for p in sources}}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=('s3k','syn6','w2lvar'),required=True);p.add_argument('--condition',choices=('full','combo','retuned_walk1_fix','retuned_walk2_fix'),default='full');p.add_argument('--profile-checks',action='store_true');p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=audit(a.family,a.condition,a.profile_checks);a.out.mkdir(exist_ok=True,parents=True);f=a.out/(a.family+'_'+a.condition+'_DISSERTATION_PRUNE_TOPOLOGY.json');f.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps({'family':a.family,'condition':a.condition,'nodes':r['actual_neural_topology']['nodes'],'nonzero_pairs':r['actual_neural_topology']['nonzero_pairs'],'FSA':{k:v for k,v in r['new_optional_FSA_static_compatibility'].items() if k in ('compatible','failure','static_pools_checked','static_pool_error_max_mV','discrete_encoder_phase_deg_at_design_band')},'historical_WALK_scores':[v['metrics'].get('kine_score') for v in r['historical_metrics_records']],'optimization_record':r['saved_optimization_record'],'output':str(f)}))
