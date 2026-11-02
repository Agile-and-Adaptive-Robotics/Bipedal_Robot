"""Distinct Syn6 single-IN/FSA + reviewed mono22/extensor-II candidate.

Electrical composition only. Historical builders/profiles remain unchanged.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse,copy,hashlib,json,sys,tempfile
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import gait_design_composition as G
import mono_central_feedback as M
from sparse_pairwise_eval import SparsePairwiseForward
SCHEMA='aarl_syn6_mono22_composition_v1'
NAME='Syn6_singleInE_canonical_FSA_mono22_extII'


def make_profile(base):
    base=Path(base).resolve();overlay=M.default_profile('syn6')
    return dict(schema=SCHEMA,name=NAME,base_profile=dict(path=str(base),sha256=M.sha(base)),
        overlay=overlay,source_sha256={str(Path(p).resolve()):M.sha(p)for p in (__file__,M.__file__,G.__file__)},
        input_contract=dict(afferent_domains_by_type_nA={k:[0.,v]for k,v in overlay['source_current_max_nA'].items()},
            other_ports='Retain source interface; finite currents required. This composition adds no reader or new physiological bound for other input ports.',
            caller_conditions=copy.deepcopy(overlay['caller_conditions'])),
        biology_scope='Reviewed22 sagittal mono-role muscles only; extensorII is an explicit engineering extension. Direct central excitatoryIb bypasses the separate IBEXC MN phase/load gate. No biarticular or independent3D completion.',
        body_training_allowed=False)


def load_profile(value):
    p=json.loads(Path(value).read_text())if isinstance(value,(str,Path))else copy.deepcopy(value)
    G.GP.exact_keys(p,('schema','name','base_profile','overlay','source_sha256','input_contract','biology_scope','body_training_allowed'),'mono22 composed profile')
    if p!=make_profile(p['base_profile']['path']):raise ValueError('Profile/source/contract changed or unused field; create a separately reviewed named profile')
    if M.sha(p['base_profile']['path'])!=p['base_profile']['sha256']:raise ValueError('Base profile changed')
    return p


def verify(n):
    q=n.syn6_mono_composition
    return dict(base_components=G.verify(n,q['base_gate']),mono_overlay=M.verify(n,P=q['P']),
        exact_complete_graph_sha256=G.GP.field_digest(G.CB.core(n)))


def validate_inputs(n,u):
    u,_=G.GP.adapt_inputs(n,u,enforce_domain=True)
    return M.validate_input_currents(n,u,P=n.syn6_mono_composition['P'])


def operator_checks(n):
    c=G.CB.core(n);rows=[];maximum=0.;zero=np.zeros(len(n.inputs))
    for group in n.mono_central_feedback_ledger['groups']:
        post=group['edges'][0]['post']
        for fraction in (0.,.25,.5,1.):
            on=SparsePairwiseForward(copy.deepcopy(c));offcore=copy.deepcopy(c)
            for e in group['edges']:offcore.g_max_non[e['post'],e['pre']]=0.
            off=SparsePairwiseForward(offcore);on.c.V[post]=off.c.V[post]=2.5
            expected=0.
            for e in group['edges']:
                voltage=e['source_max_voltage_mV']*fraction;on.c.V[e['pre']]=off.c.V[e['pre']]=voltage
                expected+=e['gmax_uS']*np.clip(voltage/5.,0,1)*(8.-2.5)
            on.forward(zero);off.forward(zero)
            actual=float((on.c.V[post]-off.c.V[post])/c.time_factor_membrane[post]);error=abs(actual-expected);maximum=max(maximum,error)
            if error>1e-10:raise AssertionError('Actual composed central current differs from contract')
            rows.append(dict(target=group['target'],type=group['type'],fraction=fraction,expected_nA=expected,actual_nA=actual,error_nA=error))
    return dict(pools=len(n.mono_central_feedback_ledger['groups']),tests=rows,maximum_error_nA=maximum)


def build_named_profile(value):
    p=load_profile(value);n,base_report=G.build_named_profile(p['base_profile']['path']);import params as P
    c=G.CB.core(n);initial=G.states(c);before=G.GP.field_digest(c)
    fields={k:getattr(c,k).copy()for k in G.GP.FIELDS+G.GP.ION_FIELDS}
    baseline_pairs=set(zip(*np.nonzero(c.g_max_non)));old_cache_rows=np.c_[n.compiled.post,n.compiled.pre,n.compiled.gmax,n.compiled.reversal,n.compiled.elo,n.compiled.elo+n.compiled.erange].copy()
    no_op=M.install(n,None)
    if not no_op['default_noop']or G.GP.field_digest(c)!=before:raise AssertionError('Default overlay changes base')
    plan=M.install(n,p['overlay'],P=P);added={(r['post'],r['pre'])for g in plan['groups']for r in g['edges']}
    if added&baseline_pairs:raise AssertionError('Overlay duplicates an existing active pair')
    selected=np.zeros(c.g_max_non.shape,dtype=bool)
    for post,pre in added:selected[post,pre]=True
    numerical_changes={}
    for key,old in fields.items():
        actual=getattr(c,key)
        if key in ('g_max_non','del_e','e_lo','e_hi'):
            if not np.array_equal(actual[~selected],old[~selected]):raise AssertionError('An old operator changed: '+key)
            numerical_changes[key]=int(np.count_nonzero(actual!=old))
        elif not np.array_equal(actual,old):raise AssertionError('Other neural parameter changed: '+key)
    if any(not np.array_equal(getattr(c,k),v)for k,v in initial.items()):raise AssertionError('Dynamic state changed during composition')
    # Explicitly qualify a DISJOINT extension of the full-matrix fingerprint.
    # No old parameter is edited and no checksum guard is silently discarded.
    gp={(r['post_index'],r['pre_index'])for r in n.graded_fsa_ledger['modified_edges']}
    cb={(r['post'],r['pre'])for g in n.central_afferent_budget_design['groups']for r in g['edges']}
    if added&(gp|cb):raise AssertionError('New central map overlaps previously budgeted operators')
    final=G.GP.field_digest(c)
    n.graded_fsa_ledger['approved_disjoint_composition_extension']=dict(base_sha256=before,final_sha256=final,new_pairs=len(added),all_old_operators_exact=True,local_and_existing_central_intersection=0)
    n.graded_fsa_ledger['compiled_after_sha256']=final
    gate=base_report['composition']['coincidence_design']
    n.syn6_mono_composition=dict(profile=p,P=P,base_gate=gate)
    checks=verify(n);ops=operator_checks(n);pools=G.GP.probe_pools(n)
    runtime=G.domain_and_backend_probe(n,gate)
    # Exercise real source step interface with explicit command changes.
    for command in (0.,2.5,0.):
        for tick in range(16):
            u=n.make_inputs();u[n.input_index('DRIVE')]=command
            u=validate_inputs(n,u);n.step(u);verify(n)
    invalid=n.make_inputs();first=plan['groups'][0]['edges'][0];invalid[n.input_index(first['input_port'])]=first['source_max_current_nA']/first['input_weight']+.01
    saved_state=c.V.copy()
    try:validate_inputs(n,invalid);rejected=False
    except ValueError:rejected=True
    if not rejected or not np.array_equal(saved_state,c.V):raise AssertionError('Invalid source current not rejected before mutation')
    # Cache-corruption and duplicate-installer guards on the final graph.
    row=plan['groups'][0]['edges'][0];slot=np.flatnonzero((n.compiled.post==row['post'])&(n.compiled.pre==row['pre']))[0];old=n.compiled.gmax[slot];n.compiled.gmax[slot]*=1.01
    try:verify(n);cache_rejected=False
    except (ValueError,AssertionError):cache_rejected=True
    finally:n.compiled.gmax[slot]=old
    if not cache_rejected:raise AssertionError('Corrupted final cache accepted')
    try:M.install(n,p['overlay'],P=P);duplicate_rejected=False
    except ValueError:duplicate_rejected=True
    if not duplicate_rejected:raise AssertionError('Duplicate overlay accepted')
    for key,v in initial.items():getattr(c,key)[:]=v
    n.named_composed_gait_profile=copy.deepcopy(p)
    report=dict(profile=p,base_profile_sha256=p['base_profile']['sha256'],base_census=dict(neurons=len(c.V),inputs=len(n.inputs),active_pairs=len(baseline_pairs)),
        complete_census=dict(neurons=len(c.V),inputs=len(n.inputs),active_pairs=int(np.count_nonzero(c.g_max_non)),new_central_pairs=len(added)),
        all_old_operators_exact=True,other_neuron_ionic_membrane_input_parameters_exact=True,numerical_changes=numerical_changes,
        approved_disjoint_fingerprint_extension=n.graded_fsa_ledger['approved_disjoint_composition_extension'],
        current_operators=ops,terminal_pool_checks=dict(pools=pools['pools_checked'],maximum_error_mV=pools['all_pool_static_error_max_mV']),
        source_domains=copy.deepcopy(n.graded_fsa_profile['physical_input_domains_nA']),actual_runtime_and_backend=runtime,
        duplicate_install_rejected=duplicate_rejected,cache_corruption_rejected=cache_rejected,out_of_domain_rejected_before_mutation=rejected,
        command_persistence_checked=[0.,2.5,0.],initial_dynamic_states_restored=True,
        final_verification=verify(n),coverage=G.coverage(n,p['name']),mono_rule_matrix=plan,
        direct_central_Ib_bypass=[r for g in plan['groups']for r in g['edges']if r['type']=='Ib'],
        mechanism_scope='New direct excitatory centralIb routes are NOT phase/load-gated IBEXC MN outputs; both mechanisms coexist. ExtensorII is an engineering extension. No new biarticular circuit,4-cell change,fullmuscle or3D closure.',
        physics_steps=0,body_bindings=0,training=False)
    # Persist enforcement for actual later electrical callers.
    original_step=n.step
    def checked_step(u):
        u=validate_inputs(n,u);answer=original_step(u);verify(n);return answer
    n.step=checked_step
    return n,report


def export(n,report,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if(out/'FINAL_COMPILED_GRAPH.json').exists():raise FileExistsError('Preserve completed candidate; choose another output')
    # Reuse the reviewed complete exporter in a scoped temporary directory.
    with tempfile.TemporaryDirectory(prefix='syn6_mono_export_',dir=out)as tmp:
        graph=G.export_named_graph(n,report,Path(tmp))
    graph.pop('claimed_missing_central_Ia_II_repair',None)
    new={(r['post'],r['pre'])for g in n.mono_central_feedback_ledger['groups']for r in g['edges']}
    for e in graph['chemical_edges']:e['origin']='NEW_REVIEWED_MONO22_CENTRAL'if(e['target_index'],e['source_index'])in new else'RETAINED_COMPOSED_BASE'
    graph.update(scoped_mono22_central_feedback_added=True,full_all_muscle_central_completion=False,
        coverage=report['coverage'],mono_rule_matrix=n.mono_central_feedback_ledger,
        direct_central_Ib_bypass_count=len(report['direct_central_Ib_bypass']),mechanism_scope=report['mechanism_scope'])
    M.write(out/'FINAL_COMPILED_GRAPH.json',graph);M.write(out/'COMPOSED_VERIFICATION.json',report)
    return dict(graph_sha256=M.sha(out/'FINAL_COMPILED_GRAPH.json'),proof_sha256=M.sha(out/'COMPOSED_VERIFICATION.json'),census=report['complete_census'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--base-profile',type=Path);p.add_argument('--profile',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if(a.profile is None)==(a.base_profile is None):raise ValueError('Choose create from base or rebuild profile')
    a.out.mkdir(parents=True,exist_ok=True)
    if a.base_profile:
        profile=make_profile(a.base_profile);target=a.out/'PROFILE.json'
        if target.exists():raise FileExistsError('Profile exists; do not overwrite')
        M.write(target,profile)
    else:profile=a.profile
    n,report=build_named_profile(profile);print(json.dumps(export(n,report,a.out)))
