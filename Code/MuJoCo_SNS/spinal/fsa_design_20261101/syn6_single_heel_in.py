"""Ben's explicit Syn6 repair: retain heel->ipsi InE, remove heel->ipsi InF.

Default-off; install before acceleration and other full-graph design profiles.
No new connection, source builder change, input encoder, or physics operation.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ.setdefault(_key,'1')
import argparse,copy,hashlib,json,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import graded_profile as G

SCHEMA='aarl_syn6_single_heel_IN_v1'


def policy():
    return dict(schema=SCHEMA,family='syn6',hypothesis='heel_ipsilateral_InE_only',
                event_contract='Preserve caller encoding; pulse versus held contact must be declared by the caller.')


def install(net,value=None):
    if value is None:return dict(applied=False,default_noop=True)
    G.assert_campaign_allowed('optional_profile_implementation')
    p=json.loads(Path(value).read_text())if isinstance(value,(str,Path))else copy.deepcopy(value)
    G.exact_keys(p,('schema','family','hypothesis','event_contract'),'Syn6 single-IN policy')
    if p!=policy():raise ValueError('Only the explicitly selected single-InE policy is supported')
    if type(net).__module__.split('.')[-1]!='build_network_syn6':raise ValueError('Actual Syn6 constructor required')
    if hasattr(net,'syn6_single_heel_in'):raise ValueError('Single-IN repair already installed')
    if hasattr(net,'_heel_stance_policy_record'):raise ValueError('Another heel hypothesis is already installed')
    if hasattr(net,'graded_fsa_ledger'):raise ValueError('Install before graded/full-graph design profiles, then reverify their composition')
    c=G.raw_compiled(net)
    if net.compiled is not c or c.spiking or c.electrical or c.electrical_rectified:
        raise ValueError('Actual graded dense graph required before accelerator construction')
    names=G.population_names(net);ix={name:i for i,name in enumerate(names)}
    required=[kind+'_'+s for s in('r','l')for kind in('HEEL','InE','InF','RG_E','RG_F')]
    if any(name not in ix for name in required):raise ValueError('Missing source cells; none are invented')
    changed=[];keep=[]
    for side in('r','l'):
        src=ix['HEEL_'+side]
        for target in('r','l'):
            for phase in('E','F'):
                if c.g_max_non[ix['RG_'+phase+'_'+target],src]!=0.:
                    raise ValueError('Existing direct heel-to-RG path would stack; separate composition review required')
            if target!=side:
                for kind in('InE','InF'):
                    if c.g_max_non[ix[kind+'_'+target],src]!=0.:
                        raise ValueError('Existing crossed heel-to-IN policy would stack')
        for relay,rg in(('InE','RG_F'),('InF','RG_E')):
            if c.g_max_non[ix[rg+'_'+side],ix[relay+'_'+side]]<=0 or c.del_e[ix[rg+'_'+side],ix[relay+'_'+side]]>=0:
                raise ValueError('Source inhibitory relay endpoint/reversal differs')
        rows=[]
        for relay in('InE','InF'):
            dst=ix[relay+'_'+side]
            fields={name:float(getattr(c,name)[dst,src])for name in('g_max_non','del_e','e_lo','e_hi')}
            if fields['g_max_non']<=0 or fields['del_e']<=0 or fields['e_hi']<=fields['e_lo']:
                raise ValueError('Expected original positive heel-to-both-IN source policy is absent')
            rows.append(dict(source=names[src],target=names[dst],source_index=src,target_index=dst,**fields))
        keep.append(rows[0]);changed.append(rows[1])
    if any(keep[0][k]!=keep[1][k]or changed[0][k]!=changed[1][k]for k in('g_max_non','del_e','e_lo','e_hi')):
        raise ValueError('Source left/right heel fields differ; do not silently fit or mirror')
    before=G.field_digest(c)
    arrays={key:v.copy()for key,v in c.__dict__.items()if isinstance(v,np.ndarray)}
    for row in changed:c.g_max_non[row['target_index'],row['source_index']]=0.
    allowed={(r['target_index'],r['source_index'])for r in changed}
    for key,old in arrays.items():
        new=getattr(c,key)
        if key=='g_max_non':
            if set(zip(*np.nonzero(old!=new)))!=allowed:raise AssertionError('Unexpected conductance delta')
        elif not np.array_equal(old,new,equal_nan=True):raise AssertionError('Out-of-scope compiled field changed: '+key)
    record=dict(schema=SCHEMA,policy=p,applied=True,changed=changed,retained=keep,
                before_sha256=before,after_sha256=G.field_digest(c),
                only_two_conductances_zeroed=True,all_other_compiled_arrays_exact=True,
                direct_heel_RG_added=False,heel_PF_paths_preserved=True,
                encoding_and_clocks_preserved=True,physics_changed=False,
                sign_scope='InE inhibits RG-F in the declared operating regime; currents follow g*s*(Erev-Vpost), not a fixed sign at all voltages.')
    net.syn6_single_heel_in=record
    return record


def verify(net):
    record=getattr(net,'syn6_single_heel_in',None)
    if record is None:return dict(applied=False)
    c=G.raw_compiled(net)
    for row in record['changed']:
        if c.g_max_non[row['target_index'],row['source_index']]!=0.:
            raise ValueError('Offending heel-to-InF branch was restored')
    for row in record['retained']:
        for key in('g_max_non','del_e','e_lo','e_hi'):
            if float(getattr(c,key)[row['target_index'],row['source_index']])!=row[key]:
                raise ValueError('Retained heel-to-InE source fields changed')
    wrapper=net.compiled
    if hasattr(wrapper,'pre'):
        for row in record['changed']:
            slots=(wrapper.pre==row['source_index'])&(wrapper.post==row['target_index'])
            if np.any(wrapper.gmax[slots]!=0):raise ValueError('Sparse wrapper retained a nonzero offending branch')
    return dict(applied=True,offending_dense_and_sparse_current_exact_zero=True,
                retained_source_fields_exact=True)


def source_fixture():
    # Capture the real full-recipe graph at the existing audit's callback boundary.
    # The callback applies no profile and advances no neurons.
    import gait_dissertation_topology_audit as D
    captured=[];original=D.static_profile_check
    def capture(network,family,params):
        captured.append(network);return {'requested':False,'network_capture_only':True}
    D.static_profile_check=capture
    try:recipe=D.audit('syn6','full',profile_checks=True)
    finally:D.static_profile_check=original
    if len(captured)!=1:raise AssertionError('Source recipe was not captured exactly once')
    return captured[0],recipe


def proof(out):
    from sparse_pairwise_eval import SparsePairwiseForward
    from heel_stance_policy import isolated_static_operator
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if(out/'RESULTS.json').exists():raise FileExistsError('Preserve completed proof; use another explicit output folder')
    source,recipe=source_fixture();before=G.field_digest(source.compiled)
    noop=install(source,None)
    if not noop['default_noop']or before!=G.field_digest(source.compiled):raise AssertionError('Default changed the source')
    n=copy.deepcopy(source);ledger=install(n,policy());c=n.compiled
    names=G.population_names(n);ix={name:i for i,name in enumerate(names)}
    operators=isolated_static_operator(n,.1,2.5)
    for row in operators:
        side=row['heel_source'][-1]
        for path in row['two_hop']:
            expected=path['target']=='RG_F_'+side
            if(expected and path['delta_target_current_nA']>=0)or(not expected and path['delta_target_current_nA']!=0):
                raise AssertionError('Unexpected isolated relay sign or target')
    current_tests=[]
    for row in ledger['changed']:
        pre,post=row['source_index'],row['target_index']
        for vpre in(0.,2.5,5.,10.):
            for vpost in(-2.5,0.,2.5,8.,10.):
                current=c.g_max_non[post,pre]*G.activation_fraction(c,post,pre,vpre)*(c.del_e[post,pre]-vpost)
                if current!=0.:raise AssertionError('Zero branch emits current')
                current_tests.append(dict(source=row['source'],target=row['target'],pre_mV=vpre,post_mV=vpost,current_nA=float(current)))
    # One actual dense update isolates the retained versus removed heel edges.
    on=copy.deepcopy(c);off=copy.deepcopy(c);on.g_max_non[:]=0.;off.g_max_non[:]=0.
    for row in ledger['retained']:
        pre,post=row['source_index'],row['target_index'];on.g_max_non[post,pre]=row['g_max_non']
        on.V[pre]=off.V[pre]=5.;on.V[post]=off.V[post]=2.5
    for row in ledger['changed']:on.V[row['target_index']]=off.V[row['target_index']]=2.5
    u=np.zeros(len(n.inputs));on.forward(u);off.forward(u);actual=[]
    for row in ledger['retained']+ledger['changed']:
        post=row['target_index'];measured=float((on.V[post]-off.V[post])/on.time_factor_membrane[post])
        expected=row['g_max_non']*(row['del_e']-2.5)if row in ledger['retained']else 0.
        if abs(measured-expected)>1e-12:raise AssertionError('Actual dense current operator mismatch')
        actual.append(dict(source=row['source'],target=row['target'],expected_nA=expected,actual_nA=measured))
    dense=copy.deepcopy(c);sparse=SparsePairwiseForward(copy.deepcopy(c));trace=[];error=0.
    for step in range(128):
        u=n.make_inputs();u[n.input_index('DRIVE')]=2.5
        if 16<=step<48:u[n.input_index('HEEL_c_r')]=.1
        if 64<=step<96:u[n.input_index('HEEL_c_l')]=.1
        dense.forward(u);sparse.forward(u)
        error=max(error,float(np.max(abs(dense.V-sparse.V))))
        if not np.isfinite(dense.V).all()or not np.isfinite(sparse.V).all():raise AssertionError('Nonfinite electrical fixture')
        trace.append(sparse.V.copy())
    if error!=0.:raise AssertionError('Dense/sparse traces differ')
    n.compiled=sparse;final=verify(n)
    # Negative guard: a separate direct reset must not be silently combined.
    poisoned=copy.deepcopy(source);poisoned.compiled.g_max_non[ix['RG_F_l'],ix['HEEL_r']]=.01
    poison_before=G.field_digest(poisoned.compiled)
    try:install(poisoned,policy())
    except ValueError:pass
    else:raise AssertionError('Stacked direct reset was accepted')
    if poison_before!=G.field_digest(poisoned.compiled):raise AssertionError('Rejected composition mutated the graph')
    if before!=G.field_digest(source.compiled):raise AssertionError('Historical source changed')
    pre=sparse.pre;post=sparse.post
    np.savez_compressed(out/'COMPILED_GRAPH.npz',names=np.array(names),pre=pre,post=post,gmax=sparse.gmax,reversal=sparse.reversal,elo=sparse.elo,erange=sparse.erange,c_m=c.c_m,g_m=c.g_m,V_rest=c.V_rest,i_b=c.i_b)
    trace=np.array(trace);np.savez_compressed(out/'NEURAL_PARITY.npz',voltage=trace,dt_s=c.dt,names=np.array(names))
    result=dict(schema=SCHEMA,fixture='Current-code reconstruction of dissertation Syn6 full recipe; repaired experimental copy only',
                source_recipe=recipe['recipe'],source_hashes=recipe['source_and_artifact_sha256'],helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                source_nodes=len(names),source_nonzero_pairs=int(np.count_nonzero(source.compiled.g_max_non)),
                final_nonzero_pairs=len(sparse.gmax),ledger=ledger,default_none_exact=True,historical_source_exact=True,
                actual_dense_isolated_current_tests=actual,zero_branch_grid=current_tests,two_hop_static=operators,
                neural_steps=128,physics_steps=0,dense_sparse_max_error_mV=error,final_verification=final,
                trace_sha256=hashlib.sha256(trace.tobytes()).hexdigest(),stacked_direct_path_rejected_before_mutation=True,
                event_encoder_changed=False,body_or_training_claim=False,
                limitation='Electrical current/backend proof; no full-network phase-response, locomotion, or combined FSA profile acceptance.')
    (out/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({key:result[key]for key in('source_nodes','source_nonzero_pairs','final_nonzero_pairs','neural_steps','physics_steps','dense_sparse_max_error_mV','trace_sha256')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    proof(parser.parse_args().out)
