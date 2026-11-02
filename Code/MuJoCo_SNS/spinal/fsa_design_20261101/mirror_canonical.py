"""Explicit right-canonical graded syn6 candidate; no learned factor alignment.

Copy actual right-side electrical parameters and target-owned connections.
Never copy voltage, gates, initialization, time, or phase. Shared node
parameters/internal edges stay fixed; shared targets use right-side sources.
The original K6 file, builder, and default compiled model remain unchanged.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
if __name__=='__main__':
    for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np

HERE=Path(__file__).resolve().parent
VECTOR_PARAMETERS=('c_m','g_m','V_rest','i_b','time_factor_membrane','theta_0','theta_leak','m','tau_theta','time_factor_threshold','V_reset')
ION_PARAMETERS=('g_ion','e_ion','pow_a','k_a','slope_a','e_a','pow_b','k_b','slope_b','e_b','tau_max_b','pow_c','k_c','slope_c','e_c','tau_max_c')
MATRIX_PARAMETERS=('g_max_non','del_e','e_lo','e_hi')
DYNAMIC_FIELDS=('V','V_last','V_0','b_gate','b_gate_last','c_gate','c_gate_last','theta','theta_last','spikes','g_spike','spike_buffer')
ALIASES={'BILAT_R':'BAL_LAT_R','BILAT_L':'BAL_LAT_L'}


def sibling(name):
    for suffix,other in (('_r','_l'),('_l','_r'),('_R','_L'),('_L','_R')):
        if name.endswith(suffix):return name[:-2]+other
    return name


def side(name):
    if name.endswith(('_r','_R')):return 'r'
    if name.endswith(('_l','_L')):return 'l'
    return None


def binding(names,*,input_ports=False):
    normalized=[ALIASES.get(n,n) if input_ports else n for n in names]
    if len(set(normalized))!=len(normalized):raise ValueError('Ambiguous normalized bilateral names')
    lookup={n:i for i,n in enumerate(normalized)};missing=[n for n in normalized if sibling(n) not in lookup]
    if missing:raise ValueError(f'Unmatched bilateral binding; no invented node/port: {missing}')
    permutation=np.array([lookup[sibling(n)] for n in normalized],int)
    if not np.array_equal(permutation[permutation],np.arange(len(names))):raise ValueError('Bilateral permutation is not an involution')
    right=np.array([i for i,n in enumerate(normalized) if side(n)=='r'],int);left=permutation[right];mid=np.array([i for i,n in enumerate(normalized) if side(n) is None],int)
    return dict(names=normalized,permutation=permutation,right=right,left=left,midline=mid,aliases=[dict(source=a,normalized=b) for a,b in zip(names,normalized) if a!=b])


def default_policy():
    return dict(schema='aarl_right_canonical_graded_v1',canonical_side='r',edge_ownership='right_destination_then_right_source_for_shared_targets',dynamic_state_policy='preserve_all',midline_policy='freeze_parameters_and_internal_edges',input_alias_policy='BILAT_R_L_to_BAL_LAT_R_L')


def load_policy(value):
    p=json.loads(Path(value).read_text()) if isinstance(value,(str,Path)) else copy.deepcopy(value)
    if p!=default_policy():raise ValueError('Unknown/unused canonical policy fields or unsupported policy choice')
    return p


def project_matrix(a,rows,columns):
    out=a.copy();out[rows['left'],:]=a[rows['right'],:][:,columns['permutation']]
    out[np.ix_(rows['midline'],columns['left'])]=a[np.ix_(rows['midline'],columns['right'])]
    return out


def parameter_errors(c,nodes,ports):
    p=nodes['permutation'];q=ports['permutation'];errors={}
    for key in VECTOR_PARAMETERS:
        if hasattr(c,key):errors[key]=float(np.max(abs(getattr(c,key)-getattr(c,key)[p])))
    for key in ION_PARAMETERS:
        if hasattr(c,key):errors[key]=float(np.max(abs(getattr(c,key)-getattr(c,key)[...,p])))
    for key in MATRIX_PARAMETERS:errors[key]=float(np.max(abs(getattr(c,key)-getattr(c,key)[np.ix_(p,p)])))
    errors['input_connectivity']=float(np.max(abs(c.input_connectivity-c.input_connectivity[np.ix_(p,q)])))
    return errors


def apply_canonical_right(spinal,policy=None):
    if policy is None:return dict(applied=False,reason='Default-off; no model field changed')
    try:
        from .graded_profile import population_names,field_digest,assert_campaign_allowed
    except ImportError:
        from graded_profile import population_names,field_digest,assert_campaign_allowed
    assert_campaign_allowed('optional_profile_implementation');chosen=load_policy(policy);c=spinal.compiled
    if hasattr(c,'c'):raise ValueError('Project dense compiled topology before constructing accelerator edge caches')
    if c.spiking or c.electrical or c.electrical_rectified:raise ValueError('This candidate covers graded chemical syn6; do not silently transfer to spiking/electrical models')
    if getattr(spinal,'graded_fsa_family',None)!='syn6':raise ValueError('Canonical-right candidate requires the actual syn6 family')
    if c.output_voltage_connectivity.shape[0]:raise ValueError('Explicit output-port binding required for nonempty output maps')
    names=population_names(spinal);nodes=binding(names);ports=binding(spinal.inputs,input_ports=True);right,left,mid=nodes['right'],nodes['left'],nodes['midline']
    for act,mi in spinal.muscles.items():
        other=sibling(act)
        if other not in spinal.muscles or spinal.muscles[other].base!=mi.base or spinal.muscles[other].groups!=mi.groups:raise ValueError(f'Unverified anatomical muscle homologue {act}')
    for r,l in zip(right,left):
        if type(spinal.net.populations[r]['type']) is not type(spinal.net.populations[l]['type']):raise ValueError(f'Neuron-type mismatch: {names[r]}/{names[l]}')
    states={k:getattr(c,k).copy() for k in DYNAMIC_FIELDS if hasattr(c,k)};identities={k:id(getattr(c,k)) for k in states};before=field_digest(c)
    old={k:getattr(c,k).copy() for k in VECTOR_PARAMETERS+ION_PARAMETERS+MATRIX_PARAMETERS+('input_connectivity',) if hasattr(c,k)}
    projected={};changed=[]
    for key in VECTOR_PARAMETERS:
        if key not in old:continue
        if old[key].shape!=(len(names),):raise ValueError(f'Unexpected neuron-parameter shape {key}')
        projected[key]=old[key].copy();projected[key][left]=old[key][right]
    for key in ION_PARAMETERS:
        if key not in old:continue
        if old[key].shape[-1]!=len(names):raise ValueError(f'Unexpected ionic parameter axis {key}')
        projected[key]=old[key].copy();projected[key][...,left]=old[key][...,right]
    for key in MATRIX_PARAMETERS:projected[key]=project_matrix(old[key],nodes,nodes)
    projected['input_connectivity']=project_matrix(old['input_connectivity'],nodes,ports)
    for key,a in projected.items():
        if key in MATRIX_PARAMETERS:continue
        indices=np.argwhere(a!=old[key]);changed.append(dict(field=key,entries_changed=int(len(indices))))
    g0,g1=old['g_max_non'],projected['g_max_non'];edges=[]
    for dst,src in zip(*np.nonzero(g0!=g1)):
        edges.append(dict(source=names[src],target=names[dst],source_index=int(src),target_index=int(dst),old_g_uS=float(g0[dst,src]),new_g_uS=float(g1[dst,src]),action='added' if g0[dst,src]==0 else 'removed' if g1[dst,src]==0 else 'replaced',ownership='right_destination_mirror' if dst in left else 'shared_target_right_source_mirror'))
    # All edits are static parameters. No initialization or phase state is copied.
    for key,a in projected.items():getattr(c,key)[:]=a
    try:
        errors=parameter_errors(c,nodes,ports)
        if max(errors.values())!=0.:raise AssertionError('Canonical parameters are not exactly permutation-symmetric')
        for key,a in states.items():
            if not np.array_equal(getattr(c,key),a) or id(getattr(c,key))!=identities[key]:raise AssertionError(f'Dynamic/initial state changed: {key}')
        for key in VECTOR_PARAMETERS:
            if key in old and not np.array_equal(getattr(c,key)[mid],old[key][mid]):raise AssertionError('Shared intrinsic parameters changed')
        for key in MATRIX_PARAMETERS:
            if not np.array_equal(getattr(c,key)[np.ix_(mid,mid)],old[key][np.ix_(mid,mid)]):raise AssertionError('Shared-to-shared matrix changed')
    except Exception:
        for key,a in old.items():getattr(c,key)[:]=a
        raise
    metadata=[]
    for key in ('syn_families','syn_tau_mult','syn_df_channel','syn_mixed'):
        value=getattr(spinal,key,None)
        if isinstance(value,dict) and 'r' in value and 'l' in value:
            metadata.append(dict(field=key,old_left=copy.deepcopy(value['l']),canonical_right=copy.deepcopy(value['r'])));value['l']=copy.deepcopy(value['r'])
    ledger=dict(schema=chosen['schema'],policy=chosen,engineering_candidate=True,interpretation='The left electrical function is now defined as the physical-name mirror of the actual right function. This is not recovery/alignment of an independently learned left NMF decomposition.',source_before_sha256=before,projected_sha256=field_digest(c),actual_neurons=len(names),bilateral_pairs=len(right),midline_nodes_frozen=[names[i] for i in mid],unmatched_bindings=[],state_and_initialization_arrays_unchanged=list(states),neuron_parameter_changes=changed,parameter_symmetry_errors=errors,edge_changes=edges,edges_added=sum(e['action']=='added' for e in edges),edges_removed=sum(e['action']=='removed' for e in edges),edges_replaced=sum(e['action']=='replaced' for e in edges),runtime_metadata_copied_from_right=metadata,builder_left_NMF_metadata='Original net.syn_W/syn_g_stats and builder graph remain historical source metadata; use final compiled candidate ledger.',input_alias_normalization=ports['aliases'],input_permutation=ports['permutation'].tolist(),neuron_permutation=nodes['permutation'].tolist(),training_allowed=False,whole_body_result=False)
    spinal.canonical_right_ledger=ledger
    return ledger


def verify_canonical(spinal):
    try:from .graded_profile import raw_compiled,population_names
    except ImportError:from graded_profile import raw_compiled,population_names
    c=raw_compiled(spinal);nodes=binding(population_names(spinal));ports=binding(spinal.inputs,input_ports=True);errors=parameter_errors(c,nodes,ports)
    if max(errors.values())!=0:raise ValueError(f'Canonical symmetry invalidated by later transform: {errors}')
    return dict(all_static_parameters_exactly_symmetric=True,maximum_error=0.,fields_checked=list(errors))


def swap_test_state(c,p):
    for key in ('V','V_last','b_gate','b_gate_last','c_gate','c_gate_last'):
        if hasattr(c,key):getattr(c,key)[:]=getattr(c,key)[...,p].copy()


def probe(model,out,repeat=False):
    sys.path[:0]=[str(HERE.parent),str(HERE.parent/'w2l_cpg'),str(HERE.parent.parent/'robot_torso'/'bpa_neural_coordination')]
    os.environ['AARL_NET']='syn6'
    import mujoco,build_network as B,gait_standing_family as F
    from graded_profile import default_profile,apply_profile,verify_profile,field_digest,population_names
    from sparse_pairwise_eval import SparsePairwiseForward
    m=mujoco.MjModel.from_binary_path(str(model));acts=[m.actuator(i).name for i in range(m.nu)];nmf=HERE.parent/'synergy_basis.npz';nmf_hash=hashlib.sha256(nmf.read_bytes()).hexdigest()
    net=B.build(acts);baseline=field_digest(net.compiled);rng=np.random.default_rng(1101);h=hashlib.sha256()
    for step in range(128):
        x=rng.uniform(0.,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;net.compiled.forward(x);h.update(net.compiled.V.tobytes())
    if h.hexdigest()!='0283aae4d70531aad60a8816611af7c00dc76d7419a646cc8947ccddb12456a4':raise AssertionError('Default syn6 neural anchor changed')
    net=B.build(acts);c=net.compiled;rng=np.random.default_rng(1103)
    c.V[:]=rng.uniform(0.,2.,len(c.V));c.V_last[:]=rng.uniform(0.,2.,len(c.V))
    for key in ('b_gate','b_gate_last','c_gate','c_gate_last'):
        if hasattr(c,key):getattr(c,key)[:]=rng.uniform(.2,.8,getattr(c,key).shape)
    untouched=field_digest(c);disabled=apply_canonical_right(net)
    if field_digest(c)!=untouched:raise AssertionError('Default-off policy modified fields')
    ledger=apply_canonical_right(net,default_policy());before_tie=field_digest(c);tie=F.tie_weights(net)
    if tie or field_digest(c)!=before_tie:raise AssertionError('Historical mean projection still changes canonical candidate')
    # Build accelerator after the topology change; apply FSA to its actual caches.
    net.compiled=SparsePairwiseForward(c);fsa=apply_profile(net,default_profile('syn6'),family='syn6');verify_profile(net);verify_canonical(net)
    if F.tie_weights(net):raise AssertionError('FSA profile broke canonical bilateral equality')
    verify_profile(net)
    p=np.array(ledger['neuron_permutation']);q=np.array(ledger['input_permutation']);names=population_names(net)
    alias_ports=[{'BAL_LAT_R':'BILAT_R','BAL_LAT_L':'BILAT_L'}.get(n,n) for n in net.inputs]
    if not np.array_equal(binding(alias_ports,input_ports=True)['permutation'],q):raise AssertionError('BILAT input alias permutation mismatch')
    a=copy.deepcopy(c);b=copy.deepcopy(c);swap_test_state(b,p);rng=np.random.default_rng(1104);digest=hashlib.sha256();errors=[]
    for step in range(128):
        x=rng.uniform(0.,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;a.forward(x);b.forward(x[q]);digest.update(a.V.tobytes())
        error=max(float(np.max(abs(b.V-a.V[p]))),float(np.max(abs(b.b_gate-a.b_gate[...,p]))),float(np.max(abs(b.c_gate-a.c_gate[...,p]))))
        errors.append(error)
        if not np.isfinite(a.V).all() or error>1e-9:raise AssertionError(f'Permutation equivariance failed at {step}: {error}')
    # Reproduce the exact same initial states/inputs through the fast candidate.
    fast=SparsePairwiseForward(copy.deepcopy(c));rng=np.random.default_rng(1104);fast_hash=hashlib.sha256()
    for step in range(128):
        x=rng.uniform(0.,.2,len(net.inputs));x[net.input_index('DRIVE')]=1.2;fast.forward(x);fast_hash.update(fast.V.tobytes())
    if fast_hash.hexdigest()!=digest.hexdigest():raise AssertionError('Canonical profile dense/accelerated trace differs')
    if hashlib.sha256(nmf.read_bytes()).hexdigest()!=nmf_hash:raise AssertionError('Source K6 file changed')
    result=dict(scope='Engineered canonical-right graded electrical candidate; no body integration,Hunt,curriculum,NMF refit,or initial-phase reset',source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),default_neural_anchor_sha256=h.hexdigest(),default_compiled_sha256=baseline,policy=ledger,fsa_consumption=verify_profile(net),fsa_encoder_tau_s=fsa['encoder_tau_s'],historical_tie_changes=0,all_static_parameters=verify_canonical(net),input_BILAT_alias_test=True,one_step_permutation_error=errors[0],maximum128step_permutation_error=max(errors),permutation_tolerance=1e-9,dense_accelerated128_exact=True,candidate_neural_state_sha256=digest.hexdigest(),source_NMF_sha256=nmf_hash,source_NMF_unchanged=True,training_allowed=False,coverage_verified=False)
    out.mkdir(parents=True,exist_ok=True)
    if repeat:
        saved=json.loads((out/'CANONICAL_RIGHT.json').read_text())
        if saved!=result:raise AssertionError('Fresh-process canonical electrical evidence differs')
        (out/'FRESH_REPEAT.json').write_text(json.dumps(dict(exact=True,candidate_neural_state_sha256=digest.hexdigest(),source_code_sha256=result['source_code_sha256']),indent=2))
    else:
        (out/'CANONICAL_RIGHT.json').write_text(json.dumps(result,indent=2));(out/'POLICY.json').write_text(json.dumps(default_policy(),indent=2))
    print(json.dumps(dict(repeat=repeat,default_exact=True,added=ledger['edges_added'],removed=ledger['edges_removed'],replaced=ledger['edges_replaced'],one_step_error=errors[0],maximum128step_error=max(errors),dense_fast_exact=True)))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--model',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--repeat',action='store_true');args=parser.parse_args();probe(args.model,args.out,args.repeat)
