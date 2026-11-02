"""NEW122/136 rat candidates: LEFT biarticular module mirrored, no body/training."""
import copy,json,re,uuid,collections,sys
import numpy as np
import non_gait_contracts as N
from rat_bilateral_design import side_key

def install(c,profile=None):
    if profile is None:return c
    if len(c.names)!=profile['source_neurons']:raise ValueError('Construction-time source census mismatch')
    n=len(c.names)
    for key,value in list(c.__dict__.items()):
        if isinstance(value,np.ndarray) and value.ndim and value.shape[0]==n and key not in ('edges','source_edges'):
            c.__dict__[key]=np.concatenate([value,value[profile['clone_left_indices']].copy()],axis=0)
    c.names=list(c.names)+profile['new_names']
    if hasattr(c,'ids'):
        c.ids=dict(c.ids)
        for i,ident in enumerate(profile['new_ids']):c.ids[ident]=n+i
    c.edges=np.asarray(profile['edges'],float).copy()
    if hasattr(c,'edge_ids'):c.edge_ids=list(profile['edge_ids'])
    c.bilateral_biarticular_design_identity=profile['identity']
    return c

def build(family):
    c,source,anchor=N.actual_native_controller(family);n=len(c.names);logical=list(c.names)
    old={k:v.copy() for k,v in c.__dict__.items() if isinstance(v,np.ndarray)};old_names=list(c.names)
    physical_io=repr((c.inputs,c.outputs));groups=collections.defaultdict(dict)
    # LF is a source typo forLEFT, proven by all actual neighbours; retain native name.
    aliases=[]
    for i,name in enumerate(logical):
        if re.sub(r'\s+','',name).upper()=='LF_RFINEXT':
            neighbours=[int(e[0]) for e in c.edges if int(e[1])==i]+[int(e[1]) for e in c.edges if int(e[0])==i]
            assert neighbours and all(side_key(c.names[k])[0]=='LH' for k in neighbours)
            logical[i]='LH_RF IN EXT';aliases.append(dict(index=i,native_name=name,logical_alias=logical[i],neighbour_names=[c.names[k] for k in neighbours]))
        side,key=side_key(logical[i]);assert side is not None,repr(logical[i])
        if side in groups[key]:raise ValueError('Ambiguous role:'+key)
        groups[key][side]=i
    assert {key for key,row in groups.items() if set(row)!=set(['RH','LH'])}=={'bfpin','rfin'}
    clones=[groups[key]['LH'] for key in ('bfpin','rfin')]
    assert all(c.params[i,3]==0 for i in clones)
    new_names=['RH_BFP IN EXT','RH_RF IN EXT'];new_ids=[str(uuid.uuid5(uuid.NAMESPACE_URL,N.sha(source)+'/bilateral-biarticular-integration/'+name)) for name in new_names]
    for j,key in enumerate(('bfpin','rfin')):groups[key]['RH']=n+j
    logical+=new_names
    mirror={}
    for row in groups.values():r,l=row['RH'],row['LH'];mirror[r]=l;mirror[l]=r
    fullparams=np.concatenate([c.params,c.params[clones]],axis=0)
    source_edges=c.edges.copy();desired={};removed_direct=[]
    def add(edge):
        key=(int(edge[0]),int(edge[1]))
        if key in desired and not np.array_equal(desired[key],edge):raise ValueError('Conflicting edge construction')
        desired[key]=edge.copy()
    for e in source_edges:
        pre,post=int(e[0]),int(e[1]);ps,pk=side_key(logical[pre]);qs,qk=side_key(logical[post])
        if ps!='RH':continue
        bi_target=bool(re.match(r'^(RH|LH)[_\s-]*(BFP|RF)(?=[_\s-])',logical[post],re.I))
        if 'PF' in logical[pre].upper() and bi_target and ('MN' in logical[post].upper() or 'IA' in logical[post].upper()):
            removed_direct.append(dict(source=logical[pre],target=logical[post],fields=e[2:].tolist()));continue
        add(e);z=e.copy();z[0],z[1]=mirror[pre],mirror[post];add(z)
    module_edges=[]
    for e in source_edges:
        if int(e[0]) not in clones and int(e[1]) not in clones:continue
        add(e);z=e.copy();z[0],z[1]=mirror[int(e[0])],mirror[int(e[1])];add(z)
        module_edges.append(dict(source=c.names[int(e[0])],target=c.names[int(e[1])],fields=e[2:].tolist()))
    edges=np.array(list(desired.values()));rc_groups=collections.defaultdict(list)
    for i,e in enumerate(edges):
        pre=int(e[0]);post=int(e[1]);_,key=side_key(logical[pre])
        if 're' in key and re.search(r'\bRE\b',logical[pre],re.I):rc_groups[post].append(i)
    rc_contracts=[]
    for post,indices in rc_groups.items():
        op=fullparams[post,0]+2.5;total=sum(edges[i,2] for i in indices)
        assert total>0
        # Declared engineering inhibitory drive20mV below operating point,
        # not an inferred physiological chloride reversal or selected side.
        for i in indices:
            edges[i,2]=.5/(fullparams[post,2]*20.)*edges[i,2]/total;edges[i,5]=op-20.
        current=float(sum(edges[i,2]*(edges[i,5]-op)*fullparams[post,2] for i in indices));assert abs(current+.5)<1e-12
        rc_contracts.append(dict(post=logical[post],operating_mV=float(op),desired_fullwindow_current_nA=-.5,reversal_mV=float(op-20),measured_current_nA=current,indices=indices))
    lookup={(int(e[0]),int(e[1])):e for e in edges}
    assert all(np.array_equal(e[2:],lookup[(mirror[int(e[0])],mirror[int(e[1])])][2:]) for e in edges)
    profile=dict(identity=family+'_bilateral_biarticular_integration_'+str(n+2),source_neurons=n,source_sha256=N.sha(source),clone_left_indices=clones,new_names=new_names,new_ids=new_ids,edges=edges.tolist(),edge_ids=[str(uuid.uuid5(uuid.NAMESPACE_URL,N.sha(source)+'/biart/'+str(int(e[0]))+'/'+str(int(e[1])))) for e in edges])
    install(c,None);assert old_names==c.names and all(np.array_equal(v,getattr(c,k)) for k,v in old.items())
    install(c,profile);assert len(c.names)==n+2 and physical_io==repr((c.inputs,c.outputs))
    for key,value in old.items():
        if key=='edges':continue
        now=getattr(c,key);assert np.array_equal(value,now[:n] if value.ndim and value.shape[0]==n and key!='source_edges' else now),key
    import native_kernel as K
    max_operator_error=0.;max_rc_error=0.
    for level in [0.,.5,1.]:
        v=np.full(n+2,-60.+20*level);before=v.copy();m=np.zeros(n+2);h=np.ones(n+2);newv=v.copy();gs=v.copy();gi=v.copy();u=np.zeros(n+2)
        K.neural_steps(v,m,h,newv,gs,gi,c.params,c.edges,c.tonic,u,c.dt*1000,1)
        syn=c.params[:,2]*(gi-gs*(before-c.params[:,0]))
        max_operator_error=max(max_operator_error,max(abs(syn[r]-syn[l]) for r,l in mirror.items()))
    for contract in rc_contracts:
        post=logical.index(contract['post']);subset=c.edges[contract['indices']];v=c.params[:,0].copy();v[subset[:,0].astype(int)]=subset[:,4];v[post]=contract['operating_mV'];before=v.copy();m=np.zeros(n+2);h=np.ones(n+2);newv=v.copy();gs=v.copy();gi=v.copy();u=np.zeros(n+2)
        K.neural_steps(v,m,h,newv,gs,gi,c.params,subset,c.tonic,u,c.dt*1000,1)
        measured=c.params[post,2]*(gi[post]-gs[post]*(before[post]-c.params[post,0]));max_rc_error=max(max_rc_error,abs(measured+.5))
    assert max_operator_error<1e-10 and max_rc_error<1e-10
    # The actual adapter dispatch still reaches the same original MN/physical ports.
    adapter=N.native_adapter_consumption_checks(c);assert adapter['actual_cached_mapping_vs_source_tuple_error_nA']<1e-10
    passive,sensors=N.native_input_contract(c,source);motor_rows,motors=N.native_motor_contract(c,source);passive['neurons'].extend(motor_rows)
    consumed=N.apply_native_passive_profile(c,passive)
    conditioned={c.names.index(row['name']) for row in sensors}
    conditional=N.native_conditional_input_checks(c,sorted(set(range(n+2))-conditioned))
    ionic=N.native_ionic_mechanism_checks(c)
    assert all(row['maximum_actual_native_step_error_mV']<1e-10 for row in conditional)
    assert all(row['compiled_reference_error']<1e-10 and row['state_bounds_preserved'] and row['positive_physical_gate_times'] for row in ionic)
    return dict(passed=True,source=N.evidence(source),anchor=anchor,source_neurons=n,candidate_neurons=n+2,source_edges=len(source_edges),candidate_edges=len(edges),new_passive_cells=[dict(name=name,source_left_name=old_names[i],source_left_parameters=old['params'][i].tolist()) for name,i in zip(new_names,clones)],aliases=aliases,left_module_edges=module_edges,removed_right_direct_biPF=removed_direct,rc_contracts=rc_contracts,profile=profile,passive_role_profile=passive,passive_role_consumption=consumed,independent_sensor_contracts=sensors,motor_domains=motors,conditional_cell_checks=conditional,retained_ionic_checks=ionic,proof=dict(default_none_exact=True,all_old_array_states_and_parameters_exact_before_explicit_passive_profile=True,all_old_physical_adapter_and_MN_bindings_exact=True,complete_first_dimension_buffers_appended=True,actual_JIT_bilateral_synaptic_operator_error_nA=max_operator_error,actual_JIT_RC_budget_error_nA=max_rc_error,actual_adapter_consumption=adapter,all_edges_fullfield_mirrored=True,physics_steps=0,training_steps=0),old_rightcopy_candidate_status='Unadopted: rejected because it removed the requested LEFT biarticular integration architecture.')

if __name__=='__main__':
    suffix=sys.argv[1] if len(sys.argv)>1 else '_v1'
    assert re.fullmatch(r'_[a-zA-Z0-9_]+',suffix)
    families=[sys.argv[2]] if len(sys.argv)>2 else ['public_rat120_precursor','rat134_GA_receptor_variant']
    assert all(family in ['public_rat120_precursor','rat134_GA_receptor_variant'] for family in families)
    rows={family:build(family) for family in families}
    path=N.RAW_RECORDS/('BILATERAL_BIARTICULAR_INTEGRATION'+suffix+'.json');assert not path.exists();path.write_text(json.dumps(dict(machine='easteregg2',families=rows,source_sha256=N.sha(__file__),physics_steps=0,training_steps=0),indent=2))
    print(json.dumps({f:dict(passed=x['passed'],candidate_neurons=x['candidate_neurons'],candidate_edges=x['candidate_edges']) for f,x in rows.items()}))
