"""Explicit NEW right-canonical rat design; original XML and dynamic states untouched."""
import json,re,collections,xml.etree.ElementTree as E
import numpy as np
import non_gait_contracts as N

def side_key(name):
    match=re.match(r'^(RH|LH)[_\s-]*(.*)',name,re.I)
    if not match:return None,None
    tail=match[2];muscle=re.match(r'^(BFA|BFP|IP|VA|RF|SO|TA|GA)(?=[_\s-])',tail,re.I)
    # EXT/FLX labels differ for the same biarticular muscle; retain joint-PF phase labels.
    if muscle:tail=re.sub(r'EXT|FLX','',tail,flags=re.I)
    return match[1].upper(),re.sub(r'[\s_]','',tail).lower()

def construct(family):
    c,source,anchor=N.actual_native_controller(family);root=E.parse(source).getroot()
    module=next(m for m in root.iter('NeuralModule') if m.findtext('ModuleName')=='IntegrateFireSim');nodes=list(module.iter('Neuron'))
    ids=[n.findtext('ID') for n in nodes];bodies={b.findtext('ID'):b.findtext('Name') for b in root.iter() if b.find('ID') is not None and b.find('Name') is not None}
    groups=collections.defaultdict(dict);logical_names=list(c.names);aliases=[]
    for i,name in enumerate(logical_names):
        side,key=side_key(name)
        if side:
            if side in groups[key]:raise ValueError('Ambiguous actual source role:'+key)
            groups[key][side]=i
    preliminary={}
    for row in groups.values():
        if set(row)==set(['RH','LH']):preliminary[row['RH']]=row['LH'];preliminary[row['LH']]=row['RH']
    unmatched=[(key,row) for key,row in groups.items() if set(row)!=set(['RH','LH'])]
    neutral=[i for i,name in enumerate(logical_names) if side_key(name)[0] is None]
    for i in list(neutral):
        incoming=[int(e[0]) for e in c.edges if int(e[1])==i];outgoing=[int(e[1]) for e in c.edges if int(e[0])==i]
        if c.names[i].strip()=='LF_RF IN EXT' and incoming and outgoing and all(side_key(c.names[k])[0]=='LH' for k in incoming+outgoing):
            logical_names[i]='LH_RF IN EXT';groups['rfin']['LH']=i;neutral.remove(i)
            aliases.append(dict(source_id=ids[i],source_name=c.names[i],logical_alias=logical_names[i],proof='Complete incoming/outgoing neighbourhood isLEFT; alltargets are left RF/BFP muscle paths. Not a right homologous interneuron.',incoming=[dict(id=ids[k],name=c.names[k]) for k in incoming],outgoing=[dict(id=ids[k],name=c.names[k]) for k in outgoing]))
    unmatched=[(key,row) for key,row in groups.items() if set(row)!=set(['RH','LH'])]
    if len(unmatched)==1 and len(neutral)==1:
        key,row=unmatched[0]
        if set(row)=={'LH'}:
            left=row['LH'];right=neutral[0]
            incoming=lambda node:{int(e[0]) for e in c.edges if int(e[1])==node}
            outgoing=lambda node:{int(e[1]) for e in c.edges if int(e[0])==node}
            expected_in={preliminary.get(i,-1) for i in incoming(left)};expected_out={preliminary.get(i,-1) for i in outgoing(left)}
            if expected_in==incoming(right) and expected_out==outgoing(right) and expected_in and expected_out and -1 not in expected_in|expected_out:
                logical_names[right]='RH_BFP IN';groups[key]['RH']=right
                aliases.append(dict(source_id=ids[right],source_name=c.names[right],logical_alias=logical_names[right],left_homologue_id=ids[left],left_source_name=c.names[left],proof='Complete incoming AND outgoing neuron-ID neighbour sets exactly mirror the uniquely identified left BFP inhibitory interneuron through59/66 already paired roles; not inferred from the erroneous name.',right_incoming=[dict(id=ids[i],name=c.names[i]) for i in sorted(incoming(right))],right_outgoing=[dict(id=ids[i],name=c.names[i]) for i in sorted(outgoing(right))]))
    left_only={key:row['LH'] for key,row in groups.items() if set(row)=={'LH'} and key in ('bfpin','rfin')}
    unresolved={key:value for key,value in groups.items() if set(value)!=set(['RH','LH']) and key not in left_only}
    mappings={};pairs=[]
    for key,row in groups.items():
        if len(row)!=2:continue
        r,l=row['RH'],row['LH'];mappings[r]=l;mappings[l]=r
        physical=[]
        for a in root.iter('Adapter'):
            if a.findtext('SourceID') in (ids[r],ids[l]) or a.findtext('TargetID') in (ids[r],ids[l]):physical.append(dict(id=a.findtext('ID'),type=a.findtext('Type'),source_id=a.findtext('SourceID'),target_id=a.findtext('TargetID'),source_name=bodies.get(a.findtext('SourceID')),target_name=bodies.get(a.findtext('TargetID')),source_data=a.findtext('SourceDataType'),target_data=a.findtext('TargetDataType')))
        pairs.append(dict(role=key,right_index=r,left_index=l,right_id=ids[r],left_id=ids[l],right_name=c.names[r],left_name=c.names[l],physical_bindings=physical,passive_or_ionic_class_equal=bool((c.params[r,3]>0)==(c.params[l,3]>0))))
    if unresolved:return dict(passed=False,family=family,unresolved_roles=unresolved,source=N.evidence(source),pairs=pairs)
    assert all(row['passive_or_ionic_class_equal'] for row in pairs)
    legacy=[]
    for key,i in left_only.items():
        adapters=[a for a in root.iter('Adapter') if a.findtext('SourceID')==ids[i] or a.findtext('TargetID')==ids[i]]
        assert not adapters and c.params[i,3]==0
        legacy.append(dict(source_id=ids[i],source_name=c.names[i],logical_alias=logical_names[i],role=key,reason='LEFT-only biarticular integrator has no genuine right counterpart. Newright-canonical design removes allincidentedges; row, intrinsicparameters, anddynamicstates retained inactive. Originalsource unchanged.',incoming=[dict(id=ids[int(e[0])],name=c.names[int(e[0])],fields=e[2:].tolist()) for e in c.edges if int(e[1])==i],outgoing=[dict(id=ids[int(e[1])],name=c.names[int(e[1])],fields=e[2:].tolist()) for e in c.edges if int(e[0])==i]))
    original=c.edges.copy();new=[];source_changes=[]
    for i,e in enumerate(original):
        pre,post=int(e[0]),int(e[1]);ps,_=side_key(logical_names[pre]);qs,_=side_key(logical_names[post])
        canonical=ps=='RH' or (ps is None and qs=='RH')
        if canonical:
            new.append(e.copy());mirrored=e.copy();mirrored[0]=mappings.get(pre,pre);mirrored[1]=mappings.get(post,post);new.append(mirrored)
        elif ps is None and qs is None:new.append(e.copy())
    candidate=np.asarray(new);assert len({(int(e[0]),int(e[1])) for e in candidate})==len(candidate)
    assert not any(int(e[0]) in left_only.values() or int(e[1]) in left_only.values() for e in candidate)
    lookup={(int(e[0]),int(e[1])):e for e in candidate}
    assert all(np.array_equal(e[2:],lookup[(mappings.get(int(e[0]),int(e[0])),mappings.get(int(e[1]),int(e[1])))][2:]) for e in candidate)
    old={(int(e[0]),int(e[1])):e for e in original};changes=[]
    for key,row in lookup.items():
        if key not in old or not np.array_equal(old[key][2:],row[2:]):changes.append(dict(source=c.names[key[0]],target=c.names[key[1]],old_fields=old[key][2:].tolist() if key in old else None,new_fields=row[2:].tolist()))
    removed=[dict(source=c.names[p],target=c.names[q],old_fields=e[2:].tolist()) for (p,q),e in old.items() if (p,q) not in lookup]
    before={k:v.copy() for k,v in c.__dict__.items() if isinstance(v,np.ndarray)}
    apply(c,None);assert all(np.array_equal(v,getattr(c,k)) for k,v in before.items())
    apply(c,candidate);assert np.array_equal(c.edges,candidate)
    assert all(np.array_equal(v,getattr(c,k)) for k,v in before.items() if k!='edges')
    # Actual source JIT consumes the replaced edge array; neural electrical probe only.
    import native_kernel as K
    v=c.params[:,0].copy();m=np.zeros(len(v));h=np.ones(len(v));newv=v.copy();gs=v.copy();gi=v.copy();u=np.zeros(len(v))
    K.neural_steps(v,m,h,newv,gs,gi,c.params,c.edges,c.tonic,u,c.dt*1000,1)
    assert np.isfinite(v).all()
    return dict(passed=True,family=family,source=N.evidence(source),anchor=anchor,node_pairs=pairs,logical_aliases=aliases,inactive_left_only_legacy_cells=legacy,source_edge_count=len(original),candidate_edge_count=len(candidate),changed_or_added_edges=changes,removed_left_only_edges=removed,profile=dict(adapter='native_right_canonical_full_fields',edges=candidate.tolist(),source_sha256=N.sha(source)),proof=dict(default_none_exact=True,nonedge_dynamic_states_unchanged=True,all_candidate_edges_fullfield_mirrored=True,actual_native_JIT_consumed=True,physics_steps=0,training_steps=0),unresolved_roles={},scope='NEW engineering candidate; source anatomical muscle bindings preserved. Source EXT/FLX labels on muscle-scoped cells are not assumed different anatomical muscles; joint RG/PF phase labels are retained.')

def apply(controller,edges=None):
    if edges is None:return
    edges=np.asarray(edges,float)
    if edges.ndim!=2 or edges.shape[1]!=6 or not np.isfinite(edges).all():raise ValueError('Valid native6column edge profile required')
    controller.edges=edges.copy()

if __name__=='__main__':
    rows={family:construct(family) for family in ['public_rat120_precursor','rat134_GA_receptor_variant']}
    path=N.RAW_RECORDS/'RAT_RIGHT_CANONICAL_FULL_FIELD_v3.json';assert not path.exists();path.write_text(json.dumps(dict(machine='easteregg2',families=rows,source_sha256=N.sha(__file__),physics_steps=0,training_steps=0),indent=2))
    print(json.dumps({f:dict(passed=x['passed'],unresolved=x['unresolved_roles'],source_edges=x.get('source_edge_count'),candidate_edges=x.get('candidate_edge_count')) for f,x in rows.items()},indent=2))
