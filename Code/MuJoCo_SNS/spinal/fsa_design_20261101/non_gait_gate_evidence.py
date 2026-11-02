"""Design-only bilateral accounting and evidence packaging; never releases training."""
from pathlib import Path
import json, re, hashlib
import numpy as np
import non_gait_contracts as N
import xml.etree.ElementTree as ET

HERE=Path(__file__).parent

def mirrored(name):
    for left,right in [('LH','RH'),('Left','Right'),('left','right')]:
        if left in name:return name.replace(left,right)
        if right in name:return name.replace(right,left)
    match=re.search(r'(^|[_ ])([lrLR])(?=[_ ]|$)',name)
    if match:return name[:match.start(2)]+({'l':'r','r':'l','L':'R','R':'L'}[match[2]])+name[match.end(2):]
    return None

def right_name(name):
    return 'RH' in name or 'Right' in name or 'right' in name or bool(re.search(r'(^|[_ ])[rR](?=[_ ]|$)',name))

def edge_account(names,edges):
    names=list(names)
    canonical=lambda name:re.sub(r'\s+','',name)
    keys=[canonical(name) for name in names]
    index={key:i for i,key in enumerate(keys) if keys.count(key)==1}
    lookup={}
    for i,row in enumerate(edges):lookup.setdefault((int(row[0]),int(row[1])),[]).append(i)
    pairs=[];unmatched=[];ineligible=[];seen=set()
    for i,row in enumerate(edges):
        pre,post=names[int(row[0])],names[int(row[1])]
        if not (right_name(pre) or right_name(post)):continue
        mp,mq=mirrored(pre) or pre,mirrored(post) or post
        candidates=lookup.get((index.get(canonical(mp),-1),index.get(canonical(mq),-1)),[])
        if len(candidates)!=1:unmatched.append(dict(source=pre,target=post,mirror_count=len(candidates)));continue
        j=candidates[0]
        if i in seen or j in seen:continue
        seen.update((i,j))
        if not np.array_equal(row[3:],edges[j,3:]):
            ineligible.append(dict(source=pre,target=post,right_edge=i,left_edge=j,right_fields=row[2:].tolist(),left_fields=edges[j,2:].tolist(),reason='Different source Erev/presynaptic window; no automatic homology assertion'));continue
        pairs.append(dict(right_edge=i,left_edge=j,right_source=pre,right_target=post,right_source_g=float(row[2]),left_source_g=float(edges[j,2]),source_gain_equal=bool(row[2]==edges[j,2])))
    return dict(eligible_pairs=pairs,source_asymmetric_pair_count=sum(not row['source_gain_equal'] for row in pairs),unmatched_right_edges=unmatched,ineligible_window_or_reversal=ineligible,
                source_asymmetry_policy='Baseline untouched. Explicit new candidate uses right canonical gain copied to verified left edge; unmatched/non-equivalent edges frozen, never independently fit.',name_mapping='Whitespace-insensitive, unique source names only; no anatomical/role-based guessed relabeling.',
                membrane_and_sensor_geometry_policy='Source cell/sensor/physical-route differences remain preserved and separately declared; edge eligibility alone is not whole-cell or anatomical equivalence.')

def apply_tied_gains(edges,pairs,canonical_values=None):
    """Explicit opt-in right-canonical application, no default alteration."""
    if canonical_values is None:return
    if len(canonical_values)!=len(pairs):raise ValueError('One value per verified pair')
    values=np.asarray(canonical_values,float)
    if not np.isfinite(values).all() or np.any(values<0):raise ValueError('Finite nonnegative conductances required')
    for pair,value in zip(pairs,values):edges[pair['right_edge'],2]=value;edges[pair['left_edge'],2]=value

def proof(account,edges):
    unchanged=edges.copy();apply_tied_gains(edges,account['eligible_pairs']);assert np.array_equal(edges,unchanged)
    pairs=account['eligible_pairs'];values=[row['right_source_g'] for row in pairs]
    if not pairs:raise AssertionError('Bilateral source has no recognized eligible pairs; inspect naming/eligibility, do not claim vacuous symmetry.')
    apply_tied_gains(edges,pairs,values)
    assert all(edges[p['right_edge'],2]==edges[p['left_edge'],2] for p in pairs)
    changed={p['right_edge'] for p in pairs}|{p['left_edge'] for p in pairs}
    assert all(np.array_equal(row,unchanged[i]) for i,row in enumerate(edges) if i not in changed)
    assert np.array_equal(edges[:,3:],unchanged[:,3:])
    return dict(default_none_exact=True,explicit_pair_values_exact=True,unmatched_edges_exact=True,windows_reversals_exact=True,body_steps=0,fit_steps=0)

def inspect_aliases(c,source):
    """A source-ID-bound logical alias requires its full functional neighbourhood."""
    names=list(c.names);root=ET.parse(source).getroot();nodes=list(next(m for m in root.iter('NeuralModule') if m.findtext('ModuleName')=='IntegrateFireSim').iter('Neuron'));id_names={n.findtext('ID'):n.findtext('Name') for n in root.iter() if n.find('ID') is not None and n.find('Name') is not None}
    found=[]
    for i,name in enumerate(names):
        if 'Knee MN ext RE' not in name:continue
        incoming=[dict(source=names[int(e[0])],source_id=nodes[int(e[0])].findtext('ID'),fields=e[2:].tolist()) for e in c.edges if int(e[1])==i]
        outgoing=[dict(target=names[int(e[1])],target_id=nodes[int(e[1])].findtext('ID'),fields=e[2:].tolist()) for e in c.edges if int(e[0])==i]
        source_ids={v['source_id'] for v in incoming}
        bound=[dict(motor=id_names.get(a.findtext('SourceID')),muscle=id_names.get(a.findtext('TargetID')),physical_target_id=a.findtext('TargetID')) for a in root.iter('Adapter') if a.findtext('Type')=='NodeToPhysical' and a.findtext('SourceID') in source_ids]
        right=any(v['source']=='R Knee MN ext' for v in incoming) and any(v['target']=='R Knee MN ext' for v in outgoing)
        wrong=any(v['source']=='L Knee MN ext' for v in incoming) or any(v['target']=='L Knee MN ext' for v in outgoing)
        found.append(dict(index=i,source_id=nodes[i].findtext('ID'),source_name=name,incoming=incoming,outgoing=outgoing,upstream_motor_physical_bindings=bound,right_function_uniquely_supported=right and not wrong,logical_alias='R Knee MN ext RE' if right and not wrong else name))
        if right and not wrong:names[i]='R Knee MN ext RE'
    return names,found

def file(path,role):
    p=Path(path);return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),role=role)

def main():
    reports={}
    for family in ['W2L108_graded','W2L704_hybrid_spiking','rat134_GA_receptor_variant']:
        c,source,anchor=N.actual_native_controller(family);edges=c.edges.copy()
        names,aliases=inspect_aliases(c,source) if family.startswith('W2L') else (c.names,[])
        account=edge_account(names,edges);account['logical_alias_evidence']=aliases;account['application_proof']=proof(account,edges);account['anchor']=anchor
        account['candidate_not_historical_source']=True
        reports[family]=account
    net=N.cad_network();c=net.compiled;names=[p['name'] for p in net.net.populations]
    post,pre=np.nonzero(c.g_max_non);edges=np.array([[a,b,c.g_max_non[b,a],c.e_lo[b,a],c.e_hi[b,a],c.del_e[b,a]] for b,a in zip(post,pre)])
    account=edge_account(names,edges);account['application_proof']=proof(account,edges)
    account['scope']='Dense compiled edge-array proposal; actual historical compiled baseline unchanged. Existing six pressure-pair weights are already shared bilateral values; hybrid allocator targets depend on actual state, not independent learned L/R parameters.'
    reports['CAD262_hybrid']=account;reports['CAD262_MN_pressure_pair']=account
    for family in ['SNS_Simscape_graded_reflexes','SNS_Simscape_spiking_reflexes','reduced_Deng14']:
        reports[family]=dict(applicable=False,reason='Unilateral antagonist primitive/mechanism; flexor/extensor cells are functionally different, not left/right homologues. No bilateral fitting variables exist.',blanket_bilateral_body_claim=False)
    destination=N.RAW_RECORDS/'NON_GAIT_BILATERAL_DESIGN_v4.json'
    if destination.exists():raise FileExistsError(destination)
    destination.write_text(json.dumps(dict(machine='easteregg2',families=reports,physics_steps=0,training_steps=0,source_sha256=N.sha(__file__)),indent=2))
    for family,account in reports.items():
        ledger_path=N.LEDGERS/(family+'.json');ledger=json.loads(ledger_path.read_text())
        ledger['bilateral_design_accounting']=dict(evidence=N.evidence(destination),applicable=account.get('applicable',True),eligible_edge_pairs=len(account.get('eligible_pairs',[])),source_asymmetric_pair_count=account.get('source_asymmetric_pair_count',0),unmatched_right_edges=len(account.get('unmatched_right_edges',[])),ineligible_edges=len(account.get('ineligible_window_or_reversal',[])),source_untouched=True,blanket_source_symmetry_claim=False,policy=account.get('source_asymmetry_policy',account.get('reason')))
        ledger_path.write_text(json.dumps(ledger,indent=2))
    print(json.dumps({key:dict(eligible_pairs=len(value.get('eligible_pairs',[])),source_asymmetric=value.get('source_asymmetric_pair_count'),unmatched=len(value.get('unmatched_right_edges',[])),ineligible=len(value.get('ineligible_window_or_reversal',[]))) for key,value in reports.items()},indent=2))

if __name__=='__main__':main()
