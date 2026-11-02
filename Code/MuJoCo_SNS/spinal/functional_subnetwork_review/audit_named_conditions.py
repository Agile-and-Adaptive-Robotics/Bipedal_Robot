"""Audit Ben's named anchors without loading neural or physical models."""
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path, PureWindowsPath
import xml.etree.ElementTree as E
import zipfile
from audit_saved_graphs import ROOT, OUT, table, digest

INTAKE = ROOT / 'tmp/functional_subnetwork_remote/named_conditions'

def serialized(x):
    return E.tostring(x, encoding='unicode')

def xml_fields(x):
    row = {}
    def visit(e, path):
        if len(e):
            for child in e:
                visit(child, path + '/' + child.tag if path else child.tag)
        elif e.text:
            row[path] = e.text
        for k, v in e.attrib.items():
            row[path + '/@' + k] = v
    visit(x, '')
    return row

def drawing_checks(root):
    result = []
    for page in root.findall('.//NervousSystem//Node'):
        raw = page.findtext('DiagramXml')
        if not raw:
            continue
        d = E.fromstring(raw)
        flow = d if d.tag == 'AddFlow' else d.find('.//AddFlow')
        entries = [x for x in flow if x.tag in ('Node', 'Link')]
        errors = []
        logical = {x.findtext('ID'): x for x in page.findall('Links/Link')}
        for x in entries:
            if x.tag != 'Link':
                continue
            handles = [int(x.attrib[k]) for k in ('Org', 'Dst')]
            if any(i < 0 or i >= len(entries) or entries[i].tag != 'Node' for i in handles):
                errors.append({'link': x.findtext('Tag'), 'invalid_handles': handles})
                continue
            l = logical.get(x.findtext('Tag'))
            if l is not None and [entries[i].findtext('Tag') for i in handles] != [l.findtext('OriginID'), l.findtext('DestinationID')]:
                errors.append({'link': x.findtext('Tag'), 'logical_endpoint_mismatch': True})
        for tag in ('Node', 'Link'):
            if int(flow.attrib[tag + 's']) != sum(x.tag == tag for x in entries):
                errors.append({'bad_count': tag})
        result.append({'page_id': page.findtext('ID'), 'page': page.findtext('Text'), 'errors': errors})
    return result

def rat_audit():
    folder = OUT / 'rat134'
    folder.mkdir(exist_ok=True)
    path = INTAKE / 'rat134/Rat134_ReceptorGeometry_Standalone.asim'
    expected = 'e3206341b7f89f3cb97584359ca83f136b3c7881b9fc0cf6d6ecaf9853606888'
    assert digest(path) == expected
    root = E.parse(path).getroot()
    nodes, edges, adapters = list(root.iter('Neuron')), list(root.iter('Connexion')), list(root.iter('Adapter'))
    node_map = {x.findtext('ID'): x for x in nodes}
    types = {x.findtext('ID'): x for x in root.iter('SynapseType')}
    bodies = {x.findtext('ID'): x for x in root.iter('RigidBody')}
    assert len(nodes) == 134 and len(edges) == 274
    rows = []
    for x in edges:
        r = xml_fields(x)
        a, b = x.findtext('SourceID'), x.findtext('TargetID')
        assert a in node_map and b in node_map
        r.update(source_name=node_map[a].findtext('Name'), target_name=node_map[b].findtext('Name'))
        t = types[x.findtext('SynapseTypeID')]
        r.update({'type/' + k: v for k, v in xml_fields(t).items()})
        r['conductance_policy'] = 'SynAmp effective; connexion G retained but not used as effective conductance'
        r['current_sign_rule'] = 'g*s*(Equil-Vpost); label alone does not fix sign'
        rows.append(r)
    table(folder / 'nodes.csv', [xml_fields(x) for x in nodes])
    table(folder / 'edges.csv', rows)
    table(folder / 'types.csv', [xml_fields(x) for x in types.values()])
    inputs = []
    for x in adapters:
        r = xml_fields(x)
        for role in ('Source', 'Target'):
            i = x.findtext(role + 'ID')
            obj = node_map.get(i)
            if obj is None:
                obj = bodies.get(i)
            r[role + 'Name'] = obj.findtext('Name') if obj is not None else 'unresolved'
        inputs.append(r)
    table(folder / 'adapters.csv', inputs)
    table(folder / 'physical_receptor_fields.csv', [xml_fields(x) for x in bodies.values() if 'Receptor' in (x.findtext('Type') or '')])
    comparator = json.loads((INTAKE / 'rat134_comparator/CONSUMED_FIELD_LEDGER.json').read_text())
    table(folder / 'refinement12_comparator_only.csv', comparator)
    project = E.parse(INTAKE / 'rat134/Rat134_ReceptorGeometry.aproj').getroot()
    ids = [x.findtext('ID') for x in project.iter() if x.tag not in ('InLinks', 'OutLinks') and x.find('ID') is not None and len(x)>1]
    # Preserve original/native graph and separate runtime overlay, not a refitted condition.
    runtime = INTAKE / 'rat134_scripts/rat134_receptor_runtime.py'
    ast.parse(runtime.read_text())
    result = {'standalone_sha256': expected, 'neurons': len(nodes), 'edges': len(edges), 'adapters': len(adapters), 'enabled_physical_inputs': sum(x.findtext('Type')=='PhysicalToNode' and x.findtext('Enabled')!='False' for x in adapters), 'duplicate_project_definition_IDs': [k for k,v in Counter(ids).items() if v>1], 'saved_pages': drawing_checks(project), 'runtime_entry': 'variant_controller()', 'runtime_source_sha256': digest(runtime), 'runtime_overlay': 'Replaces source Hill bank/parameters and physical input polynomial/limits from curated asim; do not use rat134_core.controller() alone.', 'reference_condition': json.loads((INTAKE/'rat134/CONDITION_DECISIONS.json').read_text()), 'comparator': 'refinement12 nominal trial3 only; do not install gains; failed cadence and amplitude holdouts', 'execution': 'none; static audit only'}
    (folder / 'STATIC_VALIDATION.json').write_text(json.dumps(result,indent=2))
    return result

def main():
    retrieval=json.loads((INTAKE/'RETRIEVAL_MANIFEST.json').read_text())
    mapping={r['path'].lower(): INTAKE/r['archive_path'] for r in retrieval if 'archive_path' in r}
    checks=[{'path':r['path'],'matches':digest(INTAKE/r['archive_path'])==r['sha256']} for r in retrieval if 'archive_path' in r]
    assert all(r['matches'] for r in checks)
    extra=ROOT/'tmp/functional_subnetwork_remote/dependency_evidence'
    extra_manifest=json.loads((extra/'MANIFEST.json').read_text())
    for r in extra_manifest:
        if 'archive_path' in r:
            p=extra/r['archive_path']
            assert digest(p)==r['sha256']
            mapping[r['path'].lower()]=p
    source_parity=[]
    for p in (INTAKE/'non_gait_ledgers').glob('*.json'):
        evidence=json.loads(p.read_text()).get('source_evidence',[])
        for r in evidence:
            saved=mapping.get(r['path'].lower())
            record={'ledger':p.name,'source':r['path'],'expected_sha256':r.get('sha256'),'retrieved_status':'not_retrieved' if saved is None else 'match' if digest(saved)==r.get('sha256') else 'mismatch'}
            if '/bipedal_robot/' in r['path'].replace('\\','/').lower():
                rel=r['path'].replace('\\','/').split('Bipedal_Robot/',1)[-1]
                live=ROOT/rel
                record['EB475WS4_checkout_status']='missing' if not live.exists() else 'match' if digest(live)==r.get('sha256') else 'mismatch; isolated evidence retained'
            source_parity.append(record)
    table(OUT/'NAMED_SOURCE_PARITY.csv',source_parity)
    rat_manifest=[]
    for r in json.loads((INTAKE/'rat134/MANIFEST.json').read_text())['files']:
        p=INTAKE/'rat134'/r['path']
        if not p.exists():p=INTAKE/'rat134_scripts'/r['path'].removeprefix('scripts/').removeprefix('scripts\\')
        if not p.exists():
            remote='D:\\Users\\Ben Bolen\\NeuroCampaign_20261031\\subnetworks\\Rat132_bounded_diagnostic\\Curated134_ReceptorGeometry\\'+r['path'].replace('/','\\')
            p=mapping.get(remote.lower(),p)
        rat_manifest.append({'path':r['path'],'status':'not_retrieved' if not p.exists() else 'match' if digest(p)==r['sha256'] else 'mismatch'})
    snapshots=OUT/'Gait_TimingZeroSFA/source_snapshot'
    snapshots.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(INTAKE/'gait_v4/SOURCE_SNAPSHOT.zip') as z:
        for name in z.namelist():
            assert Path(name).name==name and name not in ('..','.')
            (snapshots/name).write_bytes(z.read(name))
    candidate=json.loads((INTAKE/'gait_v4/TIMING_CANDIDATE.json').read_text())
    assert candidate['controller_class']=='TimingZeroSFA'
    assert candidate['parameters']=={'rg_adapt_jump_nA':0.,'rg_readout_C_rate':5.,'pf_rate_scale':1.,'rg_h_tau_scale':0.9496883348174533}
    assert (snapshots/'TIMING_CANDIDATE.json').read_bytes()==(INTAKE/'gait_v4/TIMING_CANDIDATE.json').read_bytes()
    syntax=[]
    for p in snapshots.glob('*.py'):
        ast.parse(p.read_text(encoding='utf-8-sig'))
        syntax.append({'file':p.name,'sha256':digest(p),'syntax':'pass'})
    dependencies=[]
    for r in candidate['dependencies']:
        p=snapshots/PureWindowsPath(r['path']).name
        dependencies.append({'file':p.name,'matches':p.exists() and digest(p)==r['sha256']})
    manifest_checks=[]
    for r in json.loads((INTAKE/'gait_v4/FILE_MANIFEST.json').read_text())['files']:
        p=mapping.get(r['path'].lower())
        manifest_checks.append({'path':r['path'],'status':'not_retrieved' if p is None else 'match' if digest(p)==r['sha256'] else 'mismatch'})
    rat=rat_audit()
    result={'retrieval_checks':checks,'extra_dependency_files_verified':len([r for r in extra_manifest if 'archive_path' in r]),'non_gait_source_parity':source_parity,'rat_manifest_checks':rat_manifest,'gait_snapshot_syntax':syntax,'gait_frozen_dependency_checks':dependencies,'gait_file_manifest_checks':manifest_checks,'rat134':rat,'training_allowed':False,'simulations_executed':0}
    (OUT/'NAMED_CONDITION_VALIDATION.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'retrieved_files_verified':len(checks),'gait_dependencies_match':all(x['matches'] for x in dependencies),'gait_manifest':dict(Counter(x['status'] for x in manifest_checks)),'rat_neurons':rat['neurons'],'rat_edges':rat['edges'],'rat_pages':len(rat['saved_pages']),'rat_page_issues':sum(len(x['errors']) for x in rat['saved_pages'])}))

if __name__=='__main__':
    main()
