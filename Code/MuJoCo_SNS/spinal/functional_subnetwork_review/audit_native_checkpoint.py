"""Read-only XML and drawing audit of the retrieved W2L108 checkpoint."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import subprocess
import xml.etree.ElementTree as E
from audit_saved_graphs import ROOT, OUT, table, digest

SRC = ROOT / 'tmp/functional_subnetwork_remote/W2L108_verified'
DEST = OUT / 'W2L108'
DEST.mkdir(exist_ok=True)
manifest = json.loads((SRC / 'MANIFEST.json').read_text())
hash_checks = [{'path': k, 'matches': (SRC / k).exists() and digest(SRC / k) == v} for k, v in manifest['files_sha256'].items()]
project = SRC / 'editable/W2L_VirtualWalker_Newer.aproj'
root = E.parse(project).getroot()
# InLinks/OutLinks contain ID references, not object definitions.
definition_ids = [x.findtext('ID') for x in root.iter() if x.tag not in ('InLinks', 'OutLinks') and x.find('ID') is not None and len(x) > 1]
duplicates = [k for k, v in Counter(definition_ids).items() if v > 1]
pages = []
for page in root.findall('.//NervousSystem//Node'):
    text = page.findtext('DiagramXml')
    if not text:
        continue
    drawing = E.fromstring(text)
    flow = drawing if drawing.tag == 'AddFlow' else drawing.find('.//AddFlow')
    entries = [x for x in flow if x.tag in ('Node', 'Link')]
    issues = []
    logical = {x.findtext('ID'): x for x in page.findall('Links/Link')}
    for x in entries:
        if x.tag != 'Link':
            continue
        endpoints = []
        for attr in ('Org', 'Dst'):
            i = int(x.attrib[attr])
            if i < 0 or i >= len(entries) or entries[i].tag != 'Node':
                issues.append('invalid handle ' + str(x.attrib))
            else:
                endpoints.append(entries[i].findtext('Tag'))
        link = logical.get(x.findtext('Tag'))
        if link is not None and endpoints != [link.findtext('OriginID'), link.findtext('DestinationID')]:
            issues.append('logical endpoint mismatch ' + str(x.findtext('Tag')))
    counts = {tag: sum(x.tag == tag for x in entries) for tag in ('Node', 'Link')}
    for tag in counts:
        if int(flow.attrib[tag + 's']) != counts[tag]:
            issues.append('header count ' + tag)
    pages.append({'page': page.findtext('Text'), **counts, 'issues': issues})

native = E.parse(SRC / 'editable/W2L_VirtualWalker_Newer_Standalone.asim').getroot()
def payload(x):
    return E.tostring(x, encoding='unicode')
neurons = native.findall('.//Neuron')
edges = native.findall('.//Connexion')
adapters = native.findall('.//Adapter')
types = [x for x in native.iter() if x.find('SynAmp') is not None]
table(DEST / 'native_nodes.csv', [{'id': x.findtext('ID'), 'name': x.findtext('Name'), 'all_native_fields_xml': payload(x)} for x in neurons])
by_id = {x.findtext('ID'): x for x in neurons}
type_by_id = {x.findtext('ID'): x for x in types}
edge_rows = []
for x in edges:
    t = type_by_id.get(x.findtext('SynapseTypeID'))
    row = {'id': x.findtext('ID'), 'source_id': x.findtext('SourceID'), 'target_id': x.findtext('TargetID'), 'source_name': by_id[x.findtext('SourceID')].findtext('Name'), 'target_name': by_id[x.findtext('TargetID')].findtext('Name'), 'all_native_fields_xml': payload(x)}
    if t is not None:
        row['effective_type_xml'] = payload(t)
    edge_rows.append(row)
table(DEST / 'native_edges.csv', edge_rows)
table(DEST / 'native_adapters.csv', [{'id': x.findtext('ID'), 'all_native_fields_xml': payload(x)} for x in adapters])
table(DEST / 'native_types.csv', [{'id': x.findtext('ID'), 'all_native_fields_xml': payload(x)} for x in types])
git = Path('C:/Users/Ben Bolen/AppData/Local/GitHubDesktop/app-3.6.7/resources/app/git/cmd/git.exe')
reference = subprocess.check_output([str(git), 'show', '23b0e785:Neuromechanical_Models/Biped_2xCPG_wSubs/Biped_2xCPG_wSubs.aproj'], cwd=ROOT)
(DEST / 'historical_reference.aproj').write_bytes(reference)
result = {'source': str(project), 'source_sha256': digest(project), 'retrieval_manifest': hash_checks, 'duplicate_project_IDs': duplicates, 'pages': pages, 'native_neurons': len(neurons), 'native_edges': len(edges), 'native_adapters': len(adapters), 'historical_reference_revision': '23b0e785', 'historical_reference_sha256': hashlib.sha256(reference).hexdigest(), 'native_load_and_simulation': 'not run; training/physics gate preserved', 'parameter_parity': 'historical evidence retained in source verification JSON; fresh full semantic comparison pending'}
(DEST / 'STATIC_NATIVE_VALIDATION.json').write_text(json.dumps(result, indent=2))
print(json.dumps({'pages': len(pages), 'page_issues': sum(len(p['issues']) for p in pages), 'duplicate_IDs': len(duplicates), 'manifest_mismatches': sum(not x['matches'] for x in hash_checks), 'neurons': len(neurons), 'edges': len(edges)}))
