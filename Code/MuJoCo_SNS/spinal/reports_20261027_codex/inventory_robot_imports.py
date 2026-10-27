"""Read-only inventory of candidate lower-body imports; no CAD inference."""
import collections
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
SIMSCAPE=REPO/'Code/Matlab/SNS_Simscape'
results=[]
for name in ('mdl_humanoid_lower_ah001_imported.slx',
             'mdl_humanoid_lower_ah001_imported_R2025a.slx'):
    path=SIMSCAPE/name
    blocks=[]
    with zipfile.ZipFile(path) as archive:
        for member in archive.namelist():
            if member.endswith('.xml') and member.startswith('simulink/'):
                root=ET.fromstring(archive.read(member))
                for block in root.iter('Block'):
                    properties={p.get('Name'):p.text for p in block.findall('P')}
                    sources=[properties.get(k,'') or '' for k in ('SourceBlock','ReferenceBlock')]
                    if any('/Joints/' in s or 'Joint' in s for s in sources):
                        blocks.append(dict(name=block.get('Name'),properties=properties,
                                           archive_member=member))
    counts=collections.Counter((b['properties'].get('SourceBlock') or
                                b['properties'].get('ReferenceBlock')) for b in blocks)
    results.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        joint_reference_counts=dict(counts),joint_blocks=blocks,
                        scope='Static SLX block inventory; not a compiled mobility test'))
(HERE/'robot_import_inventory.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
for result in results:
    print(result['path'],result['joint_reference_counts'])
