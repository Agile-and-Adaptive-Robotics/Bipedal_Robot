"""Materialize one author-source rat project folder; no circuit edits."""
from pathlib import Path
import shutil,json,hashlib,xml.etree.ElementTree as E
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
SRC=HERE/'runs/deng_source_audit_20261031';DEST=REPO/'Neuromechanical_Models/Deng_Biarticular_Rat'
DEST.mkdir(exist_ok=True)
for name in ('Biarticular_net.aproj','Biarticular_net_Standalone.asim'):
 shutil.copy2(SRC/name,DEST/name)
shutil.copy2(SRC/'Biarticular_net_Standalone.asim',DEST/'Biarticular_net.asim')
(DEST/'.gitignore').write_text('runs/\nbuild/\n*.log\n*.npz\n*.txt\n__pycache__/\n')
prov=json.loads((SRC/'PUBLIC_SOURCE_PROVENANCE.json').read_text())
(DEST/'SOURCE_PROVENANCE.json').write_text(json.dumps(dict(status='Public 2020 precursor; not verified final-2023 release',source_repository='https://github.com/vipamol/Rat-SNS-design-and-tuning',sources=prov,project_path=str(DEST/'Biarticular_net.aproj'),standalone_alias='Biarticular_net.asim is a byte-identical alias of the author standalone export.'),indent=2))
ap=E.parse(DEST/'Biarticular_net.aproj').getroot();sim=E.parse(DEST/'Biarticular_net_Standalone.asim').getroot()
def name(e):return e.findtext('Name') or e.findtext('Text')
apb={e.findtext('ID'):e for e in ap.iter('RigidBody')}
rows=[]
for b in sim.iter('RigidBody'):
 a=apb.get(b.findtext('ID'))
 if a is None:continue
 vals={}
 for field in ('Mass','Length','Width','Height','MaximumTension','Kse','Kpe','B'):
  n=a.find(field)
  if n is not None and n.get('Actual') is not None and b.findtext(field) is not None:
   vals[field]=dict(aproj_actual=n.get('Actual'),asim=b.findtext(field))
 if vals:rows.append(dict(name=name(b),type=b.findtext('Type'),values=vals))
(DEST/'runs').mkdir(exist_ok=True)
(DEST/'runs/SOURCE_UNIT_COMPARISON.json').write_text(json.dumps(rows,indent=2))
print('PROJECT',DEST/'Biarticular_net.aproj');print('UNIT_EXAMPLES',json.dumps(rows[:5]))
for field in ('FileName','Filename','ChartFile','Texture','MeshFile','CollisionMeshFile'):
 for e in ap.iter(field):
  if e.text:print('DEPENDENCY',field,e.text)
for e in sim.iter('DataChart'):print('CHART',e.findtext('Name'),e.findtext('OutputFilename'),e.findtext('StartTime'),e.findtext('EndTime'))
