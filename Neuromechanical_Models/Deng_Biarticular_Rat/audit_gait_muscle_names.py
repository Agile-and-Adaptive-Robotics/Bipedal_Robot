"""Read-only name/route inventory for rat-to-Gait transfer; no guessed assignments."""
from pathlib import Path
import xml.etree.ElementTree as E
import json,hashlib
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];BASE=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody'
source=HERE/'Biarticular_net_Standalone.asim';r=E.parse(source).getroot();neurons={n.findtext('ID'):n.findtext('Name') for n in r.iter('Neuron')};bodies={n.findtext('ID'):n for n in r.iter('RigidBody')};rat=[]
for a in r.iter('Adapter'):
 b=bodies.get(a.findtext('TargetID'))
 if a.findtext('Type')=='NodeToPhysical' and b is not None and b.findtext('Type')=='LinearHillMuscle':rat.append(dict(muscle=b.findtext('Name'),MN=neurons.get(a.findtext('SourceID')),target_data=a.findtext('TargetDataType'),adapter_gain={k:a.find('Gain').findtext(k) for k in ('Type','A','B','C','D')}))
targets=[]
for p in [BASE/'gait2392_simbody.osim',BASE/'gait2392_robotbody.osim',BASE/'gait2392_robot.osim',BASE/'mjc/gait2392_simbody/gait2392_simbody_cvt3.xml']:
 if not p.exists():continue
 root=E.parse(p).getroot()
 if p.suffix=='.osim':
  elements=[e for e in root.iter() if e.tag.endswith('Muscle') and e.get('name')];names=[e.get('name') for e in elements]
  details=[dict(name=e.get('name'),kind=e.tag,bodies=sorted({x.text for x in e.iter() if x.tag in ('body','socket_parent_frame') and x.text}),max_isometric_force=e.findtext('max_isometric_force')) for e in elements]
 else:
  names=[e.get('name') for e in root.findall('actuator/*')];details=[dict(name=e.get('name'),kind=e.tag,tendon=e.get('tendon')) for e in root.findall('actuator/*')]
 targets.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),count=len(names),names=names,details=details))
out=dict(source=str(source),rat_outputs=rat,targets=targets,scope='Names and source route bodies only. No homology claim, conversion, or innervation applied.')
(HERE/'runs/RAT_GAIT_NAME_AUDIT.json').write_text(json.dumps(out,indent=2))
print(json.dumps(dict(rat_outputs=rat,targets=[{k:v for k,v in x.items() if k not in ('details',)} for x in targets]),indent=2))
