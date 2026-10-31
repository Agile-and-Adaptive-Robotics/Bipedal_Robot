from pathlib import Path
from lxml import etree as E
import json,hashlib
R=Path(__file__).resolve().parents[3];O=R/'Neuromechanical_Models/W2L_VirtualWalker_20261031';rows=[]
expected20={'W2L_VirtualWalker_Newer_SIContact.aproj':'34a416cfce1bfd310889383de597b194f392d0acfce4f1fd85562d0f7608111d','W2L_VirtualWalker_Newer_SIContact_Standalone.asim':'171c86a3b5b678a963961ad7755c0e61a72dd030bb6ada42269d500647f451cd'}
for tag in ('SIContact','GUI1nA'):
 name='W2L_VirtualWalker_Newer_'+tag;W=R/('Neuromechanical_Models/'+name+'_20261031');files={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [W/(name+'.aproj'),W/(name+'_Standalone.asim'),W/(name+'.asim'),W/'NATIVE_SMOKE.json',W/'BUILD_SYNC.json']}
 for f,sha in files.items():
  if Path(f).name in expected20:assert sha==expected20[Path(f).name]
 p=E.parse(str(W/(name+'.aproj'))).getroot();n=E.parse(str(W/(name+'_Standalone.asim'))).getroot();gain=[]
 for a in p.findall('.//NervousSystem//Node'):
  if a.findtext('DataTypeID')=='ContactCount':gain.append({'adapter':a.findtext('ID'),'C_attributes':dict(a.find('Gain/C').attrib),'D_attributes':dict(a.find('Gain/D').attrib)})
 tonic={a.findtext('Name'):float(a.findtext('TonicStimulus')) for a in n.iter('Neuron') if a.findtext('Name') in ('L RG ext','R RG ext')}
 rows.append({'variant':tag,'hashes':files,'GUI_contact_coefficients':gain,'tonic_nA':tonic,'results':json.loads((W/'NATIVE_SMOKE.json').read_text())})
(O/'CONTACT_CANDIDATE_SNAPSHOTS.json').write_text(json.dumps({'machine':'easteregg2','tested_20nA_model_bytes_preserved':True,'no_candidate_claimed_as_user_intent':True,'candidates':rows},indent=2));print('Recorded candidate snapshots; original tested20nA model hashes unchanged.')
