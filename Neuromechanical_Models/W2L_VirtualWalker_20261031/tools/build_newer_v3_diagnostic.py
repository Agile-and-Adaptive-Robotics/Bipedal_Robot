from pathlib import Path
from historical_build_guard import require_historical_regeneration
require_historical_regeneration()
from lxml import etree as E
import copy,json,uuid,hashlib,shutil,subprocess
R=Path(__file__).resolve().parents[3];TOOLS=Path(__file__).parent;S=R/'Neuromechanical_Models/W2L_VirtualWalker_Newer_GUI1nA_20261031';NAME='W2L_VirtualWalker_NewerV3Guided';W=R/('Neuromechanical_Models/'+NAME+'_20261031')
subprocess.run([r'D:\Anaconda\envs\myo\python.exe',str(TOOLS/'build_guided_native.py'),'--historical-regeneration','--source-folder',str(S),'--source-name','W2L_VirtualWalker_Newer_GUI1nA','--variant-id','NewerV3Guided','--clearance','.03','--rail-drop','.02','--neural-clock','2'],check=True,capture_output=True)
EDGES=['cafe0027-0000-4000-8000-900000000027','cafe0029-0000-4000-8000-900000000029'];oldtype='9baf7f2a-bac2-4c80-b398-9dce2ab81f03';zeroid=str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-newer-v3-diagnostic-output-zero-20261031'));P=E.XMLParser(strip_cdata=False)
def put(n,k,v):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 z.text=str(v);return z
def zero_scaled(n,k):
 z=put(n,k,'');z.attrib.update({'Value':'0','Scale':'micro','Actual':'0'})
changes=[]
for gui,ext in [(True,'.aproj'),(False,'_Standalone.asim')]:
 f=W/(NAME+ext);tree=E.parse(str(f),P);root=tree.getroot();nodes=root.findall('.//NervousSystem//Node') if gui else list(root.iter('Neuron'))
 for n in nodes:
  if n.findtext('Text' if gui else 'Name') in ('L RG ext','R RG ext'):
   if gui:put(n,'TonicStimulus','').attrib.update({'Value':'0','Scale':'nano','Actual':'0'})
   else:put(n,'TonicStimulus','0')
 types=root.findall('.//SynapseTypes/Link') if gui else list(root.iter('SynapseType'));src=next(t for t in types if t.findtext('ID')==oldtype);z=copy.deepcopy(src);put(z,'ID',zeroid);put(z,'Name','DIAGNOSTIC V3 return off; not final architecture')
 if gui:zero_scaled(z,'MaxSynapticConductance')
 else:put(z,'SynAmp','0')
 src.getparent().append(z)
 links=root.findall('.//NervousSystem//Links/Link') if gui else list(root.iter('Connexion'));matched=[]
 for l in links:
  if l.findtext('ID') not in EDGES:continue
  assert l.findtext('SynapticTypeID' if gui else 'SynapseTypeID')==oldtype;put(l,'SynapticTypeID' if gui else 'SynapseTypeID',zeroid)
  if gui:zero_scaled(l,'SynapticConductance');put(l,'UserText','DIAGNOSTIC disabled V3 return, g=0')
  else:put(l,'G','0')
  matched.append(l.findtext('ID'))
 assert sorted(matched)==sorted(EDGES)
 changes.append({'gui':gui,'edges_disabled':matched,'new_zero_type':zeroid,'original_shared_type_retained':oldtype})
 tree.write(str(f),encoding='utf-8',xml_declaration=True)
shutil.copy2(W/(NAME+'_Standalone.asim'),W/(NAME+'.asim'))
reference=R/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031/runs/W2L_VirtualWalker_Newer_GUI1nA/newer_v3_diag_v1_ground_guides'
for filename in ['PARAMETERS.json','V3_DIAGNOSTIC.json']:
 if (reference/filename).exists():shutil.copy2(reference/filename,W/('MUJOCO_REFERENCE_'+filename))
report={'machine':'easteregg2','scope':'Separately named diagnostic ablation, NOT final V3 architecture or literature replica','source':str(S),'graph':'GUI-authoritative96neurons/204connections;2zeroeffectiveedgesretained','contact_gain_C_A_per_contact':1e-9,'contact_offset_D_A':0,'RG_ext_tonic_nA_LR':[0,0],'c1_paths':'retained','neural_clock':2,'clearance_m':.03,'horizontal_rail_offset_m':-.02,'harness':'none; physicalguidecarriageonly','perturbation':'none','changes':changes,'mujoCo_reference':str(reference),'source_hashes':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [S/'W2L_VirtualWalker_Newer_GUI1nA.aproj',S/'W2L_VirtualWalker_Newer_GUI1nA_Standalone.asim']}}
(W/'V3_DIAGNOSTIC_BUILD.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
