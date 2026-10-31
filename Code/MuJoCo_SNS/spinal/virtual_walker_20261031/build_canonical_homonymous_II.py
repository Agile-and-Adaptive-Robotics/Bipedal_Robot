"""One canonical editable NEWER:12 independent homonymous II relays, no family copy."""
from pathlib import Path
import sys,json,copy,uuid,hashlib,zipfile,importlib.util
from lxml import etree as E
import numpy as np
from import_virtual_walker import HERE,REPO,mujoco,safe
from native_kernel import NativeController
folder=REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031';project=folder/'W2L_VirtualWalker_Newer.aproj';runtime=folder/'W2L_VirtualWalker_Newer_Standalone.asim';raw=HERE/'runs/canonical_12II';raw.mkdir(exist_ok=True);parser=E.XMLParser(strip_cdata=False)
def uid(s):return str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-canonical-newer-homonymousII-20261031/'+s))
def put(n,k,v):
 x=n.find(k)
 if x is None:x=E.SubElement(n,k)
 x.text=str(v);return x
def scaled(n,k,actual,scale='None',value=None):
 x=put(n,k,'');x.attrib.update({'Actual':str(actual),'Scale':scale,'Value':str(actual if value is None else value)});return x
def val(n,k):return float(n.find(k).get('Actual'))
def reroll(n,label):
 for k,x in enumerate(n.iter('ID')):x.text=uid(label+'/id/'+str(k))
pt=E.parse(str(project),parser);st=E.parse(str(runtime),parser);p=pt.getroot();s=st.getroot();before=dict(project_sha256=hashlib.sha256(project.read_bytes()).hexdigest(),standalone_sha256=hashlib.sha256(runtime.read_bytes()).hexdigest())
with zipfile.ZipFile(raw/'BEFORE_12II_CANONICAL.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in [project,runtime,folder/'W2L_VirtualWalker_Newer.asim']:z.write(f,f.name)
assert zipfile.ZipFile(raw/'BEFORE_12II_CANONICAL.zip').testzip() is None
nodes={n.findtext('ID'):n for n in p.findall('.//NervousSystem//Node')};assert len([n for n in nodes.values() if '.Neurons.' in n.findtext('ClassName','')])==96,'Builder is one-time; do not augment twice'
base=HERE/'runs/W2L_VirtualWalker_Newer_GUI1nA';meta=json.loads((base/'IMPORT_AUDIT.json').read_text());bank=meta['muscles'];m=mujoco.MjModel.from_xml_path(str(base/'W2L_VirtualWalker_Newer_GUI1nA.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);jointnames=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R'];signatures={}
for b in bank:
 tid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+safe(b['name']));length=float(d.ten_length[tid]);arms=[]
 for name in jointnames:
  j=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name);q=m.jnt_qposadr[j];d.qpos[:]=m.qpos0;d.qpos[q]=1e-5;mujoco.mj_forward(m,d);plus=float(d.ten_length[tid]);d.qpos[q]=-1e-5;mujoco.mj_forward(m,d);minus=float(d.ten_length[tid]);arms.append(-(plus-minus)/2e-5)
 d.qpos[:]=m.qpos0;mujoco.mj_forward(m,d);signatures[b['id']]=np.r_[length,arms]
core=NativeController(runtime,bank);mn_by_muscle={bank[target]['id']:next(i for i,j in core.ids.items() if j==pre) for pre,target,coeff in core.outputs};sensors=[b for b in bank if b['type']=='LinearHillStretchReceptor'];muscles=[b for b in bank if b['type']=='LinearHillMuscle'];mapping=[]
for sensor in sensors:
 candidates=sorted([(float(np.max(abs(signatures[sensor['id']]-signatures[a['id']]))),a) for a in muscles],key=lambda x:x[0]);error,muscle=candidates[0];assert error<1e-6,(sensor['name'],'No exact actuator-route match',error)
 arm=signatures[sensor['id']][1:];ix=int(np.argmax(abs(arm)));joint=jointnames[ix];side=joint[-1];part=joint.split('_')[0];action=('dorsiflexor' if arm[ix]>0 else 'plantarflexor') if part=='ankle' else ('flexor' if arm[ix]>0 else 'extensor');mnid=mn_by_muscle[muscle['id']];mn=nodes[mnid];physical=next((n for n in nodes.values() if n.findtext('LinkedBodyPartID')==sensor['id']),None);mapping.append(dict(side=side,joint=part,action=action,sensor=sensor['name'],sensor_id=sensor['id'],sensor_GUI_id=None if physical is None else physical.findtext('ID'),actuator=muscle['name'],actuator_id=muscle['id'],MN_id=mnid,MN_name=mn.findtext('Text'),route_signature_max_error=error,moment_arm_m=float(arm[ix]),II_gain=sensor['ii_discharge'],L0_m=sensor['resting_length']))
assert len(mapping)==12 and len({(x['side'],x['joint'],x['action']) for x in mapping})==12
types=p.find('.//SynapseTypes')
for t in types:
 if t.findtext('ID')=='6b6ead54-53be-4d79-babe-5937a67057f8':put(t,'Name','Legacy MN to RE Excite')
typ=copy.deepcopy(next(n for n in types if n.findtext('ClassName','').endswith('NonSpikingChemical')));typeid=uid('homonymous-II-MN-type');put(typ,'ID',typeid);put(typ,'Name','Homonymous Type II excitation, weak');scaled(typ,'MaxSynapticConductance',.005e-6,'micro',.005);scaled(typ,'EquilibriumPotential',-.010,'milli',-10);scaled(typ,'PreSynapticThreshold',-.060,'milli',-60);scaled(typ,'PreSynapticSaturationLevel',-.040,'milli',-40);types.append(typ)
neuron_template=next(n for n in nodes.values() if n.findtext('Text','').endswith(' II'));adapter_template=next(n for n in nodes.values() if n.findtext('ClassName','').endswith('PhysicalToNodeAdapter'));syn_template=next(l for l in p.findall('.//NervousSystem//Links/Link') if l.findtext('ClassName')=='IntegrateFireGUI.DataObjects.Behavior.Synapse');adapter_link_template=next(l for l in p.findall('.//NervousSystem//Links/Link') if l.findtext('ClassName','').endswith('Adapter'))
for row in mapping:
 key=row['side']+'/'+row['joint']+'/'+row['action'];name=f"{row['side']} {row['joint'].title()} {row['action']} II relay";nid=uid(key+'/relay');aid=uid(key+'/adapter');sid=uid(key+'/MN-synapse');lin=uid(key+'/physical-adapter-link');lout=uid(key+'/adapter-relay-link');home=nodes[row['MN_id']].getparent();links=home.getparent().find('Links');assert links is not None
 if row['sensor_GUI_id'] is None:
  physical=copy.deepcopy(next(n for n in nodes.values() if n.findtext('ClassName','').endswith('StretchReceptor')));reroll(physical,key+'/physical-sensor');pid=uid(key+'/physical-sensor');put(physical,'ID',pid);put(physical,'Text',row['sensor']);put(physical,'LinkedBodyPartID',row['sensor_id']);physical.find('InLinks').clear();physical.find('OutLinks').clear();home.append(physical);nodes[pid]=physical;row['sensor_GUI_id']=pid
 n=copy.deepcopy(neuron_template);reroll(n,key+'/relay');put(n,'ID',nid);put(n,'Text',name);n.find('InLinks').clear();n.find('OutLinks').clear();E.SubElement(n.find('InLinks'),'ID').text=lout;E.SubElement(n.find('OutLinks'),'ID').text=sid;put(n,'ToolTip','Independent homonymous Type II excitatory relay; physical route/action verified. Legacy phase-II pathway is separate.');n.find('Location').attrib.update({'x':'900','y':'500'});home.append(n)
 a=copy.deepcopy(adapter_template);reroll(a,key+'/adapter');put(a,'ID',aid);put(a,'Text',name+' adapter');put(a,'OriginID',row['sensor_GUI_id']);put(a,'DestinationID',nid);put(a,'DataTypeID','II');put(a,'TargetDataTypeID','ExternalCurrent');a.find('InLinks').clear();a.find('OutLinks').clear();E.SubElement(a.find('InLinks'),'ID').text=lin;E.SubElement(a.find('OutLinks'),'ID').text=lout;gain=a.find('Gain');C=20e-9/(row['II_gain']*row['L0_m']);
 for k,v in zip('ABCD',[0.,0.,C,0.]):scaled(gain,k,v)
 a.find('Location').attrib.update({'x':'740','y':'500'});home.append(a)
 syn=copy.deepcopy(syn_template);put(syn,'ID',sid);put(syn,'OriginID',nid);put(syn,'DestinationID',row['MN_id']);put(syn,'SynapticTypeID',typeid);scaled(syn,'SynapticConductance',.005e-6,'micro',.005);put(syn,'UserText','Type II homonymous excitation, .005 uS');links.append(syn);E.SubElement(nodes[row['MN_id']].find('InLinks'),'ID').text=sid
 for ident,origin,target in [(lin,row['sensor_GUI_id'],aid),(lout,aid,nid)]:
  link=copy.deepcopy(adapter_link_template);put(link,'ID',ident);put(link,'OriginID',origin);put(link,'DestinationID',target);links.append(link)
 E.SubElement(nodes[row['sensor_GUI_id']].find('OutLinks'),'ID').text=lin;row.update(relay_id=nid,relay_name=name,adapter_id=aid,synapse_id=sid,home_MN=row['MN_name'],adapter_C_A_per_discharge=C,adapter_D_A=0.,MN_gmax_uS=.005)
# Actual runnable canonical calibration, not runtime-only overrides.
for n in nodes.values():
 if n.findtext('DataTypeID')=='ContactCount' and n.findtext('ClassName','').endswith('PhysicalToNodeAdapter'):scaled(n.find('Gain'),'C',1e-9,'nano',1);scaled(n.find('Gain'),'D',0.,'nano',0)
 if n.findtext('Text') in ('L RG ext','R RG ext'):scaled(n,'TonicStimulus',0.,'nano',0)
v3type=copy.deepcopy(next(n for n in types if n.findtext('ID')=='9baf7f2a-bac2-4c80-b398-9dce2ab81f03'));v3id=uid('positive-V3-type');put(v3type,'ID',v3id);put(v3type,'Name','V3 contralateral RG ext IN, weak positive');scaled(v3type,'MaxSynapticConductance',.02749e-6,'micro',.02749);types.append(v3type)
for link in p.findall('.//NervousSystem//Links/Link'):
 if link.findtext('ID') in ('cafe0027-0000-4000-8000-900000000027','cafe0029-0000-4000-8000-900000000029'):put(link,'SynapticTypeID',v3id);scaled(link,'SynapticConductance',.02749e-6,'micro',.02749)
for n in p.findall('.//NervousSystem//Node'):
 if '.Neurons.' in n.findtext('ClassName',''):
  for path in ('TimeConstant','CaActivation/TimeConstant','CaDeactivation/TimeConstant'):
   z=n.find(path)
   if z is not None:z.set('Value',str(float(z.get('Value'))*2));z.set('Actual',str(float(z.get('Actual'))*2))
# Reuse only the previously verified physical guide payload, never its old neural graph.
archive=HERE/'runs/archives/w2l_generated_copies_before_cleanup_20261031.zip'
with zipfile.ZipFile(archive) as z:
 gui_member=next(x for x in z.namelist() if x.endswith('/W2L_VirtualWalker_NewerV3Guided.aproj'));gp=E.fromstring(z.read(gui_member),parser);gs=E.parse(str(HERE/'runs/source_inputs/W2L_VirtualWalker_NewerV3Guided_Standalone.asim'),parser).getroot()
for dest,src in [(p,gp),(s,gs)]:
 for name in ('Left support rail','Right support rail'):
  old=next(x for x in dest.iter('RigidBody') if x.findtext('Name')==name);new=next(x for x in src.iter('RigidBody') if x.findtext('Name')==name);old.getparent().replace(old,copy.deepcopy(new))
pt.write(str(project),encoding='utf-8',xml_declaration=True);st.write(str(runtime),encoding='utf-8',xml_declaration=True)
spec=importlib.util.spec_from_file_location('safe_canonical_export',str(REPO/'tmp/w2l_lh_style_20261031/export_current_canonical.py'));module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);export=module.synchronize(project,runtime,[runtime,folder/'W2L_VirtualWalker_Newer.asim'])
ledger=dict(machine='easteregg2',before=before,after=dict(project_sha256=hashlib.sha256(project.read_bytes()).hexdigest(),standalone_sha256=hashlib.sha256(runtime.read_bytes()).hexdigest()),mapping=mapping,legacy_phase_II='All four prior phase channels and their existing adapter/neuronal edges retained unchanged.',new_II_relay_normalization='I(A)=20e−9*II_discharge/(II_gain*L0);D0. New engineering choice, separate from legacy phase gains.',calibration=dict(contact_A_per_count=1e-9,RGext_tonic_nA_LR=[0,0],V3_output_gmax_uS=.02749,source_neuron_timeconstants_scale=2,guide_clearance_m=.03,rail_height_offset_m=-.02),export=export)
(raw/'BUILD_LEDGER.json').write_text(json.dumps(ledger,indent=2));(folder/'CANONICAL_WORKING_PROFILE.json').write_text(json.dumps(dict(neural_clock_scale=2,guide_clearance_m=.03,rail_height_offset_m=-.02,II_homonymous_relays=12,II_MN_gmax_uS=.005,source_actual_parameters=True),indent=2));print(json.dumps(dict(mapping=mapping,after=ledger['after']),indent=2))
