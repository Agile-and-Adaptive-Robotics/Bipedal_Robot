from pathlib import Path
import xml.etree.ElementTree as E,copy,json,uuid,shutil
from run_comparison import REPO,OUT,LEG,NEW
BASE=REPO/'Neuromechanical_Models/W2L_Treadmill_FrozenRoot/W2L_Treadmill_FrozenRoot.aproj'
def put(n,k,v):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 c.text=str(v);return c
def num(n,k,v):
 c=put(n,k,'');c.attrib.update(Value=str(v),Scale='None',Actual=str(v));return c
for model in ['legacy','bilateral']:
 for condition in ['air','treadmill','drop','partial50']:
  name=model+'_'+condition;folder=OUT/'recording_projects'/name;folder.mkdir(parents=True,exist_ok=True)
  r=E.parse(BASE if model=='legacy' else NEW.with_name('Walker_2_Layer_CPG_BilateralRG.aproj')).getroot()
  old=E.parse(BASE).getroot();env=r.find('.//Environment');legacyenv=old.find('.//Environment')
  structs=env.find('Structures');env.remove(structs);env.append(copy.deepcopy(legacyenv.find('Structures')))
  ground=next(s for s in env.find('Structures') if s.findtext('Name')=='Ground')
  num(ground.find('LocalPosition'),'Y',-2)
  put(ground,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,0,-20,0,1')
  gb=ground.find('RigidBody');num(gb.find('LocalPosition'),'Y',-2)
  lm=gb.find('LocalMatrix')
  if lm is not None:
   values=lm.text.split(',');values[13]='-20';lm.text=','.join(values)
  org=env.find('.//Organism');root=org.find('RigidBody');org.remove(root);root=copy.deepcopy(legacyenv.find('.//Organism/RigidBody'));org.append(root)
  put(root,'Freeze','True' if condition in ['air','treadmill'] else 'False')
  anchor=next(b for b in env.iter('RigidBody') if b.findtext('Name')=='WalkingPath');put(anchor,'IsVisible','True')
  path=next(s for s in env.find('Structures') if s.findtext('Name')=='Path')
  for k,v in zip('XYZ',[-3.454,-.3,0]):num(path.find('LocalPosition'),k,v)
  put(path,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,-34.54,-3,0,1')
  belt=next(b for b in env.iter('RigidBody') if b.findtext('Name')=='Treadmill Belt');put(belt,'Texture',str(LEG.parent/'Treadmill_checker.bmp'))
  if condition=='air':
   path=next(s for s in env.find('Structures') if s.findtext('Name')=='Path')
   num(path.find('LocalPosition'),'Y',-2);put(path,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,-34.54,-20,0,1');num(anchor.find('LocalPosition'),'Y',-2);put(anchor,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,-34.54,-20,0,1')
  stimuli=r.find('.//Stimuli')
  for st in list(stimuli):
   if st.findtext('ClassName','').endswith('MotorVelocity'):stimuli.remove(st)
  motor=copy.deepcopy(next(st for st in old.find('.//Stimuli') if st.findtext('ClassName','').endswith('MotorVelocity')))
  num(motor,'Velocity',.1 if condition=='treadmill' else 0);put(motor,'Equation',.1 if condition=='treadmill' else 0);stimuli.append(motor)
  if condition=='partial50':
   meta=json.loads((OUT/model/condition/'baseline/protocol.json').read_text());st=E.SubElement(stimuli,'Stimulus')
   for k,v in dict(AssemblyFile='AnimatGUI.dll',ClassName='AnimatGUI.DataObjects.ExternalStimuli.Force',Name='50 percent declared bodyweight support',ID=str(uuid.uuid4()),StructureID=org.findtext('ID'),PartID=root.findtext('ID'),AlwaysActive='False',Enabled='True',ValueType='Constant',Equation='0').items():put(st,k,v)
   for k,v in dict(StartTime=0,EndTime=10,StepInterval=0,PositionX=0,PositionY=0,PositionZ=0,ForceX=0,ForceY=meta['support_force_N'],ForceZ=0,TorqueX=0,TorqueY=0,TorqueZ=0).items():num(st,k,v)
  num(r.find('Simulation'),'SimEndTime',10.1);put(r,'ProjectName',name);put(r,'SimulationFile',name+'.asim')
  E.ElementTree(r).write(folder/(name+'.aproj'),encoding='utf-8',xml_declaration=True)
  sim=E.parse(OUT/model/condition/'baseline/trial.asim');sr=sim.getroot();put(sr,'ProjectPath',str(folder)+'\\');put(sr,'PlaybackControlMode','1')
  put(next(b for b in sr.iter('RigidBody') if b.findtext('Name')=='WalkingPath'),'IsVisible','True')
  sim.write(folder/(name+'.asim'),encoding='utf-8',xml_declaration=True)
  for f in BASE.parent.glob('*.aform'):shutil.copyfile(f,folder/f.name)
print(OUT/'recording_projects')

