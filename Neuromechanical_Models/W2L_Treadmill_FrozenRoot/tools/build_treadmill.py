from pathlib import Path
import xml.etree.ElementTree as E,copy,uuid,shutil,json
src=Path('Neuromechanical_Models/Walker_2_Layer_CPG');dst=Path('Neuromechanical_Models/W2L_Treadmill_FrozenRoot');dst.mkdir(exist_ok=True)
ids={k:str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-w2l-treadmill-'+k)) for k in ['belt','joint','stim']}
def put(n,k,v):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 c.text=str(v);return c
def num(n,k,v,gui):
 c=put(n,k,'' if gui else v)
 if gui:c.attrib.update(Value=str(v),Scale='None',Actual=str(v))
def pos(n,k,xyz,gui):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 if gui:
  for a,v in zip('XYZ',xyz):num(c,a,v,True)
 else:c.attrib.update(zip('xyz',map(str,xyz)))
for gui in [False,True]:
 tree=E.parse(src/('Walker_2_Layer_CPG.aproj' if gui else 'Walker_2_Layer_CPG_Standalone_modern.asim'));r=tree.getroot()
 anchor=next(n for n in r.iter('RigidBody') if n.findtext('Name')=='WalkingPath')
 root=next(n for n in r.iter('RigidBody') if n.findtext('Name')=='Root');put(root,'Freeze','True')
 ground=next(n for n in r.iter('Structure') if n.findtext('Name')=='Ground')
 pos(ground,'LocalPosition' if gui else 'Position',(0,-2,0),gui)
 if gui:put(ground,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,0,-20,0,1')
 gb=ground.find('RigidBody');pos(gb,'LocalPosition' if gui else 'Position',(0,-2,0),gui)
 if gui:
  lm=gb.find('LocalMatrix')
  if lm is not None:
   vals=lm.text.split(',');vals[13]='-20';lm.text=','.join(vals)
 structure=next(n for n in r.iter('Structure') if n.findtext('Name')=='Path')
 pos(structure,'LocalPosition' if gui else 'Position',(-3.454,-.3,0),gui)
 if gui:put(structure,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,-34.54,-3,0,1')
 pos(anchor,'LocalPosition' if gui else 'Position',(-3.454,-.3,0),gui)
 if gui:put(anchor,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,-34.54,-3,0,1')
 belt=copy.deepcopy(anchor)
 for child in list(belt):
  if child.tag in ['ChildBodies','Joint']:belt.remove(child)
 for n in belt.iter('ID'):n.text=str(uuid.uuid5(uuid.NAMESPACE_URL,'treadmill-clone-'+str(n.text)))
 put(belt,'Texture','Treadmill_checker.bmp');put(belt,'ID',ids['belt']);put(belt,'Name','Treadmill Belt');put(belt,'Freeze','False');put(belt,'IsVisible','True');put(belt,'IsCollisionObject','True')
 pos(belt,'LocalPosition' if gui else 'Position',(0,.26,0),gui)
 for k,v in [('Length',4),('Width',1),('Height',.09),('Mass',4000)]:num(belt,k,v,gui)
 # Preserve original floor height initially; chart foot positions before setting contact height.
 j=copy.deepcopy(next(r.iter('Joint')))
 for n in j.iter('ID'):n.text=str(uuid.uuid5(uuid.NAMESPACE_URL,'treadmill-clone-'+str(n.text)))
 put(j,'ID',ids['joint']);put(j,'Name','Treadmill Slide');put(j,'Type','Prismatic');put(j,'PartType','AnimatGUI.DataObjects.Physical.Joints.Prismatic')
 pos(j,'LocalPosition' if gui else 'Position',(0,0,0),gui);pos(j,'Rotation',(0,0,0),gui)
 if gui:
  for x in [j,belt]:
   lm=x.find('LocalMatrix')
   if lm is not None:x.remove(lm)
 if gui:put(belt,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,0,2.6,0,1');put(j,'LocalMatrix','1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1')
 for k,v in [('EnableLimits','True'),('EnableMotor','True'),('MotorType','VelocityControl' if gui else 0)]:put(j,k,v)
 num(j,'MaxForce',10000,gui);num(j,'MaxVelocity',.2,gui)
 for k,v in [('LowerLimit',-2),('UpperLimit',2)]:num(j.find(k),'LimitPos',v,gui)
 belt.insert(0,j)
 children=anchor.find('ChildBodies')
 if children is None:children=E.SubElement(anchor,'ChildBodies')
 children.append(belt);put(anchor,'IsCollisionObject','True');put(anchor,'IsVisible','True')
 st=E.Element('Stimulus')
 fields={'Name':'Backward belt drive','ID':ids['stim'],'AlwaysActive':'False','Enabled':'True','Equation':'0.1','DisableMotorWhenDone':'False'}
 for k,v in fields.items():put(st,k,v)
 for k,v in [('StartTime',.5),('EndTime',10)]:num(st,k,v,gui)
 if gui:
  for k,v in {'AssemblyFile':'AnimatGUI.dll','ClassName':'AnimatGUI.DataObjects.ExternalStimuli.MotorVelocity','StructureID':'f18d21af-4e36-4ae6-93ad-2ce71937dba7','PartID':ids['joint'],'ValueType':'Constant'}.items():put(st,k,v)
  num(st,'StepInterval',0,True);num(st,'Velocity',.1,True);r.find('.//Stimuli').append(st)
  num(r.find('Simulation'),'SimEndTime',10.1,True);put(r,'ProjectName','W2L_Treadmill_FrozenRoot');put(r,'SimulationFile','W2L_Treadmill_FrozenRoot.asim')
 else:
  for k,v in {'Type':'MotorVelocity','StructureID':'f18d21af-4e36-4ae6-93ad-2ce71937dba7','JointID':ids['joint'],'TargetID':'Velocity'}.items():put(st,k,v)
  r.find('ExternalStimuli').append(st);put(r,'ProjectPath',str(dst.resolve())+'\\');put(r,'SimEndTime','10.1');put(r,'PlaybackControlMode','0')
  charts=r.find('DataCharts');charts.clear();chart=E.SubElement(charts,'DataChart')
  for k,v in {'ID':str(uuid.uuid4()),'Type':'FileChart','Name':'Treadmill verification','OutputFilename':'Treadmill verification.txt','SetStartEndTime':'True','Enabled':'True','CollectTimeWindow':'-1','AlwaysActive':'False','StartTime':'0','EndTime':'10','CollectInterval':'.01'}.items():put(chart,k,v)
  cols=E.SubElement(chart,'DataColumns')
  targets=[(ids['belt'],'WorldPositionX','belt_x'),(ids['belt'],'WorldPositionY','belt_y'),(root.findtext('ID'),'WorldPositionY','root_y'),(ids['joint'],'JointPosition','slide_pos'),(ids['joint'],'JointActualVelocity','slide_velocity'),(ids['joint'],'JointSetVelocity','slide_target')]
  for b in r.iter('RigidBody'):
   if 'foot' in b.findtext('Name','').lower() or 'toe' in b.findtext('Name','').lower():targets.append((b.findtext('ID'),'WorldPositionY',b.findtext('Name')+'_y'))
  for b in r.iter('RigidBody'):
   if b.findtext('IsContactSensor')=='True':targets.append((b.findtext('ID'),'ContactCount',b.findtext('Name')+'_contacts'))
  for i,(target,typ,name) in enumerate(targets):
   c=E.SubElement(cols,'DataColumn')
   for k,v in {'ID':str(uuid.uuid4()),'ColumnName':name,'Type':'DataColumn','TargetID':target,'DataType':typ}.items():put(c,k,v)
  (dst/'diagnostic_columns.json').write_text(json.dumps(targets))
 path=dst/('W2L_Treadmill_FrozenRoot.aproj' if gui else 'W2L_Treadmill_FrozenRoot_Standalone.asim');tree.write(path,encoding='utf-8',xml_declaration=True)
 if not gui:shutil.copyfile(path,dst/'W2L_Treadmill_FrozenRoot.asim')
for f in src.glob('*.aform'):shutil.copyfile(f,dst/f.name)
print(dst.resolve())










