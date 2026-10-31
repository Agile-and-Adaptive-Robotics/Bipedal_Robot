"""Build a free-root physical support rig on the separately corrected W2L copy."""
from pathlib import Path
from lxml import etree as E
import copy, json, uuid, hashlib, shutil, math
import numpy as np
import sys

W=Path(__file__).resolve().parents[1]
REPO=W.parents[1]
SOURCE=W.parent/'W2L_Geometry_20261031'
SMOOTH='09bccf17-3f71-4b1b-b397-53550cb1bd6b'
NAME='W2L_VirtualWalker'
RAIL_FRICTION=.05
RAIL_MATERIAL=str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-virtual-walker-20261031-rail-material'))

def uid(label):return str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-virtual-walker-20261031-'+label))
def put(n,k,v):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 c.text=str(v)
 return c
def num(n,k,v,gui):
 c=put(n,k,'' if gui else v)
 if gui:c.attrib.update({'Value':str(v),'Scale':'None','Actual':str(v)})
def vector(n,k,v,gui):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 if gui:
  for a,x in zip('XYZ',v):num(c,a,x,True)
 else:c.attrib.update(zip('xyz',map(str,v)))
def position(n,v,gui):
 vector(n,'LocalPosition' if gui else 'Position',v,gui)
 if gui:
  # GUI stores LocalMatrix translations in decimeters, column-major.
  put(n,'LocalMatrix',','.join(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,*[x*10 for x in v],1])))
def xyz(n,k):
 c=n.find(k)
 return np.zeros(3) if c is None else np.array([float(c.get(a,'0')) for a in 'xyz'])
def rot(v):
 x,y,z=v;cx,sx=math.cos(x),math.sin(x);cy,sy=math.cos(y),math.sin(y);cz,sz=math.cos(z),math.sin(z)
 return np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])@np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])@np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
def worlds(root):
 out={}
 def visit(b,p=np.zeros(3),R=np.eye(3)):
  p=p+R@xyz(b,'Position');R=R@rot(xyz(b,'Rotation'));out[b.findtext('ID')]=(p,R,b)
  for child in b.findall('ChildBodies/RigidBody'):visit(child,p,R)
 visit(root)
 return out
def bottom(p,R,b):
 size=np.array([float(b.findtext(k)) for k in ['Length','Height','Width']])
 return float(p[1]-np.abs(R[1])@size/2)

def configure(src,gui,meta):
 tree=E.parse(str(src));r=tree.getroot()
 root=next(b for b in r.iter('RigidBody') if b.findtext('Name')=='Root')
 path=next(b for b in r.iter('RigidBody') if b.findtext('Name')=='WalkingPath')
 materials=r.find('.//MaterialTypes');material=copy.deepcopy(next(m for m in materials if m.findtext('ID')==SMOOTH))
 put(material,'ID',RAIL_MATERIAL);put(material,'Name','Virtual walker sliding contact')
 for child in list(material):
  if child.tag.startswith('FrictionLinear'):num(material,child.tag,RAIL_FRICTION,gui)
 materials.append(material)
 put(root,'Freeze','False')
 children=root.find('ChildBodies')
 if children is None:children=E.SubElement(root,'ChildBodies')
 # Ordinary child boxes without a Joint are rigid parts of the parent body;
 # the original foot/toe contact boxes use this compound-body convention.
 for sign,side in [(1,'left'),(-1,'right')]:
  wing=copy.deepcopy(root)
  for child in list(wing):
   if child.tag in ['ChildBodies','Joint']:wing.remove(child)
  put(wing,'ID',uid('wing-'+side));put(wing,'Name',side.title()+' hip wing')
  position(wing,(0,-.04,sign*.25),gui);vector(wing,'Rotation',(0,0,0),gui);vector(wing,'COM',(0,0,0),gui)
  for k,v in [('Length',.12),('Height',.025),('Width',.30),('Mass',.001),('Density',.001/(.12*.025*.30*1e6))]:num(wing,k,v,gui)
  for k,v in [('Freeze','False'),('IsContactSensor','False'),('IsCollisionObject','True'),('IsVisible','True'),('MaterialTypeID',RAIL_MATERIAL),('Description','Rigid lateral pelvis support; 0.001 g declared mass; no added joint.')]:put(wing,k,v)
  for color in ['Ambient','Diffuse','Specular']:
   wing.find(color).attrib.update({'Red':'.3','Green':'.8','Blue':'.85','Alpha':'1'})
  children.append(wing)
 structures=r.find('.//Environment/Structures')
 for sign,side in [(1,'left'),(-1,'right')]:
  structure=E.Element('Structure')
  put(structure,'ID',uid('rail-structure-'+side));put(structure,'Name',side.title()+' virtual-walker rail');put(structure,'Type','Basic');put(structure,'IsVisible','False')
  v=(meta['root_initial_xyz'][0],meta['rail_top_y']-.02,sign*.32)
  position(structure,v,gui);vector(structure,'Rotation',(0,0,0),gui)
  E.SubElement(structure,'CollisionExclusionPairs')
  rail=copy.deepcopy(path)
  for child in list(rail):
   if child.tag in ['ChildBodies','Joint']:rail.remove(child)
  put(rail,'ID',uid('rail-'+side));put(rail,'Name',side.title()+' support rail')
  position(rail,v,gui);vector(rail,'Rotation',(0,0,0),gui);vector(rail,'COM',(0,0,0),gui)
  for k,v in [('Length',4),('Height',.04),('Width',.18),('Mass',9e8)]:num(rail,k,v,gui)
  for k,v in [('Freeze','True'),('IsVisible','False'),('IsCollisionObject','True'),('IsContactSensor','False'),('MaterialTypeID',RAIL_MATERIAL),('Texture',''),('Description','Invisible, collision-enabled support rail along world x; low-friction material permits sliding.')]:put(rail,k,v)
  if gui:
   for q in rail.findall('Transparencies/*'):q.text='100'
  else:
   for k in ['GraphicsAlpha','SimulationAlpha']:put(rail,k,0)
  structure.append(rail);structures.append(structure)
 # Keep the existing moving floor and its drive. Only Root support changes.
 belt=next(b for b in r.iter('RigidBody') if b.findtext('Name')=='Treadmill Belt')
 put(belt,'Texture',str((W/'Treadmill_checker.bmp').resolve()))
 if gui:
  put(r,'ProjectName',NAME);put(r,'SimulationFile',NAME+'.asim');put(r.find('Simulation'),'ProjectPath',str(W.resolve())+'\\');num(r.find('Simulation'),'SimEndTime',10.1,True)
 else:
  put(r,'ProjectPath',str(W.resolve())+'\\');put(r,'SimEndTime',10.1);put(r,'PlaybackControlMode','1')
  windows=r.find('WindowMgr/Windows')
  if windows is not None:windows.clear()
  charts=r.find('DataCharts');charts.clear();chart=E.SubElement(charts,'DataChart')
  for k,v in dict(ID=uid('chart'),Type='FileChart',Name='Virtual walker verification',OutputFilename='Virtual walker verification.txt',SetStartEndTime='True',Enabled='True',CollectTimeWindow='-1',AlwaysActive='False',StartTime=0,EndTime=10,CollectInterval=.002).items():put(chart,k,v)
  columns=E.SubElement(chart,'DataColumns');signals=[]
  for axis in 'XYZ':signals.append((root.findtext('ID'),'WorldPosition'+axis,'root_'+axis.lower()))
  for b in root.iter('RigidBody'):
   name=b.findtext('Name','')
   if 'hip wing' in name or name in ['foot_L','foot_R','toe_L','toe_R']:
    for axis in 'XYZ':signals.append((b.findtext('ID'),'WorldPosition'+axis,name+'_'+axis.lower()))
   if b.findtext('IsContactSensor')=='True':signals.append((b.findtext('ID'),'ContactCount',name+'_contacts'))
  for j in root.iter('Joint'):
   if j.findtext('Type')=='Hinge':signals.append((j.findtext('ID'),'JointRotation',j.findtext('Name')+'_angle'))
  for neuron in r.iter('Neuron'):
   if 'MN' in neuron.findtext('Name',''):signals.append((neuron.findtext('ID'),'MembraneVoltage',neuron.findtext('Name')))
  counts={};unique_signals=[]
  for target,kind,name in signals:
   counts[name]=counts.get(name,0)+1
   unique_signals.append((target,kind,name if counts[name]==1 else name+'_'+str(counts[name])))
  signals=unique_signals
  for target,kind,name in signals:
   c=E.SubElement(columns,'DataColumn')
   for k,v in dict(ID=uid('column-'+target+'-'+kind),Type='DataColumn',ColumnName=name,TargetID=target,DataType=kind).items():put(c,k,v)
  (W/'diagnostic_columns.json').write_text(json.dumps(['TimeSlice','Time']+[s[2] for s in signals],indent=2),encoding='utf-8')
 out=W/(NAME+('.aproj' if gui else '_Standalone.asim'))
 tree.write(str(out),encoding='utf-8',xml_declaration=True)
 # Separate parse catches broken neural drawing CDATA without changing it.
 for n in r.iter('DiagramXml'):
  if n.text and n.text.strip():E.fromstring(n.text.encode())
 ids=[n.text for n in r.iter('ID')]
 if not gui:assert len(ids)==len(set(ids)), 'Duplicate object IDs'
 for label in ['wing-left','wing-right','rail-left','rail-right','rail-structure-left','rail-structure-right']:
  assert ids.count(uid(label))==1, 'Duplicate support object ID: '+label
 return out

def main():
 global W,SOURCE,NAME
 if '--newer' in sys.argv:
  W=W.parent/'W2L_VirtualWalker_Newer_20261031';W.mkdir(exist_ok=True)
  SOURCE=SOURCE.parent/'W2L_Geometry_Newer_20261031';NAME='W2L_VirtualWalker_Newer'
  shutil.copyfile(Path(__file__).resolve().parents[1]/'.gitignore',W/'.gitignore')
 source_asim=SOURCE/'W2L_Geometry_Standalone.asim';source_aproj=SOURCE/'W2L_Geometry.aproj'
 r=E.parse(str(source_asim));root=next(b for b in r.iter('RigidBody') if b.findtext('Name')=='Root')
 ws=worlds(root);feet=[(b.findtext('Name'),bottom(p,R,b)) for p,R,b in ws.values() if b.findtext('Name') in ['foot_L','foot_R','toe_L','toe_R'] and b.findtext('Type')=='Box']
 path=next(b for b in r.iter('RigidBody') if b.findtext('Name')=='WalkingPath');pathworlds=worlds(path)
 belt=next((p,R,b) for p,R,b in pathworlds.values() if b.findtext('Name')=='Treadmill Belt');bp,bR,bb=belt
 floor_top=float(bp[1]+np.abs(bR[1])@np.array([float(bb.findtext(k)) for k in ['Length','Height','Width']])/2)
 rp=xyz(root,'Position');rr=rot(xyz(root,'Rotation'));wingcenter=rp+rr@np.array([0,-.04,.25]);wingbottom=float(wingcenter[1]-np.abs(rr[1])@np.array([.12,.025,.30])/2)
 gap=min(y for n,y in feet)-floor_top
 if gap<-.0005:raise ValueError('Initial physical foot overlaps floor: '+str(gap))
 clearance=max(0,gap)+.012
 meta=dict(source={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source_asim,source_aproj]},axes='x anterior/travel, y superior, z lateral',root_initial_xyz=rp.tolist(),root_frozen=False,foot_sole_y=dict(feet),floor_top_y=floor_top,floor_world_rotation=[0,0,0],minimum_foot_gap=gap,wing_bottom_y=wingbottom,wing_to_rail_gap=clearance,extra_support_clearance=.012,rail_top_y=wingbottom-clearance,wing_dimensions_xyz=[.12,.025,.30],wing_centers_root_xyz=[[0,-.04,.25],[0,-.04,-.25]],wing_declared_mass_g=.001,rail_dimensions_xyz=[4,.04,.18],rail_centers_z=[.32,-.32],rail_visible=False,rail_collision_enabled=True,material=dict(id=RAIL_MATERIAL,linear_friction=RAIL_FRICTION,restitution=0),support='passive normal contact, no pelvis joint or freeze; roll/pitch/yaw remain possible',recording_speed='MatchPhysicsStep (normal speed)')
 geometry_report=SOURCE/'GEOMETRY_REPORT.json'
 if geometry_report.exists():meta['poses_deg']=json.loads(geometry_report.read_text(encoding='utf-8')).get('poses_deg',{})
 for f in SOURCE.glob('*.aform'):shutil.copyfile(f,W/f.name)
 texture=SOURCE/'Treadmill_checker.bmp'
 if texture.exists():shutil.copyfile(texture,W/texture.name)
 else:shutil.copyfile(W.parent/'W2L_Treadmill_FrozenRoot/Treadmill_checker.bmp',W/'Treadmill_checker.bmp')
 configure(source_asim,False,meta);configure(source_aproj,True,meta);shutil.copyfile(W/(NAME+'_Standalone.asim'),W/(NAME+'.asim'))
 (W/'BUILD.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
 print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
