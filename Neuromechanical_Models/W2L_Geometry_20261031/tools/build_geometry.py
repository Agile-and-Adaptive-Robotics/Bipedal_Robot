"""Reposition body centres around anatomical hinge anchors in an isolated W2L copy."""
from pathlib import Path
from lxml import etree as E
import numpy as np
import math,json,copy,shutil,uuid,hashlib,sys
W=Path(sys.argv[2]) if len(sys.argv)>2 else Path(__file__).resolve().parents[1];SRC=Path(sys.argv[1]) if len(sys.argv)>1 else W.parent/'W2L_Treadmill_FrozenRoot';W.mkdir(exist_ok=True);P=E.XMLParser(strip_cdata=False)

def put(b,k,v):
 n=b.find(k)
 if n is None:n=E.SubElement(b,k)
 n.text=str(v);return n

def vec(b,k,gui):
 n=b.find(k)
 return np.array([float(n.find(a).get('Actual')) for a in 'XYZ']) if gui else np.array([float(n.get(a)) for a in 'xyz'])
def setvec(b,k,x,gui):
 n=b.find(k)
 if n is None:n=E.SubElement(b,k)
 if gui:
  for a,v in zip('XYZ',x):
   q=n.find(a)
   if q is None:q=E.SubElement(n,a)
   q.attrib.update({'Value':format(v,'.15g'),'Scale':'None','Actual':format(v,'.15g')})
 else:n.attrib.update({a:format(v,'.15g') for a,v in zip('xyz',x)})
def Rxyz(a):
 x,y,z=a;cx,sx,cy,sy,cz,sz=math.cos(x),math.sin(x),math.cos(y),math.sin(y),math.cos(z),math.sin(z)
 return np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])@np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])@np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
def rot(b,gui):return Rxyz(np.deg2rad(vec(b,'Rotation',gui)) if gui else vec(b,'Rotation',gui))
def transform(b,p,R,gui):
 ang=np.array([math.atan2(-R[1,2],R[2,2]),math.asin(max(-1,min(1,R[0,2]))),math.atan2(-R[0,1],R[0,0])])
 setvec(b,'LocalPosition' if gui else 'Position',p,gui);setvec(b,'Rotation',np.rad2deg(ang) if gui else ang,gui)
 if gui:
  m=np.eye(4);m[:3,:3]=R.T;m[3,:3]=p/.1;put(b,'LocalMatrix',','.join(format(v,'.15g') for v in m.ravel()))
def bodymap(r):return {b.findtext('Name'):b for b in r.findall('.//RigidBody') if b.findtext('Type')=='Box'}
def world(b,gui):
 p=vec(b,'LocalPosition' if gui else 'Position',gui);R=rot(b,gui);par=b.getparent()
 while par is not None and par.tag!='RigidBody':par=par.getparent()
 if par is not None:
  pp,RR=world(par,gui);return pp+RR@p,RR@R
 return p,R

def scalar(b,k,gui):return float(b.find(k).get('Actual') if gui else b.findtext(k))
def corners(b,gui):
 p,R=world(b,gui);size=np.array([scalar(b,'Length',gui),scalar(b,'Height',gui),scalar(b,'Width',gui)])
 return np.array([p+R@(size*np.array([x,y,z])/2) for x in [-1,1] for y in [-1,1] for z in [-1,1]])
report={'source':str(SRC),'poses_deg':{'L':{'hip':-10,'knee_flexion':10},'R':{'hip':20,'knee_flexion':55}},'models':{}}
for gui in [False,True]:
 src=next(SRC.glob('*.aproj')) if gui else (next(iter(SRC.glob('*Standalone.asim')),None) or next(SRC.glob('*.asim')));t=E.parse(str(src),P);r=t.getroot();bs=bodymap(r);root=bs['Root'];initial={};changed=[];oldrefs={};limit_changes=[]
 for leg in ['L','R']:
  fem=bs['femur_'+leg];tib=bs['tibia_'+leg];joint=tib.find('Joint');pf=vec(tib,'LocalPosition' if gui else 'Position',gui)+rot(tib,gui)@vec(joint,'LocalPosition' if gui else 'Position',gui);initial[leg]={'original_proximal_tibia_anchor_in_femur_m':pf.tolist(),'femur_distal_face_in_femur_m':[0,-.21,0],'original_offset_m':(pf-np.array([0,-.21,0])).tolist(),'original_proximal_tibia_box_face_in_femur_m':(vec(tib,'LocalPosition' if gui else 'Position',gui)+rot(tib,gui)@np.array([0,0,-.21])).tolist()}
 for leg in ['L','R']:
  rf=rot(bs['femur_'+leg],gui);rt=rf@rot(bs['tibia_'+leg],gui);rfp=rt@rot(bs['foot_'+leg],gui)
  f=rf@np.array([0,-1,0]);tibv=rt@np.array([0,0,1]);fv=rfp@np.array([1,0,0]);h=math.atan2(f[0],-f[1]);shank=math.atan2(tibv[0],-tibv[1]);pitch=math.atan2(fv[1],fv[0])
  oldrefs[leg]={'hip':h,'knee':h-shank,'ankle':pitch-shank}
 transform(root,vec(root,'LocalPosition' if gui else 'Position',gui),np.eye(3),gui)
 # Local box axes: femur length is Y, tibia length is Z, foot thickness is Z.
 for leg,h,k in [('L',-10,10),('R',20,55)]:
  fem=bs['femur_'+leg];tib=bs['tibia_'+leg];foot=bs['foot_'+leg]
  Rf=Rxyz([0,0,math.radians(h)]);hip=np.array([0,-.06,.06 if leg=='L' else -.06]);qh=np.array([0,.21,0]);transform(fem,hip-Rf@qh,Rf,gui);setvec(fem.find('Joint'),'LocalPosition' if gui else 'Position',qh,gui)
  # Match the proximal tibia face and knee hinge exactly to the distal femur face.
  Rt=Rxyz([0,0,math.radians(-k)])@Rxyz([math.pi/2,0,0]);qk=np.array([0,0,-.21]);pk=np.array([0,-.21,0]);transform(tib,pk-Rt@qk,Rt,gui);setvec(tib.find('Joint'),'LocalPosition' if gui else 'Position',qk,gui)
  # Keep the feet level in the pelvis frame while preserving each foot's ankle anchor.
  ankle_deg=0;pitch_deg=h-k+ankle_deg;Rfoot=Rt.T@Rf.T@Rxyz([0,0,math.radians(pitch_deg)])@Rxyz([math.pi/2,0,0]);qa=vec(foot.find('Joint'),'LocalPosition' if gui else 'Position',gui);pa=np.array([0,0,.21]);transform(foot,pa-Rfoot@qa,Rfoot,gui)
  for kind,body,newref in [('hip',fem,math.radians(h)),('knee',tib,math.radians(k)),('ankle',foot,math.radians(ankle_deg))]:
   j=body.find('Joint');offset=oldrefs[leg][kind]-newref;oldlimits=[];newlimits=[]
   for lim in ['LowerLimit','UpperLimit']:
    q=j.find(lim+'/LimitPos');oldv=float(q.get('Actual')) if gui else float(q.text);oldrad=math.radians(oldv) if gui else oldv;newrad=oldrad+offset
    if gui:q.attrib.update({'Value':format(math.degrees(newrad),'.15g'),'Scale':'None','Actual':format(math.degrees(newrad),'.15g')})
    else:q.text=format(newrad,'.15g')
    oldlimits.append(oldrad);newlimits.append(newrad)
   limit_changes.append(dict(joint=j.findtext('Name'),old_anatomical_zero_deg=math.degrees(oldrefs[leg][kind]),new_anatomical_zero_deg=math.degrees(newref),relative_limit_shift_deg=math.degrees(offset),old_relative_limits_deg=list(map(math.degrees,oldlimits)),new_relative_limits_deg=list(map(math.degrees,newlimits)),preserved_anatomical_limits_deg=[math.degrees(oldrefs[leg][kind]+v) for v in oldlimits]))
  changed.extend([fem.findtext('Name'),tib.findtext('Name'),foot.findtext('Name')])
 floor=float(corners(bs['Treadmill Belt'],gui)[:,1].max())
 soles=[float(corners(b,gui)[:,1].min()) for b in bs.values() if b.findtext('Name') in ['foot_L','foot_R','toe_L','toe_R']]
 p=vec(root,'LocalPosition' if gui else 'Position',gui);p[1]+=(floor+.01-min(soles));transform(root,p,np.eye(3),gui)
 org=root.getparent();setvec(org,'LocalPosition' if gui else 'Position',p,gui)
 if gui:
  transform(org,p,np.eye(3),True);put(r,'ProjectName','W2L_Geometry');put(r,'SimulationFile','W2L_Geometry.asim');put(r.find('Simulation'),'ProjectPath',str(W.resolve())+'\\')
 else:put(r,'ProjectPath',str(W.resolve())+'\\');put(r,'PlaybackControlMode','1')
 put(root,'Freeze','True')
 checks=[]
 for leg in ['L','R']:
  fem=bs['femur_'+leg];tib=bs['tibia_'+leg];foot=bs['foot_'+leg];pf,Rf=world(fem,gui);pt,Rt=world(tib,gui);pa,Ra=world(foot,gui)
  fd=pf+Rf@np.array([0,-.21,0]);tp=pt+Rt@np.array([0,0,-.21]);ka=pt+Rt@vec(tib.find('Joint'),'LocalPosition' if gui else 'Position',gui)
  hip=pf+Rf@vec(fem.find('Joint'),'LocalPosition' if gui else 'Position',gui);ankle=pa+Ra@vec(foot.find('Joint'),'LocalPosition' if gui else 'Position',gui)
  checks.append(dict(leg=leg,hip_anchor_world_m=hip.tolist(),distal_femur_world_m=fd.tolist(),proximal_tibia_world_m=tp.tolist(),knee_anchor_world_m=ka.tolist(),knee_face_gap_m=float(np.linalg.norm(tp-fd)),ankle_anchor_world_m=ankle.tolist(),foot_sole_y_m=float(corners(foot,gui)[:,1].min())))
 for d in r.iter('DiagramXml'):
  if d.text:E.fromstring(d.text.encode())
 dest=W/('W2L_Geometry.aproj' if gui else 'W2L_Geometry_Standalone.asim');t.write(str(dest),encoding='utf-8',xml_declaration=True)
 if not gui:shutil.copyfile(dest,W/'W2L_Geometry.asim')
 report['models']['aproj' if gui else 'asim']=dict(original=initial,checks=checks,root_position_m=p.tolist(),floor_top_y_m=floor,minimum_solid_foot_to_floor_gap_m=.01,foot_clearance_definition='minimum rotated corner of solid foot/toe boxes; contact sensor boxes excluded',hinge_limit_changes=limit_changes,diagram_xml_valid=True,changed_bodies=changed,neural_payload_sha256=hashlib.sha256(E.tostring(r.find('.//NervousSystem'))).hexdigest())
for f in SRC.glob('*.aform'):shutil.copyfile(f,W/f.name)
shutil.copyfile(W.parent/'W2L_Treadmill_FrozenRoot/Treadmill_checker.bmp',W/'Treadmill_checker.bmp')
(W/'GEOMETRY_REPORT.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report['models']['asim'],indent=2))



