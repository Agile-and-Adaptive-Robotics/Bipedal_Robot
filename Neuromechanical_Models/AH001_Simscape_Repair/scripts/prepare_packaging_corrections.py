"""Read-only source extraction plus explicit author-directed display corrections.

Does not alter survey MATs, CSVs, force predictions, or optimization results.
"""
from pathlib import Path
import json,base64,math,csv
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
import trimesh
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent;REP=HERE.parent;REPO=REP.parents[1]
OUT=REP/'packaging/output';RT=REPO/'Code/MuJoCo_SNS/robot_torso'
AUD=RT.parent/'spinal/reports_20261027_robot_control'
D=json.loads((OUT/'scene_data.json').read_text())
saved={r['site']:r for r in json.loads((REP/'packaging/build/saved_routes.json').read_text())}
def apply(T,p):
 T=np.asarray(T);return np.asarray(p)@T[:3,:3].T+T[:3,3]
def pose(R=np.eye(3),t=np.zeros(3)):
 T=np.eye(4);T[:3,:3]=R;T[:3,3]=t;return T
tree={}
def walk(e,parent=None):
 for b in e.findall('body'):
  q=np.fromstring(b.get('quat','1 0 0 0'),sep=' ');j=b.find('joint')
  tree[b.get('name')]=dict(parent=parent,T=pose(Rotation.from_quat(np.r_[q[1:],q[0]]).as_matrix(),np.fromstring(b.get('pos','0 0 0'),sep=' ')),axis=np.fromstring(j.get('axis','0 0 1'),sep=' ') if j is not None else np.array([0.,0,1]))
  walk(b,b.get('name'))
walk(ET.parse(RT/'ah001_neural.xml').getroot().find('worldbody'))
def fk(q):
 W={}
 for name,b in tree.items():
  T=b['T']@pose(Rotation.from_rotvec(b['axis']*q.get(name,0)).as_matrix())
  W[name]=(W[b['parent']] if b['parent'] else np.eye(4))@T
 return W
Wref=fk({});components=json.loads((AUD/'robot_mass_geometry.json').read_text())['components']
pc=next(c for c in components if '/02_01_PE_' in c['name']);CR=np.array(pc['transform'])[:9].reshape(3,3)
axes={};anchors={};neutral={}
for side in ('r','l'):
 def long_axis(token):
  c=next(c for c in components if token in c['name']);u=(CR@np.array(c['transform'])[:9].reshape(3,3).T)[:,1]
  return u if u[1]<0 else -u
 f=long_axis('03_01_FB_R' if side=='r' else '03_03_FB_L');t=long_axis('05_01_TI_R' if side=='r' else '05_02_TI_')
 axes[side]=(f,t);neutral[side]=math.atan2(np.cross(f,t)[2],f[:2]@t[:2])
for c in json.loads((RT/'ah001_loop_closures.json').read_text()):
 a=np.array(c['anchor_body1_m']);anchors[c['body1'][-1]]=(a,apply(np.linalg.inv(Wref[c['body2']]),apply(Wref[c['body1']],a)))

# The original feet were level in the already-flexed CAD pose. Register the
# ankle orientation so they are level at anatomical knee zero instead.
for side in ('r','l'):
 tree[f'foot_{side}']['T'][:3,:3]=Rotation.from_rotvec([0,0,neutral[side]]).as_matrix()@tree[f'foot_{side}']['T'][:3,:3]

poses=[];previous={s:np.zeros(3) for s in ('r','l')};checks=[]
for target in np.linspace(0,-120,25):
 q={};actual={}
 for side in ('r','l'):
  a,b=anchors[side];uf,ut=axes[side]
  def residual(x):
   W=fk({f'back_{side}':x[0],f'shank_{side}':x[1],f'front_{side}':x[2]})
   dp=apply(W[f'front_{side}'],a)-apply(W[f'shank_{side}'],b)
   delta=tree[f'back_{side}']['axis'][2]*x[0]+tree[f'shank_{side}']['axis'][2]*x[1]
   return np.r_[dp,delta+neutral[side]-math.radians(target)]
  sol=least_squares(residual,previous[side],xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=600)
  previous[side]=sol.x;err=np.linalg.norm(residual(sol.x)[:3])
  if err>1e-6:raise RuntimeError(f'Four-bar cannot reach {side} {target}: residual {err}')
  q.update({f'back_{side}':sol.x[0],f'shank_{side}':sol.x[1],f'front_{side}':sol.x[2]})
  W=fk(q);u=W[f'thigh_{side}'][:3,:3]@uf;v=W[f'shank_{side}'][:3,:3]@ut
  angle=np.degrees(np.arctan2(np.cross(u,v)[2],u[:2]@v[:2]));actual[side]=float(angle)
  checks.append(dict(side=side,requested_deg=float(target),measured_bone_angle_deg=float(angle),closure_error_m=float(err)))
  assert abs(angle-target)<1e-6 and -120.000001<=angle<=.000001
 poses.append(dict(knee_deg=float(target),anatomical_knee_deg=actual,W={k:v.tolist() for k,v in fk(q).items()}))
D['poses']=poses

# Register all calcn/toe survey points through one coherent foot frame. The
# old per-route root chains placed the same Achilles point 135 mm apart.
rs={r['site']:r for r in D['routes']}
soleus=rs['soleus_r'];flex=rs['flex_dig_r']
Fcalcn=Wref[soleus['offset_body']]@np.array(soleus['neutral_chain'])[2]
Ftoe=Wref[flex['offset_body']]@np.array(flex['neutral_chain'])[2]
calcn_to_toe=np.linalg.inv(Fcalcn)@Ftoe
footmeshes=[m for m in D['meshes'] if m['body']=='foot_r']
toe_tri=[]
all_tri=[]
for m in footmeshes:
 v=np.frombuffer(base64.b64decode(m['v']),dtype='<f4').reshape(-1,3).astype(float);f=np.frombuffer(base64.b64decode(m['f']),dtype='<u4').reshape(-1,3)
 tri=v[f];all_tri.append(tri)
 if '/Toe-' in m['name']:toe_tri.append(tri)
all_tri=np.concatenate(all_tri);toe_tri=np.concatenate(toe_tri)
def source_calcn_point(r,i):
 sp=saved[r['site']]['spec'];i=i%len(sp['points']);point=np.array(sp['points'][i]);cross=sp['cross'];row=i+1
 fi=0 if row<cross[0] else 2 if cross[0]==cross[1] or row>=cross[1] else 1
 fr=sp['frames'][fi]
 if fr=='toes':return apply(calcn_to_toe,point)
 if fr=='calcnA':
  ankle=Wref[r['offset_body']]@np.array(r['neutral_chain'])[1]
  return apply(np.linalg.inv(Fcalcn)@ankle,point)
 return point
toe_routes=[rs[s] for s in ('flex_dig_r','flex_hal_r','ext_dig_r','ext_hal_r')]
src=[source_calcn_point(r,-1) for r in toe_routes]
# Translation-only anatomical registration, using distal toe extent and
# dorsal/plantar surfaces. No scaling or arbitrary rotation of the anatomy.
dx=float(toe_tri[:,:,0].max()-np.mean([v[0] for v in src]))
dy=float(np.mean([(toe_tri[:,:,1].min() if r['site'].startswith('flex') else toe_tri[:,:,1].max())-v[1] for r,v in zip(toe_routes,src)]))
F=pose(t=np.array([dx,dy,0.]));foot_log=[]
for r in D['routes']:
 for i,b in enumerate(r['bodies']):
  if b!='foot_r':continue
  point=source_calcn_point(r,i);nominal=apply(F,point);registered=nominal.copy()
  # Explicit attachment registration to the CAD surface for the terminal pin.
  # Retain the original source coordinate and snap distance in the audit.
  if i==len(r['local'])-1:
   tri=toe_tri if saved[r['site']]['spec']['frames'][-1]=='toes' else all_tri
   closest=trimesh.triangles.closest_point(tri,np.repeat(nominal[None,:],len(tri),axis=0))
   registered=closest[np.argmin(np.linalg.norm(closest-nominal,axis=1))]
  foot_log.append(dict(site=r['site'],point_index=i+1,source_calcn_m=point.tolist(),registered_foot_local_m=registered.tolist(),terminal=i==len(r['local'])-1,surface_adjustment_mm=float(np.linalg.norm(registered-nominal)*1000)))
  r['local'][i]=registered.tolist()
 for i in range(len(r['local'])):
  assert np.isfinite(r['local'][i]).all()

# Register calf survey coordinates between their knee datum and the same ankle
# datum used by the foot. The survey tibia is 430 mm; CAD has a different length.
# Keep a record of this display-only affine registration rather than allowing
# guides from that survey to extend below the CAD ankle.
ankle_in_calcn=apply(np.linalg.inv(Fcalcn)@(Wref[soleus['offset_body']]@np.array(soleus['neutral_chain'])[1]),[0,0,0])
ankle_local=apply(np.linalg.inv(Wref['shank_r']),apply(Wref['foot_r'],apply(F,ankle_in_calcn)))
calf_affine=np.eye(3);calf_affine[:,1]=ankle_local/(-.430)
for r in D['routes']:
 if not any('calcn' in f for f in r['frames']):continue
 for i,b in enumerate(r['bodies']):
  if b=='shank_r':r['local'][i]=(calf_affine@np.array(r['local'][i])).tolist()
D['calf_display_registration']=dict(matrix=calf_affine.tolist(),source_ankle_y_m=-.430,target_ankle_local_m=ankle_local.tolist())

configuration=[]
for r in D['routes']:
 old=r['diameter'];oldcount=r['n'];r['source_diameter_m']=old;r['diameter']=.020
 r['minimum']=r['rest']*.75
 # Keep the saved offset direction and layout, but use the inflated 20 mm OD
 # plus the existing 5 mm bundle margin to establish centre spacing.
 offsets=np.array(r['offsets']);axis=offsets[-1]-offsets[0] if len(offsets)>1 else np.array([0.,0,1.])
 axis=axis/max(np.linalg.norm(axis),1e-12)
 if len(offsets)>1:r['offsets']=(((np.arange(oldcount)-(oldcount-1)/2)*.043)[:,None]*axis).tolist()
 assert len(r['offsets'])==oldcount
 r['split_bowden']=bool(any('knee' in d for d in r['dofs']) and r['rest']>.50)
 r['series_length_fractions']=[.5,.5] if r['split_bowden'] else [1.]
 configuration.append(dict(site=r['site'],nBPA=oldcount,nominal_diameter_mm=20,travel_gain=r['gain'],source_diameter_mm=old*1000,split_bowden=r['split_bowden'],rest_length_total_mm=r['rest']*1000))
D['correction_manifest']=dict(source='Author directives in this session; original survey results retained.',all_replacements_20mm=True,count_rule='nBPA is unchanged; never multiply it by four or by G.',foot_registration=dict(calcn_to_CADfoot=F.tolist(),calcn_to_toes=calcn_to_toe.tolist(),method='Anatomical translation using toe extent and dorsal/plantar surfaces; terminal pins registered on CAD surfaces. Preview registration, not new measured hardware data.'),bowden_rule='For knee-crossing routes with saved resting BPA length > 0.50 m: two series sections, equal halves of the total saved length, and central Bowden tendon.',knee_rule='Physical femur/tibia sagittal axes; exact four-bar closure; both sides 0 to -120 deg.',camera_rule='180 deg orbit about Blender world Z over full animation.',background='Scene world and viewport white; no global user preference change.')
(OUT/'scene_data_corrected.json').write_text(json.dumps(D,separators=(',',':')))
for name,rows in [('corrected_configurations.csv',configuration),('corrected_knee_verification.csv',checks),('foot_attachment_registration.csv',foot_log)]:
 with (OUT/name).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print(json.dumps(dict(routes=len(D['routes']),parallel_counts=sum(r['n'] for r in D['routes']),all_diameters_mm=sorted({r['diameter']*1000 for r in D['routes']}),knee_range=[min(x['measured_bone_angle_deg'] for x in checks),max(x['measured_bone_angle_deg'] for x in checks)],max_loop_error_m=max(x['closure_error_m'] for x in checks),foot_transform_translation_m=[dx,dy,0],split_routes=[r['site'] for r in D['routes'] if r['split_bowden']]),indent=2))
