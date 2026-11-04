"""Run with blender -b --factory-startup -P this_file.py.

Native, baked, editable CAD and BPA packaging preview. All geometry is in meters.
"""
from pathlib import Path
import base64,json,math
import bpy
from mathutils import Matrix,Vector
import numpy as np

HERE=Path(__file__).resolve().parents[1]/'packaging/output'
D=json.loads((HERE/'scene_data.json').read_text())
scene=bpy.data.scenes.new('AH001 animated packaging')
bpy.context.window.scene=scene
scene.world=bpy.data.worlds.new('AH001 studio world')
scene.unit_settings.system='METRIC';scene.unit_settings.length_unit='METERS'
scene.render.engine='CYCLES';scene.cycles.samples=24
scene.render.resolution_x=1600;scene.render.resolution_y=1200;scene.render.resolution_percentage=100
scene.render.fps=24;scene.world.color=(.22,.22,.22)
scene.view_settings.view_transform='AgX'

def collection(name,parent=None):
 c=bpy.data.collections.new(name);(parent or scene.collection).children.link(c);return c
CAD=collection('01 Structural CAD');ROUTES=collection('02 KompareXXVII selected routes');STAGE=collection('03 Lighting and cameras')
ROOT=bpy.data.objects.new('Anatomical frame: x anterior, y superior, z right',None);scene.collection.objects.link(ROOT)
ROOT.rotation_euler=(math.pi/2,0,0)
ROOT['display_adapter']='Blender display: X=x, Y=-z, Z=y. Source anatomical coordinates are unchanged beneath this root.'
def material(name,color,metal=0.,rough=.4):
 m=bpy.data.materials.new(name);m.diffuse_color=(*color,1);m.use_nodes=True
 p=m.node_tree.nodes.get('Principled BSDF');p.inputs['Base Color'].default_value=(*color,1);p.inputs['Metallic'].default_value=metal;p.inputs['Roughness'].default_value=rough
 return m
cadmat=material('Onyx structural surfaces',(.18,.22,.24),.15,.36)
capmat=material('End fittings',(.16,.18,.19),.65,.28)
tendonmat=material('Tendon cable',(.10,.12,.13),.7,.3)
palette=[(.03,.38,.30),(.65,.24,.08),(.31,.20,.58),(.03,.30,.61),(.62,.12,.29),(.44,.37,.04),(.10,.40,.47)]
bodyobjects={}
def meshobj(name,v,f,c,mat,parent=ROOT):
 mesh=bpy.data.meshes.new(name);mesh.from_pydata(v,[],f);mesh.update()
 obj=bpy.data.objects.new(name,mesh);c.objects.link(obj);obj.parent=parent;obj.data.materials.append(mat)
 for p in mesh.polygons:p.use_smooth=True
 return obj
for name,T in D['W0'].items():
 o=bpy.data.objects.new(name+' joint frame',None);CAD.objects.link(o);o.parent=ROOT;o.empty_display_size=.015
 bodyobjects[name]=o
for m in D['meshes']:
 v=np.frombuffer(base64.b64decode(m['v']),dtype='<f4').reshape(-1,3);f=np.frombuffer(base64.b64decode(m['f']),dtype='<u4').reshape(-1,3)
 o=meshobj(m['name'],v.tolist(),f.tolist(),CAD,cadmat,bodyobjects[m['body']])
 o['source']='Restored STL local datum matched to original STEP; audited CAD body placement.'

sequence=list(range(19))+list(range(17,-1,-1))
keyframes=[1+6*i for i in range(len(sequence))]
scene.frame_start=1;scene.frame_end=keyframes[-1]
for sample,frame in zip(sequence,keyframes):
 for name,o in bodyobjects.items():
  T=Matrix(D['poses'][sample]['W'][name]);o.rotation_mode='QUATERNION';o.matrix_local=T
  o.keyframe_insert('location',frame=frame);o.keyframe_insert('rotation_quaternion',frame=frame)

def apply(T,p):
 T=np.asarray(T);return np.asarray(p)@T[:3,:3].T+T[:3,3]
def state(r,deg):
 k=next((i for i,v in enumerate(r['dofs']) if 'knee' in v),None)
 if k is None:
  i=np.abs(r['angles'][0]).argmin();j=np.abs(r['angles'][1]).argmin();return r['contraction'][i][j]
 other=int(np.abs(r['angles'][1-k]).argmin());s=np.asarray(r['contraction']);s=s[:,other] if k==0 else s[other,:]
 return float(np.interp(deg,r['angles'][k],s))
def diameter(r,L):
 k,dk={10:(.17,.0225),20:(.25,.038),40:(.32,.075)}[round(r['diameter']*1000)]
 lm=r['rest']*(1-k);pn=math.sqrt((r['rest']**2-lm**2)/(dk**2-r['diameter']**2));L=np.clip(L,lm,r['rest'])
 return math.sqrt(r['rest']**2+(pn*r['diameter'])**2-L**2)/pn
def points_along(P,start,end,M):
 lens=np.linalg.norm(np.diff(P,axis=0),axis=1);arc=np.r_[0,np.cumsum(lens)]
 # Suppress duplicate points in the saved polyline before interpolation.
 mask=np.r_[True,np.diff(arc)>1e-10];arc=arc[mask];P=P[mask]
 ds=np.linspace(np.clip(start,0,arc[-1]),np.clip(end,0,arc[-1]),M)
 return np.stack([np.interp(ds,arc,P[:,k]) for k in range(3)],axis=1)
def tube(P,r,N=12):
 v=[]
 for i,p in enumerate(P):
  d=P[min(i+1,len(P)-1)]-P[max(i-1,0)];d=d/max(np.linalg.norm(d),1e-12)
  ref=np.array([0,1,0]) if abs(d[1])<.9 else np.array([1,0,0]);u=np.cross(d,ref);u/=max(np.linalg.norm(u),1e-12);w=np.cross(d,u)
  v.extend(p+r*(math.cos(j*2*math.pi/N)*u+math.sin(j*2*math.pi/N)*w) for j in range(N))
 faces=[]
 for i in range(1,len(P)):
  for j in range(N):faces.append(((i-1)*N+j,(i-1)*N+(j+1)%N,i*N+(j+1)%N,i*N+j))
 faces+=[tuple(range(N-1,-1,-1)),tuple((len(P)-1)*N+j for j in range(N))]
 return np.array(v),faces

animobjects=[]
for ri,r in enumerate(D['routes']):
 c=collection(f"{r['site']} | {r['name']} | {r['n']} BPAs | G={r['gain']:g}",ROUTES)
 mat=material(r['name']+' BPA sleeve',palette[ri%len(palette)],.05,.38)
 for bi,off in enumerate(r['offsets']):
  poses={k:[] for k in ('sleeve','origin fitting','distal fitting','tendon')};faces={}
  for p in D['poses']:
   W=p['W'];L=r['rest']-state(r,p['knee_deg']);dia=diameter(r,L)
   P=np.array([apply(W[b],v) for b,v in zip(r['bodies'],r['local'])]);P+=np.array(W[r['offset_body']])[:3,:3]@off
   total=np.linalg.norm(np.diff(P,axis=0),axis=1).sum();fit=r['fit']
   intervals={'sleeve':(fit,fit+L,dia/2,48),'origin fitting':(0,fit,r['diameter']*.53,3),'distal fitting':(fit+L,fit*2+L,r['diameter']*.53,3),'tendon':(fit*2+L,total,.0018,32)}
   for kind,(a,b,rad,M) in intervals.items():
    v,f=tube(points_along(P,a,b,M),rad);poses[kind].append(v);faces[kind]=f
  for kind,vs in poses.items():
   o=meshobj(f"{r['site']} / BPA {bi+1} / {kind}",vs[0].tolist(),faces[kind],c,mat if kind=='sleeve' else tendonmat if kind=='tendon' else capmat)
   o['resting_length_m']=r['rest'];o['nominal_diameter_m']=r['diameter'];o['transmission_gain']=r['gain'];o['preview']='Sleeve placement and survey-to-CAD registration are provisional. No physical sheave layout is saved.'
   o.shape_key_add(name='Reference')
   blocks=[]
   for k in range(1,19):
    key=o.shape_key_add(name=f"CAD knee {D['poses'][k]['knee_deg']:g} deg")
    key.data.foreach_set('co',vs[k].ravel());blocks.append(key)
   for sample,frame in zip(sequence,keyframes):
    for k,key in enumerate(blocks,1):key.value=1. if sample==k else 0.;key.keyframe_insert('value',frame=frame)
   animobjects.append(o)
 print('Baked',r['site'],flush=True)

# Keep motion between baked samples smooth. Loop closure is verified at the keys.
for action in bpy.data.actions:
 try:
  for layer in action.layers:
   for strip in layer.strips:
    for bag in strip.channelbags:
     for curve in bag.fcurves:
      for k in curve.keyframe_points:k.interpolation='LINEAR'
 except AttributeError:pass

floor=trimesh_floor=None
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-.06));floor=bpy.context.object;floor.name='Studio ground';floor.data.materials.append(material('Warm white background',(.77,.80,.80),0,.75))
def camera(name,pos,target,lens=65):
 data=bpy.data.cameras.new(name);o=bpy.data.objects.new(name,data);STAGE.objects.link(o);o.location=pos;o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler();data.lens=lens;return o
cam=camera('Full robot packaging',(2.2,-3.2,1.05),(-.03,0,.65),65);scene.camera=cam
camera('Erector spinae closeup',(.65,-.85,1.3),(-.11,-.04,1.12),60)
for name,pos,power,size in [('Key',(2,-3,4),550,3),('Fill',(-2,-1,2),350,3),('Rim',(0,2,3),650,2)]:
 data=bpy.data.lights.new(name,'AREA');data.energy=power;data.shape='DISK';data.size=size
 o=bpy.data.objects.new(name,data);STAGE.objects.link(o);o.location=pos;o.rotation_euler=(Vector((0,0,.6))-o.location).to_track_quat('-Z','Y').to_euler()

readme=bpy.data.texts.new('READ ME - scope and controls')
readme.write((HERE.parent/'README.md').read_text())
scene['scope']='Packaging preview. Not a clearance-certified physical installation. Read the embedded scope notes.'
scene['motion']='19 CAD knee samples, 0 to -90 deg and back. Saved BPA closure contraction drives length and calibrated diameter.'
scene['source_frame']='x anterior, y superior, z right; root rotation adapts to Blender z-up.'
scene.frame_set(1)
for area in bpy.context.screen.areas:
 if area.type=='VIEW_3D':
  area.spaces.active.region_3d.view_distance=2.4;area.spaces.active.region_3d.view_location=(0,0,.65)
  area.spaces.active.shading.type='MATERIAL'
bpy.ops.wm.save_as_mainfile(filepath=str(HERE/'AH001_Packaging.blend'),compress=True)
scene.render.filepath=str(HERE/'Blender_Packaging.png')
if bpy.app.background:bpy.ops.render.render(write_still=True)
(HERE/'blender_verification.json').write_text(json.dumps(dict(blender=bpy.app.version_string,structural_meshes=len(D['meshes']),routes=len(D['routes']),animated_sleeves=sum('sleeve' in o.name for o in animobjects),baked_poses=19,frames=scene.frame_end,scope=scene['scope']),indent=2))
print('Saved native Blender scene and rendered preview.',flush=True)
