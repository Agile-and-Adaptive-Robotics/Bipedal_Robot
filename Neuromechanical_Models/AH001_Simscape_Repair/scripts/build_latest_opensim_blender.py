"""Build Blender from API-verified OpenSim coordinates, without CAD remapping."""
from pathlib import Path
import json,base64,math,sys,colorsys
import bpy,numpy as np
from mathutils import Matrix,Vector
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/audit_20261105'
REPAIRED='--repaired' in sys.argv
PREFIX='Repaired' if REPAIRED else 'Latest'
D=json.loads((OUT/('repaired_verified_scene.json' if REPAIRED else 'opensim_verified_scene.json')).read_text())
def material(name,color):
 m=bpy.data.materials.new(name);m.diffuse_color=(*color,1);m.use_nodes=True
 bs=m.node_tree.nodes.get('Principled BSDF');bs.inputs['Base Color'].default_value=(*color,1);bs.inputs['Roughness'].default_value=.6;return m
bone=material('Original OpenSim geometry',(.5,.52,.53));blue=material('Changed right-side paths BLUE',(0,.12,.85));gray=material('Original robotbody reference paths',(.32,.32,.32))
def tube(P,r,N=8):
 P=np.array(P);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))];keep=np.r_[True,np.diff(arc)>1e-10];P=P[keep];arc=arc[keep]
 ds=np.linspace(0,arc[-1],64);P=np.stack([np.interp(ds,arc,P[:,k]) for k in range(3)],axis=1)
 vv=[];ff=[]
 for i,p in enumerate(P):
  u=P[min(i+1,63)]-P[max(0,i-1)];u/=max(np.linalg.norm(u),1e-12)
  ref=np.array([0,1,0]) if abs(u[1])<.9 else np.array([1,0,0]);v=np.cross(u,ref);v/=max(np.linalg.norm(v),1e-12);w=np.cross(u,v)
  vv.extend(p+r*(np.cos(j*2*np.pi/N)*v+np.sin(j*2*np.pi/N)*w) for j in range(N))
 for i in range(63):
  for j in range(N):ff.append((i*N+j,i*N+(j+1)%N,(i+1)*N+(j+1)%N,(i+1)*N+j))
 return np.array(vv),ff
def mesh(name,v,f,col,mat,parent):
 me=bpy.data.meshes.new(name);me.from_pydata(v,[],f);me.update();ob=bpy.data.objects.new(name,me);col.objects.link(ob);ob.parent=parent;ob.data.materials.append(mat)
 for p in me.polygons:p.use_smooth=True
 return ob
def collection(scene,name):
 c=bpy.data.collections.new(name);scene.collection.children.link(c);return c
N=len(D['poses']);sequence=list(range(N))+list(range(N-2,-1,-1));frames=[1+6*i for i in range(len(sequence))]
for branch in ('compliant','rigid'):
 s=bpy.data.scenes.new(PREFIX.upper()+' '+branch.upper()+' routes - candidate audit');bpy.context.window.scene=s
 s.unit_settings.system='METRIC';s.frame_start=1;s.frame_end=frames[-1];s.render.fps=24
 s.world=bpy.data.worlds.new('White '+branch);s.world.use_nodes=True;s.world.color=(1,1,1);s.world.node_tree.nodes['Background'].inputs[0].default_value=(1,1,1,1)
 s.world.node_tree.nodes['Background'].inputs[1].default_value=.7
 s.render.film_transparent=True
 tree=bpy.data.node_groups.new('White composite '+branch,'CompositorNodeTree');s.compositing_node_group=tree
 tree.interface.new_socket(name='Image',in_out='OUTPUT',socket_type='NodeSocketColor')
 layers=tree.nodes.new('CompositorNodeRLayers');white=tree.nodes.new('CompositorNodeRGB');white.outputs[0].default_value=(1,1,1,1)
 over=tree.nodes.new('CompositorNodeAlphaOver');output=tree.nodes.new('NodeGroupOutput')
 tree.links.new(white.outputs[0],over.inputs['Background']);tree.links.new(layers.outputs['Image'],over.inputs['Foreground']);tree.links.new(over.outputs[0],output.inputs[0])
 s.render.engine='CYCLES';s.cycles.samples=16;s.render.threads_mode='FIXED';s.render.threads=4
 s.render.resolution_x=1600;s.render.resolution_y=1200;s.render.resolution_percentage=100;s.view_settings.view_transform='Standard'
 try:s.view_settings.look='None'
 except TypeError:pass
 root=bpy.data.objects.new(branch+' OpenSim Y-up to Blender Z-up',None);s.collection.objects.link(root);root.rotation_euler.x=math.pi/2
 structural=collection(s,'01 ORIGINAL OpenSim geometry and body frames');original=collection(s,'02 ORIGINAL right-side paths GRAY');changed=collection(s,'03 OPTIMIZED candidate paths BLUE')
 body={}
 for name in D['poses'][0]['W']:
  ob=bpy.data.objects.new(branch+' '+name,None);structural.objects.link(ob);ob.parent=root;ob.rotation_mode='QUATERNION';body[name]=ob
  for sample,frame in zip(sequence,frames):
   ob.matrix_local=Matrix(D['poses'][sample]['W'][name]);ob.keyframe_insert('location',frame=frame);ob.keyframe_insert('rotation_quaternion',frame=frame)
 for m in D['meshes']:
  v=np.frombuffer(base64.b64decode(m['vertices']),dtype='<f4').reshape(-1,3);f=np.frombuffer(base64.b64decode(m['faces']),dtype='<u4').reshape(-1,3)
  mesh(branch+' '+m['name'],v.tolist(),f.tolist(),structural,bone,body[m['body']])
 def animated_route(name,paths,col,mat,radius):
  samples=[tube(p,radius)[0] for p in paths];_,faces=tube(paths[0],radius)
  ob=mesh(name,samples[0].tolist(),faces,col,mat,root);ob.shape_key_add(name='Knee 0')
  keys=[]
  for k in range(1,N):
   key=ob.shape_key_add(name=f'Knee {-5*k}');key.data.foreach_set('co',samples[k].ravel());keys.append(key)
  for sample,frame in zip(sequence,frames):
   for k,key in enumerate(keys,1):key.value=float(k==sample);key.keyframe_insert('value',frame=frame)
  return ob
 for r in D['baselineRoutes']:animated_route('REFERENCE '+r['name'],r['poses'],original,gray,.0007)
 for r in D['routes']:
  if r['branch']!=branch:continue
  paths=[]
  for pose in D['poses']:
   paths.append([(np.array(pose['W'][b])[:3,:3]@p+np.array(pose['W'][b])[:3,3]).tolist() for b,p in zip(r['assignedBodies'],r['recordedLoc'])])
  ob=animated_route('OPT '+r['group']+' | '+r['status'],paths,changed,blue,.0017)
  for key in ('nBPA','gain','family','nominal_diameter_mm','fractionMet','recordedViolation','status','seedSignatureMatches'):ob[key]=r[key]
  ob['scope']='Recorded primary-pass path. Mechanical packaging is not fabricated from a route polyline.'
 if REPAIRED:
  packaging=collection(s,'04 SOLVED full-length BPA packaging CANDIDATES')
  fittings=collection(s,'05 FITTING and tackle ENVELOPES - reeving unresolved')
  fitmat=material('Rigid envelope '+branch,(.24,.27,.28))
  def cylinder_points(cap,NV=12):
   a=np.asarray(cap['p1']);b=np.asarray(cap['p2']);u=b-a;u=u/max(np.linalg.norm(u),1e-12)
   ref=np.array([0,1,0]) if abs(u[1])<.9 else np.array([1,0,0]);v=np.cross(u,ref);v=v/max(np.linalg.norm(v),1e-12);w=np.cross(u,v)
   r=cap['r'];vv=[p+r*(math.cos(j*2*math.pi/NV)*v+math.sin(j*2*math.pi/NV)*w) for p in (a,b) for j in range(NV)]
   faces=[tuple(range(NV-1,-1,-1)),tuple(range(NV,2*NV))]+[(j,(j+1)%NV,(j+1)%NV+NV,j+NV) for j in range(NV)]
   return np.asarray(vv),faces
  groups=[r for r in D['routes'] if r['branch']==branch]
  for gi,r in enumerate(groups):
   mat=material(r['group']+' solved BPA '+branch,colorsys.hsv_to_rgb(gi/len(groups),.55,.62))
   samples=r['sleeves'];samples=[x if isinstance(x,list) else [x] for x in samples]
   assert len(samples)==N and len({len(x) for x in samples})==1
   for ci in range(len(samples[0])):
    cap=samples[0][ci];vertices,faces=cylinder_points(cap)
    name=f"{r['group']} {r['nominal_diameter_mm']}mm channel{cap['branch']} section{cap['section']} | {r['status']}"
    ob=mesh(name,vertices.tolist(),faces,packaging,mat,root);ob.shape_key_add(name='Knee 0');keys=[]
    for k in range(1,N):
     points,_=cylinder_points(samples[k][ci]);key=ob.shape_key_add(name=f'Knee {-5*k}');key.data.foreach_set('co',points.ravel());keys.append(key)
    for sample,frame in zip(sequence,frames):
     for k,key in enumerate(keys,1):key.value=float(k==sample);key.keyframe_insert('value',frame=frame)
    for key in ('nBPA','gain','nominal_diameter_mm','fractionMet','recordedViolation','status'):ob[key]=r[key]
    ob['routing_group']=r['group'];ob['channel_index']=cap['branch'];ob['section_index']=cap['section']
    ob['length_m_by_knee_pose']=r['kneeLength_m'];ob['contraction_by_knee_pose']=r['kneeContraction']
    ob['diameter_m_by_knee_pose']=[2*x[ci]['r'] for x in samples]
    ob['scope']='Solved full pressure-vessel length; route-tangent pinned placement candidate. Inspect collisions, not approved hardware.'
   # End fitting spheres are conservative envelopes, not invented catalog parts.
   for fi in range(len(r['fittings'][0]['r'])):
    f0=r['fittings'][0];rr=f0['r'][fi];center=np.asarray(f0['p'][fi]);
    cap={'p1':center+np.array([0,-rr,0]),'p2':center+np.array([0,rr,0]),'r':rr}
    # Use low-poly UV spheres and animate their centers; tackle stays wireframe.
    bpy.ops.mesh.primitive_uv_sphere_add(segments=12,ring_count=6,radius=rr)
    ob=bpy.context.object;ob.name=r['group']+' '+f0['kind'][fi]+' ENVELOPE '+branch
    for c in list(ob.users_collection):c.objects.unlink(ob)
    fittings.objects.link(ob);ob.parent=root;ob.data.materials.append(fitmat)
    if f0['kind'][fi]=='rigid_tackle':ob.display_type='WIRE';ob.hide_render=True
    for sample,frame in zip(sequence,frames):ob.location=r['fittings'][sample]['p'][fi];ob.keyframe_insert('location',frame=frame)
 # Fixed camera first, with optional orbit across playback.
 target=Vector((-.02,0,.85));pivot=bpy.data.objects.new('180 degree camera orbit '+branch,None);s.collection.objects.link(pivot);pivot.location=target
 cd=bpy.data.cameras.new('Camera '+branch);cam=bpy.data.objects.new('Camera '+branch,cd);s.collection.objects.link(cam);cam.location=(2.2,-3.3,1.15);cam.rotation_euler=(target-cam.location).to_track_quat('-Z','Y').to_euler();cd.type='ORTHO';cd.ortho_scale=2.05
 cd.ortho_scale=2.8
 world=cam.matrix_basis.copy();cam.parent=pivot;cam.matrix_basis=Matrix.Translation(-target)@world;s.camera=cam
 for frame,angle in [(1,0),(frames[-1],math.pi)]:pivot.rotation_euler.z=angle;pivot.keyframe_insert('rotation_euler',index=2,frame=frame)
 if REPAIRED:
  close_data=bpy.data.cameras.new('Packaging closeup '+branch);close=bpy.data.objects.new('Packaging closeup '+branch,close_data);s.collection.objects.link(close)
  close.location=(1.35,-2.0,.68);close_target=Vector((0,-.08,.42));close.rotation_euler=(close_target-close.location).to_track_quat('-Z','Y').to_euler();close_data.type='ORTHO';close_data.ortho_scale=1.50
 for name,pos,energy in [('Key',(2,-3,4),180),('Fill',(-2,1,3),90)]:
  ld=bpy.data.lights.new(name+branch,'AREA');ld.energy=energy;ld.size=3;lo=bpy.data.objects.new(name+branch,ld);s.collection.objects.link(lo);lo.location=pos;lo.rotation_euler=(target-lo.location).to_track_quat('-Z','Y').to_euler()
 s['scope']='Verified OpenSim route comparison, not a certified packaged robot.';s['limitations']=json.dumps(D['limitations']);s.frame_set(1)
for a in bpy.data.actions:
 for layer in a.layers:
  for strip in layer.strips:
   for bag in strip.channelbags:
    for f in bag.fcurves:
     for k in f.keyframe_points:k.interpolation='LINEAR'
bpy.context.window.scene=bpy.data.scenes[PREFIX.upper()+' COMPLIANT routes - candidate audit']
for startup_scene in list(bpy.data.scenes):
 if not startup_scene.name.startswith(PREFIX.upper()+' '):bpy.data.scenes.remove(startup_scene)
for unused in list(bpy.data.objects):
 if not unused.users_scene:bpy.data.objects.remove(unused)
for area in bpy.context.screen.areas:
 if area.type=='VIEW_3D':
  area.spaces.active.region_3d.view_perspective='CAMERA';area.spaces.active.shading.color_type='MATERIAL'
  area.spaces.active.shading.background_type='WORLD'
bpy.data.texts.new('READ ME - RESULTS ARE NOT ALL FEASIBLE').write('\n'.join(D['limitations'])+'\nGray: original robotbody paths. Blue: saved optimized primary-pass paths. Inspect object properties for actual count, gain, status and meet fraction.')
path=OUT/(PREFIX+'_Routes_OpenSim_Matched.blend');bpy.ops.wm.save_as_mainfile(filepath=str(path),compress=True)
s=bpy.context.scene;s.render.filepath=str(OUT/(PREFIX+'_Routes_OpenSim_Matched.png'));bpy.ops.render.render(write_still=True)
if REPAIRED:
 original_camera=s.camera;s.camera=bpy.data.objects['Packaging closeup compliant'];s.render.filepath=str(OUT/'Repaired_Packaging_Closeup.png');bpy.ops.render.render(write_still=True);s.camera=original_camera
print('Saved',path)
