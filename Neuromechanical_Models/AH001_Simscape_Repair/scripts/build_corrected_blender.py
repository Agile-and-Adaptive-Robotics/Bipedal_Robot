"""Build the author-corrected packaging scene, preserving earlier .blend files."""
from pathlib import Path
import json,base64,math,sys
import bpy,numpy as np
from mathutils import Matrix,Vector
SCRIPTS=Path(__file__).resolve().parent;sys.path.insert(0,str(SCRIPTS))
from packaging_geometry import parts,contraction,diameter
OUT=SCRIPTS.parent/'packaging/output';D=json.loads((OUT/'scene_data_corrected.json').read_text())
scene=bpy.data.scenes.new('AH001 corrected 20mm and Bowden')
bpy.context.window.scene=scene
scene.unit_settings.system='METRIC';scene.unit_settings.length_unit='METERS'
scene.world=bpy.data.worlds.new('AH001 white scene world');scene.world.color=(1,1,1);scene.world.use_nodes=True
nodes=scene.world.node_tree.nodes;nodes.clear();links=scene.world.node_tree.links
light=nodes.new('ShaderNodeLightPath');camera_bg=nodes.new('ShaderNodeBackground');camera_bg.inputs['Color'].default_value=(1,1,1,1);camera_bg.inputs['Strength'].default_value=1
ambient=nodes.new('ShaderNodeBackground');ambient.inputs['Color'].default_value=(1,1,1,1);ambient.inputs['Strength'].default_value=.3
mix=nodes.new('ShaderNodeMixShader');wo=nodes.new('ShaderNodeOutputWorld')
links.new(light.outputs['Is Camera Ray'],mix.inputs[0]);links.new(ambient.outputs[0],mix.inputs[1]);links.new(camera_bg.outputs[0],mix.inputs[2]);links.new(mix.outputs[0],wo.inputs['Surface'])
scene.render.engine='CYCLES';scene.cycles.samples=16;scene.render.threads_mode='FIXED';scene.render.threads=4
scene.render.resolution_x=1600;scene.render.resolution_y=1200;scene.render.resolution_percentage=100
scene.render.fps=24;scene.render.film_transparent=False
scene.view_settings.view_transform='Standard';scene.view_settings.look='Medium High Contrast' if 'Medium High Contrast' in [x.identifier for x in scene.view_settings.bl_rna.properties['look'].enum_items] else scene.view_settings.look

def coll(name,parent=None):
 c=bpy.data.collections.new(name);(parent or scene.collection).children.link(c);return c
CAD=coll('01 CAD and registered feet');ROUTES=coll('02 Individual 20mm BPA channels');STUDIO=coll('03 Camera orbit and lights')
root=bpy.data.objects.new('Source frame x anterior y superior z right',None);scene.collection.objects.link(root);root.rotation_euler=(math.pi/2,0,0)
root['display_adapter']='Blender world X=x, Y=-z, Z=y. Anatomical coordinates retained below this root.'
def material(name,color,metal=0,rough=.4):
 m=bpy.data.materials.new(name);m.diffuse_color=(*color,1);m.use_nodes=True;p=m.node_tree.nodes.get('Principled BSDF')
 p.inputs['Base Color'].default_value=(*color,1);p.inputs['Metallic'].default_value=metal;p.inputs['Roughness'].default_value=rough;return m
gray=material('Structural CAD',(.22,.27,.29),.10,.4);cap=material('20mm end fittings',(.17,.19,.20),.65,.3)
cable=material('Artificial tendon',(.07,.08,.09),.5,.4);bowden=material('Bowden housing',(.035,.045,.05),0,.65);pinmat=material('Pinned attachments',(.52,.38,.13),.7,.32)
palette=[(.03,.36,.29),(.66,.24,.10),(.35,.22,.55),(.05,.33,.57),(.60,.15,.31),(.47,.40,.06),(.12,.42,.47)]
def mesh(name,v,f,c,mat,parent=root):
 me=bpy.data.meshes.new(name);me.from_pydata(v,[],f);me.update();o=bpy.data.objects.new(name,me);c.objects.link(o);o.parent=parent;o.data.materials.append(mat)
 for p in me.polygons:p.use_smooth=True
 return o
N=len(D['poses']);seq=list(range(N))+list(range(N-2,-1,-1));frames=[1+6*i for i in range(len(seq))]
scene.frame_start=1;scene.frame_end=frames[-1]
bod={}
for name in D['W0']:
 o=bpy.data.objects.new(name+' corrected frame',None);CAD.objects.link(o);o.parent=root;o.rotation_mode='QUATERNION';o.empty_display_size=.01;bod[name]=o
for sample,frame in zip(seq,frames):
 for name,o in bod.items():
  o.matrix_local=Matrix(D['poses'][sample]['W'][name]);o.keyframe_insert('location',frame=frame);o.keyframe_insert('rotation_quaternion',frame=frame)
for m in D['meshes']:
 v=np.frombuffer(base64.b64decode(m['v']),dtype='<f4').reshape(-1,3);f=np.frombuffer(base64.b64decode(m['f']),dtype='<u4').reshape(-1,3)
 o=mesh(m['name'],v.tolist(),f.tolist(),CAD,gray,bod[m['body']]);o['source']='Restored CAD mesh; survey foot insertions share the registered foot frame.'

counts=[];calibration=[]
for ri,r in enumerate(D['routes']):
 c=coll(f"{r['site']} | {r['n']} parallel 20mm channels | G {r['gain']:g}",ROUTES)
 m=material(r['name']+' 20mm sleeves',palette[ri%len(palette)],.02,.43)
 samples=[]
 for p in D['poses']:
  objects,metrics=parts(r,p['W'],p['knee_deg']);samples.append(objects)
  calibration.append(dict(site=r['site'],knee_deg=p['knee_deg'],**metrics))
 conflict=any(x['straight_span_conflict'] for x in calibration if x['site']==r['site'])
 c['straight_span_conflict']=conflict
 c['placement_note']='Calf sleeves end at first foot guide; proximal leader shows any mount relocation. Survey lengths retained.' if any('calcn' in f for f in r['frames']) else 'Straight sleeves; routed tendons and split Bowden sections.'
 if conflict:c.name+=' | CHECK AVAILABLE LENGTH'
 for name,geo in samples[0].items():
  o=mesh(r['site']+' / '+name,geo['vertices'].tolist(),geo['faces'],c,{'sleeve':m,'fitting':cap,'cable':cable,'bowden':bowden,'pin':pinmat}[geo['kind']])
  o['nBPA_parallel']=r['n'];o['nominal_diameter_m']=.020;o['travel_gain']=r['gain'];o['component_kind']=geo['kind'];o['split_series_sections']=r['split_bowden'];o['saved_total_resting_length_m']=r['rest']
  o['provenance']='Author-directed display configuration. Original force/optimization results are unchanged.'
  o['check_available_straight_span']=conflict
  o.shape_key_add(name='Full extension')
  keys=[]
  for k in range(1,N):
   key=o.shape_key_add(name=f"Anatomical knee {D['poses'][k]['knee_deg']:g}");key.data.foreach_set('co',samples[k][name]['vertices'].ravel());keys.append(key)
  for sample,frame in zip(seq,frames):
   for k,key in enumerate(keys,1):key.value=float(k==sample);key.keyframe_insert('value',frame=frame)
 counts.append(dict(site=r['site'],parallel_channels=r['n'],sleeve_sections=r['n']*(2 if r['split_bowden'] else 1),diameter_mm=20,split_bowden=r['split_bowden']))
 print('Corrected',r['site'],flush=True)

target=Vector((-.03,0,.65));pivot=bpy.data.objects.new('Camera orbit - world Z - 180 degrees',None);STUDIO.objects.link(pivot);pivot.location=target;pivot.rotation_mode='XYZ'
camdata=bpy.data.cameras.new('Packaging orbit camera');cam=bpy.data.objects.new('Packaging orbit camera',camdata);STUDIO.objects.link(cam)
cam.location=(2.2,-3.2,1.05);cam.rotation_euler=(target-cam.location).to_track_quat('-Z','Y').to_euler();camdata.type='ORTHO';camdata.ortho_scale=2.2
world=cam.matrix_basis.copy();cam.parent=pivot;cam.matrix_parent_inverse=Matrix.Identity(4);cam.matrix_basis=Matrix.Translation(-target)@world;scene.camera=cam
pivot.rotation_euler.z=0;pivot.keyframe_insert('rotation_euler',index=2,frame=scene.frame_start)
pivot.rotation_euler.z=math.pi;pivot.keyframe_insert('rotation_euler',index=2,frame=scene.frame_end)
for name,pos,energy,size in [('Key',(2,-3,4),220,3),('Fill',(-2,-1,2),130,3),('Rim',(0,2,3),200,2)]:
 light=bpy.data.lights.new(name,'AREA');light.energy=energy;light.size=size;o=bpy.data.objects.new(name,light);STUDIO.objects.link(o);o.location=pos;o.rotation_euler=(target-o.location).to_track_quat('-Z','Y').to_euler()
# Limit interpolation edits to actions belonging to the newly built scene.
actions=set()
for o in scene.objects:
 if o.animation_data and o.animation_data.action:actions.add(o.animation_data.action)
 if o.type=='MESH' and o.data.shape_keys and o.data.shape_keys.animation_data:actions.add(o.data.shape_keys.animation_data.action)
for action in actions:
 if action is None:continue
 for layer in action.layers:
  for strip in layer.strips:
   for bag in strip.channelbags:
    for f in bag.fcurves:
     for k in f.keyframe_points:k.interpolation='LINEAR'

scene['scope']='Corrected packaging preview: all20, nBPA unchanged, two series sections for long knee routes, central Bowden, registered foot insertions.'
scene['author_directive_manifest']=json.dumps(D['correction_manifest'])
scene['knee_range_deg']='0 to -120, measured from femur and tibia axes; both sides; no hyperextension.'
scene['original_scene_preserved']='This builder creates a separate scene and saves to a new file.'
text=bpy.data.texts.new('CORRECTIONS AND LIMITATIONS');text.write(json.dumps(D['correction_manifest'],indent=2))

# Native evaluated cross-section measurements and camera-orbit validation.
verify=dict(blender=bpy.app.version_string,routes=len(counts),parallel_channels=sum(r['parallel_channels'] for r in counts),sleeve_sections=sum(r['sleeve_sections'] for r in counts),nominal_diameters_mm=[20],frames=scene.frame_end,configurations=counts,diameter_measurements=[],camera_samples=[])
test=next(o for o in scene.objects if o.name.startswith('vas_int_r / BPA 1/proximal sleeve'))
for frame in (scene.frame_start,frames[N-1],scene.frame_end):
 scene.frame_set(frame);bpy.context.view_layer.update();ob=test.evaluated_get(bpy.context.evaluated_depsgraph_get());ring=np.array([v.co[:] for v in ob.data.vertices[:16]]);diam=2*np.linalg.norm(ring-ring.mean(0),axis=1).mean()*1000
 verify['diameter_measurements'].append(dict(frame=frame,evaluated_diameter_mm=float(diam)))
 verify['camera_samples'].append(dict(frame=frame,orbit_deg=math.degrees(pivot.rotation_euler.z),camera_world_location=list(cam.matrix_world.translation)))
assert all(19.999<=r['evaluated_diameter_mm']<=38.001 for r in verify['diameter_measurements'])
assert abs(verify['camera_samples'][-1]['orbit_deg']-180)<1e-3
assert max(r['evaluated_diameter_mm'] for r in verify['diameter_measurements'])>min(r['evaluated_diameter_mm'] for r in verify['diameter_measurements'])+.01
scene.frame_set(1)
for area in bpy.context.screen.areas:
 if area.type=='VIEW_3D':
  space=area.spaces.active;space.shading.type='SOLID';space.shading.color_type='MATERIAL';space.shading.background_type='WORLD';space.overlay.show_overlays=False;space.show_gizmo=False
  space.region_3d.view_perspective='CAMERA'
scene.render.filepath=str(OUT/'Blender_Packaging_Corrected.png')
path=OUT/'AH001_Packaging_Corrected.blend'
bpy.ops.wm.save_as_mainfile(filepath=str(path),compress=True)
(OUT/'corrected_blender_verification.json').write_text(json.dumps(verify,indent=2));(OUT/'corrected_calibration_samples.json').write_text(json.dumps(calibration,indent=2))
if bpy.app.background:bpy.ops.render.render(write_still=True)
print('Saved',path,flush=True)
