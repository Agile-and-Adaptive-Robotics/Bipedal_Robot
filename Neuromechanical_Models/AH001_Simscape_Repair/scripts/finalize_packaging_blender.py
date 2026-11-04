from pathlib import Path
import json,bpy,numpy as np
here=Path(__file__).resolve().parents[1]/'packaging/output'
s=bpy.context.scene
o=next(o for o in s.objects if o.name.startswith('vas_int_r / BPA 1 / sleeve'))
s.frame_set(1);bpy.context.view_layer.update();a=np.array([v.co[:] for v in o.evaluated_get(bpy.context.evaluated_depsgraph_get()).data.vertices])
s.frame_set(109);bpy.context.view_layer.update();b=np.array([v.co[:] for v in o.evaluated_get(bpy.context.evaluated_depsgraph_get()).data.vertices])
assert np.max(np.linalg.norm(a-b,axis=1))>.001
qa=json.loads((here/'blender_verification.json').read_text());qa['evaluated_vastus_vertex_motion_m']=float(np.max(np.linalg.norm(a-b,axis=1)));qa['feet_basis']='Foot assembly axes mapped through the recorded distal mount, matching foot_geometry.py.'
(here/'blender_verification.json').write_text(json.dumps(qa,indent=2))
s.frame_set(1)
floor=s.objects.get('Studio ground')
if floor:floor.hide_set(True)
for obj in s.objects:obj.select_set(False)
bpy.context.view_layer.objects.active=o
area=bpy.context.area;area.type='VIEW_3D';v=area.spaces.active
v.shading.type='SOLID';v.shading.color_type='MATERIAL';v.overlay.show_overlays=False;v.show_gizmo=False
v.region_3d.view_location=(0,0,.65);v.region_3d.view_rotation=s.camera.matrix_world.to_quaternion();v.region_3d.view_distance=2.4;v.region_3d.view_perspective='PERSP'
s.render.threads_mode='FIXED';s.render.threads=4;s.cycles.samples=16
s.render.resolution_x=1200;s.render.resolution_y=900
s.render.filepath=str(here/'Blender_Packaging.png')
bpy.ops.wm.save_as_mainfile(filepath=str(here/'AH001_Packaging.blend'),compress=True)
bpy.ops.render.render(write_still=True)
