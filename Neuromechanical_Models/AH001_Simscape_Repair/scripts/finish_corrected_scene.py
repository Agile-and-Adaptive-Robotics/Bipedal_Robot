"""Set a flat white rendered background and final camera framing."""
from pathlib import Path
import bpy,json,sys,numpy as np
s=bpy.context.scene
s.view_settings.view_transform='Standard'
try:s.view_settings.look='None'
except TypeError:pass
s.view_settings.exposure=0;s.view_settings.gamma=1;s.render.dither_intensity=0
# Composite the transparent environment over constant white. The actual world
# remains white as well, and the viewport uses its white world colour.
s.render.film_transparent=True
tree=bpy.data.node_groups.new('AH001 white background','CompositorNodeTree')
s.compositing_node_group=tree
tree.interface.new_socket(name='Image',in_out='OUTPUT',socket_type='NodeSocketColor')
nodes=tree.nodes;links=tree.links
layers=nodes.new('CompositorNodeRLayers');white=nodes.new('CompositorNodeRGB');white.outputs[0].default_value=(1,1,1,1)
over=nodes.new('CompositorNodeAlphaOver');over.inputs['Factor'].default_value=1
out=nodes.new('NodeGroupOutput');links.new(white.outputs[0],over.inputs['Background']);links.new(layers.outputs['Image'],over.inputs['Foreground']);links.new(over.outputs[0],out.inputs[0])
s.frame_set(1);s.camera.data.ortho_scale=2.2
# Check actual evaluated Blender sleeve vertices, independently of the
# shape-key authoring loop, after reopening the saved file.
sys.path.insert(0,str(Path(__file__).parent))
from packaging_geometry import contraction
d=json.loads((Path(bpy.data.filepath).parent/'scene_data_corrected.json').read_text())
report=[]
for frame,angle in [(1,0),(73,-60),(145,-120)]:
 s.frame_set(frame);bpy.context.view_layer.update();dg=bpy.context.evaluated_depsgraph_get()
 for r in d['routes']:
  for o in s.objects:
   if not o.name.startswith(r['site']+' / ') or o.get('component_kind')!='sleeve':continue
   ob=o.evaluated_get(dg);v=np.array([x.co[:] for x in ob.data.vertices]).reshape(-1,16,3)
   centers=v.mean(axis=1);actual=float(np.linalg.norm(centers[-1]-centers[0]))
   expected=(r['rest']-contraction(r,angle))/(2 if r['split_bowden'] else 1)
   assert abs(actual-expected)<1e-5,(o.name,frame,actual,expected)
   assert np.linalg.norm(np.cross(centers-centers[0],(centers[-1]-centers[0])/actual),axis=1).max()<1e-5
   report.append(dict(object=o.name,frame=frame,length_m=actual,diameter_m=float(2*np.linalg.norm(v[0]-centers[0],axis=1).mean())))
(Path(bpy.data.filepath).parent/'evaluated_straight_sleeves.json').write_text(json.dumps(report,indent=2))
s.frame_set(1)
bpy.ops.wm.save_as_mainfile(filepath=bpy.data.filepath,compress=True)
s.render.filepath=str(Path(bpy.data.filepath).parent/'Blender_Packaging_Corrected.png');bpy.ops.render.render(write_still=True)
