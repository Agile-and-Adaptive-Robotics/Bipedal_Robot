from pathlib import Path
import bpy
s=bpy.context.scene;s.frame_set(145)
s.cycles.samples=8;s.render.resolution_percentage=65
# A fixed camera makes the full-flexion silhouette easier to compare.
p=s.objects['Camera orbit - world Z - 180 degrees']
if p.animation_data:p.animation_data.action=None
p.rotation_euler.z=0
s.render.filepath=str(Path(bpy.data.filepath).parent/'Blender_Packaging_Flexion_Check.png')
bpy.ops.render.render(write_still=True)
# This check does not save changes to the animation or presentation settings.
