"""Frame the complete feet in the saved diagnostic closeup camera."""
from pathlib import Path
import bpy
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/audit_20261105'
for branch in ('compliant','rigid'):
 camera=bpy.data.objects['Packaging closeup '+branch];target=Vector((0,-.08,.42))
 camera.rotation_euler=(target-camera.location).to_track_quat('-Z','Y').to_euler();camera.data.ortho_scale=1.50
s=bpy.data.scenes['REPAIRED COMPLIANT routes - candidate audit'];bpy.context.window.scene=s;s.frame_set(1)
bpy.ops.wm.save_as_mainfile(filepath=str(OUT/'Repaired_Routes_OpenSim_Matched.blend'),compress=True)
s.camera=bpy.data.objects['Packaging closeup compliant'];s.render.filepath=str(OUT/'Repaired_Packaging_Closeup.png');bpy.ops.render.render(write_still=True)
