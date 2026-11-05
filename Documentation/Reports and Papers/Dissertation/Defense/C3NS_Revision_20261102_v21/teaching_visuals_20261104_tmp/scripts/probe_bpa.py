"""probe_bpa.py — QA utility: render the saved v1 blend with the cable and
all text hidden (isolates the force path). Writes into the repo build dir."""
import bpy
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
bpy.data.objects["cable"].hide_render = True
for ob in bpy.data.objects:
    if ob.type == "FONT" or ob.name.startswith("lead_"):
        ob.hide_render = True
sc = bpy.context.scene
sc.frame_set(66)
out = HERE / "build" / "visual1" / "diag_forcepath.png"
sc.render.filepath = str(out)
bpy.ops.render.render(write_still=True)
print("DIAG_FORCEPATH_DONE", out)
