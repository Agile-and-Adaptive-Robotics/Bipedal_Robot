"""Measure evaluated mesh rings, independently of the diameter-law function."""
import bpy,json,numpy as np
from pathlib import Path
s=bpy.context.scene;before=s.frame_current
o=next(o for o in s.objects if o.name.startswith('vas_int_r / BPA 1 / sleeve'))
rows=[]
for frame in (1,109):
 s.frame_set(frame);bpy.context.view_layer.update()
 v=o.evaluated_get(bpy.context.evaluated_depsgraph_get()).data.vertices
 ring=np.array([v[i].co[:] for i in range(12)]);center=ring.mean(axis=0)
 rows.append(dict(frame=frame,measured_diameter_mm=float(2*np.linalg.norm(ring-center,axis=1).mean()*1000)))
s.frame_set(before)
out=Path(__file__).resolve().parents[1]/'packaging/output/measured_blender_diameters.json'
out.write_text(json.dumps(rows,indent=2));print(rows)
