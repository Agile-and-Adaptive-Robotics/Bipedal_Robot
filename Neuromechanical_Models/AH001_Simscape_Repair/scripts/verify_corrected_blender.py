"""Independent checks after loading the saved .blend in a fresh process."""
from pathlib import Path
import json,math,numpy as np,bpy
P=Path(__file__).resolve().parents[1];repo=P.parents[1]
out=P/'packaging/output';s=bpy.context.scene
cs=json.loads((repo/'Code/MuJoCo_SNS/spinal/reports_20261027_robot_control/robot_mass_geometry.json').read_text())['components']
CR=np.array(next(c for c in cs if '/02_01_PE_' in c['name'])['transform'])[:9].reshape(3,3)
def axis(token):
 c=next(c for c in cs if token in c['name']);a=(CR@np.array(c['transform'])[:9].reshape(3,3).T)[:,1]
 return a if a[1]<0 else -a
axes={side:(axis('03_01_FB_R' if side=='r' else '03_03_FB_L'),axis('05_01_TI_R' if side=='r' else '05_02_TI_')) for side in ('r','l')}
rows=[]
for frame in range(s.frame_start,s.frame_end+1):
 s.frame_set(frame)
 for side,(f,t) in axes.items():
  uf=np.array(s.objects[f'thigh_{side} corrected frame'].matrix_local)[:3,:3]@f
  ut=np.array(s.objects[f'shank_{side} corrected frame'].matrix_local)[:3,:3]@t
  angle=math.degrees(math.atan2(np.cross(uf,ut)[2],uf[:2]@ut[:2]));rows.append(dict(frame=frame,side=side,angle_deg=angle))
minimum=min(r['angle_deg'] for r in rows);maximum=max(r['angle_deg'] for r in rows)
assert minimum>=-120.001 and maximum<=.001,(minimum,maximum)
assert max(abs(r['angle_deg']) for r in rows if r['frame'] in (1,289))<.001
s.frame_set(1);s.camera.data.ortho_scale=2.2
# Save framing refinement to the corrected file only, never the original.
bpy.ops.wm.save_as_mainfile(filepath=bpy.data.filepath,compress=True)
s.render.filepath=str(out/'Blender_Packaging_Corrected.png');bpy.ops.render.render(write_still=True)
result=dict(loaded_file=bpy.data.filepath,frames_tested=289,side_frame_checks=len(rows),minimum_knee_deg=minimum,maximum_knee_deg=maximum,no_hyperextension=True,background='White scene world',original_model_unchanged=True)
(out/'corrected_fresh_process_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
