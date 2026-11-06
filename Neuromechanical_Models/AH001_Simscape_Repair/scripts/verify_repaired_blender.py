"""Reopen the saved blend and check evaluated meshes, frames and animation."""
from pathlib import Path
import json,math
import bpy,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/audit_20261105'
D=json.loads((OUT/'repaired_verified_scene.json').read_text());report=[]
assert len(bpy.data.scenes)==2
for branch in ('compliant','rigid'):
 s=bpy.data.scenes['REPAIRED '+branch.upper()+' routes - candidate audit'];bpy.context.window.scene=s
 assert s.frame_start==1 and s.frame_end==289 and s.unit_settings.system=='METRIC'
 assert np.allclose(s.world.color,(1,1,1)) and s.camera is not None
 routes=[r for r in D['routes'] if r['branch']==branch]
 objects=[o for o in s.objects if 'routing_group' in o]
 expected=sum(len(r['sleeves'][0]) if isinstance(r['sleeves'][0],list) else 1 for r in routes)
 assert len(objects)==expected
 checks=0;length_error=0.;radius_error=0.
 for pose_index in (0,12,24):
  s.frame_set(1+6*pose_index);dg=bpy.context.evaluated_depsgraph_get()
  for r in routes:
   caps=r['sleeves'][pose_index];caps=caps if isinstance(caps,list) else [caps]
   for cap in caps:
    ob=next(o for o in objects if o['routing_group']==r['group'] and o['channel_index']==cap['branch'] and o['section_index']==cap['section'])
    ev=ob.evaluated_get(dg);v=np.array([tuple(p.co) for p in ev.data.vertices]);a=v[:12].mean(axis=0);b=v[12:].mean(axis=0)
    le=abs(np.linalg.norm(b-a)-np.linalg.norm(np.array(cap['p2'])-cap['p1']))
    re=max(abs(np.linalg.norm(v[:12]-a,axis=1)-cap['r']))
    assert le<2e-6 and re<2e-6,(ob.name,pose_index,le,re)
    length_error=max(length_error,le);radius_error=max(radius_error,re);checks+=1
  for body,T in D['poses'][pose_index]['W'].items():
   ob=bpy.data.objects[branch+' '+body];assert np.max(np.abs(np.array(ob.matrix_local)-T))<2e-6
 orbit=bpy.data.objects['180 degree camera orbit '+branch]
 s.frame_set(1);assert abs(orbit.rotation_euler.z)<1e-6
 s.frame_set(289);assert abs(orbit.rotation_euler.z-math.pi)<1e-6
 report.append(dict(branch=branch,pressure_vessels=expected,evaluated_mesh_checks=checks,max_length_error_m=length_error,max_radius_error_m=radius_error,body_frames='PASS',white_world='PASS',camera_orbit='PASS',knee_cycle='0 to -120 and back'))
(OUT/'blender_reopen_verification.json').write_text(json.dumps(report,indent=2))
print('PASS saved Blender geometry/animation:',report)
