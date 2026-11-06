"""Export recorded campaign routes as blue OpenSim comparison overlays."""
from pathlib import Path
import json,csv,base64,xml.etree.ElementTree as ET,sys
import numpy as np
import opensim as osim
import vtk
from vtk.util.numpy_support import vtk_to_numpy
ROOT=Path(__file__).resolve().parents[3]
CAMP=ROOT/'Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103'
OUT=CAMP/'audit_20261105';OUT.mkdir(exist_ok=True)
source=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim'
REPAIRED='--repaired' in sys.argv
TAG='repaired' if REPAIRED else 'latest'
records=json.loads((OUT/('reoptimization/repaired_routes.json' if REPAIRED else 'recorded_routes.json')).read_text())
if REPAIRED:assert all(r['nominal_diameter_mm']==20 for r in records if r['group']=='vas_int'),'Combined vasti must use 20 mm BPAs; run fr_apply_vasti_rule_20261105 first.'
states=records if REPAIRED else json.loads((OUT/'mechanical_states.json').read_text())
state_lookup={(r['group'],r['branch']):r for r in states}
osim.Logger.setLevelString('error')
GEOM=Path('C:/OpenSim 4.3/Geometry');osim.ModelVisualizer.addDirToGeometrySearchPaths(str(GEOM))
def vec(v):return np.array([v.get(i) for i in range(3)])
def transform(body,s):
 t=body.getTransformInGround(s);R=t.R();T=np.eye(4)
 T[:3,:3]=[[R.get(i,j) for j in range(3)] for i in range(3)];T[:3,3]=vec(t.p());return T
scene=dict(source=str(source),routes=[],baselineRoutes=[],meshes=[],poses=[],limitations=[])
audit=[]
for branch in ('compliant','rigid'):
 model=osim.Model(str(source));model.setName(TAG+'_routes_'+branch+'_BLUE_COMPARISON')
 for side in ('r','l'):
  knee=model.updCoordinateSet().get('knee_angle_'+side)
  knee.setRangeMin(float(np.radians(-120)));knee.setRangeMax(0);knee.set_clamped(True)
 mus=model.updMuscles()
 for i in range(mus.getSize()):
  gp=mus.get(i).updGeometryPath();gp.setDefaultColor(osim.Vec3(.55,.55,.55))
 for rec in records:
  if rec['branch']!=branch:continue
  bodies=rec['assignedBodies'];P=np.array(rec['recordedLoc']).reshape(-1,3)
  act=osim.PathActuator();act.setName('OPT_'+rec['group']+'_'+branch+'_'+rec['status'])
  act.set_appliesForce(False);act.setOptimalForce(1)
  for i,(b,p) in enumerate(zip(bodies,P)):
   act.addNewPathPoint('opt_'+str(i+1),model.updBodySet().get(b),osim.Vec3(*map(float,p)))
  act.updGeometryPath().setDefaultColor(osim.Vec3(0,.15,1))
  model.addForce(act)
 model.finalizeConnections();s=model.initSystem();model.realizePosition(s)
 target=OUT/('gait2392_robotbody_'+TAG+'_'+branch+'_blue.osim');model.printToXML(str(target))
 # Reopen using OpenSim, verifying both serialization and exact frame use.
 check=osim.Model(str(target));sc=check.initSystem();check.realizePosition(sc)
 for rec in records:
  if rec['branch']!=branch:continue
  name='OPT_'+rec['group']+'_'+branch+'_'+rec['status'];act=osim.PathActuator.safeDownCast(check.getForceSet().get(name))
  pps=act.getGeometryPath().getPathPointSet();world=[]
  for i,(b,p) in enumerate(zip(rec['assignedBodies'],rec['recordedLoc'])):
   point=pps.get(i);T=transform(check.getBodySet().get(b),sc);expected=T[:3,:3]@p+T[:3,3]
   actual=vec(point.getLocationInGround(sc));assert np.linalg.norm(actual-expected)<1e-9
   world.append(actual.tolist())
  color=vec(act.getGeometryPath().getDefaultColor());assert np.allclose(color,[0,.15,1])
  audit.append(dict(group=rec['group'],branch=branch,status=rec['status'],nBPA=rec['nBPA'],gain=rec['gain'],recorded_meet_fraction=rec['fractionMet'],frame_validation='PASS',blue_color_validation='PASS'))
  sm=state_lookup[(rec['group'],branch)]
  scene['routes'].append(dict(**rec,mechanical=sm,neutralWorld=world))
 print('Validated OpenSim model',target)
 if branch=='compliant':
  bodyset=check.getBodySet()
  for a in np.linspace(0,-120,25):
   for side in ('r','l'):check.updCoordinateSet().get('knee_angle_'+side).setValue(sc,float(np.radians(a)),False)
   check.realizePosition(sc)
   scene['poses'].append(dict(knee_deg=float(a),W={bodyset.get(i).getName():transform(bodyset.get(i),sc).tolist() for i in range(bodyset.getSize())}))
   original=check.getMuscles()
   for mi in range(original.getSize()):
    muscle=original.get(mi)
    if not muscle.getName().endswith('_r'):continue
    path=muscle.getGeometryPath().getCurrentPath(sc)
    points=[vec(path.get(j).getLocationInGround(sc)).tolist() for j in range(path.getSize())]
    if a==0:scene['baselineRoutes'].append(dict(name=muscle.getName(),poses=[points]))
    else:next(r for r in scene['baselineRoutes'] if r['name']==muscle.getName())['poses'].append(points)
  # Read geometry from the SAME OpenSim model, including declared scales.
  tree=ET.parse(target)
  for body in tree.findall('.//BodySet/objects/Body'):
   for mesh in body.findall('./attached_geometry/Mesh'):
    fn=mesh.findtext('mesh_file');p=GEOM/fn
    if not p.exists():
     scene['limitations'].append('Missing geometry '+str(p));continue
    reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(p));reader.Update()
    tri=vtk.vtkTriangleFilter();tri.SetInputData(reader.GetOutput());tri.Update();poly=tri.GetOutput()
    v=vtk_to_numpy(poly.GetPoints().GetData()).astype(float)
    scale=np.fromstring(mesh.findtext('scale_factors','1 1 1'),sep=' ');v*=scale
    faces=vtk_to_numpy(poly.GetPolys().GetData()).reshape(-1,4)[:,1:]
    scene['meshes'].append(dict(name=body.get('name')+'/'+mesh.get('name'),body=body.get('name'),vertices=base64.b64encode(v.astype('<f4').tobytes()).decode(),faces=base64.b64encode(faces.astype('<u4').tobytes()).decode()))
scene['limitations']+=['Blue overlays reproduce the saved primary-pass mechanics rows, including its frozen/dropped-point approximations; they are not an approval of those approximations.','Saved candidate statuses and unmet torque targets are retained.','No guessed CAD-to-OpenSim foot or attachment registration is used.']
if REPAIRED:scene['limitations']+=['Full-length BPA vessels are route-tangent packaging candidates, not qualified placements.','Cam equilibrium, independent series-section equilibrium, physical reeving and complete robot-CAD collision remain unresolved.','Knee-only animation changes the knee from 0 to -120 degrees; muscles acting only on stationary joints need not contract.']
if REPAIRED:
 matlab=json.loads((OUT/'reoptimization/matlab_body_frames.json').read_text());errors=[]
 for native,approx in zip(scene['poses'],matlab):
  assert native['knee_deg']==approx['knee_deg']
  for body,T in approx['W'].items():
   A=np.asarray(T);B=np.asarray(native['W'][body]);errors.append(dict(body=body,knee_deg=native['knee_deg'],translation_error_mm=float(1000*np.linalg.norm(A[:3,3]-B[:3,3])),rotation_matrix_error=float(np.max(np.abs(A[:3,:3]-B[:3,:3])))))
 (OUT/'matlab_vs_native_body_frames.json').write_text(json.dumps(errors,indent=2))
 worst=max(e['translation_error_mm'] for e in errors)
 scene['limitations'].append(f'MATLAB natural-spline joint translation differs from native OpenSim by up to {worst:.4f} mm on this knee grid. Blue path local coordinates and displayed bone frames are native OpenSim; world-space vessel candidates use the MATLAB mechanics frames.')
scene_name='repaired_verified_scene.json' if REPAIRED else 'opensim_verified_scene.json'
(OUT/scene_name).write_text(json.dumps(scene,separators=(',',':')))
with (OUT/(TAG+'_overlay_verification.csv' if REPAIRED else 'opensim_overlay_verification.csv')).open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=audit[0]);w.writeheader();w.writerows(audit)
# GUI motion file in degrees. Model overlays remain display-only forces.
with (OUT/'knee_comparison.mot').open('w') as f:
 angles=list(np.linspace(0,-120,25))+list(np.linspace(-115,0,24))
 f.write('knee_comparison\nversion=1\nnRows=49\nnColumns=3\ninDegrees=yes\nendheader\ntime\tknee_angle_r\tknee_angle_l\n')
 for i,a in enumerate(angles):f.write(f'{i/4:.6f}\t{a:.6f}\t{a:.6f}\n')
print('Verified',len(audit),'blue overlays;',len(scene['meshes']),'OpenSim meshes.')
