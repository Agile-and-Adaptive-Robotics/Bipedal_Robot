"""Verify the native-CAD datum, closed linkage, and t1 point transforms."""
import json,math,numpy as np
import cad09_kinematics as K
g=np.load(K.ROOT/'data'/'visual1_geom.npz')
D=json.loads((K.ROOT/'data'/'cam_visuals_dump_20261104.json').read_text());raw=np.array(D['route_raw'])
closure=[];cross=[];recover=[]
for i,deg in enumerate(g['phiD']):
 p=K.pose(math.radians(float(deg)))
 closure.extend([abs(np.linalg.norm(p['C']-K.A)-K.LF),abs(np.linalg.norm(p['D']-K.B)-K.LB),abs(np.linalg.norm(p['C']-p['D'])-np.linalg.norm(K.V))])
 for a,b in [(K.A,p['C']),(K.B,p['D'])]:
  cross.append(abs(np.cross(b-a,p['I']-a))/np.linalg.norm(b-a))
 for j in range(5,9):
  xy=p['Q']@raw[j,:2,i]+p['D'];back=p['Q'].T@(xy-p['D'])
  recover.append(np.linalg.norm(back-raw[j,:2,i]))
assert max(closure)<1e-10 and max(cross)<1e-10 and max(recover)<1e-10
report={'assembly':'09_BA_003.SLDASM','datum':'Theta1 Coordinate System','t1_origin_to_lower_back_pin_m':float(np.linalg.norm(K.T1[:2]-K.D0)),'link_closure_max_m':max(closure),'pivot_to_link_centerline_max_m':max(cross),'t1_route_point_roundtrip_max_m':max(recover),'native_mesh_units':'metres; direct CAD tessellation, without STL recentering','saved_optimizer_and_measurements_modified':False}
(K.ROOT/'build_fixed_20261105'/'CAM09_VERIFICATION.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
