"""Bounded source-derived sleeve packaging screen; no hardware identification."""
import os,sys,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);sys.path.insert(0,str(ROBOT))
import numpy as np,mujoco
from scipy.spatial import ConvexHull
from dynamic_routes import DynamicRoutes
OUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mount_subset')

def inflated_hull_intersection(a,b,planes,r):
 """Conservative capsule rejection: intersection with an outer offset polytope.

 Empty interval certifies separation; nonempty does NOT prove actual collision.
 """
 value=planes[:,:3]@a+planes[:,3]-r;slope=planes[:,:3]@(b-a);lo=0.;hi=1.
 for v,s in zip(value,slope):
  if abs(s)<1e-14:
   if v>0:return False
  elif s>0:hi=min(hi,-v/s)
  else:lo=max(lo,-v/s)
  if lo>hi:return False
 return True

def main():
 OUT.mkdir(exist_ok=True);modelpath=ROBOT/'ah001_neural.xml';routepath=ROBOT/'robotbody_routes.json';source=ROBOT/'runs/verified65_exact_fast.npz';history=HERE/'runs/cad_parameter_audit_20261031/HISTORICAL_ROUTING_CANDIDATES.json'
 m=mujoco.MjModel.from_xml_path(str(modelpath));d=mujoco.MjData(m);config=json.loads(routepath.read_text());dyn=DynamicRoutes(m,config);z=np.load(source);q=z['qpos'];sample_indices=np.linspace(0,len(q)-1,17,dtype=int)
 hist={r['actuator_id']:r for r in json.loads(history.read_text())['candidates'] if r['configuration']=='single'}
 envelope_path=OUT.parent/'route_envelopes/CAD_ROUTE_ENVELOPE_REPORT.json';envelope=json.loads(envelope_path.read_text());bounds={r['name']:r for r in envelope['routes']}
 for p in (modelpath,routepath,source):
  if hashlib.sha256(p.read_bytes()).hexdigest()!=envelope['source_sha256'][p.name]:raise ValueError('Cached envelope/source mismatch: '+str(p))
 candidates=[]
 for route,tid,ids in zip(dyn.routes,dyn.tendon_ids,dyn.site_ids):
  assign=route.get('actuator_assignment')
  if not assign or assign['actuator_id_per_side'] not in hist:continue
  h=hist[assign['actuator_id_per_side']];candidates.append(dict(name=route['name'],tid=tid,ids=ids,rest=h['rest_m'],fit=h['fitting_parameter_m'],history=h,lengths=[],first_span=[],colliders=set(),origins=[],directions=[]))
 hulls={}
 for g in range(m.ngeom):
  if m.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
  mid=m.geom_dataid[g];start=m.mesh_vertadr[mid];vertices=m.mesh_vert[start:start+m.mesh_vertnum[mid]]
  hulls[g]=ConvexHull(vertices).equations
 # Source historical context explicitly uses a fully inflated20mm BPA radius38.5/2mm.
 radius=.0385/2
 for index in sample_indices:
  d.qpos[:]=q[index];dyn.update(d)
  for c in candidates:
   route_i=dyn.tendon_ids.index(c['tid']);ids=dyn.site_ids[route_i][dyn.active[route_i]];points=d.site_xpos[ids];delta=points[1]-points[0];span=np.linalg.norm(delta);unit=delta/span
   c['lengths'].append(float(d.ten_length[c['tid']]));c['first_span'].append(float(span));c['origins'].append(points[0].tolist());c['directions'].append(unit.tolist())
   a=points[0]+c['fit']*unit;b=a+c['rest']*unit
   for g,planes in hulls.items():
    rotation=d.geom_xmat[g].reshape(3,3);local_a=rotation.T@(a-d.geom_xpos[g]);local_b=rotation.T@(b-d.geom_xpos[g])
    if inflated_hull_intersection(local_a,local_b,planes,radius):c['colliders'].add(m.geom(g).name)
 rows=[]
 for c in candidates:
  lo=bounds[c['name']]['route_min_m'];hi=bounds[c['name']]['route_max_m'];fixed=hi-c['rest'];cable=fixed-2*c['fit'];strain=(hi-lo)/c['rest'];span_ok=min(c['first_span'])>=c['rest']+2*c['fit'];length_ok=cable>=0 and strain<=.25;clear=not c['colliders']
  rows.append(dict(name=c['name'],source_design=c['history'],resting_sleeve_m=c['rest'],nominal_diameter_mm=20,fully_inflated_radius_m=radius,two_fittings_each_m=c['fit'],proposed_fixed_nonmuscle_m=fixed,proposed_cable_m=cable,max_zero_stretch_contraction_fraction=strain,route_length_range_m=[lo,hi],minimum_proximal_straight_span_m=min(c['first_span']),length_budget_pass=length_ok,proximal_span_pass=span_ok,conservative_mesh_separation_certified=clear,outer_hull_intersections=sorted(c['colliders']),candidate_pass=bool(length_ok and span_ok and clear),proximal_axis_origins_world_m=c['origins'],proximal_axis_unit_directions_world=c['directions']))
 report=dict(machine='easteregg2',source_motion=str(source),sampled_pose_indices=sample_indices.tolist(),candidate_count=len(rows),certified_candidate_count=sum(r['candidate_pass'] for r in rows),routes=rows,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (modelpath,routepath,source,history,Path(__file__))},assumptions=['Historical single20mm sleeve resting length and two25mm fittings are candidate design inputs, not mounted-CAD measurements.','Cable length is derived from current sampled route maximum minus historical sleeve rest and fittings; old cable length is not silently carried across.','Sleeve assumed to pivot from the actual proximal path eye and occupy its first straight span. Guide/bracket/eye hardware, attachment loads, and swept motion between sampled poses remain unverified.','Radius19.25mm is the historical Routing27 fully-inflated envelope assumption, not new hardware identification.','Actual CAD mesh convex hulls are conservative obstacles; an intersection is a conservative rejection, not a proven mesh collision. Only empty segment/offset-hull intervals certify separation at the17sampled poses.','No force/compliance/dynamic-pressure parameter assignment or whole-robot neural-standing claim. No cable/other-BPA/guide packaging validation.'],moving_route_interface='Use complete DynamicRoutes.length_jacobian and -J.T@positive_tension with all native reference gains disabled, OR use native negative tendon force plus extra moving-site correction exactly once. Finite-difference/power parity must precede moving-CAD pressure trials.')
 (OUT/'MOUNT_SUBSET_SCREEN.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(candidates=len(rows),certified=[r['name'] for r in rows if r['candidate_pass']],length_only=[r['name'] for r in rows if r['length_budget_pass']],span_pass=[r['name'] for r in rows if r['proximal_span_pass']]),indent=2))

if __name__=='__main__':main()
