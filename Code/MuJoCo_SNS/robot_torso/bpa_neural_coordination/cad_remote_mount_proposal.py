"""Source-CAD external sleeve placement, with unverified Bowden cable explicitly separate."""
import sys,os,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent
sys.path.insert(0,str(ROBOT));os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco
from scipy.spatial import ConvexHull
from cad_mount_subset import inflated_hull_intersection
OUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mount_subset')

def boxplanes(size):return np.c_[np.r_[np.eye(3),-np.eye(3)],-np.tile(size,2)]

def main():
 modelpath=ROBOT/'ah001_neural.xml';routepath=ROBOT/'robotbody_routes.json';motion=ROBOT/'runs/verified65_exact_fast.npz';m=mujoco.MjModel.from_xml_path(str(modelpath));d=mujoco.MjData(m);config=json.loads(routepath.read_text());routes={r['name']:r for r in config['routes']};q=np.load(motion)['qpos'];indices=np.linspace(0,len(q)-1,17,dtype=int)
 envelope_path=OUT.parent/'route_envelopes/CAD_ROUTE_ENVELOPE_REPORT.json';envelope=json.loads(envelope_path.read_text());bounds={r['name']:r for r in envelope['routes']}
 for p in (modelpath,routepath,motion):
  if hashlib.sha256(p.read_bytes()).hexdigest()!=envelope['source_sha256'][p.name]:raise ValueError('Geometry envelope mismatch')
 hulls={}
 for g in range(m.ngeom):
  kind=m.geom_type[g]
  if kind==mujoco.mjtGeom.mjGEOM_PLANE:continue
  if kind==mujoco.mjtGeom.mjGEOM_MESH:
   mid=m.geom_dataid[g];v=m.mesh_vert[m.mesh_vertadr[mid]:m.mesh_vertadr[mid]+m.mesh_vertnum[mid]];hulls[g]=ConvexHull(v).equations
  elif kind==mujoco.mjtGeom.mjGEOM_SPHERE:hulls[g]=boxplanes(np.repeat(m.geom_size[g,0],3))
  elif kind==mujoco.mjtGeom.mjGEOM_BOX:hulls[g]=boxplanes(m.geom_size[g])
  elif kind in (mujoco.mjtGeom.mjGEOM_CAPSULE,mujoco.mjtGeom.mjGEOM_CYLINDER):hulls[g]=boxplanes(np.array([m.geom_size[g,0]]*2+[m.geom_size[g,0]+m.geom_size[g,1]]))
  else:raise ValueError('Unhandled obstacle geometry '+str(kind))
 d.qpos[:]=q[0];mujoco.mj_forward(m,d);rows=[];radius=.01925;fit=.025;clearance=.002
 for side in ('r','l'):
  body=m.body('thigh_'+side).id;rot=d.xmat[body].reshape(3,3);local_vertices=[]
  for g in np.flatnonzero(m.geom_bodyid==body):
   if m.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
   mid=m.geom_dataid[g];v=m.mesh_vert[m.mesh_vertadr[mid]:m.mesh_vertadr[mid]+m.mesh_vertnum[mid]];world=v@d.geom_xmat[g].reshape(3,3).T+d.geom_xpos[g];local_vertices.append((world-d.xpos[body])@rot)
  vertices=np.vstack(local_vertices);lo=vertices.min(0);hi=vertices.max(0)
  for base,direction in (('bifemsh',-1),('vas_med',1)):
   name=base+'_'+side;route=routes[name]
   if route['points'][0]['body']!='thigh_'+side or route['points'][0]['motion_functions']:raise ValueError('Need fixed source outlet on same thigh body')
   excursion=bounds[name]['route_excursion_m'];paired_excursion=max(bounds[base+'_'+s]['route_excursion_m'] for s in ('r','l'));rest=paired_excursion/.225;span=rest+2*fit;possible=span<=hi[1]-lo[1]-2*clearance
   row=dict(name=name,body='thigh_'+side,source_route_eye_body_local_m=route['points'][0]['position_m'],thigh_mesh_bounds_local_m=[lo.tolist(),hi.tolist()],route_excursion_m=excursion,derived_sleeve_rest_m=rest,nominal_diameter_mm=20,source_inflated_radius_m=radius,fittings_each_m=fit,geometric_max_contraction_fraction=.225,rest_formula='401-pose route excursion /0.225; retains0.025 strain-fraction margin below declared0.25 maximum before force-dependent series extension',straight_mount_span_m=span,available_thigh_axis_m=float(hi[1]-lo[1]),longitudinal_packaging_pass=bool(possible),alternatives=[])
   if possible:
    for extra in np.arange(0,.080001,.01):
     centre=(lo+hi)/2;centre[0]=(lo[0]-radius-clearance-extra) if direction<0 else (hi[0]+radius+clearance+extra);eye=centre.copy();eye[1]+=span/2;guide=eye.copy();guide[1]-=span;a=eye.copy();a[1]-=fit;b=a.copy();b[1]-=rest;hits=set()
     for index in indices:
      d.qpos[:]=q[index];mujoco.mj_forward(m,d);R=d.xmat[body].reshape(3,3);wa=R@a+d.xpos[body];wb=R@b+d.xpos[body]
      for g,planes in hulls.items():
       gr=d.geom_xmat[g].reshape(3,3);ga=gr.T@(wa-d.geom_xpos[g]);gb=gr.T@(wb-d.geom_xpos[g])
       if inflated_hull_intersection(ga,gb,planes,radius):hits.add(m.geom(g).name)
     attempt=dict(additional_surface_offset_m=float(extra),proximal_fixed_eye_local_m=eye.tolist(),distal_fixed_housing_guide_local_m=guide.tolist(),sleeve_axis_local_m=[a.tolist(),b.tolist()],outer_hull_intersections=sorted(hits),mesh_separation_certified=not hits)
     row['alternatives'].append(attempt)
     if not hits:
      row['selected_mount']=attempt;row['minimum_fixed_housing_endpoint_distance_m']=float(np.linalg.norm(guide-np.array(route['points'][0]['position_m'])));break
   row['sleeve_placement_certified']=bool(row.get('selected_mount'));rows.append(row)
 result=dict(machine='easteregg2',model=str(modelpath),source_motion=str(motion),sampled_pose_indices=indices.tolist(),routes=rows,certified_sleeve_placements=[r['name'] for r in rows if r['sleeve_placement_certified']],fully_calibrated_or_mounted_actuator_assignments=0,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (modelpath,routepath,motion,envelope_path,Path(__file__))},scope=['Engineering sleeve-mount proposal derived from actual thigh mesh bounds and source-route excursions, not an existing mounted mechanism.','Remote sleeve eye and housing guide are fixed on the same thigh; Bowden output would be anchored at the original fixed same-body path eye. Constant housing/remote-span length does not change the original route Jacobian in the ideal lossless model.','Only sleeve envelopes have a conservative separation certificate at17sampled whole-body poses. Cables, bend radius, fittings, eye/bracket design, attachment loads, manufacturing, BPAs against each other, and motion between poses still need verification.','No exact Bowden housing length, friction, series compliance, pressure response, sleeve/endcap mass, or calibrated force capacity is invented. These remain unassigned; no pressure simulation or full robot standing is claimed.','Inflated38.5mm diameter is the historical Routing27 envelope assumption; additional2mm clearance is an explicit engineering assumption. Length budgets use the401-pose source envelope.'])
 result['bilateral_sleeve_dimensions_tied']=True
 for row in rows:row['rest_formula']='Maximum401-pose excursion of either leg /0.225; paired sleeves share length/diameter, with0.025 geometric strain-fraction margin before series extension.'
 (OUT/'PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(certified=result['certified_sleeve_placements'],longitudinal_pass=[r['name'] for r in rows if r['longitudinal_packaging_pass']]),indent=2))

if __name__=='__main__':main()
