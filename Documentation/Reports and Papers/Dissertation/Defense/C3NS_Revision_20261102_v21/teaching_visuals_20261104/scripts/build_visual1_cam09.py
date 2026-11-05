"""Cam teaching animation from authoritative 09_BA_003 native CAD geometry.

All tibia-mounted route points are measured from the assembly's named
Theta1 Coordinate System. Native component placements are retained. The
crossed linkage is solved geometrically, with no dynamic or force refit.
"""
import sys,math,json,argparse
from pathlib import Path
import numpy as np
import bpy
sys.path.insert(0,str(Path(__file__).resolve().parent))
import render_common as rc
import cad09_kinematics as K
def build_schedule(phi):
    u=np.arange(120)/119
    s=1-np.abs(2*u-1)
    s=s*s*(3-2*s)
    return float(phi.min())+(float(phi.max())-float(phi.min()))*s

def r3(q):
 Q=np.eye(3);Q[:2,:2]=K.rot(q);return Q

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stills-only',action='store_true')
 argv=sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else []
 args=ap.parse_args(argv)
 scene,_=rc.setup_scene(res=(1440,600),samples=24);rc.add_sun()
 g=np.load(rc.HERE/'data'/'visual1_geom.npz')
 source=json.loads((rc.HERE/'data'/'cam_visuals_dump_20261104.json').read_text())
 raw=np.array(source['route_raw']);active=np.array(source['route_active'],dtype=bool)
 phis=build_schedule(g['phiD']);frames=list(range(len(phis)))
 indices=[int(np.argmin(abs(g['phiD']-p))) for p in phis]
 states=[K.pose(math.radians(float(g['phiD'][i]))) for i in indices]
 home=int(np.argmin(abs(g['phiD'])))
 t1_home=g['pivot'][home]-r3(float(g['tibia_rot'][home]))@g['t1toICR'][home]
 offset=t1_home-K.T1
 def in_femur(xy,z=0):return np.r_[xy,K.T1[2]+z]+offset
 t1=[in_femur(p['D']) for p in states]
 pivots=[in_femur(p['I'],.035) for p in states]

 # Render native parts directly; no STL recentering, Hough fit, or old assembly.
 for row in K.ASSEMBLY['components']:
  name=row['name'];v=np.load(K.DATA/row['native_mesh'])['vertices']
  tr=np.array(row['transform']);world=v@tr[:9].reshape(3,3)+tr[9:12]+offset
  is_link=name.startswith(('04_05','04_06'))
  fixed=name.startswith(('03_01','04_01'))
  if name.startswith('TibiaClip'):color=(.22,.25,.29);alpha=.85
  elif is_link:color=(.12,.14,.17);alpha=1.
  elif fixed:color=(.48,.49,.50);alpha=.30
  else:color=(.45,.47,.49);alpha=.25
  mat=rc.make_material(name+'_mat',color,alpha=alpha)
  if is_link:
   A=K.A if name.startswith('04_06') else K.B
   end=K.C0 if name.startswith('04_06') else K.D0
   pivot=in_femur(A);world-=pivot
   target='C' if name.startswith('04_06') else 'D'
   a0=math.atan2(end[1]-A[1],end[0]-A[0])
   angles=[math.atan2(p[target][1]-A[1],p[target][0]-A[0])-a0 for p in states]
   ob=rc._mesh_from_stl(world,np.arange(len(world)).reshape(-1,3),name,tuple(pivot),mat)
   rc.keyframe_linear(ob,'rotation_euler',frames,angles,index=2)
  elif fixed:
   rc._mesh_from_stl(world,np.arange(len(world)).reshape(-1,3),name,(0,0,0),mat)
  else:
   world-=t1_home
   ob=rc._mesh_from_stl(world,np.arange(len(world)).reshape(-1,3),name,(0,0,0),mat)
   rc.keyframe_linear(ob,'location',frames,[tuple(p) for p in t1])
   rc.keyframe_linear(ob,'rotation_euler',frames,[p['q'] for p in states],index=2)

 routes=[];wrap_paths=[];force_lines=[];arm_lines=[];point_paths=[[] for _ in range(4)]
 for k,(p,i) in enumerate(zip(states,indices)):
  points=raw[:,:,i].copy()
  for j in range(5,9):
   points[j]=r3(p['q'])@raw[j,:,i]+t1[k]
   point_paths[j-5].append(tuple(points[j]))
  pts=[tuple(points[j]) for j in range(9) if active[j,i] or j in (0,8)]
  pts += [pts[-1]]*(9-len(pts));routes.append(pts)
  wrap_paths.append([tuple(points[j]) for j in (5,6,7)])
  direction=points[4]-points[5];direction/=np.linalg.norm(direction)
  origin=pivots[k].copy();origin[2]=0
  foot=points[5]+np.dot(origin-points[5],direction)*direction
  arm_lines.append([tuple(origin+np.array([0,0,.04])),tuple(foot+np.array([0,0,.04]))])
  force_lines.append([tuple(points[5]+np.array([0,0,.045])),tuple(points[5]+.075*direction+np.array([0,0,.045]))])
 cable=rc.poly_curve('optimized_route',routes[0],.0035,material=rc.make_material('route_ink',rc.CABLE_BLUE))
 rc.keyframe_polyline(cable,frames,routes)
 envelope=rc.poly_curve('BPA_design_envelope',wrap_paths[0],float(source['geo']['bpaRadius']),
                        material=rc.make_material('BPA_envelope_ink',rc.CABLE_BLUE,alpha=.24))
 rc.keyframe_polyline(envelope,frames,wrap_paths)
 arm=rc.poly_curve('moment_arm',arm_lines[0],.0027,material=rc.make_material('arm_ink',rc.ARM_GOLD))
 rc.keyframe_polyline(arm,frames,arm_lines)
 force_mat=rc.make_material('force_ink',(.58,.05,.35))
 fs,fh=rc.arrow('force',*force_lines[0],.0025,force_mat)
 def sphere(name,positions,radius,color):
  bpy.ops.mesh.primitive_uv_sphere_add(radius=radius,segments=16,ring_count=8)
  ob=bpy.context.active_object;ob.name=name;rc.set_material(ob,rc.make_material(name+'_mat',color))
  rc.keyframe_linear(ob,'location',frames,[tuple(p) for p in positions]);return ob
 sphere('instantaneous_pivot',pivots,.0045,rc.ARM_GOLD)
 sphere('t1_origin',[p+np.array([0,0,.04]) for p in t1],.0035,(.35,.05,.35))
 for j,points in enumerate(point_paths):sphere('t1_route_point_'+str(j+6),points,.003,(.05,.05,.06))

 labels=[('Pivot: link intersection',(-.13,-.405,.10),(.12,.12,.13)),
         ('t1 reference',(-.15,-.48,.10),(.35,.05,.35)),
         ('Shaped cam',(.19,-.54,.10),(.12,.12,.13)),
         ('Pulling force',(.19,-.345,.10),(.58,.05,.35)),
         ('Moment arm',(.19,-.41,.10),rc.ARM_GOLD)]
 for text,loc,color in labels:
  ob=rc.text_label(text,text,.021,color=color,align='CENTER');ob.location=loc
 lead=rc.poly_curve('pivot_leader',[(-.08,-.405,.07),tuple(pivots[0])],.0006,material=rc.make_material('leader_ink',(.16,.16,.18)))
 rc.keyframe_polyline(lead,frames,[[(-.08,-.405,.07),tuple(p)] for p in pivots])
 tl=rc.poly_curve('t1_leader',[(-.10,-.48,.07),tuple(t1[0]+np.array([0,0,.04]))],.0006,material=rc.make_material('t1_leader',(.35,.05,.35)))
 rc.keyframe_polyline(tl,frames,[[(-.10,-.48,.07),tuple(p+np.array([0,0,.04]))] for p in t1])
 cam=rc.look_at((.04,-.43,1.2),(.04,-.43,0),lens=50);cam.data.type='ORTHO';cam.data.ortho_scale=.78
 readout=rc.text_label('','angle_readout',.025,align='CENTER');readout.location=(.04,-.285,.11)
 def update(k):
  scene.frame_set(k);readout.data.body=f'Knee angle: {g["phiD"][indices[k]]:+.0f} deg'
  rc.aim_arrow(fs,fh,*force_lines[k])
 rc.save_blend(rc.OUTPUT/'scenes'/'visual1_cam09.blend')
 if args.stills_only:
  k=next(k for k in range(18,60) if phis[k]>=-30)
  scene.render.resolution_x=1440;scene.render.resolution_y=600
  update(k);rc.render_frame(scene,k,str(rc.OUTPUT/'visual1'/'still_labeled.png'))
  for ob in bpy.data.objects:
   if ob.type=='FONT' or 'leader' in ob.name:ob.hide_render=True
  rc.render_frame(scene,k,str(rc.OUTPUT/'visual1'/'still_unlabeled.png'))
  print('CAM09_STILLS_DONE');return
 folder=rc.OUTPUT/'visual1'/'frames';folder.mkdir(parents=True,exist_ok=True)
 for k in frames:
  update(k);rc.render_frame(scene,k,str(folder/f'f{k:04d}.png'))
 print('CAM09_DONE')

if __name__=='__main__':main()
