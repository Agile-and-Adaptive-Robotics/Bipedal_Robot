"""Render recorded poses without correcting them, and audit actual foot clearance."""
import os,sys,json,subprocess
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
H=Path(__file__).resolve().parent;sys.path.insert(0,str(H.parent))
import numpy as np,mujoco,bsolve_ik as B
from PIL import Image,ImageDraw,ImageFont
p=Path(sys.argv[1]).resolve();r=json.loads(p.with_suffix('.result.json').read_text());assert r['assessment']['complete'] and not r.get('teacher_diagnostic',False)
subprocess.run([sys.executable,str(H/'audit_foot_clearance.py'),str(p)],check=True)
q=np.load(p.with_suffix('.clearance.npz'));rows=[]
for side in (0,1):
 air=(q['height_m'][:,side]>.001)&(q['load_bw'][:,side]<.025);bounds=np.r_[0,np.flatnonzero(air[1:]!=air[:-1])+1,len(air)];swings=[]
 for a,b in zip(bounds[:-1],bounds[1:]):
  if air[a] and b-a>=4 and q['height_m'][a:b,side].max()>=.01:swings.append(dict(start_s=float(q['t'][a]),end_s=float(q['t'][b-1]+.02),duration_s=float((b-a)*.02),peak_m=float(q['height_m'][a:b,side].max())))
 rows.append(swings)
report={'feet':rows,'criterion':'Force below .025BW, whole visual foot mesh >1mm for at least80ms, peak >=10mm; engineering diagnostic.'};p.with_suffix('.sustained_clearance.json').write_text(json.dumps(report,indent=2));print('SUSTAINED',json.dumps(report),flush=True)
z=np.load(p,allow_pickle=True);m,d=B.build_air_model(False);m.vis.global_.offwidth=900;m.vis.global_.offheight=800
plane=m.geom('ground-plane').id;m.geom_matid[plane]=-1;m.geom_rgba[plane]=[.89,.9,.91,1.]
for b in range(m.nbody):
 name=m.body(b).name or ''
 if name.endswith(('_r','_l')):
  for g in range(m.body_geomadr[b],m.body_geomadr[b]+m.body_geomnum[b]):m.geom_rgba[g]=[.27,.47,.67,1] if name.endswith('_r') else [.93,.47,.20,1]
renderer=mujoco.Renderer(m,height=600,width=440);opt=mujoco.MjvOption();opt.flags[mujoco.mjtVisFlag.mjVIS_TENDON]=False
font=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',17);frames=[];ids=np.flatnonzero((z['t']>=5)&(z['t']<13))[::50]
support=r['ground_contacts']['mean_harness_vertical_bw'] if 'ground_contacts' in r else .768
for i in ids:
 d.qpos[:]=z['qfull'][i];mujoco.mj_forward(m,d);views=[]
 for azimuth in (180.,140.):
  camera=mujoco.MjvCamera();camera.lookat[:]=[d.xpos[m.body('pelvis').id,0],d.xpos[m.body('pelvis').id,1],.85];camera.distance=2.8;camera.azimuth=azimuth;camera.elevation=-6
  renderer.update_scene(d,camera=camera,scene_option=opt)
  # Draw a visual-only 100mm floor grid; recorded poses remain untouched.
  for offset in np.linspace(-1.,1.,21):
   for axis in (0,1):
    g=renderer.scene.geoms[renderer.scene.ngeom];pos=np.array([d.xpos[m.body('pelvis').id,0],d.xpos[m.body('pelvis').id,1],d.geom_xpos[plane,2]+.0004]);pos[1-axis]+=offset
    size=np.array([.0006,.0006,.0003]);size[axis]=1.
    mujoco.mjv_initGeom(g,mujoco.mjtGeom.mjGEOM_BOX,size,pos,np.eye(3).ravel(),np.array([.5,.52,.54,1.],dtype=np.float32));renderer.scene.ngeom+=1
  views.append(Image.fromarray(renderer.render().copy()))
 f=Image.new('RGB',(880,680),'white');f.paste(views[0],(0,80));f.paste(views[1],(440,80));draw=ImageDraw.Draw(f)
 draw.text((10,5),'Gait2392 | neural supported ground stepping | 1x playback',font=font,fill='black')
 draw.text((10,29),f'Measured assistance {support:.0%} BW | 0 deg knee ceiling; 2 deg protective margin | grid 100 mm',font=font,fill='black')
 draw.text((10,53),f'Anterior and oblique views | preliminary stage-4 result | t={z["t"][i]:.1f}s',font=font,fill='black');frames.append(f)
del renderer
out=p.parent/'preview'/p.stem;out.mkdir(parents=True,exist_ok=True);frames[0].save(out/'supported_ground.gif',save_all=True,append_images=frames[1:],duration=100,loop=0);frames[len(frames)//2].save(out/'supported_ground.png')
(out/'provenance.json').write_text(json.dumps(dict(recording=str(p),frame_times=z['t'][ids].tolist(),scope='Recorded qpos, no corrected poses or teacher; actual measured support; earlier-stage regression pending.'),indent=2));print(out,flush=True)
