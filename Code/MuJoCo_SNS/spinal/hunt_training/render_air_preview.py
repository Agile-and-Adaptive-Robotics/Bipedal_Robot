"""Render recorded full states; no resimulation or model-state correction."""
import os,sys,argparse,json
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H.parent)]
import numpy as np,mujoco,bsolve_ik as B
from PIL import Image,ImageDraw,ImageFont
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('recording',type=Path);args=parser.parse_args()
source=args.recording.resolve();z=np.load(source,allow_pickle=True);model,data=B.build_air_model(False)
assert z['qfull'].shape[1]==model.nq
out=source.parent/'preview';out.mkdir(exist_ok=True)
model.vis.global_.offwidth=800;model.vis.global_.offheight=800
renderer=mujoco.Renderer(model,height=600,width=440);opt=mujoco.MjvOption();opt.flags[mujoco.mjtVisFlag.mjVIS_TENDON]=False
for b in range(model.nbody):
    name=model.body(b).name or ''
    if name.endswith(('_r','_l')):
        color=[.27,.47,.67,1.] if name.endswith('_r') else [.93,.47,.20,1.]
        for geom in range(model.body_geomadr[b],model.body_geomadr[b]+model.body_geomnum[b]):model.geom_rgba[geom]=color
font=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',18)
indices=np.flatnonzero((z['t']>=5.)&(z['t']<13.))[::25];frames=[]
for k in indices:
    data.qpos[:]=z['qfull'][k];mujoco.mj_forward(model,data);views=[]
    for azimuth in (0.,40.):
        cam=mujoco.MjvCamera();cam.lookat[:]=[data.xpos[model.body('pelvis').id,0],data.xpos[model.body('pelvis').id,1],.85]
        cam.distance=2.8;cam.azimuth=azimuth;cam.elevation=-6.
        renderer.update_scene(data,camera=cam,scene_option=opt);views.append(Image.fromarray(renderer.render().copy()))
    frame=Image.new('RGB',(880,650),'white');frame.paste(views[0],(0,50));frame.paste(views[1],(440,50));draw=ImageDraw.Draw(frame)
    draw.text((12,7),'Gait2392: neural air stepping | full vertical support',fill='black',font=font)
    draw.text((12,29),'Anterior',fill='black',font=font);draw.text((452,29),f'Oblique | t = {z["t"][k]:.2f} s',fill='black',font=font)
    frames.append(frame)
renderer.close();frames[0].save(out/'stage2_air_neural.gif',save_all=True,append_images=frames[1:],duration=50,loop=0)
frames[len(frames)//2].save(out/'stage2_air_neural.png')
(out/'provenance.json').write_text(json.dumps(dict(source=str(source),scope='Recorded full qpos; no corrected poses, resimulation, or walking claim. Nominal run; basin checks pending.',frames=len(frames)),indent=2))
print(out,flush=True)
