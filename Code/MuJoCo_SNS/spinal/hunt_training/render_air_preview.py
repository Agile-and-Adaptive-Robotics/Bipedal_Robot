"""Render recorded full states; no resimulation or model-state correction."""
import os,sys,argparse,json
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H.parent)]
import numpy as np,mujoco,bsolve_ik as B
from PIL import Image,ImageDraw,ImageFont
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('recording',type=Path)
parser.add_argument('--start',type=float,default=5.);parser.add_argument('--end',type=float,default=13.)
parser.add_argument('--fps',type=float,default=12.);parser.add_argument('--playback-speed',type=float,default=1.)
args=parser.parse_args()
if args.fps<=0 or args.playback_speed<=0 or args.end<=args.start:parser.error('Use positive fps/speed and end greater than start')
source=args.recording.resolve()
metadata_path=source.with_suffix('.result.json')
if metadata_path.exists():
    metadata=json.loads(metadata_path.read_text())
    if metadata.get('stage') not in (1,2):parser.error('This renderer labels air stepping; use an air-stage recording')
    if not metadata.get('assessment',{}).get('complete',False):parser.error('Incomplete recording; no curriculum preview generated')
z=np.load(source,allow_pickle=True);model,data=B.build_air_model(False)
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
sample_step=float(np.median(np.diff(z['t'])))
stride=max(1,round(1./(args.fps*sample_step)))
indices=np.flatnonzero((z['t']>=args.start)&(z['t']<args.end))[::stride];frames=[]
if len(indices)<2:parser.error('At least two frames are required in the selected time range')
frame_times=z['t'][indices]
deltas=np.r_[np.diff(frame_times),np.median(np.diff(frame_times))]/args.playback_speed
if np.min(deltas)<.01-1e-9:parser.error('GIF frames need at least 10 ms; lower fps or playback speed')
# Quantize the cumulative clock, avoiding drift from rounding every frame.
ticks=np.rint(np.r_[0.,np.cumsum(deltas)]*100).astype(int)
durations=(np.diff(ticks)*10).tolist()
for k in indices:
    data.qpos[:]=z['qfull'][k];mujoco.mj_forward(model,data);views=[]
    # Native world +x is anterior; MuJoCo azimuth 180 views from +x.
    for azimuth in (180.,140.):
        cam=mujoco.MjvCamera();cam.lookat[:]=[data.xpos[model.body('pelvis').id,0],data.xpos[model.body('pelvis').id,1],.85]
        cam.distance=2.8;cam.azimuth=azimuth;cam.elevation=-6.
        renderer.update_scene(data,camera=cam,scene_option=opt);views.append(Image.fromarray(renderer.render().copy()))
    frame=Image.new('RGB',(880,650),'white');frame.paste(views[0],(0,50));frame.paste(views[1],(440,50));draw=ImageDraw.Draw(frame)
    draw.text((12,7),f'Gait2392: neural air stepping | full vertical support | {args.playback_speed:g}x playback',fill='black',font=font)
    draw.text((12,29),'Anterior',fill='black',font=font);draw.text((452,29),f'Oblique | t = {z["t"][k]:.2f} s',fill='black',font=font)
    frames.append(frame)
del renderer  # MuJoCo 2.3.7 releases its renderer through destruction.
frames[0].save(out/'stage2_air_neural.gif',save_all=True,append_images=frames[1:],duration=durations,loop=0)
frames[len(frames)//2].save(out/'stage2_air_neural.png')
(out/'provenance.json').write_text(json.dumps(dict(source=str(source),result_metadata=str(metadata_path),scope='Recorded full qpos; no corrected poses, resimulation, or walking claim. See TRAINING_RECORD.md for robustness evidence.',camera_native_world_azimuth_deg=[180.,140.],frames=len(frames),frame_times_s=frame_times.tolist(),playback_speed=args.playback_speed,gif_duration_ms=sum(durations)),indent=2))
print(out,flush=True)
