"""Offscreen replay only: y-up camera, visible support, no image synthesis."""
import argparse,json,hashlib
from pathlib import Path
import xml.etree.ElementTree as E
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from import_virtual_walker import mujoco,fmt
def render(folder,start=2.,stop=4.7,title=None,guidance_caption=None,detail_caption=None,still_time=None):
 folder=Path(folder);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());tree=E.parse(folder/'model.xml');root=tree.getroot();world=root.find('worldbody');body=world.find("body[@name='Root']")
 eye=np.array([2.,.1,2.]);target=np.array([0.,-.45,0.]);forward=target-eye;forward/=np.linalg.norm(forward);right=np.cross(forward,np.array([0.,1.,0.]));right/=np.linalg.norm(right);up=np.cross(right,forward)
 E.SubElement(body,'camera',name='review_oblique',pos=fmt(eye),xyaxes=fmt(np.r_[right,up]),mode='trackcom',fovy='35')
 stat=root.find('statistic')
 if stat is None:stat=E.SubElement(root,'statistic')
 stat.set('extent','3')
 asset=root.find('asset')
 if asset is None:asset=E.SubElement(root,'asset')
 E.SubElement(asset,'texture',type='skybox',builtin='gradient',rgb1='.95 .96 .98',rgb2='.88 .9 .94',width='256',height='256')
 E.SubElement(asset,'texture',name='belt_checker',type='2d',builtin='checker',rgb1='.38 .4 .43',rgb2='.7 .72 .75',width='128',height='128')
 E.SubElement(asset,'material',name='belt_material',texture='belt_checker',texrepeat='12 3',reflectance='.05')
 visual=root.find('visual');globalv=visual.find('global');globalv.set('offwidth','900');globalv.set('offheight','700')
 E.SubElement(visual,'headlight',ambient='.4 .4 .4',diffuse='.7 .7 .7',specular='.1 .1 .1')
 for g in root.findall('.//geom'):
  name=g.get('name','')
  if name=='geom_Root':g.set('rgba','.2 .25 .3 1')
  if name in ('geom_WalkingPath','geom_Treadmill_Belt'):g.set('rgba','.65 .67 .7 1')
  if name=='geom_Treadmill_Belt':g.set('material','belt_material');g.set('rgba','1 1 1 1')
 m=mujoco.MjModel.from_xml_string(E.tostring(root,encoding='unicode'));d=mujoco.MjData(m);renderer=mujoco.Renderer(m,height=700,width=900)
 try:font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',19)
 except OSError:font=ImageFont.load_default()
 frames=[];indices=np.flatnonzero((z['t']>=start)&(z['t']<=stop))[::max(1,round(.08/r['physics_dt']))]
 for k in indices:
  d.qpos[:]=z['q'][k];mujoco.mj_forward(m,d);renderer.update_scene(d,camera='review_oblique')
  if not frames:
   debug=dict(root=d.xpos[1].tolist(),camera_position=d.cam_xpos.tolist(),camera_rotation=d.cam_xmat.tolist(),first_geom_positions=d.geom_xpos[:4].tolist(),extent=float(m.stat.extent),znear=float(m.vis.map.znear),ngeom=int(renderer.scene.ngeom));(folder/'RENDER_CAMERA.json').write_text(json.dumps(debug,indent=2))
  profile=r.get('source_profile',{});clearance=r.get('guide_clearance_m')
  if clearance is None and r.get('source_embedded_guides'):clearance=profile.get('guide_clearance_m')
  guidance=f'physical guides ({1000*clearance:g} mm slack)' if clearance is not None else 'continuous sagittal harness' if r.get('sagittal_harness_scale',0) else 'source rails only'
  circuit=title or (f'W2L {len(z["neuron_names"])}-neuron SNS + Hill muscles')
  source_clock=r.get('source_neural_clock_scale',profile.get('neural_clock_scale',1.));total_clock=source_clock*r.get('neural_time_scale',1.)
  guidance_line=guidance_caption or f'Rail-supported treadmill stepping; {guidance}'
  detail_line=detail_caption or f'Neural time constants x{total_clock:g}; mean net assistance {100*r.get("mean_net_vertical_assistance_BW",r.get("mean_rail_support_BW",0)):.2f}% BW'
  frame=Image.fromarray(renderer.render());draw=ImageDraw.Draw(frame);draw.rectangle([0,0,900,79],fill='white');draw.text((14,5),f'{circuit} | t={z["t"][k]:.2f} s',fill='black',font=font);draw.text((14,29),guidance_line,fill='black',font=font);draw.text((14,52),detail_line,fill='black',font=font);frames.append(frame)
 still_index=len(frames)//2 if still_time is None else int(np.argmin(abs(z['t'][indices]-still_time)))
 frames[0].save(folder/'SUPPORTED_STEPPING_REPLAY.gif',save_all=True,append_images=frames[1:],duration=80,loop=0);frames[still_index].save(folder/'SUPPORTED_STEPPING_STILL.png');del renderer
 provenance={'trajectory_sha256':hashlib.sha256((folder/'trajectory.npz').read_bytes()).hexdigest(),'model_sha256':hashlib.sha256((folder/'model.xml').read_bytes()).hexdigest(),'source_sha256':r.get('source_sha256'),'neuron_count':len(z['neuron_names']),'start_s':float(z['t'][indices[0]]),'stop_s':float(z['t'][indices[-1]]),'frame_count':len(frames),'frame_duration_ms':80,'title':circuit,'guidance_caption':guidance_line,'detail_caption':detail_line,'source_neural_clock_scale':source_clock,'runtime_neural_time_scale':r.get('neural_time_scale',1.),'source_embedded_guides':r.get('source_embedded_guides',False),'guide_clearance_m':clearance,'still_time_s':float(z['t'][indices[still_index]]),'scope':'Replay of recorded qpos only; no integration, optimization, or physical-model change. Visual material and camera changes affect rendering only.'}
 (folder/'RENDER_PROVENANCE.json').write_text(json.dumps(provenance,indent=2))
 print(str(folder/'SUPPORTED_STEPPING_REPLAY.gif'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('folder');p.add_argument('--start',type=float,default=2.);p.add_argument('--stop',type=float,default=4.7);p.add_argument('--title');p.add_argument('--guidance-caption');p.add_argument('--detail-caption');p.add_argument('--still-time',type=float);a=p.parse_args();render(a.folder,a.start,a.stop,a.title,a.guidance_caption,a.detail_caption,a.still_time)
