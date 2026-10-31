"""Offscreen replay only: y-up camera, visible support, no image synthesis."""
import argparse,json
from pathlib import Path
import xml.etree.ElementTree as E
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from import_virtual_walker import mujoco,fmt
def render(folder,start=2.,stop=4.7):
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
  guidance=f'physical guides ({1000*r["guide_clearance_m"]:g} mm slack)' if r.get('guide_clearance_m') is not None else 'continuous sagittal harness' if r.get('sagittal_harness_scale',0) else 'source rails only'
  circuit=f'W2L{len(z["neuron_names"])} diagnostic: tonic {r.get("symmetric_rg_ext_tonic_nA")}, V3 output {r["v3_output_scale"]:g}' if r.get('v3_output_scale',1)!=1 else 'W2L source circuit + Hill muscles'
  frame=Image.fromarray(renderer.render());draw=ImageDraw.Draw(frame);draw.rectangle([0,0,900,79],fill='white');draw.text((14,5),f'{circuit} | t={z["t"][k]:.2f} s',fill='black',font=font);draw.text((14,29),f'Rail-supported treadmill stepping; {guidance}',fill='black',font=font);draw.text((14,52),f'Neural time scale {r["neural_time_scale"]:g}; mean vertical assistance {100*r.get("mean_net_vertical_assistance_BW",r.get("mean_rail_support_BW",0)):.1f}% BW',fill='black',font=font);frames.append(frame)
 frames[0].save(folder/'SUPPORTED_STEPPING_REPLAY.gif',save_all=True,append_images=frames[1:],duration=80,loop=0);frames[len(frames)//2].save(folder/'SUPPORTED_STEPPING_STILL.png');del renderer
 print(str(folder/'SUPPORTED_STEPPING_REPLAY.gif'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('folder');p.add_argument('--start',type=float,default=2.);p.add_argument('--stop',type=float,default=4.7);a=p.parse_args();render(a.folder,a.start,a.stop)
