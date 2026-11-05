from pathlib import Path
import xml.etree.ElementTree as E
import json,hashlib,math
import numpy as np
import mujoco
from PIL import Image,ImageDraw,ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]
REPO=Path(__file__).resolve().parents[7]
W=BASE/'build_fixed_20261105'/'visual4'
W.mkdir(parents=True,exist_ok=True)
import shutil
shutil.copy2(BASE/'data'/'human_robot_rom'/'torque_source.csv',W/'torque_source.csv')
SOURCE=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
table=np.genfromtxt(W/'torque_source.csv',delimiter=',',names=True)
qlo=float(table['knee_deg'][0]);qhi=float(table['knee_deg'][-1])
positive=np.flatnonzero(np.abs(table['robot_Nm'])>1e-10)
blo=float(table['knee_deg'][positive[0]]);bhi=float(table['knee_deg'][positive[-1]])
assert np.all(np.diff(positive)==1)
root=E.parse(SOURCE).getroot()
for mesh in root.findall('./asset/mesh'):
 mesh.set('file',str(SOURCE.parent/mesh.get('file')))
for tex in root.findall('./asset/texture'):
 if tex.get('type')=='skybox':tex.set('rgb1','1 1 1');tex.set('rgb2','1 1 1')
pel=root.find('.//body[@name="pelvis"]');quat=np.fromstring(pel.get('quat'),sep=' ');quat/=np.linalg.norm(quat)
R=np.empty(9);mujoco.mju_quat2Mat(R,quat);R=R.reshape(3,3)
target=np.array([-.195,.515,.09]);distance=2.25
position=R@(target+np.array([0,0,distance]))
E.SubElement(root.find('worldbody'),'camera',name='OpenSim_plusZ',pos=' '.join(map(str,position)),
             xyaxes=' '.join(map(str,np.r_[R[:,0],R[:,1]])),fovy='34',ipd='0')
model=mujoco.MjModel.from_xml_string(E.tostring(root,encoding='unicode'))
model.vis.global_.offwidth=480;model.vis.global_.offheight=780
model.opt.gravity[:]=0
data=mujoco.MjData(model);renderer=mujoco.Renderer(model,height=780,width=480)
option=mujoco.MjvOption();option.flags[mujoco.mjtVisFlag.mjVIS_TENDON]=True
option.sitegroup[:]=0
keep={'pelvis','femur_r','tibia_r','talus_r','calcn_r','toes_r'}
for k in range(model.ngeom):
 body=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_BODY,model.geom_bodyid[k])
 if body not in keep:model.geom_rgba[k,3]=0
 else:model.geom_rgba[k]=[.66,.68,.70,1]
selected=[]
for k in range(model.ntendon):
 name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_TENDON,k) or ''
 model.tendon_rgba[k]=[.1,.1,.1,0]
 if 'bifemlh_r' in name:selected.append(k)
assert len(selected)==1,selected
model.tendon_width[:]=.004
cam=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_CAMERA,'OpenSim_plusZ')
def assign(name,v):
 j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
 assert j>=0,name
 data.qpos[model.jnt_qposadr[j]]=v
def pose(deg):
 data.qpos[:]=model.qpos0
 assign('pelvis_ty',.96);assign('hip_flexion_r',0);assign('knee_angle_r',math.radians(deg))
 for _ in range(2):
  for k in range(model.neq):
   if model.eq_type[k]!=mujoco.mjtEq.mjEQ_JOINT:continue
   j1,j2=model.eq_obj1id[k],model.eq_obj2id[k]
   x=data.qpos[model.jnt_qposadr[j2]]-model.qpos0[model.jnt_qposadr[j2]] if j2>=0 else 0
   data.qpos[model.jnt_qposadr[j1]]=model.qpos0[model.jnt_qposadr[j1]]+np.polynomial.polynomial.polyval(x,model.eq_data[k,:5])
 mujoco.mj_forward(model,data)
 outward=R.T@data.cam_xmat[cam].reshape(3,3)[:,2]
 assert np.linalg.norm(outward-[0,0,1])<1e-8,outward
 assert abs(data.qpos[model.jnt_qposadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'hip_flexion_r')]])<1e-12
 a=model.body('femur_r').id;b=model.body('tibia_r').id
 assert .37<np.linalg.norm(data.xpos[b]-data.xpos[a])<.44
def render(deg,color):
 pose(deg);model.tendon_rgba[selected[0]]=[*color,1]
 renderer.update_scene(data,camera='OpenSim_plusZ',scene_option=option)
 renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
 renderer.scene.flags[mujoco.mjtRndFlag.mjRND_REFLECTION]=False
 return Image.fromarray(renderer.render())

# Exact, unchanged graph arrays; a new compact layout fits between the models.
plt.rcParams.update({'font.family':'Arial','font.size':13})
fig,ax=plt.subplots(figsize=(4.7,3.35),dpi=220)
fig.subplots_adjust(left=.17,right=.98,bottom=.22,top=.94)
ax.plot(table['knee_deg'],abs(table['human_Nm']),'k--',lw=2.4,label='Human muscle')
ax.plot(table['knee_deg'],abs(table['robot_Nm']),color='#CD34B5',lw=2.4,label='BPA route')
ax.axvspan(blo,bhi,color='#CD34B5',alpha=.07)
ax.axvline(blo,color='#CD34B5',ls=':',lw=1.4)
ax.set(xlim=(-120,10),ylim=(0,60),xlabel='Knee angle (degrees)',ylabel='Torque magnitude (N·m)')
ax.set_xticks([-120,-80,-40,0]);ax.spines[['top','right']].set_visible(False)
ax.legend(loc='upper right',frameon=False,fontsize=12)
fig.savefig(W/'torque_comparison.png',facecolor='white');plt.close(fig)

N=120;u=np.linspace(0,1,61);down=qhi+(qlo-qhi)*(.5-.5*np.cos(np.pi*u));q=np.r_[down,down[-2:0:-1]]
qb=np.clip(q,blo,bhi)
preview=Image.new('RGB',(960,780),'white')
preview.paste(render(qlo,(.08,.08,.10)),(0,0));preview.paste(render(blo,(.804,.204,.710)),(480,0));preview.save(W/'models_preview.png')
bounds=[]
for angle in np.linspace(qlo,qhi,21):
 pose(float(angle))
 for gid in range(model.ngeom):
  if model.geom_rgba[gid,3]<=0 or model.geom_type[gid]!=mujoco.mjtGeom.mjGEOM_MESH:continue
  mid=model.geom_dataid[gid];start=model.mesh_vertadr[mid];count=model.mesh_vertnum[mid]
  verts=model.mesh_vert[start:start+count]@data.geom_xmat[gid].reshape(3,3).T+data.geom_xpos[gid]
  native=verts@R
  bounds.append([native.min(0),native.max(0)])
print('Visible native bounds:',np.min(np.array(bounds)[:,0,:],axis=0),np.max(np.array(bounds)[:,1,:],axis=0),flush=True)
manifest={'source_opensim_model':'gait2392_simbody.osim','render_model':str(SOURCE),'render_model_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'method':'Prescribed coordinate playback of the OpenSim Gait2392 model converted by MyoConverter; no dynamics integration.',
          'camera_side_in_opensim_frame':[0,0,1],'camera_up_in_opensim_frame':[0,1,0],'hip_degrees':0,
          'human_range_degrees':[qlo,qhi],'bpa_range_degrees':[blo,bhi],'last_zero_torque_sample_degrees':float(table['knee_deg'][positive[0]-1]),
          'right_motion':'Clamp synchronized human angle to the nonzero-BPA-torque sample interval.','bodyweight_used':False,
          'source_curve_sha256':hashlib.sha256((W/'torque_source.csv').read_bytes()).hexdigest(),'frames':N,'duration_ms':8000}
(W/'rom_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
if '--preview' in __import__('sys').argv:print(json.dumps(manifest,indent=2));raise SystemExit

left=[];right=[]
for k,(a,b) in enumerate(zip(q,qb)):
 left.append(render(float(a),(.08,.08,.10)));right.append(render(float(b),(.804,.204,.710)))
 if k%30==0:print('Rendered matched pair',k,flush=True)
durations=[(60,70,70)[k%3] for k in range(N)]
for name,images in [('human_full_rom.gif',left),('bpa_limited_rom.gif',right)]:
 pal=images[len(images)//4].quantize(colors=256,method=Image.MEDIANCUT)
 ims=[im.quantize(palette=pal,dither=Image.FLOYDSTEINBERG) for im in images]
 ims[0].save(W/name,save_all=True,append_images=ims[1:],loop=0,duration=durations,disposal=2,optimize=False)

# One composite preview guarantees a synchronized view outside PowerPoint too.
graph_source=Image.open(W/'torque_comparison.png').convert('RGB')
font=ImageFont.truetype(r'C:\Windows\Fonts\arialbd.ttf',36)
small=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',30)
paired=[]
for k,(a,b) in enumerate(zip(left,right)):
 graph=graph_source.copy();pen=ImageDraw.Draw(graph)
 for angle,values,color in [(q[k],table['human_Nm'],'#111111'),(qb[k],table['robot_Nm'],'#CD34B5')]:
  torque=abs(float(np.interp(angle,table['knee_deg'],values)))
  x=(.17+.81*(angle+120)/130)*graph.width
  y=(1-(.22+.72*torque/60))*graph.height
  pen.ellipse((x-13,y-13,x+13,y+13),fill=color,outline='white',width=3)
 graph.thumbnail((1000,760),Image.LANCZOS)
 canvas=Image.new('RGB',(1960,780),'white');canvas.paste(a,(0,0));canvas.paste(b,(1480,0));canvas.paste(graph,(480+(1000-graph.width)//2,(780-graph.height)//2))
 paired.append(canvas)
pal=paired[N//4].quantize(colors=256,method=Image.MEDIANCUT)
ims=[im.quantize(palette=pal,dither=Image.FLOYDSTEINBERG) for im in paired]
ims[0].save(W/'human_robot_rom_comparison.gif',save_all=True,append_images=ims[1:],loop=0,duration=durations,disposal=2,optimize=False)
left[0].save(W/'human_poster.png');right[0].save(W/'bpa_poster.png')
np.savez_compressed(W/'prescribed_angles.npz',human=q,bpa=qb,duration_ms=durations)
if hasattr(renderer,'close'):renderer.close()
print('MATCHED_ROM_GIFS_COMPLETE',flush=True)
