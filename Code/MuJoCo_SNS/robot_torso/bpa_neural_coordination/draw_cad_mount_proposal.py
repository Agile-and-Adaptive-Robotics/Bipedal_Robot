"""Actual CAD rendering with candidate sleeve envelopes; no simulated standing claim."""
import os,sys,json,xml.etree.ElementTree as ET
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent;OUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mount_subset')

def point_segment(p,a,b):
 delta=b-a;t=np.clip((p-a)@delta/(delta@delta),0,1);return float(np.linalg.norm(p-(a+t*delta)))

def segment_distance(a,b,c,d):
 u=b-a;v=d-c;w=a-c;A=u@u;B=u@v;C=v@v;D=u@w;E=v@w;den=A*C-B*B
 values=[point_segment(a,c,d),point_segment(b,c,d),point_segment(c,a,b),point_segment(d,a,b)]
 if den>1e-15:
  s=(B*E-C*D)/den;t=(A*E-B*D)/den
  if 0<=s<=1 and 0<=t<=1:values.append(float(np.linalg.norm(w+s*u-t*v)))
 return min(values)

def main():
 proposal=json.loads((OUT/'PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json').read_text());tree=ET.parse(ROBOT/'ah001_neural.xml');root=tree.getroot();root.find('compiler').set('meshdir',str(ROBOT/'meshes'));ET.SubElement(root.find('asset'),'texture',dict(name='proposal_white_background',type='skybox',builtin='flat',rgb1='1 1 1',rgb2='1 1 1',width='32',height='32'));m=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'));d=mujoco.MjData(m);q=np.load(proposal['source_motion'])['qpos'];mounts=[]
 for row in proposal['routes']:
  if not row['sleeve_placement_certified']:continue
  mounts.append(dict(name=row['name'],body=m.body(row['body']).id,local=np.array(row['selected_mount']['sleeve_axis_local_m']),radius=row['source_inflated_radius_m'],color=np.array([.835,.369,0,1.]) if row['name'].startswith('bifem') else np.array([0,.447,.698,1.])))
 minima={}
 for index in proposal['sampled_pose_indices']:
  d.qpos[:]=q[index];mujoco.mj_forward(m,d);world=[]
  for mount in mounts:world.append(mount['local']@d.xmat[mount['body']].reshape(3,3).T+d.xpos[mount['body']])
  for i in range(len(mounts)):
   for j in range(i):
    key=mounts[j]['name']+' : '+mounts[i]['name'];gap=segment_distance(*world[j],*world[i])-mounts[j]['radius']-mounts[i]['radius'];minima[key]=min(minima.get(key,float('inf')),gap)
 pair_report=dict(sampled_poses=17,surface_gap_m=minima,all_sleeve_pairs_clear=all(v>0 for v in minima.values()),adopted_simultaneous_sleeve_subset=['vas_med_r','vas_med_l'],rejected_simultaneous_flexor_pair='bifemsh_r/l overlap by18.07mm during the sampled motion; body-only separation was insufficient',scope='Mutual sleeve capsule envelopes only, at saved poses. Fittings, cables, and swept motion remain unverified.')
 (OUT/'SLEEVE_PAIR_CLEARANCE.json').write_text(json.dumps(pair_report,indent=2),encoding='utf-8')
 mounts=[r for r in mounts if r['name'] in pair_report['adopted_simultaneous_sleeve_subset']]
 d.qpos[:]=q[0];mujoco.mj_forward(m,d)
 for g in range(m.ngeom):
  if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH:m.geom_rgba[g]=[.55,.58,.6,1.]
 option=mujoco.MjvOption();option.flags[mujoco.mjtVisFlag.mjVIS_TENDON]=False;renderer=mujoco.Renderer(m,width=1000,height=900);images=[]
 for position in ([2.55,.7,0],[2.2,1.05,1.8]):
  renderer.update_scene(d,scene_option=option);scene=renderer.scene
  for mount in mounts:
   points=mount['local']@d.xmat[mount['body']].reshape(3,3).T+d.xpos[mount['body']];delta=points[1]-points[0];length=np.linalg.norm(delta);axis=delta/length;ref=np.array([1.,0,0]);ref=ref-axis*(ref@axis);ref/=np.linalg.norm(ref);rotation=np.column_stack((ref,np.cross(axis,ref),axis))
   mujoco.mjv_initGeom(scene.geoms[scene.ngeom],mujoco.mjtGeom.mjGEOM_CAPSULE,np.array([mount['radius'],length/2,0]),points.mean(0),rotation.ravel(),mount['color']);scene.ngeom+=1
  pos=np.array(position);forward=np.array([0,.62,0])-pos;forward/=np.linalg.norm(forward);right=np.cross(forward,[0,1,0]);right/=np.linalg.norm(right);up=np.cross(right,forward)
  for camera in scene.camera:camera.pos[:]=pos;camera.forward[:]=forward;camera.up[:]=up
  images.append(renderer.render().copy())
 if hasattr(renderer,'close'):renderer.close()
 del renderer
 plt.rcParams.update({'font.family':'Arial','font.size':17});fig,axes=plt.subplots(1,2,figsize=(11,6.2))
 for ax,image,title in zip(axes,images,('(A) Anterior','(B) Oblique')):ax.imshow(image);ax.axis('off');ax.set_title(title,loc='left')
 fig.legend(handles=[Line2D([],[],color=mounts[0]['color'][:3],lw=8,label='Vastus-medialis-route sleeve envelope')],loc='lower center',ncol=1,frameon=False,fontsize=17);fig.subplots_adjust(left=.01,right=.99,top=.94,bottom=.1,wspace=.01)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/('CAD_SLEEVE_MOUNT_PROPOSAL.'+ext),dpi=220)
 plt.close(fig);(OUT/'CAD_SLEEVE_MOUNT_PROPOSAL.alt.txt').write_text('Actual AH001 CAD geometry at the saved initial pose, with two blue conservative sleeve envelopes placed from measured mesh bounds for the existing vas_med route pair. Paired sleeves share dimensions. The initially considered flexor pair is excluded because its envelopes intersect during sampled motion. This is an external remote-sleeve mounting proposal, not a pressure simulation or built robot. Cables, Bowden housing, fittings, brackets, and loads remain unverified. The body pose is from an older hybrid/reference-tension recording and is used only as geometry.',encoding='utf-8');print(json.dumps(pair_report,indent=2))

if __name__=='__main__':main()
