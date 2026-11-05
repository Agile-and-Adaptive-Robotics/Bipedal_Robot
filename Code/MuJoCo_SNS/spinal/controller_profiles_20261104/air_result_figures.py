"""Plot and render saved air-probe states using the actual compiled model."""
import argparse,hashlib,json,os
from pathlib import Path
for key in('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline',type=Path);p.add_argument('--repeat-source',type=Path);p.add_argument('--gif',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 capture=a.run/'capture';source=capture/'aligned_trace.npz';z=np.load(source);r=json.loads((a.run/'body.result.json').read_text());meta=json.loads((capture/'ALIGNED_TRACE_METADATA.json').read_text());cfg=json.loads((a.run/'body.json').read_text());family=cfg['controller_family']
 if a.repeat_source:
  previous=np.load(a.repeat_source/'capture/aligned_trace.npz');comparison={}
  for key in('t_state','qpos','qvel','neural_voltage','foot_force_anatomical_N','harness_wrench_anatomical'):
   same_shape=previous[key].shape==z[key].shape
   comparison[key]=dict(exact=bool(same_shape and np.array_equal(previous[key],z[key],equal_nan=True)),maximum_absolute_difference=float(np.nanmax(abs(previous[key]-z[key])))if same_shape else None)
  (a.output/'FRESH_REPEAT_AUDIT.json').write_text(json.dumps(dict(previous=str(a.repeat_source),repeat=str(a.run),arrays=comparison,all_exact=all(v['exact']for v in comparison.values())),indent=2))
 if cfg['stage']not in(1,2)or not r['assessment']['complete']:raise ValueError('Complete air recording required')
 t=z['t_state'];tn=z['t_neural'];names=z['neural_names'].tolist();jn=z['joint_names'].tolist();keep=(t>=5)&(t<11);kn=(tn>=5)&(tn<11)
 plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none'})
 fig,ax=plt.subplots(7,1,figsize=(8.4,10),sharex=True,constrained_layout=True)
 cells=[['RG_F_r','RG_E_r'],['PF_S1_r','PF_S2_r']if family=='syn6'else['PF_F1_r','PF_E1_r'],['MN_tib_ant_r','MN_soleus_r']]
 for row,selected in enumerate(cells):
  for i,name in enumerate(selected):
   if name in names:ax[row].plot(tn[kn],z['neural_voltage'][kn,names.index(name)],'-'if i==0 else'--',color=('#4477AA','#EE6677')[i],label=name)
 for row,joint in enumerate(('hip_flexion','knee_angle','ankle_angle'),3):
  for side,color,style in(('r','#4477AA','-'),('l','#EE6677','--')):
   q=np.degrees(z['qpos'][:,int(z['joint_qpos_addresses'][jn.index(joint+'_'+side)])]);q=-q if joint=='knee_angle'else q
   ax[row].plot(t[keep],q[keep],style,color=color,label=side.upper())
 bw=float(meta['body_weight_N'])
 for side,color,style in((0,'#4477AA','-'),(1,'#EE6677','--')):ax[6].plot(t[keep],z['foot_force_anatomical_N'][keep,side,1]/bw,style,color=color,label=('Right','Left')[side]+' foot')
 for row,label in enumerate(('RG (mV)','PF (mV)','MN (mV)','Hip (deg)','Knee flexion\n(deg)','Ankle (deg)','Foot vertical\nforce / BW')):
  ax[row].set_ylabel(label);ax[row].grid(alpha=.2);ax[row].legend(loc='upper right',fontsize=8,ncol=2);ax[row].set_xlim(5,11)
 ax[6].set_ylim(-.01,.05);ax[-1].set_xlabel('Recorded time (s)')
 fig.suptitle(f'Gait2392 {family}: stage {cfg["stage"]}, suspended air stepping\nActual neural and mechanical traces; feet unloaded',fontsize=12)
 for ext in('png','pdf'):fig.savefig(a.output/('LAYERS_JOINTS_LOAD.'+ext),dpi=170)
 plt.close(fig)
 if a.baseline:
  b=np.load(a.baseline/'capture/aligned_trace.npz');fig,axes=plt.subplots(2,1,figsize=(8.4,4.4),sharex=True,constrained_layout=True)
  from scipy.signal import periodogram
  spectra=[]
  for row,side in enumerate(('r','l')):
   for data,label,color,style in((b,'Before','#AA3377','--'),(z,'Candidate','#228833','-')):
    bt=data['t_state'];mask=(bt>=5)&(bt<11);ids=data['joint_names'].tolist();col=int(data['joint_qpos_addresses'][ids.index('ankle_angle_'+side)])
    axes[row].plot(bt[mask],np.degrees(data['qpos'][mask,col]),style,color=color,lw=1,label=label)
    audit=(bt>=5)&(bt<14);frequency,power=periodogram(np.degrees(data['qpos'][audit,col]),fs=1/np.median(np.diff(bt)),window='hann');df=frequency[1]-frequency[0]
    spectra.append(dict(side=side,condition=label,power_above3Hz_fraction=float(power[frequency>3].sum()/max(power[1:].sum(),1e-30)),high_frequency_RMS_deg=float(np.sqrt(power[frequency>3].sum()*df)),scope='5-14s Hann periodogram; high-frequency RMS and fraction, not proof of oscillator mechanism'))
   axes[row].set_ylabel(side.upper()+' ankle (deg)');axes[row].grid(alpha=.2);axes[row].legend(ncol=2)
  axes[-1].set_xlabel('Recorded time (s)');fig.suptitle('syn6 ankle waveform: matched mechanics, mirrored neural changes')
  for ext in('png','pdf'):fig.savefig(a.output/('ANKLE_COMPARISON.'+ext),dpi=170)
  plt.close(fig)
  (a.output/'ANKLE_SPECTRA.json').write_text(json.dumps(spectra,indent=2))
 if a.gif:
  import mujoco
  from PIL import Image,ImageDraw,ImageFont
  modelpath=capture/meta['model_file']
  if hashlib.sha256(modelpath.read_bytes()).hexdigest()!=meta['model_sha256']:raise ValueError('Compiled model hash mismatch')
  m=mujoco.MjModel.from_binary_path(str(modelpath));d=mujoco.MjData(m);m.vis.global_.offwidth=900;m.vis.global_.offheight=800
  # Air-stage floor collision is disabled. Hide its visual plane as well;
  # otherwise saved poses appear to penetrate a supporting ground surface.
  if not r['measured']['ground_enabled']:
   for geom in range(m.ngeom):
    if m.geom_type[geom]==mujoco.mjtGeom.mjGEOM_PLANE:m.geom_rgba[geom,3]=0.
  for body in range(m.nbody):
   name=m.body(body).name or''
   if name.endswith(('_r','_l')):
    for geom in range(m.body_geomadr[body],m.body_geomadr[body]+m.body_geomnum[body]):m.geom_rgba[geom]=[.27,.47,.67,1.]if name.endswith('_r')else[.93,.47,.20,1.]
  renderer=mujoco.Renderer(m,height=540,width=400);opt=mujoco.MjvOption();opt.flags[mujoco.mjtVisFlag.mjVIS_TENDON]=False;font=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',17);frames=[];stride=max(1,round(1/(12*np.median(np.diff(t)))));indices=np.flatnonzero(keep)[::stride]
  for i in indices:
   d.qpos[:]=z['qpos'][i];d.qvel[:]=z['qvel'][i];mujoco.mj_forward(m,d);frame=Image.new('RGB',(800,600),'white')
   for column,azimuth in enumerate((180.,140.)):
    cam=mujoco.MjvCamera();cam.lookat[:]=[d.xpos[m.body('pelvis').id,0],d.xpos[m.body('pelvis').id,1],.85];cam.distance=2.8;cam.azimuth=azimuth;cam.elevation=-6;renderer.update_scene(d,camera=cam,scene_option=opt);frame.paste(Image.fromarray(renderer.render().copy()),(400*column,60))
   draw=ImageDraw.Draw(frame);draw.text((10,4),f'Gait2392 {family} | suspended air stepping | 1x playback',font=font,fill='black');draw.text((10,28),f'Ground contact disabled | feet unloaded | t={t[i]:.2f} s',font=font,fill='black');frames.append(frame)
  del renderer
  ticks=np.rint(np.r_[0.,np.cumsum(np.r_[np.diff(t[indices]),np.median(np.diff(t[indices]))])]*100).astype(int)
  frames[0].save(a.output/'AIR_STEPPING.gif',save_all=True,append_images=frames[1:],duration=(np.diff(ticks)*10).tolist(),loop=0);frames[len(frames)//2].save(a.output/'AIR_STEPPING.png')
 (a.output/'PROVENANCE.json').write_text(json.dumps(dict(run=str(a.run),capture_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),baseline=str(a.baseline)if a.baseline else None,scope='Saved neural and physical samples; actual compiled model for rendering; no resimulation or pose corrections; air only'),indent=2))
if __name__=='__main__':main()
