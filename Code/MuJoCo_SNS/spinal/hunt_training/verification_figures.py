"""Few-cycle figures and pose-only replay from one actual verification recording.

No neuronal or physics integration is performed. Missing channels remain missing.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np

def read_run(source):
 p=Path(source).resolve()
 if p.is_dir():p=p/'aligned_trace.npz'
 z=np.load(p,allow_pickle=False)
 if 't_state'in z:
  meta=json.loads((p.parent/'ALIGNED_TRACE_METADATA.json').read_text(encoding='utf-8'))
  out=dict(source=p,metadata=meta,t=z['t_state'],tn=z['t_neural'],names=z['neural_names'].tolist(),
   voltage=z['neural_voltage'],joint_names=z['joint_names'].tolist(),
   q=z['qpos'][:,z['joint_qpos_addresses']],foot=z['foot_force_anatomical_N'],
   regional=z['regional_vertical_N'],harness=z['harness_wrench_anatomical'],
   support=z['explicit_support_anatomical_N'],tf=z['force_interval_end'],stim=z['stimulus_nA'],
   ts=z['t_neural_input'],archive=z)
 else:
  # Historical sidecar adapter: preserve the logged logical clocks and state
  # their pre/post-step limitation instead of silently time-shifting samples.
  meta={'schema':'legacy_hunt_sidecars','voltage_label':'model voltage (mV; recorded normalization)',
   'native_frame':'legacy_z_up','body_weight_N':None,
   'sample_phase':'Legacy sidecars use logical timestamps. Neural/force samples precede physics; qpos is post-step. Do not use for substep latency inference.',
   'regional_sensor_definition':'No regional sensor trace unless a separate measured sidecar exists.'}
  names=z['neuro_names'].tolist();V=z['neuro'].copy();tn=z['t'];npz=p.with_suffix('.neurons.npz')
  if npz.exists():
   n=np.load(npz,allow_pickle=False);assert np.array_equal(tn,n['t'])
   for i,name in enumerate(n['cells'].tolist()):
    if name not in names:names.append(name);V=np.column_stack((V,n['voltage'][:,i]))
  ass=np.load(p.with_suffix('.assistance.npz'),allow_pickle=False);a=ass['samples'];cols=ass['columns'].tolist();col=lambda n:a[:,cols.index(n)]
  nativeF=np.array([[col(s+'_GRFx_N'),col(s+'_GRFy_N'),col(s+'_GRFz_N')]for s in('right','left')]).transpose(2,0,1)
  rot=np.array([[1,0,0],[0,0,1],[0,-1,0]])
  foot=nativeF@rot.T;w=np.column_stack([col('harness_F'+s+'_N')for s in('x','y','z')])@rot.T
  torque=np.column_stack([col('harness_M'+s+'_Nm')for s in('x','y','z')])@rot.T
  support=np.zeros((len(a),3));support[:,1]=col('explicit_support_N');regional=np.full((len(a),2,2),np.nan)
  rp=p.with_suffix('.regional_contacts.npz')
  if rp.exists():
   rz=np.load(rp,allow_pickle=False);assert np.array_equal(rz['t'],a[:,0]);regional=rz['force_N'].reshape(-1,2,2)
   meta['regional_sensor_definition']='Legacy hind/forefoot contact-normal magnitudes; vertical only on a horizontal plane.'
  out=dict(source=p,metadata=meta,t=z['t'],tn=tn,names=names,voltage=V,joint_names=z['key_joints'].tolist(),
   q=np.radians(z['q']),foot=foot,regional=regional,harness=np.column_stack((w,torque)),support=support,
   tf=a[:,0],stim=None,ts=tn,archive=z)
 if not np.all(np.diff(out['t'])>0):raise ValueError('Nonmonotonic/truncated state record; explicit repair/valid-length metadata required')
 return out

def choose_window(run,start=None,end=None,cycles=3):
 from scipy.signal import find_peaks
 lo=max(float(run[k][0])for k in('t','tn','tf'));hi=min(float(run[k][-1])for k in('t','tn','tf'))
 end=hi if end is None else min(float(end),hi);period=None
 if 'RG_E_r'in run['names']and'RG_F_r'in run['names']:
  y=run['voltage'][:,run['names'].index('RG_E_r')]-run['voltage'][:,run['names'].index('RG_F_r')]
  peaks,_=find_peaks(y,prominence=max(.2,.15*np.ptp(y)),distance=max(1,round(.25/np.median(np.diff(run['tn'])))))
  times=run['tn'][peaks];times=times[(times>=lo)&(times<=end)]
  if len(times)>=3:period=float(np.median(np.diff(times[-5:])))
 start=max(lo,end-(cycles*period if period else 3.))if start is None else max(float(start),lo)
 if not start<end:raise ValueError('Empty common verification window')
 return start,end,period

def plot(run,output,label,start,end,period):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.backends.backend_pdf import PdfPages
 plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
 output=Path(output);output.mkdir(parents=True,exist_ok=True);missing=[];colors=['#0072B2','#D55E00','#009E73','#CC79A7']
 figures=[];neural,axes=plt.subplots(5,2,figsize=(6.5,10.2),sharex=True)
 for si,s in enumerate(('r','l')):
  choices=[['RG_E_'+s,'RG_F_'+s],
   [x+s for x in('PF_BI_JOINT_RECT_','PF_BI_JOINT_HAM_','PF_IN_BI_PAIR_RECT_','PF_IN_BI_PAIR_HAM_')],
   ['MN_iliacus_'+s,'MN_bifemlh_'+s],
   [next((n for n in('MN_vas_int_'+s,'MN_vas_lat_'+s)if n in run['names']),'MN_vas_int_'+s),'MN_bifemsh_'+s],
   ['MN_soleus_'+s,'MN_tib_ant_'+s]]
  for row,names in enumerate(choices):
   ax=axes[row,si];found=0
   for j,name in enumerate(names):
    if name not in run['names']:missing.append(name);continue
    mask=(run['tn']>=start)&(run['tn']<=end)
    display=name.removesuffix('_'+s).removeprefix('MN_')
    display={'RG_E':'E','RG_F':'F','PF_BI_JOINT_RECT':'RF PF','PF_BI_JOINT_HAM':'HAM PF','PF_IN_BI_PAIR_RECT':'RF relay','PF_IN_BI_PAIR_HAM':'HAM relay'}.get(display,display)
    ax.plot(run['tn'][mask],run['voltage'][mask,run['names'].index(name)],color=colors[j],ls='--'if j>=2 else'-',lw=1.,label=display);found+=1
   ax.set_title(('Right','Left')[si]+' '+('RG','RF/HAM four-cell PF','Hip MNs','Knee MNs','Ankle MNs')[row],loc='left',fontsize=10)
   ax.set_ylabel('Model mV');ax.grid(alpha=.2)
   if found:ax.legend(loc='upper right',ncol=2,fontsize=10,framealpha=.8)
   else:ax.text(.5,.5,'Not recorded',transform=ax.transAxes,ha='center',va='center');ax.set_yticks([])
   ax.set_xlim(start,end)
 for ax in axes[-1]:ax.set_xlabel('Recorded time (s)')
 neural.suptitle(label+'\n'+run['metadata']['voltage_label'],fontsize=12);neural.tight_layout(rect=(0,0,1,.95));figures.append(('ALIGNED_NEURAL',neural))
 mechanics,axes=plt.subplots(5,2,figsize=(6.5,10.2),sharex=True);jointstats={};BW=run['metadata'].get('body_weight_N');unit='BW'if BW else'N';scale=BW or 1.
 for si,s in enumerate(('r','l')):
  for row,joint in enumerate(('hip_flexion','knee_angle','ankle_angle')):
   ax=axes[row,si];name=joint+'_'+s;mask=(run['t']>=start)&(run['t']<=end)
   if name not in run['joint_names']:missing.append(name);ax.text(.5,.5,'Not recorded',transform=ax.transAxes,ha='center');continue
   angle=np.degrees(run['q'][:,run['joint_names'].index(name)]);shown=-angle if joint=='knee_angle'else angle
   ax.plot(run['t'][mask],shown[mask],color=colors[si],lw=1.2)
   ax.set_title(('Right','Left')[si]+' '+('hip coordinate','knee flexion (+flexion)','ankle coordinate')[row],loc='left',fontsize=10);ax.set_ylabel('Angle (deg)');ax.grid(alpha=.2)
   jointstats[name]={'native_min_deg':float(angle.min()),'native_max_deg':float(angle.max()),'window_min_deg':float(angle[mask].min()),'window_max_deg':float(angle[mask].max())}
   if joint=='knee_angle':
    ax.axhline(0,color='.35',ls=':',lw=.8);jointstats[name]['strict_hyperextension_observed']=bool(np.any(angle>0.));jointstats[name]['maximum_hyperextension_deg']=max(0.,float(angle.max()))
  ax=axes[3,si];mask=(run['tf']>=start)&(run['tf']<=end);found=False
  for region in range(2):
   y=run['regional'][:,si,region]
   if np.isfinite(y).any():ax.plot(run['tf'][mask],y[mask]/scale,color=colors[region],ls=('-','--')[region],label=('Hindfoot','Forefoot')[region],lw=1.);found=True
  ax.plot(run['tf'][mask],run['foot'][mask,si,1]/scale,color='.2',lw=.9,label='Whole-foot superior')
  if not found:ax.text(.02,.08,'Regional forces not recorded',transform=ax.transAxes,fontsize=10,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
  ax.set_title(('Right','Left')[si]+' foot loading',loc='left',fontsize=10);ax.set_ylabel('Force ('+unit+')');ax.legend(ncol=1,fontsize=10);ax.grid(alpha=.2)
 mask=(run['tf']>=start)&(run['tf']<=end)
 for j in range(3):axes[4,0].plot(run['tf'][mask],run['harness'][mask,j]/scale,color=colors[j],label=('Anterior x','Superior y','Right z')[j],lw=1.,ls=('-','--',':')[j])
 axes[4,0].set_title('Measured pelvis restraint force',loc='left',fontsize=10);axes[4,0].set_ylabel('Force ('+unit+')');axes[4,0].legend(fontsize=10,ncol=1)
 for si,s in enumerate(('r','l')):axes[4,1].plot(run['tf'][mask],run['foot'][mask,si,1]/scale,color=colors[si],label=('Right foot','Left foot')[si],lw=1.)
 axes[4,1].plot(run['tf'][mask],run['harness'][mask,1]/scale,color='.4',ls='--',label='Restraint superior')
 if np.isfinite(run['support']).any():axes[4,1].plot(run['tf'][mask],run['support'][mask,1]/scale,color='.2',ls=':',label='Explicit support')
 axes[4,1].set_title('Foot loads and assistance',loc='left',fontsize=10);axes[4,1].set_ylabel('Force ('+unit+')');axes[4,1].legend(fontsize=10,ncol=1)
 for ax in axes[-1]:ax.set_xlabel('Interval end (s)'if run['metadata']['schema']=='aarl_aligned_verification_v1'else'Logical time (s)');ax.grid(alpha=.2)
 for ax in axes.ravel():ax.set_xlim(start,end)
 mechanics.suptitle(label+'\nAnatomical display: x anterior, y superior, z right',fontsize=12);mechanics.tight_layout(rect=(0,0,1,.95));figures.append(('ALIGNED_MECHANICS',mechanics))
 if run['stim']is not None and np.isfinite(run['stim']).any()and np.nanmax(np.abs(run['stim']))>0:
  fig,ax=plt.subplots(figsize=(11,2.8));mask=(run['ts']>=start)&(run['ts']<=end)
  for j in range(run['stim'].shape[1]):
   if np.isfinite(run['stim'][:,j]).any()and np.nanmax(abs(run['stim'][:,j]))>0:ax.plot(run['ts'][mask],run['stim'][mask,j],label=run['names'][j])
  ax.set(xlabel='Input application time (s)',ylabel='Applied current (nA)',title='Actually logged stimuli');ax.legend(fontsize=10,ncol=3);fig.tight_layout();figures.append(('ACTUAL_STIMULI',fig))
 with PdfPages(output/'ALIGNED_VERIFICATION.pdf')as pdf:
  for name,fig in figures:fig.savefig(output/(name+'.png'),dpi=180,bbox_inches='tight',pad_inches=.16);fig.savefig(output/(name+'.pdf'),bbox_inches='tight',pad_inches=.16);pdf.savefig(fig,bbox_inches='tight',pad_inches=.16);plt.close(fig)
 report={'source':str(run['source']),'source_sha256':hashlib.sha256(run['source'].read_bytes()).hexdigest(),'label':label,
  'window_s':[start,end],'estimated_RG_period_s':period,'missing_display_channels':sorted(set(missing)),
  'metadata':run['metadata'],'neural_channels':run['names'],'joint_measurements':jointstats,'all_recorded_neural_values_finite':bool(np.isfinite(run['voltage']).all()),
  'acceptance_policy':'Figures report actual recording; no curriculum acceptance is inferred.',
  'stimulus_policy':'Only explicitly logged nonzero stimulus arrays are plotted; none fabricated.'}
 (output/'FIGURE_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def render(run,output,label,start,end,fps=20,still_only=False,presentation=None):
 import mujoco
 from PIL import Image,ImageDraw,ImageFont
 meta=run['metadata'];source=run['source'];modelpath=source.parent/meta.get('model_file','model.mjb')
 if meta['schema']!='aarl_aligned_verification_v1'or not modelpath.exists():raise ValueError('Replay requires the exact archived compiled model; legacy source-model rebuilding is refused')
 if hashlib.sha256(modelpath.read_bytes()).hexdigest()!=meta['model_sha256']:raise ValueError('Compiled model hash mismatch')
 if meta['mujoco_version']!=mujoco.__version__:raise ValueError('Compiled replay requires its recorded MuJoCo version')
 if presentation is not None:
  if presentation['trace_sha256']!=hashlib.sha256(source.read_bytes()).hexdigest()or presentation['compiled_model_sha256']!=meta['model_sha256']:
   raise ValueError('Presentation support measurements do not belong to this exact saved trace/model')
 m=mujoco.MjModel.from_binary_path(str(modelpath));d=mujoco.MjData(m);z=run['archive'];target=np.array([(start+end)/2])if still_only else np.arange(start,end,1/fps);ids=np.unique(np.searchsorted(run['t'],target));ids=ids[ids<len(run['t'])]
 if len(ids)>200:raise ValueError('Choose a shorter few-cycle replay (at most200frames)')
 m.vis.global_.offwidth=960;m.vis.global_.offheight=640;renderer=mujoco.Renderer(m,height=520,width=440);frames=[];presentation_frames=[];font=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',16)
 if presentation is not None:
  clean_title=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',18);clean_text=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',16)
  assistance=100*presentation['support']['net_external_superior']['mean_fraction_bodyweight'];horizontal=presentation['horizontal']['restraint']['rms_magnitude_N']
 display=np.eye(3)if meta['native_frame']=='legacy_z_up'else np.array([[1,0,0],[0,0,-1],[0,1,0]])
 hide_ground=meta.get('provenance',{}).get('configuration',{}).get('stage')in(1,2)
 for i in ids:
  d.qpos[:]=z['qpos'][i];d.qvel[:]=z['qvel'][i];d.act[:]=z['activation'][i];d.ctrl[:]=z['control'][i];d.time=z['physics_clock_s'][i]
  mujoco.mj_forward(m,d) # Kinematic/force evaluation for saved pose only; never mj_step.
  pelvis=display@d.xpos[m.body('pelvis').id];views=[]
  for azimuth in(0.,135.):
   cam=mujoco.MjvCamera();cam.lookat[:]=pelvis;cam.lookat[2]=max(.7,pelvis[2]-.15);cam.distance=2.8;cam.azimuth=azimuth;cam.elevation=-5.
   renderer.update_scene(d,camera=cam)
   if hide_ground:
    for g in renderer.scene.geoms[:renderer.scene.ngeom]:
     if int(g.type)==int(mujoco.mjtGeom.mjGEOM_PLANE):g.rgba[3]=0.
   if meta['native_frame']!='legacy_z_up':
    for g in renderer.scene.geoms[:renderer.scene.ngeom]:g.pos[:]=display@g.pos;g.mat[:]=(display@g.mat.reshape(3,3)).ravel()
   views.append(Image.fromarray(renderer.render().copy()))
  img=Image.new('RGB',(880,600),'white');img.paste(views[0],(0,80));img.paste(views[1],(440,80));draw=ImageDraw.Draw(img)
  draw.text((10,5),label[:100]+' | recorded replay',font=font,fill='black');draw.text((10,30),f't={run["t"][i]:.3f}s | 1x playback | x anterior, y superior, z right',font=font,fill='black')
  draw.text((10,54),'Air: ground visual hidden; saved states, no physics integration'if hide_ground else'Saved joint states; no added support, pose correction, or physics integration',font=font,fill='black');frames.append(img)
  if presentation is not None:
   clean=Image.new('RGB',(880,600),'white');clean.paste(views[0],(0,80));clean.paste(views[1],(440,80));banner=ImageDraw.Draw(clean)
   banner.text((10,5),'Guided Gait2392 | FSA Hunt/PT | supported ground steps',font=clean_title,fill='black')
   banner.text((10,31),f'1x playback | Vertical assistance {assistance:.1f}% BW | Horizontal restraint {horizontal:.1f} N RMS',font=clean_text,fill='black')
   banner.text((10,55),f'Recorded t={run["t"][i]:.3f}s | x anterior, y superior, z right',font=clean_text,fill='black');presentation_frames.append(clean)
 if hasattr(renderer,'close'):renderer.close()
 else:del renderer
 output=Path(output);output.mkdir(parents=True,exist_ok=True)
 if not still_only:frames[0].save(output/'RECORDED_REPLAY.gif',save_all=True,append_images=frames[1:],duration=round(1000/fps),loop=0)
 frames[len(frames)//2].save(output/'RECORDED_STILL.png')
 if presentation_frames:
  if not still_only:presentation_frames[0].save(output/'PRESENTATION_REPLAY.gif',save_all=True,append_images=presentation_frames[1:],duration=round(1000/fps),loop=0)
  presentation_frames[len(presentation_frames)//2].save(output/'PRESENTATION_STILL.png')
 provenance={'recording':str(source),'model':str(modelpath),'model_sha256':meta['model_sha256'],'frame_state_times_s':run['t'][ids].tolist(),'fps':None if still_only else fps,'still_only':still_only,'physics_steps':0,'forward_evaluations':'Saved-state rendering only','display_rotation_native_to_z_up':display.tolist(),'ground_visual_hidden_for_air':hide_ground,'ground_display_scope':'Only rendered plane alpha changed; archived model, saved poses, and physical contacts unchanged.'}
 if presentation is not None:provenance['presentation_copy']={'files':['PRESENTATION_STILL.png']if still_only else['PRESENTATION_REPLAY.gif','PRESENTATION_STILL.png'],'same_rendered_views_and_frame_states_as_archival_copy':True,'banner_vertical_metric':'Net signed prescribed plus restraint superior force, mean over5–13s','banner_horizontal_metric':'RMS magnitude of measured anterior/right restraint force over5–13s','full_support_decomposition':presentation}
 (output/'REPLAY_PROVENANCE.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')

def main():
 p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True);p.add_argument('--label',required=True);p.add_argument('--start',type=float);p.add_argument('--end',type=float);p.add_argument('--cycles',type=int,default=3);p.add_argument('--render',action='store_true');p.add_argument('--still-only',action='store_true');p.add_argument('--presentation',action='store_true');p.add_argument('--fps',type=int,default=20);a=p.parse_args();r=read_run(a.source);start,end,period=choose_window(r,a.start,a.end,a.cycles);report={'missing_display_channels':[]}if a.still_only else plot(r,a.output,a.label,start,end,period)
 presentation=json.loads((Path(a.output)/'SUPPORT_DECOMPOSITION_5_TO_13.json').read_text())if a.presentation else None
 if a.render or a.still_only:render(r,a.output,a.label,start,end,a.fps,still_only=a.still_only,presentation=presentation)
 print(json.dumps({'output':str(Path(a.output).resolve()),'window_s':[start,end],'missing_channels':report['missing_display_channels'],'replay':a.render},indent=2))
if __name__=='__main__':main()
