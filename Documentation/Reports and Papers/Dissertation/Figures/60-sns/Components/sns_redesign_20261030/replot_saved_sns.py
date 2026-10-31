from pathlib import Path
import numpy as np,scipy.io as sio, json, hashlib, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'Testing_Data').exists()) if any((p/'Testing_Data').exists()for p in Path(__file__).resolve().parents) else Path(r'D:/Github/Bipedal_Robot')
D=ROOT/'Documentation/Reports and Papers/Dissertation'; S=ROOT/'Code/MuJoCo_SNS/spinal'
W=Path(__file__).resolve().parent;O=W/'output';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.labelsize':12,'xtick.labelsize':11.5,'ytick.labelsize':11.5,'legend.fontsize':11.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12.5})
EXT='#D55E00';FLX='#0072B2';GRAY='#777777';TEAL='#009E73'
manifest=[]
def save(fig,name,alt):
 floor=16.2 if name=='deng_neuromechanical_stimulus' else 14.5 if name=='spiking_ground_neuromechanics' else 13.5 if name=='w2l_feedback_verified' else 11
 for ax in fig.axes:
  ax.tick_params(labelsize=floor)
  ax.xaxis.label.set_fontsize(max(floor,ax.xaxis.label.get_fontsize()));ax.yaxis.label.set_fontsize(max(floor,ax.yaxis.label.get_fontsize()))
  ax.title.set_fontsize(max(floor,ax.title.get_fontsize()))
 for obj in [getattr(fig,'_suptitle',None),getattr(fig,'_supxlabel',None)]:
  if obj is not None:obj.set_fontsize(max(floor,obj.get_fontsize()))
 fig.canvas.draw()
 fig.savefig(O/(name+'.pdf'));fig.savefig(O/(name+'.png'),dpi=140);plt.close(fig)
 (O/(name+'.alt.txt')).write_text(alt,encoding='utf8');manifest.append({'name':name,'alt':alt})
def style(ax):ax.grid(alpha=.18);ax.tick_params(length=3)
# 6.3 three complete driven cycles, two paired half-centers per leg.
z=np.load(S/'campaigns/20260930/w2l_numpy_ref.npz',allow_pickle=True);m=sio.loadmat(ROOT/'Code/Matlab/SNS_Simscape/results/SNS_W2L_CPG_run_20260930.mat',squeeze_me=True)
t=z['t'];ts=m['t'];wn=(t>=2)&(t<=5);ws=(ts>=2)&(ts<=5)
fig,axs=plt.subplots(2,1,figsize=(6.5,4.5),sharex=True,layout='constrained')
for a,side,prefix in zip(axs,['Left','Right'],['L','R']):
 for suffix,c,mt in [('ext',EXT,'E'),('flx',FLX,'F')]:
  n=list(z['names']).index(f'{prefix} RG {suffix}')
  a.plot(t[wn],z['V'][wn,n],color=c,lw=1.8,label=('Extensor' if suffix=='ext' else'Flexor')+' (NumPy)')
  a.plot(ts[ws],m[f'v{prefix}{mt}'][ws],color=c,lw=1.1,ls='--',label=('Extensor' if suffix=='ext' else'Flexor')+' (Simulink)')
 a.set_ylabel('RG voltage\n(mV)');a.set_title(side+' rhythm generator',loc='left',fontweight='bold');a.set_xlim(2,5);style(a)
fig.suptitle('Extensor: orange; flexor: blue | NumPy solid, Simulink dashed',fontsize=10.8);axs[-1].set_xlabel('Time (s)')
save(fig,'sns_w2l_traces','Three recorded cycles, 2-5 s, of paired extensor/flexor RG voltages in the left and right rhythm generators. NumPy solid and Simulink dashed traces use the same original time axis; no alignment shift. Colors distinguish antagonists.')
#6.4 same saved physical trials; only actual logged RG-extensor channels are available.
a=np.load(S/'reports_20261027_codex/panel_control.npz');b=np.load(S/'reports_20261027_codex/panel_del_all_feedback.npz')
fig,axs=plt.subplots(5,2,figsize=(6.5,8.0),sharex=True,layout='constrained')
rows=[('rgE',0,1,'RG extensor\nreadout (mV)'),('q',0,3,'Hip angle\n(deg)'),('q',1,4,'Knee angle\n(deg)'),('q',2,5,'Ankle angle\n(deg)'),('heelF',0,1,'Heel normal\nforce (N)')]
for j,(zz,title,phase) in enumerate([(a,'Feedback intact','0.498 cycles'),(b,'Feedback disabled','0.223 cycles')]):
 tt=zz['t'];w=(tt>=12)&(tt<=16)
 for i,(key,il,ir,lab) in enumerate(rows):
  ax=axs[i,j];ax.plot(tt[w],zz[key][w,il],color=FLX,lw=1.3,label='Left');ax.plot(tt[w],zz[key][w,ir],color=EXT,lw=1.3,ls='--',label='Right');style(ax);ax.set_xlim(12,16)
  if j==0:ax.set_ylabel(lab)
 axs[0,j].set_title(title+'\nContact phase '+phase,fontsize=11.5,fontweight='bold')
 axs[-1,j].set_xlabel('Time (s)')
for i in range(5):
 lo=min(axs[i,0].get_ylim()[0],axs[i,1].get_ylim()[0]);hi=max(axs[i,0].get_ylim()[1],axs[i,1].get_ylim()[1]);axs[i,0].set_ylim(lo,hi);axs[i,1].set_ylim(lo,hi)
fig.suptitle('Left: blue solid | Right: orange dashed',fontsize=11)
fig.supxlabel('Full-trial contact-phase acceptance: 0.35-0.65 cycles',fontsize=11)
save(fig,'w2l_feedback_verified','Matched 12-16 s portions of the recorded harness-supported W2L feedback-intact and sensory-feedback-disabled runs. Rows align bilateral RG extensor readout, hip, knee, and ankle angles, and heel normal forces. Left blue solid, right orange dashed; row limits match. Full-trial phase statistics remain 0.498 and 0.223 cycles; these statistics were not recomputed from the display window.')
#D6 replace neural-only illustration by distinct clearly labeled saved ground comparison.
a=np.load(S/'spinal_run_spkbase_s3k.npz');b=np.load(S/'spinal_run_spkbase_spiking.npz')
fig,axs=plt.subplots(6,2,figsize=(6.5,8.7),sharex=True,layout='constrained')
for j,(zz,title) in enumerate([(a,'Graded controller'),(b,'Untuned spiking controller')]):
 tt=zz['t'];ww=(tt>=10)&(tt<=14);jn=list(zz['key_joints']);nn=list(zz['neuro_names']);an=list(zz['key_acts'])
 series=[(zz['neuro'][:,nn.index('RG_E_r')],zz['neuro'][:,nn.index('RG_F_r')],'RG readout\n(mV)'),(zz['neuro'][:,nn.index('PF_E1_r')],zz['neuro'][:,nn.index('PF_F1_r')],'PF readout\n(mV)'),(zz['act'][:,an.index('vas_lat_r')],zz['act'][:,an.index('semimem_r')],'Muscle\nactivation')]
 for part in ['hip_flexion','knee_angle','ankle_angle']:
  series.append((zz['q'][:,jn.index(part+'_r')],zz['q'][:,jn.index(part+'_l')],part.split('_')[0].capitalize()+' angle\n(deg)'))
 for i,(x,y,label) in enumerate(series):
  ax=axs[i,j];ax.plot(tt[ww],x[ww],color=EXT,lw=1.2,label='Extensor' if i<3 else 'Right');ax.plot(tt[ww],y[ww],color=FLX,lw=1.2,ls='--',label='Flexor' if i<3 else 'Left');style(ax);ax.set_xlim(10,14)
  if j==0:ax.set_ylabel(label)
 axs[0,j].set_title(title,fontweight='bold');axs[-1,j].set_xlabel('Time (s)')
fig.suptitle('Neural: orange extensor, blue flexor\nJoints: orange right, blue left',fontsize=11)
for i in range(6):
 lo=min(x.get_ylim()[0]for x in axs[i]);hi=max(x.get_ylim()[1]for x in axs[i]);[x.set_ylim(lo,hi)for x in axs[i]]
for text_obj in fig.findobj(matplotlib.text.Text):
 if text_obj.get_text():text_obj.set_fontsize(max(13,text_obj.get_fontsize()))
save(fig,'spiking_ground_neuromechanics','Saved 10-14 s ground-evaluation windows at identical selected parameters. Columns show graded and untuned spiking variants; rows show right RG/PF extensor-flexor readouts, right vastus-lateralis and semimembranosus muscle activation, and bilateral hip/knee/ankle angles. Neural readouts are the archived controller outputs, not individual spiking membrane events. Muscle activation is not labeled MN voltage. This is a supported coupled ground trial, not the constant-drive network-only period test.')
#Deng-derived recorded stimulus: exactly recorded V and joint q, no interpolation or new integration.
p=D/'Archive/flow_and_validation_20261028/stimulus_long_controls';base=np.load(p/'baseline.npz');zz=np.load(p/'PF_plus50_2cycles.npz');rows=json.loads((p/'stimulus_results.json').read_text());rr=next(x for x in rows if x['condition']=='PF_plus50_2cycles');tt=zz['t']-rr['start_s'];ww=(tt>=-rr['baseline_period_s'])&(tt<=3*rr['baseline_period_s']);dur=rr['duration_s']
fig,axs=plt.subplots(7,1,figsize=(6.5,9.0),sharex=True,layout='constrained')
for i,(ie,iff,label) in enumerate([(0,1,'RG voltage\n(mV)'),(4,5,'Hip PF\nvoltage\n(mV)'),(12,13,'Hip MN\nvoltage\n(mV)')]):
 axs[i].plot(tt[ww],zz['V'][ww,ie],color=EXT,lw=1.25,label='Extensor');axs[i].plot(tt[ww],zz['V'][ww,iff],color=FLX,ls='--',lw=1.25,label='Flexor');axs[i].set_ylabel(label)
for i,lab in enumerate(['Hip\nangle\n(deg)','Knee\nangle\n(deg)','Ankle\nangle\n(deg)'],3):
 j=i-3;axs[i].plot(tt[ww],np.rad2deg(base['q'][ww,j]),color=GRAY,ls='--',lw=1.35,label='No stimulus');axs[i].plot(tt[ww],np.rad2deg(zz['q'][ww,j]),color=TEAL,lw=1.35,label='Stimulus');axs[i].set_ylabel(lab)
axs[6].plot(tt[ww],np.where((tt[ww]>=0)&(tt[ww]<dur),rr['current_nA'],0),color='black',lw=1.25);axs[6].set_ylabel('PF current\n(nA)');axs[6].set_xlabel('Time from stimulus onset (s)')
for ax in axs:ax.axvspan(0,dur,color='#EECC66',alpha=.22);style(ax);ax.set_xlim(tt[ww][0],tt[ww][-1])
fig.suptitle('Neural: orange extensor, blue flexor\nAngles: green stimulus, gray dashed baseline',fontsize=11)
save(fig,'deng_neuromechanical_stimulus','Recorded reduced Deng-derived circuit with suspended limb, PF extensor stimulus +50 nA for two baseline periods (3.876 s). Common time axis spans one baseline period before to one period after stimulus. Antagonist RG, hip PF, hip MN voltages above hip/knee/ankle angles and applied current. Kinematics compare actual matched stimulus and no-stimulus saved runs; no new simulation. Shading is exact stimulus interval; times relative to onset at 10.1765 s. This is a translated reduced circuit, not an original rat-paper trace.')
(O/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8');print('Created',len(manifest),'publication figures in',O)
