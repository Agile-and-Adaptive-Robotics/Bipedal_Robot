from pathlib import Path
import json,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
R=Path(r'D:\GitHub\Bipedal_Robot')
W=R/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031'
RUN=W/'runs/W2L_VirtualWalker_Newer'
OUT=R/'Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/current_W2L_results'
OUT.mkdir(exist_ok=True,parents=True)
BLUE='#0077BB';ORANGE='#EE7733';GRAY='#aaaaaa'
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'axes.titlesize':10,'svg.fonttype':'none','pdf.fonttype':42})
base=np.load(RUN/'canonical108_air_baseline_v2/trajectory.npz',allow_pickle=True)
source=json.loads((RUN/'canonical108_air_baseline_v2/RESULT.json').read_text())['source_sha256']
provenance={'machine':'easteregg2','neuron_count':108,'chemical_edges':216,'source_sha256':source,'staging_only':True,'runs':[],'figures':[]}
def load(tag):
 p=RUN/tag;r=json.loads((p/'RESULT.json').read_text());z=np.load(p/'trajectory.npz',allow_pickle=True)
 assert r['neurons']==108 and r['nonspiking_edges']==216 and r['source_sha256']==source and len(z['neuron_names'])==108
 provenance['runs'].append({'folder':str(p),'trajectory_sha256':hashlib.sha256((p/'trajectory.npz').read_bytes()).hexdigest(),'result_sha256':hashlib.sha256((p/'RESULT.json').read_bytes()).hexdigest()})
 return p,r,z
def save(fig,stem,window):
 fs=[x.get_fontsize()for x in fig.findobj(matplotlib.text.Text)if x.get_text()];assert min(fs)>=10,(stem,min(fs))
 for ext in('png','pdf','svg'):fig.savefig(OUT/(stem+'.'+ext),dpi=300)
 provenance['figures'].append({'stem':stem,'width_in':fig.get_size_inches()[0],'height_in':fig.get_size_inches()[1],'minimum_text_pt':min(fs),'window_s':window})
 plt.close(fig)
fig,axes=plt.subplots(7,3,figsize=(6.5,7.1),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[.65,1,1,1,.85,.85,.85]})
labels=['Current\n(nA)','RG\n(mV)','Hip PF\n(mV)','Hip MN\n(mV)','Hip\n(deg)','Knee\n(deg)','Ankle\n(deg)']
for c,(layer,target)in enumerate((('RG','L RG ext'),('PF','L Hip PF ext'),('MN','L Hip MN ext'))):
 p,r,z=load('canonical108_air_'+layer+'_minus10nA_100ms_v1');t=z['t'];sel=(t>=3.4)&(t<=7.8);tt=t[sel];names=z['neuron_names'].tolist();bn=base['neuron_names'].tolist();pulse=r['neural_pulse'];on=pulse['start_s'];off=on+pulse['duration_s']
 for j,ax in enumerate(axes[:,c]):
  ax.spines[['top','right']].set_visible(False);ax.grid(color='#dddddd',lw=.4);ax.axvspan(on,off,color='#E69F00',alpha=.22,zorder=0);ax.set_xlim(3.4,7.8);ax.yaxis.set_major_locator(MaxNLocator(nbins=2));ax.set_ylabel(labels[j]if c==0 else'')
 axes[0,c].set_title('('+chr(65+c)+') '+layer+' stimulus',fontweight='bold')
 axes[0,c].step(tt,z['neural_stimulus_nA'][sel,names.index(target)],where='post',color='#333333',lw=1.2)
 for row,prefix in((1,'L RG '),(2,'L Hip PF '),(3,'L Hip MN ')):
  for tok,color,style in(('ext',BLUE,'-'),('flx',ORANGE,'--')):
   n=prefix+tok;axes[row,c].plot(tt,np.interp(tt,base['t'],base['neuro'][:,bn.index(n)]),style,color=GRAY,lw=1.15,zorder=1);axes[row,c].plot(tt,z['neuro'][sel,names.index(n)],style,color=color,lw=1.1,zorder=2)
 for q,row in enumerate((4,5,6)):
  axes[row,c].plot(tt,np.interp(tt,base['t'],base['angles'][:,q]),color=GRAY,lw=1.15)
  axes[row,c].plot(tt,z['angles'][sel,q],color='#222222',lw=1.1)
 axes[-1,c].set_xticks([4,5,6,7]);axes[-1,c].set_xlabel('Time (s)')
for row in range(7):
 low=min(ax.get_ylim()[0]for ax in axes[row]);high=max(ax.get_ylim()[1]for ax in axes[row])
 for ax in axes[row]:ax.set_ylim(low,high)
for ax in axes[0]:ax.set_ylim(-11,1);ax.set_yticks([-10,0])
fig.legend(handles=[Line2D([],[],color=BLUE,label='Extensor'),Line2D([],[],color=ORANGE,ls='--',label='Flexor'),Line2D([],[],color=GRAY,label='Unperturbed')],loc='outside upper center',ncol=3,frameon=False,handlelength=1.4,columnspacing=1.)
save(fig,'W2L108_left_RG_PF_MN_stimuli',[3.4,7.8])
p,r,z=load('canonical108_ground_rustuna_winner_fresh_v1');t=z['t'];sel=(t>=4)&(t<=8);tt=t[sel];bw=r['robot_mass_kg']*9.81
fig,axes=plt.subplots(6,1,figsize=(6.5,7.),sharex=True,layout='constrained')
for i,ax in enumerate(axes):ax.spines[['top','right']].set_visible(False);ax.grid(color='#dddddd',lw=.45);ax.yaxis.set_major_locator(MaxNLocator(nbins=3));ax.set_xlim(4,8)
for side,color,style,label in((0,BLUE,'-','Left'),(1,ORANGE,'--','Right')):
 for row,q in enumerate((0,1,2)):axes[row].plot(tt,z['angles'][sel,q+3*side],style,color=color,lw=1.3,label=label)
 axes[3].plot(tt,z['foot_vertical'][sel,side]/bw*100,style,color=color,lw=1.1)
 axes[4].plot(tt,np.min(z['foot_clearance'][sel,2*side:2*side+2],axis=1)*1000,style,color=color,lw=1.3)
net=sum(z[name][:,1]for name in('rail_force_world','guide_force_world','harness_force'))/bw*100
axes[5].plot(tt,net[sel],color='#222222',lw=1.1)
for ax,label in zip(axes,['Hip (deg)','Knee (deg)','Ankle (deg)','Foot load\n(% BW)','Whole-foot\nclearance (mm)','Net assistance\n(% BW)']):ax.set_ylabel(label)
axes[0].legend(ncol=2,frameon=False,loc='upper right');axes[-1].set_xlabel('Time (s)');axes[4].axhline(0,color='#777777',lw=.7)
save(fig,'W2L108_tuned_ground_motion_load',[4,8])
rows=[]
for tag,label in [('canonical108_ground_baseline_v1','Toe-down initial pose'),('canonical108_ground_nearflat_source_phase_v1','Near-flat IK pose'),('canonical108_ground_nearflat_ankle_crossed_v1','IK + ankle phase correction'),('canonical108_ground_rustuna_winner_fresh_v1','IK + phase correction + fitting')]:
 p,r,z=load(tag);a=json.loads((p/'ASSESSMENT.json').read_text());rows.append({'condition':label,'run':str(p),'swings_L':a['legs']['L']['unloaded_clear_swings'],'swings_R':a['legs']['R']['unloaded_clear_swings'],'alternation_percent':a['alternating_transition_fraction']*100,'mean_net_assistance_percent_BW':r['mean_net_vertical_assistance_BW']*100,'knee_min_L_deg':r['anatomical_knee_min_deg'][0],'knee_min_R_deg':r['anatomical_knee_min_deg'][1],'duration_s':r['duration_completed_s']})
(OUT/'NUMERICAL_RESULTS.json').write_text(json.dumps(rows,indent=2))
tab=['\\begin{table}[htbp]','\\centering','\\caption{Changes in supported stepping after initial-pose correction, ankle phase correction, and bounded controller fitting. Swing counts and alternation are evaluated from 1 to 10 s; assistance is the mean net vertical support force over the full 10 s, normalized by model body weight.}','\\label{tab:w2l108_ground_changes}','\\begin{tabular}{@{}lrrr@{}}','\\toprule','Condition & Swings L/R & Alternation (\\%) & Assistance (\\% BW) \\\\','\\midrule']
for row in rows:tab.append(f"{row['condition']} & {row['swings_L']}/{row['swings_R']} & {row['alternation_percent']:.2f} & {row['mean_net_assistance_percent_BW']:.2f} \\\\")
tab+=['\\bottomrule','\\end{tabular}','\\end{table}']
table_text='\n'.join(tab)+'\n'
table_text=table_text.replace('\\begin{tabular}{@{}lrrr@{}}','\\small\n\\begin{tabularx}{\\textwidth}{@{}Xrrr@{}}').replace('Condition & Swings L/R & Alternation (\\%) & Assistance (\\% BW)','Condition & \\shortstack{Swings\\\\L/R} & \\shortstack{Alternation\\\\(\\%)} & \\shortstack{Assistance\\\\(\\% BW)}').replace('\\end{tabular}','\\end{tabularx}')
(OUT/'RESULTS_TABLE.tex').write_text(table_text)
best=json.loads((W/'runs/rustuna108_ground_home_v1/BEST.json').read_text());trials=[json.loads(x)for x in(W/'runs/rustuna108_ground_home_v1/TRIALS.jsonl').read_text().splitlines()];assert len(trials)==24 and all(x['diagnostics']['physical_gate_passed']for x in trials)
provenance['optimizer']={'best_file':str(W/'runs/rustuna108_ground_home_v1/BEST.json'),'trials_file':str(W/'runs/rustuna108_ground_home_v1/TRIALS.jsonl'),'accepted':len(trials),'total':len(trials),'params':best['trial']['params'],'objective':'Post-1-s upward assistance + foot-force and relative-belt-speed penalties, subject to swing/unloading/alternation/knee/root-height acceptance. Not a metabolic objective.'}
(OUT/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2));print(OUT);print(json.dumps(rows,indent=2))
