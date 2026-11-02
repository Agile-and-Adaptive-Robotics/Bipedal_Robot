"""Plot saved records only. No controller imports, model construction, or integration."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent
OUT=HERE/'assets/neural/slide_panels'
CAMPAIGN=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h')
REPO=Path(r'D:\GitHub\Bipedal_Robot')
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7']
plt.rcParams.update({'font.family':'Arial','font.size':14,'axes.labelsize':14,
 'xtick.labelsize':14,'ytick.labelsize':14,'legend.fontsize':13,
 'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
records=[]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(fig,name,sources,scope,extra=None):
    for ext in ('svg','png','pdf'):
        fig.savefig(OUT/(name+'.'+ext),dpi=220)
    records.append({'name':name,'sources':[{'path':str(p),'sha256':digest(p)}for p in sources],
        'scope':scope,'width_height_ratio':fig.get_figwidth()/fig.get_figheight(),
        'physics_steps':0,'neural_steps':0,'amplitude_normalization':False,
        'interpolation':False,'extra':extra or {},
        'outputs':{ext:digest(OUT/(name+'.'+ext))for ext in ('svg','png','pdf')}})
    plt.close(fig)

def ylim(ax,values,pad=.07):
    lo,hi=float(np.min(values)),float(np.max(values));margin=pad*max(hi-lo,.01)
    ax.set_ylim(lo-margin,hi+margin)

aligned=CAMPAIGN/'figures/ground9_rc_s4/aligned_trace.npz'
data=np.load(aligned,allow_pickle=False);names=data['neural_names'].tolist()
t=data['t_state'];mask=(t>=7.663)&(t<=13.)
groups=[(['RG_E_r','RG_F_r'],['Extensor','Flexor'],[0,1],['-','-'],'Rhythm generation'),
 (['PF_BI_JOINT_RECT_r','PF_BI_JOINT_HAM_r','PF_IN_BI_PAIR_RECT_r','PF_IN_BI_PAIR_HAM_r'],
  ['RF PF','HAM PF','RF relay','HAM relay'],[0,1,2,3],['-','-','--','--'],'RF/HAM pattern formation'),
 (['MN_bifemlh_r'],['Long-head MN'],[1],['-'],'Biceps femoris long head MN')]
fig,axes=plt.subplots(1,3,figsize=(12.5,4.0))
fig.subplots_adjust(left=.065,right=.99,bottom=.22,top=.64,wspace=.34)
for i,(cells,labels,cs,ls,title) in enumerate(groups):
    ax=axes[i]
    for c,label,color,style in zip(cells,labels,cs,ls):
        ax.plot(t[mask],data['neural_voltage'][mask,names.index(c)],color=COLORS[color],ls=style,lw=1.5,label=label)
    # Preserve the original row scale, including both sides and both source hip MNs.
    scale_cells=cells if i<2 else ['MN_iliacus_r','MN_bifemlh_r']
    both=scale_cells+[n.removesuffix('_r')+'_l'for n in scale_cells]
    ylim(ax,data['neural_voltage'][mask][:,[names.index(n)for n in both]])
    ax.set(xlim=(7.663,13.),xticks=[8,10,12],xlabel='Recorded time (s)',ylabel='Model mV')
    ax.grid(alpha=.17)
    pos=ax.get_position();center=(pos.x0+pos.x1)/2
    fig.text(center,.96,title if i<2 else 'Biceps femoris long head\nMN',ha='center',va='top',fontsize=15.5,weight='bold')
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.10),ncol=2 if i==1 else 1,frameon=False,
       handlelength=1.8,handletextpad=.45,columnspacing=.9,borderaxespad=0.)
save(fig,'NEURAL_RIGHT_HAM_HIERARCHY',[aligned],
 'Actual right RG and RF/HAM PF plus their saved downstream biceps-femoris-long-head MN. All share the same post-update aligned clocks. Third panel retains the original bilateral hip-MN row y scale.',
 {'window_s':[7.663,13.],'cells':[c for group in groups for c in group[0]],'clocks':'t_state; original aligned post-update state clock'})

# Contact comparison: the exact data arrays archived by its original figure generator.
contact=CAMPAIGN/'figures/stage5_contact_strength_review/FIGURE_DATA.npz'
z=np.load(contact,allow_pickle=False)
fig=plt.figure(figsize=(10,5.625));axes=[]
for col,x in enumerate((.10,.57)):
    axes.append([fig.add_axes([x,.515,.38,.22]),fig.add_axes([x,.13,.38,.20])])
for col in range(2):
    for ri,row in enumerate((1,6)):
        ax=axes[col][ri];tt=z[f'case{col}_row{row}_t'];yy=z[f'case{col}_row{row}_values']
        for j in range(2):ax.plot(tt,yy[:,j],color=COLORS[j],ls=('-','--')[j],lw=1.5)
        shared=np.concatenate([z[f'case{k}_row{row}_values'].ravel()for k in range(2)])
        ylim(ax,shared,.08)
        ax.set_xlim(7.666,13);ax.set_xticks([8,10,12]);ax.grid(alpha=.18)
        if ri==0:ax.tick_params(labelbottom=False)
        else:ax.set_xlabel('Recorded time (s)',fontsize=16)
        if col==0:ax.set_ylabel('Model mV' if ri==0 else 'Force (N)',fontsize=16)
        else:ax.tick_params(labelleft=False)
fig.text(.5,.98,'Ankle pattern-formation activity',ha='center',va='top',fontsize=19,weight='bold')
fig.text(.29,.905,'Quarter strength',ha='center',va='top',fontsize=18)
fig.text(.76,.905,'Original strength',ha='center',va='top',fontsize=18)
handles=[Line2D([],[],color=COLORS[0],lw=1.5,label='Plantarflexor'),Line2D([],[],color=COLORS[1],lw=1.5,ls='--',label='Dorsiflexor')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.858),ncol=2,frameon=False,fontsize=16,handlelength=2.2)
fig.text(.5,.462,'Vertical foot loading',ha='center',va='top',fontsize=19,weight='bold')
handles=[Line2D([],[],color=COLORS[0],lw=1.5,label='Right foot'),Line2D([],[],color=COLORS[1],lw=1.5,ls='--',label='Left foot')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.413),ncol=2,frameon=False,fontsize=16,handlelength=2.2)
save(fig,'CONTACT_PF_AND_FOOT_LOADING_WIDE',[contact],
 'Original figure arrays, row 1 right ankle PF and row 6 bilateral vertical foot force. All legends outside plot axes; original amplitudes, row scales, colors, dashed styles, and native timestamps preserved.',
 {'window_s':[7.666,13.],'clocks':'PF: original 50 ms pre-update t; force: original force interval ends','conditions':['quarter contact-to-PF strength','original strength']})

meta_path=CAMPAIGN/'figures/ground9_rc_s4/ALIGNED_TRACE_METADATA.json'
meta=json.loads(meta_path.read_text());bw=meta['body_weight_N'];tf=data['force_interval_end'];m=(tf>=7.663)&(tf<=13.)
fig,ax=plt.subplots(figsize=(8,4.5));fig.subplots_adjust(left=.135,right=.98,bottom=.205,top=.68)
ys=[data['foot_force_anatomical_N'][:,0,1]/bw,data['foot_force_anatomical_N'][:,1,1]/bw,
    data['harness_wrench_anatomical'][:,1]/bw,data['explicit_support_anatomical_N'][:,1]/bw]
labels=['Right foot','Left foot','Vertical restraint','Prescribed lift']
for y,label,color,style in zip(ys,labels,[COLORS[0],COLORS[1],'.4','.2'],['-','-','--',':']):
    ax.plot(tf[m],y[m],label=label,color=color,ls=style,lw=1.6)
ax.set(xlim=(7.663,13.),xticks=[8,10,12],xlabel='Interval end (s)',ylabel='Force (BW)');ax.grid(alpha=.2)
fig.text(.55,.98,'Foot loads and assistance',ha='center',va='top',fontsize=19,weight='bold')
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.90),ncol=2,frameon=False,fontsize=15,handlelength=2.1,columnspacing=1.4)
save(fig,'FOOT_LOADS_AND_ASSISTANCE_WIDE',[aligned,meta_path],
 'Same source arrays, original BW denominator, line styles, colors, and force-interval clocks as ALIGNED_MECHANICS. Legend moved outside axes; no change to amplitudes.',
 {'window_s':[7.663,13.],'body_weight_N':bw,'clocks':'force_interval_end'})

# The complete RG stimulus data, replacing the old image whose bottom tick labels were cropped.
stimdir=REPO/'Code/MuJoCo_SNS/spinal/reports_20261028_publication_flow/stimulus_long_controls'
record=next(r for r in json.loads((stimdir/'stimulus_results.json').read_text())if r['condition']=='RG_plus50_2cycles')
basefile=stimdir/'baseline.npz';stimfile=stimdir/'RG_plus50_2cycles.npz'
base=np.load(basefile,allow_pickle=False);stim=np.load(stimfile,allow_pickle=False)
tt=stim['t']-record['start_s'];m=(tt>=-2)&(tt<=9)
assert np.array_equal(stim['t'],base['t'])
assert np.array_equal(stim['V'][tt<0],base['V'][tt<0])
fig,ax=plt.subplots(figsize=(8,4.5));fig.subplots_adjust(left=.16,right=.98,bottom=.23,top=.77)
ax.plot(tt[m],base['V'][m,0],'--',color='.35',lw=1.6,label='Baseline')
ax.plot(tt[m],stim['V'][m,0],color=COLORS[0],lw=1.6,label='Stimulus')
ax.axvspan(0,record['duration_s'],color='#EECC66',alpha=.35)
ax.set(xlabel='Time from stimulus onset (s)',ylabel='RG extensor (mV)',xticks=[-2,0,2,4,6,8]);ax.grid(alpha=.15)
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.98),ncol=2,frameon=False,fontsize=17,handlelength=2.4)
save(fig,'DENG_RG_STIMULUS_COMPLETE',[basefile,stimfile,stimdir/'stimulus_results.json'],
 'Actual reduced Deng14 saved baseline and +50 nA RG stimulus, unchanged voltages and source-relative time. Complete time ticks and external Baseline/Stimulus key. Shading is the actual stimulus interval.',
 {'stimulus_start_absolute_s':record['start_s'],'duration_s':record['duration_s'],'source_relative_window_s':[-2,9]})

(OUT/'RENDER_FIX_PROVENANCE.json').write_text(json.dumps(records,indent=2),encoding='utf8')
print(json.dumps([{'name':r['name'],'aspect':r['width_height_ratio']}for r in records],indent=2))
