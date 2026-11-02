"""Projection-sized versions of existing recorded figures. No simulation."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent;OUT=HERE/'assets/neural/slide_panels'
CAMPAIGN=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h')
REPO=Path(r'D:\GitHub\Bipedal_Robot')
BLUE,ORANGE='#0072B2','#D55E00'
plt.rcParams.update({'font.family':'Arial','font.size':18,'axes.labelsize':18,
 'xtick.labelsize':18,'ytick.labelsize':18,'legend.fontsize':18,
 'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
records=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(fig,name,sources,scope,extra):
    for ext in ('svg','png','pdf'):fig.savefig(OUT/(name+'.'+ext),dpi=220)
    records.append({'name':name,'sources':[{'path':str(p),'sha256':sha(p)}for p in sources],
        'scope':scope,'figure_inches':list(fig.get_size_inches()),'aspect':fig.get_figwidth()/fig.get_figheight(),
        'physics_steps':0,'neural_steps':0,'amplitude_normalization':False,'interpolation':False,
        'details':extra,'outputs':{ext:sha(OUT/(name+'.'+ext))for ext in ('svg','png','pdf')}})
    plt.close(fig)

stimdir=REPO/'Code/MuJoCo_SNS/spinal/reports_20261028_publication_flow/stimulus_long_controls'
summary=stimdir/'stimulus_results.json'
record=next(r for r in json.loads(summary.read_text())if r['condition']=='RG_plus50_2cycles')
basefile=stimdir/'baseline.npz';stimfile=stimdir/'RG_plus50_2cycles.npz'
base=np.load(basefile,allow_pickle=False);stim=np.load(stimfile,allow_pickle=False)
tt=stim['t']-record['start_s'];mask=(tt>=-2)&(tt<=9)
assert np.array_equal(stim['t'],base['t'])
assert np.array_equal(stim['V'][tt<0],base['V'][tt<0])
fig,ax=plt.subplots(figsize=(6.4,3.8));fig.subplots_adjust(left=.18,right=.985,bottom=.24,top=.76)
ax.plot(tt[mask],base['V'][mask,0],'--',color='.35',lw=1.6,label='Baseline')
ax.plot(tt[mask],stim['V'][mask,0],color=BLUE,lw=1.6,label='Stimulus')
ax.axvspan(0,record['duration_s'],color='#EECC66',alpha=.35)
ax.set(xlabel='Time from stimulus onset (s)',ylabel='RG extensor (mV)',xticks=[-2,0,2,4,6,8]);ax.grid(alpha=.15)
ax.xaxis.label.set_size(16);ax.yaxis.label.set_size(16);ax.tick_params(labelsize=16)
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.56,.995),ncol=2,frameon=False,fontsize=16,handlelength=2.2,columnspacing=1.)
save(fig,'DENG_RG_STIMULUS_PROJECTION',[basefile,stimfile,summary],
 'Same actual saved reduced Deng14 RG extensor data, baseline, and +50 nA stimulus as the complete derivative. Full external key and time ticks.',
 {'window_s':[-2,9],'stimulus_duration_s':record['duration_s'],'font_pt':16,'target_frame':[525,268]})

contact=CAMPAIGN/'figures/stage5_contact_strength_review/FIGURE_DATA.npz';z=np.load(contact,allow_pickle=False)
fig=plt.figure(figsize=(10,4));axes=[]
for x in (.235,.655):axes.append([fig.add_axes([x,.52,.325,.235]),fig.add_axes([x,.17,.325,.185])])
for col in range(2):
    for ri,row in enumerate((1,6)):
        ax=axes[col][ri];t=z[f'case{col}_row{row}_t'];y=z[f'case{col}_row{row}_values']
        for j in range(2):ax.plot(t,y[:,j],color=(BLUE,ORANGE)[j],ls=('-','--')[j],lw=1.5)
        vals=np.concatenate([z[f'case{k}_row{row}_values'].ravel()for k in range(2)])
        lo,hi=float(vals.min()),float(vals.max());pad=.08*max(hi-lo,.01);ax.set_ylim(lo-pad,hi+pad)
        ax.set_xlim(7.666,13.);ax.set_xticks([8,10,12]);ax.grid(alpha=.18);ax.tick_params(labelsize=18)
        if ri==0:ax.set_yticks([-2.5,0,2.5]);ax.tick_params(labelbottom=False)
        else:ax.set_yticks([0,1000])
        if col==1:ax.tick_params(labelleft=False)
fig.text(.3975,.99,'Quarter strength',ha='center',va='top',fontsize=19,weight='bold')
fig.text(.8175,.99,'Original strength',ha='center',va='top',fontsize=19,weight='bold')
fig.text(.015,.64,'Ankle pattern-\nformation activity\n(model mV)',ha='left',va='center',fontsize=16,linespacing=1.15)
fig.text(.015,.255,'Vertical foot\nloading (N)',ha='left',va='center',fontsize=16,linespacing=1.15)
handles=[Line2D([],[],color=BLUE,lw=1.5,label='Plantarflexor'),Line2D([],[],color=ORANGE,lw=1.5,ls='--',label='Dorsiflexor')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.61,.90),ncol=2,frameon=False,fontsize=18,handlelength=1.9,columnspacing=1.2)
handles=[Line2D([],[],color=BLUE,lw=1.5,label='Right foot'),Line2D([],[],color=ORANGE,lw=1.5,ls='--',label='Left foot')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.61,.465),ncol=2,frameon=False,fontsize=18,handlelength=1.9,columnspacing=1.2)
fig.text(.61,.008,'Recorded time (s)',ha='center',va='bottom',fontsize=18)
save(fig,'CONTACT_PF_AND_FOOT_LOADING_PROJECTION',[contact],
 'Same archived contact-comparison arrays, original row scales, amplitudes, colors, line styles, and PF versus force clocks. All legends external.',
 {'window_s':[7.666,13.],'axes_and_legend_font_pt':18,'target_frame':[870,333],'rows':['Ankle pattern-formation activity','Vertical foot loading']})

aligned=CAMPAIGN/'figures/ground9_rc_s4/aligned_trace.npz';meta_path=aligned.parent/'ALIGNED_TRACE_METADATA.json'
data=np.load(aligned,allow_pickle=False);meta=json.loads(meta_path.read_text());bw=meta['body_weight_N'];t=data['force_interval_end'];mask=(t>=7.663)&(t<=13.)
fig,ax=plt.subplots(figsize=(5.2,4.2));fig.subplots_adjust(left=.205,right=.985,bottom=.235,top=.61)
ys=[data['foot_force_anatomical_N'][:,0,1]/bw,data['foot_force_anatomical_N'][:,1,1]/bw,
 data['harness_wrench_anatomical'][:,1]/bw,data['explicit_support_anatomical_N'][:,1]/bw]
for y,label,color,style in zip(ys,['Right foot','Left foot','Restraint','Lift'],[BLUE,ORANGE,'.4','.2'],['-','-','--',':']):
    ax.plot(t[mask],y[mask],color=color,ls=style,label=label,lw=1.6)
ax.set(xlim=(7.663,13.),xticks=[8,10,12],xlabel='Time (s)',ylabel='Force (BW)');ax.grid(alpha=.2);ax.tick_params(labelsize=18)
fig.text(.57,.985,'Vertical forces',ha='center',va='top',fontsize=19,weight='bold')
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.57,.88),ncol=2,frameon=False,fontsize=17,handlelength=1.65,columnspacing=.8,handletextpad=.5)
save(fig,'FOOT_LOADS_AND_ASSISTANCE_PROJECTION',[aligned,meta_path],
 'Same aligned source forces and body-weight denominator. Restraint means signed superior harness force; lift is prescribed superior assistance. Times are original force-interval ends.',
 {'window_s':[7.663,13.],'body_weight_N':bw,'axes_font_pt':18,'legend_font_pt':17,'target_frame':[380,295]})

path=REPO/'Code/MuJoCo_SNS/spinal/reports_20261027_codex/perturbation_results.json'
items=json.loads(path.read_text());deletions=[x for x in items if x['job']['condition'].startswith('del_')]
labels=[x['job']['condition'].replace('del_','').replace('_',' ')for x in deletions]
left=[x['metrics']['cycL']for x in deletions];right=[x['metrics']['cycR']for x in deletions]
assert len(labels)==7
fig,ax=plt.subplots(figsize=(10.5,4.5));fig.subplots_adjust(left=.23,right=.985,bottom=.22,top=.80)
y=np.arange(len(labels))
ax.barh(y-.16,left,height=.3,color='#4477AA',label='Left')
ax.barh(y+.16,right,height=.3,color='white',edgecolor='#EE6677',hatch='///',label='Right')
ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(0,37)
ax.set_xlabel('Merged heel-contact events, 4.5–30 s',fontsize=18);ax.tick_params(labelsize=18)
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.61,.985),ncol=2,frameon=False,fontsize=18,handlelength=2.)
save(fig,'W2L_PATHWAY_REMOVALS_PROJECTION',[path],
 'Exact saved frozen W2L pathway-removal counts from original plot_results.py. Same order, blue Left bars, striped red Right bars, and event definition. Legend outside axes. Counts are contact events, not physiological strides.',
 {'conditions':labels,'left_counts':left,'right_counts':right,'axes_and_legend_font_pt':18,'original_generator':str(path.parent/'plot_results.py')})

(OUT/'PROJECTION_RENDER_PROVENANCE.json').write_text(json.dumps(records,indent=2),encoding='utf8')
print(json.dumps([{'name':r['name'],'aspect':r['aspect']}for r in records],indent=2))
