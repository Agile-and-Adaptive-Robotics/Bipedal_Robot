import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,
                     'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,
                     'pdf.fonttype':42,'ps.fonttype':42})
BLUE,RED,TEAL,PURPLE='#4477AA','#EE6677','#228833','#AA3377'


def traces(ax,name,kind='q'):
    with np.load(HERE/(name+'.npz')) as d:
        t=d['t'];w=t>=4.5
        if kind=='q':
            ax.plot(t[w],d['q'][w,1],color=BLUE,label='Left',lw=1.1)
            ax.plot(t[w],d['q'][w,4],color=RED,label='Right',ls='--',lw=1.1)
            ax.set_ylabel('Knee angle (deg)')
        elif kind=='force':
            ax.plot(t[w],d['heelF'][w,0],color=BLUE,label='Left',lw=.8)
            ax.plot(t[w],d['heelF'][w,1],color=RED,label='Right',ls='--',lw=.8)
            ax.set_ylabel('Heel normal force (N)')
        elif kind=='com':
            ax.plot(t[w],d['com'][w,0]-d['com'][w,0][0],color=TEAL,lw=1.3)
            ax.set_ylabel('Forward displacement (m)')
        ax.set_xlabel('Time (s)')
        ax.grid(alpha=.2)
        ax.spines[['top','right']].set_visible(False)


def save(fig,name,alt):
    for extension in ('png','pdf'):
        fig.savefig(HERE/(name+'.'+extension),dpi=180,bbox_inches='tight')
    (HERE/(name+'.alt.txt')).write_text(alt,encoding='utf-8')
    plt.close(fig)


fig,axes=plt.subplots(3,2,figsize=(12,9),layout='constrained',sharex=True)
for col,(name,title) in enumerate((('panel_control','Feedback intact'),
                                  ('panel_del_all_feedback','All sensory feedback disabled'))):
    for row,kind in enumerate(('q','force','com')):
        traces(axes[row,col],name,kind)
    axes[0,col].set_title(title)
    axes[0,col].legend(loc='upper right',ncol=2)
fig.suptitle('Frozen W2L conductances: sensory-feedback comparison\nHarness-supported model; feedforward support setting 0.25',fontsize=12)
save(fig,'feedback_comparison','Three rows compare left and right knee angles, heel normal forces, and forward displacement with sensory feedback intact or fully disabled. All conductances, drive, and support settings are frozen. The disabled arm fails the alternation gate. The support harness remains active in both arms. Knee flexion is negative.')

fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
for ax,name,title in zip(axes[0],('control','stim_mild','stim_strong'),
                         ('(A) Control','(B) RG drive reduction: 1.5 nA','(C) RG drive reduction: 3 nA')):
    traces(ax,'panel_'+name)
    ax.set_title(title)
    ax.set_xlim(5,16)
    if name!='control':
        ax.axvspan(8,9,color='#CCBB44',alpha=.3)
axes[0,0].legend(ncol=2)
for ax,name,title in zip(axes[1,:2],('push_fwd','push_back'),
                         ('(D) Forward push: 40 N, 0.15 s','(E) Backward push: 40 N, 0.15 s')):
    traces(ax,'panel_'+name,'com')
    ax.set_title(title)
    ax.axvspan(8,8.15,color='#CCBB44',alpha=.4)
items=json.loads((HERE/'perturbation_results.json').read_text())
deletions=[x for x in items if x['job']['condition'].startswith('del_')]
labels=[x['job']['condition'].replace('del_','').replace('_',' ') for x in deletions]
y=np.arange(len(labels))
axes[1,2].barh(y-.16,[x['metrics']['cycL'] for x in deletions],height=.3,color=BLUE,label='Left')
axes[1,2].barh(y+.16,[x['metrics']['cycR'] for x in deletions],height=.3,color='white',edgecolor=RED,hatch='///',label='Right')
axes[1,2].set_yticks(y,labels);axes[1,2].invert_yaxis()
axes[1,2].set_xlabel('Merged heel-contact events, 4.5–30 s')
axes[1,2].set_title('(F) Frozen pathway removals');axes[1,2].legend(ncol=2)
fig.suptitle('W2L stimulus and pathway-removal experiments\nTrial 458, 56 verified conductance dimensions; harness-supported',fontsize=12)
save(fig,'w2l_perturbation_panels','Six panels show control knee angles, two rhythm-generator drive-reduction pulses, forward and backward pushes, and independent left and right contact-event counts after pathway deletions. The stimulus reduces tonic RG extensor drive; it is not the Deng paper original excitation protocol. All sensory feedback removal and commissural removal fail the existing heel-contact stepping gate. The model retains the support harness. Contact events are not asserted to be physiological strides. Two legacy ankle-to-RG/PF keys have no matching edges and were inert; all 56 effective conductance dimensions were verified.')

validation=json.loads((HERE/'validation_results.json').read_text())
summary={kind:[x['metrics']['survive'] for x in validation if x['job']['kind']==kind]
         for kind in ('repeat','basin','frozen_support_ladder')}
print('Validation', {k:(sum(v),len(v)) for k,v in summary.items()})
print('Figures saved',HERE)
