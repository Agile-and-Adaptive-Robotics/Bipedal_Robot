"""Plot recorded neural-control evidence without classifying contact chatter as gait."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
def main(run):
    data=np.load(run);t=data['time']
    plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.titlesize':12,'axes.labelsize':11,'legend.fontsize':10})
    fig,axes=plt.subplots(2,2,figsize=(11,7.5),constrained_layout=True)
    styles=[dict(color='#4477AA',linestyle='-',label='Right'),dict(color='#EE6677',linestyle='--',label='Left')]
    for i,style in enumerate(styles):
        axes[0,0].plot(t,np.rad2deg(data['knee_coordinate'][:,i]),**style)
        axes[0,1].plot(t,data['rg'][:,2*i+1]-data['rg'][:,2*i],**style)
        axes[1,0].plot(t,data['contact_loads'][:,i].sum(axis=1),**style)
    axes[0,0].set(title='A  Knee motion',ylabel='Change from CAD reference (deg)')
    axes[0,1].set(title='B  Neural rhythm signal',ylabel='RG flexor minus extensor (mV)')
    axes[1,0].set(title='C  Ground contact',ylabel='Vertical foot force (N)')
    axes[1,1].plot(t,100*data['rig_force'][:,1]/(75*9.81),color='#228833',label='Actual vertical support')
    axes[1,1].set(title='D  External support',ylabel='Body weight supported (%)')
    for ax in axes.flat:
        ax.set_xlabel('Time (s)');ax.grid(alpha=.2);ax.legend(loc='best')
    fig.suptitle('AH001 hybrid SNS control, 75 kg, with vertical and angular support')
    stem=Path(run).stem
    fig.savefig(HERE/(stem+'_evidence.png'),dpi=180)
    fig.savefig(HERE/(stem+'_evidence.pdf'))
    plt.close(fig)
    (HERE/(stem+'_evidence.alt.txt')).write_text('Four time-series panels show knee changes from the CAD reference, RG flexor-minus-extensor potentials, vertical ground forces for each foot, and measured external vertical support. Right is solid blue, left is dashed red. Vertical and angular pelvis assistance are active. The figure does not by itself establish walking.')
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run',type=Path)
    main(parser.parse_args().run)
