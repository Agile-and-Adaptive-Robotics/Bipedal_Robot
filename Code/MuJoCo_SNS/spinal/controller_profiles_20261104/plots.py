"""Short-window neural replay plots; reference kinematics are not body results."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    z=np.load(a.input/'BEST_TRACE.npz');cfg=json.loads((a.input/'MANIFEST.json').read_text());final=json.loads((a.input/'FINAL.json').read_text())
    if 'configuration'in cfg:cfg=cfg['configuration']
    elif 'fit'in cfg:
        cfg=json.loads((Path(cfg['fit'])/'MANIFEST.json').read_text())['configuration'];cfg={**cfg,'phase':'D'}
    t=z['t'];keep=(t>=t[-1]-3*cfg['target_period_s']);x=t[keep]-t[keep][0];a.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none'})
    fig,ax=plt.subplots(4,1,figsize=(8.3,8),sharex=True,constrained_layout=True)
    for col in (1,3):ax[0].plot(x,z['rg'][keep,col],'-'if col==1 else'--',label=str(z['rg_names'][col]))
    names=z['pf_names'].tolist();selected=[]
    for fragment in ('HIP-F_r','HIP-E_r','BI_JOINT_RECT_r','BI_JOINT_HAM_r','PF_F1_r','PF_E1_r','PF_S1_r','PF_S2_r'):
        selected += [i for i,n in enumerate(names)if n.endswith(fragment)and i not in selected]
    if len(selected)<2:selected=np.argsort(np.ptp(z['pf'][keep],axis=0))[-4:].tolist()
    for i in selected[:4]:ax[1].plot(x,z['pf'][keep,i],label=names[i])
    acts=z['acts'].tolist()
    for base,style in (('iliacus_r','-'),('vas_int_r','--'),('bifemsh_r','-.'),('tib_ant_r',':')):
        i=acts.index(base);ax[2].plot(x,z['mn'][keep,i],style,label=base);ax[3].plot(x,z['desired'][keep,i],style,label=base)
    for y in ax:y.grid(alpha=.2);y.legend(loc='upper right',fontsize=8,ncol=2)
    ax[0].set_ylabel('RG (mV)');ax[1].set_ylabel('PF (mV)');ax[2].set_ylabel('MN / 5 mV');ax[3].set_ylabel('Target activation');ax[-1].set_xlabel('Time in displayed window (s)')
    fig.suptitle(cfg['family']+' — offline block '+cfg['phase']+(' passed'if final['accepted']else' failed')+'\nExpected afferents replayed; no body integration')
    for ext in ('png','pdf','svg'):fig.savefig(a.output/(cfg['family']+'_'+cfg['phase']+'_NEURAL.'+ext),dpi=180)
    plt.close(fig)
if __name__=='__main__':main()
