"""Scientific review plot from recorded neural/mechanical states of one trial."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from import_virtual_walker import HERE
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':10})
def draw(folder,start=2.,stop=4.):
 folder=Path(folder);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());names=z['neuron_names'].tolist();mask=(z['t']>=start)&(z['t']<=stop);t=z['t'][mask];fig,axes=plt.subplots(7,1,figsize=(6.6,9),sharex=True,layout='constrained');blue,orange='#0072B2','#D55E00'
 for ax,prefix,label in zip(axes[:3],['L RG','L Knee PF','L Knee MN'],['RG voltage\n(mV)','Knee PF voltage\n(mV)','Knee MN voltage\n(mV)']):
  for tail,c,ls in [('ext',blue,'-'),('flx',orange,'--')]:ax.plot(t,z['neuro'][mask,names.index(prefix+' '+tail)],color=c,ls=ls,lw=1.05,label='Extensor' if tail=='ext' else 'Flexor')
  ax.set_ylabel(label)
 for ax,k,label in zip(axes[3:6],[0,1,2],['Hip angle\n(deg)','Knee flexion\n(deg)','Ankle angle\n(deg)']):
  for side,offset,c,ls in [('Left',0,'.15','-'),('Right',3,'#AA3377','--')]:ax.plot(t,z['angles'][mask,k+offset],color=c,ls=ls,lw=1.1,label=side)
  ax.set_ylabel(label)
 for side,k,c,ls in [('Left',0,'.15','-'),('Right',1,'#AA3377','--')]:axes[6].plot(t,z['foot_vertical'][mask,k],color=c,ls=ls,lw=1.1,label=side)
 axes[6].set_ylabel('Foot vertical\nforce (N)');axes[-1].set_xlabel('Time (s)')
 for ax in axes:ax.grid(alpha=.17);ax.set_xlim(start,stop)
 h,l=axes[0].get_legend_handles_labels();h2,l2=axes[3].get_legend_handles_labels();fig.legend(h+h2,l+l2,ncol=4,loc='outside upper center',frameon=False)
 path=folder/'SHORT_NEURAL_MECHANICAL.png';fig.savefig(path,dpi=220);plt.close(fig);return path
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('folder');p.add_argument('--start',type=float,default=2.);p.add_argument('--stop',type=float,default=4.);a=p.parse_args();print(draw(a.folder,a.start,a.stop))
