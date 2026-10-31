"""Readable question-slide views of the saved native AnimatLab trials."""
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_comparison import OUT
plt.rcParams.update({'font.family':'Arial','font.size':14,'axes.spines.top':False,'axes.spines.right':False})
def load(m,c,n='baseline'):
 z=np.load(OUT/m/c/n/'signals.npz');return z['data'],list(z['labels'])
def save(f,n):
 f.tight_layout();f.savefig(OUT/'figures'/(n+'.png'),dpi=200,bbox_inches='tight',pad_inches=.1);f.savefig(OUT/'figures'/(n+'.pdf'),bbox_inches='tight',pad_inches=.1);plt.close(f)
f,ax=plt.subplots(1,2,figsize=(9.2,3.7),sharey=True)
for i,m in enumerate(['legacy','bilateral']):
 for c,color,style,label in [('air','#999999','-.','Root fixed'),('drop','#4477AA','-','Free'),('partial50','#EE6677','--','50% upward force')]:
  a,l=load(m,c);ax[i].plot(a[:,1],a[:,l.index('root_y')],color=color,ls=style,label=label,lw=1.8)
 ax[i].set_title('Original controller' if m=='legacy' else 'II + contact + bilateral RG');ax[i].set_xlabel('Time (s)');ax[i].grid(alpha=.2)
ax[0].set_ylabel('Pelvis height (m)');ax[1].legend(fontsize=12,loc='upper right');save(f,'defense_support')
f,ax=plt.subplots(2,1,figsize=(9.2,3.9),sharex=True,sharey=True)
for i,m in enumerate(['legacy','bilateral']):
 for n,color,style,label in [('baseline','#999999','-','No pulse'),('deng2019_PF','#4477AA','-','PF pulse'),('deng2019_RG','#EE6677','--','RG pulse')]:
  a,l=load(m,'air',n);ax[i].plot(a[:,1],1000*a[:,l.index('L Hip MN ext')],color=color,ls=style,lw=1.5,label=label)
 ax[i].axvspan(2,2.1,color='#CCBB44',alpha=.3);ax[i].set_xlim(1.5,3.6);ax[i].set_ylabel('L hip MN (mV)');ax[i].set_title('Original controller' if m=='legacy' else 'II + contact + bilateral RG',fontsize=14)
ax[0].legend(ncol=3,fontsize=12,loc='upper right');ax[-1].set_xlabel('Time (s)');save(f,'defense_deng2019')
with (OUT/'phase_response.csv').open() as h:rows=list(csv.DictReader(h))
f,ax=plt.subplots(2,2,figsize=(9.2,4.2),sharex=True,sharey=True)
for i,c in enumerate(['air','treadmill']):
 for j,sign in enumerate([-1,1]):
  for m,color,style,marker,label in [('legacy','#4477AA','-','o','Original'),('bilateral','#EE6677','--','s','Newer')]:
   r=[r for r in rows if r['model']==m and r['condition']==c and float(r['current_nA'])*sign>0]
   r.sort(key=lambda x:float(x['phase']));ax[i,j].plot([float(x['phase']) for x in r],[float(x['phase_advance_first_cycles']) for x in r],color=color,ls=style,marker=marker,ms=4,label=label)
  ax[i,j].axhline(0,color='#BBBBBB',lw=.7);ax[i,j].set_title(c+' | '+('-0.5 nA' if sign<0 else '+0.5 nA'),fontsize=14);ax[i,j].grid(alpha=.15)
  if j==0:ax[i,j].set_ylabel('Advance / period')
  if i==1:ax[i,j].set_xlabel('Stimulus phase')
ax[0,0].legend(fontsize=12,loc='upper left');save(f,'defense_deng2022')
