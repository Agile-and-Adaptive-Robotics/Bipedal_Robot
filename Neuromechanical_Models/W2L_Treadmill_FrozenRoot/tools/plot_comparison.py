from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_comparison import OUT,onsets
plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
def data(m,c,n='baseline'):
 z=np.load(OUT/m/c/n/'signals.npz');return z['data'],list(z['labels'])
def save(fig,name):
 fig.tight_layout();fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight',pad_inches=.12);fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.12);plt.close(fig)
summary=[]
for m in ['legacy','bilateral']:
 for c in ['air','treadmill','drop','partial50']:
  a,l=data(m,c);t=a[:,1]
  row=dict(model=m,condition=c,root_initial_y=float(a[0,l.index('root_y')]),root_min_y=float(a[:,l.index('root_y')].min()),root_final_y=float(a[-1,l.index('root_y')]),belt_displacement_m=float(a[-1,l.index('belt_x')]-a[0,l.index('belt_x')]))
  falls=np.where(a[:,l.index('root_y')]<.5)[0]
  row['first_root_below_0_5_m_s']=float(t[falls[0]]) if len(falls) else None
  if c=='treadmill':row['mean_belt_speed_driven_m_s']=row['belt_displacement_m']/9.5
  for side in ['L','R']:
   cols=[i for i,k in enumerate(l) if k.startswith('contact_') and ('_'+side) in k]
   row[side+'_contact_fraction']=float(np.mean(np.any(a[:,cols]>0,axis=1)))
  summary.append(row)
  f,ax=plt.subplots(3,2,figsize=(12,8),sharex=True)
  for j,joint in enumerate(['Hip','Knee','Ank']):
   for s,side in enumerate(['L','R']):
    for term,color,style in [('ext','#4477AA','-'),('flx','#EE6677','--')]:
     key=f'{side} {joint} MN {term}'
     if key in l:ax[j,s].plot(t,1000*a[:,l.index(key)],color=color,ls=style,lw=.9,label=term)
    ax[j,s].set_title(side+' '+joint);ax[j,s].set_ylabel('MN voltage (mV)');ax[j,s].legend(loc='upper right',fontsize=8)
  for x in ax[-1]:x.set_xlabel('Time (s)')
  f.suptitle(f'{m}: {c} | saved controller, no retuning',y=1.01);save(f,f'{m}_{c}_MN')
  f,ax=plt.subplots(3,1,figsize=(11,7),sharex=True)
  for j,joint in enumerate(['hip','knee','ankle']):
   for side,color,style in [('L','#4477AA','-'),('R','#EE6677','--')]:
    key=f'angle_{joint}_{side}'
    if key in l:ax[j].plot(t,np.rad2deg(a[:,l.index(key)]),color=color,ls=style,lw=.9,label=side)
   ax[j].set_ylabel(joint+' (deg)');ax[j].legend(loc='upper right');ax[j].grid(alpha=.2)
  ax[-1].set_xlabel('Time (s)');f.suptitle(f'{m}: {c} | native AnimatLab joint convention');save(f,f'{m}_{c}_joints')
 for c in ['air','treadmill','drop','partial50']:
  base,l=data(m,c)
  f,ax=plt.subplots(3,2,figsize=(12,8),sharex=True)
  for col,trial in enumerate(['deng2019_PF','deng2019_RG']):
   a,la=data(m,c,trial)
   for row,key in enumerate(['L RG ext','L Hip PF ext','L Hip MN ext']):
    ax[row,col].plot(base[:,1],1000*base[:,l.index(key)],color='#BBBBBB',lw=1,label='baseline')
    ax[row,col].plot(a[:,1],1000*a[:,la.index(key)],color='#4477AA',lw=.9,label=trial)
    ax[row,col].axvspan(2,2.1,color='#CCBB44',alpha=.4);ax[row,col].set_xlim(1.3,4.5);ax[row,col].set_ylabel(key+' (mV)')
   ax[0,col].legend(fontsize=8);ax[-1,col].set_xlabel('Time (s)')
  f.suptitle(f'{m}: {c} | -10 nA, 2.0–2.1 s; Deng 2019 protocol transferred to W2L',y=1.01);save(f,f'{m}_{c}_deng2019')
(OUT/'baseline_summary.json').write_text(json.dumps(summary,indent=2))
prc=[]
for plan in json.loads((OUT/'phase_plan.json').read_text()):
 if 'period' not in plan:continue
 m,c,T=plan['model'],plan['condition'],plan['period'];b,lb=data(m,c);bt=onsets(b,lb.index('L RG ext'))
 for path in (OUT/m/c).glob('deng2022_*'):
  meta=json.loads((path/'protocol.json').read_text());start=meta['pulse'][1];a,l=data(m,c,path.name);at=onsets(a,l.index('L RG ext'));k=int(np.sum(bt<=start))
  first=(bt[k]-at[k])/T if k<len(bt) and k<len(at) else None
  late=(bt[k+2]-at[k+2])/T if k+2<len(bt) and k+2<len(at) else None
  prc.append(dict(model=m,condition=c,phase=float(path.name.split('phase')[1]),current_nA=meta['pulse'][3]*1e9,duration_s=meta['pulse'][2],phase_advance_first_cycles=first,phase_advance_third_cycles=late,baseline_period_s=T,stationary=plan['stationary']))
with (OUT/'phase_response.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(prc[0]));w.writeheader();w.writerows(prc)
f,ax=plt.subplots(2,2,figsize=(11,8),sharex=True)
for row,c in enumerate(['air','treadmill']):
 for col,sgn in enumerate([-.5,.5]):
  for m,color,marker,style in [('legacy','#4477AA','o','-'),('bilateral','#EE6677','s','--')]:
   q=sorted([r for r in prc if r['model']==m and r['condition']==c and r['current_nA']==sgn],key=lambda r:r['phase'])
   ax[row,col].plot([r['phase'] for r in q],[r['phase_advance_first_cycles'] for r in q],color=color,marker=marker,ls=style,label=m)
  ax[row,col].axhline(0,color='.7',lw=.7);ax[row,col].set_title(f'{c}: {sgn:+g} nA');ax[row,col].set_ylabel('First-cycle advance / baseline period');ax[row,col].legend()
for x in ax[-1]:x.set_xlabel('Stimulus phase relative to L RG ext onset')
f.suptitle('Deng 2022 pulse protocol transferred to W2L: duration = 5% of measured period',y=1.01);save(f,'deng2022_phase_response')
f,ax=plt.subplots(1,2,figsize=(12,4),sharey=True)
for j,m in enumerate(['legacy','bilateral']):
 for c,color,style in [('air','#BBBBBB',':'),('treadmill','#228833','-.'),('drop','#4477AA','-'),('partial50','#EE6677','--')]:
  a,l=data(m,c);ax[j].plot(a[:,1],a[:,l.index('root_y')],color=color,ls=style,label=c)
 ax[j].set_title(m);ax[j].set_xlabel('Time (s)');ax[j].legend();ax[j].grid(alpha=.2)
ax[0].set_ylabel('Pelvis/root height (m)');save(f,'support_comparison')
print(json.dumps(summary,indent=2))
