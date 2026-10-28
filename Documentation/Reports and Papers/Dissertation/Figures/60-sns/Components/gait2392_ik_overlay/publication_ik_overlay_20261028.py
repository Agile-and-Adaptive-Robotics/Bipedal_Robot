from pathlib import Path
import os,sys,json
R=Path('D:/GitHub/Bipedal_Robot');S=R/'Code/MuJoCo_SNS/spinal';O=S/'reports_20261028_publication_flow/human_ik';sys.path.insert(0,str(S))
import numpy as np
import kine_ref as K
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'axes.titlesize':13})
z=np.load(O/'human_fresh.npz',allow_pickle=True);n=int(np.sum(z['t']>0))+1;t=z['t'][:n];q=z['q'][:n];contact=z['contact'][:n];neuro=z['neuro'][:n];ref=K.load_reference();grid=K.GRID
fig,axes=plt.subplots(4,2,figsize=(7.2,8.2),gridspec_kw={'height_ratios':[1,1,1,.65]});metrics={};columns={'cycle_percent':grid};m=t>=2
for side,col in [('r',0),('l',1)]:
 mean,periods,duty,cfrac,used,ncyc,on=K.sim_side(t,q,neuro,2.,side,contact)
 assert used and ncyc>0,'No contact-defined model cycles'
 cols=(3,4,5) if side=='r' else (8,9,10);cycs,ps,raw=K._cycles_from_onsets(t[m],t[m],on,cols,q[m]);pcontact=[]
 for start,end in zip(on[:-1],on[1:]):
  if not (.35<=end-start<=2.5):continue
  w=(t>=start)&(t<end)
  if w.sum()<8:continue
  pcontact.append(np.interp(grid,(t[w]-start)/(end-start)*100,(contact[w,col]>50).astype(float)))
 cp=np.mean(pcontact,axis=0);columns[side+'_model_stance_probability']=cp
 # Reconstruct reference stance from the exact GRF-defined cycle used by load_reference.
 ti,names,vals=K._read_mot(K.IK_MOT);tg,gn,gv=K._read_mot(K.GRF_MOT);force=gv[:,gn.index('ground_force_vy' if side=='r' else '1_ground_force_vy')-1];starts=K._loading_onsets(tg,force);b,e=K._first_pair_in_ik(starts,ti);cpr=(np.interp(grid,(tg-b)/(e-b)*100,force)>50).astype(float);columns[side+'_reference_stance']=cpr
 metrics[side]={'cycles':ncyc,'model_period_mean_s':float(np.mean(periods)),'model_contact_duty':duty,'reference_period_s':ref['T_'+side],'reference_duty':ref['duty_'+side]}
 for row,j in enumerate(['hip','knee','ankle']):
  ax=axes[row,col];mu=cycs[j];sd=np.std(raw[j],axis=0);rr=ref[side][j];ax.plot(grid,mu,color='#0072B2',lw=1.6,label='Model mean');ax.fill_between(grid,mu-sd,mu+sd,color='#0072B2',alpha=.14,lw=0);ax.plot(grid,rr,'--',color='black',lw=1.5,label='OpenSim IK');ax.set_title(('Right' if side=='r' else 'Left')+' '+j);ax.set_ylabel('Angle (deg)');ax.grid(alpha=.2);ax.set_xlim(0,100)
  columns[f'{side}_{j}_model_mean_deg']=mu;columns[f'{side}_{j}_model_sd_deg']=sd;columns[f'{side}_{j}_reference_deg']=rr;metrics[side]['rmse_'+j+'_deg']=float(np.sqrt(np.mean((mu-rr)**2)))
 ax=axes[3,col];ax.plot(grid,cp,color='#0072B2',lw=1.5);ax.plot(grid,cpr,'--',color='black',lw=1.5);ax.set_ylim(-.05,1.05);ax.set_xlim(0,100);ax.set_ylabel('Stance\nfraction');ax.set_xlabel('Gait cycle (%)');ax.grid(alpha=.2)
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False);fig.subplots_adjust(top=.94,bottom=.075,left=.12,right=.98,hspace=.5,wspace=.38);fig.savefig(O/'gait2392_ik_overlay.pdf');fig.savefig(O/'gait2392_ik_overlay.png',dpi=180);plt.close(fig)
np.savetxt(O/'gait2392_ik_overlay.csv',np.column_stack(list(columns.values())),delimiter=',',header=','.join(columns),comments='')
(O/'overlay_metrics.json').write_text(json.dumps(metrics,indent=2));(O/'overlay_methods.txt').write_text('Fresh saved round8_best3 replay, no retuning. Bilateral hip/knee/ankle angles in OpenSim sign convention; knee flexion negative. Walk starts 2 s. Contact >50 N, each leg own loading onset, cycles .35-2.5 s retained per kine_ref. No phase shift, mean removal, amplitude scaling, or fitted alignment. Curves mean +/- one cycle SD. Bottom strips: fraction of retained model cycles in stance at each phase, and reference stance from paired GRF. Each side normalized to its own cycle; strips do not establish interleg phase equivalence. Reference is subject01_walk1_ik.mot with subject01_walk1_grf.mot.\n');print(json.dumps(metrics),flush=True)
