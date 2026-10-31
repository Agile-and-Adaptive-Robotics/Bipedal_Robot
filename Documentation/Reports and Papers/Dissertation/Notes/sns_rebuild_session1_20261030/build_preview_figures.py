"""Review-only short-window figures; original data and active figures stay unchanged."""
import os,json,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());O=Path(__file__).parent/'previews';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':13,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':13,'pdf.fonttype':42})
blue,orange='#0072B2','#D55E00';manifest=[]
p=R/'Code/MuJoCo_SNS/spinal/campaigns/20260930/w2l_numpy_ref.npz';q=R/'Code/Matlab/SNS_Simscape/results/SNS_W2L_CPG_run_20260930.mat'
a=np.load(p,allow_pickle=True);b=loadmat(q);names=a['names'].tolist();t=a['t'];ts=b['t'].ravel();fig,axes=plt.subplots(2,1,figsize=(6.5,4.8),sharex=True,layout='constrained')
for ax,side,letter in zip(axes,('L','R'),'AB'):
 for family,short,color in (('ext','E',blue),('flx','F',orange)):
  mask=(t>=2)&(t<=5);m=(ts>=2)&(ts<=5)
  ax.plot(t[mask],a['V'][mask,names.index(side+' RG '+family)],color=color,lw=1.35,label=family.title()+' NumPy')
  ax.plot(ts[m],b['v'+side+short].ravel()[m],color=color,lw=1.05,ls='--',label=family.title()+' Simulink')
 ax.set_title('('+letter+') '+('Left' if side=='L' else 'Right')+' rhythm generator',loc='left',fontsize=13)
 ax.set_ylabel('Voltage relative\nto rest (mV)');ax.grid(alpha=.17);ax.set_xlim(2,5)
fig.legend(*axes[0].get_legend_handles_labels(),ncol=2,loc='outside upper center',frameon=False);axes[-1].set_xlabel('Time (s)');fig.savefig(O/'fig6_3_compact_engine_overlay.png',dpi=220);plt.close(fig)
manifest.append(dict(figure='fig6_3_compact_engine_overlay.png',sources=[str(p),str(q)],window_s=[2,5],scope='Same existing traces; no time shift, normalization, or new simulation. Grouped E/F by leg; 3-second view.'))
S=R/'Neuromechanical_Models/Walker_2_Layer_CPG'
def chart(name):
 p=S/name;header=p.open().readline().strip().split('\t');z=np.loadtxt(p,skiprows=1);manifest.append(dict(source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return header,z
nr,rg=chart('Rhythm Generator.txt');npf,pf=chart('L Hip PF.txt');nm,mn=chart('L Hip Sensory Motor.txt');nj,joints=chart('L Angles.txt')
print('MN_COLUMNS',nm,'JOINT_COLUMNS',nj)
diagnostic={name:{n:[float(z[:,i].min()),float(z[:,i].max())] for i,n in enumerate(names[2:],2)} for name,names,z in [('rg',nr,rg),('pf',npf,pf),('mn',nm,mn),('joints',nj,joints)]}
(O/'animatlab_column_ranges.json').write_text(json.dumps(diagnostic,indent=2))
import xml.etree.ElementTree as ET
model=ET.parse(S/'Walker_2_Layer_CPG_Standalone_modern.asim').getroot();types={e.findtext('ColumnName'):e.findtext('DataType') for e in model.iter('DataColumn')}
assert all(types[k]=='JointRotation' for k in ('hip_L','knee_L','ankle_L'))
fig,axes=plt.subplots(6,1,figsize=(6.5,8.3),sharex=True,layout='constrained')
for ax,names,z,channels,label in [(axes[0],nr,rg,['L RG ext','L RG flx'],'RG\n(mV)'),(axes[1],npf,pf,['L Hip PF ext','L Hip PF flx'],'Hip PF\n(mV)'),(axes[2],nm,mn,['L Hip MN ext','L Hip MN flx'],'Hip MN\n(mV)')]:
 mask=(z[:,1]>=1)&(z[:,1]<=2.5)
 for key,color,style in zip(channels,(blue,orange),('-','--')):ax.plot(z[mask,1],1000*z[mask,names.index(key)],color=color,ls=style,lw=1.1,label=key)
 ax.set_ylabel(label)
for ax,key in zip(axes[3:],('hip_L','knee_L','ankle_L')):
 mask=(joints[:,1]>=1)&(joints[:,1]<=2.5);ax.plot(joints[mask,1],np.degrees(joints[mask,nj.index(key)]),color='.2',lw=1.1);ax.set_ylabel(key.split('_')[0].title()+' angle\n(deg)')
for ax in axes:
 ax.grid(alpha=.17);ax.set_xlim(1,2.5);ax.tick_params(labelsize=16);ax.xaxis.label.set_fontsize(16);ax.yaxis.label.set_fontsize(16)
fig.legend(*axes[0].get_legend_handles_labels(),ncol=2,loc='outside upper center',frameon=False,prop={'size':16});axes[-1].set_xlabel('Time (s)',fontsize=16);fig.savefig(O/'animatlab_short_neural_mechanical_preview.png',dpi=220);plt.close(fig)
manifest.append(dict(figure='animatlab_short_neural_mechanical_preview.png',window_s=[1,2.5],joint_units='Converted from JointRotation radians to degrees; original joint signs retained',scope='Existing recorded left-leg channels; not a new simulation or matched platform-release experiment.'))
sp=R/'Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal4_gate2_traces.npz'
if sp.exists():
 z=np.load(sp);t=z['t'];mask=(t>=2)&(t<=4.5);fig,axes=plt.subplots(2,1,figsize=(6.5,3.7),sharex=True,layout='constrained')
 axes[0].plot(t[mask],z['v_ns'][mask],color=blue,lw=1.1);axes[0].set_ylabel('RG E–F readout\n(mV)');axes[0].set_title('(A) Graded network',loc='left',fontsize=13)
 for level,key,color in ((1,'cnt_e',blue),(0,'cnt_f',orange)):
  events=np.flatnonzero(np.diff(z[key])>0)+1;events=events[(t[events]>=2)&(t[events]<=4.5)];axes[1].vlines(t[events],level-.3,level+.3,color=color,lw=.6)
 axes[1].set_yticks([0,1],['RG F','RG E']);axes[1].set_ylim(-.6,1.6);axes[1].set_ylabel('Spike events');axes[1].set_title('(B) Untuned spiking variant',loc='left',fontsize=13)
 for ax in axes:ax.set_xlim(2,4.5);ax.grid(alpha=.17)
 axes[-1].set_xlabel('Time (s)');fig.savefig(O/'figD6_short_neural_preview.png',dpi=220)
 if '--pdf' in __import__('sys').argv:fig.savefig(O/'spiking_mirror_rhythm_short.pdf')
 plt.close(fig)
 manifest.append(dict(figure='figD6_short_neural_preview.png',source=str(sp),sha256=hashlib.sha256(sp.read_bytes()).hexdigest(),window_s=[2,4.5],scope='Neural-only record. Spike events from positive cumulative-count increments; no unrelated joint-angle traces appended.'))
(O/'PREVIEW_PROVENANCE.json').write_text(json.dumps(manifest,indent=2))
print('Created review-only NumPy/Simulink preview and audited AnimatLab channels')
