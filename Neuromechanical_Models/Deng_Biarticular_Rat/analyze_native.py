"""Short-window plots and honest diagnostics from the bounded native rat run."""
import os,json,sys,hashlib
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as E
HERE=Path(__file__).resolve().parent;OUT=HERE/'runs/native_source_smoke'
path=OUT/'source_baseline.txt'
names=[' '.join(s.strip().split()) for s in path.open().readline().strip().split('\t')]
z=np.loadtxt(path,skiprows=1);z=z[np.any(z[:,2:]!=0.,axis=1)];z=z[(z[:,1]>=0)&(z[:,1]<=5.001)]
assert np.all(np.isfinite(z)) and np.all(np.diff(z[:,1])>=0)
def signal(name):return z[:,names.index(name)]
t=z[:,1];use=(t>=1.)
periods={}
for name in ('LH_RG EXT','LH_RG FLX','RH_RG EXT','RH_RG FLX'):
 s=signal(name)*1000;peaks,_=find_peaks(s[use],prominence=.5,distance=max(1,int(.10/np.median(np.diff(t)))))
 ts=t[use][peaks];periods[name]=dict(peaks=len(ts),period_s=float(np.median(np.diff(ts))) if len(ts)>1 else None,range_mV=float(np.ptp(s[use])))
period=periods['LH_RG EXT']['period_s'];end=4.8;start=max(1.,end-3*period) if period and periods['LH_RG EXT']['peaks']>=3 else 1.
m=(t>=start)&(t<=end)
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False})
blue,orange='#0072B2','#D55E00';fig,axes=plt.subplots(6,1,figsize=(6.5,8.1),sharex=True,layout='constrained')
for ax,channels,label in [(axes[0],('LH_RG EXT','LH_RG FLX'),'RG voltage\n(mV)'),(axes[1],('LH_Hip PF EXT','LH_Hip PF FLX'),'Hip PF voltage\n(mV)'),(axes[2],('LH_BFA MN EXT','LH_IP MN FLX'),'Hip MN voltage\n(mV)')]:
 for name,color,ls in zip(channels,(blue,orange),('-','--')):ax.plot(t[m],1000*signal(name)[m],color=color,ls=ls,lw=1.15)
 ax.set_ylabel(label)
for ax,key in zip(axes[3:],('LH_HipZ angle','LH_Knee angle','LH_AnkleZ angle')):
 ax.plot(t[m],np.degrees(signal(key)[m]),color='.2',lw=1.1);ax.set_ylabel(key.split('_')[1].split()[0]+' coordinate\n(deg)')
for ax in axes:ax.grid(alpha=.17);ax.set_xlim(start,end)
axes[0].plot([],[],color=blue,label='Extensor');axes[0].plot([],[],color=orange,ls='--',label='Flexor')
fig.legend(*axes[0].get_legend_handles_labels(),ncol=2,loc='outside upper center',frameon=False)
axes[0].set_title('Left leg: native suspended baseline',loc='left',fontsize=11)
axes[-1].set_xlabel('Time (s)');fig.savefig(HERE/'Native_baseline_neural_joint.png',dpi=220);plt.close(fig)
fig,axes=plt.subplots(4,1,figsize=(6.5,5.7),sharex=True,layout='constrained')
axes[0].set_title('Left-leg biarticular MN activity: native suspended baseline',loc='left',fontsize=11)
for ax,name in zip(axes[:3],('LH_BFP MN EXT','LH_RF MN FLX','LH_GA MN FLX')):
 ax.plot(t[m],1000*signal(name)[m],color=blue,lw=1.15);ax.set_ylabel(name.split()[0].replace('LH_','')+' MN\n(mV)')
axes[-1].plot(t[m],signal('Root height')[m],color='.2');axes[-1].set_ylabel('Root y\n(m)')
for ax in axes:ax.grid(alpha=.17);ax.set_xlim(start,end)
axes[-1].set_xlabel('Time (s)');fig.savefig(HERE/'Native_biarticular_activity.png',dpi=220);plt.close(fig)
fig,axes=plt.subplots(2,1,figsize=(6.5,4.4),sharex=True,sharey=True,layout='constrained')
for ax,side in zip(axes,('LH','RH')):
 for suffix,color,ls in zip(('EXT','FLX'),(blue,orange),('-','--')):ax.plot(t[m],signal(side+'_RG '+suffix)[m]*1000,color=color,ls=ls,lw=1.2,label=suffix)
 ax.set_ylabel(side+' RG\n(mV)');ax.set_ylim(-64,-54);ax.grid(alpha=.17);ax.set_xlim(start,end)
axes[0].legend(loc='upper right',frameon=False,ncol=2);axes[-1].set_xlabel('Time (s)')
fig.savefig(HERE/'Native_bilateral_RG_diagnostic.png',dpi=220);plt.close(fig)
root=E.parse(HERE/'Biarticular_net_Standalone.asim').getroot()
stimuli=[dict(name=e.findtext('Name'),type=e.findtext('Type'),enabled=e.findtext('Enabled'),start=e.findtext('StartTime'),end=e.findtext('EndTime'),target=e.findtext('TargetNodeID')) for e in root.findall('.//ExternalStimuli/Stimulus')]
metrics=dict(machine='easteregg2',status='Native AnimatLab baseline completed; no tuning/ground-walking claim',duration_recorded_s=float(t[-1]),valid_samples=len(t),
 source='Public 2020 Biarticular_net precursor, not verified final2023',physical_muscles=sum(e.findtext('Type')=='LinearHillMuscle' for e in root.iter('RigidBody')),
 neurons=len(list(root.iter('Neuron'))),root_frozen=root.find('.//Organism/RigidBody/Freeze').text=='True',
 root_y_range_m=[float(signal('Root height').min()),float(signal('Root height').max())],
 rhythm_diagnostics=periods,plot_window_s=[start,end],
 joint_ranges_deg={name:[float(np.degrees(signal(name)[use]).min()),float(np.degrees(signal(name)[use]).max())] for name in names if name.endswith('angle')},
 stimuli=stimuli,knee_guard_note='Source joint limits retained. Recorded angles are native coordinates from source initial pose, not yet anatomically calibrated extension-neutral angles.',
 raw_trace=str(path),raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(HERE/'NATIVE_BASELINE_RESULTS.json').write_text(json.dumps(metrics,indent=2))
np.savez_compressed(OUT/'native_trace.npz',trace=z,names=np.array(names))
print(json.dumps(metrics,indent=2))
if __name__=='__main__':pass
