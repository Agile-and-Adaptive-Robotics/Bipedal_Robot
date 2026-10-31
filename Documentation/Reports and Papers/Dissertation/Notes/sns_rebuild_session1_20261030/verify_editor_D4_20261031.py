"""Read saved trajectories; write coordination evidence and a figure, never simulate."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = Path(r'D:\GitHub\Bipedal_Robot')
D = R / 'Documentation/Reports and Papers/Dissertation'
N = Path(__file__).parent
O = N / 'editor_D4_verification'
O.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
protected = [D/'ProofFinal/chapters/60-sns.tex', D/'ProofFinal/chapters/95-AppendixD.tex', D/'ProofFinal/main.tex', D/'ProofFinal/thesis.bib', D/'ProofFinal/chapters/20-muscles.tex']
before = {str(p): sha(p) for p in protected}
legacy = R/'Code/MuJoCo_SNS/spinal/reports_20261028_publication_flow/stimulus_long_controls'
b = np.load(legacy/'baseline.npz')
p = np.load(legacy/'PF_plus50_2cycles.npz')
rows = json.loads((legacy/'stimulus_results.json').read_text())
row = next(x for x in rows if x['condition']=='PF_plus50_2cycles')
t = b['t']
def crossings(time, signal, threshold):
    return time[np.flatnonzero((signal[1:]>threshold)&(signal[:-1]<=threshold))+1]
events = {str(c):crossings(t,b['V'][:,c],-45).tolist() for c in (0,1,4,5,12,13)}
e = np.array(events['0']); f = np.array(events['1'])
combined = np.sort(np.r_[e[e>8],f[f>8]])
on = row['start_s']
display = (t>=on-2)&(t<=on+6)
xml = R/'Code/MuJoCo_SNS/spinal/w2l_mujoco/w2l_air.xml'
ankle_range = np.rad2deg(np.fromstring(ET.parse(xml).find(".//joint[@name='ankle_L']").get('range'),sep=' '))
legacy_metrics = {}
for name,z in [('baseline',b),('PF_plus50_2cycles',p)]:
    q = np.rad2deg(z['q'][:,2]); w=q[display]
    legacy_metrics[name] = dict(ankle_min_max_deg=[float(q.min()),float(q.max())],display_zero_upcrossings=int(np.sum((w[1:]>0)&(w[:-1]<=0))),display_duration_s=float(t[display][-1]-t[display][0]),fraction_outside_nominal_limits_by_over_2deg=float(np.mean((q<ankle_range[0]-2)|(q>ankle_range[1]+2))))

W = R/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031'
run = W/'runs/W2L_VirtualWalker_Newer'
base_folder = run/'canonical108_air_baseline_v2'
pf_folder = run/'canonical108_air_PF_minus10nA_100ms_v1'
base = np.load(base_folder/'trajectory.npz')
stim = np.load(pf_folder/'trajectory.npz')
result = json.loads((pf_folder/'RESULT.json').read_text())
source = result['source_sha256']
assert len(stim['neuron_names'])==108 and result['neurons']==108 and result['nonspiking_edges']==216
tt=stim['t']; pulse=result['neural_pulse']; onset=pulse['start_s']; end=onset+pulse['duration_s']
pre=tt<onset
assert np.array_equal(base['t'],tt)
assert np.array_equal(base['neuro'][pre],stim['neuro'][pre])
assert np.array_equal(base['angles'][pre],stim['angles'][pre])
model=ET.parse(pf_folder/'model.xml')
offsets=np.array([-10.,10.,0.,20.,55.,0.])
names=stim['joint_names'].tolist()
bounds=np.array([np.rad2deg(np.fromstring(model.find(".//joint[@name='"+name+"']").get('range'),sep=' '))+offset for name,offset in zip(names,offsets)])
native_xml=R/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim'
assert sha(native_xml)==source
window=(tt>=3.4)&(tt<=7.8)
canonical={}
for label,z in [('baseline',base),('PF_minus10nA_100ms',stim)]:
    canonical[label]=dict(hip_min_max_deg=np.stack([z['angles'][:,[0,3]].min(0),z['angles'][:,[0,3]].max(0)]).tolist(),hip_min_max_deg_display=np.stack([z['angles'][window][:,[0,3]].min(0),z['angles'][window][:,[0,3]].max(0)]).tolist(),hip_fraction_within_half_degree_upper_limit_full12s=[float(np.mean(z['angles'][:,j]>=bounds[j,1]-.5))for j in (0,3)],hip_fraction_within_half_degree_upper_limit_display=[float(np.mean(z['angles'][window,j]>=bounds[j,1]-.5))for j in (0,3)],knee_min_deg=z['angles'][:,[1,4]].min(0).tolist(),ankle_min_max_deg=np.stack([z['angles'][:,[2,5]].min(0),z['angles'][:,[2,5]].max(0)]).tolist())

plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'axes.titlesize':10,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(7,2,figsize=(6.5,7.2),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[.6,1,1,1,.9,.9,.9]})
labels=['Current\n(nA)','RG\n(mV)','Hip PF\n(mV)','Hip MN\n(mV)','Hip\n(deg)','Knee\n(deg)','Ankle\n(deg)']
ns=stim['neuron_names'].tolist()
blue='#0077BB'; orange='#EE7733'; gray='#aaaaaa'
for side,col in [('L',0),('R',1)]:
    for j,ax in enumerate(axes[:,col]):
        ax.spines[['top','right']].set_visible(False);ax.grid(color='#dddddd',lw=.4);ax.axvspan(onset,end,color='#E69F00',alpha=.22);ax.set_xlim(3.4,7.8);ax.set_ylabel(labels[j]if col==0 else '')
    axes[0,col].set_title('('+chr(65+col)+') '+('Left'if col==0 else'Right')+' leg',fontweight='bold')
    target=ns.index(side+' Hip PF ext')
    axes[0,col].step(tt[window],stim['neural_stimulus_nA'][window,target],where='post',color='#222222',lw=1.1)
    axes[0,col].set_ylim(-11,1);axes[0,col].set_yticks([-10,0])
    for rowi,prefix in [(1,'RG'),(2,'Hip PF'),(3,'Hip MN')]:
        for token,color,style in [('ext',blue,'-'),('flx',orange,'--')]:
            k=ns.index(side+' '+prefix+' '+token)
            axes[rowi,col].plot(tt[window],base['neuro'][window,k],style,color=gray,lw=1.1)
            axes[rowi,col].plot(tt[window],stim['neuro'][window,k],style,color=color,lw=1.0)
    for q,rowi in enumerate((4,5,6)):
        j=q+col*3
        axes[rowi,col].plot(tt[window],base['angles'][window,j],color=gray,lw=1.2)
        axes[rowi,col].plot(tt[window],stim['angles'][window,j],color='#222222',lw=1.0)
    for bound in bounds[col*3]:axes[4,col].axhline(bound,color='#884488',ls=':',lw=.8)
    axes[4,col].set_ylim(-23,25);axes[4,col].set_yticks([-18,0,20])
    axes[5,col].set_ylim(-2,60);axes[5,col].set_yticks([0,30,60])
    axes[6,col].set_ylim(-16,9);axes[6,col].set_yticks([-12,0,6])
    axes[-1,col].set_xlabel('Time (s)');axes[-1,col].set_xticks([4,5,6,7])
for j in (1,2,3):
    lo=min(x.get_ylim()[0]for x in axes[j]); hi=max(x.get_ylim()[1]for x in axes[j])
    for ax in axes[j]:ax.set_ylim(lo,hi)
fig.legend(handles=[Line2D([],[],color=blue,label='Extensor'),Line2D([],[],color=orange,ls='--',label='Flexor'),Line2D([],[],color=gray,label='Baseline'),Line2D([],[],color='#884488',ls=':',label='Hip limits')],loc='outside upper center',ncol=4,frameon=False,handlelength=1.1,columnspacing=.75)
stem=O/'W2L108_PF_BILATERAL_LIMITS'
for ext in ('pdf','png','svg'):fig.savefig(stem.with_suffix('.'+ext),dpi=240)
plt.close(fig)

inputs=[N/'EDITOR_D4_CURRENT.pdf',Path(r'C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf'),legacy/'baseline.npz',legacy/'PF_plus50_2cycles.npz',legacy/'stimulus_results.json',xml,W/'run_virtual_walker.py',native_xml,base_folder/'trajectory.npz',pf_folder/'trajectory.npz',pf_folder/'RESULT.json',pf_folder/'model.xml']
for oldname,current in [('60-sns_remote.tex','60-sns.tex'),('95-AppendixD_remote.tex','95-AppendixD.tex')]:
    old=D/'Archive/session1_sns_20261031/inputs_eb475ws4'/oldname
    new=D/'ProofFinal/chapters'/current
    diff=''.join(difflib.unified_diff(old.read_text(encoding='utf-8').splitlines(True),new.read_text(encoding='utf-8').splitlines(True),fromfile=str(old),tofile=str(new)))
    (O/('V5_VS_CAPTURED_REMOTE_'+current+'.diff')).write_text(diff,encoding='utf-8')
geometry=O/'human_geometry'
geometry.mkdir(exist_ok=True)
for filename in ['Gait2392_human_reference_oblique.png','Gait2392_human_reference_two_view.png','PROVENANCE.json']:
    shutil.copy2(R/'tmp/human_reference_static_20261031'/filename,geometry/filename)
after={str(p):sha(p)for p in protected}
assert before==after,'Working manuscript changed during verification; investigate concurrent edits.'
array_hashes={str(folder/(name+'.npz')):{k:dict(shape=list(z[k].shape),dtype=str(z[k].dtype),sha256=hashlib.sha256(np.ascontiguousarray(z[k]).tobytes()).hexdigest())for k in z.files}for folder,name,z in [(legacy,'baseline',b),(legacy,'PF_plus50_2cycles',p),(base_folder,'trajectory',base),(pf_folder,'trajectory',stim)]}
audit=dict(machine='easteregg2',simulation_reruns=0,working_manuscript_unchanged=True,working_source_hashes=after,input_hashes={str(p):sha(p)for p in inputs},saved_array_hashes=array_hashes,legacy=dict(channel_order=['RG_HC_ext','RG_HC_flx','RG_IN_ext','RG_IN_flx','PF_hip_HC_ext','PF_hip_HC_flx','PF_ka_HC_ext','PF_ka_HC_flx','PF_hip_IN_ext','PF_hip_IN_flx','PF_ka_IN_ext','PF_ka_IN_flx','MN_hip_ext','MN_hip_flx'],record_dt_s=float(np.median(np.diff(t))),onsets_at_minus45mV=events,same_cell_period_after8s=float(np.median(np.diff(e[e>8]))),alternating_E_F_interval_after8s=float(np.median(np.diff(combined))),nominal_ankle_range_deg=ankle_range.tolist(),metrics=legacy_metrics,saved_stimulus=row),canonical108=dict(source_sha256=source,pulse=pulse,anatomical_offset_deg=offsets.tolist(),joint_names=names,nominal_anatomical_bounds_deg=bounds.tolist(),metrics=canonical,prestimulus_neural_and_angle_arrays_identical=True,plot_window_s=[3.4,7.8],plot_has_no_smoothing_or_time_warp=True,minimum_font_pt=10,plot_native_inches=[6.5,7.2]),preview_hashes={str(stem.with_suffix('.'+ext)):sha(stem.with_suffix('.'+ext))for ext in ('pdf','png','svg')})
(O/'D4_VERIFICATION.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({'legacy':audit['legacy'],'canonical108':audit['canonical108'],'preview':str(stem.with_suffix('.png'))},indent=2))
