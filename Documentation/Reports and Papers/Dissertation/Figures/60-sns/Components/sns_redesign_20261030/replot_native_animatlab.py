from pathlib import Path
import numpy as np,hashlib,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());D=R/'Documentation/Reports and Papers/Dissertation';P=R/'Neuromechanical_Models/Biped_2xCPG_wSubs/DataTool_8.txt';W=Path(__file__).resolve().parent/'output';W.mkdir(exist_ok=True)
names=P.open().readline().strip().split('\t');z=np.loadtxt(P,skiprows=1);t=z[:,1];w=(t>=2)&(t<=3.5)
assert np.isfinite(z[w]).all();assert np.allclose(np.diff(t[w]),.0002,atol=1e-8)
plt.rcParams.update({'font.family':'Arial','font.size':13,'axes.labelsize':13,'xtick.labelsize':13,'ytick.labelsize':13,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(6,1,figsize=(6.5,7.8),sharex=True,layout='constrained');blue='#0072B2';orange='#D55E00'
for i,(chans,scale,label) in enumerate([(['RG_L_stance','RG_R_stance'],1000,'Stance RG\nvoltage (mV)'),(['PF_L_hip_ext'],1000,'Left hip PF\next voltage\n(mV)'),(['MN_L_hip_ext'],1000,'Left hip MN\next voltage\n(mV)'),(['hip_L_deg','hip_R_deg'],1,'Hip angle\n(deg)'),(['knee_L_deg','knee_R_deg'],1,'Knee angle\n(deg)'),(['ankle_L_deg','ankle_R_deg'],1,'Ankle angle\n(deg)')]):
 for k,ch in enumerate(chans):axs[i].plot(t[w],scale*z[w,names.index(ch)],color=[blue,orange][k],ls=['-','--'][k],lw=1.1)
 axs[i].set_ylabel(label);axs[i].set_xlim(2,3.5);axs[i].grid(alpha=.15)
axs[-1].set_xlabel('Time (s)');fig.suptitle('Recorded AnimatLab phase-1 model\nLeft: blue solid | Right: orange dashed',fontsize=13)
name='animatlab_phase1_joint_angles';fig.savefig(W/(name+'.pdf'));fig.savefig(W/(name+'.png'),dpi=150);plt.close(fig)
alt='Native AnimatLab phase-1 recording, 2-3.5 s, without resampling or filtering. Rows show left/right stance RG voltages, left hip extensor PF voltage, left hip extensor MN voltage, and bilateral hip, knee, and ankle angles. Neural voltages converted from V to mV, joint channels already degrees. All within-window samples, including exact zero neural samples, are retained. No ground-force measurement or proof of walking is added.'
(W/(name+'.alt.txt')).write_text(alt,encoding='utf8')
for dest in ['60-sns','95-AppendixD']:
 folder=D/'Figures'/dest;dep=folder/'Deprecated/sns_redesign_20261030';dep.mkdir(parents=True,exist_ok=True)
 for ext in ['pdf','png','alt.txt']:
  target=folder/(name+'.'+ext)
  if target.exists() and not (dep/target.name).exists():shutil.copy2(target,dep/target.name)
  shutil.copy2(W/(name+'.'+ext),target)
info={'source':str(P.relative_to(R)),'sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'display_s':[2,3.5],'samples':int(w.sum()),'dt_s':.0002,'neural_zero_samples_retained':int(np.sum(z[w,names.index('RG_L_stance')]==0)),'model':'walk new new tester added 2 axis_phase1.asim','scope':'instrumented phase-1 configuration; not asserted original shared-RG model','units':'voltage V to mV; original recorded JointRotationDeg values unchanged'}
(W/'native_phase1_manifest.json').write_text(json.dumps(info,indent=2));print(json.dumps(info))
