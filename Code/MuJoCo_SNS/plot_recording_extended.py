"""Select named saved signals and plot them without rerunning a model."""
import argparse,fnmatch,json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('recording',type=Path);p.add_argument('--kind',choices=['joints','neurons','controls','activation','pressure','foot_forces','harness','harness_moments'],default='joints')
p.add_argument('--signals',nargs='+');p.add_argument('--indices',type=int,nargs='+');p.add_argument('--list',action='store_true')
p.add_argument('--out',type=Path);p.add_argument('--start',type=float,default=0.);p.add_argument('--end',type=float,default=float('inf'))
a=p.parse_args();z=np.load(a.recording,allow_pickle=True)
is_bpa=all(k in z for k in ('q','p','f','h','rgE'))
if a.kind=='joints':
    values=z['q'];unit='Angle (deg)'
    if 'key_joints' in z:names=z['key_joints'].tolist()
    elif is_bpa and values.shape[1]==6:names=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R']
    else:p.error('This joint recording needs an explicit coordinate-name adapter')
elif a.kind=='neurons':
    if 'voltage' in z:values=z['voltage'];names=z['cells' if 'cells' in z else 'names'].tolist()
    elif 'neuro' in z:values=z['neuro'];names=z['neuro_names'].tolist()
    elif is_bpa:values=z['rgE'][:,None];names=['RG_E_L']
    else:p.error('No named neural log is present')
    unit='Membrane/readout (mV)'
elif a.kind=='activation':
    if not ('passive_force' in z and 'acts' in z):p.error('Native activation plotting requires a mechanics diagnostic; main-run act is a command subset.')
    values=z['act'];names=z['acts'].tolist();unit='Native activation (0–1)'
elif a.kind=='controls':
    if 'acts' in z:values=z['ctrl'];names=z['acts'].tolist()
    elif 'key_acts' in z:values=z['act'];names=z['key_acts'].tolist()
    else:p.error('Commands were not recorded; the BPA compact log has --kind pressure instead')
    unit='Excitation command (0–1)'
elif a.kind=='pressure':
    if not is_bpa:p.error('Pressure requires a BPA compact recording')
    values=z['p'];unit='Pressure (kPa)';manifest=a.recording.with_suffix('.result.json')
    report=json.loads(manifest.read_text()) if manifest.is_file() else {}
    names=report.get('pressure_channels',[f'pressure_{i}' for i in range(values.shape[1])])
elif a.kind=='foot_forces' and is_bpa:
    values=z['f'];names=['foot_L_normal','toe_L_normal','foot_R_normal','toe_R_normal'];unit='Contact normal force (N)'
else:
    if not all(k in z for k in ('samples','columns')):p.error('Use the Hunt .assistance.npz recording for harness or three-axis foot forces')
    columns=z['columns'].tolist();prefixes={'foot_forces':('right_GRF','left_GRF'),'harness':('harness_F',),'harness_moments':('harness_M',)}[a.kind]
    selected=[i for i,name in enumerate(columns) if name.startswith(prefixes)]
    values=z['samples'][:,selected];names=[columns[i] for i in selected]
    unit='Moment (N·m)' if a.kind=='harness_moments' else 'Force (N)'
if values.ndim!=2 or values.shape[1]!=len(names):p.error('Saved names do not match the requested signal matrix')
if not names:p.error('No channels of the requested kind were recorded')
if a.list:
    print('\n'.join(f'{i}: {name}' for i,name in enumerate(names)));raise SystemExit(0)
if not a.out:p.error('--out is required unless --list is used')
indices=list(a.indices or [])
for pattern in a.signals or []:
    matches=[i for i,name in enumerate(names) if fnmatch.fnmatchcase(name,pattern)]
    if not matches:p.error('No saved channel matches '+pattern+'; use --list or record the missing neurons explicitly')
    indices.extend(matches)
indices=list(dict.fromkeys(indices or list(range(min(6,len(names))))))
if any(i<0 or i>=len(names) for i in indices):p.error('Channel index outside the saved matrix')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':11})
t=z['t'] if 't' in z else z['samples'][:,z['columns'].tolist().index('t')]
keep=(t>=a.start)&(t<a.end)
if not keep.any():p.error('No samples lie in the requested time range')
fig,ax=plt.subplots(figsize=(10,5));colors=['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377'];styles=['-','--',':','-.'];groups={}
for j,i in enumerate(indices):
    name=names[i];side=name[-2:].lower() if name.lower().endswith(('_r','_l')) else '';base=name[:-2] if side else name
    groups.setdefault(base,len(groups));group=groups[base]
    ax.plot(t[keep],values[keep,i],label=name,color=colors[group%len(colors)],linestyle='--' if side=='_l' else '-' if side else styles[group%len(styles)],
        marker='s' if side=='_l' else 'o' if side else ['o','s','^','d'][group%4],markersize=3,markevery=max(1,int(keep.sum())//14))
ax.set(xlabel='Time (s)',ylabel=unit);ax.grid(alpha=.25);ax.legend(fontsize=10,ncol=2);fig.tight_layout()
a.out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.out,dpi=180);plt.close(fig);print(a.out.resolve())
