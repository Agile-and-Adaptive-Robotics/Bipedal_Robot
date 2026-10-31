"""Fresh-process replay, scoped summary, and scientific plot for the fixture."""
import time
PROCESS_START=time.perf_counter()
import json,hashlib,sys,shutil,subprocess
from pathlib import Path
import numpy as np
import run_fixture as fixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

O=fixture.OUT
IMPORT_S=time.perf_counter()-PROCESS_START
results=json.loads((O/'RESULTS.json').read_text())
cfg=results['training_winner']['configuration']
plant=next(r['plant'] for r in results['results'] if r['case']=='nominal')
tick=time.perf_counter()
r,tr,_=fixture.run(cfg,plant,save=O/'fresh_process_nominal')
fresh_run_s=time.perf_counter()-tick
original=np.load(O/'nominal_feedback.npz')['trace']
fresh_difference=float(np.max(np.abs(original-tr)))
assert fresh_difference==0.
native=np.load(O/'nominal_native.npz')['trace'];fast=np.load(O/'nominal_repeat.npz')['trace']
parity=dict(q_rad=float(np.max(np.abs(native[:,1]-fast[:,1]))),neurons_mV=float(np.max(np.abs(native[:,4:14]-fast[:,4:14]))),
    pressure_kpa=float(np.max(np.abs(native[:,20:24]-fast[:,20:24]))),force_N=float(np.max(np.abs(native[:,24:28]-fast[:,24:28]))))
rows=[]
for case in ('nominal','heldout_slow_weak','heldout_fast_strong','heldout_push'):
    on=next(x for x in results['results'] if x['case']==case and x['feedback'])
    off=next(x for x in results['results'] if x['case']==case and not x['feedback'])
    rows.append(dict(case=case,feedback_on_rmse_deg=on['rms_angle_error_deg'],feedback_off_rmse_deg=off['rms_angle_error_deg'],
        reduction_percent=100*(1-on['rms_angle_error_deg']/off['rms_angle_error_deg']),
        no_hyperextension=on['no_hyperextension'] and off['no_hyperextension']))
summary=dict(scope='Fixed-base one-joint actuator-interface fixture; imposed PF waveform, not autonomous CPG or walking.',
    comparisons=rows,fresh_process_trace_max_difference=fresh_difference,fresh_process_import_s=IMPORT_S,
    first_run_in_fresh_process_s=fresh_run_s,
    first_run_note='Existing on-disk Numba cache used; includes first dispatch/cache loading and fixture setup; not an uncached compilation benchmark.',
    native_numba_parity=parity)
(O/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
snap=O/'source_snapshots';snap.mkdir(exist_ok=True)
for p in (Path(fixture.__file__),Path(fixture.__file__).with_name('actuator_interface.py'),Path(__file__)):
    content=p.read_bytes();(snap/(hashlib.sha256(content).hexdigest()+p.suffix)).write_bytes(content)

plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
blue,orange,gray='#0072B2','#D55E00','.5'
fig,axes=plt.subplots(5,2,figsize=(9.2,9.5),sharex=True,layout='constrained')
for c,case,title in [(0,'nominal','Nominal pressure response'),(1,'heldout_slow_weak','Held-out slower/weaker actuators')]:
    z=np.load(O/(case+'_feedback.npz'))['trace'];off=np.load(O/(case+'_ablated.npz'))['trace'];m=(z[:,0]>=2)&(z[:,0]<=8);t=z[m,0]
    axes[0,c].plot(t,np.degrees(z[m,3]),color='black',ls=':',label='Desired');axes[0,c].plot(t,np.degrees(z[m,1]),color=blue,label='Feedback on');axes[0,c].plot(t,np.degrees(off[m,1]),color=orange,ls='--',label='Feedback off')
    axes[1,c].plot(t,z[m,4],color=blue,label='Flexor');axes[1,c].plot(t,z[m,5],color=orange,ls='--',label='Extensor')
    axes[2,c].plot(t,z[m,12],color=blue,label='Flexor');axes[2,c].plot(t,z[m,13],color=orange,ls='--',label='Extensor')
    axes[3,c].plot(t,z[m,20],color=gray,ls=':',label='Command');axes[3,c].plot(t,z[m,22],color=blue,label='Actual')
    axes[4,c].plot(t,z[m,24],color=gray,ls=':',label='Requested');axes[4,c].plot(t,z[m,26],color=blue,label='Delivered')
    axes[0,c].set_title('('+('A' if c==0 else 'B')+') '+title.replace(' slower/','\nslower/'),loc='left',fontsize=11)
    for row in range(5):
        ax=axes[row,c];ax.grid(alpha=.18);ax.set_xlim(2,8)
        if c==0:ax.set_ylabel(['Joint angle\n(deg)','PF voltage\n(mV)','MN voltage\n(mV)','Flexor pressure\n(kPa)','Flexor tension\n(N)'][row])
    axes[-1,c].set_xlabel('Time (s)')
for row in range(5):
    axes[row,1].legend(*axes[row,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(1.01,.5),frameon=False,fontsize=10)
fig.savefig(O/'neural_bpa_fixture_preview.png',dpi=220);plt.close(fig)

net=json.loads((O/'NETWORK.json').read_text())
lines=['digraph BPA_fixture {','graph [rankdir=LR, bgcolor="white", label="One-joint SNS–BPA fixture: prescribed PF input; no RG or torque allocator", labelloc=t, fontsize=16];','node [shape=ellipse,fontname="Arial",fontsize=11];','edge [fontname="Arial",fontsize=9];']
for cluster,names,label in [('pf',fixture.NAMES[:2],'Prescribed pattern input'),('sens',fixture.NAMES[2:6],'Actual position/velocity encoders'),('mn',fixture.NAMES[6:],'Motor and reciprocal inhibitory neurons')]:
    lines.append('subgraph cluster_'+cluster+' { label="'+label+'";')
    for name in names:lines.append('"'+name+'" [style=filled,fillcolor="'+('#d5eef7' if cluster=='pf' else '#fff0cc' if cluster=='sens' else '#e0ebd5')+'"];')
    lines.append('}')
for e in net['connections']:
    inh=e['reversal_mV']<0
    lines.append('"'+e['source']+'" -> "'+e['destination']+'" [label="'+format(e['conductance_uS'],'.3g')+' µS",color="'+('#b2182b' if inh else '#2166ac')+'",arrowhead="'+('tee' if inh else 'normal')+'"];')
lines+=['actuator [shape=box,label="Command conversion → pressure state → BPA force\nMoment arms → fixed-base MuJoCo joint"];','MN_F -> actuator [style=dashed,label="normal output map"];','MN_E -> actuator [style=dashed];','actuator -> II_positive [style=dashed,label="joint state"];','actuator -> II_negative [style=dashed];','actuator -> Ia_positive [style=dashed];','actuator -> Ia_negative [style=dashed];','}']
(O/'network.dot').write_text('\n'.join(lines),encoding='utf-8')
dot=next((p for p in [Path(sys.executable).parent/'Library/bin/dot.exe',Path('C:/Program Files/Graphviz/bin/dot.exe')] if p.exists()),None)
if dot is not None:
    subprocess.run([str(dot),'-Tsvg',str(O/'network.dot'),'-o',str(O/'network.svg')],check=True)
    subprocess.run([str(dot),'-Tpng',str(O/'network.dot'),'-o',str(O/'network.png')],check=True)
print(json.dumps(summary,indent=2))
