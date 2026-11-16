"""Replot the accepted stage-2 aligned air recording, without changing samples."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
SOURCE = Path('D:/Users/Ben Bolen/ControllerTraining_20261104/hunt_pose_air_validation_v3/stage2/capture')
OUT = HERE / 'air_figure'
OUT.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
z = np.load(SOURCE / 'aligned_trace.npz', allow_pickle=False)
meta = json.loads((SOURCE / 'ALIGNED_TRACE_METADATA.json').read_text())
t = z['t_state']
assert np.array_equal(t, z['t_neural'])
window = (t >= 7) & (t <= 11)
names = z['neural_names'].tolist()
jn = z['joint_names'].tolist()
addresses = dict(zip(jn, z['joint_qpos_addresses'].tolist()))
plt.rcParams.update({'font.family':'Arial', 'font.size':10, 'axes.labelsize':10,
    'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,
    'pdf.fonttype':42,'path.simplify':False,'axes.linewidth':1.0})
fig, axes = plt.subplots(4, 1, figsize=(14/2.54, 15/2.54), sharex=True)
fig.subplots_adjust(left=.19,right=.975,top=.93,bottom=.08,hspace=.72)
blue,red = '#0077BB','#EE6677'
def line(ax, y, label, color, style='-'):
    ax.plot(t[window], y[window], label=label, color=color, ls=style, lw=1.15)
line(axes[0],z['neural_voltage'][:,names.index('RG_F_r')],'Right',blue)
line(axes[0],z['neural_voltage'][:,names.index('RG_F_l')],'Left',red,'--')
axes[0].set_ylabel('RG-F voltage\n(mV)')
axes[0].legend(loc='upper center',bbox_to_anchor=(.5,1.34),ncol=2,frameon=False)
line(axes[1],z['neural_voltage'][:,names.index('MN_vas_int_r')],'Vastus intermedius',blue)
line(axes[1],z['neural_voltage'][:,names.index('MN_bifemsh_r')],'Biceps femoris SH',red,'--')
axes[1].set_ylabel('Right MN voltage\n(mV)')
axes[1].legend(loc='upper center',bbox_to_anchor=(.5,1.50),ncol=1,frameon=False,handlelength=2.5)
for side,color,style in [('r',blue,'-'),('l',red,'--')]:
    line(axes[2],np.rad2deg(z['qpos'][:,addresses['knee_angle_'+side]]),side,color,style)
axes[2].set_ylabel('Knee angle\n(deg)')
line(axes[3],z['harness_wrench_anatomical'][:,1],'Suspension',blue)
line(axes[3],z['explicit_support_anatomical_N'][:,1],'Explicit support',red,'--')
axes[3].set_ylabel('Superior force\n(N)')
axes[3].legend(loc='upper center',bbox_to_anchor=(.5,1.32),ncol=2,frameon=False)
axes[3].set_xlabel('Time (s)')
for i,ax in enumerate(axes):
    ax.set_xlim(7,11)
    ax.grid(True,alpha=.2,lw=.5)
    ax.spines[['top','right']].set_visible(False)
    ax.text(-.22,1.08,chr(65+i),transform=ax.transAxes,fontweight='bold')
fig.savefig(OUT/'biw2lvar_air_aligned_authoritative.pdf')
fig.savefig(OUT/'biw2lvar_air_aligned_authoritative.png',dpi=180)
plt.close(fig)
(OUT/'biw2lvar_air_aligned_authoritative.alt.txt').write_text(
    'Four aligned panels from the accepted BiW2Lvar+ stage-2 suspended-air recording, '
    '7 to 11 seconds. A: bilateral RG flexor membrane voltages. B: right vastus '
    'intermedius and biceps femoris short-head motor-neuron voltages. C: bilateral '
    'native flexion-negative knee angles. D: measured superior suspension force and '
    'separate zero explicit support force, both in newtons. All saved samples are retained.\n')
prov = {'source_files':{str(p):sha(p) for p in [SOURCE/'aligned_trace.npz',
    SOURCE/'ALIGNED_TRACE_METADATA.json',SOURCE/'model.mjb',SOURCE.parent/'body.json']},
    'plotted_window_s':[7,11], 'sample_selector':'7 <= t_state <= 11',
    'number_plotted':int(window.sum()),'saved_samples_total':len(t),
    'first_last_plotted_s':[float(t[window][0]),float(t[window][-1])],
    'native_joint_sign':'Flexion negative; no sign reversal',
    'transformations':['qpos radians to degrees using numpy.rad2deg',
        'anatomical superior force column 1; already transformed by recording capture'],
    'resampling':False,'smoothing':False,'decimation':False,'phase_shift':False,
    'amplitude_scaling':False,'force_normalization':'None; raw newtons',
    'legend_and_channel_map':{'A':['RG_F_r','RG_F_l'],
        'B':['MN_vas_int_r','MN_bifemsh_r'],'C':['knee_angle_r','knee_angle_l'],
        'D':['harness_wrench_anatomical[:,1]','explicit_support_anatomical_N[:,1]']},
    'coordinate_addresses':addresses,'minimum_font_pt':10,'width_cm':14,
    'scoring_window_note':'Display only. Adopted support statistic uses legacy result arrays after 5 s; see AIR_EVIDENCE.md.',
    'compiled_model_sha256_matches_metadata':sha(SOURCE/'model.mjb')==meta['model_sha256'],
    'output_hashes':{p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file()}}
(OUT/'AIR_FIGURE_PROVENANCE.json').write_text(json.dumps(prov,indent=2)+'\n')
print(json.dumps({'samples':prov['number_plotted'],'window':prov['first_last_plotted_s'],
    'out':str(OUT),'model_matches':prov['compiled_model_sha256_matches_metadata']}))
