"""Plot original syn6 failure directly from raw arrays on the original time clock.

No model integration, averaging, smoothing, channel alignment, or decimation.
"""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=Path(r'D:\GitHub\Bipedal_Robot')
RAW=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
FIT=RAW/'waveform_edge_D_20261104/syn6'
SCORER=RAW/'source_snapshot_pf_design_20261105/pf_timing_audit.py'
FSA=ROOT/'Code/MuJoCo_SNS/spinal/fsa_results/fsa_backsolve.npz'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    manifest=json.loads((FIT/'MANIFEST.json').read_text())
    parent=Path(manifest['fit'])
    configuration=json.loads((parent/'MANIFEST.json').read_text())['configuration']
    tree=ast.parse(SCORER.read_text())
    target_fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='primitive_target')
    env=dict(np=np,left_phase_offset=lambda:0.)
    exec(compile(ast.Module(body=[target_fn],type_ignores=[]),str(SCORER),'exec'),env)
    period=configuration['target_period_s']
    start,end=2*period,4*period
    with np.load(FIT/'BEST_TRACE.npz') as raw, np.load(FSA) as fsa:
        keep=(raw['t']>=start)&(raw['t']<=end)
        time=raw['t'][keep]
        names=raw['pf_names'].tolist()
        actual={channel:np.clip(raw['pf'][keep,names.index(f'PF_S{channel}_r')]/5.,0.,1.)
                for channel in [1,2,5]}
        wanted={channel:env['primitive_target'](fsa,configuration,'r',time,channel-1)
                for channel in [1,2,5]}
    plt.rcParams.update({'font.family':'Arial','font.size':10.,'axes.labelsize':10.,
                         'axes.titlesize':10.,'xtick.labelsize':10.,'ytick.labelsize':10.,
                         'legend.fontsize':10.,'pdf.fonttype':42,'ps.fonttype':42,
                         'axes.linewidth':.8})
    fig,axes=plt.subplots(3,1,figsize=(14/2.54,11.7/2.54),sharex=True,sharey=True)
    for index,(ax,channel) in enumerate(zip(axes,[1,2,5])):
        ax.plot(time,actual[channel],color='#4477AA',lw=1.4,label='Actual PF transmission')
        ax.plot(time,wanted[channel],color='#222222',lw=1.3,ls=(0,(4,2)),label='Intended primitive')
        ax.axhline(.2,color='#888888',lw=.7,ls=(0,(1,2)))
        ax.set_xlim(start,end);ax.set_ylim(-.06,1.08);ax.set_yticks([0,.5,1.])
        ax.text(.015,.93,f'({chr(65+index)}) Right S{channel}',transform=ax.transAxes,
                ha='left',va='top',fontsize=10,fontweight='bold',
                bbox=dict(facecolor='white',edgecolor='none',alpha=.87,pad=.5))
        ax.spines[['right','top']].set_visible(False)
        ax.grid(axis='y',lw=.4,color='#dddddd')
        ax.tick_params(axis='both',length=3)
    axes[-1].set_xlabel('Original calibration time (s)')
    axes[-1].set_xticks([2.5,3.,3.5,4.,4.5])
    fig.supylabel('Normalized PF transmission',x=.015,fontsize=10)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.54,.995),
               ncol=1,frameon=False,handlelength=3.,borderaxespad=0.,labelspacing=.25)
    fig.subplots_adjust(left=.14,right=.985,bottom=.115,top=.84,hspace=.20)
    pdf=HERE/'PF_timing_authoritative.pdf'
    png=HERE/'PF_timing_authoritative.png'
    fig.savefig(pdf)
    fig.savefig(png,dpi=180)
    plt.close(fig)
    alt=('Three aligned panels compare original syn6 right PF channels S1, S2, and S5 '
         'against intended primitive transmission on the unchanged 2.466–4.932 s calibration '
         'clock. Actual transmission is solid blue; intended primitive is dashed black. '
         'The dotted gray line is the 0.2 active threshold. S1 and S2 show timing differences; '
         'S5 remains active. This raw two-cycle display is distinct from the nine-cycle '
         'averaged full audit establishing zero of twelve channel passes.')
    (HERE/'PF_timing_authoritative.alt.txt').write_text(alt+'\n',encoding='utf-8')
    provenance=dict(figure=str(pdf),figure_sha256=sha(pdf),png_sha256=sha(png),
                    size_cm=[14.,11.7],minimum_text_pt=10.,raw_source=str(FIT/'BEST_TRACE.npz'),
                    raw_sha256=sha(FIT/'BEST_TRACE.npz'),source_files={str(p):sha(p) for p in
                    [SCORER,FSA,FIT/'MANIFEST.json',parent/'MANIFEST.json',Path(__file__)]},
                    requested_display_window_s=[start,end],actual_display_sample_window_s=[float(time[0]),float(time[-1])],
                    samples_per_channel=int(len(time)),side='right',channels=[1,2,5],
                    transformation='clip(raw PF membrane voltage / 5 mV,0,1), matching original outgoing transmission; '
                    'intended target evaluated from the raw FSA right primitive using original target_start/period interpolation',
                    target_start_s=configuration['target_start_s'],target_period_s=period,
                    source_array='t and pf from BEST_TRACE.npz; time and r_source from fsa_backsolve.npz',
                    no_smoothing=True,no_averaging=True,no_decimation=True,no_phase_shift=True,
                    scoring_window_s=[2*period,(10+199/200)*period],
                    scoring_method='Full0/12 conclusion remains original nine-cycle mean audit, not recomputed from displayed two cycles.',
                    implementation='Matplotlib vector PDF with embedded TrueType text; exact140 mm width, no tight-bbox scaling.')
    (HERE/'PF_figure_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(json.dumps({k:provenance[k] for k in ['size_cm','samples_per_channel','actual_display_sample_window_s','figure_sha256']}))

if __name__=='__main__':
    main()
