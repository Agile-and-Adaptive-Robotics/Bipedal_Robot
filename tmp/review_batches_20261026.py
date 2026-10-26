from pathlib import Path
import shutil, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(r'D:\Github\Bipedal_Robot')
base=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023'
inbox=base/'REVIEW_INBOX'
out=inbox/'REVIEWED_20261026'
patch=out/'patches'; figs=out/'figures'; evidence=out/'evidence'
for p in (patch,figs,evidence): p.mkdir(parents=True,exist_ok=True)
v=base/'Revision_E/review/targeted_reconciliation'
egg=inbox/'easteregg2/20261026_0730_sns-torque-evidence'
eb=inbox/'EB475WS4/20261026_0700_sns-demos-ablation'
deng=inbox/'easteregg2_2/20261026_deng-original-replication'
for name in ['P01_replacement.tex','P01_companion_table.tex','P04_figure61_block.tex','P05_controller_history.tex','P06_appendix_conversion_section.tex','P07_conversion_summary.tex']:
    shutil.copy2(v/name,patch/name)
table=(patch/'P01_companion_table.tex').read_text(encoding='utf-8')
rows=table.splitlines()
a=next(i for i,s in enumerate(rows) if s.startswith('SNS-Simscape library'))
b=next(i for i,s in enumerate(rows) if s.startswith('Gait2392/SNS'))
rows[a],rows[b]=rows[b],rows[a]
(patch/'P01_companion_table.tex').write_text('\n'.join(rows)+'\n',encoding='utf-8')

def write(name,s): (patch/name).write_text(s.strip()+'\n',encoding='utf-8')
write('EB_feedback_ablation.tex',r'''
At fixed controller and drive settings, removal of contact and Ib force-feedback pathways reduced stepping in the supported W2L model from 33 left and 28 right cycles to one left and no right cycles. Forward displacement decreased from \qty{1.16}{\meter} to \qty{0.09}{\meter}. These measures cover the \qty{25.499}{\second} scored window after settling within each \qty{30}{\second} simulation. The harness configuration was unchanged, although its measured support contribution differed between conditions. Ia and II pathways remained active. This comparison demonstrates the contribution of the removed pathways at this operating point; it uses the standard muscle-actuated W2L model.
\begin{figure}[t]
\centering
\includegraphics[width=0.85\textwidth]{ablation_pair_metrics_20261026.pdf}
\caption[Contact and force-feedback ablation]{Contact and Ib force-feedback ablation at fixed controller and drive settings in the supported W2L model. Counts and displacement refer to the post-settling scored window. The retained Ia and II pathways distinguish this test from complete deafferentation.}
\label{fig:sns_ablation_pair}
\end{figure}
''')
write('EG_controlled_tests.tex',r'''
\subsection{Controlled Perturbation and Pathway Tests}\label{sec:w2l_perturbation}
The experiment drivers apply neural stimuli, external forces, pathway deletions, and body-mass changes to a fixed controller configuration. Replaying all stored controller fields gives 17 left and 16 right strides in the control condition. The corresponding counts are 16/16 after crossed-heel deletion, 13/13 after Ia-chain deletion, 16/14 after Renshaw-cell deletion, 17/15 after heel-cascade deletion, and 13/12 after commissural deletion. Removing Ib load feedback stops stepping in both legs. Thus, Ib feedback is necessary for stepping at this tested operating point; the other deletions retain stepping with varying counts.

The neural stimuli and external-force tests also retain bilateral stepping. Counts alone do not establish balance recovery: the stronger stimulus and both push conditions exceed the \ang{45} tilt criterion. Across ten body-mass conditions, stride counts span 32--38 per leg and the lowest recorded pelvis height is \qty{0.8981}{\meter}. Nine conditions exceed the tilt criterion; one meets the full survival gate. These comparisons demonstrate controlled changes in neural and mechanical conditions while distinguishing continued stepping from postural stability.

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{w2l_perturbation_fully_applied.pdf}
\caption[Controlled W2L perturbation and pathway tests]{Supported W2L tests with all stored controller fields applied. The control produces 17/16 left/right strides. Ib deletion stops stepping; commissural deletion retains 13/12 strides. Stimulus and push trials retain stepping, although several exceed the tilt criterion. The controller replay does not constitute a new optimization with the corrected parameter forwarding.}
\label{fig:w2l_perturbation}
\end{figure}
''')
pulley=(egg/'patches/P04_torque_pulley_section.tex').read_text(encoding='utf-8')
pulley=pulley.replace('the tackle halves the insertion-side force','the two-part reversed tackle halves the insertion-side force')
write('EG_pulley_discussion.tex',pulley)
write('EG_extensor_merged_figure.tex',r'''
\begin{figure}[p]
\centering
\makebox[\textwidth][l]{\textbf{(A)}}\par
\includegraphics[width=\textwidth]{Ext10mm_52cm-eps-converted-to.pdf}\par\medskip
\makebox[\textwidth][l]{\textbf{(B)}}\par
\includegraphics[width=0.52\textwidth]{KneeExt_10mm_passthru.png}
\caption{Measured versus expected isometric torque using a $\diameter$\qty{10}{mm}, $l_{0}=\qty{51.8}{\cm}$ extensor BPA (A), and the corresponding routing configuration (B). In (A), blue dots show measured results, the dashed black line is the original prediction, the solid magenta line is the optimized prediction, and yellow dots show the hybrid calculation. In (B), the BPA passes between the femoral condyles close to the instantaneous center of rotation.}
\label{fig:Ext10mm_52cm}\label{fig:Ext_passthru}
\end{figure}
''')
write('EG_physical_figures_captionless_OPTION.tex',r'''
\begin{figure}[p]
\centering
\refstepcounter{figure}\label{fig:extensor_physical_final}
\includegraphics[width=0.97\textwidth]{extensor_20mm_physical.pdf}
\end{figure}
\begin{figure}[p]
\centering
\refstepcounter{figure}\label{fig:flexor_physical_final}
\includegraphics[width=0.97\textwidth]{flexor_20mm_physical.pdf}
\end{figure}
''')
write('DENG_original_engine.tex',r'''
The archived rat-hindlimb model of \citet{deng_neuromechanical_2019} provides an original-engine demonstration of the experimental workflow. The model was exported to standalone AnimatLab simulation form and instrumented with neural-voltage charts. In a \qty{10.1}{\second} run without external neural drive, rhythm-generator extensor and flexor populations alternate with an approximately \qty{0.23}{\second} period, while the two legs fire in phase. The baseline retains the archived model's mechanical torque input.

A \qty{2}{\nano\ampere} current was applied to the left rhythm-generator flexor from \SIrange{2}{2.1}{\second}, using the stimulation timing reported by \citet{deng_biarticular_2023}. A flexor burst begins at \qty{2.000}{\second}, compared with \qty{2.433}{\second} in the baseline. The post-stimulus period is \qty{0.224}{\second}, compared with \qty{0.221}{\second} in the baseline. Applying the same current and timing to the knee pattern-formation flexor produces a post-stimulus rhythm-generator phase difference of \qty{-2.3 \pm 2.4}{\percent} of the baseline period. An additional, longer knee-pattern-formation stimulus from \SIrange{2}{2.5}{\second} illustrates the motor-output response in Figure~\ref{fig:deng_original_perturbation}D. These are observations from the archived antagonist model in its original engine, rather than a quantitative validation of the Simulink translation or a reproduction of the later biarticular model's results.
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{deng_original_reviewed.pdf}
\caption[Archived Deng-model neural perturbations]{Neural perturbations in the archived antagonist model, simulated in its original AnimatLab engine. Solid lines show the baseline; dashed lines show stimulated conditions. (A, B)~Left rhythm-generator flexor and extensor responses to rhythm-generator flexor stimulation. (C)~Left rhythm-generator extensor response to knee pattern-formation flexor stimulation. The current is \qty{2}{\nano\ampere} from \SIrange{2}{2.1}{\second} in A--C. (D)~Left knee extensor motoneuron response to an extended knee pattern-formation flexor stimulus from \SIrange{2}{2.5}{\second}. Shading marks each panel's actual stimulus window.}
\label{fig:deng_original_perturbation}
\end{figure}
''')
(out/'deng_biarticular_2023.bib').write_text('''@incollection{deng_biarticular_2023,
 author = {Deng, Kaiyu and Hunt, Alexander J. and Chiel, Hillel J. and Quinn, Roger D.},
 title = {Biarticular Muscles Improve the Stability of a Neuromechanical Model of the Rat Hindlimb},
 booktitle = {Biomimetic and Biohybrid Systems},
 series = {Lecture Notes in Computer Science},
 volume = {14158}, pages = {20--37}, year = {2023},
 publisher = {Springer}, address = {Cham},
 doi = {10.1007/978-3-031-39504-8_2},
 url = {https://link.springer.com/chapter/10.1007/978-3-031-39504-8_2}
}
''',encoding='utf-8')
for name in ['ablation_pair_metrics_20261026.pdf','fig_kneereflex_v4.pdf','fig_beer_v4.pdf']:
    shutil.copy2(eb/'figures'/name,figs/name)
for name in ['extensor_20mm_physical.pdf','flexor_20mm_physical.pdf','FlxGrp_shared_axes.eps']:
    shutil.copy2(egg/'figures'/name,figs/name)
for name in ['w2l_perturbation_fully_applied.pdf']:
    shutil.copy2(base/'Revision_E/figs/60-sns'/name,figs/name)
for name in ['Ext10mm_52cm-eps-converted-to.pdf','KneeExt_10mm_passthru.png']:
    shutil.copy2(base/'Revision_E/figs/30-torque'/name,figs/name)
shutil.copytree(root/'tmp/batch_deng_data',evidence/'deng_raw',dirs_exist_ok=True)

plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10})
fig,axs=plt.subplots(4,1,figsize=(6.4,8.3),sharex=True)
spec=[('run_P1b_rgflexor_2p1',3,'RG flexor: RG stimulus',2.1),('run_P1b_rgflexor_2p1',2,'RG extensor: RG stimulus',2.1),('run_P2b_kneepfflexor_2p1',2,'RG extensor: knee PF stimulus',2.1),('run_P2_kneepfflexor',6,'Knee extensor MN: extended knee PF stimulus',2.5)]
def load(run):
    a=np.genfromtxt(evidence/'deng_raw'/run/'deng_rg_mn.txt',delimiter='\t',skip_header=1)
    a=a[np.isfinite(a).all(axis=1)]
    return a[(a[:,1]>=1.6)&(a[:,1]<=3.6)]
b=load('run_B0_baseline')
for i,(ax,(run,col,title,end)) in enumerate(zip(axs,spec)):
    a=load(run)
    ax.plot(b[:,1],b[:,col]*1000,color='#4477AA',lw=1,label='Baseline')
    ax.plot(a[:,1],a[:,col]*1000,color='#EE6677',ls='--',lw=1,label='Stimulated')
    ax.axvspan(2,end,color='0.9',zorder=0)
    ax.set_title(title,loc='left',fontsize=10,pad=7)
    ax.text(-.19,1.12,'('+chr(65+i)+')',transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax.set_ylabel('Voltage (mV)')
    if i<3:
        values=np.r_[b[:,col],a[:,col]]*1000
        margin=max(0.5,(values.max()-values.min())*.12)
        ax.set_ylim(values.min()-margin,values.max()+margin)
    else:
        ax.set_ylim(-105,5)
    ax.spines[['top','right']].set_visible(False)
axs[0].legend(loc='upper center',bbox_to_anchor=(.5,1.75),ncol=2,frameon=False,fontsize=10)
axs[-1].set_xlabel('Time (s)'); axs[-1].set_xlim(1.6,3.6)
fig.subplots_adjust(left=.20,right=.98,top=.88,bottom=.07,hspace=.75)
fig.savefig(figs/'deng_original_reviewed.pdf'); fig.savefig(figs/'deng_original_reviewed.png',dpi=180)
(figs/'deng_original_reviewed.alt.txt').write_text('Four neural-voltage panels compare baseline solid traces and stimulated dashed traces in the original AnimatLab Deng antagonist model. A and B show rhythm-generator flexor and extensor responses to RG flexor stimulation; C shows RG extensor response to knee pattern-formation flexor stimulation. A-C use 2.0-2.1 s stimuli. D shows knee extensor motor-neuron response during extended 2.0-2.5 s knee PF stimulation. Shaded regions show those distinct stimulus windows.',encoding='utf-8')
manifest=[]
for source in [eb,egg,deng]:
    for p in source.rglob('*'):
        if p.is_file(): manifest.append({'file':str(p.relative_to(inbox)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(out/'incoming_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(out)
