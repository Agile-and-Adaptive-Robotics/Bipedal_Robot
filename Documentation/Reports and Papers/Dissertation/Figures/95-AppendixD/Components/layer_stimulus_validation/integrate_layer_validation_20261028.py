from pathlib import Path
import json,shutil,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('D:/Github/Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';A=D/'Archive/flow_and_validation_20261028';P=D/'ProofFinal/chapters';F=D/'Figures'
before=A/'sources_before';before.mkdir(exist_ok=True)
for n in ['10-introduction','15-background','50-database','60-sns','70-overall','80-conclusion','95-AppendixD']:
 shutil.copy2(P/(n+'.tex'),before/(n+'.tex'))
protected={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'20-muscles.tex',P/'30-torque.tex',P/'94-AppendixC.tex',R/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex']}
(A/'protected_hashes.json').write_text(json.dumps(protected,indent=2))
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
rows=json.loads((A/'stimulus_controls/stimulus_results.json').read_text());base=np.load(A/'stimulus_controls/baseline.npz')
fig,axs=plt.subplots(3,3,figsize=(8.4,7.1),sharex=True)
for i,(name,label) in enumerate([('RG_plus50_1','RG: +50 nA, 1 s'),('PF_plus50_0.5','PF: +50 nA, 0.5 s'),('MNflex_minus500_0.5','MN: −500 nA, 0.5 s')]):
 row=next(r for r in rows if r['condition']==name);z=np.load(A/'stimulus_controls'/f'{name}.npz');t=z['t']-row['start_s'];mask=(t>=-2)&(t<=6)
 for j,(key,col,scale,ylabel) in enumerate([('V',0,1,'RG extensor, mV'),('motor',1,1,'Flexor MN activation'),('q',0,180/np.pi,'Hip angle, deg')]):
  ax=axs[i,j];ax.plot(t[mask],base[key][mask,col]*scale,'--',color='.35',lw=1.4,label='Baseline');ax.plot(t[mask],z[key][mask,col]*scale,color='#0072B2',lw=1.25,label='Stimulus');ax.axvspan(0,row['duration_s'],color='#EECC66',alpha=.35);ax.set_ylabel(ylabel);ax.grid(alpha=.15)
 axs[i,0].set_title(label,loc='left',fontweight='bold')
for ax in axs[-1]:ax.set_xlabel('Time from stimulus onset, s')
axs[0,2].legend(loc='upper right',fontsize=9);fig.tight_layout()
dest=F/'95-AppendixD';dest.mkdir(exist_ok=True)
fig.savefig(dest/'layer_stimulus_validation.pdf',bbox_inches='tight');fig.savefig(A/'layer_stimulus_validation.png',dpi=160);plt.close(fig)
comp=dest/'Components/layer_stimulus_validation';comp.mkdir(parents=True,exist_ok=True)
shutil.copy2(__file__,comp/Path(__file__).name)
with (comp/'stimulus_results.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(comp/'alt-text.txt').write_text('Three matched baseline/stimulus comparisons show RG voltage, flexor motor activation, and suspended-model hip angle. RG stimulation shifts subsequent rhythm; PF and direct MN interventions suppress output without shifting the RG clock. Gold shading marks stimulus duration.')
shutil.copy2(A/'human_ik/gait2392_ik_overlay.pdf',F/'60-sns/gait2392_ik_overlay.pdf')
ikcomp=F/'60-sns/Components/gait2392_ik_overlay';ikcomp.mkdir(parents=True,exist_ok=True)
for n in ['gait2392_ik_overlay.csv','overlay_metrics.json','overlay_methods.txt']:shutil.copy2(A/'human_ik'/n,ikcomp/n)
shutil.copy2(R/'tmp/publication_ik_overlay_20261028.py',ikcomp/'publication_ik_overlay_20261028.py')
(ikcomp/'alt-text.txt').write_text('Right and left hip, knee, and ankle cycles compare the supported human model with Gait2392 inverse kinematics. Mean model motion and cycle variability show substantial hip and knee timing differences. Bottom strips compare foot loading fractions with reference stance and swing.')
changes=[]
def edit(n,fn):
 p=P/(n+'.tex');s=p.read_text(encoding='utf-8');old=s;s=fn(s);p.write_text(s,encoding='utf-8');changes.append({'file':n,'changed':s!=old})
def sub(s,a,b):
 if a not in s:raise ValueError('Missing anchor: '+a[:100])
 return s.replace(a,b)
def sns(s):
 for a,b in [('selected functional-verification and proof-of-concept demonstrations of the resulting tools','tests of the circuits, interfaces, and coupled models'),('Selected neural ports','Selected neural implementations'),('one tonic-drive Simulink realization reached a fixed point','one Simulink realization reached a steady, nonoscillating state under constant neural drive'),('Under tonic drive, one Simulink circuit realization reached a fixed point','Under constant neural drive, one Simulink circuit realization reached a steady, nonoscillating state'),('a numerical failure sentinel','a reserved numerical failure value')]:s=sub(s,a,b)
 s=sub(s,'a pruned version, a six-synergy version, and a W2L-layout version','a version with selected connections removed, a version organized around six grouped muscle-activation patterns, and a version adapting the W2L controller organization')
 s=sub(s,'Free network gains were searched with a tree-structured Parzen estimator, and each candidate was evaluated through simulation.','Free network gains, which scale connection strengths and sensory contributions, were searched with a tree-structured Parzen estimator. This sequential search uses earlier simulation results to propose new parameter combinations; each candidate is then evaluated through simulation.')
 a=s.index('The library assembles Ia stretch');b=s.index('\\begin{figure}[p]',a)
 s=s[:a]+r'''The library assembles Ia stretch excitation, reciprocal inhibition, and Ib autogenic inhibition in a reduced-order antagonist knee model (Figure~\ref{fig:sns_reflex_current}C). The knee comparison starts from zero muscle activation under sustained neural drive and a constant \qty{0.5}{\newton\meter} flexion load. It tests startup behavior, rather than a brief tendon-tap response. At the end of the five-second window, the Simulink and NumPy implementations reach \ang{44.4} and \ang{43.0}, respectively; the maximum angular difference is \ang{3.6}. The implementations show similar startup behavior with a measurable numerical discrepancy; this comparison does not validate transient reflex recovery.

A second demonstration tests feedback in an arm holding a beverage container as it is loaded (Figure~\ref{fig:sns_beverage_current}). With the reflex enabled, arm sag was \ang{2.00}, compared with \ang{7.85} when it was disabled. As the forearm sags, biceps stretch feedback excites the biceps motoneuron and reciprocally inhibits the triceps motoneuron. The separate Ib pathways provide autogenic inhibition of each muscle. The reduced sag demonstrates a corrective effect of the implemented sensory pathways in the coupled neural and actuator blocks.

'''+s[b:]
 s=sub(s,'\\subsection{Reflex Response}','\\subsection{Coupled Block-Model Tests}')
 marker='\\subsection{Conversion Checks and Muscle-Force Agreement}'
 add=r'''\subsection{Layer-Specific Stimulation}\label{sec:layer_stimulus_main}

A missing motor burst can arise from a disturbance of the rhythm generator or from suppression farther downstream. I applied stimuli separately to the RG, PF, and MN layers of a reduced Deng-derived circuit \citep{deng_neuromechanical_2019,nourse_sns_2023}. A resetting deletion suppresses expected motor activity and shifts the subsequent rhythm; a non-resetting deletion suppresses activity while preserving the original schedule. Matched baseline and stimulus runs distinguish the rhythm-generator clock from the motor output.

A one-second RG stimulus shifted both RG and motor timing by 0.307 cycles. A shorter PF stimulus and direct MN inhibition suppressed flexor output without a lasting phase shift. A longer PF stimulus shifted motor timing by 0.208 cycles while leaving the RG clock unchanged. Thus, a motor phase shift alone did not establish resetting of the rhythm generator. The smaller-current sweep produced no observable deletions. Appendix~\ref{sec:app_layer_stimulus} gives the protocols, all stronger-current controls, and representative neural and joint responses. These are controlled tests of a translated circuit coupled to a suspended model, rather than a replication of the original rat's locomotion.

'''
 s=sub(s,marker,add+marker)
 # Interface checks precede the larger coupled controller comparison.
 a=s.index('\\subsection{Closed-Loop Model Coupling}');b=s.index('\\subsection{Feedback Contribution',a);bridge=s[a:b];s=s[:a]+s[b:];s=sub(s,'\\subsection{Testing Alternative Controllers on Gait2392}',bridge+'\\subsection{Testing Alternative Controllers on Gait2392}')
 a=s.index('\\subsection{Feedback Contribution')
 add=r'''A replay of the saved human-model controller produced ten detected cycles per leg. Figure~\ref{fig:gait2392_ik_overlay} compares the joint trajectories directly with the Gait2392 inverse-kinematics reference and its associated ground-reaction forces. The model's mean stance fractions were 0.71 and 0.65 on the right and left, compared with 0.61 in the reference. Substantial hip and knee timing differences remain visible despite the detected cycles. The overlay therefore complements the scalar gait score by showing which parts of the movement disagree. Appendix~\ref{sec:app_ik_overlay} reports the angle errors and alignment method.

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.75\textheight,keepaspectratio]{gait2392_ik_overlay.pdf}
\caption[Gait2392 model motion compared with inverse kinematics]{Supported human-model joint motion compared with the Gait2392 inverse-kinematics reference. Model curves show the mean and one standard deviation across ten cycles per leg. Each leg is aligned to its own foot-loading onset; no angle offset, amplitude adjustment, or fitted phase shift is applied. Bottom strips compare model contact probability across cycles with reference stance and swing from ground-reaction forces, using a \qty{50}{\newton} threshold. Separate leg alignment does not establish bilateral phase agreement.}
\label{fig:gait2392_ik_overlay}
\end{figure}

'''
 s=s[:a]+add+s[a:]
 s=sub(s,'Forward progression and joint cycling persisted, but the feedback-disabled condition failed the contact-phase criterion','A phase of 0.5 cycles indicates alternating foot loading; the supported-stepping criterion accepts 0.35--0.65 cycles. Forward progression and joint cycling persisted, but the feedback-disabled condition failed this criterion')
 return s
edit('60-sns',sns)
def appendix(s):
 s=s.replace('failure sentinel','reserved failure value').replace('The sentinel itself','That value itself')
 s=sub(s,'The systematic pattern-formation and motoneuron stimulation series was not completed.','The systematic pattern-formation and motoneuron stimulation series for these full Gait2392 networks was not completed. The reduced-circuit tests in Section~\\ref{sec:app_layer_stimulus} address a separate, smaller model.')
 add=r'''\subsection{RG, PF, and MN Stimulus Tests}\label{sec:app_layer_stimulus}

The fourteen-neuron Deng-derived implementation follows the circuit used in SNS-Toolbox demonstrations \citep{deng_neuromechanical_2019,nourse_sns_2023}. A \qty{0.1}{\milli\second} Euler step reproduced the existing neural integrator's complete 24-second baseline voltage trace exactly in the tested process. The baseline period was \qty{1.938}{\second}. Neural output drove the existing suspended, unilateral W2L mechanical model without afferent return. Hip commands came from the MN outputs; knee and ankle commands used the existing PF-to-activation mapping. This configuration tests transmitted circuit effects, not the original rat mechanics or bilateral walking.

The initial sweep applied $+2$ and $+10$ nA to each layer for 0.1 and 0.5 seconds, together with two direct MN inhibitory controls. None met the deletion criterion. A subsequent stronger-current control series retained the same circuit and mechanical parameters (Table~\ref{tab:layer_stimulus}). RG and PF stimuli began one-quarter cycle after a baseline RG onset; MN inhibition began 0.05 seconds into a baseline motor burst. These currents are mechanistic controls, not experimentally calibrated physiological stimuli.

A deletion required at least \qty{20}{\milli\second} of baseline motor activity in the stimulus window and suppression of at least half of that activity, using an activation threshold of 0.5. Phase residuals compare the final four recovered onsets with the nearest baseline onset, normalized by the baseline period. Residuals at or below 0.05 cycles were classified as non-resetting, and those at or above 0.15 as shifted; intermediate values remained indeterminate. RG and motor residuals were assessed separately. In particular, the longer PF pulse shifted the motor pattern without resetting the RG clock.

\begin{table}[htbp]
\centering\small
\caption[Layer-specific stimulation controls]{Stronger-current controls in the reduced Deng-derived circuit. Phase residuals are fractions of the baseline cycle. ``No deletion'' means the specified suppression criterion was not met.}
\label{tab:layer_stimulus}
\begin{tabularx}{\textwidth}{@{}lrrcc>{\raggedright\arraybackslash}X@{}}
\toprule
Target & nA & s & RG phase & MN phase & Outcome \\
\midrule
RG & $+50$ & 0.5 & 0.052 & 0.052 & Deletion; phase indeterminate \\
RG & $+50$ & 1.0 & 0.307 & 0.307 & Resetting deletion \\
PF & $+50$ & 0.5 & 0.000 & 0.000 & Non-resetting deletion \\
PF & $+50$ & 1.0 & 0.000 & 0.208 & Motor pattern shifted; RG preserved \\
MN & $-100$ & 0.5 & 0.000 & 0.000 & No deletion \\
MN & $-500$ & 0.5 & 0.000 & 0.000 & Non-resetting deletion \\
\bottomrule
\end{tabularx}
\end{table}

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.74\textheight,keepaspectratio]{layer_stimulus_validation.pdf}
\caption[Layer-specific stimuli and downstream responses]{Matched baseline and stimulus responses in the reduced Deng-derived circuit coupled to a suspended model. Rows show RG stimulation, PF stimulation, and direct flexor-MN inhibition; columns show RG voltage, flexor-MN activation, and hip angle. Shading marks the applied stimulus. The RG intervention shifts subsequent timing, whereas the displayed PF and MN interventions preserve the RG schedule. Joint motion shows the downstream mechanical response, without establishing stance or walking.}
\label{fig:layer_stimulus_validation}
\end{figure}

'''
 s=sub(s,'\\subsection{Archived AnimatLab Models}',add+'\\subsection{Archived AnimatLab Models}')
 add=r'''\subsection{Direct Kinematics Overlay}\label{sec:app_ik_overlay}

The saved human-model controller reproduced its cycle-shape score of $-145.07$ in an independent process. The direct overlay in Figure~\ref{fig:gait2392_ik_overlay} uses the Gait2392 \texttt{subject01\_walk1\_ik.mot} trajectory and paired ground-reaction-force file. Model cycles begin when vertical foot loading crosses \qty{50}{\newton}; reference stance is defined using the same force threshold. Ten cycles per side, with durations between 0.35 and 2.5 seconds, were retained using the existing cycle detector. Each leg was aligned to its own onset and resampled over the gait cycle, without fitting phase, subtracting the mean angle, or scaling amplitude.

Mean cycle periods were \qty{1.285}{\second} on the right and \qty{1.254}{\second} on the left, compared with \qty{1.233}{\second} and \qty{1.213}{\second} in the reference. Model stance fractions were 0.710 and 0.649, compared with 0.610 in the reference. Table~\ref{tab:ik_absolute_rmse} reports absolute-angle errors between mean cycles. These errors retain posture offsets and are distinct from the mean-removed cycle-shape errors used in tuning. The hip and knee mismatches show why cycle counts and a scalar score must be accompanied by the trajectories and contact timing.

\begin{table}[htbp]
\centering
\caption[Absolute-angle errors against Gait2392 inverse kinematics]{Unshifted absolute-angle RMSE between mean model and reference cycles, in degrees.}
\label{tab:ik_absolute_rmse}
\begin{tabular}{lrr}
\toprule
Joint & Right & Left \\
\midrule
Hip & 28.03 & 24.65 \\
Knee & 27.44 & 28.45 \\
Ankle & 8.71 & 7.38 \\
\bottomrule
\end{tabular}
\end{table}

'''
 s=sub(s,'\\subsection{Pathway-Removal Tests}',add+'\\subsection{Pathway-Removal Tests}')
 return s
edit('95-AppendixD',appendix)
def overall(s):
 s=sub(s,'What those same evaluations localize as missing defines the work remaining here:','These evaluations identify the remaining tasks:')
 a=s.index('The verification series comes first');b=s.index('The existing integrations provide',a)
 s=s[:a]+r'''Further verification should apply specified stimuli to the rhythm-generator, pattern-formation, and motoneuron layers, while recording both neural timing and mechanical responses. The reduced-circuit tests in Section~\ref{sec:app_layer_stimulus} distinguish RG resetting from suppression of downstream motor output. Extending these tests to the full Gait2392 networks will require matched stimulus protocols and separate assessment of the RG clock and motor pattern. Physiological and modeling studies provide reference phenomena for these comparisons \citep{selionov_tonic_2009,rybak_modelling_2006,markin_afferent_2010,shevtsova_modeling_2016}.

The next experiments should stimulate, remove, or add specified pathways while retaining the other controller and mechanical settings. Existing removal tests establish a starting point (Section~\ref{sec:walker_variant_results}); adding a candidate sensory or commissural pathway would test whether it improves a defined response. Matched comparisons of spiking and non-spiking controllers should use the same tuning and evaluation conditions. The observed differences in their limb coordination motivate that comparison (Section~\ref{sec:disc_spiking}).

'''+s[b:]
 s=s.replace('\\subsection{Extending the Sensorium}','\\subsection{Adding Sensory and Higher-Level Control}').replace('Once stepping is stable, three supraspinal contributions will be added, and a fourth is in line behind them','Four planned extensions address sensory estimation, coordination, and action selection').replace('as an adaptive gateway between higher-level locomotion commands and the CPG','to adapt commands and coordination between higher-level control and the CPG').replace('a cerebellar adaptive gateway between high-level commands and the pattern generator','a cerebellar model for adapting commands and coordination')
 s=s.replace('hand-drive withdrawn entirely','externally supplied neural drive removed').replace('at weaned drive','with reduced neural drive')
 return s
edit('70-overall',overall)
def conclusion(s):
 s=s.replace('elucidated the relationship','characterized the relationship').replace('They establish an operational experimental workflow; unsupported walking and general controller transfer remain open.','These tests provide a repeatable framework for comparing neural and mechanical hypotheses. General controller transfer and independent balance remain open.')
 a=s.index('Together, these contributions connect');s=s[:a]+r'''These actuator and joint models, the sensory-afferent application, and the neuromechanical tools make design assumptions testable: force against measurement, routing against torque requirements, and sensory pathways against controlled changes in movement. The future studies in Chapter~\ref{ch:overall} extend these tests toward dynamic actuation, independent balance, and hardware. This work provides an experimental foundation for biomimetic robot design and for studying how neural control interacts with mechanical structure.
'''
 return s
edit('80-conclusion',conclusion)
(A/'integration_manifest.json').write_text(json.dumps(changes,indent=2))
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
print('Integrated figures and prose; protected sources unchanged.')
