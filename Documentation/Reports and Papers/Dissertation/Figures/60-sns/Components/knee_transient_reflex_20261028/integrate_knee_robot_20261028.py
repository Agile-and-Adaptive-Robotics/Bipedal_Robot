from pathlib import Path
import shutil,json,hashlib
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitz
R=Path('D:/Github/Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';A=D/'Archive/knee_and_robot_update_20261028';K=A/'knee_transient_reflex_20261028';F=D/'Figures';P=D/'ProofFinal/chapters'
before=A/'sources_before';before.mkdir(exist_ok=True)
for n in ['60-sns','70-overall','95-AppendixD']:shutil.copy2(P/(n+'.tex'),before/(n+'.tex'))
protected={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'20-muscles.tex',P/'30-torque.tex',P/'94-AppendixC.tex',R/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex']};(A/'protected_hashes.json').write_text(json.dumps(protected,indent=2))
m=loadmat(K/'replay_summary.mat',squeeze_me=True);t=m['t'];C=m['C'];Q=m['P'];fine=m['F'];d=-(Q[:,0]-C[:,0])*180/np.pi;sel=t>=30;idx=np.flatnonzero(sel)[np.argmax(d[sel])];baseline=np.mean(C[(t>=25)&(t<30),0]);bad=np.flatnonzero((t>=30.03)&((abs((Q[:,0]-baseline)*180/np.pi)>.01)|(abs(Q[:,1])>1e-4)))
metrics={'paper_extension_deg':float(d[idx]),'peak_latency_s':float(t[idx]-30),'recovery_s':float(t[bad[-1]+1]-30),'step_error_deg':float(np.max(abs(Q[:,0]-fine[:,0]))*180/np.pi),'pre_pulse_control_equal':bool(np.array_equal(C[t<30],Q[t<30]))};assert metrics['pre_pulse_control_equal'];assert abs(metrics['paper_extension_deg']-1.090792830529)<1e-9;assert abs(metrics['recovery_s']-6.355)<1e-9;(A/'independent_knee_array_check.json').write_text(json.dumps(metrics,indent=2))
comp=F/'60-sns/Components/knee_transient_reflex_20261028';comp.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':13,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(2,1,figsize=(6.53,3.9),sharex=True,gridspec_kw={'height_ratios':[.7,2.3]});ix=(t>=29)&(t<=40);x=t[ix]-30
ax[0].plot(x,Q[ix,8],color='#0072B2',lw=1.5);ax[0].set_ylabel('Ia input, nA');ax[0].set_yticks([0,2]);ax[0].set_ylim(-.1,2.3)
ax[1].plot(x,-Q[ix,0]*180/np.pi,color='#0072B2',lw=1.6,label='Stimulus');ax[1].plot(x,-C[ix,0]*180/np.pi,'--',color='black',lw=1.4,label='No stimulus');ax[1].set_ylabel('Knee angle, deg');ax[1].set_xlabel('Time from pulse onset, s');ax[1].legend(loc='upper right',fontsize=11);ax[1].text(.98,.42,'Extension positive',transform=ax[1].transAxes,ha='right',fontsize=11)
for a in ax:a.set_xlim(-1,10);a.grid(alpha=.2)
fig.tight_layout();fig.savefig(comp/'transient_panel.pdf');plt.close(fig)
old=F/'60-sns/knee_reflex_current_layout.pdf';shutil.copy2(old,A/'previous_knee_layout.pdf');src=fitz.open(old);panel=fitz.open(comp/'transient_panel.pdf');out=fitz.open();page=out.new_page(width=470,height=740);page.show_pdf_page(fitz.Rect(0,20,470,300),panel,0);page.insert_text((8,14),'A  Transient sensory-current response',fontsize=12);page.show_pdf_page(fitz.Rect(0,305,470,740),src,0,clip=fitz.Rect(0,305,470,740));out.save(comp/'knee_transient_current_layout.pdf');out[0].get_pixmap(matrix=fitz.Matrix(1.4,1.4)).save(A/'knee_transient_current_layout.png');out.close();src.close();shutil.copy2(comp/'knee_transient_current_layout.pdf',F/'60-sns/knee_transient_current_layout.pdf')
# The previous startup layout is superseded, but preserved for audit.
dep=F/'60-sns/Deprecated';dep.mkdir(exist_ok=True);old.rename(dep/'knee_reflex_current_layout.pdf')
shutil.copy2(__file__,comp/Path(__file__).name);shutil.copy2(K/'run_transient.m',comp/'run_transient.m');shutil.copy2(K/'replay_summary.json',comp/'replay_summary.json');(comp/'alt-text.txt').write_text('A 2 nA, 30 ms extensor Ia current pulse evokes positive knee extension of about 1.09 degrees and recovery to the matched control. Other panels retain the antagonist-BPA model and its six sensory connections. A modified graded-response configuration is used.')
ac=F/'95-AppendixD/Components/knee_transient_reflex_20261028';ac.mkdir(exist_ok=True);shutil.copy2(K/'replay_response.pdf',ac/'replay_response.pdf');shutil.copy2(K/'replay_response.pdf',F/'95-AppendixD/knee_transient_response.pdf')
asset=A/'AH001_20261028';prov=json.loads((asset/'asset_provenance.json').read_text());img=asset/'ah001_abductor_control_still.png';assert hashlib.sha256(img.read_bytes()).hexdigest()==prov['outputs'][img.name]
shutil.copy2(img,F/'60-sns/ah001_supported_control.png');rc=F/'60-sns/Components/ah001_supported_control';rc.mkdir(exist_ok=True)
for n in ['asset_provenance.json','PLACEMENT_RECOMMENDATION.md','ah001_abductor_control_still.alt.txt']:shutil.copy2(asset/n,rc/n)
def edit(n,fn):
 p=P/(n+'.tex');s=p.read_text(encoding='utf8');p.write_text(fn(s),encoding='utf8')
def sns(s):
 a=s.index('The library assembles Ia stretch');b=s.index('A second demonstration tests feedback',a)
 s=s[:a]+r'''The library assembles Ia stretch excitation, reciprocal inhibition, and Ib autogenic inhibition in a reduced-order antagonist knee model (Figure~\ref{fig:sns_reflex_current}C). After the model reaches a steady operating point, a \qty{2}{\nano\ampere}, \qty{30}{\milli\second} current pulse at the extensor Ia sensory-neuron input evokes \ang{1.09} of extension, peaking \qty{111}{\milli\second} after pulse onset. The knee then returns toward the matched no-stimulus trajectory, meeting the specified recovery thresholds after \qty{6.355}{\second}. The displayed configuration broadens the MN and Ia synaptic saturation ramps; the original narrow-ramp settings maintained oscillations. This is a transient test of the implemented feedback and actuator equations, rather than a physiologically calibrated tendon-tap model. Appendix~\ref{sec:app_knee_transient} gives the configuration and numerical checks.

'''+s[b:]
 s=s.replace('{knee_reflex_current_layout.pdf}','{knee_transient_current_layout.pdf}').replace('Knee-model implementation comparison: startup kinematics under constant load and sustained drive (A), antagonist-BPA leg model (B), and all six implemented feedback connections (C).','Transient knee response to a 2 nA, 30 ms sensory-current pulse after settling (A), antagonist-BPA leg model (B), and the six implemented feedback connections (C). Flexion is negative. The graded-response configuration uses broadened MN and Ia synaptic saturation ramps; its response is compared with a matched no-stimulus control.')
 anchor='No harness-free walking was demonstrated.'
 s=s.replace(anchor,anchor+r'''

A newer CAD-model comparison tests whether additional abductor control prevents the legs from crossing (Figure~\ref{fig:ah001_supported_control}). Under the same assistance settings, muscle-mediated separation feedback changed reversed ankle ordering from \qty{67}{\percent} of sampled states to zero. The model uses a layered SNS for timing together with an engineering force allocator. Its MN target-current compensation cancels incoming PF and reflex currents, so this comparison demonstrates supported model integration and a separation-control effect, rather than full neural pattern control. Appendix~\ref{sec:app_cad_supported_stance} gives the matched measures and foot-mount sensitivity.

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.55\textheight,keepaspectratio]{ah001_supported_control.png}
\caption[CAD-derived robot under supported hybrid control]{Anterior (A) and oblique (B) views of the same recorded CAD-model state. Right and left legs are blue and gold; selected hip-abducting routes are green. The 75 kg model uses 54 reference-tension actuators with hybrid neural and postural control. Vertical assistance supplies 80\% of body weight in this interval; lateral and angular assistance remain active. Attachments, the welded foot mount, and passive knee stops are provisional.}
\label{fig:ah001_supported_control}
\end{figure}
''')
 return s
edit('60-sns',sns)
def app(s):
 marker=r'\subsection{A Published Neural Mechanism as a Verification Target}'
 add=r'''\subsection{Transient Knee Sensory-Current Response}\label{sec:app_knee_transient}

The knee block model was tested from a settled operating point using a rectangular current pulse added after the spindle transducer and before the extensor Ia sensory-neuron membrane input. Stimulus and no-stimulus cases start from the same initial conditions and run continuously for \qty{50}{\second}; the \qty{2}{\nano\ampere} pulse acts from 30.000 to \qty{30.030}{\second}. The original configuration did not settle: its knee angle spanned \ang{0.712} during the 25--30 s baseline window. Its former startup comparison is therefore not used to demonstrate transient reflex recovery.

The graded-response comparison changes exactly six parameters: the two MN output ramp widths from 1 to \qty{10}{\milli\volt}, and the four Ia synaptic ramp widths from 0.5 to \qty{10}{\milli\volt}. Equations, thresholds, conductances, wiring, mechanics, tonic drives, and activation dynamics are retained. This reduces loop gain as an explicit engineering comparison, rather than fitting physiological parameters. The simple actuator-force law is distinct from the characterized BPA model of Chapter~\ref{ch:muscles}. The flexor activation is zero at baseline, so the experiment does not establish suppression of a tonically active antagonist.

The modified configuration passes baseline settling tests for angle, velocity, all six membrane states, and both activations. Its paper-coordinate baseline is \ang{-42.6392}; positive excursions indicate extension. Extension peaks at \ang{1.0908} after \qty{111}{\milli\second}. Recovery requires displacement within \ang{0.01} of baseline and absolute angular velocity below \qty{1e-4}{\radian\per\second}; both remain within those bounds from \qty{6.355}{\second} after stimulus onset. The late residual is approximately \ang{9.3e-7}. A fresh-process replay reproduces the resampled pulse arrays exactly. Halving the ODE45 maximum step from 1 to \qty{0.5}{\milli\second}, at unchanged tolerances, changes the angle trace by at most \ang{5.21e-5}. Comparisons use linear interpolation onto a 1 ms grid, without smoothing. These checks establish numerical recovery in the model, not physiological reflex fidelity.

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.76\textheight,keepaspectratio]{knee_transient_response.pdf}
\caption[Transient knee sensory-current test]{Stimulus, knee angle, angular velocity, MN voltages, activations, and actuator forces in the graded-response knee model. The pulse enters the extensor Ia sensory neuron after settling. Knee angle is compared with the no-stimulus control; extensor and flexor signals are distinguished in the other paired panels. The display applies only the sign conversion to flexion-negative coordinates.}
\label{fig:app_knee_transient}
\end{figure}

'''
 s=s.replace(marker,add+marker)
 a=s.index('The model completed a \\qty{20}{\\second} supported crouched-stance test.');b=s.index('The W2L replay metrics',a)
 s=s[:b]+r'''A subsequent matched 40-second comparison retained the same geometry, initial pose, reference-tension capacities, and external assistance in both conditions. During 15--40 s, vertical assistance supplied 80\% of body weight. The additional control combines muscle-mediated separation feedback (160 N/m, capped at 20 N) with a 6\% activation floor on currently abducting routes, compared with the 5\% global floor (Table~\ref{tab:cad_abductor}). Full recorded trajectories reproduce bit for bit between fresh reference and accelerated runs. This is numerical replay within the specified model, rather than a robustness test.

\begin{table}[htbp]
\centering\small
\caption[Matched CAD-model separation-control comparison]{Matched CAD-model comparison over 15--40 s. Vertical assistance remains 80\% of body weight in both conditions.}
\label{tab:cad_abductor}
\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}Xrr@{}}
\toprule
Measure & Control off & Control on \\
\midrule
Reversed ankle ordering (\%) & 67 & 0 \\
Minimum ankle separation (mm) & $-38.3$ & 142.4 \\
Median ankle separation (mm) & $-6.2$ & 159.0 \\
Mean right TFL activation (\%) & 22.3 & 81.8 \\
Mean left TFL activation (\%) & 99.6 & 99.8 \\
Maximum sampled knee-loop error (mm) & 4.485 & 0.387 \\
Maximum sampled interleg penetration (mm) & 2.700 & 0 \\
\bottomrule
\end{tabularx}
\end{table}

The 262-neuron controller receives 230 inputs and drives 54 reference-tension actuators. Its MN target-current compensation cancels incoming PF and reflex synaptic currents to track the force allocator's targets. Consequently, the effective comparison concerns RG timing, engineered force allocation, and mechanical response, rather than normal PF/reflex authority over activation. Toe-out and small post-ramp knee excursions remain, and the stepping acceptance tests fail.

The foot attachment is an unverified assembly datum. A separate hypothesis translates each whole foot 60 mm toward the proximal shank, retaining tibial geometry and mass. With otherwise unchanged hybrid control, median ankle separation changes from 159.0 to \qty{57.3}{\milli\meter}, and sampled interleg penetration reaches \qty{1.316}{\milli\meter}. A separate W2L-organization timing probe has 544 neurons and 230 inputs; four named ankle indices alias knee-PF cells and are not additional neurons. Its shifted-foot condition passes bilateral unloading and lifting counts but fails repeated alternating foot apices and forward progression. The inherited MN compensation remains active, so this is not a full PF-driven controller comparison.

Route sites can lie outside CAD surfaces and transmit force through their named parent bodies, equivalent to ideal rigid offset brackets or guides. Bracket mass, compliance, packaging, and collision-aware wrapping remain unmodeled. These attachment assumptions are separate from verification of the source moving-point functions and complete path derivatives. The older 20-second stance test and newer 80\%-assistance comparison are different configurations, not a matched support-reduction experiment.

'''+s[b:]
 return s
edit('95-AppendixD',app)
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
print(json.dumps(metrics,indent=2))
