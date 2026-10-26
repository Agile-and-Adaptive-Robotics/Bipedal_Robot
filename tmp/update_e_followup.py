from pathlib import Path
import json,shutil,re,difflib,hashlib
R=Path(r'D:\Github\Bipedal_Robot'); E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review';F=V/'easteregg2_followup'
p=E/'chapters/60-sns.tex';s=p.read_text(encoding='utf-8')
assert not (V/'60-sns_before_evidence_followup.tex').exists()
shutil.copy2(p,V/'60-sns_before_evidence_followup.tex')
for n in ('REPORT.md','ZCODE_HANDOFF.md'):
 shutil.copy2(E/n,V/('before_followup_'+n))
shutil.copy2(V/'self_scoring_matrix.json',V/'self_scoring_matrix_before_followup.json')
mark=r'\subsection{Matched Neural Implementations}'
i=s.index(mark)
s=s[:i]+r'''\subsection{Published Neural Behaviors}\label{sec:sns_published_behaviors}

An activity-based reduction of the two-level circuit of \citet{rybak_modelling_2006} was tested against selected published neural behaviors. At the published drive settings, increasing extensor drive produced a longer extensor phase, while increasing flexor drive reversed that relationship. The deletion tests distinguished the two circuit levels: a PF-level pulse suppressed flexor MN bursts without resetting the RG schedule, whereas an RG-level pulse shifted the subsequent rhythm. The first returning flexor burst differed from the extrapolated schedule by 0.001 periods after the PF perturbation and 0.416 periods after the RG perturbation. All 13 specified checks passed, including deterministic repeat runs.

This result reproduces selected qualitative behaviors of the published neural model. The implementation uses activity-based populations in place of the original spiking formulation, calibrated global conductances, and adapted pulse amplitudes and integration settings. It therefore tests the two-level explanation of resetting and non-resetting deletions without claiming an exact numerical reproduction of the original model or any walking mechanics. The Deng-derived and Di Russo-derived tests below have different scopes and should not inherit this result.

'''+s[i:]
s=s.replace('The agreement is specific to this circuit and input protocol.', 'The agreement is specific to this circuit and input protocol. Together with the contact-removal test, it checks the contact-driven design principle; it is not a quantitative comparison with the original Li experiment.')
s=s.replace('Its period was \\qty{1.938}{\\second}, compared with approximately \\qty{1.94}{\\second} in the NumPy reference.', 'Its period was \\qty{1.938}{\\second}, compared with approximately \\qty{1.94}{\\second} in our NumPy reference, which was itself translated from the Nourse/Deng parameter tables.')
s=s.replace('It does not demonstrate pointwise voltage agreement:', 'The original Deng AnimatLab project was extracted for connectivity and parameters but was not run for this comparison. The published Deng perturbation experiment therefore remains untested in translation. This comparison also does not demonstrate pointwise voltage agreement:')
marker='The Simulink air-stepping demonstration connects'
i=s.index(marker)
s=s[:i]+r'''A separate W2L BPA-actuation test replaced all 12 muscle actuators with the Festo force model. Both knees completed at least three air-stepping excursions, with periods of approximately \qty{1.54}{\second}, and the model subsequently maintained height and foot loading with retained support. Peak pressure was \qty{372}{\kPa}, below the \qty{620}{\kPa} ceiling. This is direct evidence that the controller, BPA-force calculation, and mechanical model operated together. The retained harness supplied substantial support, and no matched Hill-muscle versus BPA comparison was performed. This W2L test is distinct from the converted 92-actuator human model used in the bridge tests.

'''+s[i:]
start=s.index('The coupled W2L evaluations address')
end=s.index('A separate controller derived',start)
s=s[:start]+r'''The coupled W2L tests ask how changing the body or selected neural pathways affects a fixed controller. Their interpretation required separating an initial parameter-loading error from a later corrected replay. The original search stored 61 parameter fields, but 55 individual pathway fields did not reach the model. Only six searched fields affected those evaluations; the remaining pathway strengths stayed at their anchors. The original mass-change and support-reduction results therefore belong to that effective controller, not to a fully optimized set of individual conductances.

After correcting the runner to apply the stored individual fields, I evaluated the resulting candidate with all settings held fixed (Figure~\ref{fig:w2l_perturbation}). This candidate had not been optimized under the corrected parameter mapping. In the \qty{16}{\second} control record, the left and right heels completed 17 and 16 contact cycles. Removing the crossed heel pathway produced 16/16 cycles, the Ia chain 13/13, the Renshaw pathways 16/14, the heel cascade 17/15, and commissural coupling 13/12. Removing Ib load feedback produced no heel-contact cycles on either side. Thus commissural removal reduced cycling but did not abolish it at this operating point. Among the tested removals, Ib feedback was required for the recorded contact cycling in this configuration. Continued cycling after other removals does not by itself establish unchanged phase, posture, or gait quality.

One-second RG stimuli and \qty{40}{\newton}, \qty{0.15}{\second} fore-aft pushes also left bilateral contact cycling present. Nevertheless, the stronger stimulus reached approximately \ang{62.8} pelvis tilt, and the forward and backward pushes reached approximately \ang{49.6} and \ang{47.7}. These runs do not satisfy a \ang{45} posture criterion merely because cycling persists. Recovery of the original gait requires a separately defined phase, posture, and timing measure.

The corrected candidate was also tested with thigh or distal-leg masses and inertias changed by $\pm10$ and $\pm20$ percent, and whole-leg properties changed by $\pm10$ percent. All ten conditions retained bilateral contact cycling, but nine exceeded the \ang{45} tilt gate. Only one condition passed the recorded survival criterion. Minimum pelvis height across these records was approximately \qty{0.898}{\metre}. This is evidence of persistent cycling with inadequate posture robustness, rather than ten successful robust-walking tests.

All these evaluations retained a virtual harness with vertical height feedback and constraints on lateral motion, roll, and yaw. The nominal weight-support setting was \qty{25}{\percent}; the realized force also depended on the feedback terms. Earlier searches reduced the nominal setting from 70 to 10 percent with re-tuning between settings, but those records precede full application of the individual fields. Their failed local conductance probes, 9/12 at nominal 25 percent and 7/12 at nominal 10 percent, also describe the earlier effective controller. They cannot certify the robustness of the corrected candidate. These are W2L results, not demonstrations of the BPA-routed Gait2392 robot model or harness-free walking.

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.56\textheight,keepaspectratio]{w2l_perturbation_fully_applied.pdf}
\caption[W2L perturbations with the stored conductances applied]{W2L records after applying the stored individual conductances, with retained harness support. (A--C) Knee motion in control and one-second RG-stimulus conditions. (D, E) Center-of-mass position and pelvis tilt after forward and backward pushes at \qty{8}{\second}. (F) Independently counted left and right heel-contact cycles after pathway removals. The control has 17/16 cycles; commissural removal retains 13/12, while Ib removal gives 0/0. Continued cycling does not establish posture recovery: both push conditions exceed \ang{45} tilt. The figure's ``fell'' field records the height-based fall detector, not passage of the tilt criterion.}
\label{fig:w2l_perturbation}
\end{figure}

'''+s[end:]
mark=r'\section{Capabilities and Remaining Limits}'
i=s.index(mark)
s=s[:i]+r'''\subsection{Feedback and Reduced Hand Drive}

The human-model feedback-tuning campaign provides a complementary test of sustained activity. With the selected feedback parameters fixed, bilateral stepping-like motion persisted as the added hand-drive currents were reduced to zero, with 8 and 11 contact cycles at zero scale. The original feedback configuration lost bilateral activity below 0.85 of those added currents. However, the descending drive remained active, the mechanics retained support, and the motion was in place. This comparison tests the effect of feedback re-tuning on dependence on added drive; it is not a fixed-controller feedback-intact versus feedback-disabled experiment. A matched feedback-ablation pair on the walking BPA model remains absent.

'''+s[i:]
s=s.replace('A brief-pulse test produced self-sustained neural bursts, short-window comparisons checked network translation, and the recorded mass-change tests retained supported progression in the W2L model.', 'The Rybak-derived neural model reproduced selected phase-resetting behaviors, a brief-pulse test produced self-sustained bursts in the Deng-derived circuit, and short-window comparisons checked network translation. The W2L tests connected pathway removals to bilateral contact cycling while exposing posture failures under perturbation.')
s=s.replace('and parameter loading and deletion scoring require correction before the corresponding optimization and pathway claims can be accepted.', 'and applying previously unused parameter fields changed the W2L operating point and the interpretation of pathway removals. Its persistent cycling under mass changes did not satisfy the posture criterion in nine of ten conditions.')
p.write_text(s,encoding='utf-8')
base=(V/'60-sns_revision_D.tex').read_text(encoding='utf-8')
(V/'changes.diff').write_text(''.join(difflib.unified_diff(base.splitlines(True),s.splitlines(True),fromfile='Revision_D/chapters/60-sns.tex',tofile='Revision_E/chapters/60-sns.tex')),encoding='utf-8')
(V/'followup_changes.diff').write_text(''.join(difflib.unified_diff((V/'60-sns_before_evidence_followup.tex').read_text(encoding='utf-8').splitlines(True),s.splitlines(True),fromfile='E_before_followup',tofile='E_updated')),encoding='utf-8')
matrix=json.loads((V/'self_scoring_matrix.json').read_text())
anchors={'E04':'A separate Simulink','E06':'The corrected candidate was also','E07':'The coupled W2L tests ask','E08':'After correcting the runner','E10':'The resulting contribution'}
for row in matrix:
 anchor=anchors.get(row['id'],row['text'][:45]); idx=next(i for i,l in enumerate(s.splitlines(),1) if anchor in l);row['line']=idx;row['text']=s.splitlines()[idx-1]
 if row['id']=='E06':row.update(change='Distinguished old six-field optimization from fully applied candidate; 9/10 mass cases fail tilt gate.',limit='Persistent contact cycling is not posture robustness; candidate has not been optimized after mapping correction.')
 if row['id']=='E07':row.update(change='Preserved original 55/61 loading failure as history; inspected newer remote runner and saved outputs.',limit='Evaluator still filters non-KV keys; verify fresh-process campaign/verification drivers before redo.')
 if row['id']=='E08':row.update(change='Independently recounted remote NPZ logs; restored corrected figure with 17/16 control, 13/12 commissural, and 0/0 Ib.',limit='Claims restricted to this operating point; push and strong-stimulus posture failures are explicit.')
 if row['id']=='E10':row.update(change='Libraries were never supplied; Ben does not want to request them. Avoided unsupported novelty claim.',limit='Unavailable comparison is a documented limitation, not an action awaiting Ben.')
for id_,anchor,concern in [('E11','An activity-based reduction','Selected Rybak published neural behaviors reproduced with documented reductions'),('E12','The human-model feedback-tuning','Human wean removes added hand drive, not all descending drive or all feedback')]:
 idx=next(i for i,l in enumerate(s.splitlines(),1) if anchor in l);matrix.append(dict(id=id_,concern=concern,section='Verification of the Implemented Interfaces',file='chapters/60-sns.tex',line=idx,text=s.splitlines()[idx-1],change=concern,limit='Scoped evidence; human/advisor acceptance pending.',draft_score=2,acceptance_score=0))
(V/'self_scoring_matrix.json').write_text(json.dumps(matrix,indent=2),encoding='utf-8')
print('Updated E with remote records; preserved earlier draft and report.')
