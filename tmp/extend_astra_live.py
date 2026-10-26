from pathlib import Path
import json
base=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\LIVE_ASTRA_20261026')
patches=json.loads((base/'patches.json').read_text(encoding='utf-8'))
bg=(base/'before/15-background.tex').read_text(encoding='utf-8')
old=next(p for p in bg.split('\n\n') if p.startswith('The third line, in which'))
new=r'''The third line, in which this dissertation sits, uses synthetic nervous systems whose neurons and synapses explicitly represent candidate physiological circuits. SNS-Toolbox supports conductance-based networks at speeds suitable for selected real-time robotic applications \citep{nourse_sns_2023}. Related work includes motor-module bipedal walkers \citep{ichimura_pathological_2019}, two-layer rhythm-generator and pattern-formation control of a musculoskeletal macaque model \citep{nishizaki_gait_2026}, and CPG control of muscle-driven robotic hind legs \citep{fukuoka_autonomous_2021}. These studies establish the value of coupled neural and mechanical models. For the humanoid design studied here, the remaining methodological question is how to implement and verify candidate spinal circuits with human musculotendon mechanics and custom BPA actuation, while distinguishing neural effects from conversion or coupling errors. Chapter~\ref{ch:sns} addresses that question through specific model implementations, graphical neural blocks, and controlled verification and feedback tests.'''
patches.append(dict(file='15-background.tex',old=old,new=new,reason='Replace unsupported field-wide novelty claim with the scoped methodological gap; retain cited prior work.'))
(base/'after/15-background.tex').write_text(bg.replace(old,new),encoding='utf-8')
app=(base/'after/95-AppendixD_LIVE.tex').read_text(encoding='utf-8')
old=r'\subsection{Archived AnimatLab Models}\label{sec:animatlab_validation}'
new=r'''\subsection{A Published Neural Mechanism as a Verification Target}\label{sec:rybak_verification}

The reduced Python implementation of the two-level circuit of \citet{rybak_modelling_2006} passed 13 checks covering drive-dependent activity, repeatability, and the distinction between resetting and non-resetting deletions. Stimulating the pattern-formation extensor population suppressed flexor motor output while the rhythm generator continued, so flexor activity resumed at its expected phase. Stimulating the rhythm-generator extensor population also suppressed flexor output but shifted the subsequent rhythm phase. This comparison tests the functional separation between rhythm generation and pattern formation: similar missing motor bursts can arise through different neural mechanisms. It is a neural-only verification of the translated circuit, independent of a mechanical model and separate from the Simulink demonstrations.

\subsection{Archived AnimatLab Models}\label{sec:animatlab_validation}'''
assert app.count(old)==1
patches.append(dict(file='95-AppendixD.tex',old=old,new=new,reason='Substantiate the Rybak table entry with the verified resetting/non-resetting mechanism; identify Python and neural-only scope.'))
(base/'after/95-AppendixD.tex').write_text(app.replace(old,new),encoding='utf-8')
(base/'patches.json').write_text(json.dumps(patches,ensure_ascii=False,indent=2),encoding='utf-8')
print(len(patches))
