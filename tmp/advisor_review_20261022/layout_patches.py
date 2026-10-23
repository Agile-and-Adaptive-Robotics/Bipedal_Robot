from pathlib import Path
import json
H=Path(__file__).resolve().parent;D=H.parent.parent/'Documentation/Reports and Papers/Dissertation'
caps={
'60-sns.tex':r'\caption[Simulink SNS library blocks and neuron implementation]{Core SNS-Simscape blocks (top) and graded-neuron signal flow (bottom). Synaptic, leak, and applied currents drive membrane integration; threshold, inverse-slope gain, and saturation produce bounded output. Voltage feeds both current calculations. Other library blocks implement synaptic summation, diameter-specific BPA force, and biological muscle. The generic BPA block remains a placeholder; these views do not demonstrate quantitative knee regulation or completed CAD coupling.}',
'70-overall.tex':r'\caption[Construction of the braided pneumatic actuator]{Proposed layered BPA finite-element model, not to scale. (A) Wall layers, contact interfaces, internal pressure, and end constraints. Longitudinal symmetry permits in-plane expansion. At half length, impose axial symmetry, or replace it with axial motion and an optional load. (B) Illustrative exploded cutaway identifying the cap, wall materials, and longitudinal symmetry plane. Middle rubber is separated only in (A); cap details and layer spacing are schematic.}'}
ps=[]
for file,cap in caps.items():
    p=D/'ProofFinal/chapters'/file;s=p.read_text(encoding='utf-8');prefix=cap.split(']{')[0]+']'
    old=next(l for l in s.splitlines() if l.startswith(prefix));pairs=[(old,cap)]
    if file=='60-sns.tex':
        pairs.extend([(r'\includegraphics[width=0.96\textwidth,height=0.42\textheight,keepaspectratio]{sns_library_overview.pdf}',r'\includegraphics[width=0.94\textwidth,height=0.42\textheight,keepaspectratio]{sns_library_overview.pdf}'),(r'\par\medskip'+'\n'+r'\includegraphics[width=\textwidth]{simulink_neuron_detail.pdf}',r'\par\smallskip'+'\n'+r'\includegraphics[width=0.94\textwidth]{simulink_neuron_detail.pdf}')])
    else:
        pairs.append((r'\includegraphics[width=\textwidth]{bpa_construction.pdf}\par\medskip'+'\n'+r'\includegraphics[width=\textwidth]{bpa_exploded_cutaway.pdf}',r'\includegraphics[width=0.94\textwidth]{bpa_construction.pdf}\par\smallskip'+'\n'+r'\includegraphics[width=0.94\textwidth]{bpa_exploded_cutaway.pdf}'))
    for a,b in pairs:
        assert s.count(a)==1;s=s.replace(a,b,1);ps.append({'file':file,'old':a,'new':b})
    p.write_text(s,encoding='utf-8',newline='')
    target=H/'overleaf_expected'/file;t=target.read_text(encoding='utf-8')
    for e in ps:
        if e['file']==file:assert t.count(e['old'])==1;t=t.replace(e['old'],e['new'],1)
    target.write_text(t,encoding='utf-8')
(H/'layout_patches.json').write_text(json.dumps(ps),encoding='utf-8')
print(len(ps),'caption and sizing patches')
