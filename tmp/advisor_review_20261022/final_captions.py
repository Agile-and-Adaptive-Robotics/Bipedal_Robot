from pathlib import Path
import json
H=Path(__file__).resolve().parent;D=H.parent.parent/'Documentation/Reports and Papers/Dissertation'
caps={
'60-sns.tex':r'\caption[Simulink SNS library blocks and neuron implementation]{Core SNS-Simscape blocks (top) and graded-neuron implementation (bottom). Summed currents drive membrane integration; threshold, inverse-slope gain, and saturation bound the output. Voltage feeds both current calculations. The generic BPA block remains a placeholder. These views do not demonstrate quantitative knee regulation or complete CAD coupling.}',
'70-overall.tex':r'\caption[Construction of the braided pneumatic actuator]{Proposed layered BPA finite-element model, not to scale. (A) Contact interfaces, internal pressure, and end constraints. Longitudinal symmetry permits in-plane expansion; at half length, impose axial symmetry or allow axial motion with an optional load. (B) Illustrative cutaway identifying materials and longitudinal symmetry. Middle rubber is separated only in (A); cap details and layer spacing are schematic.}'}
ps=[]
for file,cap in caps.items():
    p=D/'ProofFinal/chapters'/file;s=p.read_text(encoding='utf-8');prefix=cap.split(']{')[0]+']';old=next(l for l in s.splitlines() if l.startswith(prefix))
    ps.append({'file':file,'old':old,'new':cap});s=s.replace(old,cap,1);p.write_text(s,encoding='utf-8')
    p=H/'overleaf_expected'/file;s=p.read_text(encoding='utf-8');assert s.count(old)==1;p.write_text(s.replace(old,cap,1),encoding='utf-8')
(H/'final_caption_patches.json').write_text(json.dumps(ps))
