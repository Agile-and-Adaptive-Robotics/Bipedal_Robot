from pathlib import Path
import json,re
H=Path(__file__).resolve().parent;D=H.parent.parent/'Documentation/Reports and Papers/Dissertation'
patches=[]
for file in ['10-introduction.tex','60-sns.tex','70-overall.tex']:
    old=(H/'overleaf_before/chapters'/file).read_text(encoding='utf-8');new=(D/'ProofFinal/chapters'/file).read_text(encoding='utf-8')
    pairs=[]
    starts={'10-introduction.tex':[r'\caption[Dissertation roadmap]'],'60-sns.tex':[r'\caption[Simulink SNS library blocks',r'\caption[Five-stage training'],'70-overall.tex':[r'\caption[Construction of the braided']}
    for prefix in starts[file]:pairs.append((next(x for x in old.splitlines() if x.startswith(prefix)),next(x for x in new.splitlines() if x.startswith(prefix))))
    if file=='60-sns.tex':
        pairs.extend([(r'\includegraphics[width=\textwidth,height=0.21\textheight,keepaspectratio]{simulink_neuron_detail.pdf}',r'\includegraphics[width=\textwidth]{simulink_neuron_detail.pdf}'),('which also shows the internal implementation of the non-spiking neuron block as a native Simulink export (bottom)','which also shows a functional schematic of the non-spiking neuron implementation (bottom)')])
    if file=='70-overall.tex':pairs.append(('bpa_fea_exploded_half_cutaway_v2.png','bpa_exploded_cutaway.pdf'))
    for a,b in pairs:
        assert old.count(a)==1,(file,a)
        patches.append({'file':file,'old':a,'new':b});old=old.replace(a,b,1)
    (H/'overleaf_expected').mkdir(exist_ok=True)
    (H/'overleaf_expected'/file).write_text(old,encoding='utf-8')
(H/'figure_patches.json').write_text(json.dumps(patches),encoding='utf-8')
print(len(patches),'targeted replacements')
