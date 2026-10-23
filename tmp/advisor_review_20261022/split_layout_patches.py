from pathlib import Path
import json
H=Path(__file__).resolve().parent
ps=json.loads((H/'layout_patches.json').read_text());out=[]
for e in ps:
    if '\n' not in e['old']:out.append(e);continue
    for a,b in zip(e['old'].splitlines(),e['new'].splitlines()):
        if a==r'\par\medskip':
            for p in [H/'overleaf_expected/60-sns.tex',H.parent.parent/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex']:
                s=p.read_text(encoding='utf-8');s=s.replace(r'\par\smallskip'+'\n'+r'\includegraphics[width=0.94\textwidth]{simulink_neuron_detail.pdf}',r'\par\medskip'+'\n'+r'\includegraphics[width=0.94\textwidth]{simulink_neuron_detail.pdf}');p.write_text(s,encoding='utf-8')
            continue
        out.append(dict(e,old=a,new=b))
(H/'layout_patches.json').write_text(json.dumps(out))
