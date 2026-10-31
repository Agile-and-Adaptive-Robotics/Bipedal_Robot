from pathlib import Path
import re,shutil
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());D=R/'Documentation/Reports and Papers/Dissertation';p=D/'ProofFinal/chapters/95-AppendixD.tex';backup=D/'Archive/session1_sns_20261031/before_matched_ground_figure.tex'
assert not backup.exists();shutil.copy2(p,backup);s=p.read_text(encoding='utf-8')
blocks=list(re.finditer(r'\\begin\{figure\}[\s\S]*?\\end\{figure\}',s));old_neural=[m.group() for m in blocks if r'\label{fig:spike_rhythm}' in m.group()];old_ground=[m.group() for m in blocks if r'\label{fig:spike_ground}' in m.group()];assert len(old_neural)==len(old_ground)==1
new=r'''\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{spiking_ground_neuromechanics_short_20261031.pdf}
\caption[Matched neural and mechanical responses with graded and spiking controllers]{Matched \qtyrange{10}{12.5}{\second} ground-trial window: graded controller (left) and untuned spiking controller (right), at the same selected parameters and support settings. Rows show right-side RG and PF readouts, vastus-lateralis and semimembranosus activation, and bilateral hip, knee, and ankle angles. Neural rows pair extensor (orange solid) and flexor (blue dashed) outputs; joint rows pair right (orange solid) and left (blue dashed) legs. Spiking signals are filtered readouts, and muscle activation is the actuator state.}
\label{fig:spike_rhythm}\label{fig:spike_ground}
\end{figure}'''
s=s.replace(old_neural[0],'').replace(old_ground[0],new)
s=s.replace(r'respectively (Figure~\ref{fig:spike_rhythm}).','respectively.',1)
p.write_bytes(s.replace('\n','\r\n').encode('utf-8'))
print('Reconciled EB475WS4 matched-ground figure approach using local archived arrays; both existing labels preserved.')
