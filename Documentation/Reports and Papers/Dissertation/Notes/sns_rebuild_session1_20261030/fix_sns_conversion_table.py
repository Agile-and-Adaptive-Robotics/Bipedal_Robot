from pathlib import Path
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());p=R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/95-AppendixD.tex';s=p.read_text(encoding='utf-8');start=s.index('\\subsection{AnimatLab Conversions}');stop=s.index('\\section{Supported BPA-Robot Evaluations}',start);block=s[start:stop]
assert block.count('\\begin{tabular}{@{}lcc@{}}')==1
block=block.replace('\\begin{tabular}{@{}lcc@{}}',r'\begin{tabularx}{\textwidth}{@{}p{0.22\textwidth}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}X@{}}').replace('\\end{tabular}','\\end{tabularx}').replace('\\multirow{2}{*}','\\multirow{2}{=}')
s=s[:start]+block+s[stop:];p.write_bytes(s.replace('\n','\r\n').encode('utf-8'));print('Applied EB475WS4 table-width correction without changing any recorded value.')
