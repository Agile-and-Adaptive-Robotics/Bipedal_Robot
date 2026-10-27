from pathlib import Path
import re
R=Path(r'D:\GitHub\Bipedal_Robot')
P=R/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex'
S=P.with_name('frontiers_SupplementaryMaterial.tex')
p=P.read_text(encoding='utf-8')
# These headings are unused template placeholders after the references.
p=p.replace(r'\section*{Figures}','').replace(r'\section*{Tables}','')
P.write_text(p,encoding='utf-8',newline='\n')
s=S.read_text(encoding='utf-8')
st=s.index(r'\subsection{Configuration-Survey Coordinates and Actuator Sizes}')
en=s.index(r'\end{table}',st)+len(r'\end{table}')
inventory=s[st:en].replace(r'\begin{table}[p]',r'\begin{table}[htbp]')
s=s[:st]+s[en:]
s=s.replace(r'\end{document}',r'\clearpage'+'\n'+inventory+'\n'+r'\end{document}')
S.write_text(s,encoding='utf-8',newline='\n')
