from pathlib import Path
import re,shutil
R=Path(r'D:\GitHub\Bipedal_Robot');O=R/'tmp/torque_review_20261027'
D=R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex'
B=D.with_name('15-background.tex');DB=D.with_name('50-database.tex')
d=D.read_text(encoding='utf-8')
m=next(m for m in re.finditer(r'\\begin\{table\}.*?\\end\{table\}',d,re.S) if r'\label{tab:mesh_penalties}' in m.group())
old=m.group();cap=re.search(r'\\caption\{(.*?)\}\n\\label',old,re.S).group(1)
rows=old[old.index(r'$J_{\text{worst}}$'):old.index(r'\bottomrule')]
head=r'\textbf{Penalty term} & \textbf{Weight} & \textbf{Meaning} \\'
new=r'''\begingroup
\small\singlespacing
\begin{longtable}{@{}p{0.28\textwidth}p{0.17\textwidth}p{0.49\textwidth}@{}}
\caption{'''+cap+r'''}\label{tab:mesh_penalties}\\
\toprule
'''+head+r'''
\midrule
\endfirsthead
\multicolumn{3}{l}{\tablename\ \thetable\ (continued)}\\
\toprule
'''+head+r'''
\midrule
\endhead
\midrule
\multicolumn{3}{r}{Continued on next page}\\
\endfoot
\bottomrule
\endlastfoot
'''+rows+r'''\end{longtable}
\endgroup'''
d=d[:m.start()]+new+d[m.end():]
D.write_text(d,encoding='utf-8',newline='\n')
b=B.read_text(encoding='utf-8').replace(r'\includegraphics[width=0.72\textwidth]{masters_leg_model.pdf}',r'\includegraphics[width=0.68\textwidth]{masters_leg_model.pdf}')
B.write_text(b,encoding='utf-8',newline='\n')
shutil.copy2(DB,O/'50-database_before.tex')
db=DB.read_text(encoding='utf-8')
db=db.replace(r'\begin{tabularx}{\textwidth}{@{}lllX@{}}',r'\begin{tabularx}{\textwidth}{@{}p{0.10\textwidth}p{0.20\textwidth}p{0.23\textwidth}L@{}}',1)
DB.write_text(db,encoding='utf-8',newline='\n')
print('Repaired the torque weight table, background figure height, and database table widths without changing their contents.')
