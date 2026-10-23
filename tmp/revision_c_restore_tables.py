from pathlib import Path
import re, shutil, difflib
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
S=C/'review/table_restoration_20261023';S.mkdir(exist_ok=True)
p=C/'chapters/30-torque.tex';q=C/'chapters/94-AppendixC.tex'
s=p.read_text();a=q.read_text();old_s=s;old_a=a
for f in [p,q]:
    assert not (S/f.name).exists(), 'Already applied'
    shutil.copy2(f,S/f.name)
def get_table(t,key):
    return next(m.group() for m in re.finditer(r'\\begin\{(table|sidewaystable)\}.*?\\end\{\1\}',t,re.S) if '\\label{'+key+'}' in m.group())
blocks={key:get_table(a,key) for key in ['tab:FlxPinPreGof','tab:FlxPinGof','tab:ExtPinGof']}
for key,b in blocks.items(): a=a.replace(b,'')
# Restore tables next to their corresponding result figures.
for key,fig in [('tab:FlxPinPreGof','fig:FlxPin_preopt'),('tab:FlxPinGof','fig:FlxPin_optimized')]:
    pos=s.index(r'\end{figure}',s.index(r'\label{'+fig+'}'))+len(r'\end{figure}')
    s=s[:pos]+'\n\n'+blocks[key]+s[pos:]
summary=get_table(s,'tab:torque_fit_summary')
s=s.replace(summary,blocks['tab:ExtPinGof'])
s=s.replace('Table~\\ref{tab:torque_fit_summary} summarizes the training and validation means at the adopted terms. The full per-test comparisons, including maximum residuals and the earlier force-model comparison, are retained in Appendix~\\ref{sec:app_torque_fit_details}.','Table~\\ref{tab:ExtPinGof} gives the per-test extensor fit, including training and validation means. The earlier cross-configuration fit is retained in Appendix~\\ref{sec:app_torque_fit_details}.')
xi=r'''\begin{table}[htbp]
\centering\small
\caption{Adopted transmission corrections. The stiffness pair is shared. The wrapping factor was identified on the extensor and is also used in the flexor redesign; it was not fitted to the pinned-flexor tests.}
\label{tab:AdoptedXi}
\begin{tabular}{@{}lrrrr@{}}
\toprule
Configuration & $\chi_{0}$ (mm) & $\chi_{1}$ (N/m) & $\chi_{2}$ (N/m) & $\chi_{3}$\\
\midrule
Flexor & $+3.94$ & $3.998\times10^{4}$ & $1.473\times10^{4}$ & $0.158$\\
Extensor & $-6.38$ & $3.998\times10^{4}$ & $1.473\times10^{4}$ & $0.158$\\
\bottomrule
\end{tabular}
\end{table}'''
s=s.replace(get_table(s,'tab:AdoptedXi'),xi)
a=a.replace('These tables retain the per-test errors, fold assignments, and maximum residuals underlying Chapter~\\ref{ch:torque}. The first compares three actuator-force predictions before transmission correction. The subsequent tables give the adopted corrected fits and the earlier cross-configuration results; these are distinct identification stages.','The earlier cross-configuration fit below is distinct from the adopted two-stage identification. The pre-optimization comparison and the adopted pinned-flexor and pinned-extensor per-test fits appear in Chapter~\\ref{ch:torque} (Tables~\\ref{tab:FlxPinPreGof}, \\ref{tab:FlxPinGof}, and~\\ref{tab:ExtPinGof}).')
for key,b in blocks.items():
    assert get_table(s,key)==b and ('\\label{'+key+'}') not in a
p.write_text(s,encoding='utf-8');q.write_text(a,encoding='utf-8')
(S/'changes.diff').write_text(''.join(difflib.unified_diff(old_s.splitlines(True),s.splitlines(True),fromfile='before/30-torque.tex',tofile='chapters/30-torque.tex'))+''.join(difflib.unified_diff(old_a.splitlines(True),a.splitlines(True),fromfile='before/94-AppendixC.tex',tofile='chapters/94-AppendixC.tex')),encoding='utf-8')
(S/'README.md').write_text('Restored three complete per-test GoF tables to Chapter 4, byte-for-byte, including all data and scoring conditions. Removed redundant mean-only summary (preserved in before snapshot). Reoriented adopted corrections into flexor/extensor rows and chi0–chi3 columns. chi3 is transferred from extensor to flexor redesign, not fitted on pinned-flexor data. Earlier cross-configuration table remains in Appendix C. User confirmed four columns chi0–chi3. No result numbers recalculated.\n',encoding='utf-8')
print('Restored three GoF tables and two-row correction table')
