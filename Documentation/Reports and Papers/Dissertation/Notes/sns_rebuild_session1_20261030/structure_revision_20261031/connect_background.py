from pathlib import Path
import re
p=Path(r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\structure_revision_20261031')
c=(p/'60-sns.tex').read_text(encoding='utf-8'); a=(p/'95-AppendixD.tex').read_text(encoding='utf-8')
refs=r'''The circuit diagrams in Chapter~\ref{ch:background} provide the reference architectures for these tests: the Rybak RG/PF hierarchy (Figure~\ref{fig:rybak_cpg_circuit}), the two-joint PF populations in the Shevtsova model (Figure~\ref{fig:shevtsova_spinal_circuit}), and the Deng sensory--motor organization (Figure~\ref{fig:deng_connectome}). The Di Russo reflex network and Shinohara bilateral network (Figures~\ref{fig:afferent_connectomes} and~\ref{fig:shinohara_connectome}) provide further comparisons for afferent routing and interleg coordination. The experiments here test selected implemented mechanisms from this literature.

'''
c=c.replace('Biarticular mechanics require explicit consideration when organizing the pattern-formation layer.',refs+'Biarticular mechanics require explicit consideration when organizing the pattern-formation layer.',1)
c=c.replace('The training workflow combines the four neural-training operations described by \\citet{hunt_development_2017} with my five-stage physical curriculum','The training workflow combines the four neural-training operations presented by \\citet{hunt_inrc_2019} with my five-stage physical curriculum',1)
for name,text in [('60-sns.tex',c),('95-AppendixD.tex',a)]:
 text=re.sub(r'\bpilot\b','simulation',text,flags=re.I)
 (p/name).write_text(text,encoding='utf-8')
print('Background diagram references added; pilot terminology removed.')
print('Motivation words (including captions):',len(c[c.index('\\section{Motivation'):c.index('\\section{Knowledge')].split()))
print('Gap words:',len(c[c.index('\\section{Knowledge'):c.index('\\begin{sidewaysfigure}')].split()))
