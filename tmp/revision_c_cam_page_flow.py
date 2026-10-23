from pathlib import Path
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
p=C/'chapters/94-AppendixC.tex';s=p.read_text()
s=s.replace("Figures~\\ref{fig:app_cam_AF} and~\\ref{fig:app_cam_GH} show the current cam design across the knee range. The existing vector montage is divided between pages to retain readable labels; its paths and numerical annotations are unchanged.\n\n",'')
start=s.index(r'\section{Cam-Route Pose Sequence}');end=s.index(r'\section{Initial Cam-Route Isometric Holds}',start)
s=s[:start]+s[start:end].replace(r'\begin{figure}[p]',r'\begin{figure}[htbp]')+s[end:]
p.write_text(s,encoding='utf-8')
print('Kept cam figure section with the figure and removed layout commentary from prose')
