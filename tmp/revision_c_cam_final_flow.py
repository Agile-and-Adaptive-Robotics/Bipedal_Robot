from pathlib import Path
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
p=C/'chapters/94-AppendixC.tex';s=p.read_text();s=s.replace('\\section{Cam-Route Pose Sequence}\\label{sec:app_cam_geometry}\n\n','');p.write_text(s,encoding='utf-8')
p=C/'chapters/30-torque.tex';s=p.read_text();s=s.replace(r'Panel letters follow the full pose sequence in Appendix~\ref{sec:app_cam_geometry}.',r'Panel letters follow the full sequence in Figures~\ref{fig:app_cam_AF} and~\ref{fig:app_cam_GH}.');p.write_text(s,encoding='utf-8')
print('Removed the orphan cam subsection heading; pose figures remain in supporting route section')
