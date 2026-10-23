from pathlib import Path
import re,fitz
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
d=fitz.open(C/'figs/30-torque/extensor_cam_route.pdf');rect=d[0].rect;print('CAM RECT',rect)
w,h=rect.width,rect.height
# Crop the existing vector montage at its tile boundaries; no plot values or pixels are edited.
def graphic(left,bottom,right,top,width=r'\linewidth'):
    return r'\includegraphics[width='+width+',trim='+f'{left:.3f}bp {bottom:.3f}bp {right:.3f}bp {top:.3f}bp'+r',clip]{extensor_cam_route.pdf}'
parts=[]
for left,bottom,right,top in [(0,2*h/3,2*w/3,0),(w/3,0,w/3,2*h/3),(2*w/3,0,0,2*h/3)]:
    parts.append(r'\begin{minipage}[t]{0.32\textwidth}\centering'+'\n'+graphic(left,bottom,right,top)+'\n'+r'\end{minipage}')
p=C/'chapters/30-torque.tex';s=p.read_text()
for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',s,re.S):
    if r'\label{fig:ExtensorRoute_redesigned}' not in m.group():continue
    b=r'''\begin{figure}[p]
\centering
'''+ '\hfill\n'.join(parts)+r'''
\caption[Cam-extensor route geometry]{Cam-extensor route at \ang{-120} flexion (A) and \ang{10} extension (H), with the original clearance and contact legend. Squares mark the femoral origin, diamonds the tibial insertion, and blue lines the actuator route. Panel letters follow the full pose sequence in Appendix~\ref{sec:app_cam_geometry}.}
\label{fig:ExtensorRoute_redesigned}
\end{figure}'''
    s=s[:m.start()]+b+s[m.end():];break
p.write_text(s,encoding='utf-8')
p=C/'chapters/94-AppendixC.tex';s=p.read_text()
insert=r'''\section{Cam-Route Pose Sequence}\label{sec:app_cam_geometry}

Figures~\ref{fig:app_cam_AF} and~\ref{fig:app_cam_GH} show the current cam design across the knee range. The existing vector montage is divided between pages to retain readable labels; its paths and numerical annotations are unchanged.

\begin{figure}[p]
\centering
'''+graphic(0,h/3,0,0,r'\textwidth')+r'''
\caption{Cam-extensor poses A--F, spanning \ang{-120} to \ang{-27}. Blue lines show the route; colored outlines show the clearance envelopes. The remaining poses and original legend are in Figure~\ref{fig:app_cam_GH}.}
\label{fig:app_cam_AF}
\end{figure}

\begin{figure}[p]
\centering
'''+graphic(0,0,0,2*h/3,r'\textwidth')+r'''
\caption{Cam-extensor poses G and H at \ang{-10} and \ang{10}, with the original legend, \qty{96.9}{\mm} cam-drum radius, and contact-release schedule. Panel letters continue Figure~\ref{fig:app_cam_AF}.}
\label{fig:app_cam_GH}
\end{figure}
\clearpage

'''
idx=s.index(r'\section{Initial Cam-Route Isometric Holds}');s=s[:idx]+insert+s[idx:];p.write_text(s,encoding='utf-8')
p=C/'review/04_Torque.md';s=p.read_text();s+='The portrait cam montage was too tall at full text width. The main chapter now shows its full-flexion and full-extension tiles plus the original legend. The complete unchanged vector montage is framed as two readable groups in Appendix C, with original panel letters retained.\n';p.write_text(s,encoding='utf-8')
print('Reframed existing vector tiles for readable main and appendix figures')
