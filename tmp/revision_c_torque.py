from pathlib import Path
import re, shutil, hashlib, json
import fitz
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
p=C/'chapters/30-torque.tex';s=p.read_text(encoding='utf-8-sig')
app=C/'chapters/94-AppendixC.tex';a=app.read_text(encoding='utf-8-sig')
(C/'review/before/94-AppendixC.tex').write_bytes(app.read_bytes())
a=a.replace(r'{figs/94-AppendixC/}',r'{figs/94-AppendixC/}{figs/30-torque/}')
def para_start(st):
    global s
    start=s.index(st);end=s.find('\n\n',start)
    if end<0:end=len(s)
    return s[start:end]
archive=[]
# Keep optimization diagnostics available, but out of the main argument.
for st in ['The Bi-pulley batch applied', 'The robot-body model itself carries', 'The same redesign pipeline was then extended']:
    old=para_start(st);archive.append(old);s=s.replace(old,'')
pat=r'\\begin\{table\}\[htbp\].*?\\end\{table\}'
for m in list(re.finditer(pat,s,re.S)):
    if r'\label{tab:bipulley_batch}' in m.group():
        archive.insert(1,m.group());s=s.replace(m.group(),'');break
old=para_start('The route redesign of Section~\\ref{sec:ongoing} placed')
archive.insert(0,old)
s=s.replace(old,r'''At unity transmission, the knee-flexor pulley route increased peak predicted knee-axis torque from \qty{81.23}{\newton\meter} to \qty{85.74}{\newton\meter}, while limiting off-axis torque to \qty{0.803}{\newton\meter}. This comparison measures the effect of the two optimized routing configurations; their different objective functions do not isolate the pulley alone. The pulley serves as a redirect in this configuration rather than a force multiplier. Appendix~\ref{sec:app_routing_extensions} records the optimizer diagnostics and the broader route searches.

A subsequent transmission study re-solved actuator sizes for each tackle configuration at five anatomical muscle sites. A two-BPA pair through a two-part tackle preserves net actuator force while doubling insertion travel; increasing travel can extend the angle range over which a route meets its torque requirement. The corresponding torque and margin surfaces are shown in Figures~\ref{fig:PulleySurfaces} and~\ref{fig:PulleyMargins_resolved}. These calculations use the moving-force-point routing, with the adjacent joint fixed for the biarticular evaluation planes; they do not establish simultaneous performance across both joints. Physical verification of the pulley configurations remains future work.''')
start=s.index(r'\begin{figure}[p]',s.index('The Mono/Bi-pulley algorithms extend'))
end=s.index(r'\end{figure}',start)+len(r'\end{figure}')
s=s[:start]+s[end:]
old=para_start('Physical testing of the redesigned extensor route')
s=s.replace(old,r'''Nine isometric holds have been recorded on the rebuilt cam extensor route with a $\diameter$\qty{20}{\mm} BPA, spanning knee angles from \ang{-119} to \ang{5} and gauge pressures from \qty{256}{\kPa} to \qty{488}{\kPa}. Six entered wrench reductions give knee torques of \qtyrange{28.8}{50.58}{\newton\meter}; three reductions remain pending (Table~\ref{tab:cam_holds}). These measurements document initial load-bearing tests, not a completed validation of the predicted torque envelope. The measured-versus-predicted comparison and the corresponding flexor measurements remain to be completed.''')
# Correct the old shared-route statement against the committed endpoint/route code.
s=s.replace('under the modeling assumption that both actuators share one attachment route rather than a physically mirrored pair; making the mirrored pair explicit is the next design refinement.','with a separately solved route for each actuator. The second insertion mirrors the first across the sagittal plane, while its origin preserves the first route\'s lateral endpoint offset; its wrapping and clearance are evaluated independently.')
s=s.replace('The flexor route, built on the same two-BPA arrangement sharing one attachment route,','The flexor route, built from two separately evaluated BPA paths,')
s=s.replace('These redesign solutions are the model definitions carried into the neuromechanical simulation work of Sections~\\ref{sec:preliminary_sim_methods}--\\ref{sec:preliminary_sim_results} and Chapter~\\ref{ch:futurework}.','These are mechanical design predictions. Their complete incorporation into the converted neuromechanical models remains future integration work.')
# Remove duplicate optimizer-winner history from the final-design paragraph.
s,n=re.subn(r'; the optimizer winner before that adjustment carried .*? at the same checks','',s);assert n==1
archive.append(r'The unadjusted cam optimizer solution had minimum margins of $+31.30$\,\% against the required envelope and $+42.34$\,\% against the plain human target; the adjusted design used in Chapter~\ref{ch:torque} has margins of $+30.83$\,\% and $+39.77$\,\%, respectively.')
a+='\n'+r'''\section{Supporting Route Searches and Transmission Sweeps}\label{sec:app_routing_extensions}

The following diagnostics describe earlier optimizer runs and integration status. The five-site batch table concerns the original unity-transmission search; it is distinct from the subsequent size re-solves plotted below. Objective values from different formulations are not directly comparable.

'''+ '\n\n'.join(archive)+'\n'
for name,label,caption in [
('pulley_torque_resolved.pdf','fig:PulleySurfaces',r'Predicted insertion-torque magnitude versus joint angle and supply pressure for two BPAs through a two-part tackle (force gain 1, travel gain 2), with actuator sizes re-solved. Panels: (A) soleus, (B) tibialis anterior, (C) biceps femoris short head, (D) medial gastrocnemius with ankle fixed at zero, and (E) rectus femoris evaluated at the hip. Dashed curves mark equality with the maximum-isometric OpenSim requirement; no equality curve lies within the soleus panel.'),
('pulley_margin_resolved.pdf','fig:PulleyMargins_resolved',r'Torque margins at \qty{620}{\kPa} for the single-BPA baseline, two BPAs through a two-part tackle, and four BPAs through a four-part tackle. Actuator sizes are re-solved for each configuration. A positive fractional margin exceeds the maximum-isometric OpenSim requirement; the dashed line marks zero. Panels and evaluation planes follow Figure~\ref{fig:PulleySurfaces}.')]:
    a+='\n'+r'\begin{figure}[p]'+'\n'+r'\centering'+'\n'+r'\includegraphics[width=\textwidth]{'+name+'}\n'+r'\caption{'+caption+'}\n'+r'\label{'+label+'}\n'+r'\end{figure}'+'\n'
    d=fitz.open(C/'figs/30-torque'/name);d[0].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(C/'review'/name.replace('.pdf','.png'));print(name,d[0].get_text()[:1800].encode('ascii','backslashreplace').decode())
a+='\n'+r'''\clearpage
\section{Initial Cam-Route Isometric Holds}

Table~\ref{tab:cam_holds} records the nine holds on the rebuilt cam route. Load-cell angle is measured from the tibia axis. The measured reaction-point distance is \qty{259.9}{\mm}; the reduced torque refers to the moving knee center. Dashes indicate reductions not yet entered, rather than zero torque.

\begin{table}[htbp]
\centering\small
\caption{Recorded cam-route holds. The available torque reductions do not yet constitute a measured-versus-predicted envelope comparison.}
\label{tab:cam_holds}
\begin{tabular}{@{}rrrrrr@{}}
\toprule
Hold & Load (N) & Knee angle ($^\circ$) & Pressure (kPa) & LC angle ($^\circ$) & Torque (N\,m)\\
\midrule
1 & 180.2 & -26.2 & 437 & 28 & --\\
2 & 181.5 & -30 & 420 & 26 & --\\
3 & 168 & -61 & 425 & 23 & --\\
4 & 181 & -43 & 435 & 27 & 50.58\\
5 & 116 & -119 & 256 & 10.5 & 32.8\\
6 & 100 & 5 & 480 & 21 & 28.8\\
7 & 126 & -1 & 488 & 27 & 36.86\\
8 & 173.4 & -22 & 435 & 28 & 49.95\\
9 & 140 & -86.5 & 350 & 22 & 39.71\\
\bottomrule
\end{tabular}
\end{table}
'''
p.write_text(s,encoding='utf-8');app.write_text(a,encoding='utf-8')
(C/'review/04_Torque.md').write_text('Corrected the false pulley-position caption: these are joint-angle/pressure and transmission-configuration surfaces. The companion paper plain PDFs are stale; converted its current EPS exports, whose headings explicitly state force gain 1/travel x2, using MiKTeX Ghostscript. Supporting surfaces, original batch diagnostics, adductor conversion status, and whole-body route-builder history moved to Appendix C. Main text retains the final knee design and distinguishes the old unity-transmission search from the subsequent size re-solves. Corrected shared flexor-route claim against flexorBpa2Endpoints20mm.m and predictKneeFlexor20mm.m. New cam table transcribes only workbook values from ExtTest20mm_2; cached row 16 torque was not substituted for row 17 actual torque. The earlier three direct-route holds and pressure-band comparison were omitted from cam validation: different routes/loads do not validate the predicted envelope. Final cam metrics retained; superseded optimizer metrics moved to Appendix C. Full baseline text retained in review/before.\n',encoding='utf-8')
manifest={}
for old,new in [('PulleyTorqueSurfaces_20261020.eps','pulley_torque_resolved.eps'),('PulleyMarginSurfaces_20261020.eps','pulley_margin_resolved.eps')]:
    source=R/'Documentation/Reports and Papers/Knee_Torque_Test/Figures'/old;shutil.copy2(source,C/'figs/30-torque'/new);manifest[new]={'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
(C/'review/torque_figure_sources.json').write_text(json.dumps(manifest,indent=2))
