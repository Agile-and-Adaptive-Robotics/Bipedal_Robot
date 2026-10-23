from pathlib import Path
import re,shutil,json
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/user_sweep_20261023'
def read(name):return (C/'chapters'/name).read_text(encoding='utf-8-sig')
def save(name,s):(C/'chapters'/name).write_text(s,encoding='utf-8')
def block(s,kind,key):return next(m.group() for m in re.finditer(r'\\begin\{('+kind+r')\}.*?\\end\{\1\}',s,re.S) if '\\label{'+key+'}' in m.group())
def caption(b,short,long):
 start=b.index(r'\caption');i=b.index('{',start);depth=1;j=i+1
 while depth:
  if b[j]=='{':depth+=1
  elif b[j]=='}':depth-=1
  j+=1
 return b[:start]+r'\caption['+short+']{'+long+'}'+b[j:]
# Preserve phone originals; make descriptive, exact copies.
photos=[('1000041724.jpg','15-background/festo_bpas_inflated_deflated.jpg'),('1000041706.jpg','30-torque/cam_extensor_hardware.jpg')]
for original,rel in photos:shutil.copy2(R/'Documentation/Reports and Papers/Dissertation/Figures/components'/original,C/'figs'/rel)
s=read('15-background.tex');b=block(s,'figure','fig:festo_bpa_cad')
new=r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.78\textwidth,height=0.54\textheight,keepaspectratio]{festo_bpas_inflated_deflated.jpg}
\caption[Festo BPAs with and without inflation]{Festo BPAs in the test fixture: the inflated actuator is wider and shorter than the deflated actuator beside it.}
\label{fig:festo_bpa_cad}
\end{figure}'''
s=s.replace(b,new).replace('shows the commercial actuator assembly used in the laboratory.','shows the change in actuator shape with inflation.')
s=s.replace(r'\textbf{(a)} Sagittal',r'\textbf{(A)} Sagittal').replace(r'\textbf{(b)} Exploded',r'\textbf{(B)} Exploded')
b=block(s,'figure','fig:aarl_opensim_comparison')
new=r'''\begin{figure}[htbp]
\centering
\begin{minipage}[t]{0.30\textwidth}\vspace{0pt}
\raggedright\textbf{(A)}\par\smallskip
\centering\includegraphics[width=\linewidth,height=0.40\textheight,keepaspectratio]{aarl_bpa_biped.png}
\end{minipage}\hfill
\begin{minipage}[t]{0.30\textwidth}\vspace{0pt}
\raggedright\textbf{(B)}\par\smallskip
\centering\includegraphics[width=\linewidth,height=0.40\textheight,keepaspectratio]{robot_pam_routes.pdf}
\end{minipage}\hfill
\begin{minipage}[t]{0.30\textwidth}\vspace{0pt}
\raggedright\textbf{(C)}\par\smallskip
\centering\includegraphics[width=\linewidth,height=0.40\textheight,keepaspectratio]{human_gait2392_paths.pdf}
\end{minipage}
\caption[Robot skeleton and robot/human muscle paths]{\textbf{(A)} Laboratory robot skeleton. \textbf{(B)} Robot PAM routes. \textbf{(C)} Human Gait2392 muscle paths used as design references \citep{bolen_bipedal_2019, morrow_optimization_2020}.}
\label{fig:aarl_opensim_comparison}
\end{figure}'''
s=s.replace(b,new)
s=s.replace("Figure~\\ref{fig:aarl_opensim_comparison} shows the two sides of that comparison: the laboratory's bipedal skeleton and the OpenSim muscle-routing reference developed for it.","Figure~\\ref{fig:aarl_opensim_comparison} places the laboratory skeleton beside the robot PAM routes and human muscle paths used in that comparison.")
b=block(s,'figure','fig:masters_routing_and_arms')
new=r'''\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{moment_arm_methods.pdf}
\caption[Equivalent geometric constructions of joint torque]{Equivalent fixed-axis constructions. \textbf{(A)} The cross product gives the moment of a force about the joint. \textbf{(B)} Projecting a joint-to-path vector onto the force-path normal gives the signed moment arm used in Equation~\eqref{eq:project_ma}. The joint axis points out of the page.}
\label{fig:masters_routing_and_arms}
\end{figure}'''
s=s.replace(b,new).replace('show the geometry and the torque comparison.','show the moment-arm construction and the torque comparison.')
save('15-background.tex',s)
# Chapter 3: layout only; every source plot stays byte-identical.
s=read('20-muscles.tex')
for key,height in [('fig:Fmax_fit','0.29'),('fig:F-P_contraction','0.29'),('fig:GenForceFit','0.27')]:
 b=block(s,'figure',key);n=b.replace(r'\begin{figure}[H]',r'\begin{figure}[p]')
 n=n.replace('%[tbp]','');n=n.replace(r'\subfloat',r'\subfloat')
 n=re.sub(r'\\includegraphics\[width=[^\]]+\]',r'\\includegraphics[width=0.94\\textwidth,height='+height+r'\\textheight,keepaspectratio]',n)
 n=n.replace(r'\begin{figure}[p]',r'\begin{figure}[p]'+'\n'+r'\centering')
 s=s.replace(b,n)
b=block(s,'table','tab:model_compare');data=b.split(r'\midrule',1)[1].split(r'\bottomrule',1)[0]
# Unpivot the existing 4x3 results into readable rows, without changing values.
rows=[]
raw=re.findall(r'\\textbf\{([^}]*)\}\s*&\s*(\d+)\s*&([^\n]+)',data)
diam=''
for diameter,length,values in raw:
 if diameter:diam=diameter
 vals=[v.strip() for v in values.split('&')];vals[-1]=vals[-1].split(r'\\')[0].strip()
 for model,j in [('Bolen',0),("S\'arosi",3),('Martens and Boblan',6)]:rows.append(f'{diam} & {length} & {model} & '+' & '.join(vals[j:j+3])+r'\\')
assert len(rows)==12
n=r'''\begin{table}[p]
\centering\small
\caption{BPA force-model comparison. Diameter and resting length are in mm; RMSE and maximum absolute error are in N. FVU is dimensionless.}
\label{tab:model_compare}
\begin{tabular}{@{}rrlrrr@{}}
\toprule
Diameter & Length & Model & RMSE & FVU & Max. error\\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}'''
s=s.replace(b,n);save('20-muscles.tex',s)
s=read('30-torque.tex')
old='The pinned-flexor data do not uniquely determine the stiffness pair, because a family of pairs reproduces the training tests nearly equally well, so the adopted pair was selected by a requirement external to the fit: both redesign configurations must meet the required 5\\,\\% torque-overshoot target simultaneously. '
assert old in s;s=s.replace(old,'')
s=s.replace(r'$\chi_{0}$ (mm)',r'$\chi_{0}$ (m)').replace(r'Flexor & $+3.94$',r'Flexor & $+0.00394$').replace(r'Extensor & $-6.38$',r'Extensor & $-0.00638$')
# Critical scoring conditions move out of the table caption and into result prose.
scoring=r'''For Table~\ref{tab:FlxPinPreGof}, the rigid-body and improved-characterization predictions are scored on the same points where the rigid-body strain is defined: 75 of 90, 17 of 20, 14 of 19, and 15 of 22 points in panels A--D, respectively. The hybrid calculation is scored on all measured points. These pre-optimization scores therefore use a different point set from the full-test corrected-model scores in Table~\ref{tab:FlxPinGof}.

'''
b=block(s,'sidewaystable','tab:FlxPinPreGof');s=s.replace(b,scoring+caption(b,'Pre-optimization pinned-flexor goodness of fit',r'Pre-optimization pinned-flexor goodness of fit. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. Scoring-point restrictions are stated in the text.'))
for key,title in [('tab:FlxPinGof','Optimized pinned-flexor goodness of fit'),('tab:ExtPinGof','Optimized pinned-extensor goodness of fit')]:
 b=block(s,'table',key);s=s.replace(b,caption(b,title,title+r' at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.'))
# Fold membership is already explicit in the rows; retain exclusions in nearby prose.
s=s.replace('Table~\\ref{tab:FlxPinGof} completes this comparison',r'The \qty{47.9}{\cm} flexor test was not used in identification and is evaluated separately in Appendix~\ref{app:model_geometry}. Table~\ref{tab:FlxPinGof} completes this comparison',1)
# Explain design bounds before presenting routes, using Ben's account, without inventing bounds.
anchor='The route redesign built on these values addresses'
pos=s.index(anchor);s=s[:pos]+r'''The initial route searches bounded the BPA origins and insertions near the corresponding human muscle attachments to preserve biomimetic routing. Those anatomical bounds had to be relaxed where the available BPA contraction could not cover the knee range of motion. The resulting routes then required a second tradeoff: reducing the actuator forces carried by the attachment hardware after parts failed during testing. The cam redesign changes the moment arm and contact geometry to retain the torque requirement at lower actuator force.

'''+s[pos:]
s=s.replace('The cam solution converts the razor-thin margin of the direct route into genuine design headroom. On the final geometry, which folds the optimizer winner together with a verified adjustment of the distal route point and tendon length,','The final cam geometry includes the optimized route and the verified adjustment of the distal route point and tendon length. Its')
s=s.replace('Its the minimum torque margin','Its minimum torque margin')
b=block(s,'figure','fig:ExtensorRoute_redesigned');n=b.replace(r'\caption[Cam-extensor route geometry]',r'\par\medskip'+'\n'+r'\includegraphics[width=0.46\textwidth,height=0.31\textheight,keepaspectratio]{cam_extensor_hardware.jpg}'+'\n'+r'\caption[Cam-extensor route geometry and hardware]');n=n.replace('Cam-extensor route at','Cam-extensor route (top) at').replace('Panel letters follow the full sequence','The photograph (bottom) shows the assembled two-BPA hardware. Panel letters follow the full sequence');s=s.replace(b,n)
# Reuse Figure 2.6 to explain Eq.4.2 without shifting Chapter 4 figure numbers.
s=s.replace('Moment arm is calculated in the tibial body frame.',r'Here $\mathbf{p}_{f}$ is the force-path segment vector, and $\mathbf{p}_{proj,xy}$ locates a point on that line relative to the joint, projected onto the plane normal to the joint axis. The normalized cross product in Equation~\eqref{eq:project_ma} is the in-plane normal to the force path. Its dot product with the projected position gives the signed moment arm (Figure~\ref{fig:masters_routing_and_arms}). Moment arm is calculated in the tibial body frame.')
s=s.replace(r'In Eq.~\eqref{eq:project_ma}, $\mathbf{p}_{proj,xy}$ is the free muscle segment projected onto the $XY$--plane, which is normal to the joint axis, and $\mathbf{p}_{f}$ is the projected muscle segment vector. ','')
# Place wrapping diagram in the existing deflection figure, preserving numbering.
b=block(s,'figure','fig:xiProjection');n=b.replace(r'\includegraphics[height=0.72\textheight]{xiProjection.pdf}',r'\makebox[\textwidth][l]{\textbf{(A)}}'+'\n'+r'\includegraphics[height=0.40\textheight]{xiProjection.pdf}\par\smallskip'+'\n'+r'\makebox[\textwidth][l]{\textbf{(B)}}'+'\n'+r'\includegraphics[width=0.83\textwidth,height=0.20\textheight,keepaspectratio]{bpa_wrapping_loss.pdf}')
n=caption(n,'Compliance projections and BPA wrapping geometry',r'\textbf{(A)} Bracket displacement and its projection along the force path; tendon elongation acts in series. \textbf{(B)} Wrapping radius $R$ and bend angle $|\theta_{\mathrm{wrap}}-\theta_{\mathrm{k}}|$ in the length-loss correction of Equation~\eqref{eq:deltaL}. The lost usable length is a fitted correction, not the difference between arc and chord lengths.')
s=s.replace(b,n);s=s.replace('A wrapping-loss term $\\chi_{3}$ was introduced',r'Figure~\ref{fig:xiProjection}\textbf{B} shows the wrapping geometry. A wrapping-loss term $\chi_{3}$ was introduced')
save('30-torque.tex',s)
s=read('50-database.tex').replace('lack of a structured account, available to this project,','lack of a structured database, available to this project,').replace("The collection began with my Zotero library, the laboratory's Biology library, and the papers already curated in Airtable.","The collection began with the Zotero libraries of the AARL and myself, and papers already curated in Airtable.")
start=s.index('I evaluated AI-assisted');end=s.index('\n\n',start)
s=s[:start]+r'''I evaluated AI-assisted literature identification, summarization, and categorization in work presented at Neuroscience 2024 \citep{bolen_using_2024}. The comparison included ResearchRabbit, Litmaps, txyz, SciSpace, and MirrorThink. Some tools helped identify related papers or summarize straightforward material, but categorization remained inadequate: none received a good rating, and only txyz was rated adequate. Complex models and nuanced physiological findings required substantial human correction. The tools therefore offered limited help with the central curation task of accurately assigning papers to database categories. I retained responsibility for inclusion, classification, and deletion decisions after checking the source.'''+s[end:];save('50-database.tex',s)
s=read('60-sns.tex')
s=s.replace('The controller separated rhythm generation from pattern formation and included muscle-length, velocity, and force feedback.','The initial controller used one rhythm generator (RG) to control both sides, with a separate pattern-formation layer. I added afferent feedback and heel- and toe-contact signals during the subsequent development.')
s=s.replace("The muscle-feedback gains were not systematically tuned, and the combined model did not include the Li walker's ground-contact feedback.","This suspended stepping result preceded the later feedback and bilateral-RG extensions; it does not demonstrate walking with those additions.")
s=s.replace('I constructed the right-side subsystems from the left-side organization and linked the two RGs through commissural pathways.','I constructed the right-side subsystems from the left-side organization and linked the two RGs through commissural pathways. This later model is referred to as the Walker Two-Layer (W2L) model.')
s=s.replace('The missing step was a tested link between spinal neural circuits and the human-scale and BPA-actuated mechanical models needed to study this robot.','The missing step was a tested link between spinal neural circuits and the human-scale and BPA-actuated mechanical models needed to study this robot, including the many joint degrees of freedom and 92 musculotendon actuators represented in Gait2392.')
s=s.replace('AnimatLab retains the reference controllers but is no longer maintained.','AnimatLab retains the reference controllers but is unsupported, relies on an older physics engine, and ran these simulations slowly.')
start=s.index('Compared with the earlier SNS models');end=s.index(' Figure~\\ref{fig:simulink_sns_implementation}',start)
s=s[:start]+'The library retains the conductance-based neuron and synapse formulation supplied by SNS-Toolbox; it does not introduce a new neuron model. Verification uses selected circuits under matched inputs.'+s[end:]
# Delete the known-false welded-joint account. Leave implementation questions to ZCode.
start=s.index('Two mechanical import routes were investigated.');end=s.index('\n\n',start)
s=s[:start]+'OpenSim-derived and SolidWorks-derived mechanisms were imported into Simscape Multibody. These mechanical imports are separate from the standalone reflex demonstrations.'+s[end:]
s=s.replace("Mechanism import achieved; the SolidWorks knee's welded joint and incomplete muscle integration remain limitations",'Mechanism import achieved; standalone reflex demonstrations do not establish control of the imported models')
s=s.replace('I developed the conversion, coupling, graphical-library, and controller-port tools below','I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below')
s=s.replace('MuJoCo--SNS coupling &','Models in MuJoCo and SNS-Toolbox &').replace('Controller ports and test scripts &','Matched circuit implementations &')
save('60-sns.tex',s)
# Also remove the same fit-selection narrative from Appendix C, preserving its old form in snapshot.
s=read('94-AppendixC.tex')
s=re.sub(r'The pinned-flexor data do not uniquely determine.*?The adopted pair was chosen because.*?\(extensor\)\.', 'The adopted values and their front-row provenance are recorded below.',s,flags=re.S)
save('94-AppendixC.tex',s)
print('Applied prose, table, photograph, and figure-layout changes')
