from pathlib import Path
import re, shutil, csv
import pandas as pd
R=Path(r'D:\GitHub\Bipedal_Robot'); O=R/'tmp/torque_review_20261027'
D=R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex'
A=D.with_name('94-AppendixC.tex');B=D.with_name('15-background.tex')
P=R/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex';S=P.with_name('frontiers_SupplementaryMaterial.tex')
d,a,p,s,b=[x.read_text(encoding='utf-8') for x in [D,A,P,S,B]]
shutil.copy2(B,O/'15-background_before.tex')
# Ben prefers content completeness over the display-item ceiling for now.
original=(O/'Main_before.tex').read_text(encoding='utf-8')
frameblocks=[]
for lab in ['fig:xiFrameTransform','fig:bktFrame']:
    block=next(m.group() for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',original,re.S) if '\\label{'+lab+'}' in m.group())
    frameblocks.append(block)
    s=re.sub(r'\\begin\{figure\}.*?\\end\{figure\}',lambda m:'' if '\\label{'+lab+'}' in m.group() else m.group(),s,flags=re.S)
    p=p.replace('Supplementary Figure~\\ref{'+lab+'}','Figure~\\ref{'+lab+'}')
p=p.replace('In each local frame the bracket is modeled', '\n\n'.join(frameblocks)+'\n\nIn each local frame the bracket is modeled')
s=s.replace('\\subsection{Bracket Reference Frames}','')

# Remove legacy excluded-test narrative from the chapter's result exposition.
st=d.index('Results are shown in Fig. \\ref{fig:ExtPin_group}. The initial identification')
en=d.index('\n\n',st)
d=d[:st]+r'''An initial single-bracket extensor fit yielded $\chi_0=\qty{-10.9}{\mm}$ and $\chi_3=0.123$; its historical goodness-of-fit values are retained in Table~\ref{tab:TorqueFit}. The adopted two-stage refit instead uses the shared flexor stiffness pair and the six-test fold defined in Section~\ref{sec:identification}. Figure~\ref{fig:ExtPin_group} and Table~\ref{tab:ExtPinGof} report that adopted fit.'''+d[en:]
d=d.replace('at 0.172, on the adopted fold.', 'at 0.172.')
d=d.replace('The methods, results, and discussion follow the manuscript noted above; the physical-validation section additionally reports the newer sweep results, which postdate the manuscript draft.', 'The methods, results, and discussion follow the manuscript noted above and include the physical-validation sweeps of the redesigned routes.')
d=d.replace('then the stiffness-aware identification of the correction terms', 'then the identification of the correction terms')

# costFunctionKnee.m minimizes -Ct - 10*Ca, where Ct averages all three
# torque components and Ca = mean(dot(uH,uR)-1). The previous equation had
# the alignment sign and torque normalization wrong for a maximized score.
d=d.replace('The attachment locations of the BPAs are themselves design variables.', 'The attachment locations of the BPAs are themselves design variables. Morrow and I developed the mesh-placement method as part of this work \\citep{morrow_optimization_2020}.')
d=d.replace('The score it maximizes rewards torque surplus over the human curve in the direction the human curve demands, together with alignment of the robot and human torque vectors,', 'The implementation minimizes a cost that rewards torque surplus in the required direction and penalizes misalignment of the human and robot torque vectors,')
old=r'''S = \frac{1}{N}\sum_{j=1}^{N}\left[\,w_{\mathbf{k}}\,D_{\mathbf{k},j} + \sum_{a} w_{a}\,D_{a,j}\right] + \frac{w_{u}}{N}\sum_{j=1}^{N}\left(1 - \hat{\mathbf{m}}^{\mathrm{H}}_{j}\cdot\hat{\mathbf{m}}^{\mathrm{R}}_{j}\right),'''
new=r'''C_{\mathrm{mesh}} = -\frac{1}{3N}\sum_{j=1}^{N}\left[\,w_{\mathbf{k}}\,D_{\mathbf{k},j} + \sum_{a\in\{x,y\}} w_{a}\,D_{a,j}\right] + \frac{w_{u}}{N}\sum_{j=1}^{N}\left(1 - \hat{\mathbf{m}}^{\mathrm{H}}_{j}\cdot\hat{\mathbf{m}}^{\mathrm{R}}_{j}\right).'''
assert old in d;d=d.replace(old,new)
needle='The knee-axis component carries weight $w_{\\mathbf{k}} = 5$, the off-axis components carry $w_{a} = 0.25$, and the alignment term carries $w_{u} = 10$.'
d=d.replace(needle,needle+' The factor $3N$ matches the implementation\'s average over the three torque components at each pose. A zero robot-torque vector is assigned the direction opposite the human vector for the alignment term. Table~\\ref{tab:mesh_original_weights} lists the terms of this cost.')
tbl=r'''
\begin{table}[htbp]
\caption{Terms of the compliance-unaware knee-placement cost. The cost is minimized; a surplus reduces it, while a direction mismatch increases it.}
\label{tab:mesh_original_weights}
\begin{tabularx}{\textwidth}{@{}p{0.24\textwidth}p{0.12\textwidth}X@{}}
\toprule
\textbf{Term} & \textbf{Weight} & \textbf{Meaning} \\
\midrule
Knee torque surplus & $5$ & Signed surplus $D_{\mathbf{k},j}$, summed with the off-axis terms and divided by $3N$; enters with a negative sign. \\
Off-axis torque surplus & $0.25$ each & Signed surplus $D_{x,j}$ and $D_{y,j}$, using the human torque's sign on each axis; enters with a negative sign. \\
Direction mismatch & $10$ & Mean $1-\hat{\mathbf{m}}^H_j\cdot\hat{\mathbf{m}}^R_j$; zero for aligned vectors and two for opposing vectors. \\
\bottomrule
\end{tabularx}
\end{table}
'''
d=d.replace("The dissertation's route redesign re-derives",tbl+"\nThe dissertation's route redesign re-derives")
# Place the PhD placement-paper citations in the torque chapter.
b=b.replace(' \\citep{morrow_optimization_2020}', '')
b=b.replace(', morrow_optimization_2020','')

# Make the axis-mean comparison explicit using the already saved CSV.
metric=r'''For comparison, the average of the two axis-wise mean coverages is \qty{38.39}{\percent} for the baseline and \qty{66.17}{\percent} for the updated configurations. On a complete rectangular grid, averaging row and column coverages equals the overall cell fraction exactly; it cannot reduce the influence of edge poses there. The two summaries differ only on the four hip-envelope domains with excluded cells. The two metrics and the torque-ratio score are reported together, so their different weightings remain visible.'''
for name in ['d','p']:
    t=locals()[name]
    t=t.replace('These two metrics support different conclusions and should be considered together.', 'These metrics support different conclusions and should be considered together.\n\n'+metric)
    t=t.replace('same-named single muscle of the Gait2392 model on a $5\\times5$ grid', 'same-named single muscle of the Gait2392 model \\citep{delp_opensim_2007,seth_opensim_2018} on a $5\\times5$ grid')
    t=t.replace('The comparison includes all 27 actuator sites.', 'The BPA force model is evaluated at \\qty{620}{\\kPa} \\citep{bolen_2026}. The comparison includes all 27 actuator sites.')
    locals()[name]=t
p=p.replace('\\begin{tabular}{@{}lcc@{}}\n\\toprule\n\\textbf{Routing} & \\textbf{Torque meets or exceeds human (\\%)} & \\textbf{Mean torque-ratio score}', '\\begin{tabular}{@{}lccc@{}}\n\\toprule\n\\textbf{Configuration} & \\textbf{Grid coverage (\\%)} & \\textbf{Axis mean (\\%)} & \\textbf{Torque-ratio score}')
p=p.replace('Baseline & 37.28 & $-1.93$', 'Baseline & 37.28 & 38.39 & $-1.93$')
p=p.replace('Transmission-ladder winners & 65.19 & $-3.63$', 'Updated & 65.19 & 66.17 & $-3.63$')
p=p.replace('Whole-Body Capability Against Human Torque Surfaces','Whole-Leg Capability Against Human Torque Surfaces')
p=p.replace('two-\\qty{20}{\\mm}-BPA configurations', 'configurations with two $\\diameter$\\qty{20}{\\mm} BPAs')
p=p.replace('cut held-out prediction error by two thirds or more on the pinned-knee test families', 'reduce held-out prediction error on both pinned-knee test families')
p=p.replace('A robot designer adds effects that no biological model carries:', 'The robot hardware adds further sources of error:')
p=p.replace('What these models do not provide is a torque prediction that survives contact with a physical joint.', 'The present tests evaluate errors that remain when these force and routing models are applied to a physical joint.')
p=p.replace('The terms validated here attach to the force characterization', 'The correction terms used here attach to the force characterization')
# Restoring frames raises this reference count without losing any figures.
p=p.replace('a transfer whose testing remains future work (Section~\\ref{sec:limitations})', 'a transfer evaluated in the redesigned flexor\'s isometric tests')

# Tables in Appendix C previously ran beyond the margin. Wrap descriptive
# columns rather than shrinking scientific text below the figure font floor.
def tablex(text,label,spec,changes=None):
    ms=list(re.finditer(r'\\begin\{table\}.*?\\end\{table\}',text,re.S))
    m=next(m for m in ms if '\\label{'+label+'}' in m.group())
    block=m.group()
    block=re.sub(r'\\begin\{tabular\}\{[^\n]+\}',lambda _:r'\begin{tabularx}{\textwidth}{'+spec+'}',block,count=1)
    block=block.replace('\\end{tabular}','\\end{tabularx}')
    if changes:
        for x,y in changes:block=block.replace(x,y)
    return text[:m.start()]+block+text[m.end():]
a=tablex(a,'tab:app_bracket_points',r'@{}lrrrX@{}')
a=tablex(a,'tab:app_stiffness_arrays',r'@{}p{0.30\textwidth}p{0.22\textwidth}X@{}')
a=tablex(a,'tab:app_opt_bounds',r'@{}p{0.20\textwidth}p{0.13\textwidth}Xp{0.14\textwidth}@{}')
a=tablex(a,'tab:app_gof',r'@{}p{0.24\textwidth}Xrrr@{}',[(r'\textbf{Rigid RMSE}',r'\shortstack{\textbf{Rigid}\\\textbf{RMSE}}')])
a=tablex(a,'tab:bipulley_batch',r'@{}p{0.22\textwidth}p{0.15\textwidth}X@{}')
# Permit line breaks in the long literal filenames in Appendix C prose.
a=re.sub(r'\\texttt\{([^{}]*(?:\.m|\.mat|\.log|\.csv)[^{}]*)\}',lambda m:r'\nolinkurl{'+m.group(1).replace(r'\_','_')+'}',a)

# Provide the requested per-site DoFs and sizing without running MATLAB.
rows=list(csv.DictReader((R/'tmp/paper_bridge/capability_table_20261022/BiPulley_survey2_perMuscle_20261023.csv').open(encoding='utf-8')))
assert len(rows)==27
abbr={'lumbar_extension':'lumbar ext.','lumbar_bending':'lumbar bend.','hip_flexion_r':'hip flex.','hip_adduction_r':'hip add.','hip_rotation_r':'hip rot.','knee_angle_r':'knee','ankle_angle_r':'ankle','subtalar_angle_r':'subtalar','mtp_angle_r':'MTP'}
def esc(t):return t.replace('_',r'\_')
lines=[]
for row in rows:
    pri=abbr.get(row['primary_dof'],row['primary_dof'].replace('_r','').replace('_',' '))
    sec=abbr.get(row['secondary_dof'],row['secondary_dof'].replace('_r','').replace('_',' '))
    if sec in ['','none','NA']:sec='--'
    dia=int(float(row['bpa_size_mm']));baseline='4'+r'$\times$'+'20' if dia==40 else '1'+r'$\times$'+str(dia)
    lines.append(f"{esc(row['site'])} & {esc(pri)} & {esc(sec)} & {dia} & {baseline} \\\\")
inventory=r'''\begin{table}[p]
\centering
\caption{Degrees of freedom and nominal BPA sizes in the 27-site configuration survey. Diameter follows the inherited robot sizing rule and is retained on the updated side. At nominal \qty{40}{\mm} sites, the baseline instead uses four parallel \qty{20}{\mm} BPAs; other baseline sites use one actuator of the listed diameter. The two coordinate columns identify the graded surface, rather than every coordinate a multiarticular route could cross.}
\label{tab:survey_dofs_sizes}
\begingroup\small\renewcommand{\arraystretch}{1.1}
\begin{tabular}{@{}lllll@{}}
\toprule
\textbf{Site} & \textbf{Primary DoF} & \textbf{Secondary DoF} & \textbf{Dia. (mm)} & \textbf{Baseline (mm)} \\
\midrule
'''+ '\n'.join(lines)+r'''
\bottomrule
\end{tabular}
\endgroup
\end{table}
'''
a+='\n\\section{Configuration-Survey Coordinates and Actuator Sizes}\n\n'+inventory
s=s.replace('\\end{document}', '\\subsection{Configuration-Survey Coordinates and Actuator Sizes}\n\n'+inventory+'\n\\end{document}')
d=d.replace('The comparison includes all 27 actuator sites.', 'The comparison includes all 27 actuator sites (Table~\\ref{tab:survey_dofs_sizes}).')
p=p.replace('The comparison includes all 27 actuator sites.', 'The comparison includes all 27 actuator sites (Supplementary Table~\\ref{tab:survey_dofs_sizes}).')
for path,text in zip([D,A,P,S,B],[d,a,p,s,b]):path.write_text(text,encoding='utf-8',newline='\n')
print('Main-paper display items:',len(re.findall(r'\\begin\{(?:figure|table)\}',p)))
print('Inventory rows:',len(lines))
