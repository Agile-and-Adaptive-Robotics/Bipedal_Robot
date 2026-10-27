from pathlib import Path
import re, shutil, difflib, json

ROOT = Path(r'D:\GitHub\Bipedal_Robot')
OUT = ROOT/'tmp/torque_review_20261027'
OUT.mkdir(parents=True, exist_ok=True)
D = ROOT/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex'
A = D.with_name('94-AppendixC.tex')
P = ROOT/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex'
S = P.with_name('frontiers_SupplementaryMaterial.tex')
paths = [D,A,P,S]
before = {}
for p in paths:
    before[p] = p.read_text(encoding='utf-8')
    shutil.copy2(p, OUT/(p.stem+'_before.tex'))
d,a,p,s = [before[x] for x in paths]

def replace_once(text, old, new):
    assert text.count(old)==1, (old[:90],text.count(old))
    return text.replace(old,new,1)

# Keep frame construction details available without increasing the paper's
# combined display-item count when the two explanatory schematics return.
moved=[]
for label in ['fig:xiFrameTransform','fig:bktFrame']:
    blocks=list(re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',p,re.S))
    block=next(m.group() for m in blocks if '\\label{'+label+'}' in m.group())
    moved.append(block.replace('[H]','[htbp]'))
    p=p.replace(block+'\n','',1)
    p=p.replace('Figure~\\ref{'+label+'}', 'Supplementary Figure~\\ref{'+label+'}')
    p=p.replace('Figure \\ref{'+label+'}', 'Supplementary Figure~\\ref{'+label+'}')
s=s.replace('\\end{document}', '\\subsection{Bracket Reference Frames}\n\n'+'\n\n'.join(moved)+'\n\n\\end{document}')
s=s.replace('\\usepackage{amssymb}','\\usepackage{amssymb}\n\\usepackage{mathabx}')
s=s.replace('\\textbf{$\\phi$}', '\\textbf{$\\diameter$}')
s=s.replace('Comparison of pre- and post-optimization torque prediction fits for BPAs', 'Historical single-bracket torque-prediction fits, preceding the adopted two-stage identification, for BPAs')

# Explain the vector calculation before reducing it to a scalar force balance.
start=p.index('Only the deflection component along the line of action')
end=p.index('\\begin{equation}\n\\delta_{br}',start)
p=p[:start]+r'''The local force vector is the force magnitude times its unit direction,
\begin{equation}
\mathbf{F}_{br}=F\hat{\mathbf{u}}, \qquad \hat{\mathbf{u}}=[u_x,u_y,u_z]^T.
\label{eq:bktforcevector}
\end{equation}
Hooke's law relates this force to the bracket displacement vector:
\begin{equation}
\mathbf{F}_{br}=\mathbf{K}_{br}\boldsymbol{\eta}_{br},\qquad
\boldsymbol{\eta}_{br}=\mathbf{C}_{br}\mathbf{F}_{br},\qquad
\mathbf{C}_{br}=\mathbf{K}_{br}^{-1}.
\label{eq:bktdeflection}
\end{equation}
For the insertion-bracket assignment of Equation~\eqref{eq:bktstiffness}, this gives
\begin{equation}
\boldsymbol{\eta}_{br}=F
\begin{bmatrix}u_x/\chi_1\\u_y/\chi_2\\u_z/\chi_1\end{bmatrix}.
\label{eq:bktcomponentdeflection}
\end{equation}
The displacement need not be parallel to the force. Only its projection along the local line of action contributes to the scalar shortening used in the force balance (Figure~\ref{fig:xiProjection}):
\begin{equation}
\delta_{br}=\hat{\mathbf{u}}^T\boldsymbol{\eta}_{br}.
\label{eq:bktprojecteddeflection}
\end{equation}
Substitution of Equation~\eqref{eq:bktdeflection} gives
'''+p[end:]
old='In this frame the bracket is modeled as three independent linear springs, and the identified stiffness array is y-symmetric,'
p=replace_once(p,old,'In each local frame the bracket is modeled as three independent linear springs. For the pinned-flexor insertion bracket, the identified stiffness assignment is')
old='with $\\chi_{1}$ the effective stiffness in the $x_{br}$ and $z_{br}$ directions and $\\chi_{2}$ the effective bending stiffness along the cantilever direction $y_{br}$. The symmetry makes the compliance prediction independent of the frame convention: the pitch-only and two-rotation constructions of Supplementary Figure~\\ref{fig:xiFrameTransform} differ by a rotation about $\\hat{\\mathbf{y}}_{br}$, and a diagonal stiffness array is invariant under that rotation exactly when the two stiffnesses in its plane are equal.'
new=r'''with $\chi_{1}$ assigned to $x_{br}$ and $z_{br}$, and $\chi_{2}$ to $y_{br}$. The stiffness pair is shared across configurations, but its axis assignment follows the hardware model: the pinned-flexor origin bracket uses $\operatorname{diag}(\chi_1,\chi_1,\chi_2)$, the pinned extensor and biomimetic flexor use $\operatorname{diag}(\chi_1,\chi_2,\chi_2)$, and the biomimetic extensor uses Equation~\eqref{eq:bktstiffness}. The flexor evaluator constructs both bracket frames with the first rotation alone; the pinned-extensor evaluator uses the two-rotation construction. A rotation about $y_{br}$ leaves a stiffness matrix unchanged only when its $x_{br}$ and $z_{br}$ stiffnesses are equal. Frame independence therefore applies to Equation~\eqref{eq:bktstiffness}, not to every assignment.'''
p=replace_once(p,old,new)
p=p.replace('The frame is constructed by two rotations', 'The two-rotation frame is constructed')
assumptions=r'''This reduction assumes quasi-static tensile loading, linear translational springs, and additive displacements along the force path. The direction is held fixed during the scalar solve, then updated with the displaced attachments; rotational compliance and force--moment coupling are not modeled. With a second compliant bracket, the same reduction gives
\begin{equation}
k_{eq}^{-1}=\sum_j\hat{\mathbf{u}}_j^T\mathbf{K}_{br,j}^{-1}\hat{\mathbf{u}}_j+k_{ten}^{-1},
\label{eq:two_bracket_stiffness}
\end{equation}
where the tendon term is omitted for a direct attachment. All unit directions and stiffness matrices in a projection are expressed in the same local frame.

'''
p=p.replace('The compliance displacement shortens the actuator,',assumptions+'The compliance displacement shortens the actuator,')
p=p.replace('Both ends of the actuator path can be compliant:', 'Figure~\\ref{fig:xiProjection} shows how both ends of the actuator path can be compliant:')

projection=r'''\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{FIGPATHxiProjection_review.pdf}
\caption[Compliance projections along the force path]{Local bracket displacement $\boldsymbol{\eta}_{br,j}$ and its projection $\delta_{br,j}$ along the force direction $\hat{\mathbf{u}}_j$. The inset shows that an anisotropic bracket can deflect off the force axis. The two projected bracket displacements and the tendon elongation add to the actuator shortening. Displacements are exaggerated for clarity.}
\label{fig:xiProjection}
\end{figure}
'''
wrapping=r'''\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{FIGPATHbpa_knee_wrapping_review.pdf}
\caption[BPA crossing a flexed knee]{Sagittal schematic of an extensor BPA crossing a flexed robot knee, with its braided body following an anterior contact arc between the proximal attachment and distal tendon. Dashed radii bound the contact arc; $R$ is the BPA-path wrapping radius, and $\alpha=\abs{\theta_{\mathrm{wrap}}-\theta_k}$ is the bend angle after contact begins, expressed in radians in Equation~\eqref{eq:deltaL}. The drawing illustrates the local circular-wrap approximation, not the full multi-contact route. The fitted usable-length loss $\Delta l$ is not an arc-minus-chord length.}
\label{fig:bpaWrapping}
\end{figure}
'''
p=p.replace('\\subsection{Fictive Strain and Wrapping Loss}', projection.replace('FIGPATH','Figures/')+'\n\\subsection{Fictive Strain and Wrapping Loss}')
wrap_explanation=r'''As the knee flexes, the extensor BPA contacts the anterior joint surface and bends around it (Figure~\ref{fig:bpaWrapping}). The rigid routed-length calculation does not account for the resulting change in usable contraction. We represent this effect with fictive strain: an empirical increase in the contraction seen by the force model, beyond the actuator's geometric contraction. The correction is phenomenological; it does not resolve wall stresses or the deformation of the braided actuator. A wrapping-loss term $\chi_3$ scales the associated reduction in usable length:'''
st=p.index('The extensor configuration required one further correction,')
en=p.index('\n\\begin{equation}',st)
p=p[:st]+wrap_explanation+p[en:]
p=p.replace('\\subsection{Identification Procedure}', wrapping.replace('FIGPATH','Figures/')+'\n\\subsection{Identification Procedure}')

# Separate the dissertation's composite float and match the paper's corrected
# identification description, while retaining its appendix cross-references.
st=d.index('\\begin{figure}[htbp]\n\\centering\n\\makebox[\\textwidth][l]{\\textbf{(A)}}')
en=d.index('\\end{figure}',st)+len('\\end{figure}')
d=d[:st]+projection.replace('FIGPATH','')+d[en:]
st=d.index('The extensor configuration required one further correction,')
en=d.index('\n\\begin{equation}',st)
d=d[:st]+wrap_explanation.replace('We represent','I represent')+d[en:]
d=d.replace('\\subsubsection{Stiffness-Aware Multiobjective Identification}',wrapping.replace('FIGPATH','')+'\n\\subsubsection{Identification of the Torque-Model Correction Terms}')
st=d.index('Multiobjective optimization using \\texttt{gamultiobj}')
en=d.index('The objective vector of the identification is defined explicitly.',st)
intro=r'''I identified the correction terms in two stages using \texttt{gamultiobj} in MATLAB. The flexor stage fits $(\chi_0,\chi_1,\chi_2)$ with compliance at both attachments. The adopted fold trains on the $l_0=\qty{48.5}{\cm}$ and \qty{41.0}{\cm} tests and validates on the \qty{45.7}{\cm} and \qty{40.6}{\cm} tendon tests. The extensor stage locks $(\chi_1,\chi_2)$ to that same flexor-front member and fits only $(\chi_0,\chi_3)$. Its adopted fold trains on the \qty{43.6}{\cm}, \qty{43.2}{\cm}-tendon, and \qty{45.7}{\cm} tests and validates on the \qty{40.0}{\cm} plain, \qty{40.0}{\cm}-tendon, and \qty{46.5}{\cm} tests. The extensor length offset is allowed to be negative. Appendix~\ref{app:model_geometry} records the test constants, bounds, stiffness-axis assignments, and adopted values in Tables~\ref{tab:app_flxpin_tests}, \ref{tab:app_extpin_tests}, \ref{tab:app_opt_bounds}, \ref{tab:app_stiffness_arrays}, and~\ref{tab:app_xi_values}.

An initial single-bracket identification preceded the two-bracket model. Its results are retained as historical comparisons, but the two-stage procedure defines the adopted correction terms. Each candidate is evaluated using RMSE, FVU, and maximum absolute residual. These metrics jointly measure torque error; none isolates a unique source such as level offset, waveform shape, or angular misalignment.

'''
d=d[:st]+intro+d[en:]
# Match the metric interpretation to its actual definition (including the
# evaluator's n-1 normalization), without inventing a shape-only metric.
old='The maximum absolute residual records the largest pointwise error. The solver minimizes the objective vector'
d=d.replace(old,'The maximum absolute residual records the largest pointwise error. RMSE and the maximum residual have torque units, whereas FVU is dimensionless. For a given test, FVU is proportional to squared RMSE, so the two metrics are not independent measures of shape and level. Across tests, their different scaling changes the relative weighting of the errors. The solver minimizes the objective vector')
d=d.replace('described above. The production runs', 'using held-out tests. The production runs')
d=d.replace('The root is found with a bracketed scalar solver (\\texttt{fzero} in MATLAB) on the physically admissible interval between zero deflection and the deflection that would bring the actuator to its maximum usable contraction. Termination is controlled by the bracket, not by a residual tolerance: the iteration closes in on the root and stops when the bracket has collapsed to adjacent floating-point values, so the returned $\\delta$ is converged to machine precision. Two safeguards act as secondary stop conditions.', 'The root is found with a bracketed scalar solver (\\texttt{fzero} in MATLAB) on the physically admissible interval between zero deflection and the deflection that would bring the actuator to its maximum usable contraction. The evaluator supplies no custom solver options, so termination uses the default displacement tolerance. This numerical tolerance does not establish the accuracy of the physical model. Two safeguards act as secondary stop conditions.')
st=d.index('These three principal stiffnesses are not independent')
en=d.index('\\begin{equation}\n\\mathbf{K}_{br} = \n',st)
d=d[:st]+r'''The principal stiffnesses are parameterized by the two fitted values $\chi_1$ and $\chi_2$, with axis assignments fixed by configuration (Table~\ref{tab:app_stiffness_arrays}). For the pinned-flexor insertion bracket, the assignment is $k_x=k_z=\chi_1$ and $k_y=\chi_2$. This $y$-symmetric array is invariant to a change of frame orientation about $\hat{\mathbf{y}}_{br}$, so the one- and two-rotation constructions give the same compliance for this assignment. That invariance does not hold for the other assignments. In particular, the pinned-flexor origin bracket uses $[\chi_1,\chi_1,\chi_2]$, the pinned extensor and biomimetic flexor use $[\chi_1,\chi_2,\chi_2]$, and the biomimetic extensor uses the insertion-bracket assignment. Both pinned-flexor brackets use the one-rotation frame in the evaluator of record; the pinned extensor uses the two-rotation frame. The general projected-compliance calculation below applies to each assignment, with the appropriate local matrix. For the insertion bracket,

'''+d[en:]
d=d.replace('a compliant coupling between them would break this.', 'an additional compliant coupling would require its own term in the series compliance.')
d=d.replace('Here $k_{eq}$ is the series sum of every compliant element', 'Here $k_{eq}$ is the equivalent stiffness of every compliant element')
d=d.replace('which is an approximation, and one the results justify:', 'which is a local approximation:')
d=d.replace('and the post hoc evaluation of the unused \\qty{47.9}{\\cm} test improves as well (Appendix~\\ref{app:model_geometry}).', 'on the adopted fold. The excluded-test bookkeeping is retained in Appendix~\\ref{app:model_geometry}.')
# Remove an obsolete exception clause from the historical extensor narrative,
# not the historical results or the appendix's exclusion record.
d=re.sub(r', with the exception of maximum absolute residual on the \$l_\{0\}=\\qty\{48\.0\}\{\\cm\}\$ BPA', '',d)

# Correct the empirical loss interpretation in both documents.
for name in ['p','d']:
    t=locals()[name]
    t=t.replace('exactly where a wrapped actuator has the least contraction to spare', 'where the unbent actuator would otherwise produce its largest force')
    t=t.replace('so the wrapping loss grows quadratically as the actuator approaches its resting length and vanishes at full contraction.', 'so the fitted wrapping loss grows quadratically toward resting length and vanishes at the normalized free-contraction limit.')
    locals()[name]=t

# The latest whole-leg CSV compares configurations on the same underlying
# master's routing; it does not establish an attachment-geometry improvement.
survey=r'''The capability survey grades each actuator against the same-named single muscle of the Gait2392 model on a $5\times5$ grid of two joint coordinates. The hip adduction domain is clipped by its flexion-dependent envelope. A cell meets the requirement when the BPA is taut and its torque has the human torque's sign and at least its magnitude. The coverage metric is the fraction of included grid cells meeting that condition; the aggregate is the unweighted mean of the per-actuator fractions, not a continuous-area integral or a whole-leg torque sum.

The comparison includes all 27 actuator sites. The baseline uses the master's-era attachment routing at unity transmission; its eight nominal \qty{40}{\mm} sites are evaluated with four parallel \qty{20}{\mm} BPAs, while the other sites use the survey's base configuration. The updated side selects the saved transmission-ladder winner at each site, with actuator length and tendon length re-solved per configuration. It retains the same underlying attachment routing and the original diameter specification, including \qty{40}{\mm} actuators at those eight sites. The comparison therefore combines transmission, actuator sizing, and the baseline diameter substitution; it does not isolate an effect of attachment placement or compare two all-\qty{20}{\mm} assemblies.

The mean grid coverage increases from \qty{37.28}{\percent} to \qty{65.19}{\percent}, improving at 23 sites and tying at four (Table~\ref{tab:CapabilitySurvey}). A separate sign-aware torque-ratio score averages $\min(M_R/M_H,1)$ over the graded cells, assigning zero to inactive cells. Its mean falls from $-1.93$ to $-3.63$. Negative ratios penalize torque in the wrong direction, so greater coverage can coexist with a worse mean ratio where severe wrong-sign predictions remain. These two metrics support different conclusions and should be considered together.

'''
st=p.index('The knee redesigns are graded one joint at a time')
en=p.index('\\begin{table}[H]',st)
p=p[:st]+survey+p[en:]
p=p.replace('Whole-body capability of the hand-tuned and optimized BPA routings against single-muscle human torque surfaces, mean per-actuator share of the graded domain over the 19 muscle-actuators graded on both routings.', 'Capability across 27 actuator sites against single-muscle human torque grids. Means are taken over sites. The baseline substitutes four parallel \\qty{20}{\\mm} BPAs at each nominal \\qty{40}{\\mm} site; the transmission-ladder winners retain the original diameter specifications. Negative torque-ratio scores reflect wrong-sign torque, not a torque deficit in newton meters.')
p=p.replace('\\textbf{Active force (\\%)}', '\\textbf{Mean torque-ratio score}')
p=p.replace('Hand-tuned & 66.32 & 97.05', 'Baseline & 37.28 & $-1.93$')
p=p.replace('Optimized & 75.30 & 95.23', 'Transmission-ladder winners & 65.19 & $-3.63$')
st=d.index('The placement optimization began on the rigid-body model',d.index('\\subsection{Masters-Era'))
en=d.index('\\subsection{Physical Tests',st)
d=d[:st]+survey.replace('(Table~\\ref{tab:CapabilitySurvey})','').replace('The capability survey','The configuration survey')+d[en:]
st=d.index('The extensor requirement is set by vastus medialis',d.index('\\subsection{Minimum-Actuator'))
en=d.index('\\subsection{Masters-Era',st)
d=d[:st]+r'''The vastus medialis study evaluates one to three parallel \qty{20}{\mm} BPAs and one to three tackle line parts at \qty{620}{\kPa}, with actuator and tendon lengths re-solved for each combination. The attachment geometry is fixed by the study's routing family. Both OpenSim and cable-path moment arms are tested, with and without the compliance correction. In the corrected mode the calculation uses the pulley model directly; this model does not include the wrapping-loss term $\chi_3$, so it is not the same prediction as the cam redesign.

No tested configuration meets the target at all five graded angles. The best coverage is four of five, obtained with three BPAs through a two-part tackle using the compliance correction and the cable-path moment arm. Its four passing-angle margins are $+56.3$, $+40.1$, $+35.4$, and $+1.0$~\unit{\newton\meter}. It fails at \ang{10}, where this routing family approaches zero torque. The cable path uses a ring point \qty{0.75}{\mm} off the bone surface. These results identify a limitation of the tested routing family; they do not establish a minimum actuator count capable of meeting the complete envelope. The separate cam redesign meets its own target over its design grid with two BPAs (Section~\ref{sec:route_redesign_results}).

'''+d[en:]
d=d.replace('\\subsection{Minimum-Actuator Count for the Vastus Medialis Requirement}', '\\subsection{Actuator-Count and Transmission Study for Vastus Medialis}')
d=d.replace('\\subsection{Masters-Era and Updated Design Numbers}', '\\subsection{Whole-Leg Transmission Capability}')

# Fix conspicuous journal-register defects and unsupported concluding rhetoric.
p=p.replace('the cam converts the razor-thin margin', 'the cam increases the minimum margin')
p=p.replace('The cam converts the razor-thin margin of the direct route into genuine headroom.', 'The cam increases the torque margin relative to the direct route.')
p=p.replace('The geometry confirms that the headroom is not bought with routing risk:', 'The geometric checks give')
p=p.replace('That procedure is inexpensive next to the test campaign it saves, which is the point of having it.', 'Transfer to other systems requires validation against their measured torque.')
p=p.replace('Additional testing such as isokinetic and isotonic that are frequently used in in human dynamometer testing may yield equally as important insights.', 'Isokinetic and isotonic tests would assess behavior outside the isometric conditions evaluated here.')
st=p.index('This work started with an improved braided pneumatic actuator force-length-pressure model')
en=p.index('\\section*{Resource Identification Initiative}',st)
p=p[:st]+r'''Isometric knee-torque measurements on pinned and biomimetic test stands showed that actuator force characterization alone does not determine delivered joint torque. A hybrid calculation using measured muscle length and moment arm helped separate force-model error from routing and attachment effects. The corrected torque model represents a constant length offset, effective attachment compliance, and an empirical wrapping loss. A two-stage identification shares one stiffness pair across flexor and extensor configurations and reduces prediction error on the held-out tests.

Route redesign under the corrected model produces two-\qty{20}{\mm}-BPA configurations that meet the selected human-muscle torque targets over their design grids, with the required \qty{5}{\percent} margin. The cam extensor increases predicted torque headroom while reducing the required operating pressure. Physical tests meet or exceed the respective human targets at all 15 measured flexor points and all nine measured extensor points. These measurements validate the tested isometric configurations; dynamic loading and transfer to other joints remain to be evaluated.

'''+p[en:]
# Preserve old appendix figures and exclusion bookkeeping, but fix column-vector
# projection notation and identify historical table values as historical.
a=a.replace('c_{b} = \\hat{\\mathbf{u}}\\,\\mathbf{K}_{br}^{-1}\\hat{\\mathbf{u}}^{\\top}', 'c_{b} = \\hat{\\mathbf{u}}^{\\top}\\mathbf{K}_{br}^{-1}\\hat{\\mathbf{u}}')

for path,after in zip(paths,[d,a,p,s]):
    path.write_text(after,encoding='utf-8',newline='\n')
    diff=''.join(difflib.unified_diff(before[path].splitlines(True),after.splitlines(True),fromfile=str(path)+' (session start)',tofile=str(path)))
    (OUT/(path.stem+'.diff')).write_text(diff,encoding='utf-8')
    print(path.name, len(before[path]),'->',len(after))
    oldlabels=set(re.findall(r'\\label\{([^}]+)\}',before[path]))
    newlabels=set(re.findall(r'\\label\{([^}]+)\}',after))
    print('labels removed', sorted(oldlabels-newlabels),'added',sorted(newlabels-oldlabels))
print('main display items',len(re.findall(r'\\begin\{(?:figure|table|sidewaystable)\}',p)))
