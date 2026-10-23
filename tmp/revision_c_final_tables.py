from pathlib import Path
import re
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/final_sweep_20261023'
p=C/'chapters/30-torque.tex';s=p.read_text();a_path=C/'chapters/94-AppendixC.tex';a=a_path.read_text()
def table(label):
    for m in re.finditer(r'\\begin\{table\}.*?\\end\{table\}',s,re.S):
        if r'\label{'+label+'}' in m.group():return m.group()
    raise ValueError(label)
moved=[]
for label in ['tab:FlxPinPreGof','tab:FlxPinGof','tab:ExtPinGof','tab:TorqueFit']:
    old=table(label);b=old
    if label=='tab:FlxPinPreGof':
        b=b.replace(r'\begin{table}[htbp]',r'\begin{sidewaystable}[p]'+'\n'+r'\centering\small\setlength{\tabcolsep}{4pt}').replace(r'\end{table}',r'\end{sidewaystable}')
    else:
        b=b.replace(r'\begin{table}[htbp]',r'\begin{table}[p]'+'\n'+r'\centering\small')
        if label=='tab:TorqueFit':
            b=b.replace(r'\begin{tabular}{@{}lccc@{}}',r'\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.55\textwidth}*{3}{>{\centering\arraybackslash}X}@{}}')
        else:
            b=b.replace(r'\begin{tabular}{@{}lcccccc@{}}',r'\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}}')
        b=b.replace(r'\end{tabular}',r'\end{tabularx}')
    moved.append(b)
    replacement=''
    if label=='tab:ExtPinGof':
        replacement=r'''Table~\ref{tab:torque_fit_summary} summarizes the training and validation means at the adopted terms. The full per-test comparisons, including maximum residuals and the earlier force-model comparison, are retained in Appendix~\ref{sec:app_torque_fit_details}.

\begin{table}[htbp]
\centering\small
\caption{Mean pinned-knee errors before and after transmission correction. RMSE is in N\,m; FVU is dimensionless. Per-test values and fold assignments are in Appendix~\ref{sec:app_torque_fit_details}.}
\label{tab:torque_fit_summary}
\begin{tabular}{@{}llrrrr@{}}
\toprule
 & & \multicolumn{2}{c}{RMSE} & \multicolumn{2}{c}{FVU}\\
\cmidrule(lr){3-4}\cmidrule(l){5-6}
Route & Set & Before & Corrected & Before & Corrected\\
\midrule
Flexor & Training & 8.656 & 1.624 & 2.1213 & 0.0998\\
Flexor & Validation & 7.588 & 1.490 & 1.2657 & 0.0447\\
Extensor & Training & 2.498 & 0.997 & 2.2502 & 0.5146\\
Extensor & Validation & 2.554 & 1.098 & 2.4236 & 0.3965\\
\bottomrule
\end{tabular}
\end{table}'''
    s=s.replace(old,replacement)
old=table('tab:AdoptedXi')
new=r'''\begin{table}[htbp]
\centering\small
\caption{Adopted transmission corrections. The stiffness pair is shared across flexor and extensor configurations; parameter bounds and full provenance are in Appendix~\ref{app:model_geometry}.}
\label{tab:AdoptedXi}
\begin{tabularx}{\textwidth}{@{}ll>{\raggedright\arraybackslash}X@{}}
\toprule
Parameter & Value & Configuration and role\\
\midrule
$\chi_{0}$ & \qty{+3.94}{\mm} & Flexor effective-length offset\\
$\chi_{1}$ & \qty{3.998e4}{\N\per\m} & Shared effective stiffness\\
$\chi_{2}$ & \qty{1.473e4}{\N\per\m} & Shared effective stiffness\\
$\chi_{0}$ & \qty{-6.38}{\mm} & Extensor effective-length offset\\
$\chi_{3}$ & \num{0.158} & Dimensionless wrapping-loss scale\\
\bottomrule
\end{tabularx}
\end{table}'''
s=s.replace(old,new)
start=s.index('The extensor and flexor isometric measurements were reduced');end=s.find('\n\n',start)
if end<0:end=len(s)
physical=r'''\subsection{Physical Tests of the Redesigned Routes}\label{sec:redesign_physical}

The redesigned cam-extensor route was tested at nine isometric holds, and the two-BPA flexor route at fifteen holds. Both test sets span knee angles from \ang{-119} to \ang{5}. The test-specific analysis scripts transform the measured load-cell wrench to the moving knee center and evaluate model predictions at each measured angle and pressure. Figures~\ref{fig:extensor_physical_final} and~\ref{fig:flexor_physical_final} compare these reductions with the human-reference torques and the optimized \qty{620}{\kPa} curves.

The measured torques meet or exceed the human reference over most of the tested range. Agreement with the pressure-specific prediction is uneven, particularly for the flexor at intermediate angles. The tests therefore support the redesigned routes' torque capacity while identifying prediction errors that remain to be explained. Since pressure varied between holds, the measured points do not constitute a full-pressure envelope at every angle.

\begin{figure}[p]
\centering
\includegraphics[width=0.97\textwidth]{extensor_20mm_physical.png}
\caption[Physical tests of the cam-extensor route]{Cam-extensor tests. Blue points are measured torques; yellow points are model predictions at the measured angles and pressures. The pink curve is the optimized \qty{620}{\kPa} prediction, and the black dotted curve is the human target. Point area scales with pressure.}
\label{fig:extensor_physical_final}
\end{figure}

\begin{figure}[p]
\centering
\includegraphics[width=0.97\textwidth]{flexor_20mm_physical.png}
\caption[Physical tests of the redesigned flexor route]{Two-BPA flexor tests with \qty{56}{\cm} actuators. Blue points are measured torques; yellow points are predictions at the measured angles and pressures. The dashed pink curve is the optimized \qty{620}{\kPa} prediction, and the black dotted curve is the human target. Point area scales with pressure.}
\label{fig:flexor_physical_final}
\end{figure}'''
s=s[:start]+physical+s[end:];p.write_text(s,encoding='utf-8')
a+='\n'+r'''\section{Supporting Torque-Prediction Fit Statistics}\label{sec:app_torque_fit_details}

These tables retain the per-test errors, fold assignments, and maximum residuals underlying Chapter~\ref{ch:torque}. The first compares three actuator-force predictions before transmission correction. The subsequent tables give the adopted corrected fits and the earlier cross-configuration results; these are distinct identification stages.

'''+ '\n\n'.join(moved)+'\n';a_path.write_text(a,encoding='utf-8')
p=C/'chapters/70-overall.tex';s=p.read_text().replace('; The monthly literature-scraping extension','. The monthly literature-scraping extension');p.write_text(s,encoding='utf-8')
(S/'table_moves.md').write_text('Moved tab:FlxPinPreGof, tab:FlxPinGof, tab:ExtPinGof, and tab:TorqueFit to Appendix C, retaining every numeric body entry and original label. Main Chapter 4 now has one four-row mean-error summary and the five-row adopted-correction table. Summary values are copied from the existing mean rows, not recalculated or fitted. The adopted table drops a duplicate bounds column; bounds and full axis assignments remain in existing Appendix C tables. Full original caption context is retained with the moved tables. The physical-test figures are unchanged copies of the committed MATLAB script exports, with input counts/ranges read from the scripts. No spreadsheet torque cells were used.\n',encoding='utf-8')
print('Moved four long tables, added compact summary, and inserted completed physical-test plots')
