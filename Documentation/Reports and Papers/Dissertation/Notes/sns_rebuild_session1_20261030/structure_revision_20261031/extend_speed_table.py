from pathlib import Path
import json
p=Path(r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\structure_revision_20261031')
s=(p/'60-sns.tex').read_text(encoding='utf-8')
b=s.rfind('\\begin{table}',0,s.index('\\label{tab:sns_runtime}'));e=s.index('\\end{table}',b)+len('\\end{table}')
s=s[:b]+r'''\begin{table}[htbp]
\centering
\caption[Measured neural and coupled-model acceleration]{Execution times on easteregg2 after warm-up. The Gait2392 reference uses native SNS-Toolbox; W2L uses the same ported equations in Python. Neural-only rows exclude mechanical integration. The equation-kernel row compares a Di Russo-inspired Python implementation with its compiled Numba implementation.}
\label{tab:sns_runtime}
\small
\begin{tabularx}{\textwidth}{@{}Xrrrr@{}}
\toprule
Computation & Workload & Reference (s) & Fast (s) & Speedup \\
\midrule
\multicolumn{5}{@{}l}{Coupled MuJoCo models} \\
Gait2392 balance & 8 s & 122.120 & 13.253 & $9.21\times$ \\
W2L80 suspended & 4 s & 17.208 & 0.970 & $17.75\times$ \\
\addlinespace[5pt]
\multicolumn{5}{@{}l}{Neural calculations only} \\
SNS-Toolbox, 1,027 neurons & 500 updates & 24.534 & 0.245 & $100.17\times$ \\
Same benchmark, quiet repeat & 500 updates & 14.306 & 0.253 & $56.55\times$ \\
Python equation kernel & 1 update & $1.978\times10^{-3}$ & $3.385\times10^{-6}$ & $584.33\times$ \\
\bottomrule
\end{tabularx}
\end{table}'''+s[e:]
s=s.replace('Table~\\ref{tab:sns_runtime} gives matched coupled-loop timings after warm-up, excluding loading, compilation, file output, and animation rendering.','Table~\\ref{tab:sns_runtime} distinguishes complete coupled-loop timings from neural-only timings after warm-up. Loading, compilation, file output, and animation rendering are excluded.')
s=s.replace('The saved Gait2392 positions, controls, center-of-mass, and neural arrays were identical between backends.','The 1,027-neuron benchmark accelerated fixed-conductance neural updates by 100.17 times, with a maximum voltage discrepancy of \\qty{3.55e-15}{\\milli\\volt}. A quiet-machine repeat measured 56.55 times acceleration. The separately compiled Di Russo-inspired neural equations reached 584.33 times acceleration over their Python implementation. Mechanical integration, contact, and muscle calculations account for the smaller speedups of complete coupled simulations.\n\nThe saved Gait2392 positions, controls, center-of-mass, and neural arrays were identical between backends.')
(p/'60-sns.tex').write_text(s,encoding='utf-8')
print('Speed table extended with verified high-ratio neural measurements.')
