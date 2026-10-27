from pathlib import Path
import re
R=Path(r'D:\GitHub\Bipedal_Robot')
D=R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex'
A=D.with_name('94-AppendixC.tex')
P=R/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex'
S=P.with_name('frontiers_SupplementaryMaterial.tex')
d,a,p,s=[x.read_text(encoding='utf-8') for x in [D,A,P,S]]
# Keep the optimization-weight table within the dissertation margins.
assert r'\begin{tabular}{@{}llp{8.2cm}@{}}' in d
d=d.replace(r'\begin{tabular}{@{}llp{8.2cm}@{}}',r'\begin{tabularx}{\textwidth}{@{}p{0.28\textwidth}p{0.17\textwidth}L@{}}',1)
st=d.index(r'\begin{tabularx}{\textwidth}{@{}p{0.28')
en=d.index(r'\end{tabular}',st)
d=d[:en]+d[en:].replace(r'\end{tabular}',r'\end{tabularx}',1)
# Remove an Appendix-C version of the unity-gain causal claim Ben had
# already rejected in the torque paper. Preserve the numerical run outcomes.
st=a.index('What the two algorithms buy is therefore reach and screening:')
en=a.index('\n\n',st)
a=a[:st]+'These clearance results motivate evaluating attachment sites with more space for the reeved routes. The pulley solutions remain optimizer-stage results and have not been propagated into the OpenSim or MuJoCo model trees. Physical isometric tests of the direct redesigned flexor and cam extensor are reported in Section~\\ref{sec:redesign_physical}.'+a[en:]
# Avoid a long inline force-balance expression in the appendix paragraph.
a=a.replace(r'the force balance $F_{\max}F^{*}(\varepsilon^{*}(\delta), P^{*}) = k_{\text{eff}}\,\delta$ is solved',r'the force balance of Equation~\eqref{eq:complianceforcebalance} is solved')
# Supplementary methods need the release calculation as well as its output
# schedule, as requested by Ben. Port the same derivation, changing only
# the section level and the worked-example reference.
st=a.index(r'\section{Contact and Release Angles of the Wrapped Route}')
en=a.index(r'\section{Supporting Route Searches',st)
release=a[st:en].replace(r'\section{',r'\subsection{',1)
release=release.replace(r'Table~\ref{tab:app_ext_seed}',r'Table~\ref{table:route_redesign}')
# Work out the actual local label instead of inventing a table reference.
blocks=list(re.finditer(r'\\begin\{table\}.*?\\end\{table\}',s,re.S))
route=next(m.group() for m in blocks if 'Routing via points of the \\qty{20}' in m.group())
lab=re.search(r'\\label\{([^}]+)\}',route).group(1)
release=release.replace(r'Table~\ref{table:route_redesign}',r'Table~\ref{'+lab+'}')
# The fixed 9-point representation is shared, but distinguish the appendix's
# direct-route schedule from the supplementary table's cam schedule.
release=release.replace('is the worked example, and the release schedule the calculation produces is the one quoted there.', 'gives the cam-route example; its release schedule follows the same rule on the cam geometry.')
s=s.replace(r'\end{document}',release+'\n'+r'\end{document}')
s=s.replace('The flexor design routes two BPAs on mirrored paths:', 'The flexor design uses a staggered pair of BPA paths:')
p=p.replace('its contact-release schedule is', 'its contact-release schedule is')
# Add an explicit pointer from the main-paper route method to the derivation.
p=p.replace('Three variants extend the search space.', 'The contact-release calculation is given in the Supplementary Material (Section~\\ref{app:contact_release}). Three variants extend the search space.')
p=p.replace('nine measured data measurements','nine measurements')
for x,t in zip([D,A,P,S],[d,a,p,s]):x.write_text(t,encoding='utf-8',newline='\n')
print('Restored supplementary contact-release derivation; label',lab)
