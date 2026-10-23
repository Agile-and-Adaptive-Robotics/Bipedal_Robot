from pathlib import Path
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
p=C/'chapters/70-overall.tex';s=p.read_text(encoding='utf-8')
s=s.replace('both are documented as complete in the chapter\'s methods (Section~\\ref{sec:preliminary_sim_methods})','the implemented interfaces and their limits are described in Section~\\ref{sec:preliminary_sim_methods}')
s=s.replace('production walking configurations have been tuned against a walking reference and evaluated on the ground','supported stepping configurations have been tuned against a walking reference and evaluated on the ground')
s=s.replace('The first of these tests is complete at the rhythm-generator level, where tonic drive evokes an approximately \\qty{0.9}{\\second} antiphase cycle in the isolated network (Section~\\ref{sec:preliminary_sim_methods});',r'An isolated rhythm-generator test produced an approximately \qty{0.9}{\second} antiphase cycle under tonic drive (Appendix~\ref{app:simulation_verification});')
s=s.replace("The author's aspiration for this program is",'My aspiration for this program is')
s=s.replace('excititory','excitatory')
p.write_text(s,encoding='utf-8')
p=C/'chapters/02-abstract.tex';(C/'review/before/02-abstract.tex').write_bytes(p.read_bytes());s=p.read_text(encoding='utf-8-sig')
old='The biomimetic knee driven by optimized BPA placements was shown to meet or exceed human isometric knee torque over most of the range of motion.'
# Replace only the active abstract; preserve existing commented historical drafts.
idx=s.index(old);s=s[:idx]+'The corrected torque model supports route redesigns predicted to meet human-scaled knee-torque targets, with full experimental validation of the redesigned envelopes still pending.'+s[idx+len(old):]
p.write_text(s,encoding='utf-8')
print('Aligned abstract claim and remaining discussion references')
