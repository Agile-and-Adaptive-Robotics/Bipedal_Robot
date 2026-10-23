from pathlib import Path
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
p=C/'chapters/60-sns.tex';s=p.read_text();start=s.index('I subsequently expanded');end=s.index('\n\n\\begin{figure}',start)
s=s[:start]+r'''I subsequently expanded the neural model into left- and right-side circuits with interconnected rhythm generators (RGs). I constructed the right-side subsystems from the left-side organization and linked the two RGs through commissural pathways. Further circuit schematics and implementation records are provided in Appendix~\ref{app:simulation_verification}.

Air stepping did not establish whether the bilateral coupling and sensory feedback could coordinate the legs under load. Manual tuning could not separate controller-design limitations from implementation errors or problems in mechanical coupling. Repeatable tests of the circuits, actuator interface, and closed-loop model were needed to resolve those questions.'''+s[end:];p.write_text(s,encoding='utf-8')
