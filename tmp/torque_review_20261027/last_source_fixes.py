from pathlib import Path
R=Path(r'D:\GitHub\Bipedal_Robot')
d=R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex'
p=R/'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex'
t=d.read_text(encoding='utf-8')
t=t.replace('The \\qty{47.9}{\\cm} flexor test was not used in identification and is evaluated separately in Appendix~\\ref{app:model_geometry}. ', '')
t=t.replace('The coefficient values reported with the results figures belong to the initial single-bracket identification. Re-identifying the correction terms over the full test sets with the two-stage procedure of Section~\\ref{sec:identification} yields the adopted values', 'The results figures use the adopted two-stage identification. The procedure of Section~\\ref{sec:identification} yields')
if t!=d.read_text(encoding='utf-8'):
    d.write_text(t,encoding='utf-8',newline='\n')
t=p.read_text(encoding='utf-8')
t=t.replace('The flexor redesign was a two-BPA mirrored-route solution rather than a single shared route.', 'The flexor redesign uses two staggered BPA routes rather than a single shared route.')
t=t.replace('solved from it by a mirrored-endpoint rule:', 'derived from it by an offset-preserving endpoint rule:')
if t!=p.read_text(encoding='utf-8'):
    p.write_text(t,encoding='utf-8',newline='\n')
