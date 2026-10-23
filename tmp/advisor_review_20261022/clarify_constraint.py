from pathlib import Path
import json
H=Path(__file__).resolve().parent;R=H.parent.parent
for p in [R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/70-overall.tex',H/'overleaf_expected/70-overall.tex']:
    b=p.read_bytes();assert b.count(b'impose axial symmetry')==1;p.write_bytes(b.replace(b'impose axial symmetry',b'fix axial displacement',1))
(H/'constraint_patch.json').write_text(json.dumps({'file':'70-overall.tex','old':'impose axial symmetry','new':'fix axial displacement'}))
