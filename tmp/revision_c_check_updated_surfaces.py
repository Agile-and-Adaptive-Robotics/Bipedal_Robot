from pathlib import Path
import fitz
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Knee_Torque_Test/Figures'
for name in ['PulleyTorqueSurfaces_20261020.pdf','PulleyMarginSurfaces_20261020.pdf']:
    d=fitz.open(D/name);s=d[0].get_text();print(name,s[-1050:].encode('ascii','backslashreplace').decode())
