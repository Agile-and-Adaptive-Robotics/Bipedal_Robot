from pathlib import Path
import fitz
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
for name in ['PulleyTorqueSurfaces_20261020-eps-converted-to','PulleyMarginSurfaces_20261020-eps-converted-to']:
    p=C/'figs/30-torque'/f'{name}.pdf';d=fitz.open(p);d[0].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(C/'review'/f'{name}.png');print(name,d[0].get_text()[:1500].encode('ascii','backslashreplace').decode())
