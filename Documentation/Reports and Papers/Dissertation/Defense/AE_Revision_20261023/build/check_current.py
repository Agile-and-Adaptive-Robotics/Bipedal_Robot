from pathlib import Path
import fitz
r=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C\figs\30-torque')
o=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\AE_Revision_20261023\build')
for key,name in [('corrected_current','FlxPin_optimized-eps-converted-to.pdf'),('hybrid_current','FlxPin_preopt-eps-converted-to.pdf')]:
 d=fitz.open(r/name);d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(o/f'{key}.png')
