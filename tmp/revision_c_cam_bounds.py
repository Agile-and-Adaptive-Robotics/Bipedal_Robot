from pathlib import Path
import fitz
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
d=fitz.open(C/'figs/30-torque/extensor_cam_route.pdf')
for b in d[0].get_text('blocks'):
    if (410<b[1]<475) or (190<b[1]<250):print(tuple(round(x,2) for x in b[:4]),b[4].encode('ascii','backslashreplace').decode().replace('\n',' / '))
