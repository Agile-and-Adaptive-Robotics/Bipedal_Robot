from pathlib import Path
import fitz
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');S=C/'review/user_sweep_20261023'
for rel in ['15-background/steeleknee.pdf','30-torque/bktFrame.pdf']:
 d=fitz.open(C/'figs'/rel);p=d[0];print(rel,repr(p.get_text())[:1800]);p.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(S/(Path(rel).stem+'.png'))
d=fitz.open(r'C:\Users\Ben Bolen\Zotero\storage\J8K862TH\SfN24_Bolen.pdf');d[0].get_pixmap(matrix=fitz.Matrix(.6,.6)).save(S/'poster.png')
