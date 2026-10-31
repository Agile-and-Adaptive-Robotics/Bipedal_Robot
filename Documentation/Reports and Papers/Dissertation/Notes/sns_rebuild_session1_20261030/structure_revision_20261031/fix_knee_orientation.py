from pathlib import Path
import pymupdf as f
base=Path(r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation')
p=base/'Figures/60-sns/knee_transient_current_layout.pdf';old=f.open(p)
out=base/'Notes/sns_rebuild_session1_20261030/structure_revision_20261031'
new=f.open();page=new.new_page(width=470,height=600)
for r in [f.Rect(0,0,470,433),f.Rect(0,573,470,600)]:page.show_pdf_page(r,old,0,clip=r)
page.show_pdf_page(f.Rect(227,433,437,573),old,0,clip=f.Rect(20,433,230,573))
page.show_pdf_page(f.Rect(23,433,233,573),old,0,clip=f.Rect(230,433,440,573))
new.save(out/'knee_transient_orientation.pdf');new[0].get_pixmap(matrix=f.Matrix(1,1)).save(out/'knee_orientation_preview.png')
print('Panel C corrected: flexor left, extensor right; connections and results preserved.')
