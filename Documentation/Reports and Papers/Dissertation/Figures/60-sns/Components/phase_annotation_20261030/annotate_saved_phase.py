from pathlib import Path
import pymupdf as fitz
p=Path(__file__).resolve().parent
s=fitz.open(p/'w2l_feedback_before_annotation.pdf');o=fitz.open();q=o.new_page(width=s[0].rect.width,height=s[0].rect.height+48);q.show_pdf_page(fitz.Rect(0,0,s[0].rect.width,s[0].rect.height),s,0);y=s[0].rect.height
q.insert_text((35,y+18),'Median bilateral contact phase: 0.498 cycles',fontsize=15,color=(0.13,0.45,0.2));q.insert_text((460,y+18),'Median bilateral contact phase: 0.223 cycles',fontsize=15,color=(0.7,0.15,0.15));q.insert_text((168,y+39),'Acceptance interval: 0.35-0.65 cycles; right condition outside interval',fontsize=14)
o.save(p.parents[1]/'w2l_feedback_verified.pdf')
