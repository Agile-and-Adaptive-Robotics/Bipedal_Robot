import pymupdf as f
from pathlib import Path
p=Path(r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\submission_live_20261031')
d=f.open(p/'Bolen_Dissertation_SESSION1.pdf')
print('Pages:',len(d))
need=['Figure 6.1:','Figure 6.2:','Matched coupled-loop','Gait2392 balance','hysteresis term','Figure D.4:','Native W2L']
for n in need:
 hits=[i+1 for i,x in enumerate(d) if n in x.get_text()]
 print(n,hits)
 for i in hits:
  if n in ['Figure 6.1:','Figure 6.2:','Gait2392 balance']:
   d[i-1].get_pixmap(matrix=f.Matrix(1.3,1.3)).save(p/('page_'+str(i)+'.png'))
for i in range(150,177):
 text=d[i].get_text()
 if 'Figure 6.' in text or 'Table 6.' in text:
  print(i+1, '\n'.join([x for x in text.splitlines() if 'Figure 6.' in x or 'Table 6.' in x]))
