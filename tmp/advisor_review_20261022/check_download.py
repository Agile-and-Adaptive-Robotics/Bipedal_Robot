from pathlib import Path
import fitz
H=Path(__file__).resolve().parent
p=Path(r'C:\Users\Ben Bolen\Downloads\1332d1d2-1b69-4eeb-a3a5-89ad039ebcc7.tmp')
data=p.read_bytes();print('PDF_HEADER',data[:8],'EOF',data[-30:])
if data.startswith(b'%PDF-') and b'%%EOF' in data[-100:]:
    doc=fitz.open(p);print('PAGES',len(doc))
    for num in ['6.6','7.2']:
        hits=[i for i,page in enumerate(doc) if 'Figure '+num+':' in page.get_text()]
        print(num,'pages',[i+1 for i in hits])
        for i in hits:
            page=doc[i];page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(H/f'overleaf_figure_{num}.png')
            spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
            print('YMAX',max(s['bbox'][3] for s in spans),'PAGEHEIGHT',page.rect.height)
