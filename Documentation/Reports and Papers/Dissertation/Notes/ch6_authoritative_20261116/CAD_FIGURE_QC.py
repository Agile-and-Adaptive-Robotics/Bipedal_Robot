"""Read-only vector font/scale audit for the active Chapter 6/Appendix D sources."""
from pathlib import Path
import re,json
import fitz

HERE=Path(__file__).resolve().parent
DISS=HERE.parent.parent
# book12pt letterpaper calculation: floor((11-3.5)in/14.5pt)*14.5pt+12pt.
TEXT_WIDTH_BP=6.5*72
TEXT_HEIGHT_TEX_PT=int((7.5*72.27)/14.5)*14.5+12
TEXT_HEIGHT_BP=TEXT_HEIGHT_TEX_PT*72/72.27

def dim(text):
    text=text.strip()
    if text.endswith('cm'):return float(text[:-2])*72/2.54
    if text.endswith('in'):return float(text[:-2])*72
    if text.endswith('pt'):return float(text[:-2])*72/72.27
    for unit,size in [('\\textwidth',TEXT_WIDTH_BP),('\\textheight',TEXT_HEIGHT_BP)]:
        if text.endswith(unit):return float(text[:-len(unit)]or'1')*size
    return None

rows=[]
for chapter in ['60-sns','95-AppendixD']:
    src=DISS/'Final/chapters'/f'{chapter}.tex'
    text=src.read_text()
    for match in re.finditer(r'\\includegraphics\[([^]]+)\]\{([^}]+)\}',text):
        opts,name=match.groups();path=DISS/'Figures'/chapter/name
        row={'chapter':chapter,'line':text[:match.start()].count('\n')+1,'file':str(path),'options':opts,'exists':path.exists()}
        rows.append(row)
        if path.suffix.lower()!='.pdf' or not path.exists():continue
        doc=fitz.open(path);p=doc[0];w,h=p.rect.width,p.rect.height
        ww=re.search(r'(?:^|,)width=([^,]+)',opts);hh=re.search(r'(?:^|,)height=([^,]+)',opts)
        tw=dim(ww[1])if ww else None;th=dim(hh[1])if hh else None
        factors=[x for x in [tw/w if tw else None,th/h if th else None]if x is not None]
        factor=min(factors)if factors else 1
        spans=[s for b in p.get_text('dict')['blocks']if 'lines'in b for line in b['lines']for s in line['spans']if s['text'].strip()]
        sizes=[s['size']for s in spans]
        row.update({'source_size_cm':[w*2.54/72,h*2.54/72],'applied_scale':factor,'rendered_size_cm':[w*factor*2.54/72,h*factor*2.54/72],
                    'source_min_font_bp':min(sizes)if sizes else None,'rendered_min_font_bp':min(sizes)*factor if sizes else None,
                    'small_rendered_text_examples':[{'text':s['text'],'font_bp':s['size']*factor}for s in spans if s['size']*factor<9.9][:15],
                    'unambiguous_10pt_floor':bool(sizes)and min(sizes)*factor>=9.9,
                    'available_caption_height_cm':(TEXT_HEIGHT_BP-h*factor)*2.54/72})
        pix=p.get_pixmap(matrix=fitz.Matrix(1.4*factor,1.4*factor),alpha=False)
        render=DISS.parents[2]/'tmp/ch6_authoritative_20261116/CAD_figure_renders';render.mkdir(parents=True,exist_ok=True)
        pix.save(render/(chapter+'_'+path.stem+'.png'))
staged=[]
for path in [HERE/'CAD_supported_traces_20261116.pdf',*sorted((HERE/'CAD_sns_compact').glob('*.pdf'))]:
    p=fitz.open(path)[0]
    spans=[s for b in p.get_text('dict')['blocks']if 'lines'in b for line in b['lines']for s in line['spans']if s['text'].strip()]
    staged.append({'file':str(path),'size_cm':[p.rect.width*2.54/72,p.rect.height*2.54/72],'minimum_font_bp':min(s['size']for s in spans)})
out={'textwidth_bp':TEXT_WIDTH_BP,'textheight_tex_pt':TEXT_HEIGHT_TEX_PT,'textheight_bp':TEXT_HEIGHT_BP,'dimension_source':'Canonical psutex.cls delegates to installed book.cls/bk12.clo; class overrides textwidth6.5in but not textheight. These are source-derived dimensions, not an observed successful build.','rows':rows,'staged':staged}
(HERE/'CAD_FIGURE_QC.json').write_text(json.dumps(out,indent=2)+'\n')
for row in rows:
    if 'rendered_min_font_bp'in row:
        print(Path(row['file']).name,'min=',round(row['rendered_min_font_bp']or 0,2),'height=',round(row['rendered_size_cm'][1],2),'caption_room_cm=',round(row['available_caption_height_cm'],2))
