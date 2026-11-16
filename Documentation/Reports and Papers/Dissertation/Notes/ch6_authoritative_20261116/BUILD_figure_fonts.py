"""Read PDF text spans at declared TeX placement dimensions; no PDF edits."""
from pathlib import Path
import json, re
import fitz

ROOT=Path(__file__).resolve().parents[5]
DISS=ROOT/'Documentation/Reports and Papers/Dissertation'
NOTES=Path(__file__).resolve().parent
report=[]
for name in ('60-sns','95-AppendixD'):
    text=(DISS/'Final/chapters'/f'{name}.tex').read_text(encoding='utf-8')
    for options,file in re.findall(r'\\includegraphics(?:\[([^]]*)\])?\{([^}]+)\}',text):
        path=DISS/'Figures'/name/file
        if path.suffix.lower()!='.pdf': continue
        if not path.exists():
            report.append({'path':str(path),'placement_options':options,'missing_during_edit':True})
            continue
        doc=fitz.open(path)
        page=doc[0]
        # psutex sets width=6.5in; standard 12pt letterpaper book computes
        # textheight=floor(7.5in/14.5pt)*14.5pt+12pt = 548.5 TeX points.
        values={'\\textwidth':6.5*72,'\\textheight':548.5*72/72.27}
        def dimension(expr):
            expr=expr.strip()
            for key,value in values.items():
                if expr.endswith(key):
                    prefix=expr[:-len(key)]
                    return float(prefix or '1')*value
            m=re.fullmatch(r'([0-9.]+)(cm|mm|in|pt)',expr)
            if m: return float(m[1])*{'cm':72/2.54,'mm':72/25.4,'in':72,'pt':72/72.27}[m[2]]
            return None
        scales=[]
        for axis,expr in re.findall(r'(width|height)=([^,]+)',options):
            amount=dimension(expr)
            if amount: scales.append(amount/(page.rect.width if axis=='width' else page.rect.height))
        scale=min(scales) if scales else None
        spans=[span for block in page.get_text('dict')['blocks'] if 'lines' in block
            for line in block['lines'] for span in line['spans'] if span['text'].strip()]
        fonts=sorted({round(span['size'],3) for span in spans})
        small=[{'text':span['text'],'source_pt':span['size'],'placed_pt':span['size']*scale}
               for span in spans if scale and span['size']*scale<9.999]
        report.append({'path':str(path),'placement_options':options,
                       'natural_size_bp':[page.rect.width,page.rect.height],
                       'assumed_placement_scale':scale,'source_font_sizes_pt':fonts,
                       'below_10pt_spans':small,
                       'limit':'Source-scale inference, not compiled-page verification. Math superscripts/subscripts require human review.'})
        doc.close()
(NOTES/'BUILD_FIGURE_FONT_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
for item in report:
    if item.get('missing_during_edit'):
        print(Path(item['path']).name,'not yet available during editing')
        continue
    print(Path(item['path']).name, 'scale=',item['assumed_placement_scale'],
          'min source pt=',min(item['source_font_sizes_pt']) if item['source_font_sizes_pt'] else None,
          'small spans=',len(item['below_10pt_spans']))
