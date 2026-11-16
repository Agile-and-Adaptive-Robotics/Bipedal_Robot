"""Measure actual compiled sans-serif figure spans, with raster limits stated."""
from pathlib import Path
import argparse, hashlib, json, re
import fitz

NOTES=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--revision',default='r6')
args=parser.parse_args()
record=json.loads((NOTES/f'BUILD_RENDER_{args.revision}.json').read_text())
doc=fitz.open(record['pdf'])
if hashlib.sha256(Path(record['pdf']).read_bytes()).hexdigest()!=record['pdf_sha256']:
    raise SystemExit('PDF hash mismatch')
results=[]
for target in record['pages']:
    page=doc[target['pdf_page_1based']-1]
    spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']
        if s['text'].strip() and re.search(r'Arial|DejaVu|LiberationSans|Helvetica',s['font'],re.I)]
    if not spans: continue
    results.append({'pdf_page':target['pdf_page_1based'],'printed_label':target['printed_label'],
        'span_count':len(spans),'minimum_figure_text_bp':min(s['size'] for s in spans),
        'below_9p98bp_spans':[{'text':s['text'],'font':s['font'],'size_bp':s['size'],'bbox':list(s['bbox'])}
                              for s in spans if s['size']<9.98]})
report={'pdf':record['pdf'],'pdf_sha256':record['pdf_sha256'],
    'method':'Actual transformed PDF text spans in Arial/DejaVu/LiberationSans/Helvetica, cross-checked with source figure font records.',
    'raster_limit':'Raster text has no recoverable PDF font size; assessed visually at publication size, not claimed as a numeric font-size measurement.',
    'minimum_detected_figure_text_bp':min(item['minimum_figure_text_bp'] for item in results),
    'below_floor_count':sum(len(item['below_9p98bp_spans']) for item in results),'pages':results}
(NOTES/f'BUILD_ACTUAL_FIGURE_FONT_CHECK_{args.revision}.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='pages'},indent=2))
