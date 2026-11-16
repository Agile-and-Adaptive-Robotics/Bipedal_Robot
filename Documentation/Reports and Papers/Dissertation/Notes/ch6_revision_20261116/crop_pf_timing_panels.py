"""Select three unaltered vector panels from the recorded syn6 PF audit.

The full twelve-channel original remains in the portable experiment handoff.
This composition changes neither data coordinates nor phase alignment.
"""
from pathlib import Path
from copy import deepcopy
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import (
    NameObject, DictionaryObject, DecodedStreamObject, RectangleObject,
)

dissertation=Path(__file__).resolve().parents[2]
source=dissertation/'Defense/Session1_Handoff_20261104/experiments_20261104/circuit_trace_20261105/syn6/timing_audit/PF_TIMING_AUDIT.pdf'
target=dissertation/'Figures/95-AppendixD/syn6_pf_timing_selected_20261116.pdf'
original=PdfReader(source).pages[0]
scale=792/1870
width=960*scale
panel_height=308*scale
height=3*panel_height+78
writer=PdfWriter()
page=writer.add_blank_page(width=width,height=height)
for index,top_px in enumerate((73,391,1345)):
    p=deepcopy(original)
    bottom=864-(top_px+308)*scale
    p.cropbox=RectangleObject([0,bottom,width,bottom+panel_height])
    y=34+(2-index)*panel_height
    page.merge_transformed_page(p,Transformation().translate(0,y-bottom))

# Opaque annotation only replaces the small original legend and adds a common
# 10-point legend and x label. Existing axes, ticks, waveforms, and titles remain.
font=DictionaryObject({
    NameObject('/Type'):NameObject('/Font'),
    NameObject('/Subtype'):NameObject('/Type1'),
    NameObject('/BaseFont'):NameObject('/Helvetica'),
})
resources=page['/Resources']
fonts=resources.get('/Font')
if fonts is None:
    fonts=DictionaryObject()
    resources[NameObject('/Font')]=fonts
fonts[NameObject('/AnnotationFont')]=writer._add_object(font)
top=34+3*panel_height
stream=DecodedStreamObject()
stream.set_data(f'''q
1 1 1 rg 46.5 {top-46.7:.4f} 114 26.7 re f
0 0 0 rg
BT /AnnotationFont 12 Tf 25 {height-15:.4f} Td (Selected syn6 PF timing failures: right side) Tj ET
0.29 0.48 0.70 RG 1.7 w 38 {height-33:.4f} m 62 {height-33:.4f} l S
0 0 0 rg BT /AnnotationFont 10 Tf 68 {height-36:.4f} Td (Actual transmission) Tj ET
0.15 0.15 0.15 RG [5 3] 0 d 220 {height-33:.4f} m 244 {height-33:.4f} l S
0 0 0 rg BT /AnnotationFont 10 Tf 250 {height-36:.4f} Td (Intended primitive) Tj ET
BT /AnnotationFont 10 Tf 119 7 Td (Original calibration cycle \(percent\)) Tj ET
Q
'''.encode('ascii'))
tick_commands='\n'.join(f'BT /AnnotationFont 10 Tf {x*scale-(2.8 if value==0 else 5.6 if value<100 else 8.4):.4f} 23 Td ({value}) Tj ET' for value,x in zip((0,20,40,60,80,100),(140,295,450,605,760,915)))
stream.set_data(stream.get_data()+b'\nq 0 0 0 rg\n'+tick_commands.encode('ascii')+b'\nQ\n')
contents=page.get('/Contents')
if hasattr(contents,'get_object'):
    contents=contents.get_object()
from pypdf.generic import ArrayObject
if isinstance(contents,ArrayObject):
    contents.append(writer._add_object(stream))
else:
    page[NameObject('/Contents')]=ArrayObject([page['/Contents'],writer._add_object(stream)])
with target.open('wb') as f:
    writer.write(f)
target.with_suffix('.alt.txt').write_text('Representative unaltered right-side syn6 PF transmission panels S1, S2, and S5 from the full timing audit. Solid blue curves are actual transmission and dashed black curves are the intended primitive. S1/S2 timing is displaced and S5 remains active throughout the cycle.\n',encoding='utf-8')
print(target)
