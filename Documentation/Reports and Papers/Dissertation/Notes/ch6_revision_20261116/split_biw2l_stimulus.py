"""Compose unchanged recorded vector rows into two readable landscape strips."""
from pathlib import Path
from copy import deepcopy
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import NameObject,DictionaryObject,DecodedStreamObject,RectangleObject,ArrayObject

d=Path(__file__).resolve().parents[2]
source=d/'Defense/Session1_Handoff_20261104/experiments_20261104/BiW2LvarX_active_phase_stimulus.pdf'
original=PdfReader(source).pages[0]
width=597.6
scale=.4
sets=[('biw2l_stimulus_neural_20261116.pdf',[(81,299),(318,536),(555,773),(792,1010)],'Applied currents and RG-F, hip-flexor PF, and iliacus MN voltages from the recorded separate layer-current pulses and matched baseline.'),
      ('biw2l_stimulus_mechanics_20261116.pdf',[(1028,1247),(1265,1484),(1502,1722)],'Hip, knee, and ankle responses from the same recorded baseline and layer-current pulse trials, with identical time and angle coordinates.')]
for name,ranges,alt in sets:
    writer=PdfWriter()
    heights=[(b-a)*scale for a,b in ranges]
    height=sum(heights)+27
    page=writer.add_blank_page(width=width,height=height)
    y=height
    for (a,b),h in zip(ranges,heights):
        y-=h
        bottom=720-b*scale
        p=deepcopy(original)
        p.cropbox=RectangleObject([0,bottom,width,bottom+h])
        page.merge_transformed_page(p,Transformation().translate(0,y-bottom))
    font=DictionaryObject({NameObject('/Type'):NameObject('/Font'),NameObject('/Subtype'):NameObject('/Type1'),NameObject('/BaseFont'):NameObject('/Helvetica')})
    resources=page['/Resources']
    fonts=resources.get('/Font')
    if fonts is None:
        fonts=DictionaryObject();resources[NameObject('/Font')]=fonts
    fonts[NameObject('/AnnotationFont')]=writer._add_object(font)
    cmds=['q 0 0 0 rg']
    for value,x in zip(('7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0'),(284,454,622,790,958,1126,1294,1462)):
        cmds.append(f'BT /AnnotationFont 10 Tf {x*scale-7.5:.4f} 17 Td ({value}) Tj ET')
    cmds.append('BT /AnnotationFont 10 Tf 260 3 Td (Recorded time \(s\)) Tj ET Q')
    stream=DecodedStreamObject();stream.set_data('\n'.join(cmds).encode('ascii'))
    contents=page.get('/Contents')
    if hasattr(contents,'get_object'):contents=contents.get_object()
    if isinstance(contents,ArrayObject):contents.append(writer._add_object(stream))
    else:page[NameObject('/Contents')]=ArrayObject([page['/Contents'],writer._add_object(stream)])
    target=d/'Figures/95-AppendixD'/name
    with target.open('wb') as f:writer.write(f)
    target.with_suffix('.alt.txt').write_text(alt+'\n',encoding='utf-8')
    print(target)
