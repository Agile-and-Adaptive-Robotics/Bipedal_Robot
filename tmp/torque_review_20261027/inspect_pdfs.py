from pathlib import Path
import json, subprocess, logging
logging.getLogger('pypdf').setLevel(logging.ERROR)
from pypdf import PdfReader
ROOT=Path(r'D:\GitHub\Bipedal_Robot');OUT=ROOT/'tmp/torque_review_20261027'
files={
'paper':ROOT/'Documentation/Reports and Papers/Knee_Torque_Test/build_review_20261027/Main.pdf',
'supplement':ROOT/'Documentation/Reports and Papers/Knee_Torque_Test/build_review_20261027/frontiers_SupplementaryMaterial.pdf',
'dissertation':ROOT/'Documentation/Reports and Papers/Dissertation/ProofFinal/build_r12_torque_review_20261027/main.pdf'}
terms=['Local bracket displacement','BPA crossing a flexed','Identification of the Torque','Identification Procedure','Actuator-Count and Transmission','Whole-Leg Transmission','Capability across 27','Historical single-bracket','Bracket Reference Frames','Penalty terms of the mesh','Bracket stiffness arrays','Optimizer parameterization','Degrees of freedom and nominal','Afferent classes represented','Terms of the compliance-unaware']
mapping={}
for name,path in files.items():
    doc=PdfReader(path).pages; mapping[name]={'pages':len(doc),'targets':{}}
    texts=[];pages=set()
    for i,page in enumerate(doc):
        text=page.extract_text();texts.append(f'\n\nPDF PAGE {i+1}\n{text}')
        for term in terms:
            if term in text:
                mapping[name]['targets'].setdefault(term,[]).append(i+1)
                if i>15 or name!='dissertation':pages.add(i)
    if name=='paper':pages.update(range(len(doc)))
    if name=='supplement':pages.update(range(len(doc)))
    (OUT/f'{name}_extracted.txt').write_text(''.join(texts),encoding='utf-8')
    for i in sorted(pages):
        poppler=Path(r'C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
        subprocess.run([str(poppler),'-f',str(i+1),'-l',str(i+1),'-scale-to','1250','-singlefile','-png',str(path),str(OUT/f'{name}_page_{i+1:03}')],check=True,capture_output=True)
    if name in ['paper','supplement']:
        from PIL import Image,ImageOps,ImageDraw
        nums=sorted(pages)
        for batch in range(0,len(nums),12):
            group=nums[batch:batch+12]
            sheet=Image.new('RGB',(1000,370*((len(group)+3)//4)),'#d9d9d9')
            draw=ImageDraw.Draw(sheet)
            for j,n in enumerate(group):
                im=Image.open(OUT/f'{name}_page_{n+1:03}.png');im.thumbnail((245,345))
                x=(j%4)*250;y=(j//4)*370
                sheet.paste(im,(x,y+20));draw.text((x+8,y+3),f'{name} page {n+1}',fill='black')
            sheet.save(OUT/f'{name}_contact_{batch//12+1}.png')
    print(name,len(doc),'pages',mapping[name]['targets'])
(OUT/'page_map.json').write_text(json.dumps(mapping,indent=2),encoding='utf-8')
