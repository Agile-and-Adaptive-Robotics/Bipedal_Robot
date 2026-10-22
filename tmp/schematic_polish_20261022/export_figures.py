from pathlib import Path
import json, subprocess
import win32com.client
import fitz

HERE=Path(__file__).resolve().parent
STAGE=HERE/'staged'
GS=Path(r'D:\MiKTeX\miktex\bin\x64\mgs.exe')
app=win32com.client.DispatchEx('PowerPoint.Application')
try:
    for chapter,name in json.loads((STAGE/'jobs.json').read_text()):
        src=STAGE/chapter/'components'/f'{name}.pptx'
        pdf=STAGE/chapter/f'{name}.pdf'
        doc=app.Presentations.Open(str(src),ReadOnly=True,WithWindow=False)
        try:doc.SaveAs(str(pdf),32)
        finally:doc.Close()
        eps=src.with_suffix('.eps')
        subprocess.run([str(GS),'-q','-dNOPAUSE','-dBATCH','-sDEVICE=eps2write',f'-sOutputFile={eps}',str(pdf)],check=True)
        png=HERE/'renders'/f'{name}.png';png.parent.mkdir(exist_ok=True)
        subprocess.run([str(GS),'-q','-dNOPAUSE','-dBATCH','-sDEVICE=png16m','-r150','-dTextAlphaBits=4','-dGraphicsAlphaBits=4',f'-sOutputFile={png}',str(pdf)],check=True)
        page=fitz.open(pdf)[0]
        print(name,tuple(round(x,1) for x in page.rect),flush=True)
finally:app.Quit()
