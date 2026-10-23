from pathlib import Path
import win32com.client,json
from PIL import Image,ImageDraw
root=Path(__file__).resolve().parents[1]
ppt=root/'output/Bolen_Defense_Assertion_Evidence_20261023_v6.pptx'
try:
    app=win32com.client.GetActiveObject('PowerPoint.Application');spawned=False
except Exception:
    app=win32com.client.Dispatch('PowerPoint.Application');spawned=True
pres=app.Presentations.Open(str(ppt),ReadOnly=True,Untitled=False,WithWindow=False)
out=root/'build/render';out.mkdir(exist_ok=True)
audit=[]
for i in range(1,pres.Slides.Count+1):
    slide=pres.Slides(i)
    slide.Export(str(out/f'slide_{i:02}.png'),'PNG',1600,900)
    for shape in slide.Shapes:
        try:
            if shape.HasTextFrame and shape.TextFrame.HasText:
                tf=shape.TextFrame2
                audit.append({'slide':i,'text':shape.TextFrame.TextRange.Text,'height':shape.Height,'boundHeight':tf.TextRange.BoundHeight,'width':shape.Width,'boundWidth':tf.TextRange.BoundWidth,'left':shape.Left,'top':shape.Top})
        except Exception: pass
pres.SaveAs(str(root/'output/Bolen_Defense_Assertion_Evidence_20261023_v6.pdf'),32)
pres.Close()
if spawned:app.Quit()
(root/'build/powerpoint_geometry.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
files=sorted(out.glob('slide_*.png'))
for start in range(0,len(files),12):
    canvas=Image.new('RGB',(1600,780),'#ddd');d=ImageDraw.Draw(canvas)
    for j,p in enumerate(files[start:start+12]):
        im=Image.open(p);im.thumbnail((395,225));x=j%4*400;y=j//4*260
        canvas.paste(im,(x,y+25));d.text((x+5,y+5),p.stem,fill='black')
    canvas.save(root/f'build/final_montage_{start//12+1}.jpg')
print('PowerPoint exported',len(files),'slides and PDF')





