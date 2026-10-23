from PIL import Image,ImageOps,ImageDraw
from pathlib import Path
root=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense')
out=root/'AE_Revision_20261023'/'build'
files=sorted((root/'deck_builder'/'deck_render').glob('slide_*.png'))
for start in range(0,len(files),12):
 canvas=Image.new('RGB',(1600,1020),'#dddddd'); d=ImageDraw.Draw(canvas)
 for j,p in enumerate(files[start:start+12]):
  im=Image.open(p); im.thumbnail((395,305)); x=(j%4)*400; y=(j//4)*340
  canvas.paste(im,(x,y+25)); d.text((x+5,y+5),p.stem,fill='black')
 canvas.save(out/f'original_montage_{start//12+1}.jpg')
