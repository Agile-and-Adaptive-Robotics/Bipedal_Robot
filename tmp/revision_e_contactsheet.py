from pathlib import Path
from PIL import Image,ImageDraw
E=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\Revision_E')
files=sorted((E/'review').glob('page-*.png'))
canvas=Image.new('RGB',(1600,((len(files)+3)//4)*550),'#cccccc')
for i,p in enumerate(files):
 im=Image.open(p);im.thumbnail((390,520))
 x=(i%4)*400;y=(i//4)*550
 canvas.paste(im,(x,y+20));ImageDraw.Draw(canvas).text((x+8,y+2),p.stem,fill='black')
canvas.save(E/'review/contactsheet.png')
