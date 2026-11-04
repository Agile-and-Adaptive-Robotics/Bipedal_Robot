from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
here=Path(__file__).resolve().parents[1]/'packaging'
frames=[]
for k,p in enumerate(sorted((here/'build/frames').glob('knee_*.png'))):
 im=Image.open(p).convert('RGB');d=ImageDraw.Draw(im)
 try:font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',22)
 except OSError:font=ImageFont.load_default()
 d.text((24,22),'AH001 | Vastus intermedius packaging preview',fill='#174d24',font=font)
 d.text((24,52),f'CAD shank rotation: {-5*k} degrees | diameter follows saved knee contraction',fill='#273c43',font=font)
 d.text((24,920),'Preview registration and sleeve placement; physical tackle hardware is not specified.',fill='#46565c',font=font)
 frames.append(im.quantize(colors=128))
assert len(frames)==19
seq=frames+frames[-2:0:-1]
seq[0].save(here/'output/Defense_Knee_Sweep.gif',save_all=True,append_images=seq[1:],duration=140,loop=0,optimize=False)
print('GIF:',here/'output/Defense_Knee_Sweep.gif')
