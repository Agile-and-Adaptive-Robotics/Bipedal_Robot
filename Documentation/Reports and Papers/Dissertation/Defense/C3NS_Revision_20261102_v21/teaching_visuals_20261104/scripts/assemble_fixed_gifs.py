"""Encode verified repaired renders as looping, aspect-preserving GIFs."""
from pathlib import Path
from PIL import Image
import json,hashlib
root=Path(__file__).resolve().parents[1]/'build_fixed_20261105'
jobs=[('visual1','cam_extensor_loop.gif',120),('visual2','reverse_pulley_loop.gif',106),('visual3','knee_frames_loop.gif',136)]
manifest=[]
for folder,name,count in jobs:
 frames=sorted((root/folder/'frames').glob('f*.png'))
 assert len(frames)==count,(folder,len(frames),count)
 # A shared palette sampled throughout the animation retains moving colors.
 samples=Image.new('RGB',(100*12,100),'white')
 for j in range(12):
  with Image.open(frames[round(j*(count-1)/11)]) as im:
   im=im.convert('RGB');im.thumbnail((100,100));samples.paste(im,(100*j,0))
 palette=samples.quantize(colors=256,method=Image.MEDIANCUT)
 ims=[]
 for path in frames:
  with Image.open(path) as im:
   im=im.convert('RGB');im.thumbnail((900,900),Image.LANCZOS)
   ims.append(im.quantize(palette=palette,dither=Image.FLOYDSTEINBERG))
 durations=[(60,70,70)[j%3] for j in range(count)]
 out=root/folder/name
 ims[0].save(out,save_all=True,append_images=ims[1:],duration=durations,loop=0,disposal=2,optimize=False)
 with Image.open(out) as gif:
  assert gif.n_frames>=50,(folder,gif.n_frames)
  assert gif.info.get('loop')==0
  manifest.append({'visual':folder,'path':str(out),'size':gif.size,'rendered_frames':count,'frames':gif.n_frames,'duration_ms':sum(durations),'sha256':hashlib.sha256(out.read_bytes()).hexdigest()})
 print(name,'encoded',flush=True)
 del ims
(root/'GIF_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('ALL_THREE_GIFS_VERIFIED')
