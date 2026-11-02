from PIL import Image,ImageSequence
from pathlib import Path
import json
im=Image.open(r'D:\GitHub\Bipedal_Robot\Knee_Flexor_20mm.gif')
counts=[]
for f in ImageSequence.Iterator(im):
 pix=f.convert('RGB').crop((464,124,512,160))
 counts.append(sum(max(c)-min(c)>25 for c in pix.getdata()))
result={'frames':len(counts),'mask_source_bounds':[464,124,512,160],'colored_geometry_pixels_max':max(counts),'scope':'Top-left clipped coordinate-axis fragment; femur is blue and outside this rectangle in all frames.'}
assert max(counts)==0,result
(Path(__file__).parent/'build/gif_mask_check.json').write_text(json.dumps(result,indent=2))
print(result)
