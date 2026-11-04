from pathlib import Path
import hashlib,json,subprocess
import numpy as np
from PIL import Image

W=Path(__file__).resolve().parent
meta=json.loads((W/'build/native_gif_export.json').read_text(encoding='utf-8-sig'))
assert meta['status']==3
video=Path(meta['output'])
ffmpeg=r'D:\Anaconda\envs\myo\Library\bin\ffmpeg.exe'
frames=W/'build/media_frames';frames.mkdir(exist_ok=True)
rows=[]
for slide,times in zip(meta['source_slides'],[(1,7),(15,21),(27,33),(37,40)]):
    images=[];paths=[]
    for t in times:
        p=frames/f'slide_{slide:03d}_t{t}.png'
        subprocess.run([ffmpeg,'-hide_banner','-loglevel','error','-y','-ss',str(t),'-i',str(video),'-frames:v','1','-update','1',str(p)],check=True)
        images.append(np.asarray(Image.open(p).convert('RGB'),dtype=float));paths.append(str(p))
    delta=np.abs(images[0]-images[1]);changed=float((delta.max(axis=2)>8).mean())
    rows.append({'slide':slide,'frame_seconds':times,'frames':paths,'mean_absolute_pixel_change':float(delta.mean()),'changed_pixel_fraction':changed,'motion_verified':changed>.001,'method':'Native PowerPoint video export; within-slide frame difference. Visual pose/crop inspection is recorded separately.'})
assert all(r['motion_verified'] for r in rows),rows
(W/'build/native_gif_playback.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
print(json.dumps([{'slide':r['slide'],'motion_verified':r['motion_verified'],'changed_fraction':round(r['changed_pixel_fraction'],4)} for r in rows]))
