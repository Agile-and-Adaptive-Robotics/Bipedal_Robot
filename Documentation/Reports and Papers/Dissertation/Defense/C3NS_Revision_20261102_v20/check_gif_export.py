from pathlib import Path
import subprocess,json
import numpy as np
from PIL import Image
W=Path(__file__).resolve().parent
out=W/'build/gif_playback_frames';out.mkdir(exist_ok=True)
ffmpeg=Path(r'D:\Anaconda\envs\myo\Library\bin\ffmpeg.exe')
video=W/'build/route_gif_native_playback.mp4'
times=[1,3,5,7,9,11]
report=[]
for idx,slide in enumerate([46,47,49]):
    images=[]
    for t in times:
        path=out/f'slide_{slide:03d}_t{t:02d}.png'
        subprocess.run([str(ffmpeg),'-hide_banner','-loglevel','error','-y','-ss',str(idx*12+t),'-i',str(video),'-frames:v','1',str(path)],check=True)
        images.append(np.array(Image.open(path).convert('RGB'),dtype=np.int16))
    # Exclude captions and borders; test the model region of each native slide.
    roi=(850,130,1210,600) if slide==46 else (150,130,750,620) if slide==47 else (200,130,850,620)
    x0,y0,x1,y1=roi
    changes=[int((np.max(np.abs(a[y0:y1,x0:x1]-b[y0:y1,x0:x1]),axis=2)>40).sum()) for a,b in zip(images,images[1:])]
    report.append({'slide':slide,'sample_seconds':times,'model_roi':roi,'changed_pixels_between_samples':changes,'motion_verified':min(changes)>100,'scope':'Native PowerPoint video export of original slide objects; visible movement at six times across a full source GIF cycle.'})
assert all(x['motion_verified'] for x in report),report
(W/'build/native_gif_playback.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report))
