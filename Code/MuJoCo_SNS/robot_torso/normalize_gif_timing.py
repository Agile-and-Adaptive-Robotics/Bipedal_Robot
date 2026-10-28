"""Set GIF centisecond timing without altering any encoded image pixels."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
folder=ROOT/'Documentation/Reports and Papers/Dissertation/Defense/slideshow/AH001_20261028'
path=folder/'ah001_abductor_control.gif';raw=bytearray(path.read_bytes());before=hashlib.sha256(raw).hexdigest()
assert raw[:6] in (b'GIF87a',b'GIF89a')
p=13
if raw[10]&128:p+=3*(2**((raw[10]&7)+1))
positions=[]
def blocks(index):
    while raw[index]:index+=1+raw[index]
    return index+1
while raw[p]!=0x3b:
    if raw[p]==0x21:
        if raw[p+1]==0xf9:
            assert raw[p+2]==4;positions.append(p+4)
        p=blocks(p+2)
    elif raw[p]==0x2c:
        packed=raw[p+9];p+=10
        if packed&128:p+=3*(2**((packed&7)+1))
        p+=1;p=blocks(p)
    else:raise ValueError('Unexpected GIF block')
meta=json.loads((folder/'asset_provenance.json').read_text());fps=meta['gif_fps']
assert len(positions)==meta['gif_frames']
for i,position in enumerate(positions):
    delay=round((i+1)*100/fps)-round(i*100/fps)
    raw[position:position+2]=delay.to_bytes(2,'little')
path.write_bytes(raw)
meta['outputs'][path.name]=hashlib.sha256(raw).hexdigest()
meta['gif_duration_ms']=sum(int.from_bytes(raw[i:i+2],'little')*10 for i in positions)
meta['timing_normalization']=dict(encoded_image_pixels_unchanged=True,original_sha256=before,frames=len(positions))
(folder/'asset_provenance.json').write_text(json.dumps(meta,indent=2));print(meta['gif_duration_ms'],len(positions))
