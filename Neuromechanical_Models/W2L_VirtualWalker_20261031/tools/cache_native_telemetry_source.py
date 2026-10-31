from pathlib import Path
import urllib.request,json
W=Path(__file__).resolve().parents[1]/'build/guided_source_audit';W.mkdir(parents=True,exist_ok=True)
for name in ['MovableItem','BodyPart']:
 url='https://raw.githubusercontent.com/NeuroRoboticTech/AnimatLabPublicSource/master/Libraries/AnimatSim/'+name+'.cpp'
 try:
  b=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'local-model-audit'}),timeout=12).read();(W/(name+'.cpp')).write_bytes(b);print(name,len(b))
 except Exception as e:print(name,str(e))
