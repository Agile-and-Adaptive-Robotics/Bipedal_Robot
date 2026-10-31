from pathlib import Path
import urllib.request
W=Path(__file__).resolve().parents[1]/'build/guided_source_audit'
for folder,name in [('BulletAnimatSim','BlRigidBody'),('AnimatGUI/DataObjects/Physical','RigidBody')]:
 ext='.vb' if folder.startswith('AnimatGUI') else '.cpp';url='https://raw.githubusercontent.com/NeuroRoboticTech/AnimatLabPublicSource/master/Libraries/'+folder+'/'+name+ext
 try:
  b=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'native-telemetry-audit'}),timeout=12).read();(W/(name+ext)).write_bytes(b);print(name,len(b))
 except Exception as e:print(name,str(e))
