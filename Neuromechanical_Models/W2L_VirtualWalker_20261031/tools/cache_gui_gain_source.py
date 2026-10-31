from pathlib import Path
import urllib.request,json,hashlib
W=Path(__file__).resolve().parents[1]/'build/native_contact_unit_probe_20261031';D=W/'gui_source';D.mkdir(exist_ok=True)
report=[]
for name,path in [('Polynomial','DataObjects/Gains/Polynomial.vb'),('Sigmoid','DataObjects/Gains/Sigmoid.vb'),('ScaledNumber','Framework/ScaledNumber.vb')]:
 url='https://raw.githubusercontent.com/NeuroRoboticTech/AnimatLabPublicSource/master/Libraries/AnimatGUI/'+path
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'local-unit-audit'});data=urllib.request.urlopen(req,timeout=12).read();(D/(name+'.vb')).write_bytes(data);report.append({'url':url,'file':str(D/(name+'.vb')),'sha256':hashlib.sha256(data).hexdigest()});print(name,len(data))
 except Exception as e:report.append({'url':url,'error':str(e)});print(name,str(e))
(D/'SOURCE_MANIFEST.json').write_text(json.dumps(report,indent=2))
