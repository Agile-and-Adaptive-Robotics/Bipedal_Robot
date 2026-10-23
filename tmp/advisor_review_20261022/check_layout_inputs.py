from pathlib import Path
import zipfile,difflib
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';z=zipfile.ZipFile(D/'Bolen_Dissertation.zip')
for name in ['main.tex','psutex.cls']:
    old=z.read(name).decode('utf-8');new=(D/'ProofFinal'/name).read_text(encoding='utf-8')
    print(name,'same',old.replace('\r\n','\n')==new.replace('\r\n','\n'))
    if old!=new:print(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True)))[:3500])
