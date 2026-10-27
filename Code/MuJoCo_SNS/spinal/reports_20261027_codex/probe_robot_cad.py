"""Read-only assembly inspection; no mate or material changes are saved."""
import json
import sys
from pathlib import Path
import pythoncom
import win32com.client

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
sys.path.insert(0,str(Path.home()/'.agents/skills/solidworks/scripts'))
from sw_session import connect,constants

def value(obj,name,*args):
    member=getattr(obj,name)
    if hasattr(member,'_oleobj_'):
        return member
    return member(*args) if callable(member) else member

pythoncom.CoInitialize()
sw=connect()
sw=win32com.client.dynamic.DumbDispatch(sw._oleobj_)
print('SOLIDWORKS',value(sw,'RevisionNumber'),flush=True)
assembly=REPO/'Solid_Models/00_BipedalRobot_Redesign/10_AH_001.SLDASM'
errors=win32com.client.VARIANT(pythoncom.VT_BYREF|pythoncom.VT_I4,0)
warnings=win32com.client.VARIANT(pythoncom.VT_BYREF|pythoncom.VT_I4,0)
options=constants.swOpenDocOptions_Silent|constants.swOpenDocOptions_ReadOnly
doc=sw.OpenDoc6(str(assembly),constants.swDocASSEMBLY,options,'',errors,warnings)
if doc is None:
    raise RuntimeError(f'Assembly open failed: errors={errors.value}, warnings={warnings.value}')
out=dict(source=str(assembly),open_errors=errors.value,open_warnings=warnings.value,
         version=value(sw,'RevisionNumber'),components=[])
for component in value(doc,'GetComponents',False) or []:
    transform=value(component,'Transform2')
    out['components'].append(dict(name=value(component,'Name2'),
        path=value(component,'GetPathName'),fixed=bool(value(component,'IsFixed')),
        transform=list(value(transform,'ArrayData')) if transform else None))
(HERE/'robot_cad_inventory.json').write_text(json.dumps(out,indent=2))
print('Opened read-only;',len(out['components']),'components.',flush=True)
