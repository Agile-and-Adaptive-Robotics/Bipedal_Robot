
import os, sys
sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect
import win32com.client.dynamic
from win32com.client import VARIANT
import pythoncom
parts = sys.argv[1:]
sw = connect()
swd = win32com.client.dynamic.Dispatch("SldWorks.Application")
for p in parts:
    OUT = r'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\cad\stl'
    o = os.path.join(OUT, os.path.splitext(os.path.basename(p))[0] + ".STL")
    if os.path.exists(o) and os.path.getsize(o) > 0:
        continue
    try:
        errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        doc = swd.OpenDoc6(p, 1, 3, "", errs, warns)
        if doc is None:
            print("OPEN-FAIL", os.path.basename(p), flush=True); continue
        doc.SaveAs3(o, 0, 1)
        t = doc.GetTitle
        t = t() if callable(t) else t
        swd.QuitDoc(t)
        print("ok" if os.path.exists(o) and os.path.getsize(o) > 0 else "EMPTY",
              os.path.basename(p), flush=True)
    except Exception as e:
        print("FAIL", os.path.basename(p), str(e)[:60], flush=True)
try:
    sw.ExitApp()
except Exception:
    pass
