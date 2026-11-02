#!/usr/bin/env python3
"""Try ISldwks2smApp::SaveAsSimMechXml(dispid 8) with the target XML path."""
import os
import sys
import time

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402
import pythoncom  # noqa: E402

EXPORT_DIR = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\cad\export"
XML = os.path.join(EXPORT_DIR, "ah001_full.xml")

sw = connect()
addin = sw.GetAddInObject("sldwks2sm.Sldwks2smApp")
disp = addin._oleobj_

for label, dispid, args in (("SaveAsSimMechXml(file)", 8, (XML,)),):
    print(label, "...")
    t0 = time.time()
    try:
        res = disp.Invoke(dispid, 0, pythoncom.DISPATCH_METHOD, True, *args)
        print(f"  returned {res!r} after {time.time()-t0:.1f}s")
    except Exception as e:
        print(f"  ERR after {time.time()-t0:.1f}s: {e}")

print("files:", sorted(os.listdir(EXPORT_DIR)) if os.path.isdir(EXPORT_DIR) else [])
for extra in ("ah001_full_error.txt",):
    p = os.path.join(EXPORT_DIR, extra)
    if os.path.exists(p):
        print("---", extra, "---")
        with open(p, "r", errors="replace") as f:
            print(f.read()[:6000])
