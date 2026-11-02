#!/usr/bin/env python3
"""Raw IDispatch::Invoke attempts on ISldwks2smApp::ExportCADXML (dispid 7)."""
import os
import sys
import threading
import time

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402

import pythoncom  # noqa: E402

EXPORT_DIR = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\cad\export"
XML = os.path.join(EXPORT_DIR, "ah001_full.xml")

sw = connect()
addin = sw.GetAddInObject("sldwks2sm.Sldwks2smApp")
print("addin:", addin)
disp = addin._oleobj_

results = {}


def call_raw(name, dispid, *args, timeout=90):
    box = {}

    def runner():
        try:
            pythoncom.CoInitialize()
            box["res"] = disp.Invoke(dispid, 0, pythoncom.DISPATCH_METHOD, True, *args)
        except Exception as e:
            box["err"] = str(e)

    th = threading.Thread(target=runner, daemon=True)
    th.start()
    th.join(timeout)
    if th.is_alive():
        return f"TIMEOUT after {timeout}s (still running)"
    if "err" in box:
        return "ERR " + box["err"]
    return "OK " + repr(box.get("res"))


print("ExportCADXML(filename) main-thread invoke ...")
t0 = time.time()
try:
    res = disp.Invoke(7, 0, pythoncom.DISPATCH_METHOD, True, XML)
    print("returned", repr(res), f"after {time.time()-t0:.1f}s")
except Exception as e:
    print(f"ERR after {time.time()-t0:.1f}s:", e)
print("files:", os.listdir(EXPORT_DIR) if os.path.isdir(EXPORT_DIR) else None)
err_txt = os.path.join(EXPORT_DIR, "ah001_full_error.txt")
if os.path.exists(err_txt):
    print("--- error report ---")
    with open(err_txt, "r", errors="replace") as f:
        print(f.read()[:5000])
