#!/usr/bin/env python3
"""Load the Simscape Multibody Link add-in into the running SW and enumerate
its COM interface, then attempt the CAD-XML export of the open assembly."""
import os
import sys

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402

DLL = r"C:\Program Files\MATLAB\R2025b\bin\win64\cl_sldwks2sm.dll"
REPAIR = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair"
EXPORT_DIR = os.path.join(REPAIR, "cad", "export")
os.makedirs(EXPORT_DIR, exist_ok=True)


def cop(obj, member):
    v = getattr(obj, member)
    return v() if callable(v) else v


def typeinfo_members(obj, label):
    print(f"--- {label} typeinfo ---")
    try:
        ole = obj._oleobj_
        tinfo = ole.GetTypeInfo(0)
        attr = tinfo.GetTypeAttr()
        print("interface:", tinfo.GetDocumentation(-1)[0], "funcs:", attr.cFuncs)
        for i in range(attr.cFuncs):
            fd = tinfo.GetFuncDesc(i)
            try:
                names = tinfo.GetNames(fd.memid)
                name = names[0] if isinstance(names, tuple) else names
            except Exception:
                name = f"memid{fd.memid}"
            nparams = getattr(getattr(fd, "descparams", None), "cArgs", None)
            print(f"   {name}  (memid={fd.memid}, invkind={fd.invkind}, nargs={nparams})")
    except Exception as e:
        print("typeinfo ERR", e)


sw = connect()
doc = sw.ActiveDoc
title = cop(doc, "GetTitle")
print("ActiveDoc:", title)
assert "10_AH_001" in str(title)

rc = sw.LoadAddIn(DLL)
print("LoadAddIn rc:", rc, "(0=OK, 1=error, 2=already loaded-ish)")

addin = None
for prog in ("sldwks2sm.Sldwks2smApp", "sldwks2sm.Sldwks2smApp.1",
             "SimscapeMultibodyLink.AddIn", "SimMechanicsLink.AddIn"):
    try:
        addin = sw.GetAddInObject(prog)
        print("GetAddInObject", prog, "->", addin)
        if addin is not None:
            break
    except Exception as e:
        print("GetAddInObject", prog, "ERR", e)

if addin is None:
    sys.exit("no add-in object handle")

typeinfo_members(addin, "MultibodyLink")

# dump one FUNCDESC's attributes to learn the arg-count field name
try:
    tinfo = addin._oleobj_.GetTypeInfo(0)
    attr = tinfo.GetTypeAttr()
    fd = tinfo.GetFuncDesc(7)   # ExportCADXML
    print("PyFUNCDESC attrs:", [a for a in dir(fd) if not a.startswith("_")])
    print("fd details:", {a: getattr(fd, a) for a in dir(fd) if not a.startswith("_") and not callable(getattr(fd, a))})
except Exception as e:
    print("fddesc ERR", e)

xml_path = os.path.join(EXPORT_DIR, "ah001_full.xml")
print("calling ExportCADXML() with no args ...")
try:
    res = addin.ExportCADXML()
    print("ExportCADXML() result:", res)
except Exception as e:
    print("ExportCADXML() ERR:", e)

import time
for _ in range(10):
    time.sleep(1)
print("export dir listing:", sorted(os.listdir(EXPORT_DIR))[:30])
# search typical default output locations
for base in (os.path.dirname(os.path.normpath(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\00_BipedalRobot_Redesign\10_AH_001.SLDASM")),):
    hits = [f for f in os.listdir(base) if f.lower().endswith(".xml") or "ah001" in f.lower()]
    print("asm dir xml-ish:", hits)
err_txt = os.path.join(EXPORT_DIR, "ah001_full_error.txt")
if os.path.exists(err_txt):
    with open(err_txt, "r", errors="replace") as f:
        print("--- error report ---")
        print(f.read()[:4000])
