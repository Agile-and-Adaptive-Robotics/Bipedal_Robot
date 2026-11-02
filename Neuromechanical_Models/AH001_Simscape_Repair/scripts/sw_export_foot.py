#!/usr/bin/env python3
"""Export Foot Ankle.SLDASM via Multibody Link + read material assignments.

1. Connect (launch SW).
2. Open the LEG assembly read-only first: read each unique part's material
   name + density (verify Ben's Onyx/steel ruling).
3. Close it, open Foot Ankle.SLDASM read-only, export CAD XML to
   cad/export_foot/ via SaveAsSimMechXml, dismiss the completion dialog by
   posting WM_COMMAND/IDOK.
4. Close docs; leave SW running for now (more work may follow).
"""
import json
import os
import sys
import time
import threading

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402
import win32com.client.dynamic  # noqa: E402
from win32com.client import VARIANT  # noqa: E402
import pythoncom  # noqa: E402
import win32con  # noqa: E402
import win32gui  # noqa: E402

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
LEG_ASM = os.path.join(REPO, r"Solid_Models\00_BipedalRobot_Redesign\10_AH_001.SLDASM")
FOOT_ASM = os.path.join(REPO, r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Foot Ankle.SLDASM")
REPAIR = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
EXPORT_DIR = os.path.join(REPAIR, "cad", "export_foot")
os.makedirs(EXPORT_DIR, exist_ok=True)
LOG = os.path.join(REPAIR, "logs", "sw_foot_export_20261101.json")

DLL = r"C:\Program Files\MATLAB\R2025b\bin\win64\cl_sldwks2sm.dll"


def cop(obj, member):
    v = getattr(obj, member)
    return v() if callable(v) else v


def open_asm(path, readonly=True):
    swd = win32com.client.dynamic.Dispatch("SldWorks.Application")
    opts = 1 | (2 if readonly else 0)
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = swd.OpenDoc6(path, 2, opts, "", errs, warns)
    return doc, errs.value, warns.value


def close_doc(sw, title):
    try:
        sw.QuitDoc(title)
    except Exception:
        pass


def dismiss_export_dialog(timeout=240):
    """Wait for the #32770 'Export Complete' dialog and post IDOK."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            hwnd = win32gui.FindWindow("#32770", "Simscape Multibody Link Export Complete")
            if hwnd:
                win32gui.PostMessage(hwnd, win32con.WM_COMMAND, 1, 0)
                time.sleep(1.0)
                if not win32gui.IsWindow(hwnd):
                    return True
                win32gui.PostMessage(hwnd, win32con.WM_CLOSE, 0, 0)
                time.sleep(1.0)
                return not win32gui.IsWindow(hwnd)
        except Exception:
            pass
        time.sleep(1.0)
    return False


report = {}

sw = connect()
report["sw"] = cop(sw, "RevisionNumber")

# ---- 1. materials: open each unique PART directly and read its material name
swd = win32com.client.dynamic.Dispatch("SldWorks.Application")

def part_material(part_path):
    try:
        errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        pdoc = swd.OpenDoc6(part_path, 1, 3, "", errs, warns)  # part, silent+readonly
        if pdoc is None:
            return "open-fail"
        m = pdoc.GetMaterialPropertyName2("", "")
        m = m() if callable(m) else m
        swd.QuitDoc(cop(pdoc, "GetTitle"))
        return str(m)
    except Exception as e:
        return f"err {e}"

# leg assembly component paths were captured earlier in sw_census_20261101.json
census = json.load(open(os.path.join(REPAIR, "logs", "sw_census_20261101.json")))
mats = {}
seen = set()
for c in census["components"]:
    base = c["name"].split("-")[0] if "-" in c["name"] else c["name"]
    if base in seen or not c["path"]:
        continue
    seen.add(base)
    mats[base] = part_material(c["path"])
report["leg_materials"] = mats
print(json.dumps(mats, indent=1)[:2500])

# ---- 2. foot export
rc = swd.LoadAddIn(DLL)
report["LoadAddIn"] = rc
addin = swd.GetAddInObject("sldwks2sm.Sldwks2smApp")
report["addin"] = addin is not None
if addin is None:
    with open(LOG, "w") as f:
        json.dump(report, f, indent=1)
    sys.exit("no addin")

doc, errs, warns = open_asm(FOOT_ASM, readonly=True)
report["foot_open"] = {"errors": errs, "warnings": warns}
doc = swd.ActiveDoc
title = cop(doc, "GetTitle")
report["foot_title"] = str(title)

# foot materials too
mats = {}
seen = set()
try:
    comps = doc.GetComponents(False)
    for c in comps:
        name = str(cop(c, "Name2"))
        base = name.split("-")[0] if "-" in name else name
        if base in seen:
            continue
        seen.add(base)
        try:
            pdoc = c.GetModelDoc2()
            if pdoc is None:
                continue
            mats[base] = str(pdoc.GetMaterialPropertyName2("", ""))
        except Exception:
            continue
except Exception as e:
    mats["_error"] = str(e)
report["foot_materials"] = mats

xml_path = os.path.join(EXPORT_DIR, "foot_ankle.xml")
print("exporting ->", xml_path)
disp = addin._oleobj_
result_box = {}

def runner():
    try:
        pythoncom.CoInitialize()
        result_box["res"] = disp.Invoke(8, 0, pythoncom.DISPATCH_METHOD, True, xml_path)
    except Exception as e:
        result_box["err"] = str(e)

# marshal problem: run in main thread, dismiss dialog from a helper thread
dismiss_thread = threading.Thread(target=lambda: dismiss_export_dialog(300), daemon=True)
dismiss_thread.start()
try:
    res = disp.Invoke(8, 0, pythoncom.DISPATCH_METHOD, True, xml_path)
    report["export_result"] = repr(res)
except Exception as e:
    report["export_error"] = str(e)
dismiss_thread.join(timeout=10)

time.sleep(2)
report["export_dir"] = sorted(os.listdir(EXPORT_DIR))
close_doc(swd, str(title))

with open(LOG, "w") as f:
    json.dump(report, f, indent=1)
print(json.dumps({k: v for k, v in report.items() if k != "export_dir"}, indent=1)[:3000])
print("files:", report["export_dir"][:25])
