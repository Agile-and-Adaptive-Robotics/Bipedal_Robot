#!/usr/bin/env python3
"""sw_batch_stl.py -- export every leg-assembly part + foot part to coarse STL
for the route-render visualization. Reuses the early-bound COM route."""
import json
import os
import sys

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402
import win32com.client.dynamic  # noqa: E402
import win32com.client.gencache as gc  # noqa: E402
import win32com.client.makepy as mk  # noqa: E402
from win32com.client import VARIANT  # noqa: E402
import pythoncom  # noqa: E402

mk.GenerateFromTypeLibSpec("SldWorks.Application", 0)
M = gc.GetModuleForProgID("SldWorks.Application")

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
OUT = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair\cad\stl")
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair\logs\sw_batch_stl_20261102.json")

report = {"done": [], "fail": []}
census = json.load(open(os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair\logs\sw_census_20261101.json")))
parts = []
seen = set()
for c in census["components"]:
    p = c["path"]
    if not p or not p.lower().endswith(".sldprt"):
        continue
    base = os.path.basename(p)
    if base in seen:
        continue
    seen.add(base)
    parts.append(p)
# foot parts
FOOT_DIR = os.path.join(REPO, r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)")
for f in os.listdir(FOOT_DIR):
    if f.lower().endswith(".sldprt") and f not in seen:
        parts.append(os.path.join(FOOT_DIR, f))
        seen.add(f)
print("parts to export:", len(parts))

sw = connect()
swd = win32com.client.dynamic.Dispatch("SldWorks.Application")

for p in parts:
    out = os.path.join(OUT, os.path.splitext(os.path.basename(p))[0] + ".STL")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        report["done"].append({"part": os.path.basename(p), "stl": out, "cached": True})
        continue
    try:
        errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        doc = swd.OpenDoc6(p, 1, 3, "", errs, warns)  # part, silent+readonly
        if doc is None:
            report["fail"].append({"part": os.path.basename(p), "err": "open-fail"})
            continue
        ok = doc.SaveAs3(out, 0, 1)  # current version, silent
        title = doc.GetTitle
        title = title() if callable(title) else title
        swd.QuitDoc(title)
        sz = os.path.getsize(out) if os.path.exists(out) else 0
        (report["done"] if sz > 0 else report["fail"]).append(
            {"part": os.path.basename(p), "stl": out, "bytes": sz, "saveas3": str(ok)})
    except Exception as e:
        report["fail"].append({"part": os.path.basename(p), "err": str(e)})

with open(LOG, "w") as f:
    json.dump(report, f, indent=1)
print("done:", len(report["done"]), "fail:", len(report["fail"]))
for fl in report["fail"][:8]:
    print("  FAIL", fl)
