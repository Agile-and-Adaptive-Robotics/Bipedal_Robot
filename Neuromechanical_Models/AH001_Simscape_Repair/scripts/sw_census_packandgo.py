#!/usr/bin/env python3
"""AH001 repair - SW census + Pack and Go (laptop, 2026-11-01).

Runs against the ALREADY-OPEN original assembly (kept open by prior probe):
  1. component census (name, path, suppression) -> JSON
  2. feature walk: mate features (names) -> identify four-bar knee mates
  3. Pack and Go the assembly into AH001_Simscape_Repair/cad/ (self-contained copy)
  4. registry probe for the Simscape Multibody Link SW add-in + find its DLL
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
REPAIR = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
CAD_DIR = os.path.join(REPAIR, "cad")
LOG = os.path.join(REPAIR, "logs", "sw_census_20261101.json")
os.makedirs(CAD_DIR, exist_ok=True)


def cop(obj, member):
    v = getattr(obj, member)
    return v() if callable(v) else v


sw = connect()
doc = sw.ActiveDoc
title = cop(doc, "GetTitle")
print("ActiveDoc:", title)
assert "10_AH_001" in str(title), "expected the original assembly to be open"

out = {"title": str(title), "sw": cop(sw, "RevisionNumber")}

# ---------- 1. component census
comps = doc.GetComponents(False)
census = []
suppressed = 0
for c in comps:
    try:
        name = str(cop(c, "Name2"))
    except Exception:
        name = "?"
    try:
        path = str(cop(c, "GetPathName"))
    except Exception:
        path = ""
    sup = False
    try:
        ss = cop(c, "SuppressionState")
        sup = int(ss) == 1  # swComponentSuppressed
    except Exception:
        try:
            sv = cop(c, "Suppress")
            sup = int(sv) != 0
        except Exception:
            pass
    if sup:
        suppressed += 1
    census.append({"name": name, "path": path, "suppressed": sup})
out["n_components"] = len(census)
out["n_suppressed"] = suppressed
out["n_missing_paths"] = sum(1 for c in census if not c["path"])
out["components"] = census
print(f"components={len(census)} suppressed={suppressed} missing_paths={out['n_missing_paths']}")

# ---------- 2. mate feature census
mates = []
feat = None
fm = None
try:
    fm = cop(doc, "FeatureManager")
except Exception as e:
    print("FeatureManager ERR", e)
for getter in (lambda: doc.GetMates(),
               lambda: cop(cop(fm, "GetFirstFeature"), "GetNextFeature"),
               lambda: cop(doc, "FirstFeature")):
    try:
        x = getter()
        print("mate-census route ->", type(x).__name__)
        break
    except Exception as e:
        print("mate-census route fail:", e)

guard = 0
try:
    feat = cop(fm, "GetFirstFeature")
except Exception as e:
    print("fm.GetFirstFeature ERR", e)
while feat is not None and guard < 8000:
    guard += 1
    try:
        tname = str(cop(feat, "GetTypeName2"))
    except Exception:
        tname = ""
    if "mate" in tname.lower():
        try:
            mname = str(cop(feat, "Name"))
        except Exception:
            mname = "?"
        mtype = ""
        try:
            fdef = cop(feat, "GetDefinition")
            mtype = str(cop(fdef, "Type"))
        except Exception:
            pass
        mates.append({"name": mname, "ftype": tname, "mtype": mtype})
    try:
        feat = cop(feat, "GetNextFeature")
    except Exception:
        break
out["n_mate_features"] = len(mates)
out["mate_features"] = mates
print(f"mate features={len(mates)}")
kneed = [m for m in mates if any(k in m["name"].upper() for k in ("KT", "KB", "FL_", "BL_", "KNEE", "HINGE", "PIN"))]
out["knee_related_mates"] = kneed
print("knee-related mate names (first 40):")
for m in kneed[:40]:
    print("   ", m)

# ---------- 3. Pack and Go
try:
    ext = doc.Extension
    pg = cop(ext, "PackAndGo")
    # target: single flat folder with suffix _r  (keeps names unique)
    names = pg.GetDocumentNames()
    to_names = []
    for n in names:
        to_names.append(os.path.join(CAD_DIR, os.path.basename(n)))
    ok_name = pg.SetSaveToName(False, os.path.join(CAD_DIR, "10_AH_001_packed.SLDASM"))
    ok_flatten = pg.FlattenToSingleFolder(True, CAD_DIR)
    # after SetSaveToName/flatten, set each doc's target explicitly (belt+braces)
    st = pg.SetDocumentNames(to_names)
    from win32com.client import VARIANT
    import pythoncom
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    res = ext.SaveAs3(os.path.join(CAD_DIR, "10_AH_001_packed.SLDASM"), 0, 1)  # current version, silent
    out["packandgo"] = {"GetDocumentNames": len(names), "SetSaveToName": bool(ok_name),
                        "FlattenToSingleFolder": bool(ok_flatten), "SetDocumentNames": bool(st),
                        "SaveAs3_result": res}
    print("PackAndGo:", out["packandgo"], "files now in cad/:", len(os.listdir(CAD_DIR)))
except Exception as e:
    out["packandgo_error"] = str(e)
    print("PackAndGo ERR:", e)

# ---------- 4. Multibody Link add-in discovery
addins = []
for hive in ("HKLM", "HKCU"):
    try:
        r = subprocess.run(["reg", "query", hive + r"\SOFTWARE\SolidWorks\AddIns", "/s"],
                           capture_output=True, text=True, timeout=30)
        txt = r.stdout
    except Exception:
        txt = ""
    for line in txt.splitlines():
        if any(k in line.lower() for k in ("simscape", "simmechanics", "smlink", "smls")):
            addins.append(f"{hive}: {line.strip()}")
out["registry_link_mentions"] = addins[:40]
print("registry link mentions:", len(addins))

# MATLAB link folder discovery
matlab_link = []
for base in (r"C:\Program Files\MATLAB",):
    if os.path.isdir(base):
        for rel, dirs, files in os.walk(os.path.join(base, "R2025b", "toolbox", "physmod", "sm", "link")):
            for f in files:
                if f.lower().endswith((".dll", ".addin")) and "sw" in rel.lower():
                    matlab_link.append(os.path.join(rel, f))
            if len(matlab_link) > 30:
                break
out["matlab_link_files"] = matlab_link[:30]
print("matlab link dlls:", len(matlab_link))

with open(LOG, "w") as f:
    json.dump(out, f, indent=1)
print("log:", LOG)
