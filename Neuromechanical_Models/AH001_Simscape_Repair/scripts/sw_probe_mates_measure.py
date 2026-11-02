#!/usr/bin/env python3
"""sw_probe_mates_measure.py -- one SW pass, read-only:
1. enumerate the FOOT assembly's mates (name/type/components) via any working
   API route (FirstFeature variants, GetMates, Extension) for the (A) cleanup
2. measure the LEG assembly: tibia component boxes, knob component boxes
3. measure the FOOT assembly: component boxes (esp. sole-defining parts) and
   the Ankle Exterior (housing) placement
Everything logged to JSON. No files modified.
"""
import json
import os
import sys

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402
import win32com.client.dynamic  # noqa: E402
from win32com.client import VARIANT  # noqa: E402
import pythoncom  # noqa: E402

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
LEG = os.path.join(REPO, r"Solid_Models\00_BipedalRobot_Redesign\10_AH_001.SLDASM")
FOOT = os.path.join(REPO, r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Foot Ankle.SLDASM")
OUT = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair\logs\sw_probe_mates_20261102.json")


def cop(obj, member):
    v = getattr(obj, member)
    return v() if callable(v) else v


def open_asm(swd, path):
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = swd.OpenDoc6(path, 2, 3, "", errs, warns)  # silent + readonly
    return doc, errs.value, warns.value


def try_mate_routes(doc, log):
    """Try every plausible mate-enumeration API; log what works."""
    routes = {}
    # route 1: AssemblyDoc.GetMates
    try:
        mates = doc.GetMates()
        routes["GetMates"] = str(type(mates))
        if mates:
            routes["GetMates_count"] = len(mates)
            log["mates_GetMates"] = []
            for mte in mates:
                try:
                    e = {"type": str(mte.Type)}
                    try:
                        e["name"] = str(cop(mte, "Name"))
                    except Exception:
                        pass
                    log["mates_GetMates"].append(e)
                except Exception as ex:
                    log["mates_GetMates"].append({"err": str(ex)})
    except Exception as ex:
        routes["GetMates"] = f"ERR {ex}"
    # route 2: feature walk via name variants
    feat = None
    for cand in ("FirstFeature", "IGetFirstFeature", "GetFirstFeature"):
        try:
            f = getattr(doc, cand)
            feat = f() if callable(f) else f
            routes[cand] = "OK" if feat is not None else "None"
            if feat is not None:
                routes["feat_route"] = cand
                break
        except Exception as ex:
            routes[cand] = f"ERR {ex}"
    if feat is not None:
        items = []
        guard = 0
        while feat is not None and guard < 4000:
            guard += 1
            try:
                tname = str(cop(feat, "GetTypeName2"))
            except Exception:
                tname = ""
            if "mate" in tname.lower():
                ent = {"ftype": tname}
                for nm_member in ("Name",):
                    try:
                        ent["name"] = str(cop(feat, nm_member))
                    except Exception:
                        pass
                try:
                    fd = cop(feat, "GetDefinition")
                    if fd is not None:
                        try:
                            ent["mtype"] = str(fd.Type)
                        except Exception:
                            pass
                except Exception:
                    pass
                items.append(ent)
            try:
                nxt = getattr(feat, "GetNextFeature")
                feat = nxt() if callable(nxt) else nxt
            except Exception:
                break
        routes["mate_features"] = len(items)
        log["mate_features"] = items
    # route 3: Extension.GetMates? / doc.Mates?
    for nm in ("Mates", "GetMates2", "MateGroups"):
        try:
            v = getattr(doc, nm)
            v2 = v() if callable(v) else v
            routes[nm] = "OK" if v2 is not None else "None"
        except Exception as ex:
            routes[nm] = f"ERR {ex}"
    return routes


def comp_boxes(doc, log, key):
    boxes = {}
    try:
        comps = doc.GetComponents(False)
        for c in comps:
            try:
                name = str(cop(c, "Name2"))
            except Exception:
                continue
            box = None
            for meth, args in (("GetBox", (1,)), ("GetBox", ()), ("GetBodyBox", ())):
                try:
                    b = getattr(c, meth)(*args)
                    if b:
                        box = list(b)
                        break
                except Exception:
                    continue
            boxes[name] = box
    except Exception as ex:
        boxes["_error"] = str(ex)
    log[key] = boxes


sw = connect()
swd = win32com.client.dynamic.Dispatch("SldWorks.Application")
log = {"sw": str(cop(swd, "RevisionNumber"))}

# ---- foot: mates + boxes
doc, errs, warns = open_asm(swd, FOOT)
log["foot_open"] = {"errors": errs, "warnings": warns}
doc = swd.ActiveDoc
log["foot_mate_routes"] = try_mate_routes(doc, log)
comp_boxes(doc, log, "foot_component_boxes")
try:
    log["foot_assembly_box"] = list(doc.GetBox(1))
except Exception as ex:
    log["foot_assembly_box"] = f"ERR {ex}"
swd.QuitDoc(str(cop(doc, "GetTitle")))

# ---- leg: boxes only
doc, errs, warns = open_asm(swd, LEG)
log["leg_open"] = {"errors": errs, "warnings": warns}
doc = swd.ActiveDoc
comp_boxes(doc, log, "leg_component_boxes")
try:
    log["leg_assembly_box"] = list(doc.GetBox(1))
except Exception as ex:
    log["leg_assembly_box"] = f"ERR {ex}"
swd.QuitDoc(str(cop(doc, "GetTitle")))

with open(OUT, "w") as f:
    json.dump(log, f, indent=1)
print("routes:", json.dumps(log.get("foot_mate_routes"), indent=1))
print("mate feature count:", len(log.get("mate_features", [])))
print("sample:", log.get("mate_features", [])[:6])
print("log:", OUT)
