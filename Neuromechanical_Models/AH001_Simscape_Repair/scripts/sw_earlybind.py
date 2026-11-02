#!/usr/bin/env python3
"""sw_probe2.py -- early-bound mate enumeration + component boxes (read-only).

Route that works (2026-11-02): makepy by ProgID IN-PROCESS, then
    M = gencache.GetModuleForProgID('SldWorks.Application')
    doc = M.ModelDoc2(swd.ActiveDoc._oleobj_)
gives IModelDoc2 with FirstFeature/GetNextFeature (property GetTypeName2).
"""
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

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
LEG = os.path.join(REPO, r"Solid_Models\00_BipedalRobot_Redesign\10_AH_001.SLDASM")
FOOT = os.path.join(REPO, r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Foot Ankle.SLDASM")
OUT = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair\logs\sw_probe_mates_20261102.json")

mk.GenerateFromTypeLibSpec("SldWorks.Application", 0)
M = gc.GetModuleForProgID("SldWorks.Application")

sw = connect()
swd = win32com.client.dynamic.Dispatch("SldWorks.Application")


def cop(o, m):
    v = getattr(o, m)
    return v() if callable(v) else v


def open_asm(path):
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    swd.OpenDoc6(path, 2, 3, "", errs, warns)
    return M.ModelDoc2(swd.ActiveDoc._oleobj_), errs.value, warns.value


def wrap(obj, cls="IFeature"):
    try:
        return getattr(M, cls)(obj._oleobj_)
    except Exception:
        return obj


def enumerate_mates(doc):
    mates = []
    feat = wrap(cop(doc, "FirstFeature"))
    guard = 0
    while feat is not None and guard < 6000:
        guard += 1
        try:
            tname = str(cop(feat, "GetTypeName2"))
        except Exception:
            tname = ""
        if "mate" in tname.lower():
            ent = {"ftype": tname, "name": str(cop(feat, "Name"))}
            # the Mates folder holds individual mates as SUB-features
            try:
                sub = wrap(cop(feat, "GetFirstSubFeature"))
                while sub is not None:
                    sub = wrap(sub)
                    try:
                        stn = str(cop(sub, "GetTypeName2"))
                    except Exception:
                        stn = ""
                    ment = {"ftype": stn, "name": str(cop(sub, "Name"))}
                    try:
                        fd = cop(sub, "GetDefinition")
                        if fd is not None:
                            try:
                                ment["mtype"] = int(fd.Type)
                            except Exception:
                                pass
                            n = fd.GetMateEntityCount()
                            ents2 = []
                            for i in range(n):
                                try:
                                    cmp_ = fd.MateEntity(i).ReferenceComponent
                                    ents2.append(str(cop(cmp_, "Name2")))
                                except Exception:
                                    ents2.append("?")
                            ment["components"] = ents2
                    except Exception:
                        pass
                    if ment.get("ftype"):
                        mates.append(ment)
                    try:
                        sub = wrap(cop(sub, "GetNextSubFeature"))
                    except Exception:
                        sub = None
            except Exception:
                pass
            try:
                fd = cop(feat, "GetDefinition")
                if fd is not None:
                    try:
                        ent["mtype"] = int(fd.Type)
                    except Exception:
                        pass
                    try:
                        n = fd.GetMateEntityCount()
                        ents = []
                        for i in range(n):
                            me = fd.MateEntity(i)
                            try:
                                cmp_ = me.ReferenceComponent
                                c2 = cop(cmp_, "Name2")
                                ents.append(str(c2))
                            except Exception:
                                ents.append("?")
                        ent["components"] = ents
                    except Exception as ex:
                        ent["ent_err"] = str(ex)
            except Exception:
                pass
            mates.append(ent)
        try:
            feat = wrap(cop(feat, "GetNextFeature"))
        except Exception:
            break
    return mates


def comp_boxes(swd):
    boxes = {}
    adoc = M.AssemblyDoc(swd.ActiveDoc._oleobj_)
    comps = adoc.GetComponents(False)
    for c in comps:
        try:
            name = str(cop(c, "Name2"))
        except Exception:
            continue
        box = None
        try:
            ce = wrap(c, "IComponent2")
            for style in (1, 0, 2, 3):
                try:
                    b = ce.GetBox(style)
                    if b:
                        box = [float(x) for x in b]
                        break
                except Exception:
                    continue
            if box is None:
                # body bounding box route
                bodies = c.GetBodies2(0)
                if bodies:
                    bb = wrap(bodies[0], "IBody2").GetBoundingBox("")
                    box = [float(x) for x in bb]
        except Exception:
            pass
        boxes[name] = box
    return boxes


log = {}

# ---- FOOT
doc, errs, warns = open_asm(FOOT)
log["foot_open"] = {"errors": errs, "warnings": warns}
log["foot_mates"] = enumerate_mates(doc)
log["foot_boxes"] = comp_boxes(swd)
try:
    log["foot_assembly_box"] = [float(x) for x in doc.GetBox(1)]
except Exception as ex:
    log["foot_assembly_box"] = f"ERR {ex}"
swd.QuitDoc(str(cop(doc, "GetTitle")))

# ---- LEG
doc, errs, warns = open_asm(LEG)
log["leg_open"] = {"errors": errs, "warnings": warns}
log["leg_mates_count_only"] = len(enumerate_mates(doc))
log["leg_boxes"] = comp_boxes(swd)
try:
    log["leg_assembly_box"] = [float(x) for x in doc.GetBox(1)]
except Exception as ex:
    log["leg_assembly_box"] = f"ERR {ex}"
swd.QuitDoc(str(cop(doc, "GetTitle")))

with open(OUT, "w") as f:
    json.dump(log, f, indent=1)
print("foot mates:", len(log["foot_mates"]))
for mte in log["foot_mates"][:30]:
    print("  ", mte.get("name"), mte.get("mtype"), mte.get("components"))
print("foot box sample:", {k: v for k, v in list(log["foot_boxes"].items())[:3]})
print("log:", OUT)
