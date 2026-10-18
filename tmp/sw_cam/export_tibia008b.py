"""Export 05_01_TI_R_008 part + drawing for CAM prep (bootstrap #2).

sw_session.py's gen_py cache is broken on this machine today; this uses raw
CoCreateInstance + win32com dynamic dispatch (verified working), with callable()
guards for members that flip property/method across bindings.
Only closes docs that were not already open before this script ran.
"""
import os, json, traceback
import pythoncom
from win32com.client import VARIANT
from win32com.client import dynamic as dyn

KA = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly"
PART = os.path.join(KA, "05_01_TI_R_008.SLDPRT")
DRW = os.path.join(KA, "05_01_TI_R_008.SLDDRW")
OUT = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"
os.makedirs(OUT, exist_ok=True)

# enum ints (no gencache constants available)
swDocPART, swDocASSEMBLY, swDocDRAWING = 1, 2, 3
swOpenDocOptions_Silent = 1
swSaveAsCurrentVersion, swSaveAsOptions_Silent = 0, 1

clsid = pythoncom.MakeIID("{6AF263BB-EB9F-4176-89E9-4F892EB0CA3D}")
disp = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER,
                                  pythoncom.IID_IDispatch)
sw = dyn.Dispatch(disp)
rev = sw.RevisionNumber
if callable(rev):
    rev = rev()
print("SW handle OK, version:", rev)
print("visible:", sw.Visible)

report = {"part": PART, "drawing": DRW, "sw": rev, "exports": [],
          "bbox_mm": None, "material": None, "errors": [], "notes": []}

def was_open(path):
    try:
        return sw.GetOpenDocumentByName(path) is not None
    except Exception:
        return False

def open_doc(path, dtype):
    pre = was_open(path)
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warn = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.OpenDoc6(path, dtype, swOpenDocOptions_Silent, "", errs, warn)
    # dynamic dispatch may return tuple or single object depending on binding
    if isinstance(doc, tuple):
        doc = doc[0]
    if doc is None:
        raise RuntimeError(f"open failed {path} errors={errs.value}")
    return doc, pre

try:
    drw, pre_drw = open_doc(DRW, swDocDRAWING)
    for ext in (".pdf", ".dxf"):
        p = os.path.join(OUT, "05_01_TI_R_008_drawing" + ext)
        err = drw.SaveAs3(p, swSaveAsCurrentVersion, swSaveAsOptions_Silent)
        report["exports"].append({"file": p, "err": err})
    if not pre_drw:
        sw.QuitDoc(drw.GetTitle())
except Exception as e:
    report["errors"].append("drawing: " + repr(e))
    traceback.print_exc()

try:
    part, pre_part = open_doc(PART, swDocPART)
    bodies = part.GetBodies2(0)
    if isinstance(bodies, tuple):
        body = bodies[0]
    else:
        body = bodies
    box = body.GetBodyBox()
    report["bbox_mm"] = [round(v * 1000.0, 4) for v in box]
    try:
        report["material"] = part.GetMaterialPropertyName2("")
    except Exception as e:
        report["material"] = "unavailable: " + repr(e)
    for ext in (".step", ".stl"):
        p = os.path.join(OUT, "05_01_TI_R_008" + ext)
        err = part.SaveAs3(p, swSaveAsCurrentVersion, swSaveAsOptions_Silent)
        report["exports"].append({"file": p, "err": err})
    if not pre_part:
        sw.QuitDoc(part.GetTitle())
except Exception as e:
    report["errors"].append("part: " + repr(e))
    traceback.print_exc()

print(json.dumps(report, indent=2))
with open(os.path.join(OUT, "export_report.json"), "w") as f:
    json.dump(report, f, indent=2)
