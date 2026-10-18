"""Export 05_01_TI_R_008 part + drawing for CAM prep.

Attaches to the running SLDWORKS instance (Ben is live in the GUI). Exports:
  drawing -> PDF (then rasterized separately) + DXF (2-axis CAM reference)
  part    -> STEP + STL + exact bounding box + material name
Only closes docs that were not already open before this script ran.
"""
import os, sys, json

sys.path.insert(0, r"C:\Users\Ben Bolen\.agents\skills\solidworks\scripts")
from sw_session import connect, constants
import pythoncom
from win32com.client import VARIANT

KA = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly"
PART = os.path.join(KA, "05_01_TI_R_008.SLDPRT")
DRW = os.path.join(KA, "05_01_TI_R_008.SLDDRW")
OUT = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"
os.makedirs(OUT, exist_ok=True)

sw = connect()
print("attached; SW version:", sw.RevisionNumber())

def was_open(path):
    try:
        return sw.GetOpenDocumentByName(path) is not None
    except Exception:
        return False

def open_doc(path, dtype):
    pre = was_open(path)
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = sw.OpenDoc6(path, dtype, 0, "", errs)[0]   # non-silent: native files open clean
    if doc is None:
        raise RuntimeError(f"open failed {path} errors={errs.value}")
    return doc, pre

report = {"part": PART, "drawing": DRW, "exports": [], "bbox_mm": None,
          "material": None, "errors": []}

# ---- drawing ----
try:
    drw, pre_drw = open_doc(DRW, constants.swDocDRAWING)
    for ext in (".pdf", ".dxf"):
        p = os.path.join(OUT, "05_01_TI_R_008_drawing" + ext)
        err = drw.SaveAs3(p, constants.swSaveAsCurrentVersion, constants.swSaveAsOptions_Silent)
        report["exports"].append({"file": p, "err": err})
    if not pre_drw:
        sw.QuitDoc(drw.GetTitle())
except Exception as e:
    report["errors"].append(f"drawing: {e}")

# ---- part ----
try:
    part, pre_part = open_doc(PART, constants.swDocPART)
    body = part.GetBodies2(0)[0] if hasattr(part.GetBodies2(0), "__getitem__") else part.GetBodies2(0)
    box = body.GetBodyBox()   # meters: xmin,ymin,zmin,xmax,ymax,zmax
    report["bbox_mm"] = [round(v * 1000.0, 4) for v in box]
    try:
        report["material"] = part.GetMaterialPropertyName2("")
    except Exception as e:
        report["material"] = f"unavailable: {e}"
    for ext, name in ((".step", "05_01_TI_R_008"), (".stl", "05_01_TI_R_008")):
        p = os.path.join(OUT, name + ext)
        err = part.SaveAs3(p, constants.swSaveAsCurrentVersion, constants.swSaveAsOptions_Silent)
        report["exports"].append({"file": p, "err": err})
    if not pre_part:
        sw.QuitDoc(part.GetTitle())
except Exception as e:
    report["errors"].append(f"part: {e}")

print(json.dumps(report, indent=2))
with open(os.path.join(OUT, "export_report.json"), "w") as f:
    json.dump(report, f, indent=2)
