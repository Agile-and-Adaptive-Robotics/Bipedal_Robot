"""Fix pass: export part STEP/STL + bbox (GetBodies2 needs 2 args), title guard."""
import os, json, traceback
import pythoncom
from win32com.client import VARIANT
from win32com.client import dynamic as dyn

KA = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly"
PART = os.path.join(KA, "05_01_TI_R_008.SLDPRT")
DRW = os.path.join(KA, "05_01_TI_R_008.SLDDRW")
OUT = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"

swDocPART, swDocDRAWING = 1, 3
clsid = pythoncom.MakeIID("{6AF263BB-EB9F-4176-89E9-4F892EB0CA3D}")
disp = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER, pythoncom.IID_IDispatch)
sw = dyn.Dispatch(disp)

def title_of(doc):
    t = doc.GetTitle
    return t() if callable(t) else t

def was_open(path):
    try:
        return sw.GetOpenDocumentByName(path) is not None
    except Exception:
        return False

report = {"exports": [], "bbox_mm": None, "material": None, "closed": [], "errors": []}

# close the drawing left open by the previous run (it's fully exported)
try:
    d = sw.GetOpenDocumentByName(DRW)
    if d is not None:
        t = title_of(d)
        sw.QuitDoc(t)
        report["closed"].append(t)
except Exception as e:
    report["errors"].append("close drawing: " + repr(e))

try:
    pre = was_open(PART)
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warn = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    part = sw.OpenDoc6(PART, swDocPART, 1, "", errs, warn)
    if isinstance(part, tuple):
        part = part[0]
    if part is None:
        raise RuntimeError(f"open failed errors={errs.value}")
    bodies = part.GetBodies2(0, False)
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
        err = part.SaveAs3(p, 0, 1)
        report["exports"].append({"file": p, "err": err})
    if not pre:
        sw.QuitDoc(title_of(part))
        report["closed"].append("part")
except Exception as e:
    report["errors"].append("part: " + repr(e))
    traceback.print_exc()

print(json.dumps(report, indent=2))
old = {}
try:
    with open(os.path.join(OUT, "export_report.json")) as f:
        old = json.load(f)
except Exception:
    pass
old.update({k: v for k, v in report.items() if v})
old["pass2"] = report
with open(os.path.join(OUT, "export_report.json"), "w") as f:
    json.dump(old, f, indent=2)
