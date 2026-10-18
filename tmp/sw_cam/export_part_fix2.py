"""Part export with instance-readiness retries (hidden SW instance may be booting)."""
import os, json, time, traceback
import pythoncom
from win32com.client import VARIANT
from win32com.client import dynamic as dyn

KA = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly"
PART = os.path.join(KA, "05_01_TI_R_008.SLDPRT")
OUT = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008"

clsid = pythoncom.MakeIID("{6AF263BB-EB9F-4176-89E9-4F892EB0CA3D}")

def get_sw(timeout=180):
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            disp = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER,
                                              pythoncom.IID_IDispatch)
            sw = dyn.Dispatch(disp)
            if hasattr(sw, "OpenDoc6"):
                return sw
            print("  got handle without OpenDoc6 yet; waiting...")
        except Exception as e:
            print("  CoCreate:", e)
        time.sleep(10)
    raise RuntimeError("SW never became ready")

report = {"exports": [], "bbox_mm": None, "material": None, "errors": [], "closed": []}
sw = get_sw()
print("SW ready; visible:", sw.Visible)

def title_of(doc):
    t = doc.GetTitle
    return t() if callable(t) else t

def was_open(path):
    try:
        return sw.GetOpenDocumentByName(path) is not None
    except Exception:
        return False

try:
    pre = was_open(PART)
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warn = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    part = sw.OpenDoc6(PART, 1, 1, "", errs, warn)   # swDocPART, silent
    if isinstance(part, tuple):
        part = part[0]
    if part is None:
        raise RuntimeError(f"open failed errors={errs.value}")
    time.sleep(2)
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
try:
    with open(os.path.join(OUT, "export_report.json")) as f:
        old = json.load(f)
except Exception:
    old = {}
old["pass2"] = report
if report["bbox_mm"]:
    old["bbox_mm"] = report["bbox_mm"]
if report["material"]:
    old["material"] = report["material"]
old.setdefault("exports", [])
old["exports"] += report["exports"]
with open(os.path.join(OUT, "export_report.json"), "w") as f:
    json.dump(old, f, indent=2)
