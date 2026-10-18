"""Query exact geometry of 05_01_TI_R_008 via SW API: circular edges (holes) + planar faces (steps)."""
import os, json, time, math, traceback
import pythoncom
from win32com.client import VARIANT
from win32com.client import dynamic as dyn

PART = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\05_01_TI_R_008.SLDPRT"
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
            print("  waiting for OpenDoc6...")
        except Exception as e:
            print("  CoCreate:", e)
        time.sleep(10)
    raise RuntimeError("SW not ready")

sw = get_sw()

def title_of(doc):
    t = doc.GetTitle
    return t() if callable(t) else t

errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
warn = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
part = sw.OpenDoc6(PART, 1, 1, "", errs, warn)
if isinstance(part, tuple):
    part = part[0]
time.sleep(2)

report = {"circles": [], "planes": [], "bbox_mm": None, "material": None, "errors": []}

try:
    bodies = part.GetBodies2(0, False)
    body = bodies[0] if isinstance(bodies, tuple) else bodies
    box = body.GetBodyBox()
    report["bbox_mm"] = [round(v * 1000.0, 4) for v in box]

    # --- circular edges: exact holes ---
    edges = body.GetEdges()
    if not isinstance(edges, tuple):
        edges = (edges,)
    seen = set()
    for e in edges:
        try:
            curve = e.GetCurve
            curve = curve() if callable(curve) else curve
            is_circ = curve.IsCircle
            is_circ = is_circ() if callable(is_circ) else is_circ
            if not is_circ:
                continue
            cp = curve.CircleParams          # [x,y,z, i,j,k, radius, ...]
            cp = cp() if callable(cp) else cp
            center = tuple(round(v * 1000.0, 4) for v in cp[0:3])
            axis = tuple(round(v, 4) for v in cp[3:6])
            radius = round(cp[6] * 1000.0, 4)
            try:
                s, t0v, e0, t1v = e.GetEndPointsOrDefault()
            except Exception:
                s = None
            key = (center, axis, radius)
            if key in seen:
                continue
            seen.add(key)
            # endpoints give the extent along the circle's plane placement
            if s is not None:
                ptA = tuple(round(v * 1000.0, 3) for v in s[0:3])
                ptB = tuple(round(v * 1000.0, 3) for v in e0[0:3])
            else:
                ptA = ptB = None
            report["circles"].append({"center": center, "axis": axis,
                                      "r_mm": radius, "d_mm": round(2 * radius, 4),
                                      "ptA": ptA, "ptB": ptB})
        except Exception:
            pass

    # --- planar faces: step levels ---
    faces = body.GetFaces()
    if not isinstance(faces, tuple):
        faces = (faces,)
    pseen = {}
    for f in faces:
        try:
            surf = f.GetSurface
            surf = surf() if callable(surf) else surf
            isp = surf.IsPlane
            isp = isp() if callable(isp) else isp
            if not isp:
                continue
            pp = surf.PlaneParams          # [i,j,k (normal), x,y,z (root)]
            pp = pp() if callable(pp) else pp
            n = tuple(round(v, 4) for v in pp[0:3])
            root = tuple(round(v * 1000.0, 3) for v in pp[3:6])
            fb = f.GetBox
            fb = fb() if callable(fb) else fb
            corners = [round(v * 1000.0, 3) for v in fb]
            key = (n, root)
            area_proxy = (corners[3] - corners[0]) * (corners[4] - corners[1])
            if key not in pseen or area_proxy > pseen[key]:
                pseen[key] = area_proxy
        except Exception:
            pass
    report["planes"] = [{"normal": k[0], "root": k[1], "area_proxy": round(v, 1)}
                        for k, v in sorted(pseen.items(), key=lambda kv: -kv[1])]

    try:
        report["material"] = part.GetMaterialPropertyName2("")
    except Exception as e:
        report["material"] = "unavailable: " + repr(e)

    sw.QuitDoc(title_of(part))
except Exception as e:
    report["errors"].append(repr(e))
    traceback.print_exc()

print(json.dumps(report, indent=1))
with open(os.path.join(OUT, "geometry_api.json"), "w") as f:
    json.dump(report, f, indent=1)
