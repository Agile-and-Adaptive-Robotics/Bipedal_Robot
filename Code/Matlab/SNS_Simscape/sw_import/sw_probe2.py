# sw_probe2.py - typed COM probe of 09_BA_003.SLDASM.
# SW's Application dispatch does not expose type info (EnsureDispatch/CastTo
# fail), so we import the makepy-generated module directly and wrap raw
# _oleobj_ in its interface classes. Read-only, never saves.
import win32com.client
import importlib
import pythoncom
import sys, os, json

ASM = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.SLDASM"
HERE = os.path.dirname(os.path.abspath(__file__))

pythoncom.CoInitialize()
# locate the generated SW API module (the real API typelib: 83A33D31-27C5-11CE)
import win32com.gen_py
gencache_dir = os.path.dirname(win32com.gen_py.__file__)
swmod_name = [f[:-3] for f in os.listdir(gencache_dir)
              if f.endswith(".py") and f.startswith("83A33D31")]
if not swmod_name:
    print("makepy module not found - run makepy on sldworks.tlb first")
    sys.exit(1)
mod = importlib.import_module("win32com.gen_py.%s" % swmod_name[0])
print("typed module:", swmod_name[0])

def typed(obj, cls_names):
    if not isinstance(cls_names, (list, tuple)):
        cls_names = [cls_names]
    for cn in cls_names:
        cls = getattr(mod, cn, None)
        if cls is not None and isinstance(obj, cls):
            return obj
    for cn in cls_names:
        cls = getattr(mod, cn, None)
        if cls is not None:
            return cls(obj._oleobj_)
    raise AttributeError("none of %s in generated module" % cls_names)

sw = win32com.client.Dispatch("SLDWorks.Application")
try:
    sw = typed(sw, "ISldWorks")
except Exception as e:
    print("ISldWorks wrap failed (%s), continuing dynamic" % e)
sw.Visible = False
print("SW:", sw.RevisionNumber())

errs, warns = 0, 0
doc = sw.OpenDoc6(ASM, 2, 0, "", errs, warns)   # 2=assembly, 0=read-write (byref entity calls refused on read-only docs)
if isinstance(doc, tuple):                       # typed byref outputs come back in a tuple
    doc, errs, warns = (list(doc) + [None, 0, 0])[:3]
if doc is None:
    print("FAILED to open, errs=", errs, "warns=", warns)
    sys.exit(1)
doc = typed(doc, ["IModelDoc2", "ModelDoc2"])
print("opened:", doc.GetTitle())
docAsm = typed(doc, "IAssemblyDoc")   # GetComponents
docFeat = doc                          # IModelDoc2: FirstFeature traversal

comps = docAsm.GetComponents(False)
out = {"components": [], "mates": []}
for c in comps:
    c = typed(c, "IComponent2")
    t = c.Transform2
    a = t.ArrayData if t is not None else None
    out["components"].append({
        "name": c.Name2,
        "fixed": c.IsFixed(),
        "transform": [round(float(x), 10) for x in a] if a is not None else None,
        "path": c.GetPathName()})
print("components:", len(out["components"]))
for c in out["components"]:
    print("  ", c["name"], "FIXED" if c["fixed"] else "", "T" if c["transform"] else "-")

# mates: each IComponent2 knows its mates (IMate = the dispatch-friendly class)
seen = {}
for c in out["components"]:
    comp = docAsm.GetComponentByName(c["name"])
    comp = typed(comp, "IComponent2")
    try:
        ml = comp.GetMates()
    except Exception as e:
        print("  GetMates failed on %s: %s" % (c["name"], e))
        continue
    if ml is None:
        continue
    for m in ml:
        m = typed(m, "IMate")            # GetEntity lives here
        mtype = typed(m, "IMate2").Type  # Type lives here
        mDyn = win32com.client.Dispatch(m._oleobj_)  # byref args need VARIANTs
        ent = []
        for i in (0, 1):
            try:
                e = m.IGetEntity(i)          # plain int; doc is read-write now
            except Exception as ex:
                print("      (IGetEntity(%d) failed: %s)" % (i, ex))
                continue
            if e is None:
                print("      (IGetEntity(%d) -> None)" % i)
                continue
            if e is None:
                print("      (GetEntity(%d) -> None)" % i)
                continue
            e = typed(e, "IEntity")
            oc = e.GetComponent()
            oc = typed(oc, "IComponent2") if oc is not None else None
            refinfo = None
            try:
                face = typed(e, "IFace2")
                surf = typed(face.GetSurface(), "ISurface")
                if surf.IsCylinder():
                    cp = surf.CylinderParams
                    refinfo = {"cyl_axis": [round(float(x), 6) for x in cp[3:6]],
                               "cyl_pt": [round(float(x), 6) for x in cp[0:3]],
                               "r": round(float(cp[6]), 6)}
            except Exception:
                refinfo = None
            ent.append({"etype": e.GetType(),
                        "comp": oc.Name2 if oc is not None else None,
                        "ref": refinfo})
        if not ent:
            continue
        key = (mtype, tuple(sorted(str(e.get("comp")) + str(e.get("etype")) for e in ent)))
        if key in seen:
            continue
        seen[key] = {"type": mtype, "entities": ent}
        print("  mate type %s" % mtype)
        for e in ent:
            print("      ", e)
out["mates"] = list(seen.values())
print("unique mates:", len(out["mates"]))

with open(os.path.join(HERE, "sw_probe2_out.json"), "w") as f:
    json.dump(out, f, indent=1, default=str)
print("wrote sw_probe2_out.json")

sw.CloseDoc(doc.GetTitle())
print("PROBE2 DONE")
