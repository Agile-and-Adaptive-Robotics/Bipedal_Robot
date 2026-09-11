# sw_probe.py - read-only probe of 09_BA_003.SLDASM via COM.
# pywin32 dynamic dispatch quirk: zero-arg SolidWorks "methods" come back as
# PROPERTIES. zcall handles both. Opens READ-ONLY, never saves.
import win32com.client
from win32com.client import gencache
import pythoncom
import sys, os, json

ASM = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.SLDASM"
HERE = os.path.dirname(os.path.abspath(__file__))


def zcall(o, name, *args):
    attr = getattr(o, name)
    if callable(attr):
        return attr(*args)
    return attr


def items(coll):
    """SW returns either a collection (Item/Count) or a plain tuple/list."""
    if isinstance(coll, (tuple, list)):
        for x in coll:
            yield x
        return
    n = coll.Count
    for i in range(1, n + 1):
        yield zcall(coll, 'Item', i)


pythoncom.CoInitialize()
sw = win32com.client.Dispatch("SLDWorks.Application")
sw.Visible = False
print("SW version:", getattr(sw, 'RevisionNumber'))

errs = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
warns = win32com.client.VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
doc = sw.OpenDoc6(ASM, 2, 1, "", errs, warns)   # 2=asm, 1=read-only
if doc is None:
    print("FAILED to open, errs=", errs.value, "warns=", warns.value)
    sys.exit(1)
print("opened:", zcall(doc, 'GetTitle'))

comps = zcall(doc, 'GetComponents', False)
out = {"components": [], "mates": []}
for c in items(comps):
    tr = None
    try:
        t = zcall(c, 'Transform2')
        if t is not None:
            tr = [round(float(x), 9) for x in t.ArrayData]
    except Exception:
        pass
    out["components"].append({
        "name": c.Name2,
        "fixed": zcall(c, 'IsFixed'),
        "transform": tr,
        "path": zcall(c, 'GetPathName')})
print("components:", len(out["components"]))

with open(os.path.join(HERE, "sw_probe_out.json"), "w") as f:
    json.dump(out, f, indent=1, default=str)
print("wrote components to sw_probe_out.json")

try:
    mates = zcall(doc, 'GetMates')
    for m in items(mates):
        ent = []
        try:
            mes = zcall(m, 'GetMateEntities')
            for me in items(mes):
                ent.append({"type": zcall(me, 'MateEntityType')})
        except Exception:
            pass
        try:
            mtype = zcall(m, 'Type')
        except Exception:
            mtype = None
        out["mates"].append({"name": zcall(m, 'Name'), "type": mtype, "entities": ent})
    print("mates:", len(out["mates"]))
except (AttributeError, pythoncom.com_error):
    # GetMates hidden from plain dispatch: traverse MateFolder features instead
    print("GetMates unavailable via dispatch; traversing MateFolder features...")
    feat = zcall(doc, 'FirstFeature')
    while feat is not None:
        try:
            tname = zcall(feat, 'GetTypeName2')
        except Exception:
            tname = ''
        if tname == 'MateFolder':
            sub = zcall(feat, 'GetSubFeatures')
            if sub is not None:
                for m in items(sub):
                    try:
                        out["mates"].append({"name": zcall(m, 'Name'),
                                             "type": zcall(m, 'Type'), "entities": []})
                    except Exception:
                        pass
        feat = zcall(feat, 'GetNextFeature')
    print("mates via features:", len(out["mates"]))

with open(os.path.join(HERE, "sw_probe_out.json"), "w") as f:
    json.dump(out, f, indent=1, default=str)
print("wrote sw_probe_out.json")

sw.CloseDoc(zcall(doc, 'GetTitle'))
print("PROBE DONE")
