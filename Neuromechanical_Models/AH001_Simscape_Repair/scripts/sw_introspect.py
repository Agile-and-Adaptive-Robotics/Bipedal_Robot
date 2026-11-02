#!/usr/bin/env python3
"""Introspect the makepy v33 wrappers for the members we need."""
import sys
sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402
import win32com.client  # noqa: E402

sw = connect()
doc = sw.ActiveDoc
print("doc type:", type(doc).__name__, "title:", doc.GetTitle() if callable(doc.GetTitle) else doc.GetTitle)

def show(obj, label, pat=None):
    try:
        ms = [m for m in dir(obj) if not m.startswith("_") and (pat is None or pat.lower() in m.lower())]
        print(f"{label}: {len(ms)} members", ms[:40] if pat else ms[:25])
    except Exception as e:
        print(label, "ERR", e)

show(doc, "ModelDoc2 (head)")
show(doc, "ModelDoc2 mass/mates", "mate")
show(doc, "ModelDoc2 mass", "mass")
ext = doc.Extension
print("Extension type:", type(ext).__name__)
show(ext, "Extension mass", "mass")
show(doc, "GetComponents-ish", "component")

# try assembly cast
for cast in ("AssemblyDoc", "IAssemblyDoc"):
    try:
        asm = win32com.client.CastTo(doc, cast)
        print("CastTo OK:", cast, type(asm).__name__)
        show(asm, f"{cast} mates", "mate")
        break
    except Exception as e:
        print("CastTo fail:", cast, e)

# component path via call_or_prop
comps = doc.GetComponents(False)
print("n comps:", len(comps))
c0 = comps[0]
print("comp0 type:", type(c0).__name__)
for member in ("GetPathName", "PathName", "Name2", "Suppress", "SuppressionState", "IsResolved"):
    try:
        v = getattr(c0, member)
        try:
            v = v()
        except TypeError:
            pass
        print(f"comp0.{member} =", str(v)[:90])
    except Exception as e:
        print(f"comp0.{member} ERR {e}")

# mass properties the old way
try:
    res = doc.GetMassProperties2(1)
    print("GetMassProperties2:", str(res)[:200])
except Exception as e:
    print("GetMassProperties2 ERR", e)
