"""Probe the Simscape Multibody Link add-in (sldwks2sm) COM interface.

Dumps the typelib members of sldwks2sm.Sldwks2smApp and tries to load the
add-in into a (possibly new) SolidWorks session, then fetch the add-in
object via ISldWorks::GetAddInObject.
"""
import sys
import winreg
import pythoncom
from win32com.client import Dispatch

TLB = "{D9AB98CD-176A-4BCB-B178-BB2D0DEE2F92}"
CLSID = "{2666BDBF-5207-4731-9976-13172BEB124F}"
DLL = r"C:\Program Files\MATLAB\R2025b\bin\win64\cl_sldwks2sm.dll"


def dump_typelib():
    """Enumerate typeinfo names + member names of the sldwks2sm typelib."""
    try:
        tlb = pythoncom.LoadTypeLib(DLL)  # typelib embedded in the add-in DLL
    except Exception as e:
        print("LoadTypeLib(DLL) failed:", e)
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, rf"TypeLib\{TLB}") as k:
            ver = winreg.EnumKey(k, 0)
        maj, mnr = (int(x) for x in ver.split("."))
        tlb = pythoncom.LoadRegTypeLib(TLB, maj, mnr, 0)
    print(f"typelib loaded: {tlb.GetTypeInfoCount()} typeinfos")
    for i in range(tlb.GetTypeInfoCount()):
        ti = tlb.GetTypeInfo(i)
        name, doc, ctx, helpfile = tlb.GetDocumentation(i)
        attr = ti.GetTypeAttr()
        kind = attr.typekind  # 3=interface 5=coclass 1=enum
        print(f"\n[{i}] {name} (typekind={kind}, funcs={attr.cFuncs})")
        if kind in (3, 4, 5):  # interface, dispatch, coclass
            for f in range(attr.cFuncs):
                fd = ti.GetFuncDesc(f)
                names = ti.GetNames(fd.memid)
                print(f"    memid=0x{fd.memid & 0xffffffff:x}: {names}")
        elif kind == 1:  # enum
            for v in range(attr.cVars):
                vd = ti.GetVarDesc(v)
                vnames = ti.GetNames(vd.memid, 1)
                print(f"    value {vnames} = {vd._.lpvarValue}")


def sw_connect():
    sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
    from sw_session import connect
    sw = connect()
    rev = sw.RevisionNumber() if callable(sw.RevisionNumber) else sw.RevisionNumber
    print("SW connected, rev:", rev)
    return sw


def load_addin(sw):
    ok = sw.LoadAddIn(DLL)
    print("LoadAddIn ->", ok)
    # try common object names
    for nm in ("sldwks2sm.Sldwks2smApp", "sldwks2sm.Sldwks2smApp.1",
               "Simscape.Multibody.Link", "sldwks2sm"):
        try:
            obj = sw.GetAddInObject(nm)
            print(f"GetAddInObject({nm!r}) ->", obj)
        except Exception as e:
            print(f"GetAddInObject({nm!r}) raised: {e}")


def standalone_create():
    try:
        obj = Dispatch("sldwks2sm.Sldwks2smApp.1")
        print("standalone Dispatch OK:", obj)
        print("dir hints: _methods_", [m for m in dir(obj) if not m.startswith("_")][:40])
    except Exception as e:
        print("standalone Dispatch failed:", e)


if __name__ == "__main__":
    dump_typelib()
    print("\n=== SW session ===")
    sw = sw_connect()
    load_addin(sw)
