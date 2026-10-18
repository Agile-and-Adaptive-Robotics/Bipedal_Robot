import pythoncom, sys
print("python:", sys.version)
print("pythoncom file:", getattr(pythoncom, "__file__", "?"))
print("CoCreateInstance type:", type(pythoncom.CoCreateInstance))
print("GetActiveObject type:", type(pythoncom.GetActiveObject))

import win32com.client.dynamic as dyn
clsid = pythoncom.MakeIID("{6AF263BB-EB9F-4176-89E9-4F892EB0CA3D}")
try:
    disp = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER, pythoncom.IID_IDispatch)
    sw = dyn.Dispatch(disp)
    print("HANDLE OK; version:", sw.RevisionNumber())
    ad = sw.ActiveDoc
    print("active doc:", ad.GetTitle() if ad else "(none)")
except Exception as e:
    import traceback; traceback.print_exc()
