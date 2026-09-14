"""Open 09_BA_003.SLDASM and drive the Simscape Multibody Link add-in's
export methods directly via COM (no GUI dialog).

Tries, in order:
  1. SaveAsSimMechXml(filename)  -- headless XML export (hope: Multibody XML)
  2. reports what landed on disk so we can identify the XML flavor.
"""
import os
import sys
import pythoncom
from win32com.client import Dispatch, VARIANT

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import constants

ASM = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.SLDASM"
OUT = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\sw_import\smlink_export_test.xml"

sw = Dispatch("SldWorks.Application")  # dynamic dispatch
sw.Visible = True
errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
doc = sw.OpenDoc6(ASM, constants.swDocASSEMBLY, constants.swOpenDocOptions_Silent,
                  "", errs, warns)
if doc is None:
    raise SystemExit(f"OpenDoc6 failed errs={errs.value} warns={warns.value}")
def zcall(x, *a):
    """SolidWorks zero-arg methods come back as pywin32 PROPERTIES."""
    return x(*a) if callable(x) else x

print("opened:", zcall(doc.GetTitle))

DLL = r"C:\Program Files\MATLAB\R2025b\bin\win64\cl_sldwks2sm.dll"
loaded = sw.LoadAddIn(DLL)
print("LoadAddIn ->", loaded)
addin = sw.GetAddInObject("sldwks2sm.Sldwks2smApp")
print("addin object:", addin)
if addin is None:
    raise SystemExit("add-in object unavailable after LoadAddIn")

if os.path.exists(OUT):
    os.remove(OUT)

try:
    r = addin.SaveAsSimMechXml(OUT)
    print("SaveAsSimMechXml ->", r)
except Exception as e:
    print("SaveAsSimMechXml raised:", repr(e))

if os.path.exists(OUT):
    sz = os.path.getsize(OUT)
    print("XML written:", OUT, sz, "bytes")
    with open(OUT, "r", errors="replace") as f:
        head = f.read(600)
    print("--- head ---")
    print(head)
else:
    print("no XML produced")
