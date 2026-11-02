#!/usr/bin/env python3
"""AH001 repair - SW COM probe (laptop DESKTOP-5Q16KE9, 2026-11-01).

1. Connect to SOLIDWORKS, report version + add-ins (looking for Simscape
   Multibody Link).
2. Open the ORIGINAL 10_AH_001.SLDASM read-only-ish (silent open, no save),
   report open errors, component count, unresolved/suppressed components,
   total mass (if obtainable), and mate inventory summary.
3. Pack and Go the assembly into the repair folder (cad/), remapping all
   references so the working copy is self-contained.
4. Reopen the packed copy and verify component count + mass match.
"""
import json
import math
import os
import sys

sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect  # noqa: E402

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
SRC_ASM = os.path.join(REPO, r"Solid_Models\00_BipedalRobot_Redesign\10_AH_001.SLDASM")
REPAIR = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
CAD_DIR = os.path.join(REPAIR, "cad")
LOG = os.path.join(REPAIR, "logs", "sw_probe_assembly_20261101.json")
os.makedirs(CAD_DIR, exist_ok=True)
os.makedirs(os.path.dirname(LOG), exist_ok=True)

report = {"source": SRC_ASM, "exists": os.path.exists(SRC_ASM)}

sw = connect()
rev = sw.RevisionNumber
report["sw_revision"] = rev() if callable(rev) else rev
report["sw_visible"] = bool(sw.Visible)

# --- add-in inventory
try:
    names = []
    for i in range(1, 50):
        try:
            t, name, loaded, active, _, _, _, _ = sw.GetAddInInfo(i)
            if t == 0 and name == "":
                break
            names.append({"type": t, "name": name, "loaded": bool(loaded), "active": bool(active)})
        except Exception:
            break
    report["addins"] = names
except Exception as e:
    report["addins_error"] = str(e)

sw_open_doc_options_silent = 1     # swOpenDocOptions_Silent
sw_open_doc_options_readonly = 2   # swOpenDocOptions_ReadOnly
swDocASSEMBLY = 2

def open_asm(path, readonly=True):
    from win32com.client import VARIANT
    import win32com.client.dynamic
    import pythoncom
    # early-bound makepy wrapper rejects byref VARIANTs (known trap) -> dynamic
    swd = win32com.client.dynamic.Dispatch("SldWorks.Application")
    opts = sw_open_doc_options_silent | (sw_open_doc_options_readonly if readonly else 0)
    errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
    doc = swd.OpenDoc6(path, swDocASSEMBLY, opts, "", errs, warns)
    return doc, errs.value, warns.value

def asm_inventory(doc, label):
    inv = {"label": label}
    comps = doc.GetComponents(False)  # all, incl. suppressed
    comps = [c for c in comps if c is not None]
    inv["n_components"] = len(comps)
    unresolved, suppressed, paths = 0, 0, []
    for c in comps:
        try:
            p = c.GetPathName() or ""
        except Exception:
            p = ""
        if not p:
            unresolved += 1
        else:
            paths.append(p)
        try:
            # 1 = swComponentSuppressed
            if c.Suppress == 1:
                suppressed += 1
        except Exception:
            pass
    inv["n_unresolved"] = unresolved
    inv["n_suppressed"] = suppressed
    # mate summary
    try:
        mates = doc.GetMates()
        inv["n_mates"] = len(mates) if mates else 0
    except Exception as e:
        inv["n_mates"] = f"err {e}"
    # total mass
    try:
        mp = doc.Extension.CreateMassProperty()
        mp.UseSystemUnits = False
        m = mp.Mass
        com = mp.CenterOfMass
        inv["total_mass_kg"] = m
        inv["com_m"] = list(com)
    except Exception as e:
        inv["mass_error"] = str(e)
    return inv, paths

def call_or_prop(obj, name):
    v = getattr(obj, name)
    return v() if callable(v) else v

doc, errs, warns = open_asm(SRC_ASM, readonly=True)
if doc is None:
    report["open"] = {"errors": errs, "warnings": warns, "doc": None}
else:
    # use the EARLY-bound doc for inventory (makepy knows GetPathName etc.)
    doc = sw.ActiveDoc
    inv, paths = asm_inventory(doc, "original")
    report["open"] = {"errors": errs, "warnings": warns, "doc": call_or_prop(doc, "GetTitle")}
    report["inventory"] = inv
    report["component_paths"] = paths

# Multibody Link add-in handle check (direct, not via GetAddInInfo)
link_progs = ["SimscapeMultibodyLink.AddIn", "SimscapeMultibodyLink.App",
              "SimMechanicsLink.AddIn", "SW2URDF.SW2URDF"]
report["addin_handles"] = {}
for prog in link_progs:
    try:
        h = sw.GetAddInObject(prog)
        report["addin_handles"][prog] = "none" if h is None else str(type(h))
    except Exception as e:
        report["addin_handles"][prog] = f"err {e}"

with open(LOG, "w") as f:
    json.dump(report, f, indent=1)

print(json.dumps({k: v for k, v in report.items() if k != "component_paths"}, indent=1))
print("component_paths written to log:", LOG)
