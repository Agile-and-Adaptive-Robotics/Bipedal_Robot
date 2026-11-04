#!/usr/bin/env python3
"""sw_stl_chunked.py -- STL export in SW chunks with restart (the cycling
batch crashed; fresh SW every CHUNK parts). Caches existing STLs."""
import os
import subprocess
import sys

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
REP = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair")
OUT = os.path.join(REP, "cad", "stl")
os.makedirs(OUT, exist_ok=True)

PARTS = [
    r"Solid_Models\00_BipedalRobot_Redesign\02_Pelvis\02_01_PE_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\03_Femur\03_01_FB_R_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\03_Femur\03_02_FH_R_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\03_Femur\03_03_FB_L_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\03_Femur\03_04_FH_L_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\04_Knee\04_01_KT_R_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\04_Knee\04_02_KB_R_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\04_Knee\04_03_KT_L_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\04_Knee\04_04_KB_L_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\04_Knee\04_05_BL_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\04_Knee\04_06_FL_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\05_Tibia\05_01_TI_R_001.SLDPRT",
    r"Solid_Models\00_BipedalRobot_Redesign\05_Tibia\05_02_TI_001.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Ankle Exterior.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Ankle Middle Joint Shell.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Ankle Interior Joint.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Ankle Bottom Joint.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\New foot - Copy - Copy.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Left mid foot.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\New Front foot.SLDPRT",
    r"Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Toe.SLDPRT",
]

CHUNK = 2
inner = r'''
import os, sys
sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect
import win32com.client.dynamic
from win32com.client import VARIANT
import pythoncom
parts = sys.argv[1:]
sw = connect()
swd = win32com.client.dynamic.Dispatch("SldWorks.Application")
for p in parts:
    out = sys.argv[0]
    o = os.path.join(OUT, os.path.splitext(os.path.basename(p))[0] + ".STL")
    if os.path.exists(o) and os.path.getsize(o) > 0:
        continue
    try:
        errs = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        warns = VARIANT(pythoncom.VT_BYREF | pythoncom.VT_I4, 0)
        doc = swd.OpenDoc6(p, 1, 3, "", errs, warns)
        if doc is None:
            print("OPEN-FAIL", os.path.basename(p), flush=True); continue
        doc.SaveAs3(o, 0, 1)
        t = doc.GetTitle
        t = t() if callable(t) else t
        swd.QuitDoc(t)
        print("ok" if os.path.exists(o) and os.path.getsize(o) > 0 else "EMPTY",
              os.path.basename(p), flush=True)
    except Exception as e:
        print("FAIL", os.path.basename(p), str(e)[:60], flush=True)
try:
    sw.ExitApp()
except Exception:
    pass
'''
inner = inner.replace("out = sys.argv[0]", "OUT = r'%s'" % OUT)
inner_path = os.path.join(REP, "scripts", "_sw_chunk_inner.py")
open(inner_path, "w").write(inner)

todo = [os.path.join(REPO, p) for p in PARTS
        if not (os.path.exists(os.path.join(OUT, os.path.splitext(os.path.basename(p))[0] + ".STL"))
                and os.path.getsize(os.path.join(OUT, os.path.splitext(os.path.basename(p))[0] + ".STL")) > 0)]
print("to export:", len(todo))
for i in range(0, len(todo), CHUNK):
    chunk = todo[i:i+CHUNK]
    print("== chunk", i // CHUNK + 1, len(chunk), "parts", flush=True)
    r = subprocess.run([sys.executable, "-u", inner_path] + chunk,
                       capture_output=True, text=True, timeout=900)
    print(r.stdout.strip())
    if r.returncode != 0:
        print("chunk rc", r.returncode, r.stderr.strip()[-200:])
print("done; STLs:", len([f for f in os.listdir(OUT) if f.endswith('.STL')]))
