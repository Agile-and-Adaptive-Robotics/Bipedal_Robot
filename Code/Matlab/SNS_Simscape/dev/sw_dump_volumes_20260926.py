# sw_dump_volumes_20260926.py — exact solid volume + bbox per humanoid skeleton part
# (Onyx mass rescale + feet/spine/muscle placement). Read-only; attaches to a
# running SOLIDWORKS if present, else launches.
import sys, os, json, math
sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect
from win32com.client import constants

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\00_BipedalRobot_Redesign"
PARTS = [
    "02_Pelvis\\02_01_PE_001.SLDPRT",
    "03_Femur\\03_01_FB_R_001.SLDPRT",
    "03_Femur\\03_02_FH_R_001.SLDPRT",
    "03_Femur\\03_03_FB_L_001.SLDPRT",
    "03_Femur\\03_04_FH_L_001.SLDPRT",
    "04_Knee\\04_01_KT_R_001.SLDPRT",
    "04_Knee\\04_02_KB_R_001.SLDPRT",
    "04_Knee\\04_03_KT_L_001.SLDPRT",
    "04_Knee\\04_04_KB_L_001.SLDPRT",
    "04_Knee\\04_05_BL_001.SLDPRT",
    "04_Knee\\04_06_FL_001.SLDPRT",
    "05_Tibia\\05_01_TI_R_001.SLDPRT",
    "05_Tibia\\05_02_TI_001.SLDPRT",
]

sw = connect()
errs = None

out = {}
for rel in PARTS:
    path = os.path.join(REPO, rel)
    name = os.path.splitext(os.path.basename(rel))[0]
    doc = None
    for pth in (path, os.path.join(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly", os.path.basename(rel))):
        if not os.path.isfile(pth):
            continue
        doc = sw.OpenDoc6(pth, constants.swDocPART, 0, "", 0, 0)
        if isinstance(doc, tuple):
            doc = doc[0]
        if doc is not None:
            break
    if doc is None:
        print("OPEN FAILED:", name)
        continue
    try:
        from win32com.client import CastTo
        pdoc = CastTo(doc, "IPartDoc")
        vol = None
        try:
            bodies = pdoc.GetBodies2(-1, False)  # -1 = all types, NO temporary bodies
        except Exception:
            bodies = None
        if bodies:
            if not isinstance(bodies, tuple):
                bodies = (bodies,)
            vol = 0.0
            nb = 0
            for b in bodies:
                mp = b.GetMassProperties(1000.0)  # SI; vol independent of density
                nb += 1
                if mp is not None:
                    arr = list(mp)
                    print("   body type=%s vol=%.9f" % (b.GetType(), float(arr[0])))
                    vol += float(arr[0])
            print("   (%d bodies)" % nb)
        bb = pdoc.GetPartBox(True)  # meters, part coords
        out[name] = {
            "volume_m3": vol,
            "bbox_m": [float(x) for x in bb] if bb else None,
        }
        print("%-16s vol=%9.6f m3  bbox=%s" % (name, vol if vol else -1,
              ["%.4f" % v for v in bb] if bb else None))
    finally:
        # leave docs open if Ben had them open (lock files present pre-run)
        lock = os.path.join(os.path.dirname(path), "~$" + os.path.basename(rel)[1:] + ".SLDPRT")
        was_open = os.path.isfile(os.path.join(os.path.dirname(path), "~$" + os.path.basename(rel)))
        if not was_open:
            sw.CloseDoc(os.path.basename(path))

dst = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev\part_volumes_20260926.json"
json.dump(out, open(dst, "w"), indent=1)
print("wrote", dst)
