# sw_massprop_probe_20260926.py — volume via IModelDocExtension::CreateMassProperty
# + assigned material name per part (cross-check: V = mass_xml / rho_material)
import sys, os, json
sys.path.insert(0, r"C:\Users\Ben\.zcode\skills\solidworks\scripts")
from sw_session import connect
from win32com.client import constants, CastTo

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
    entry = {}
    # material name (early-bound IModelDoc2 has GetMaterialPropertyName2? try both)
    mat = ""
    try:
        r = doc.GetMaterialPropertyName2("", "")
        if isinstance(r, tuple):
            mat = str(r[0])
        else:
            mat = str(r)
    except Exception as e:
        mat = "<%s>" % e.__class__.__name__
    entry["material"] = mat
    # volume via Extension mass property
    vol = None
    try:
        mprop = doc.Extension.CreateMassProperty()
        vol = float(mprop.Volume)
    except Exception as e:
        entry["massprop_error"] = repr(e)[:120]
    entry["volume_m3"] = vol
    try:
        pdoc = CastTo(doc, "IPartDoc")
        entry["bbox_m"] = [float(x) for x in pdoc.GetPartBox(True)]
    except Exception:
        pass
    out[name] = entry
    print("%-16s mat=%-28s vol=%s" % (name, mat, ("%.9f" % vol) if vol is not None else "FAILED"))

dst = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev\part_volumes_20260926.json"
json.dump(out, open(dst, "w"), indent=1)
print("wrote", dst)
