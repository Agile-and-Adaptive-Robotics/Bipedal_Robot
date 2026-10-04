import scipy.io as sio
import numpy as np

RES = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def val(v):
    if isinstance(v, np.ndarray) and not v.dtype.names:
        if v.size == 1:
            return f"{v.item():.6g}"
        if v.size <= 8:
            return "[" + ", ".join(f"{x:.6g}" for x in v.ravel()) + "]"
        return f"array{v.shape} min={np.nanmin(v):.6g} max={np.nanmax(v):.6g}"
    if isinstance(v, np.ndarray) and v.dtype.names:
        return "struct:" + ",".join(v.dtype.names)
    return repr(v)[:180]

jobs = [
    ("Bif", "Bifemsh_20mm_Result.mat",
     dict(top=["fBest", "fG", "obj", "Jseed", "exitflagG", "exitflagP", "cBest",
               "boundViolation", "collisionFeasible", "deltaLmSigned",
               "maxContractionTravel", "relativeContractionBest", "torqueMarginFraction"],
          ctx=["Xi0", "Xi1", "Xi2", "Xi3", "requiredTorqueMargin", "BPAcount",
               "targetPressure", "KMAX", "torqueScale"]),
    ("Vas", "Vas_Pam_20mm_Result.mat",
     dict(top=["fBest", "fG", "obj", "exitflagG", "cBest", "minMargin", "iWorst",
               "humanBest", "shortfallFraction", "contractionBest", "torqueMargin"],
          ctx=["Xi0", "Xi1", "Xi2", "Xi3", "requiredTorqueMargin", "BPAcount",
               "targetPressure", "KMAX", "torqueScale"])),
    ("Vas1519", "Vas_Pam_20mm_Result_20260920_1519.mat",
     dict(top=["fBest", "exitflagG", "cBest", "XiUsed"],
          ctx=["Xi0", "Xi1", "Xi2", "Xi3"])),
    ("BifPul", "Bifemsh_20mm_Result_pulley_20260926_1652.mat",
     dict(top=["fBest", "fG", "fRefined", "Jseed", "exitflagG", "exitflagP",
               "cfgBestScore", "cfgBestF", "cfgBestG", "cfgBestBPA", "cfgBPACounts",
               "cfgLineParts", "nCfg", "nB", "kCfg", "maxEvalsG", "maxEvalsP",
               "boundViolation", "collisionFeasible", "smokeMode", "cBest"],
          ctx=["Xi0", "Xi1", "Xi2", "Xi3", "requiredTorqueMargin", "BPAcount",
               "pulleyBPACount", "tackleLineParts", "pulleyRoutingMode",
               "pulleyExitIndex", "pulleyGainBounds", "optimizePulleyGain"])),
    ("BiPullBif", "BiPulley_opensim_bifemsh_r_20260926_1702.mat",
     dict(top=["f", "f0", "fBest", "fS", "ok", "pairName", "tauTarget",
               "tauArmNominal", "T", "N1", "N2", "budgetS", "budgetP", "exitflagS",
               "exitflagP", "isSmoke", "minPair", "needPair", "baseXi0", "baseXi1",
               "baseXi2", "basePressure", "Lrig", "Lchar", "summaryFile"],
          ctx=[], st=["margin", "worstC", "obj", "tauAbs"],
          spec=["name", "source", "notes", "fmax", "fmaxSource", "pulley", "cross", "kin"])),
]

for tag, fn, spec in jobs:
    print("=" * 70)
    print(tag, "->", fn)
    m = sio.loadmat(f"{RES}\\{fn}", squeeze_me=True, struct_as_record=True)
    for k in spec.get("top", []):
        if k in m:
            try:
                print(f"  {k} = {val(m[k])}")
            except Exception as e:
                print(f"  {k} ERR {type(e).__name__}: {e}")
        else:
            print(f"  {k} <absent>")
    for sub, flds in (("ctx", spec.get("ctx", [])), ("st", spec.get("st", [])),
                      ("spec", spec.get("spec", []))):
        if sub in m and flds:
            s = m[sub]
            names = s.dtype.names if hasattr(s, "dtype") and s.dtype.names else []
            for f in flds:
                if f in (names or []):
                    try:
                        print(f"  {sub}.{f} = {val(s[f])}")
                    except Exception as e:
                        print(f"  {sub}.{f} ERR {e}")
                else:
                    print(f"  {sub}.{f} <absent>")
