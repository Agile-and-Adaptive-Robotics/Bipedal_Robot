import scipy.io as sio
import numpy as np

RES = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def val(v, tag):
    if isinstance(v, np.ndarray) and not v.dtype.names:
        if v.size == 1:
            return f"{v.item():.6g}"
        if v.size <= 8:
            return "[" + ", ".join(f"{x:.6g}" for x in v.ravel()) + "]"
        return f"array{v.shape} min={np.nanmin(v):.6g} max={np.nanmax(v):.6g}"
    return repr(v)[:200]

def dump_fields(struct, fields, label):
    names = struct.dtype.names
    for f in fields:
        if f in names:
            try:
                print(f"  {label}.{f} = {val(struct[f].flat[0] if struct[f].size else struct[f], f)}")
            except Exception as e:
                print(f"  {label}.{f} ERR {e}")

jobs = [
    ("Bif", "Bifemsh_20mm_Result.mat",
     ["fBest", "fG", "obj", "Jseed", "exitflagG", "exitflagP", "cBest", "boundViolation",
      "collisionFeasible", "deltaLmSigned", "maxContractionTravel", "relativeContractionBest"]),
    ("Vas", "Vas_Pam_20mm_Result.mat",
     ["fBest", "fG", "obj", "exitflagG", "cBest", "minMargin", "iWorst", "humanBest",
      "shortfallFraction", "contractionBest"]),
    ("Vas1519", "Vas_Pam_20mm_Result_20260920_1519.mat",
     ["fBest", "exitflagG", "cBest", "XiUsed"]),
    ("BifPul", "Bifemsh_20mm_Result_pulley_20260926_1652.mat",
     ["fBest", "fG", "fRefined", "Jseed", "exitflagG", "exitflagP", "exitRefined",
      "cfgBestScore", "cfgBestF", "cfgBestG", "cfgBestBPA", "cfgBPACounts", "cfgLineParts",
      "nCfg", "nB", "kCfg", "maxEvalsG", "maxEvalsP", "boundViolation", "collisionFeasible",
      "smokeMode", "cBest"]),
    ("BiPullBif", "BiPulley_opensim_bifemsh_r_20260926_1702.mat",
     ["f", "f0", "fBest", "fS", "ok", "pairName", "tauTarget", "tauArmNominal", "T",
      "N1", "N2", "budgetS", "budgetP", "exitflagS", "exitflagP", "isSmoke", "minPair",
      "needPair", "baseXi0", "baseXi1", "baseXi2", "basePressure", "Lrig", "Lchar",
      "summaryFile"]),
]

for tag, fn, fields in jobs:
    print("=" * 70)
    print(tag, "->", fn)
    m = sio.loadmat(f"{RES}\\{fn}", squeeze_me=True, struct_as_record=True)
    top = np.zeros(1, dtype=[(k, "O") for k in m if not k.startswith("__")])
    for k in m:
        if not k.startswith("__"):
            top[k] = m[k]
    dump_fields(top, fields, tag)
    if "ctx" in m:
        ctx = m["ctx"]
        dump_fields(ctx, ["Xi0", "Xi1", "Xi2", "Xi3", "requiredTorqueMargin", "BPAcount",
                          "pulleyBPACount", "tackleLineParts", "pulleyRoutingMode",
                          "pulleyExitIndex", "targetPressure", "KMAX", "torqueScale",
                          "pulleyGainBounds", "optimizePulleyGain"], "ctx")
    if "st" in m:
        dump_fields(m["st"], ["margin", "worstC", "obj", "tauAbs"], "st")
    if "spec" in m:
        dump_fields(m["spec"], ["name", "source", "notes", "fmax", "fmaxSource", "pulley",
                                "cross", "kin"], "spec")
