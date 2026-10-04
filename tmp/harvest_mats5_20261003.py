import scipy.io as sio
import numpy as np

RES = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def val(v):
    if isinstance(v, np.ndarray) and not v.dtype.names:
        if v.size == 1:
            return "%.6g" % v.item()
        if v.size <= 8:
            return "[" + ", ".join("%.6g" % x for x in v.ravel()) + "]"
        return "array%s min=%.6g max=%.6g" % (v.shape, np.nanmin(v), np.nanmax(v))
    if isinstance(v, np.ndarray) and v.dtype.names:
        return "struct:" + ",".join(v.dtype.names)
    return repr(v)[:180]

def show(m, prefix, names):
    for k in names:
        if k not in m:
            print("  %s.%s <absent>" % (prefix, k))
            continue
        try:
            print("  %s.%s = %s" % (prefix, k, val(m[k])))
        except Exception as e:
            print("  %s.%s ERR %s" % (prefix, k, e))

def show_sub(m, sub, names):
    if sub not in m or not names:
        return
    s = m[sub]
    names_all = s.dtype.names if hasattr(s, "dtype") and s.dtype.names else []
    for f in names:
        if f not in (names_all or []):
            print("  %s.%s <absent>" % (sub, f))
            continue
        try:
            print("  %s.%s = %s" % (sub, f, val(s[f])))
        except Exception as e:
            print("  %s.%s ERR %s" % (sub, f, e))

CTX_COMMON = ["Xi0", "Xi1", "Xi2", "Xi3", "requiredTorqueMargin", "BPAcount"]

m = sio.loadmat(RES + "\\" + "Bifemsh_20mm_Result.mat", squeeze_me=True, struct_as_record=True)
print("=" * 70); print("Bifemsh_20mm_Result.mat")
show(m, "top", ["fBest", "fG", "obj", "Jseed", "exitflagG", "exitflagP", "cBest",
                "boundViolation", "collisionFeasible", "deltaLmSigned",
                "maxContractionTravel", "relativeContractionBest", "torqueMarginFraction"])
show_sub(m, "ctx", CTX_COMMON + ["targetPressure", "KMAX", "torqueScale"])

m = sio.loadmat(RES + "\\" + "Vas_Pam_20mm_Result.mat", squeeze_me=True, struct_as_record=True)
print("=" * 70); print("Vas_Pam_20mm_Result.mat")
show(m, "top", ["fBest", "fG", "obj", "exitflagG", "cBest", "minMargin", "iWorst",
                "humanBest", "shortfallFraction", "contractionBest", "torqueMargin"])
show_sub(m, "ctx", CTX_COMMON + ["targetPressure", "KMAX", "torqueScale"])

m = sio.loadmat(RES + "\\" + "Vas_Pam_20mm_Result_20260920_1519.mat", squeeze_me=True, struct_as_record=True)
print("=" * 70); print("Vas_Pam_20mm_Result_20260920_1519.mat")
show(m, "top", ["fBest", "exitflagG", "cBest", "XiUsed"])
show_sub(m, "ctx", CTX_COMMON)

m = sio.loadmat(RES + "\\" + "Bifemsh_20mm_Result_pulley_20260926_1652.mat", squeeze_me=True, struct_as_record=True)
print("=" * 70); print("Bifemsh_20mm_Result_pulley_20260926_1652.mat")
show(m, "top", ["fBest", "fG", "fRefined", "Jseed", "exitflagG", "exitflagP",
                "cfgBestScore", "cfgBestF", "cfgBestG", "cfgBestBPA", "cfgBPACounts",
                "cfgLineParts", "nCfg", "nB", "kCfg", "maxEvalsG", "maxEvalsP",
                "boundViolation", "collisionFeasible", "smokeMode", "cBest"])
show_sub(m, "ctx", CTX_COMMON + ["pulleyBPACount", "tackleLineParts", "pulleyRoutingMode",
                                 "pulleyExitIndex", "pulleyGainBounds", "optimizePulleyGain"])

m = sio.loadmat(RES + "\\" + "BiPulley_opensim_bifemsh_r_20260926_1702.mat", squeeze_me=True, struct_as_record=True)
print("=" * 70); print("BiPulley_opensim_bifemsh_r_20260926_1702.mat")
show(m, "top", ["f", "f0", "fBest", "fS", "ok", "pairName", "tauTarget",
                "tauArmNominal", "T", "N1", "N2", "budgetS", "budgetP", "exitflagS",
                "exitflagP", "isSmoke", "minPair", "needPair", "baseXi0", "baseXi1",
                "baseXi2", "basePressure", "Lrig", "Lchar", "summaryFile"])
show_sub(m, "st", ["margin", "worstC", "obj", "tauAbs"])
show_sub(m, "spec", ["name", "source", "notes", "fmax", "fmaxSource", "pulley", "cross", "kin"])
