import scipy.io as sio
import numpy as np

RES = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def g(d, k):
    v = d.get(k)
    if v is None:
        return None
    if isinstance(v, np.ndarray):
        if v.dtype.names:
            return v.flat[0]
        if v.size <= 6:
            return v.ravel().tolist()
        return v
    return v

files = {
    "Bif": "Bifemsh_20mm_Result.mat",
    "Vas": "Vas_Pam_20mm_Result.mat",
    "Vas1519": "Vas_Pam_20mm_Result_20260920_1519.mat",
    "BifPul": "Bifemsh_20mm_Result_pulley_20260926_1652.mat",
    "BiPullBif": "BiPulley_opensim_bifemsh_r_20260926_1702.mat",
}

for tag, fn in files.items():
    print("=" * 70)
    print(tag, "->", fn)
    m = sio.loadmat(f"{RES}\\{fn}", squeeze_me=True, struct_as_record=True)
    keys = [k for k in m if not k.startswith("__")]
    print("TOP KEYS:", sorted(keys))
    for k in keys:
        if any(s in k.lower() for s in ("xi", "margin", "jseed", "jbest", "score", "obj", "iter", "pulley", "fail", "ok")) and not hasattr(m[k], "dtype") or (isinstance(m[k], np.ndarray) and not m[k].dtype.names and m[k].size <= 6):
            try:
                print(f"  {k} = {m[k].ravel().tolist()}")
            except Exception:
                pass
    # ctx Xi + required margin
    ctx = m.get("ctx")
    if ctx is not None and hasattr(ctx, "_fieldnames"):
        for f in ("Xi0", "Xi1", "Xi2", "Xi3", "requiredTorqueMargin", "BPAcount",
                  "pulleyBPACount", "tackleLineParts", "pulleyRoutingMode",
                  "pulleyExitIndex", "targetPressure", "KMAX", "torqueScale"):
            if f in ctx._fieldnames:
                v = getattr(ctx, f)
                v = v.tolist() if isinstance(v, np.ndarray) and v.size <= 6 else v
                print(f"  ctx.{f} = {v}")
    # predBest extras
    pb = m.get("predBest")
    if pb is not None and hasattr(pb, "_fieldnames"):
        for f in ("ok", "pulleyConfig", "pulleyGain", "BPAcount", "bpa", "slackCount",
                  "infeasibleCount", "TorqueInsZ", "TorqueZ", "momentArm", "failReason"):
            if f in pb._fieldnames:
                v = getattr(pb, f)
                if isinstance(v, np.ndarray):
                    if v.size > 6:
                        print(f"  predBest.{f}: shape={v.shape} min={np.nanmin(v):.6g} max={np.nanmax(v):.6g}")
                        continue
                    v = v.tolist()
                print(f"  predBest.{f} = {v}")
    st = m.get("st")
    if st is not None and hasattr(st, "_fieldnames"):
        for f in st._fieldnames:
            v = getattr(st, f)
            v = v.tolist() if isinstance(v, np.ndarray) and v.size <= 6 else v
            print(f"  st.{f} = {v}")
    spec = m.get("spec")
    if spec is not None and hasattr(spec, "_fieldnames"):
        for f in ("name", "source", "notes", "fmax", "fmaxSource", "pulley", "cross"):
            if f in spec._fieldnames:
                print(f"  spec.{f} = {getattr(spec, f)}")
    tm = m.get("torqueMarginFraction")
    if tm is not None:
        print(f"  torqueMarginFraction min={np.min(tm):.6g} max={np.max(tm):.6g}")
    tm2 = m.get("torqueMargin")
    if tm2 is not None:
        print(f"  torqueMargin min={np.min(tm2):.6g} max={np.max(tm2):.6g}")
    for extra in ("minMargin", "J", "JBest", "objective", "fval", "exitflag", "iterations"):
        if extra in m and not hasattr(m[extra], "_fieldnames"):
            v = m[extra]
            v = v.tolist() if isinstance(v, np.ndarray) and v.size <= 6 else v
            print(f"  {extra} = {v}")
