# One-off verification for prose-fix round 1 (2026-10-03). Reads recorded artifacts only.
import csv, json
import numpy as np
from scipy.io import loadmat

base = r"D:\Github\Bipedal_Robot"
spinal = base + r"\Code\MuJoCo_SNS\spinal"
res = base + r"\Code\Matlab\Mesh_Optimization\Results"

# 1) gait validation CSV: class split + means for s3k
rows = []
with open(spinal + r"\gait_validation_20260930.csv") as f:
    for r in csv.DictReader(f):
        rows.append((r["variant"], r["reference"], float(r["kine_score"])))

def cls(ref):
    if ref.startswith("subject") and "_Run_" in ref:
        return "running"
    if ref.startswith("ong_"):
        return "ong_walking"
    return "other:" + ref

for variant in ("s3k", "s3kpruned", "syn6", "w2lvar"):
    sc = [s for v, ref, s in rows if v == variant]
    run = [s for v, ref, s in rows if v == variant and cls(ref) == "running"]
    ong = [s for v, ref, s in rows if v == variant and cls(ref) == "ong_walking"]
    oth = [(ref, s) for v, ref, s in rows if v == variant and cls(ref).startswith("other")]
    print(variant, "n_total", len(sc), "mean_all %.1f" % np.mean(sc),
          "| n_run", len(run), "mean_run %.1f" % np.mean(run),
          "| n_ong", len(ong), "mean_ong %.1f" % np.mean(ong),
          "| other", oth)

# 2) Case_40_motion provenance: inspect npz keys + any source metadata
d = np.load(spinal + r"\gait_refs\Case_40_motion.npz", allow_pickle=True)
print("Case_40 npz keys:", sorted(d.files))
for k in d.files:
    a = d[k]
    if a.dtype.kind in "US" or a.size <= 8:
        try:
            print("  ", k, "=", a)
        except Exception as e:
            print("  ", k, "unreadable:", e)
    else:
        print("  ", k, a.shape, a.dtype)
# compare against one ong + one run ref for T/duty context
for name in ("ong_speed_100", "subject01_Run_20002"):
    dd = np.load(spinal + r"\gait_refs\%s.npz" % name, allow_pickle=True)
    keep = {k: dd[k] for k in dd.files if k in ("T", "duty", "source", "fs", "speed")}
    print(name, {k: (v.tolist() if v.size < 5 else v.shape) for k, v in keep.items()})

# 3) BiPulley bifemsh mat: st.tauAbs NaN census
m = loadmat(res + r"\BiPulley_opensim_bifemsh_r_20260926_1702.mat", squeeze_me=False, struct_as_record=False)
print("BiPulley bifemsh top keys:", [k for k in m if not k.startswith("__")])
st = m["st"]
print("st fields:", st._fieldnames)
tau = np.asarray(st.tauAbs)
print("tauAbs shape", tau.shape, "total", tau.size, "NaN", int(np.isnan(tau).sum()),
      "finite", int(np.isfinite(tau).sum()), "max finite %.4f" % np.nanmax(tau))

# 4) pulley + direct flexor mats: torque field names + peaks
for f in ("Bifemsh_20mm_Result_pulley_20260926_1652.mat", "Bifemsh_20mm_Result.mat"):
    mm = loadmat(res + "\\" + f, squeeze_me=True, struct_as_record=False)
    keys = [k for k in mm if not k.startswith("__")]
    print(f, "top keys:", keys)
    for k in keys:
        v = mm[k]
        if hasattr(v, "_fieldnames"):
            tf = [fn for fn in v._fieldnames if "orque" in fn or "tau" in fn.lower()]
            if tf:
                print("   ", k, "torque-ish fields:", tf)
                for fn in tf:
                    arr = np.asarray(getattr(v, fn))
                    if arr.dtype.kind in "fiu":
                        print("      ", fn, arr.shape, "min %.4f max %.4f" % (np.nanmin(arr), np.nanmax(arr)))
