import scipy.io as sio
import numpy as np
import csv, os

R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def show(name, v, depth=0, maxdepth=2):
    pad = "  " * depth
    if isinstance(v, np.ndarray):
        if v.dtype.names:  # struct array
            print(f"{pad}{name}: struct array shape={v.shape} fields={v.dtype.names}")
            if depth < maxdepth:
                flat = v.ravel()
                for i, item in enumerate(flat[:3]):
                    for f in v.dtype.names:
                        show(f"{name}[{i}].{f}", item[f], depth+1, maxdepth)
        elif v.size <= 12:
            print(f"{pad}{name}: {np.array2string(np.asarray(v).ravel(), precision=6)} (shape {v.shape})")
        else:
            a = np.asarray(v).ravel()
            nanmask = np.isnan(a) if a.dtype.kind == 'f' else np.zeros(a.shape, bool)
            print(f"{pad}{name}: shape {v.shape}, n={a.size}, allNaN={bool(nanmask.all()) if a.dtype.kind=='f' else 'n/a'}, "
                  f"min={np.nanmin(a) if a.dtype.kind in 'fiu' else '?'}, max={np.nanmax(a) if a.dtype.kind in 'fiu' else '?'}")
    else:
        print(f"{pad}{name}: {type(v).__name__} {str(v)[:200]}")

def dump(fname, keys=None):
    p = os.path.join(R, fname)
    print("=" * 100)
    print("FILE:", fname)
    d = sio.loadmat(p, squeeze_me=False, struct_as_record=True)
    for k, v in d.items():
        if k.startswith("__"):
            continue
        show(k, v, 0)

for f in ["Vas_Pam_20mm_Result.mat",
          "Vas_Pam_20mm_Result_20260920_1519.mat",
          "Vas_Pam_20mm_Result_20260925.mat",
          "Bifemsh_20mm_Result.mat",
          "Bifemsh_20mm_Result_pulley_20260926_1652.mat",
          "BiPulley_opensim_bifemsh_r_20260926_1702.mat"]:
    dump(f)

print("=" * 100)
print("GAIT CSV library means (excluding the kine_score tuning row):")
rows = list(csv.DictReader(open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\gait_validation_20260930.csv")))
from collections import defaultdict
by = defaultdict(list)
mostneg = {}
for r in rows:
    v, ref, s = r["variant"], r["reference"], float(r["kine_score"])
    if ref == "kine_score":
        print(f"  {v}: tuning-ref = {s:.4f}")
        continue
    by[v].append((ref, s))
for v, lst in by.items():
    scores = [s for _, s in lst]
    ref, worst = min(lst, key=lambda t: t[1])
    print(f"  {v}: n={len(lst)} library mean = {np.mean(scores):.1f}, "
          f"rounded {np.mean(scores):.0f}, most-negative {ref} {worst:.1f} (rounded {worst:.0f}), "
          f"ratio mean/tuning = {abs(np.mean(scores))/186.4:.1f}x-scale-check")
    # count running vs walking vs other
    run = [x for x in lst if "_Run_" in x[0]]
    ong = [x for x in lst if x[0].startswith("ong_")]
    other = [x for x in lst if "_Run_" not in x[0] and not x[0].startswith("ong_")]
    print(f"    counts: run={len(run)} ong={len(ong)} other={[o[0] for o in other]}")
