"""Targeted: locate front row with X1~39978.7/X2~14733.75 in the flx noT3 pooled front;
read XiUsed from the redesign mats."""
import scipy.io as sio
import numpy as np
from pathlib import Path

TD = Path(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo")
MO = Path(r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results")

m = sio.loadmat(TD / "minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat", squeeze_me=True, struct_as_record=True)
fc = np.atleast_2d(m["filtered_results"])
ac = np.atleast_2d(m["all_candidates"])
print("filtered_results", fc.shape, "cols; all_candidates", ac.shape)
# find row with X1 ~ 39978.7 (linear) in filtered_results
for arr, name, lin_cols in [(fc, "filtered_results", (3, 4, 5)), (ac, "all_candidates", None)]:
    hits = []
    for i in range(arr.shape[0]):
        row = arr[i]
        for j in range(arr.shape[1]):
            v = row[j]
            if abs(v - 39978.7285) < 0.1:
                hits.append((i + 1, j, row))
    for h in hits[:3]:
        print(f"{name} 1-based row {h[0]} col {h[1]}:", np.array2string(h[2], precision=6, max_line_width=200))

for f in ["Bifemsh_20mm_Result.mat", "Vas_Pam_20mm_Result_20260920_1519.mat", "Vas_Pam_20mm_Result.mat"]:
    p = MO / f
    if not p.exists():
        print(f"\n== {f}: MISSING")
        continue
    r = sio.loadmat(p, squeeze_me=True, struct_as_record=True)
    print(f"\n== {f} ==")
    for k in r:
        if k.startswith("__"):
            continue
        v = r[k]
        if hasattr(v, "dtype") and v.dtype.names and "XiUsed" in v.dtype.names:
            xu = v["XiUsed"]
            print(f"   {k}.XiUsed = {xu}")
        elif k in ("LOCKSRC", "flxRow", "xi1lock", "xi2lock"):
            print(f"   {k} = {v}")
    # some store XiUsed at top level inside a nested struct; dump shallow fields of each top struct
    for k in r:
        if k.startswith("__"):
            continue
        v = r[k]
        if hasattr(v, "dtype") and v.dtype.names:
            try:
                names = [n for n in v.dtype.names]
                print(f"   [{k}] fields: {names}")
                for n in names:
                    x = v[n]
                    if np.ndim(x) == 0 or (hasattr(x, "size") and x.size <= 8):
                        print(f"      {n} = {x}")
            except Exception as e:
                print("   field dump err", e)
