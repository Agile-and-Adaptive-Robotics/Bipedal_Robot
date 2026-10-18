"""Check undated (2022-era) result mats for the initial single-bracket identification."""
import scipy.io as sio
import numpy as np

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

for f in [r"\minimizeFlxPin10_results.mat", r"\minimizeExt10mmX3_results.mat"]:
    print(f"\n===== {f} =====")
    try:
        m = sio.loadmat(BASE + f, squeeze_me=True, struct_as_record=False)
    except Exception as e:
        print("  LOAD FAILED:", e)
        continue
    for k in ["results_sort_actual", "results_sort", "filtered_results", "sol_actual", "pick"]:
        if k in m:
            v = m[k]
            if np.isscalar(v):
                print(f"  {k} = {v!r}")
            elif isinstance(v, np.ndarray) and v.dtype.kind != "O":
                print(f"  {k}: shape {v.shape}")
                if v.ndim == 2:
                    for r in range(min(3, v.shape[0])):
                        print(f"    row{r}: {np.array2string(v[r], precision=6)}")
    # pick the PICK-selected solution if present
    for k in m:
        if k.startswith("__"):
            continue
        v = m[k]
        if isinstance(v, str) and len(v) < 120:
            print(f"  {k} = {v!r}")
