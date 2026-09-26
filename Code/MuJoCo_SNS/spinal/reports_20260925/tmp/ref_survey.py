"""Read-only survey of gait_refs/*.npz per-family scalars (goal3 plan input)."""
import sys, io, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[2]
REFS = HERE / "gait_refs"
rows = []
for p in sorted(REFS.glob("*.npz")):
    z = np.load(p, allow_pickle=True)
    rows.append({
        "ref": p.stem,
        "keys": sorted(z.files)[:0] or None,  # keys dumped once below
        "duty_r": round(float(z["duty_r"]), 3),
        "duty_l": round(float(z["duty_l"]), 3),
        "T_r": round(float(z["T_r"]), 3),
        "knee_min_r": round(float(z["knee_min_r"]), 1),
        "hip_range_r": round(float(z["hip_range_r"]), 1),
        "ankle_range_r": round(float(z["ankle_range_r"]), 1),
        "lag_rl": round(float(z["lag_rl"]), 3),
        "ds": round(float(z["ds"]), 3),
        "grf_style": str(z["grf_style"]),
    })
one = np.load(REFS / "ong_speed_100.npz", allow_pickle=True)
out = {"npz_keys": sorted(one.files), "refs": rows}
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
print(json.dumps(out, indent=1))
