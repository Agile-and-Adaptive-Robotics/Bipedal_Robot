"""Recalculate the paired muscle-force statistics directly from the archived traces."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

here = Path(__file__).resolve().parent
diss = here.parents[1]
root = diss.parents[2]
spinal = root / "Code/MuJoCo_SNS/spinal"
frames = diss / "Figures/95-AppendixD/Components/force_comparison/force_arrays.npz"
summary = spinal / "muscle_force_compare.csv"
reference = root / "Solid_Models/OpenSim/Gait2392_Robotbody/ResultsBSolve/zz_bsolve_StaticOptimization_force.sto"
source = reference.read_text().splitlines()
header = next(i + 1 for i, line in enumerate(source) if line.strip().lower() == "endheader")
names = source[header].split()
data = np.loadtxt(source[header + 1:])

with np.load(frames, allow_pickle=False) as saved:
    times, muscles, mj, os = (saved[k] for k in ("t", "muscles", "mj", "osim"))
assert mj.shape == os.shape == (121, 12)
np.testing.assert_allclose(times[[0, -1]], [0.5, 2.5])
with summary.open(newline="") as f:
    rows = list(csv.reader(f))[1:]
assert len(rows) == len(muscles)
results = {}
for j, row in enumerate(rows):
    muscle = str(muscles[j])
    assert muscle == row[0]
    model, ref = mj[:, j], os[:, j]
    np.testing.assert_allclose(np.interp(times, data[:, 0], data[:, names.index(muscle)]), ref)
    baseline_error = np.sum((ref - ref.mean()) ** 2)
    model_error = np.sum((model - ref) ** 2)
    r2 = 1 - model_error / baseline_error
    corr = np.corrcoef(ref, model)[0, 1]
    ratio = model.mean() / ref.mean()
    np.testing.assert_allclose([r2, corr**2, ratio, ref.max(), model.max()],
                               [float(row[i]) for i in (1, 2, 4, 5, 6)], rtol=1e-9, atol=1e-12)
    results[muscle] = dict(reference_mean_N=ref.mean(), converted_mean_N=model.mean(),
                           mean_ratio=ratio, force_R2=r2, signed_Pearson_r=corr, shape_r2=corr**2,
                           squared_error_relative_to_reference_mean=model_error / baseline_error)
agree = [m for m, v in results.items() if .88 <= round(v["force_R2"], 2) <= .98 and
         abs(v["mean_ratio"] - 1) <= .1]
assert len(agree) == 8
output = dict(reference="OpenSim static-optimization estimate at each pose", poses=len(times),
              intervals_s=[float(times[0]), float(times[-1])], muscles=results, close_agreement=agree,
              source_sha256={p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in [frames, reference, summary, spinal / "_muscle_force_compare.py"]})
idx, = np.flatnonzero(np.isclose(times, 1.550))
j, = np.flatnonzero(muscles == "rect_fem_r")
output["rectus_example"] = dict(time_s=float(times[idx]), OpenSim_N=float(os[idx,j]), MuJoCo_N=float(mj[idx,j]))
np.testing.assert_allclose([os[idx,j], mj[idx,j]], [96.38, 415.57], atol=.005)
(here / "FORCE_CHECK.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
for m in ("rect_fem_r", "bifemsh_r"):
    print(m, results[m])
