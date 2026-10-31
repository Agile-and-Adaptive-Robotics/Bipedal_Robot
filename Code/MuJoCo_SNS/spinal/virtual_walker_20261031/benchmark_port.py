"""Matched interpreted/compiled end-to-end timing, excluding loading/JIT/save."""
import json
import numpy as np
from run_virtual_walker import run
from import_virtual_walker import HERE
kwargs=dict(mode='air',duration=4.,neural_time_scale=2.)
a=run(**kwargs,reference=True,tag='benchmark_final_reference');b=run(**kwargs,tag='benchmark_final_numba');base=HERE/'runs/W2L_VirtualWalker';za=np.load(base/'benchmark_final_reference/trajectory.npz');zb=np.load(base/'benchmark_final_numba/trajectory.npz');diff={k:float(np.max(abs(za[k]-zb[k]))) for k in ['q','angles','neuro','forces','root','rail_vertical','foot_vertical','foot_clearance']}
out=dict(machine='easteregg2',simulation_duration_s=4.,neural_time_scale=2.,reference_warmed_wall_s=a['warmed_wall_s'],numba_warmed_wall_s=b['warmed_wall_s'],acceleration=a['warmed_wall_s']/b['warmed_wall_s'],numba_sim_seconds_per_wall_second=4/b['warmed_wall_s'],max_abs_differences=diff,matched_equations=True,within1e_minus8=all(x<1e-8 for x in diff.values()),scope='Source80-neuron/150-edge circuit plus24-unit Hill/receptor bank, MuJoCo2.3.7 implicitfast2ms, native.2ms neural substeps, diagnostics/post-step forward included; excludes source/model loading, compilation/warm-up, NPZ compression and rendering. Reference calls Python scalar versions of the same equations, not SNS-Toolbox dense backend or AnimatLab. Both cases use the same transformed anatomical frame and physics.')
(HERE/'runs/BENCHMARK_FINAL.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
