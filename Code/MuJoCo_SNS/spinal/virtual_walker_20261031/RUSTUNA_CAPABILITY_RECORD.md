# Rustuna capability record

Session 1, easteregg2. Separate optimizer check; no walker optimization or mechanical perturbation study has been run.

## Isolated runtime

Rustuna requires Python 3.11 or newer according to its [getting-started guide](https://rustuna.readthedocs.io/en/latest/tutorial/getting-started/). The current MuJoCo/SNS environment is Python 3.10.21, so it cannot directly import the supported wheel. Existing `D:\Anaconda\python.exe` is Python 3.11.5. Installed the Rustuna 0.2.1 Windows prebuilt wheel only into ignored `D:\GitHub\Bipedal_Robot\tmp\rustuna_runtime_py311`, with `--only-binary=:all: --no-deps`. No conda transaction, NumPy replacement, MuJoCo environment change, or Optuna database migration was performed.

Reproduce installation:

```powershell
& 'D:\Anaconda\python.exe' -m pip install --only-binary=:all: --no-deps --target 'D:\GitHub\Bipedal_Robot\tmp\rustuna_runtime_py311' 'rustuna==0.2.1'
```

Rustuna can therefore supply suggestions from the Python 3.11 process while an isolated persistent Python 3.10 worker performs SNS/MuJoCo evaluations. This avoids importing the old scientific stack into the newer interpreter or paying a fresh simulator startup/JIT compilation on every trial. A model-specific worker and its baseline gate still require implementation and validation; this runtime check alone is not that transfer pipeline.

## Tested API and bounded benchmark

`rustuna_capability.py` performs 256 TPE trials on the four-dimensional known quadratic with optimum `[0.2,−0.4,0.6,−0.8]`. Seed 100031, 10 startup trials, and independent TPE were fixed explicitly. Both backends reproduced their complete parameter/value history exactly in a second fresh study. Histories are saved in `runs\rustuna_capability` with hashes.

| Backend | Interpreter | First study wall time | Per trial | Best objective |
| --- | --- | ---: | ---: | ---: |
| Rustuna 0.2.1 | Python 3.11.5 | 0.0559 s | 0.218 ms | 0.01043 |
| Optuna 5.0.0 | Python 3.10.21 | 1.8263 s | 7.134 ms | 0.001371 |

This is approximately 32.7× less optimizer bookkeeping time for this cheap objective, not a MuJoCo speed-up. The Python versions differ, and the suggestion histories and solution quality differ between backends. The repeated studies took 0.0422 s and 1.9615 s respectively. In a multi-second simulation, the absolute saving of several milliseconds per trial is modest; physics/neural acceleration remains the dominant opportunity.

The actual installed API differs from some latest documentation: `Study.best_params` and `Study.best_value` are absent in 0.2.1. `Study.best_trial.params` works. Seeded `TPESampler`, `suggest_float`, `create_study`, and `optimize` were verified. CMA was not tested; the official [sampler reference](https://rustuna.readthedocs.io/en/latest/api/samplers/) identifies its Python `cmaes` backend, so a Rustuna name alone does not imply a fully Rust CMA implementation. Feature availability should follow the installed version, not untested documentation examples.

Run the checks:

```powershell
& 'D:\Anaconda\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\rustuna_capability.py' --backend rustuna
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\rustuna_capability.py' --backend optuna
```

## Earlier transfer scope, now paused

The user clarified that the desired change is rat-controller transfer to Gait2392 or the CAD/BPA model, innervating only the corresponding named muscles. The proposed synthetic +10% limb-mass experiment was not executed. An unused opt-in `physics_variants.py` helper was prepared before this clarification, but it is not the user's requested transfer study and should not be wired into that study by default.

The named-muscle audit is `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\RAT_TO_GAIT_MUSCLE_MAP.md`. The public rat precursor is 120 neurons and currently has a right RG latch. The separate actual108 W2L baseline is not a rat baseline and must not be used to certify rat-to-human transfer. An accepted rat baseline, explicit group/head mappings, a source-specific persistent worker, and corresponding sensor/actuator calibration are prerequisites before a few-hour bounded transfer optimization. Existing Optuna SQLite studies remain untouched.

## Current authorized work: newer108 W2L contact feedback

The user subsequently paused rat/Gait transfer and requested the current newer108 W2L ground model, specifically testing excessive contact feedback. The optimizer driver is prepared for that target. The physics agent owns source parsing, pathway classification, ablations, physical objective, and final acceptance. No model-family copies, synthetic mass-change study, or existing Optuna database changes are part of this work.

The current first-stage budget is 64 trials or 900 seconds. Candidate variables are shared bilateral contact-pathway class gains, limited to classes actually present/consumed in the source circuit, with at most an additional shared MN recruitment gain if needed. An evaluated, verified baseline/ablation protocol must precede tuning; it need not already be a walking pass. Dragging, nonfinite behavior, knee hyperextension, and apparatus-only support without foot loading must not receive a physical stepping pass. Fresh-process winner confirmation remains a separate requirement after optimization.

## Persistent worker and checkpoint contract

- `rustuna_runtime.py`: loads the isolated 0.2.1 wheel in Python 3.11.
- `rustuna_transfer.py`: bounded ask/tell driver, complete parameter/diagnostic JSONL recording, exact seeded resume, provenance hashes, and persistent original deadline.
- `rustuna_myo_worker.py`: imports one model-specific adapter in the existing Python 3.10 interpreter, keeps its scientific imports/JIT state alive, and exchanges JSON lines.
- `rustuna_resume_selfcheck.py`: a 64-trial toy study matched exactly after a JSON checkpoint at trial 32; a corrupted checkpoint was rejected. An explicit enqueued candidate, including a zero-valued parameter, was also verified under the actual installed API.
- `rustuna_worker_selfcheck.py`: 100 known-objective requests reused the same worker process, gave exact answers, and averaged 0.128 ms round-trip overhead. This excludes physics and objective cost.

The adapter exports `handle_request(request)`. A `preflight` request receives the complete manifest and must return `ok: true` with `ready_for_optimization` and baseline/source evidence. An `evaluate` request receives `trial_number`, a unique `evaluation_id`, `evaluation_seed`, full `params`, and the manifest; it returns `ok: true`, a finite `value`, and physical diagnostics. Use the evaluation ID in artifact names so an interrupted attempt is preserved. Non-protocol output is redirected into the worker stderr log. Unknown or unused optimization keys must be rejected by the adapter.

Manifest fields: `adapter_path`, `output_dir` inside this package's ignored `runs`, `search_space` with float `low`/`high` and optional `log`, `fixed_parameters`, `seed`, `max_trials`, `wall_budget_s`, and `trial_timeout_s`. Optional `initial_candidates` contains complete sampled-parameter dictionaries, allowing exact baseline/off-gain points rather than hoping continuous TPE samples zero. Additional model-specific fields pass through unchanged to the adapter. Set the W2L source path and expected census explicitly to avoid stale-family confusion.

```powershell
& 'D:\Anaconda\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\rustuna_transfer.py' --manifest 'FULL_PATH_TO_REVIEWED_MANIFEST.json' --probe-only
& 'D:\Anaconda\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\rustuna_transfer.py' --manifest 'FULL_PATH_TO_REVIEWED_MANIFEST.json'
# Only for an interrupted run with identical configuration, adapter, and Rustuna version:
& 'D:\Anaconda\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\rustuna_transfer.py' --manifest 'FULL_PATH_TO_REVIEWED_MANIFEST.json' --resume
```

The driver cannot turn its lowest objective into a scientific pass. BEST.json is a candidate checkpoint; inspect its physical diagnostics and repeat it in a fresh worker. The false-preflight toy adapter was verified to prevent any optimization trials. No model optimization has been launched by this capability subtask.
