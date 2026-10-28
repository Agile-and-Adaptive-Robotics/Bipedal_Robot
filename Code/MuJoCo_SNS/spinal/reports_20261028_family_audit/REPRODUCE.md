# Reproduce on easteregg2

Machine: **easteregg2**. Repository: **`D:\GitHub\Bipedal_Robot`**.
Interpreter: **`D:\Anaconda\envs\myo\python.exe`**.
Do not substitute the other machines' interpreter paths.

Run in PowerShell:

```powershell
Set-Location -LiteralPath 'D:\GitHub\Bipedal_Robot'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\audit_families_20261028.py' w2lvar robot54
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\verify_preallocated_families_20261028.py' w2lvar robot54
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\family_human_trials_20261028.py' w2lvar baseline --preallocated
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\family_authority_trials_20261028.py' w2lvar --authority 0.1 --ramp --backend preallocated
```

For the inventory and neural verifier, replace `w2lvar` with `s3k` or `syn6`,
and `robot54` with `human92`, to cover all six combinations. These commands
overwrite their own named exploratory outputs, not production winners.

The optional fast evaluator is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\preallocated_exact_eval.py`.
It requires NumPy and Numba from the existing environment. It preserves fixed-tau
NaP and the NumPy synaptic reduction order. It rejects other gated backend
classes, including paper-specific voltage-dependent kinetics. For s3k, the
existing `exact_fast_eval.py` can be faster; choose from the measured results.

The human trial script supports `--reference`, no backend flag (existing
exact-fast), or `--preallocated` (new evaluator). Baseline files use distinct
suffixes, and all saved numeric arrays are compared in `RESULTS.json`.
The production `runner.py`, family builders, and parameter module are imported
without edits. The driver changes gains only inside its dedicated process.

`run_family_campaign_20261028.py` resumes the bounded human screen by skipping
existing trial JSONs. It uses three worker processes and its own logs. It does
not use or modify the existing Optuna databases. To rerun an individual case,
call the single-trial script explicitly. Gains added after the original
curricula are loaded explicitly by this experimental script and recorded in
each trial's `compiled.effective_gains`; do not assume an older `--best` loader
will apply them.

CAD `--delete-pf` removes only PF-to-MN conductances after 15 s. `--authority 0`
is the masked-current negative control, and `--authority 0.1 --ramp` exposes
10% of native synaptic current over 5–15 s. Full immediate release fails.
None of these exploratory conditions is an accepted walking controller.

The `human92` runner retains its legacy converted-model frame for baseline
replay. The CAD robot remains x anterior, y superior, z right, with y-down
gravity. These experiments do not silently rotate either model or retune
their contacts, constraints, support rig, timestep, or actuator force model.
