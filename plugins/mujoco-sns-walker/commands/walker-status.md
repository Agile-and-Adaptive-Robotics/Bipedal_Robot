---
description: One-stop walker status - optuna studies, top trials, curriculum log milestones, stray processes, newest npz files
---

# /walker-status

Print the full state of the walker tuning campaign in one read-only call. Use this
FIRST in any walker session, and after any long run finishes, instead of hand-writing
optuna/sqlite queries and log greps.

## Run

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" status.py
```

(If `CLAUDE_PLUGIN_ROOT` is not set in the shell, use the repo source path:
`<repo>\plugins\mujoco-sns-walker\scripts\walker.cmd status.py`.)

With a study name for a top-trials table with dynamic parameter columns:

```bat
...\walker.cmd status.py --study curr_s3k_nocross
```

## Output sections

1. `== optuna studies ==` - every study in `optuna_walk.db` with trial count + best value/trial.
2. top trials of `--study` (value-sorted, parameter columns chosen from what the trials
   actually carry - never stale hardcoded keys).
3. `== curriculum logs (newest) ==` - newest 3 `curriculum_*.log`/`curr_s*.log` with their
   last milestone lines (`== stage N best`, `DONE`, `FAILED`, `seeded`).
4. `== walker processes ==` - live `runner.py`/`_curriculum.py`/`run_*.bat` processes
   (empty = nothing running; check this after TaskStop to catch orphans).
5. `== newest spinal_run*.npz ==` - newest 5 run artifacts with sizes/ages.

## Failure modes

- `optuna read failed: ...` - db missing/locked (a running study holds it briefly);
  the other sections still print.
- `could not locate ...spinal` - run from the repo, or set `AARL_SPINAL`.
- Wrong machine paths - `AARL_PYTHON` env var overrides the python; `walker.cmd`
  tries EB475WS4 -> easteregg2 -> laptop candidates automatically.
