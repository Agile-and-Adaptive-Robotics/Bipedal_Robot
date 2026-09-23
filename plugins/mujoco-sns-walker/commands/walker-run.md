---
description: Launch walker runs safely - compile gate first, then env-resolved background launch with dated log and lock-free npz
---

# /walker-run

Launch `runner.py` or `_curriculum.py` WITHOUT the fragile cd/set/redirect cmd chains.
Always gate first, then launch detached, then poll - never babysit turn-by-turn.

## Step 1 - compile gate (always, before any launch)

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" run_gate.py compile
```

PASS required. A FAIL means a syntax error in runner/_curriculum/params/build_network -
fix before launching anything.

## Step 2 - launch

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" run_gate.py go <tag> runner <runner flags...>
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" run_gate.py go <tag> curriculum <stage> [n_trials]
```

Examples:

```bat
... run_gate.py go s3l runner --eval --drive 2.077
... run_gate.py go s3l curriculum 3 40
... run_gate.py go repro runner --fitted --best10 --eval
```

What `go` does: resolves the SNS python + sets CONDA_PREFIX in-process (no `set` chains),
writes `<tag>_<YYYYMMDD>.log` in the spinal folder, sets `AARL_NPZ=spinal_run_<tag>.npz`
(unique name -> no Spyder/AV lock collisions on `spinal_run.npz`), launches detached,
prints PID + log path.

Use `print` instead of `go` for a dry run of the resolution.

## Step 3 - poll

- `/walker-status` (studies table + log milestones + process list), or
- tail the printed log; `findstr /c:"== stage" /c:"DONE" /c:"FAILED" <log>` for milestones.

Curriculum runs take hours (40 trials x ~1-3 min eval). Wait with sleep/ping between
checks; do not spin. After a study finishes, `/walker-probe` the winner.

## Rules

- Study names live in `_curriculum.py` main() (`{1: curr_s1..., 2: curr_s2b..., 3: curr_s3k_nocross}`);
  a NEW round requires Ben's naming/edit first - the plugin never edits wiring.
- Seeding a new stage-3 round from a previous winner: `study_manage.py seed <study> <trial>`
  BEFORE `run_gate.py go ...`.
- NEVER purge/delete studies without `/study-manage` archive + explicit confirmation.
