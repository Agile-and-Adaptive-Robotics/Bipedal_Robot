---
description: Optuna study lifecycle - show top trials, seed stage3 json from a winner, archive trials to dated json, purge (archive-first, confirmed)
---

# /study-manage

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" study_manage.py <action> <study> [options]
```

## Actions

- `show <study> [-n K]` - top K trials with full params (read-only).
- `seed <study> --trial N` - rewrite `curriculum_stage3.json` from trial N so the next
  `_curriculum.py 3` round chains from that winner. Prints the pf_gain/contra_swing
  shorthand like `_update_stage3.py` did.
- `archive <study...>` - dump every trial (number/value/params/state) to
  `curriculum_<study>_archive_<YYYYMMDD>.json`. NO delete. Do this before purging or
  before any risky db operation.
- `purge <study...> --yes` - archive FIRST (refuses if the study is absent), then
  `optuna.delete_study`. DESTRUCTIVE: never run without Ben's explicit go, and always
  state the archive path in the report.

## Standing rules

- Exploited/failed-study cleanups: archive under a name that says WHY
  (rename after the fact is fine), e.g. the 2026-09-20 exploit archives.
- Renaming a stage-3 study for a new round is a `_curriculum.py` source edit
  (Ben's wiring) - this tool only manages the db side.
- After purge, `/walker-status` to confirm the study list.
