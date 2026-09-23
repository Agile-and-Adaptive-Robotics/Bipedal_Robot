---
description: Probe a curriculum study winner - runs one eval and reports the full per-leg gait metric set (duty, cycles, frozen flags, kine)
---

# /walker-probe

Run the generalized s3k-style probe on an optuna study: load the best trial (or a given
one), merge `curriculum_stage3.json` + v10 BASE_MUL + trial params, run one in-process
`runner --eval`, and print the full per-leg metric set from `kine_ref.compare`.

## Run

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" probe_study.py --study <name> [options]
```

Options:
- `--trial N` - probe a specific trial instead of the best
- `--no-stage3` - skip merging curriculum_stage3.json params
- `--set key=val` - float override, repeatable (e.g. `--set renshaw=0.0 --set pm_add=0.2`)
- `--renshaw X` - BASE_MUL renshaw (default 0.5, the probe family convention)
- `--npz NAME` / `--out FILE` - output names (defaults are study-scoped)

## Reading the output

- `frozen_l`/`frozen_r` True = that leg completed no cycles (the s3g-s3k failure mode).
- `duty_*` vs human reference ~0.61; `n_cycles_*` >= 3 = real gait; `lag_rl` ~ antiphase.
- `kine_score` - higher (less negative) is better; comparable ONLY within the same
  objective version.
- Writes `<study>_probe_out.json` in the spinal folder for later diffing.

Takes ~1-3 min (one 16 s eval). If `nan=True` appears, the eval diverged - check drive
and ky_scale before trusting anything.

## Compare studies

Run it per study and tabulate `kine_score, duty_r/l, n_cycles_r/l, frozen_r/l` -
that comparison table is what Ben wants to see, never a single-leg metric alone.
