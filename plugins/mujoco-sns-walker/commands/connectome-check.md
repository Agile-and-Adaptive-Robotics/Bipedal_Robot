---
description: Validate a connectome spec JSON against the runner's real consumption contract before a run silently applies 0 rules
---

# /connectome-check

Validate `connectome_gains.json` (or any spec file) against what runner.py actually
reads (the `spec.get("rules", {})` consumer). Catches the silent failure classes:

- block-editor export (nodes/synapses top-level) saved as connectome_gains.json ->
  parses fine, applies 0 rules                    [ERROR, the known trap]
- `gain_key` not in `params.G` -> rule silently skipped            [ERROR]
- two enabled rules sharing one gain_key (ia_recip & ii_inh both ->
  `ia_to_antagonist`) -> dict order decides                       [WARN]
- disabled rules are SKIPPED not zeroed (runner `continue`s - the
  CONNECTOME.md "disabled -> 0.0" note is aspirational)           [INFO]
- any enabled rule with hops>=1 sets the single global `full_rules`
  (no per-rule direct/IN mixing exists)                           [INFO]
- `gain` sign is discarded (runner takes abs)                     [INFO]

## Run

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" connectome_check.py [--file X.json]
```

Default file: `<spinal>/connectome_gains.json`. Exit 1 = fix before running; 0 = OK.
"FILE ABSENT" is a valid state (no spec applied - params defaults run as-is).

ALWAYS run this after Ben exports a new spec from the editor and BEFORE launching a
run that should honor it - the only runtime trace of a bad spec is a single startup
line `connectome_gains.json applied (0 rules, full_rules=0.0)`.
