---
description: Four-place wiring audit for a new params.G knob - params default, runner JSON-RULE loader, curriculum KEYS/suggest/seed - read-only
---

# /knob-check

The recurring trap: a new gain knob wired in one place but missing in another silently
misbehaves (the v7 renshaw omission cost a silent 0.14 kine delta; a missing seed key
gets SAMPLED by optuna - the stage-2 lesson). Audit all four places, read-only:

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" knob_check.py <knob> [<knob> ...]
```

Checked:
1. `params.py` - default entry in the G dict
2. `runner.py` - consumed by the `--best` JSON loader (the JSON RULE tuple or a
   dedicated `if key in best` branch)
3. `_curriculum.py` - present in the stage KEYS tuple
4. `_curriculum.py` - has a `suggest_*` call AND a seed-dict entry

Exit 1 = incomplete. **Wiring fixes are Ben's call** - this command reports, it never
edits. Run it for every knob a new study round introduces BEFORE launching the study,
and again before recording a winner json as reproducible.

Related standing rule (JSON RULE): every knob a study uses must also be SAVED in the
winner json params and have a loader branch - a `--bestX` reproduction that silently
drops a knob is the classic false-regression.
