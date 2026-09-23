---
description: Serve the connectome editors over localhost (fixes the block editor's file:// template load) and open both in the browser
---

# /connectome-gui

Open Ben's connectome editors reliably. The block editor loads its w2l template via a
synchronous XHR to a relative path - browser-dependent on `file://`, a plain same-origin
GET over `http://localhost` - so serve the spinal folder instead of double-clicking the
html.

## Run (background)

```bat
"%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" serve_editor.py
```

Run it as a background task. It prints and opens:

- rules editor: `http://127.0.0.1:8765/connectome_editor.html`
  (exports `connectome_gains.json` - the ONLY spec runner.py consumes)
- block editor: `http://127.0.0.1:8765/connectome_block_editor.html`
  (design/documentation tool - its export is NOT consumed by the runner)

Options: `--port N`, `--no-open` (just print URLs).

NOTE for Ben: localStorage is per-origin - bookmarks made under `file://` will not
appear under `localhost` (and vice versa). Export/import JSON to move tabs.

## After Ben edits and exports

1. He saves `connectome_gains.json` into the spinal folder.
2. Run `/connectome-check` on it (validators catch the 0-rules silent drop).
3. Only then launch runs that should honor it (`/walker-run`).

Stop the server by killing the background task when done.
