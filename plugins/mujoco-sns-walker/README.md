# mujoco-sns-walker (ZCode plugin)

Workflow tooling for the AARL MuJoCo + SNS-Toolbox spinal walker
(`Code\MuJoCo_SNS\spinal`). Built 2026-09-23 from the recurring pain points of
the walker debugging/tuning sessions. Source lives here (repo-tracked, text
only); ZCode installs a copy from the local dev marketplace.

## Commands

| Command | Purpose |
|---|---|
| `/walker-status` | optuna studies + top trials + log milestones + processes + npz |
| `/walker-run` | compile gate → env-resolved detached launch, dated log, lock-free npz |
| `/walker-probe` | generalized s3k-style winner probe (full per-leg gait metrics) |
| `/study-manage` | show / seed / archive / purge (purge = archive-first + `--yes`) |
| `/connectome-check` | validate connectome_gains.json vs runner's real contract |
| `/connectome-gui` | serve + open both connectome editors over localhost |
| `/circuit-vision` | paper-figure ↔ our-circuit vision diff, figure audits, Ben markups |
| `/supervisor-gate` | spawn SUPERVISOR_AGENT.md at milestones M1/M2/M3 |
| `/knob-check` | four-place wiring audit for a new params.G knob (read-only) |

## MCP servers (`.mcp.json`)

Declares `zai-mcp-server` (vision), `web-search-prime`, `web-reader` under the
same names as the user-scope config — user entries override (no duplicates);
on machines without them the plugin provides the servers itself. Auth travels
via the `${ZAI_API_KEY}` template → set the `ZAI_API_KEY` user environment
variable (done on EB475WS4; `setx ZAI_API_KEY <key>` elsewhere). **Never
commit a real key into this folder.** `zai-mcp-server` needs Node ≥ 22 on PATH
(EB475WS4 + easteregg2: `D:\NodeJS\node-v22.23.2-win-x64`, user PATH —
easteregg2's installed 2026-09-23; laptop: not yet).

## Hooks

- `UserPromptSubmit` — off-peak/keep-running phrases inject a standing
  autonomy directive (keep iterating, don't wait for go-ahead).
- `PreCompact` — Windows toast + beep when session context compacts
  (the "start a fresh chat" nudge).

PowerShell hook scripts must stay PURE ASCII (PS 5.1 misparses BOM-less
UTF-8 CJK literals — build Chinese phrases from `[char]` codepoints).

## Scripts

All invoked via `scripts\walker.cmd <script.py> [args]`, which resolves the
SNS env python per machine (EB475WS4 → easteregg2 → laptop; `AARL_PYTHON` /
`AARL_SPINAL` override) and sets `CONDA_PREFIX`. Study commands
(status/study-manage/probe) additionally need `optuna` in that env (5.0.0 on
EB475WS4 + easteregg2; the shared `optuna_walk.db` is schema 12 — keep the
optuna major version matched across machines).

## Maintenance

After source edits: bump `.zcode-plugin/plugin.json` version + the
`marketplace.json` entry, then in ZCode: marketplace sources → refresh
`dev-bipedal-robot-0c995be2` → Personal → the plugin → Update.
