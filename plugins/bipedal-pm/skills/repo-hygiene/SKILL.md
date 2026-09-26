---
name: repo-hygiene
description: Repo file-structure audit and cleaning for Bipedal_Robot - stray files, naming consistency, large-file relocation, folder moves with reference-impact analysis. Use when Ben asks to organize or clean up files/folders, asks what a stray file or folder is, wants naming conventions enforced, or wants large files handled. Read-only by default; nothing moves or deletes without Ben's explicit per-item approval.
---

# repo-hygiene

Repo: Bipedal_Robot. Machine guard: check `hostname` FIRST and read only that machine's
bullet in AGENTS.md (paths differ per machine). Repo-root resolution: `BIPEDAL_REPO` env
-> cwd -> `D:\Github\Bipedal_Robot` / `D:\GitHub\Bipedal_Robot` /
`C:\Users\Ben\Documents\GitHub\Bipedal_Robot`.

## Hard rules
1. **Audit-only by default.** Never delete, move, or rename without Ben's explicit yes
   on the SPECIFIC item. "Clean it up" is not consent for a specific irreversible act.
2. **This repo is reference-dense.** Hardcoded paths live in AGENTS.md (directory map),
   skills/plugin configs, `.vscode\`, `HowToRunCode.md`, MATLAB/Python scripts,
   `.aproj`/`.asim` XML, optuna db pathing. ANY move/rename runs the reference-impact
   protocol (Step 3) first.
3. **Large files**: propose moving OUT of the repo (`D:\temp\<name>\` like
   gait_lib_staging, or Ben's archive drive), or gitignore-if-regenerable. Never
   silently delete data.
4. **Mass naming standardization = POST-DEFENSE** (dissertation deadline is priority 1).
   Flag new violations; do not churn proven paths this week.

## Step 1 - audit
```
powershell -NoProfile -ExecutionPolicy Bypass -File "<scripts>\hygiene_audit.ps1" -Mode audit
```
Scripts location: installed plugin -> `%CLAUDE_PLUGIN_ROOT%\scripts\hygiene_audit.ps1`;
source checkout -> `<repo>\plugins\bipedal-pm\scripts\hygiene_audit.ps1`.
Read-only. Reports: unexpected repo-root entries (vs the AGENTS.md canonical map, with
origin hints), suspicious names repo-wide, zero-byte files, files > 50 MB, advisories.

## Step 2 - per-item disposition with Ben
For each finding: inspect it (size, mtime, first lines for text; recurse count for
dirs), give the best-evidence origin guess, and one recommendation:
KEEP+RELOCATE / DELETE / ARCHIVE to D:\temp / GITIGNORE / BEN'S CALL. Show evidence,
let Ben rule, record the ruling (AGENTS.md or the digest file).

## Step 3 - move protocol (only after Ben approves a specific move)
1. Impact grep:
   `powershell ... hygiene_audit.ps1 -Mode refs -Name "<name>"`
   (text extensions; skips .git, junctions, >20 MB files; flags AGENTS.md and .vscode hits).
2. Enumerate EVERY hit class and write them into the plan: AGENTS.md directory map
   (ALWAYS), `.gitignore` patterns, `.vscode` settings, skills/plugin configs,
   MATLAB/Python scripts, README/DESIGN docs, `.aproj`/`.asim` XML.
3. Apply the move (PowerShell `Move-Item -LiteralPath ...`; GitHub Desktop shows the
   rename; Ben commits).
4. Fix references in the SAME pass; update the AGENTS.md directory-map entry.
5. Smoke-verify one affected script (import/run) before declaring done.

## Known candidates (2026-09-26 baseline audit - re-verify before acting)
| item | evidence | proposed disposition |
|---|---|---|
| root file `0` | 26-node/46-edge connectome layout dump, 2026-09-24, edge TARGETS missing | KEEP, rename + relocate to `Code\MuJoCo_SNS\spinal\` - Ben rules |
| root dir `-p` | empty | delete (Ben's yes) |
| `.DS_Store` | macOS artifact | delete |
| `MUJOCO_LOG.TXT`, `opensim.log` | run logs at root | delete (opensim.log already gitignored) |
| `slprj` (root) | Simulink cache; gitignore only covers SNS_Simscape paths | gitignore + delete |
| `_git_*.bat` (4) | one-off git helpers | archive |
| `spring_series.m` | stray at root | find owner tree or archive |
| `temp`, `.zcode_tmp` | temp dirs | review, then clean |
| `CHATGPT_REPORT.md.bak_handoff_20260909` | dated backup | archive |
| `Jeffrey's practice pcb` | personal folder | Ben's call |
| `Code\MuJoCo_SNS` placement | Ben floated Neuromechanical_Models or a Python folder | NOT before defense - heavy entanglement (AGENTS.md map, HowToRunCode.md, MCP config, optuna db, plugin skills, dozens of scripts); schedule the full refs-protocol move post-defense |

## Naming conventions (advisory)
Current dominant conventions per tree: MATLAB trees keep their historical names
(PascalCase-ish with spaces in a few legacy folders), Python code is snake_case,
Solid_Models uses part codes (`NN_NN_XX_###`). Flag NEW violations; no mass renames
pre-defense. `Documentation\Project_Digests\hygiene_<YYYYMMDD>.md` holds the tracked
report when Ben wants one saved.
