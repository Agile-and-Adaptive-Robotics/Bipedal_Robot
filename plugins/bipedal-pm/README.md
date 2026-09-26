# bipedal-pm

Project-management plugin for the Bipedal_Robot repo (AARL, PSU). Four jobs, all in
service of priority 1: the dissertation (Sept 2026).

| command / skill | job |
|---|---|
| `/diss` (dissertation-focus) | dissertation status, placeholder audits, figure integrity, ranked task list |
| `/digest` (project-digest) | compile readmes/reports/git/figures into a dated digest + dissertation-disposition verdicts |
| `/handoff` (chatgpt-handoff) | write CHATGPT_HANDOFF.md, read/verify CHATGPT_REPORT.md |
| `/hygiene` (repo-hygiene) | read-only file-structure audit; reference-safe move proposals (nothing moves without Ben's explicit approval) |

SessionStart hook prints a short banner: dissertation priority/deadline, latest
dissertation edit, CHATGPT file ages, known root strays.

## Scripts
- `hooks/scripts/session_banner.ps1` - read-only, fail-silent, always exit 0.
- `scripts/hygiene_audit.ps1` - `audit` mode (unexpected root entries vs the AGENTS.md
  canonical map, suspicious names, zero-byte files, >50 MB files) and `refs -Name X`
  mode (reference-impact grep before any move/rename). Never modifies anything.

## Notes
- Repo-root resolution: `BIPEDAL_REPO` env -> cwd -> the three known machine paths
  (mirrors mujoco-sns-walker's wcommon.py). Installed copies live in the host plugin
  cache and must never assume the repo path.
- Optional `Documentation\Reports and Papers\Dissertation\Notes\DEFENSE_DATE.txt`
  (single date line) makes the banner show a countdown.
- Vision QA of dissertation figures reuses the project's existing zai-mcp-server tools
  (already provided by the mujoco-sns-walker plugin) and/or visual-judge subagents; this
  plugin does not duplicate an MCP server for it.
