---
name: project-digest
description: Compile project knowledge for Bipedal_Robot - readmes, reports, summaries, git history, figures - into a dated management digest with a dissertation-disposition verdict for every item. Use when Ben asks for a status digest, an end-of-week compilation, "what happened across the project", or what should go into the dissertation.
---

# project-digest

Goal: one dated snapshot Ben can act on, always in service of PRIORITY 1 (dissertation,
Sept 2026). Repo-root resolution: `BIPEDAL_REPO` env -> cwd -> known machine candidates.

## Inputs (sweep in this order; skip gracefully if missing)
1. **Repo state**: git log since the last digest (try `git` on PATH; else GitHub Desktop
   git at `C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-*\resources\app\git\cmd\git.exe`,
   `--oneline --since=...`); plus `Documentation\Project_Digests\` mtimes for the digest
   cadence.
2. **Standing summaries**: `AGENTS.md` (canonical memory - read the NEWEST dated
   sections), `README.md`, `PROJECT_INSTRUCTIONS.md`, `Code\MuJoCo_SNS\README.md`,
   `Code\MuJoCo_SNS\spinal\DESIGN.md`, `SADb_audit\README.md`.
3. **Reports**: `spinal\reports_20260923\EXEC_SUMMARY_20260923.md`,
   `spinal\reports_20260924\overnight_report_20260924.md`, `CHATGPT_REPORT.md`,
   Dissertation `summary.md` + `CPG_DISSERTATION_UPDATE_NOTES.md` + `Notes\figure_todo.md`.
4. **Artifacts**: figures (`Dissertation\CPG_airstepping_figs\`, `Figures\Aim1|Aim2\`,
   `ProofFinal\figs\`), result params/npz/mat, newest first by mtime.
5. **Fan-out**: for breadth, dispatch read-only Explore subagents, one per area
   (Code\Matlab incl. SNS_Simscape; Code\MuJoCo_SNS; Neuromechanical_Models;
   SADb_audit + Documentation). Each returns <=10 bullets: what changed recently,
   current status/numbers, open items. Aggregate; don't let them dump file lists.

## Output -> Documentation\Project_Digests\digest_<YYYY-MM-DD>.md (create folder on demand)
(a) Since last digest (commits + notable new/changed files).
(b) Workstream state, one short paragraph each: Xi-correction / Mesh_Optimization,
    MuJoCo-SNS spinal walker, SNS_Simscape/Simulink, AnimatLab, SolidWorks/CAD,
    SADb, Dissertation - carry the latest NUMBER or verdict from the cited report.
(c) Artifact inventory (figures/reports/params, dated).
(d) **Dissertation disposition table** - every item gets exactly one verdict:
    INCLUDE (target section/figure slot) | ALREADY IN (cite chapter/label) |
    SKIP (why) | LATER (post-defense).
(e) Blockers needing Ben (AGENTS.md "FOR BEN"/pending items, verbatim pointers).
(f) Ranked next actions by dissertation impact.

## Rules
- Never invent numbers - quote them from the cited report/script output and cite the
  file path for every claim.
- Conflicts between reports are FLAGGED, not silently resolved (e.g., stale
  `--fitted --bestN` references vs current-physics baselines).
- Digest <= 400 lines. Present the disposition table in chat too, not only in the file.
