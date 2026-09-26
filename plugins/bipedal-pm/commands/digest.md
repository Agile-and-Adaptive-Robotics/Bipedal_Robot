---
description: Compile a dated project digest - readmes, reports, git, figures - with dissertation-disposition verdicts
---

# /digest

Run the bipedal-pm **project-digest** skill: sweep the standing summaries, report
folders, git log since the last digest, and artifact mtimes (fan out read-only Explore
subagents for the big areas), then write
`Documentation\Project_Digests\digest_<YYYY-MM-DD>.md` with the workstream states and
the dissertation-disposition table (INCLUDE / ALREADY IN / SKIP / LATER). Present the
disposition table in chat, not just the file.

$ARGUMENTS may set the lookback window (default: since the last digest file).
