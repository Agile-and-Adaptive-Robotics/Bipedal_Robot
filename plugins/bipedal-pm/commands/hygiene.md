---
description: Read-only repo hygiene audit (strays, naming, large files) and reference-safe move planning
---

# /hygiene

Run the bipedal-pm **repo-hygiene** skill.

- No arguments: run the full read-only audit (`hygiene_audit.ps1 -Mode audit`) and
  present findings as an evidence table with a recommended disposition per item
  (KEEP+RELOCATE / DELETE / ARCHIVE / GITIGNORE / BEN'S CALL). Nothing is modified.
- With a name (e.g. `/hygiene MuJoCo_SNS`): run `-Mode refs -Name <name>` and report
  the reference-impact analysis plus the full move protocol (fix classes, AGENTS.md
  map update, smoke test) as a plan for Ben's approval.
