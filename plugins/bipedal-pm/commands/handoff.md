---
description: Write/update CHATGPT_HANDOFF.md or read/verify CHATGPT_REPORT.md
---

# /handoff

Run the bipedal-pm **chatgpt-handoff** skill. Mode from $ARGUMENTS:

- `write <topic>` (default) - update `CHATGPT_HANDOFF.md` for the named topic following
  its existing structure: date + machine + next step, no reorganization, no clobbering
  other agents' sections.
- `verify` - read `CHATGPT_REPORT.md`, extract claims/numbers/files, spot-verify against
  the repo and AGENTS.md, and deliver the confirmed / unverified / contradicted lists.
- `catchup` - summary of what the other AI reports plus what needs Ben's ruling.
