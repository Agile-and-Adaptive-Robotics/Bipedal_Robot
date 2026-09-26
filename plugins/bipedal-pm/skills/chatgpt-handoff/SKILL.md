---
name: chatgpt-handoff
description: Write or update CHATGPT_HANDOFF.md (the brief for other AI assistants when ZCode is unavailable) and read/verify CHATGPT_REPORT.md (their report back). Use when Ben asks to write a handoff, catch up on what another AI (ChatGPT/Codex) did, or verify a ChatGPT report's claims.
---

# chatgpt-handoff

Files: repo-root `CHATGPT_HANDOFF.md` and `CHATGPT_REPORT.md`. Standing convention:
both stay current whenever work is handed off across assistants. `AGENTS.md` is the
canonical machine/session memory - the handoff POINTS at it, never duplicates it.

## WRITE (update the handoff)
1. Read the current file first; follow its existing section structure - do not
   reorganize or reformat it.
2. Update only the sections the work touched; prepend or amend the relevant Active Work
   brief. Every entry carries: date + machine name (`hostname`) + the next concrete
   step, in one paragraph.
3. Never delete or rewrite another agent's sections; stale entries get at most a dated
   `(superseded <date> - see <section>)` note.
4. Machine-specific instructions must name the machine explicitly (EB475WS4 /
   easteregg2 / DESKTOP-5Q16KE9) - paths differ per machine (see AGENTS.md machine map).
5. Never put credentials, API keys, or license keys in these files.
6. Ben commits via GitHub Desktop; do not commit yourself.

## READ (verify a report)
1. Read `CHATGPT_REPORT.md` (or a `.bak` copy) end to end. Extract: claims, numbers,
   files created/modified, and requested next steps.
2. Spot-verify each load-bearing claim: does the named file exist with a plausible
   mtime? does the cited number actually appear in the referenced report/script? does
   anything CONTRADICT AGENTS.md facts (machine paths, bit-exactness contracts, Ben's
   standing rulings)? Produce three lists: confirmed / unverified / contradicted.
3. Fold durable outcomes into the right home: repo facts -> AGENTS.md (Ben commits),
   project status -> a /digest entry, open work -> the handoff's Active Work section.
4. Report back to Ben: the acceptance list plus anything he must rule on. Unverified
   claims stay labeled unverified - do not launder them into AGENTS.md.
