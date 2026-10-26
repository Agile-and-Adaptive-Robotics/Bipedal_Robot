# Context backup manifest — 2026-10-24 (staged by ZCode on EB475WS4)

Ben's laptop is down; he is on mobile. This bundle insures the dissertation
revision context that exists ONLY as uncommitted files on EB475WS4.

## Contents

- `Revision_D/` — Ben's prose working copy (changed chapters: Introduction,
  Background, Database, SNS), copied WITHOUT its regenerable `build_*`
  directories. Includes `chapters/`, `figs/`, `review/` (before-copies,
  `changes.json`, `changes.diff`, `conservation.json`), `main.tex`,
  `thesis.bib`, `BolenFrontiers22.bib`, `Frontiers-Harvard.bst`,
  `psutex.cls`, `REPORT.md`, `ZCODE_HANDOFF.md` (the handoff prompt of
  record), `Revision_D_local.pdf`. Provenance:
  `Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\Revision_D`.
- `CHATGPT_REPORT.md` — the live advisor-thread matrix (Cxx/Txx/Rxx),
  copied from the repo root of EB475WS4.
- `connectome_audit/` — ChatGPT's 2026-10-23 connectome/optimization audit
  (`AUDIT.md`, `ZCODE_HANDOFF.md`, `inventory.json`), independently
  machine-verified on EB475WS4 on 2026-10-24 (accurate except AUDIT.md
  finding 2's "four phase cells and six NMF synergies" sentence, which
  misattributes syn6's six NMF channels to the Di Russo draft).

## Restore rules

- Extract OUTSIDE the repository (e.g.
  `D:\GitHub\_zcode_backups\diss_context_20261024\`) so the working tree
  stays clean. Ben commits; agents do not.
- `thesis.bib` traveled inside a binary zip: byte-safe. Never round-trip it
  through default-encoding PowerShell; Latin-1 byte-safe edits only.
- Overleaf remains read-only: no chapter uploads from this copy without a
  fresh live-source comparison.
- `build_*` LaTeX dirs were excluded; rebuild fresh (MiKTeX present on both
  machines at `D:\MiKTeX\miktex\bin\x64`).
