# Vibe-edit queue — Ben's catalogue-and-execute list

Opened 2026-10-28 (late). Ben gives "vibe edits" as he spots them; agents catalogue
them here verbatim, propose the exact wording, and execute only on his go (the batch
approval counts as the go for entries marked QUEUED below). Status: QUEUED | DONE
(kept, with before/after) | DROPPED (Ben rejected).

## Standing patterns behind the queue (Ben + advisor, 2026-10-28)

- A request for more detail means bounds and settings, not a paragraph.
- Introduce the concept or frame the results and discussion; do not spit out reams of
  data and figures without framing.
- Compression must keep the one-sentence engineering motivation. Never a bare "cam
  solver" mention; instead, for example: "Due to concerns about material fatigue and
  failure, I decided to reduce the force exerted by the BPAs and increase the moment
  arm to get the required torque."
- Do not signal incomplete work; frame open items as opportunities for researchers to
  build on.
- No step-by-step or corrective-action stories; the committee has about a week to read
  the dissertation.

## Entries

### VE-01 — 70-overall.tex: end-of-intro "remaining work" framing — DONE (batch 1)

Ben (2026-10-28): at the end of the intro paragraph of Chapter 7 (overall), wording
like "remaining work to be done" indicates incomplete work. "Instead we'd frame it as
opportunities for researchers to build work off of." Note: Ben's own 26 Oct tracked
change already converted one instance ("remaining tasks" -> "tasks researchers can now
take"); this entry covers the end-of-intro instance.

Executed 2026-10-29 (vibe-editor). One passage, `chapters/70-overall.tex` line 4.

- Before: "This chapter considers how the actuator models, torque corrections, sensory-afferent resource, and simulation tools support biomimetic robot design, and identifies the experiments and integration work still required."
- After: "This chapter considers how the actuator models, torque corrections, sensory-afferent resource, and simulation tools support biomimetic robot design, and identifies the experiments and integration work on which researchers can now build."
- Rationale: minimal clause-level reframe of the end-of-intro sentence so open items read as opportunities for researchers to build on, not incomplete work; first clause, numbers, labels, and figures untouched. Ben's 26 Oct tracked change ("tasks researchers can now take", line 63) verified untouched.

### VE-02 — 94-AppendixC.tex: "results are bounded" justification — DONE (batch 1)

Ben (2026-10-28): the agent-written line about the 48 cm test ("results are bounded")
is wrong. The 42 cm and 48 cm tests were not used because the brackets shifted in the
bolt holes and the bracket had to be redesigned — but we do NOT tell the audience that
story. Fix: minimal, honest exclusion statement (those tests were excluded from
identification; training/validation folds unchanged), no "bounded" claim, no
bracket-shift narrative.

Executed 2026-10-29 (vibe-editor). Two sentence-level passages, both inside this
entry's scope, `chapters/94-AppendixC.tex`:

1. Fold paragraph (line 264):
   - Before: "Tests 3, 4, and 9 (42\,cm, 42\,cm with tendon, and 48\,cm) were not used because the attachment bracket shifted during those early tests."
   - After: "Tests 3, 4, and 9 (42\,cm, 42\,cm with tendon, and 48\,cm) were excluded from identification."
   - Rationale: minimal, honest exclusion statement; the bracket-shift/redesign story is removed. The fold structure in the same paragraph is unchanged and matches the surrounding tables exactly (flexor: training 1 and 5, validation 2 and 4, test 3 not used; extensor: training 5, 6, and 7, validation 1, 2, and 8; Tables~`tab:app_gof` and~`tab:app_gof_ext_all`).
2. Results paragraph (line 308):
   - Before: "The 43\,cm-tendon extensor test retains an FVU above unity despite improvement over its baseline. These outcomes bound the accuracy of the shared correction model."
   - After: "The 43\,cm-tendon extensor test retains an FVU above unity despite improvement over its baseline."
   - Rationale: deletes the "results are bounded" claim Ben rejected. The factual observations of the unused 48\,cm test and the 43\,cm-tendon test, and their Table~`tab:app_gof_ext_all` callout, are retained.

Deleted material contained no numbers (nothing to flag verbatim); no labels, tables,
figures, equations, or recorded values were altered.


## Status 2026-10-29 (overnight pass)

- VE-01 (remaining work -> opportunities framing): DONE, 70-overall line 4 area + new isokinetic paragraph line 57.
- VE-02 (bounded -> excluded-tests truth): DONE, 94-AppendixC line 264.
- Round-5 ledger entry appended to SLIM_BATCH1_LEDGER_20261028.md (isokinetic move, FlxGrp rename, ch.6 sweep, Overleaf sync record).
- Ben-only items outstanding: (1) drag FlxGrp.pdf into Overleaf figs/30-torque (live compile shows missing-file until then); (2) the open advisor comment threads incl. motor-schematic square position; (3) torque-paper Overleaf push blocked by browser-channel failure tonight - local paper Main.tex (33 pp) is ahead of live (32 pp).



## Torque paper sync verified 2026-10-29 (evening)

- Live Overleaf paper (project 622e9536b8ea1e0a7f167cf5) extracted and compared against
  local Knee_Torque_Test\ProofFinal: Main.tex and frontiers_SupplementaryMaterial.tex are
  character-identical (95,356 and 23,058 chars after CRLF normalization; zero diff lines).
- The earlier "local 33 pp ahead of live 32 pp" staging note is superseded - it was based on
  checkpoint page counts, not a live read; the 10-28 18:29 local edits did reach Overleaf.
- No paper push required. Remaining Ben items: the four figure drags for the dissertation
  (reverse_pulley_biomimetic.pdf required; FlxGrp.pdf, Ext10mm_52cm_current_publication.pdf,
  Redesigned_flexor_pressure_publication.pdf replacements).

