# BEN_TODO — running tally of things only Ben can do

Living list; agent updates it after each round. Items are numbered for reference in chat.
Machine: EB475WS4. Last updated: 2026-10-17.

## Do now

1. **Overleaf finish (upload dragged 2026-10-17 — two steps left, ~30 seconds in your
   Chrome):** (a) delete the five obsolete chapter files from the Overleaf project tree —
   20-methods.tex, 30-results.tex, 40-ankle.tex, 50-futurework.tex, 60-conclusion.tex
   (they are inert: the new main.tex does not include them, but they clutter the tree);
   (b) click Recompile and eyeball the first pages (the agent's local gate says 264 pages,
   0 errors, 0 undefined references, 0 missing figures; the in-app browser is not logged
   into Overleaf, so this verification is yours). Old figs entries from the pre-restructure
   layout can be deleted whenever; they are unreferenced.

2. **Overleaf review UI: disposition the 15 pending tracked changes (T01–T15 in
   [CHATGPT_REPORT.md](D:/Github/Bipedal_Robot/CHATGPT_REPORT.md)).** The pushed files already
   carry their content applied, so accepting in the UI reconciles the record. If you REJECT
   any, Overleaf's file reverts and re-diverges from the local copy — tell the agent and it
   re-syncs.

## Content decisions (small, one ruling each)

3. **Airtable screenshot** for the database chapter: drop your capture at
   Figures/50-database/airtable_base.png (Papers grid with the linked fields visible).
   Caption is drafted; the agent wires the figure in once the file exists.

4. **Database numbers framing:** ch. 6 carries the retained corpus (561 papers / 523 notes /
   86 afferent-classified / 258 translation notes) while the Discussion chapter carries the
   audited corpus (943 / 833 / 95 / 363). Both are correct and explicitly scoped, but pick
   which framing leads, or bless both as-is.

5. **"Faster than real time" (background, simulation-tooling paragraph)** contradicts the SNS
   chapter's honest "full-body co-simulation is not yet running." One-sentence fix ready;
   awaiting your nod since it weakens a claim.

6. **Tracked-change T03/T04 wording:** the advisor's accepted edit was "I used these tools to
   convert..."; the restructure pass reverted it to the passive. Restore the advisor's wording
   (one line) or ratify the passive.

7. **Figure QC verdicts (from the campaign, nothing auto-changed):** cam winner figures 1–3
   FAILED print legibility (curve-style collisions), winner 4 PASSED; route-check figures
   FAILED (sub-legible at print scale); the knee GIF first frame is unusable as a print still.
   Regenerate the failed ones with your scripts if you want them in the dissertation.

8. **Master's leg-model figure is intrinsically low-res** (the 2019 source itself is low-res;
   feet cut off in the original). Regenerable from code for print-crisp if wanted.

9. **24 inferred comment anchors:** the detached advisor comments were placed by inferred
   anchors (the ledger marks them "inferred"). A skim of
   [comment_ledger.json](D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20260915/comment_ledger.json)
   … correct path: Notes/revision_campaign_20261015/comment_ledger.json — confirm the
   placements so the inference retires.

## Housekeeping

10. **Commit the working tree** (GitHub Desktop): the whole revision campaign — restructured
    chapters, moved figures, thesis.bib, campaign notes — is uncommitted.

11. **Build-directory cleanup:** ProofFinal now holds build_wf, build_final, build_final2,
    build_r2, build_r3 output folders (only the latest, build_r3, matters; the rest are
    deletable, though build_wf/main.pdf was locked by an open viewer).

12. Optional, defense-adjacent: deposit the balance-platform demonstration videos (for the
    record, if the quantitative video analysis ever gets added) and the SNS slideshow GIF
    supplementary-media decision.

## Notation (for reading agent reports)

- **C01–C46** = the 46 OPEN reviewer comment threads, in order, in
  [CHATGPT_REPORT.md](D:/Github/Bipedal_Robot/CHATGPT_REPORT.md) ("#### C12 — ajh26 — ...").
  Example: C12 is the advisor's chapter-restructure comment.
- **T01–T15** = the 15 PENDING tracked changes (accept/reject in Overleaf's review UI).
  Example: T02 is "requirement" → "improvement" in the background.
- **R01–R32** = the 32 already-resolved comment threads (archive section of the same file).
