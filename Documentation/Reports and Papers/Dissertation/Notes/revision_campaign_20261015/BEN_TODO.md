# BEN_TODO — running tally of things only Ben can do

Living list; agent updates it after each round. Items are numbered for reference in chat.
Machine: EB475WS4. Last updated: 2026-10-18 (writing-quality campaign session).

## Do now

1. **Commit the working tree (GitHub Desktop).** The whole 2026-10-17/18 revision wave is
   uncommitted: the journal-register chapter passes, the Tables 4.2/4.3 corrections,
   the publication-note rescope, the abstract/acknowledgments naming fixes, the cleaned
   AGENTS.md, the regenerated figures, and the campaign notes
   (RULINGS_20261017_session.md). Commit BEFORE further pulls. (tmp\ is disposable by
   convention — the diss_baseline_20261017 snapshot in there is not for committing.)

2. **Overleaf re-push (drag-and-drop, your Chrome):** drag the changed chapter .tex
   files to the project root and the new figure files into the project's
   figs/30-torque/ folder (the agent will list the exact changed set in the final
   report once the verification build lands). Then in the Overleaf UI: the five
   obsolete chapter files (20-methods, 30-results, 40-ankle, 50-futurework,
   60-conclusion) may still need deleting if the earlier cleanup did not cover them,
   and click Recompile.

3. **Overleaf review UI: disposition the 15 pending tracked changes (T01–T15).** The
   pushed files already carry their content applied; accepting in the UI reconciles
   the record. Any REJECT re-diverges the trees — tell the agent and it re-syncs.

## Figure decisions (nothing auto-wired; verdicts only)

4. `Figures/30-torque/` now holds the regenerated candidates from this campaign, none
   wired into the .tex: `Vas_Pam_20mm_cam_winner_{1,2,3}_r2.png` and
   `route_check_run8{,b}_r2.png` (restyled to your figure spec: Tol palette, distinct
   markers AND line styles, 10 pt floor), each with a `.alt.txt` sidecar, plus
   `opensim_robot_isometric.png` — a 300 dpi isometric render of the converted robot
   model (candidate replacement/companion for the low-res master's-era leg figure).
   Say which to wire and where, or discard. Cam winner 4 already passed. The knee GIF
   first frame remains unusable as a print still.

5. **Airtable screenshot** for the database chapter (carried): drop your capture at
   Figures/50-database/airtable_base.png; caption is drafted and the agent wires it.

## Content rulings (small, one ruling each)

6. **Voice:** you predicted the advisor will want more first-person "I". It was NOT
   mass-applied (style contract forbids it, and a full-voice change days before
   submission is high risk). The advisor's accepted T03/T04 wording ("I used these
   tools to convert...") IS restored in the background as the single sanctioned "I".
   A light first-person pass over the introduction, discussion, and conclusion is a
   one-hour job on your nod — ideally after the advisor sees the current draft.

7. **Database numbers framing** (carried): Chapter 5 carries the retained corpus
   (561/523/86/258), the Discussion the audited corpus (943/833/95/363). Both are
   explicitly scoped; bless both or pick which leads.

8. **Graphical abstract placement** (writer flag): fig:bpa_graphical_abstract is
   Chapter 3's graphical abstract but sits in the Background chapter. Left in place;
   moving it is your call.

9. **24 inferred comment anchors** (carried): skim comment_ledger.json and confirm
   the inferred placements so the inference retires.

## Housekeeping

10. **Pinned-knee mats (your stated plan):** the clean-named records of record —
    `minimizeFlxPin10_results.mat` and `minimizeExt10mmX3_results.mat` — are now what
    the dissertation and every script cite. When you delete the other pinned-knee
    result mats, keep `FlxPinBPASet.mat` and `ExtPinBPASet.mat` (test data, still
    cited by Appendix C).

11. Old build dirs (build_wf, build_final, build_final2, build_r2 through the latest
    minus one) are deletable; keep the newest build of record.

12. Defense-deck wiring (15+2 new assets) and the GIF supplementary-media decision
    (carried).

## Notation (for reading agent reports)

- **C01–C46** = the 46 OPEN reviewer comment threads, in order, in
  [CHATGPT_REPORT.md](D:/Github/Bipedal_Robot/CHATGPT_REPORT.md).
- **T01–T15** = the 15 PENDING tracked changes (accept/reject in Overleaf's review UI).
- **R01–R32** = the 32 already-resolved comment threads (archive section of the same file).
- Session rulings of record (tables 4.2/4.3 folds, publication note, clean mats, χ4,
  SADb naming, "pool" wording, wean fold-in):
  [RULINGS_20261017_session.md](D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/RULINGS_20261017_session.md).
