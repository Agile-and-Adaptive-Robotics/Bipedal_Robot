# Round-4 review of Revision C — 2026-10-23

Verdict: **adopt as the round-4 base.** Revision C executes the round-3 direction
faithfully: the appendix relocations are real moves with provenance preserved, not
deletions or pointer-swaps; the build is clean; and the riskiest new claims check
out under independent verification. A short correction list follows, plus the
gap-identification expansion Ben asked for. Scored per the summary.md ledger.

## Independently verified (this review)

- **Build**: `review/build.log` (5 passes) — 0 errors, 0 undefined
  references/citations. `Revision_C.pdf` = 254 pages.
- **Label survival**: every externally referenced round-2 label resolves —
  Revision C carries them as dual labels on the relocated sections
  (`sec:sns_tuning_methods`+`sec:sns_curriculum`, `sec:spiking_methods`+
  `sec:spiking_results`+`sec:disc_spiking`, `sec:gait_library_results`,
  `sec:walker_variant_results`, `sec:sns_controller_comparison`,
  `sec:crossengine_methods`, `sec:animatlab_validation` in Appendix D;
  `tab:sadb_afferents`/`tab:sadb_pathways`, `fig:sadb_clusters` in Appendix E).
- **Numeric provenance**: `review/numeric_change_review.json` (176 entries)
  records every numeral-bearing line removed or moved, including my round-2
  torque numbers (28.8–50.6 N·m, 259.9 mm) and the basin-gate values.
- **New caption numbers**: the Appendix E citation-view caption claims 57
  in-corpus citers and 13 cited papers for Grillner (1985). Verified exactly
  against `SADb_audit/export/sadb_cites.json` (57 records list
  recFA1PeRSPTv9POu; its reference list has 13 entries). **add score** —
  new quantitative claims in captions are checkable and check out.
- **Assets**: the new app screenshots (`sadb_pivot_app.png`,
  `sadb_clusters_app.png`), `extensor_cam_route.pdf`,
  `pulley_torque_resolved.pdf`/`pulley_margin_resolved.pdf`, and
  `design_workflow.pdf` all exist in `figs/`.
- **Round-3 checklist, all done**: Appendix B-equivalent relocation real
  (Appendices D and E); the spiking overclaim narrowed ("preserved a continuous
  motor-activation interface while changing rhythm timing and, in some runs,
  limb coordination" — and Appendix D states outright "This construction does
  not require the network dynamics or resulting movement to remain unchanged");
  the 55.6% arithmetic fixed ("The spiking period is 55.6% longer", stated
  twice); the −25 sentinel explained correctly in both the tuning text and the
  ablation caption ("numerically greater than the genuine-gait scores"); the
  bibliometric-viewer history and "queries that no bibliography can" gone; the
  app/database distinction now exact ("The application distributes paper
  records and derived indexes, rather than the complete five-table curation
  base"; "the non-paper tables have not yet been exported as complete record
  sets into the application").

## Scores for Revision C (labeled examples)

- **add score** — the walking-status sentence is now the best in the program's
  history: "These two studies distinguish reduced neural assistance from reduced
  mechanical assistance; neither demonstrates harness-free walking."
- **add score** — SNS main chapter at 1,244 words with five figures and the
  capability table still establishes every tool; the scoping line "A successful
  interface test does not establish walking performance" in the table caption.
- **add score** — the requirements narrowing kept the advisor's correction
  ("These requirements concerned integration for this robot project") while
  Chapter 7's opening now carries the advisor's outcome-first framing ("The
  outcome that organizes this discussion is that improved first-pass design...
  is now achievable").
- **add score** — the SADb chapter states the store/distribution split Ben
  ruled, scoped to what the app contains today.
- **minus score** — the completed-study drafting assumption is documented only
  in `review/06_Completed_study_wording.md`; no FLAG marker exists at any claim
  site in the .tex (grep confirms zero "FLAG" in chapters/). The summary.md
  precedent (FLAG-BEN blocks at the claim, gating release) exists precisely so
  a drafting assumption cannot silently survive into the defended PDF.
- **minus score** — one person regression: "The author's experience with
  \qty{10}{\mm} actuators" (10-introduction, motivation paragraph 2) against
  the 2026-10-22 first-person ruling. (The Appendix B grep hit is a false
  positive: "the authoritative model state".)
- **minus score (mine, not C's)** — round 2 presented the Excel "Torque actual"
  rows as final adjoint-reduced torques; Ben's clarification that rows 16–17
  are sanity checks means my 28.8–50.6 N·m paragraph asserted presentation
  status for sanity data. Revision C removed it correctly. Logged in the
  summary.md ledger this session.

## Correction list for round 4 (small, focused)

1. **Flag the drafting assumption at the claim sites** (the summary.md FLAG-BEN
   pattern): a `% FLAG-BEN:` comment block at the 30-torque physical-validation
   paragraph and at the Conclusion sentence "The isometric tests demonstrate
   that the redesigned configurations meet or exceed the human reference
   torques over most of the tested range of motion." The abstract's outcome
   sentence is pre-existing baseline text the advisor has seen; flag it too or
   leave per Ben's call.
2. **"The author's experience" → "My experience"** (10-introduction ¶2).
3. **Gap explicitness (Ben's request)** — see the next section.
4. **Ben's call**: the spiking basin-gate pass (1% jitter, 12/12, 0.89× swing)
   was dropped from Appendix D (removal logged). It was the robustness
   counterpart to the timing results and the hardware-transfer argument;
   consider restoring one sentence in Appendix D's MuJoCo spiking subsection.
5. **Advisor note**: INTRO-31's outcome-first sentence no longer sits in the
   introduction's Chapter-7 description (it moved to Chapter 7's own opening,
   which serves the intent better); worth a word to ajh26 so he sees his edit
   honored.
6. Bookkeeping only: Revision C's round-2 word-count baseline (9,958 torque /
   6,865 SNS) differs from my measured finals (9,154 / 6,794) — their counter
   was modified ("escaped-percent and end-of-file fixes") and their starting
   commit may include small post-session edits; the deltas are measurement
   artifacts, not content disputes. Their final counts are the ones of record
   for the C files.

## The knowledge-gap expansion (Ben's standing request)

The actuator and torque gaps are explicit and concrete ("blocked rational
actuator selection"; "rested on unverified calculations"). The two neural-side
gaps were softened into *requirements* during the compression, which is exactly
where the request bites:

- Intro, neural paragraph: "sensory-feedback design draws on physiological
  measurements distributed across species, preparations, and primary papers...
  A structured resource can make those measurements easier to locate" — this
  says what a resource would do, not what was missing. Proposed round-4 wording:
  the afferent properties a feedback designer needs are scattered across a
  century of primary literature in cat, mouse, rat, and human preparations;
  reviews organize slices of it, but none can be filtered by afferent class,
  preparation, or pathway, and no structured, searchable resource existed from
  which a robot designer could draw them with their experimental context.
- Same paragraph, second half: "adapting and verifying their interfaces for a
  BPA-driven biped remains a project-specific task" — again a task, not a gap.
  Proposed: no assembled pipeline coupled a BPA-routed robot model to a
  biologically grounded neural controller in a maintained simulation
  environment; the lineage's contact-driven controllers existed only in an
  unmaintained engine and had never been coupled to a human-scale
  musculoskeletal model or exercised under load.
- The four-item summary paragraph should return to the "no X existed" form
  (round 2 had it) rather than "four connected needs."
- Optional one-sentence gap statements at the heads of Chapters 5 and 6 —
  WITHOUT re-broadening the SNS requirements section, which the advisor
  deliberately narrowed ("huh?" on the original breadth). The intro is the
  right home for the field-level gaps; the chapter keeps the project-scoped
  requirements.
