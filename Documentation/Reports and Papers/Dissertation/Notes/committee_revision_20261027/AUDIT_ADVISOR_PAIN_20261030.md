# Advisor-pain audit, consolidated report — 2026-10-30

Source of truth: the 2026-10-29 23:52 Overleaf ZIPs (snapshots and SHA-256s in
Archive/audit_20261030/SOURCE_MANIFEST.json; dissertation builds clean at 306 pp,
0 errors; paper confirmed 33 pp by Ben's recompile — the 38-page ZIP log was stale).
Rubric: penalty weights mined from the ajh26 comment corpus (integrity-insult class
-15, tracked changes -10/file, "huh?" clarity -6, non-sequitur -6, repetition -5,
wall >450 words -4, awkward -3, float defects -5). Full rubric:
tmp/audit_20261030/RUBRIC_ADVISOR_PAIN_20261030.md.

IMPORTANT interpretation: scores are a pattern-density index calibrated to advisor
pain, NOT a quality grade. Every unit below 70 has a bounded, surgical fix list and
NO proposed number changes anywhere.

## Scoreboard

| Unit | Score | Wall count | Worst pattern | Fix scale |
|---|---|---|---|---|
| Ch 4 (30-torque) | 24 | 0 | first-use jargon (4) + lineage block | ~6 clusters, under a day |
| Ch 5 (50-database) | 29 | 1 (506 w) | ONE surviving I1 integrity-pattern line (L50) | ~9 sentences |
| Ch 6 (60-sns) | 36 | 1 (501 w) | jargon (4) + tools list x4 | gloss-and-delete pass |
| Ch 7 (70-overall) | 27 | 1 (901 w, sec 7.1) | ch.-re-list wall + 3 caption defects | split/trim + captions |
| App A | 95 | 0 | compliance-array dup with C | pointer-ize |
| App B | 83 | 0 | c1/V3 unglossed; AnimatLab catalogue dup with D | small |
| App C | 65 | 1 (583 w) | process narration + table verbatim-dup from ch 30 | bounded |
| App D | 45 | 0 | SIX figures with no callouts; 3 orphaned sections | mechanical |
| App E | 95 | 0 | tables lack \ref callouts | 2 lines |
| Paper Main | 3 | 4 (884/770/548/529 w) | walls + repetition + callout order | one revision pass |
| Paper Supplement | 44 | 0 | mid-doc bib; 5 uncited items; 2 contradictions vs main | bounded |

Appendices overall 77 (mean) / 66 (size-weighted). Paper combined 24.

## What is verified CLEAN everywhere

- T1: ZERO tracked-change artifacts in any audited file (the -10 class never fired).
- Numbers: every prose number reconciles with its table in paper and diss ch 4/5;
  no notation violations (\diameter used, no "plant", Modern Robotics forms).
- I1: exactly ONE survivor dissertation-wide (50-database L50); chapters 4, 6, 7 and
  all appendices are clean of the integrity-hazard class.
- All appendices referenced from main chapters; no orphan appendix at chapter level.
- Motor-schematic squares (open advisor thread) verified absent from ch 7's five figures.

## Cross-cutting pattern classes (the machine-learned part)

1. INTEGRITY-PATTERN (P1, fix first): 50-database.tex L50 "I retained responsibility
   for inclusion, classification, and deletion decisions after checking the source."
   Same muscle as the sentence that drew "question your scientific integrity".
   Rewrite: fold the fact into the methods sentence, drop the verdict framing.
2. ADDED-THEN-REMOVED NARRATION (his 5.9 comment), four survivors:
   ch 5 L16 (Zotero groups then narrowed), ch 6 L145 (narrow-ramp settings), App D L60
   (former startup comparison), ch 4 L342-348 (placement-search lineage block, also a
   spine-anatomy non-sequitur). All become: state criteria/final state once.
3. FIRST-USE JARGON (his "huh?" class), ~14 dissertation-wide: W2L (ch 6 L4), FVU
   (ch 4 L314), l_620 (ch 4 L83), Hunt model (ch 4 L470, also uncited), Mono/Bi-pulley
   algorithms + insertion-side measure (ch 4 L749 - causes apparent 85.56-vs-85.74
   contradiction), production configuration / recovered (ch 6 L196/L209), isobaric
   (ch 7 L37), pivot views (ch 7 L14), c1 and V3 (App B L201), OpenAlex + attached-PDF
   indicators (ch 5 L8/L22), rho collision (supplement L242).
4. WALLS (default remedy DELETE): ch 7 sec 7.1 (901 w), App C L385-396 (583 w),
   ch 5 results opener L20-28 (506 w), ch 6 sec 6.4.1 L71-93 (501 w), paper x4.
   Chapter 4 has ZERO walls - the compression campaign held.
5. FLOAT DEFECTS: App D six figures with zero callouts (the exact "figure?" pattern);
   ch 4 L476 results paragraph missing opening callout; ch 4 Table 4.11 dual-unit cells
   (+5.02% next to 0.0501 fraction); ch 7 two figures with captions BAKED INTO ARTWORK
   atop LaTeX captions (pipeline.pdf, design_workflow.pdf - need regeneration);
   ch 7 L47 caption sentence contradicts rendered panels (BPA layer separation).
6. REPETITION: tools list x4 in ch 6 (L4, L39, Table 6.1, L63); capability taxonomy x3
   (L67/137/234); FVU values x3 in ch 4 (L606/663 + table); abstract-to-intro
   contribution restatement in paper; App C tab:TorqueFit verbatim duplicate of ch 30.

## Ben-decision flags (NONE-flag-only; never auto-edit)

1. Ch 4 Table 4.11: harmonize % vs fraction columns? touches displayed values.
2. Ch 4 L66: tape-measure resolution/uncertainty - he rejected the uncertainty task;
   flag stands for the defense Q&A only.
3. Ch 5 L22: "Within the application's Papers category, containing 943 entries" -
   redundant but carries a number of record.
4. Ch 7 L47: verify BPA construction figure against vector art (caption may be
   reversed between panels A and B).
5. Paper supplement L82: vastus INTERMEDIUS table vs main text grading extensor vs
   vastus MEDIALIS - which is right?
6. Paper supplement L156: tendon 25.00 mm vs main's "effective 5.00 mm tendon
   parameter" - terminology or real discrepancy?
7. App A L198 vs App C L227: provenance drift (frame line 306 vs 305).

## Paper submission blockers (beyond prose)

- 13 figures + 4 tables = 17 floats vs Frontiers' 15-float limit (template note at
  Main L729): merge or promote two floats to the supplement.
- Supplement places \bibliography mid-document; move to end.
- Build warning: fig:extensor_misalign multiply defined - circular \externaldocument
  Main<->Supplement; verify build chain before submission.
- Citation order: Figure 11 first cited at L410 before Figure 8.
- Five supplement items never specifically cited from main (blanket pointer only).
- Duplicate bib entry hunt_modeling_2017-1 (uncited; safe to drop).

## Suggested pass order (defense 2026-11-06)

P1 ch 5 L50 integrity sentence. P2 the four process-narration sites. P3 the ~14
jargon glosses (mechanical). P4 walls (delete-first). P5 App D callouts + App C
dedup. P6 ch 7 caption repairs (needs figure regeneration for the two baked-caption
figures - MATLAB). P7 paper revision pass + submission blockers (post-defense OK
if time is short, but the float limit and supplement bib are quick).

Per-unit full findings, wall tables, and three-highest-risk sentences live in the
six agent reports (this session's agent outputs) and the scripts/wall measurements
in tmp/audit_20261030/. Per standing rules: no edits were made to the snapshot;
local ProofFinal was only synced for the ch 6 title (see below).

## Session notes

- Ch 6 title synced local<-live: "Tools for Neuromechanical Simulation"
  (live edit made after the 10-29 morning push; no cross-references break;
  rebuild verified 306 pp, 0 errors, identical counts).
- Ben's figure drags of 10-29 all landed (reverse_pulley_biomimetic.pdf et al.
  present in the live ZIP; compile complete).
- Paper 33 pp confirmed live by Ben's recompile; the 38-page figure in the ZIP's
  bundled Main.log was stale. Flag retracted.


---

## Round 2 (2026-10-30, later): Chapter 2 expansion — 15-background.tex

Agent audit complete. Score: **0/100 (floored; raw -128)** — but the deficit is
concentrated in cheap mechanical fixes; I1/I2/T1/X1 all clean.

Headline findings:
- Terminology first-use failures (his "huh?" class): BPA used at L8, defined at L50
  (42-line gap); PAM (L123) and RoM never introduced; Ia/Ib/II afferents used at L166,
  defined at L173 (a section late); fictive locomotion, V3, Renshaw, locomotor
  primitives never glossed.
- Cross-chapter duplication with the Introduction: the BPA paragraph at ch2 L50 is
  verbatim intro L8; the three-differences list at ch2 L61 duplicates intro L8; the
  ethics enumeration at ch2 L270 duplicates intro L6.
- Misdirected cross-reference: L244 points cat-hindlimb models at sec:bg_compneuro
  (74-word section containing neither citation; should be bg_neural/bg_sensory).
- Reproduced-figure license gaps: Kengoro (L19), Rybak Fig 1 (L216), Shevtsova
  Fig 5.8 (L226) lack the "Reproduced from ... under CC BY 4.0" pattern the other
  five carry.
- One wall: L284-290, 524 words closing the ethics section (DELETE ~100 weakest words,
  horse aside first).
- Duplicate callout sentences L34/36; double-captioned moment-arm figure caption (L148)
  points ch2 readers at a ch4 equation by number.
- Open verification: figs/15-background not in ZIP extract, so caption-vs-panel
  rendering unchecked; render pass against live ProofFinal pending.

Full findings/wall/terminology tables in the session agent report; scripts at
tmp/audit_20261030/scripts/wall_measure_15background.py.



---

## Round 3 (2026-10-30 evening): Ben's findings incorporated (verified against the 2026-10-29 ZIP snapshot)

Verdict codes: CONFIRMED (verified at cited lines), PARTIAL (exists but needs Ben's
wording/ruling), DIRECTIVE (new content or MATLAB work; no verification needed beyond
location). All are folded into the fix campaign; none touch numbers except where a
number is the subject (flagged).

### Dissertation chapter 4 (30-torque.tex)

| # | Ben finding | Verdict | Location / notes |
|---|---|---|---|
| B1 | Table 4.3: FLEX and EXT weights on separate rows | CONFIRMED | Cost tables at L361 (flexor terms) and L418ff; mixed cells at L425-426 ("$10^7$ FLEX, $10^4$ EXT"). Fix: split FLEX/EXT into separate rows (or visual equivalent) in both weight tables. |
| B2 | subsubsection "Constrained Knee-Route Redesign" has no number; section end too long | CONFIRMED | L374 subsubsection renders unnumbered (class depth). Options: promote to subsection or leave unnumbered by design - Ben pick. End-of-section text L440-460 overlaps audit's L461-468 run; trim with the L342-348 lineage cut. |
| B3 | Fig 4.10 legend northeast in panel B, display order like other torque figures | DIRECTIVE (MATLAB) | Ext10mm_52cm_current_publication.pdf regeneration; legend position 'northeast' panel B; series draw order consistent with ExtPin/FlxPin figures. |
| B4 | Table 4.10: margin as percent, drop plus sign, add flexor worst-case angle, designed-vs-required wording | CONFIRMED | Margin table L674-685: "Worst-case angle & \ang{-20.2} & --" (L682, flexor cell empty - value needed from Mesh_Optimization results or Ben); "+5.02 %" cell (L680) and 0.0501 fraction column - unify to percent, remove plus. "Minimum margin" (designed) vs "Required margin" (5% spec) relabel per Ben's reading - CONFIRM that reading is correct per campaign records. |
| B5 | Pressure headroom dual nature; min operating pressure ~200 kPa (manifold open/close) | DIRECTIVE (content) | Belongs at the cam/redesign pressure discussion (L687ff, 269-507 kPa range). The ~200 kPa floor is Ben's operator observation - needs either his data pointer or "in my operating experience" attribution per the number-provenance rule. |
| B6 | Fig 4.14 "Prediction at..." points before measured | DIRECTIVE (MATLAB) | Redesigned_flexor_pressure_publication.pdf generator; same file used in paper Fig 13 panel B - ONE regeneration fixes both. |
| B7 | Remove R_mount calculation/discussion; cut section to ~2 paragraphs; "RoM" not "angle range" | CONFIRMED | L763-777 (equations eq:reverse_pulley_mount, sqrt(5)F20 discussion) delete; Mono/Bi-pulley subsection L749-798 compresses to ~2 paragraphs; sweep for "angle range" -> "RoM" after the C1 fix list defines RoM. NOTE: this partially reverses the 10-29 directives pass which ADDED the full vector equation - Ben's new ruling supersedes. |
| B8 | Explain briefly what "graded" means | PARTIAL | L792 defines the grading inline (sign + magnitude vs same-named muscle on 5x5 grid). Definition exists but is dense; add a one-clause plain gloss at first use. ("graded range" as a phrase does not appear; the grading-over-range scheme at L792 is the target.) |
| B9 | 4x 20mm BPAs in parallel = one 40mm, in theory | DIRECTIVE (content) | Natural home: L792ff (the baseline substitution sentence names exactly this substitution) - add one sentence stating the theoretical force equivalence. |

### Dissertation chapter 5 (50-database.tex)

| # | Ben finding | Verdict | Location / notes |
|---|---|---|---|
| B10 | Fig 5.2 too much white space at bottom, crop | DIRECTIVE (image) | sadb_pivot_app.png (L41); crop bottom margin, keep 0.96 width. |
| B11 | 5.3.2's future-improvements content belongs in ch 7 | CONFIRMED (paper 4.3.2) | Mapping: Ben's "5.3.2" = paper Physical Validation subsubsection (audit-numbered 5.3.2, Main.tex L604-628). The hedge tails at L610, L612, L627 are the future-flavored lines; the L627 segue he praised stays. Diss ch 7 already carries the isokinetic/testing opportunities framing - port essence there if not already covered (it largely is, at 70-overall L57). |
| B12 | Monthly scraper is built (factual) | CONFIRMED gap | No scraper/automation mention found in ch5 text; wherever the text implies manual updates (or omits the scraper), state that a monthly scraper is built and running. Exact sentence TBD at fix time with Ben's one-line description of what it does. |
| B13 | Pictures closer to front | CONFIRMED feasible | Both floats sit at Results (L30-44); intro L4-16 is figureless. Move sadb_bubble callout+float reference into the opening overview so Fig 5.1 lands by p.2 of the chapter. |
| B14 | Advisor comments today dispositioned | SATISFIED by this list | Ben's findings ARE the dispositions for today's comments; ignored/later items per his rulings (FVU stays ignored; motor-schematic thread he handles). Browser panel re-verification deferred (channel down). |

### Appendix C (94-AppendixC.tex)

| # | Ben finding | Verdict | Location / notes |
|---|---|---|---|
| B15 | Table C13 redundant with text | CONFIRMED | C13 = L456 configuration-survey coordinates table (sec C.10). Overlaps ch4 L792ff survey prose + Table tab:survey_dofs_sizes pointer. Keep the table as the data of record OR fold to supplement - Ben pick; default delete-per-redundancy. |
| B16 | C.9.1 duplicates ch 4 | CONFIRMED | C.9.1 "Goodness-of-Fit Metrics" L429-434 re-defines RMSE/FVU/Max already given in ch4 (first expansion currently in a ch4 table caption L486 - consolidate BOTH by defining once in ch4 prose and pointing from C.9.1). |
| B17 | C.10 jargonistic; C14 lacks master's-project mention | CONFIRMED | C.10 L452-497 + Table C14 (L502). Zero master's/Morrow/bolen_determination_2019 mentions in AppC (grep clean) while ch4 L792 credits the master's model. Fix: plain-language rewrite + one attribution sentence. |
| B18 | Tables C5/C6 to meters (C1/C4 already meters) | CONFIRMED | C5 = L120 (bio extensor inputs), C6 = L133 (measured dimensions; mm values 561/560/828/826 per measured-dimensions record). Convert to m for consistency with C1 (L23) and C4 (L105). NUMBER FORMAT change only - values unchanged, Ben has ruled. |
| B19 | C.4 repeats chapter equations | CONFIRMED (mild) | C.4 "Stiffness Arrays and Compliance" L205-257 contains 1 equation block duplicating ch4's compliance relation; replace with pointer. |

### Paper (Main.tex)

| # | Ben finding | Verdict | Location / notes |
|---|---|---|---|
| B20 | Intro->Methods segue: "the force model must see beyond" wrong register; ignores 27-actuator results | CONFIRMED | L101 exact sentence: "the excess contraction the force model must see beyond the actuator's physical contraction". Rewrite in neutral register; rebalance contribution framing to headline the redesigned-route + transmission-survey outcomes (27 sites, L643ff) alongside the corrected model. |
| B21 | Dissertation abstract to one page; propose rewording | PROPOSAL below | Current 606 words / ~2.3 pp double-spaced. Proposal drafted (this file, next section). |
| B22 | Old pulley schematic -> replace | CONFIRMED | Paper L401 includes Figures/reverse_pulley_principle.pdf (old, carries superseded F_BPA(n-1) labeling). Replace with reverse_pulley_biomimetic.pdf (regenerated with F_R label per B31) - file must be added to paper project. |
| B23 | Pulley discussion after meet/exceed section | CONFIRMED | Figure block currently in Methods L401-404; pulley discussion text already at L633 (Whole-Leg 4.4). Move figure+its methods text to 4.4, after the meet/exceed results (4.3.2). |
| B24 | Morrow 2020 cost function: mention, not spell out | CONFIRMED | 4.3.1 (L581-603) spells out J_flx/J_ext weight lists (align block L582-583+); compress to one sentence + citation of the original formulation (Morrow 2020 = master's cost; cite bolen_determination_2019/morrow key per Ben's bib preference). Also resolves audit's 770-word wall. |
| B25 | Fig 11 hybrid removed (plot, caption, discussion) | CONFIRMED | Diss version already regenerated (10-29); paper still has hybrid in figure + L499ff discussion + caption L505. Replace figure file with the regenerated no-hybrid PDF + scrub caption/prose. |
| B26 | 4.3.1 too long | CONFIRMED | Same section as B24; the audit's L592-603 770-word wall. Fix via B24 + provenance-clause deletions. |
| B27 | 4.3.2 important results first | CONFIRMED | L606 opens with setup, callout mid-paragraph, headline result at L610. Reorder: figure callout first sentence, then meet/exceed headline, then details (same pattern as diss ch4 round-2 fix). |
| B28 | Pulley schematic belongs in 4.4 | CONFIRMED | Same as B23; 4.4 = Whole-Leg Capability (L629), which is where the transmission content lives. |
| B29 | Fig 13: "Cam extensor" awkward; redraw both panels same linewidths/fonts; panel B prediction-before-measured; MFA 0.75, MEC none | CONFIRMED | L616-625: panels are two separately generated PDFs (Ext20_current_publication + Redesigned_flexor_pressure_publication). "Cam extensor" appears in caption L623 AND prose. Regenerate as matched pair (single script/style); rename panel in caption+prose (e.g. "cam-assisted extensor route"). |
| B30 | 4.4 cost function given there or in supplement + supp table | DIRECTIVE | The transmission-search cost currently only summarized; add the cost function either in 4.4 or as a new supplement section with its own table. Needs Ben's choice + the cost expression source (PULLEY_SWEEP directive records). |
| B31 | New schematic R -> F_R | DIRECTIVE (MATLAB) | Regeneration spec for reverse_pulley_biomimetic.pdf: reaction-force label F_R (not R / R_c). |

### Supplement (frontiers_SupplementaryMaterial.tex)

| # | Ben finding | Verdict | Location / notes |
|---|---|---|---|
| B32 | Table S1 redundant with text | CONFIRMED | S1 = tab:FlxPinPreGof (L47-66), pre-optimization GOF stats duplicating in-text pre-optimization discussion + table; candidate: fold into the main table or delete. |
| B33 | Fig S2 (unused 42 cm pinned extensor) not journal-worthy | CONFIRMED cut | fig:extensor_misalign (L186-191) shows the excluded 42 cm test; the exclusion is already truthfully stated in text. Cut figure + its referencing sentence. |

### One-page abstract proposal (B21; for Ben's go)

Current: 606 words, ~2.3 double-spaced pages. Target: ~250 words, one page.

PROPOSED (235 words):

"Biomimetic robots offer a platform for studying neuromechanical control through
experiments that would be impractical or unethical with living subjects. This
dissertation advances a long-running effort at Portland State University to build a
trunk-down humanoid robot actuated by pneumatic artificial muscles and controlled by a
synthetic nervous system. First, the isometric force of Festo braided pneumatic
actuators (BPAs) was characterized across resting lengths, contractions, and
pressures, producing new maximum-force equations and a three-coefficient surface that
predicts isometric force from relative strain and relative pressure. Second, the
actuator model was carried to the joint level through isometric knee-torque
experiments on a pinned knee and a biomimetic knee with a migrating center of
rotation. A hybrid calculation isolated the prediction errors, and multiobjective
identification of length-offset, compliance, and wrapping-loss correction terms
substantially improved torque prediction across configurations. Redesigned flexor and
cam-extensor routes, each driven by two BPAs, met or exceeded selected human isometric
torque references at every measured angle; these discrete measurements do not
establish dynamic performance. Third, supporting tools were developed: a curated
sensory-afferent literature database with a browser application, and neuromechanical
simulation tools that couple synthetic nervous systems to human and robot
biomechanical models, verified through staged tests separating implementation
agreement, sensory contribution, and mechanical response."

Every claim is already established in the current text; nothing new added; the SNS
verification clause reflects the ch-6 framing without new claims.

### Session notes for this round

- Ben's item mapping resolved: his "5.3.2" = paper Discussion subsubsection "Physical
  Validation Against Human Torque" (audit numbering 5.3.2; Frontiers numbering 4.3.2).
- Paper figure numbering: Fig 11 = biomimetic-extensor (Ext10mm); Fig 13 =
  physical-validation pair (Ext20 + redesigned flexor). Verified by position.
- B7 partially supersedes the 10-29 reverse-pulley directives (full vector equation
  now to be cut); recorded so the ledger does not treat this as a conflict.
- The browser review panel could not be re-read tonight (channel down); B14 satisfied
  by treating Ben's list as the disposition record, per his framing.



---

## Round 4 (2026-10-30 night): Evaluation of the other audit (output/pdf/Dissertation_and_Torque_Paper_Audit_20261030.pdf)

Overall: this is a strong, complementary audit. Where it overlaps my rubric audit it
agrees (ch 7 walls, ch 6 repetition, App C history, hybrid carryover, pivot screenshot,
figures-to-front). Its distinctive strength is PROVENANCE: it chased definitions,
saved records, rendered font sizes, code, and bibliographic metadata — catching five
classes my pattern-density audit did not. Its weaknesses: no chapter 2 coverage, no
terminology first-use sweep, no tracked-change scan of sources, and two structural
recommendations that must yield to Ben's explicit rulings. Adopt: all five P1 items
and most P2s. The merged master list is at the end of this section.

### Per-item response (every item of the other audit)

| ID | Their finding | My response / merit |
|---|---|---|
| P1 A1 | Appendix C says identified flexor uses two-rotation frames, then pitch-only; App A + Table C9 repeat the latter | ADOPT — highest-value find. My appendix agent caught only the A/C line-305-vs-306 drift. Resolution exactly as they say: align saved-fit vs current-evaluator descriptions from the dated record; no refit. |
| P1 T01 | +30.83%/+39.77% margin pair violates m_h=1.05(1+m_req)-1 if same domain | ADOPT as provenance check. Likely resolution: the two minima are taken over DIFFERENT curves/envelopes (required envelope vs plain target) so the argmins differ; document domains explicitly, change NO number. |
| P1 T02 | Paper defines r, F as magnitudes but introduces signed Mz = rF | ADOPT — my paper agent missed it. Fix: magnitude form as ch4 L78-82 already uses. |
| P1 T03 | Supplement has two unreconciled coordinate sets for the named biomimetic flexor | ADOPT — deeper than my agents' vastus/tendon flags (three distinct supplement-consistency items now: coordinates, vastus int/med, tendon 25-vs-5). Needs frame/version mapping; Ben's data. |
| P1 SNS-S01 | Gait-score description ("zero perfect, mean-removed") conflicts with current kine_ref.py (mean-angle term + bilateral bonus 40); tables show -145..-10,040 with no equation | ADOPT — defense-question hazard my prose audit could not catch. Document the archived scorer per result family, equation or term table, sentinel separately. No reruns. |
| P1 T04 | Bundled paper PDFs stale (38/5 pp vs 33/10 source) | ADOPT — consistent with Ben's recompile confirmation; exclude stale PDFs from any share package. |
| SADb S1 | Old criticized Table screenshot still the appendix demo | ADOPT — matches the known open advisor item and my report; needs Ben's current filtered capture. |
| SADb S2 | "38 records without grounding text" conflates empty notes with absent sources | ADOPT — 5.13-pattern class, cheap fix after checking the saved audit; complements (does not replace) my ch5 L50 integrity-sentence finding. Both stay. |
| SADb S3 | Pivot screenshot unreadable at final size | ADOPT — extends Ben's B10 (crop): may need column subset or native table; agrees the crop alone may not suffice. |
| SADb S4 | Bubble-map callout CONTENT too small | ADOPT (moderate) — pairs with the standing 10-pt floor rule; selective detail crops. |
| SADb S5 | Counts repeated x5; first visual p.124 | ADOPT — matches my ch5 R1s + Ben's B13 verbatim. |
| SADb S6 | No retrievable release/snapshot identifier for the app | ADOPT — one availability statement; Ben decides release venue/identifier. NEW reproducibility catch. |
| SNS-S02 | Duty caption contradicts zero-cycle row (Table D5) | ADOPT — cheap caption fix after definition check. |
| SNS-S03 | Old startup comparison presented as "matched/settled" | ADOPT — extends my AppD L60 finding; rename to final-window stats; state convention. |
| SNS-S04 | Fig 7.1 workflow lacks the second acceptance gate; failed loop leads straight to hardware | ADOPT — sharp catch; figure edit (add return route) or prose narrowing. |
| SNS-S05 | Figs 6.3/6.4B labels too small; B tiny in its frame | ADOPT — adds two items to the MATLAB regeneration queue. |
| SNS-S06 | Ch 7 walls (7:8/12/65/75) | ADOPT — agrees my ch7 audit (901-word 7.1); same delete-first remedy. |
| SNS-S07 | "Movement is not verification" repeated | ADOPT — = my ch6 taxonomy-x3 cluster. |
| SNS-S08 | Configuration map table (model/controller/support/result) | ADOPT — excellent, cheap defense protection; new appendix table. |
| SNS-S09 | "Accurately enough to commit" overclaims | ADOPT — one-sentence scope fix to first-pass static sizing. |
| SNS-S10 | Fig 6.8 lacks the phase annotation that IS the result (0.498->0.223 vs 0.35-0.65 band) | ADOPT — annotate median phase + criterion band, matched axis limits; MATLAB task. |
| SNS-S11 | Debug bookkeeping in App D ("dead", aliases, "zero hand drive") | ADOPT — merges with my AppD C1/A1 items. |
| T05 | Bib still calls the Actuators 15(5):230 paper unpublished | ADOPT — verified vs DOI; fix the paper bib entry. |
| T06 | Paper Fig 11 hybrid carryover | ADOPT — = Ben's B25 exactly; both manuscripts get the regenerated no-hybrid figure. |
| T07 | Legend text 5-6 pt at final size (Figs 9, 10, S1) | ADOPT — measured vector-text sizes; merges with Ben's B3/B29 redraw specs and the 10-pt floor. |
| T08 | Move cost definitions to Methods, hardware results to Results | PARTIAL — the 4.3.1 compression and cost-to-Methods/supplement parts match Ben's B24/B26/B30 and are adopted; the full move-hardware-evidence-to-Results restructure EXCEEDS Ben's directives (his diss keeps interpretation in Discussion) — hold for Ben's explicit choice. |
| T09 | Dense numeric prose -> compact tables | ADOPT — matches my walls + Ben's table preference. |
| T10 | Tendon 25 mm vs effective 5 mm undefined relationship | ADOPT — = Ben's flag; their remedy (define physical vs nominal vs effective parameter, preserve numbers) is right. |
| T11 | Author-list mismatch (7 in diss note vs 4 in paper) + title drift | ADOPT — human decision; correctly not an agent call. |
| T12 | Abstract overstates extensor-only first measurements; flexor has 15 points | ADOPT — wording fix; state both routes tested. |
| T13 | Equivalent-algebra repetition (Hooke chain) ch4 + paper | ADOPT (moderate) — compress substitutions inline, keep advisor-requested mechanics. |
| T14 | Fig 4.12 missing axis labels/units on A/H panels | ADOPT — adds to regeneration queue (my visual pass did not cover this figure). |
| T15 | Supplement inventory sentence promises a table that is not there | ADOPT — = my agent's uncited-items cluster; cheap fix. NOTE CONFLICT: they recommend KEEPING table S1; Ben's B32 says S1 seems redundant. Ben rules — S1 cut or folded per his go. |
| A2 | Fig B1 caption letters vs actual rows | ADOPT — correct the caption to (A)/(B)/(C) rows. |
| A3 | AppC C1/C2 legibility; C2 delete candidate | ADOPT — matches my agent's C1 letterless-panel finding; C2 crop-or-delete accepted. |
| A4 | AppC history/constants repetition | ADOPT — = my agent's I2/R1 wall cluster; same remedies. |
| A5 | AppB GUI-acceptance bookkeeping | ADOPT (moderate). |
| Rubric | Their assistant-failure rubric (40/30/25/20/12/10/8/6/2) | RESPECT — different object than mine (scores the assistant, not the manuscript); their calibration ledger is honest and complements my pain-weighted manuscript scores. Keep both. |
| Frontiers check | NO 15-float cap for Original Research; 12,000-word max instead | ACCEPTED CORRECTION OF MY AUDIT — my paper blocker "17 floats vs 15 limit" is retracted pending article-type confirmation; the 12k-word limit and ordered-figure/axis-unit rules replace it. |

### Corrections to my own audit from this comparison

1. Float-limit blocker retracted (see above).
2. Five catch-classes added to the master list: margin-relation provenance (T01),
   sign-definition (T02), coordinate/version mapping (T03+A1), score definition
   (SNS-S01), release identifier (SADb S6).
3. The "AI-tool rating methodology" open item is WITHDRAWN per Ben's clarification
   (his judgment, no documented method) — the other audit's withdrawal is correct and
   supersedes the CH45_POLISH note.
4. Their supplement-figure S2 ("slippage photograph") and Ben's "S2 = unused 42 cm
   extensor test" may be the same figure (extensor_misalign); Ben's cut ruling wins;
   verify identity at fix time.

### Where my audit stands (items the other audit did not cover)

Chapter 2 entirely (0/100 unit: terminology first-use failures, intro duplication,
license gaps, misdirected ref); the dissertation-wide first-use jargon sweep (14
sites); the tracked-change artifact scan (clean); the ch5 L50 integrity-pattern
sentence; the added-then-removed narration class (4 sites); reproduced-figure
licensing; and the wall measurements per float gap in chs 4/5/6.

### Merged master fix list (union, ordered for the campaign)

P1 provenance (before any prose pass): A1 frames, T01 margins, T03 coordinates,
SNS-S01 score, T02 sign, T05 bib, T11 authors (Ben).
Then Ben's directives (B1-B33) + my P1-P3 (integrity sentence, narrations, jargon)
+ adopted P2s (S01 scope sentence, S10 phase annotation, S04 workflow gate, S02/S03
AppD captions, S06 release statement, T12-T15, A2-A5), with MATLAB regeneration
queue consolidated: Figs 4.10/4.12/4.14, paper Figs 9/10/11/13 + S1, new pulley
schematic (F_R), Figs 6.3/6.4/6.8, 7.1 workflow, 5.2 crop, SADb captures.



---

## Round 5 (2026-10-30, final): Evaluation of the follow-up audit (Dissertation_Audit_Followup_20261030.pdf)

The follow-up is a parallel execution of the same three ordered tasks plus an
adversarial cross-review of MY audit (both frozen versions of it). Its quality is
high and its criticisms of my Round 3 are mostly correct. I arbitrated every
factual dispute against primary evidence (compiled .aux table numbers, file hashes,
PDF text extraction, and chapter structure). Results below.

### Arbitration: follow-up RIGHT, my Round 3 corrected (5 identity errors owned)

| Dispute | Evidence | Verdict |
|---|---|---|
| Ben's "5.3.2" mapping | 50-database.tex HAS subsections: 5.2.1 AI-Assisted Triage (L48), 5.3.1 Uses of the Database (L53), 5.3.2 Further Improvements (L61) | FOLLOW-UP RIGHT. 5.3.2 = diss ch5 "Further Improvements"; its content moves to ch 7. My paper-4.3.2 mapping retracted. |
| Table C13 identity | Build aux: C13 = tab:TorqueFit (L435) = biomimetic GoF duplicate of ch30's tab:BioTorqueGof; C14 = survey table (L456); C15 = capability (L502) | FOLLOW-UP RIGHT (matches my own appendix agent's duplicate finding, which I failed to connect). My B15/B17 line mappings were off by one; corrected. |
| Tables C5/C6 identity | Aux: C5 = bio flexor inputs (L120), C6 = bio extensor inputs (L133), C7 = measured dimensions | My B18 line identities wrong (I attached measured-dimensions to C6). Conversion directive stands; targets are the two bio-input tables. |
| Fig 4.14 state | Hashes: diss asset 0721D2DD/22,774 B (regenerated 10-29) vs paper 79814B2A/24,389 B | FOLLOW-UP RIGHT. Diss Fig 4.14 already prediction-first; the PAPER's copy is stale. B6 applies to the paper only; fix = copy the regenerated PDF into the paper project. |
| Abstract word count | Machine count of my proposal: 206 words | FOLLOW-UP RIGHT; my "235" label wrong. Their two-page compile claim is dubious (206 words double-spaced is ~1 page) - settle with a compile proof at implementation; trimming to ~160 is the safe target either way. |
| T01 argmin explanation | Their affine rebuttal is mathematically correct: pointwise min preserves an affine relation | MY EXPLANATION WRONG, accepted. The discrepancy requires genuinely different domains/targets/definitions from the saved records, not different argmins. |
| Callout-first rule | Ben's own B27 = claim-first | FOLLOW-UP RIGHT that "open with figure callout" was overstated; correct rule = claim first, callout early. |
| design_workflow baked caption | pdftotext of the PDF: only diagram-box text, no caption sentence | FOLLOW-UP RIGHT; my ch7 agent's second baked-caption F1 retracted. |

### Arbitration: follow-up overcorrected (2 items)

1. **pipeline.pdf double-caption IS REAL.** Text extracted verbatim from the figure:
   "Computer-in-the-loop: test embedded control against the simulated body before
   transfer to the physical robot." sits inside pipeline.pdf beneath the LaTeX
   caption at 70-overall L91. Their blanket "double-caption claims unsupported" is
   wrong for this figure; only the design_workflow half of my agent's claim was
   wrong. Pipeline keeps its FIX-CAPTION (regenerate without the baked line).
2. **Ch5 L50 is not the excluded AI comparison.** The retained-responsibility
   sentence is fixable without touching the tool-rating sentences; Ben's own 2:57
   reply sanctioned rewriting that paragraph. Finding stands, fix scoped to the
   responsibility sentence; ratings text untouched per the exclusion. Their R4-X04
   concern is noted but conflates the two.

### Their new contributions (verified where checkable)

- Ch2 science findings CONFIRMED at source: DEAs "kept them out of legged robots"
  (L30) vs DEAnsect (Science Robotics 2019) counterexample; asano_design_2017 cited
  among "humanoids" BPA examples (L61 tail) while L8 correctly describes those robots
  as brushless-motor-driven - a real citation-consistency defect. Geyer/Markin
  taxonomy, ethics-prohibition scoping, and accuracy-language findings are
  well-evidenced primary-source checks (not re-verified line-by-line here; their
  citations check out).
- 42-card current-state disposition ledger; 159-word one-page abstract alternative;
  generator-agnostic figure advice; ch2 float-interruption (B2-L01) and Fig 2.7
  0.62-width findings; Morrow-cost attribution correction (the spelled-out paper
  costs are the current constrained costs, not Morrow's original).

### Net changes to the master fix list

1. B11 retargeted: diss ch5 5.3.2 "Further Improvements" content moves to ch7
   (paper 4.3.2 keeps only the results-first reorder).
2. B15 retargeted to the tab:TorqueFit duplicate (L435, delete/fold; ch30 is the
   table of record); B17's table is the C14 survey table (L456) with master's
   attribution added.
3. B18 retargeted to the C5/C6 bio-input tables (meters conversion).
4. B6 retargeted to the paper asset only (copy regenerated PDF; no MATLAB rerun).
5. Pipeline.pdf regeneration retained; design_workflow caption fix dropped; the
   design_workflow second-gate fix (their SNS-S04) confirmed by the same text
   extraction ("Meets targets? ... Closed-loop verification -> Commit to hardware"
   with no return path).
6. Ch2 adds six science findings + float-placement/sizing items to the queue.
7. Abstract: two candidate proposals (206-word mine, corrected; 159-word theirs);
   final wording and length settled by compile proof at implementation.

### Overall verdict

This follow-up is the strongest of the three audit documents: it corrected five
real identity errors in my Round 3 with primary evidence, added genuinely new
literature-accuracy findings my pattern audit could not produce, and its one
systematic weakness is occasional overcorrection (dismissing the real pipeline
double-caption; conflating the L50 integrity sentence with the excluded AI
comparison; a dubious two-page abstract claim). With this round, the three audits
now converge on a single evidence-checked master list. No further audit passes are
needed; the next step is implementation on Ben's go.

