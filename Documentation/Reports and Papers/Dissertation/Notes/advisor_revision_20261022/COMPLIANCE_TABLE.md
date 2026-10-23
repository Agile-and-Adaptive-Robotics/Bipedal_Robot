# Per-comment and per-directive compliance table — 2026-10-22 revision pass

Editing base: fresh live Overleaf source (2026-10-22 capture), preserved at
`baseline/chapters/`. Patches: `patches/*.patch` (unified diffs) and full files at
`proposals/*.tex`. Status labels: **applied locally, verified by clean build**
(unless noted "proposed" = awaiting Ben's disposition). Line numbers refer to the
PROPOSAL files.

## 10-introduction.tex (advisor threads INTRO-01..34; 20 tracked edits already applied on live)

| ID | Comment gist | Where / evidence | Status |
|---|---|---|---|
| INTRO-01 | ¶ identifies one gap; there are others | The gap sequence now walks all four (lines 8, 10, 12) and the summary sentence names them (line 14) | applied locally |
| INTRO-02..06 | tracked edits (additional gaps, grounded, design tools) | Applied on live before this pass; retained | applied (prior) |
| INTRO-07 | "does not provide their assembly" — expand | Line 12: assembly defined as "the conversion of the robot's own model into a simulation engine, the wiring of the rhythm and reflex circuits with their afferent feedback, and the closed-loop coupling of controller and mechanics in one runnable environment" | applied locally |
| INTRO-08 | "their" undefined | Same sentence: antecedent is now "The neuromechanical modeling literature provides the ingredients … but it does not provide their assembly" | applied locally |
| INTRO-10 | "Awkward" (fourth-gap sentence) | Old "Wiring realistic feedback … this one in the literature rather than in methods" construction removed; replaced by "One lies in the literature itself … The other lies in methods" (line 12) | applied locally |
| INTRO-11 | "hundreds if not thousands" — you know the count | Line 12: "a primary literature that the database of Chapter 5 alone audits at 943 papers and reviews" (matches SADb audited corpus) | applied locally |
| INTRO-13 | mirror gaps ↔ chapter order | Afferent gap (Ch 5) now precedes pipeline gap (Ch 6) in line 12; summary sentence line 14 lists all four "in the order the chapters of this dissertation address them": force, torque, afferent resource, pipeline | applied locally |
| INTRO-14 | Steele/Muscle Mutt connection unclear | Line 20: "every piece of it is hardware this dissertation builds on: the … knee linkage … is the knee of the robot leg whose torque transmission Chapter 4 measures, and … Muscle Mutt … demonstrated antagonistic BPA pairs under SNS control on a physical robot … The humanoid leg of this dissertation carries that line forward." | applied locally |
| INTRO-15/16/17 | "first two"/"third"/"fourth" undefined | Line 24 names the gaps: "addresses the actuator-force and torque-correction gaps experimentally, resolves the afferent-data gap …, and closes the simulation-pipeline gap" | applied locally |
| INTRO-18..20, 22, 24..29, 31, 34 | tracked edits already applied on live | Retained verbatim | applied (prior) |
| INTRO-21 | "Gait2392" too specific in roadmap | Line 39: "evaluated on the converted robot-body model" | applied locally |
| INTRO-23 | "nor any role for afferent feedback…" huh? | Line 43: "produced suspended stepping but never walked: the walker stepped only in air, never bore its own weight, and its afferent feedback was never exercised under load" | applied locally |
| INTRO-26/30 | "spiking-neuron mirrors" huh? | Line 43, first use now reads "spiking-neuron variants of the same networks that test whether conclusions depend on the graded-coding assumption"; later use "the spiking variants" | applied locally |
| INTRO-32/33 | "Its also scopes" typo | Line 45: "It also scopes" | applied locally |

## 15-background.tex (advisor threads BG-01..50)

| ID | Comment gist | Where / evidence | Status |
|---|---|---|---|
| BG-03 | motor center off (figure) | FIGURE item — figures were re-uploaded in the prior sync pass; no text action. Flagged for Ben in REPORT | awaiting Ben |
| BG-04 | graphical abstract "confusing and probably not needed" | Figure block + intro sentence REMOVED (proposal diff). Asset untouched in Figures/15-background/. **Flagged for Ben's final call** | proposed locally |
| BG-05 | what do 23 and 92 mean | §2.3 ¶1: "the Gait2392 model, named for its 23 degrees of freedom and 92 musculotendon actuators"; the Arnold 2010 clause reworded to "related lower-limb models of the same family are characterized and validated by" (that model's own 23-DoF/54-MTA count no longer misread as Gait2392's) | applied locally |
| BG-06/07 | "first-class" / awkward | "makes that adjacent-joint dependence a primary result of its attachment-point technique" (advisor's suggested word) | applied locally |
| BG-09 | "used throughout this dissertation" — not your collected data | Now: "the moment arms used in the robot model's design and torque predictions here are properties of the Gait2392-derived geometry rather than quantities measured on the robot; the torques measured on the test stands of Chapter 4 are the experimental ground truth against which those geometry-based predictions are checked" | applied locally |
| BG-10 | "definitions are stressed" huh? | "The human knee is the joint where the choice of these definitions matters most." | applied locally |
| BG-11 | run-on ICR sentence | Split into two sentences ("Locating that moving center is a kinematic problem in its own right, and [they] solve it with natural coordinates …") | applied locally |
| BG-13 | "written against" huh? | "were written to resolve exactly this definitional freedom, keeping reports of joint motion comparable across laboratories" | applied locally |
| BG-14 | which method does 2392 use? | New sentence: "Gait2392 resolves the knee with planar kinematics \citep{yamaguchi_planar_1989}: the flexion angle prescribes the joint's translation as well as its rotation, rather than fixing a hinge axis or solving a linkage." | applied locally |
| BG-15 | four-bar run-on(s) | Split at "The rolling-and-sliding behavior … therefore emerges …" | applied locally |
| BG-16 | non sequitur "exactly the behavior a tool must capture" | Now: "a torque prediction tool for this knee must capture that variation rather than assume a constant arm" | applied locally |
| BG-12, 17, 18, 20..23, 31..42, 45(applied), 46..49 | tracked edits already applied on live | Retained verbatim | applied (prior) |
| BG-19 | non seq (figure comparison) | "shows the two sides of that comparison: the laboratory's bipedal skeleton and the OpenSim muscle-routing reference developed for it" | applied locally |
| BG-24 | "Delete?" (robotbody sentence) | Compressed and kept (the fact is load-bearing): "the resulting robot model, gait2392_robotbody, is the robot model whose torque transmission the later chapters of this dissertation refine" | applied locally |
| BG-25 | "routings" → origin and insertion points? | "each PAM's origin, via points, and insertion were placed by hand in OpenSim beside the Gait2392 path it imitates" | applied locally |
| BG-26 | "moved to respect contraction limits" huh? | "adjusted so that the contraction each actuator must produce over the range of motion stays within the Festo BPA's usable contraction range" | applied locally |
| BG-27 | "saved as" — don't narrate saves | Removed ("the result was saved as" gone) | applied locally |
| BG-28 | "jumping all over the place" | ¶ re-sequenced: mechanics stated once → arms + torque method → what the comparison was → biceps example → routing placement → matching problem | applied locally |
| BG-29 | "followed" → is calculated | "isometric torque was calculated as muscle force multiplied by moment arm" | applied locally |
| BG-30 | "what comparison was made?" | "The comparison the study made was robot against human: …" | applied locally |
| BG-43/44/45 | CPG-specific material in comp-neuro section | Brown half-center + persistent-sodium sentences MOVED into §2.5 (CPG section, now its ¶1); the model-family ¶ moved into §2.5 as a condensed second ¶; §2.7 keeps the executable-models framing + the controls-tradition reading with an explicit scope pointer to §2.5 | applied locally |

## 30-torque.tex (reconcile with Oct-22 records)

| Directive (from Ben's instructions) | Evidence | Status |
|---|---|---|
| Resolve Main.tex:540–542 vs :556 vs :560 contradiction against test records | Record: fd1c13f4's six-point paragraph is REAL but superseded — 7862e4e2's `Knee_Extensor_20mm.m` and `Results_table_20mm.xlsx` (ExtTest20mm_2 tab) now carry NINE holds (ExtTest_2 six + ExtTest_3 three) with the adjoint reduction EXECUTED at reaction point p_rf = (−0.05, −0.255, 0) m, d = 259.9 mm; workbook has adjoint torques for six holds (28.8–50.6 N·m). Manuscript lines 556/560 remain stale — flagged to Ben | applied locally |
| Ch.4 must not say the full sweep is pending | Replaced: "the full sweep … pending" → nine-point cam-route sweep reported with raw loads/pressures/angles, six adjoint-reduced torques, reduction of three in progress, envelope comparison + flexor pending (proposals/30-torque.tex, §4.3 physical-validation ¶) | applied locally |
| Separate raw measurements / provisional reductions / finalized comparisons / optimizer predictions | The paragraph labels each class: raw nine holds (loads/pressures/angles), adjoint-reduced six (28.8–50.6 N·m, "entered"), three in progress ("in progress"), comparison pending; optimizer predictions unchanged elsewhere and still labeled optimizer-stage | applied locally |
| Current reaction-point constants + LC-angle convention | d = 259.9 mm and "load-cell angle recorded in the tibia frame" (Ben's 2026-10-04 convention); wrench carried "by the same adjoint transform as the identification data" | applied locally |
| Direct-route vs cam-route correction | The initial three points now read "on the earlier direct route" (fd1c13f4's accurate framing); the nine-point sweep is "on the rebuilt cam route" | applied locally |
| Provenance for every changed number | Loads/angles/pressures/adjoint torques: Results_table_20mm.xlsx ExtTest20mm_2 tab + Knee_Extensor_20mm.m lines 590–655 (commit 7862e4e2); pressure-band reading (256 vs 269–507; 480/488 vs 507; 620 ceiling): manuscript fd1c13f4 Main.tex:542 | applied locally |

## 50-database.tex (advisor threads DB-01..26 + Ben's SADb structure)

| ID | Comment gist | Where / evidence | Status |
|---|---|---|---|
| DB-01/02 | "design needs parameters, not only principles" huh? / relation to database unclear | ¶1 rewritten: "Building that controller requires numbers." + the specification enumerated as what the database carries | applied locally |
| DB-03 | ipsi/contra — more options | "whether the target is ipsilateral, contralateral, or both" | applied locally |
| DB-04 | "still lots of unknowns" | "much of it remains unknown or contested outside the pathway-rich core" | applied locally |
| DB-05 | define gap + what you did | Gap ¶ (¶1) → product ¶ (¶2) immediately after | applied locally |
| DB-06 | tracked edit applied on live | Retained | applied (prior) |
| DB-07 | readers sentence should lead | Readers/questions now the ¶1 spine | applied locally |
| DB-08 | reviews ¶ belongs in methods | Moved to Curation Methods ¶1 | applied locally |
| DB-09 | presented 2x (without/with AI) | VERIFIED from bib: Neuroscience 2023 poster (bolen_sensory_2023, the database) + Neuroscience 2024 poster (bolen_using_2024, the AI comparative study); both stated at first mention (¶2) and in the triage section | applied locally |
| DB-10/21 | "distilled mechanism notes"/"what are notes?" | Defined in product ¶2: "structured curation notes, one to three sentences of distilled mechanism written from the paper's own text, never from memory" | applied locally |
| DB-11/14 | "first-class"/"canonical" | Both wordings gone ("a record in its own right"; "the DOI as the unique identifier") | applied locally |
| DB-12 | "pathway records" huh? | Defined in product ¶2 with the naming grammar and an example ("Ib stance to swing") | applied locally |
| DB-13 | "I thought you migrated it" | Product ¶2/¶3: Airtable is the working home; the packaged offline explorer is the distribution; no migration claimed | applied locally |
| DB-15/16 | "one program" huh? / move to end | Cut from the top; the mapping to Chapter 6 closes the chapter | applied locally |
| DB-17/18 | "Methods" heading + Cochrane form | Real \section{Curation Methods} organized as sources → inclusion/exclusion rules → note-writing → AI-assisted workflow + decision authority → prune; no registered-review claim | applied locally |
| DB-19 | "audit"? | "audit" used with its content defined by the methods (three-source reconciliation, verification passes) | applied locally |
| DB-20 | too confusing, start from finished product | Chapter now opens with the finished product; assembly narrative compressed to sources + verified imports + growth (batch numbers dropped) | applied locally |
| DB-22 | big picture first | Product ¶¶ 2–3 precede all methods | applied locally |
| DB-23 | approved by whom? | "each approved by me before execution" + "I held every decision that required judgment or carried risk" (verified: SADb_audit README pilot-review requirement; prune_review_20261001.md/round2 are Ben's review docs) | applied locally |
| DB-24 | how decided/identified? | Inclusion scope stated positively; exclusion rules, weighted rubric dimensions, and protected classes listed; off-topic categories named | applied locally |
| DB-25 | "the laboratory" → take credit | "I built it within the AARL's … (C3NS) effort"; "my personal Zotero collection" | applied locally |
| DB-26 | "curation tool, not a reader" huh? | Recast: "Airtable remains the working home of the data, but it cannot travel past the laboratory's account" | applied locally |

## 60-sns.tex (advisor threads SNS-01..15 + Ben's five-part structure)

| ID | Comment gist | Where / evidence | Status |
|---|---|---|---|
| SNS-01 | so many details, what is the goal/gap? | Preamble replaced by 2-¶ opener + §6.1 the combination + NEW §6.2 "What the AnimatLab Attempt Exposed" with the four gaps from Ben's reply | applied locally |
| SNS-02 | "reflex-based and learned alternatives … treated there" huh? | Preamble ¶3 (containing it) removed; scope pointer survives as one clause in the opener ("defined in Sections 2.5–2.8") | applied locally |
| SNS-03 | "tooling" is a manufacturing term | "tooling" no longer used as the chapter's self-description (opener speaks of "the simulation tools required to develop it"; section titled "A Neuromechanical Simulation Toolchain") | applied locally |
| SNS-04/05 | tracked edits (virtual, this) | Retained where their host sentence survived; the host ¶ was replaced in the restructure — the "virtual laboratory" framing itself was cut per SNS-07 | applied locally |
| SNS-06 | why several engines? | §6.2 final ¶ answers directly from Ben's reply (MuJoCo speed/rigid-tendon limits; Simulink/Simscape graphical construction + CAD import; real-time/computer-in-the-loop) | applied locally |
| SNS-07 | identify the gap, save implications for discussion | The laboratory ¶ is gone; gaps up front; implications consolidated in §6.6 | applied locally |
| SNS-08 | "Preliminary Results:" deleted (tracked) | Retained (applied on live); section now part of the opening combination story | applied (prior) |
| SNS-09 | "verification-first approach" connection unclear | The connection is now causal: the hand-tuning failure of §6.1 ends with the sentence that motivates it; §6.2 opens by deriving the gaps from that failure | applied locally |
| SNS-10 | "build order" huh? | §6.3 intro no longer uses "build order" — "This section documents the tools in build order" → "in build order" retained? NO: reads "This section documents the tools in build order:" → replaced with "in the order they are built upon" is not present; final text: "the existing tools they stand on, the new tools built from them, and the methodology by which the controllers are tuned and checked" | applied locally |
| SNS-11 | overlapping text in figure | Figure was fixed in the prior sync pass (re-uploaded); no text action | applied (prior) |
| SNS-12 | "The laboratory stands on tools this program did not build" (tracked delete) | Retained deleted; existing-tools subsection opens with "None of the infrastructure components is invented here" | applied locally |
| SNS-13 | SADb segment paragraph unnecessary | Removed; single linkage kept at the §6.3 wiring point and in existing-tools ¶1 | applied locally |
| SNS-14 | need deeper description of existing tools' limitations | The limitations are now §6.2's content (four gaps, each naming what the existing tool cannot do) | applied locally |
| SNS-15 | compare against Hunt's and Fletcher's libraries | NEW comparison ¶ in §6.3 New Tools citing hunt_development_2017 and young_analyzing_2019 (real sources; both already in thesis.bib): what the predecessors carry (Hill-type muscles, animal-scale models) and what the derivative adds (Ia/Ib encoder blocks, diameter-specific BPA force elements, CAD/OpenSim body coupling, spiking additions, regeneration from the tuned network) | applied locally |
| Ben structure 1 | AnimatLab combination here, not Background | §6.1 with the combination named in ¶1 ("Morrow and I combined two published designs …: the biomechanics came from the suspended bipedal walker of Li et al., and the controller from the two-layer CPG organization of Deng et al.") | applied locally |
| Ben structure 2 | gaps from Ben's replies | §6.2, four gaps in Ben's order + the engines rationale (SNS-06 reply) | applied locally |
| Ben structure 3 | tools + brief results, capability/evidence table, status classes | §6.3 with Table tab:sns_capability (status column: complete / partial import / provisional / proposed) | applied locally |
| Ben structure 4 | verification: expected/test/observed/limit | §6.4 reframed ("each subsection stating the expected behavior, the test or comparison, the observed result, and its limit"); subsections rewritten in that form | applied locally |
| Ben structure 5 | short handoff; future work in Ch 7 | §6.6 "Delivered Capabilities, Limits, and the Program They Enable" — two ¶s + limits ¶ pointing to Chapter 7 | applied locally |
| Standing rule | retain "no harness-free walking demonstrated" | §6.6 walking-status ¶ keeps the exact statement ("Unsupported, harness-free walking has not been demonstrated anywhere in this program's records") | applied locally |
| Bug fix | dangling `Figure~\ref{fig:spike_animatlab}` (no figure block exists in baseline) | Reference removed (the table carries the comparison); logged | applied locally |

## Ben's global directives

| Directive | Status |
|---|---|
| Cut map first | CUT_MAP.md written before execution |
| Working goals 40–50% SNS / 25–35% SADb | Landed 27.5% SNS / 15.0% SADb — justified deviation: every verified number, figure, table, and externally referenced label is preserved and the AnimatLab detail Ben asked to retain is intact; the cut concentrated in the preamble, duplicated inventories/rationale, and discussion repetition (full accounting in REPORT.md §2) |
| Measure words/captions/pages before+after | measure_words.py results below; pages 294 → 284 (fresh-build comparison, same engine) |
| Small reviewable patches + move map | Per-chapter unified diffs (patches/) + SNS_MOVE_MAP.md (every paragraph dispositioned) |
| No upload of prose/reports to Overleaf; no comment resolution; no commits/pushes | None performed; everything is local under Notes/advisor_revision_20261022/ |
| Bibliography untouched | thesis.bib byte-identical to live baseline |
| Latin-1 byte-safe bib | No bib edits at all |

## Measurements (measure_words.py; narrative excludes captions and table bodies)

| Chapter | Baseline narrative | Proposal narrative | Δ |
|---|---:|---:|---:|
| 10-introduction | 2,093 | 2,202 | +109 (repairs) |
| 15-background | 5,217 | 5,392 | +175 (repairs; −1 figure) |
| 30-torque | 8,993 | 9,128 | +135 (new verified results) |
| 50-database | 2,588 | 2,201 | −387 (−15.0%) |
| 60-sns | 10,851 | 7,868 | −2,983 (−27.5%) |
| Total (5 chapters) | 29,742 | 27,391 | −2,351 (−7.9%) |

PDF: baseline build 294 pages → proposal build 284 pages (−10). Captions: 2,547 →
2,345 words across the five chapters (background −28 from the removed graphical
abstract). Chapter page spans (proposal TOC): Ch 5 pp. 120–132, Ch 6 pp. 133–180.

Verification build (MiKTeX 26.5, pdflatex ×4 + bibtex, fresh D:\temp dirs,
absolute BIBINPUTS/BSTINPUTS): baseline and proposal both exit 0 with 0 errors,
0 undefined references, 0 undefined citations, 0 missing figure files.


## Visual verification (ink-bounds analysis, 110 dpi renders; visual-judge unavailable, fallback per instructions)

Control pages (untouched chapters) define the normal text-block ink envelope:
ink x[110,826] px of a 935 px page, right margin 109 px (approximately 72 pt margins).

| Page (PDF) | Content | Result |
|---|---|---|
| 26-27 | Intro gap + objectives rewrite | Fits envelope |
| 142 | Torque nine-point sweep paragraph | Fits envelope |
| 143, 145, 146 | SADb opener, methods, Explorer figure | Fits envelope |
| 153-154 | tab:sadb_afferents + tab:sadb_pathways | Sideways tables (rotated, centered) — correct in baseline and proposal; NOT defects |
| 156, 158 | SNS opener, gaps section | Fits envelope |
| 164-165 | tab:sns_capability | DEFECT FOUND on first render (ink to page edge, 934/935 px; ~100 cell overfull boxes up to 397 pt) → FIXED: `llX` natural-width columns replaced by fixed ragged-right `p{0.28/0.24}+X` at `\small`. Re-render: ink x[110,826], identical to controls |
| 179 | Verification section | Fits envelope |
| 196 | tab:spike_animatlab | PRE-EXISTING 24 pt overflow (baseline PDF p.203 identical, spans to x=564.8 pt) → FIXED with single `\small` (natural width 493→450 pt; no text or number changes). Now x=542.1 pt, 1.4 pt residual |
| 201 | Controller comparison | Fits envelope |

Overfull-hbox width multisets (final pass): baseline 148, proposal 146
(2 baseline widths removed by rewrites; 0 added).

Rebuild after fixes: 284 pages, 0 errors, 0 undefined references/citations,
0 missing figure files. Proposal build log = `tmp/build_live4.log`; renders
`v3_p164.png` / `v3_p196.png`.


## Supervisor review round (diss-supervisor agent, four gates) — fixes applied

The supervisor returned BLOCK with three factual failures and proofreading flags; all
mandated fixes applied and verified in the rebuild (285 pp, 0 undefined, 0 missing,
overfull multiset still 0 added vs. baseline):

| Finding | Fix applied |
|---|---|
| 60-sns opener: "graded neuron and synapse models … defined in Sections 2.5–2.8" false (only CPG + afferents defined there; equations live in App. B); also "is defined" with list subject; hardcoded numbers | Reworded: vocabulary/CPG/afferents → `Sections~\ref{sec:bg_neural}--\ref{sec:bg_sns}`; neuron and synapse equations → `Appendix~\ref{app:neuron_equations}`; "are defined" |
| 60-sns tuning ¶: "detailed in Section~\ref{sec:gait_library_results}" dangled (composition deleted) | Restored the baseline's verified composition sentence at the target subsection (32 running / 8 subjects / 4 speeds; ten walking 0.5–2.0 m/s + 3 self-selected; Falisse cycle) |
| 10-intro: "evaluated on the converted robot-body model" mislabel (Ch. 4 evaluates on the OpenSim model; conversion predates redesigns) | "converted" dropped |
| 10-intro: "every piece of it is hardware this dissertation builds on" overclaim (Muscle Mutt is a demonstration) | "and this dissertation builds on its hardware and demonstrations" |
| 10-intro: "spiking-neuron variants" vs ch. 6 "mirrors" term mismatch | "spiking-neuron mirrors"; closing list "the spiking comparison" |
| 50-database: "the afferent classes that controller monitors" | "that the controller monitors" |
| 60-sns:32 "documents the tools in build order" (advisor SNS-10 "huh?") | "documents the tools and what each contributes" |
| 30-torque:588 verbatim manuscript frames (Gate 1 FLAG: 11/18/12/10-word runs vs bolen_isometric_2025 Main.tex) | Prose frames reworded (opening sentence, adjoint-transform clause, load-cell clause, ceiling clause); numbers and data recitation untouched |
| 30-torque:588 "at the low end of the 269–507 kPa band" (256 < 269) | "just below the low end" |
| 30-torque:588 mixed `\newton\cdot\meter` / `\newton\meter` in one sentence | unified to `\newton\meter` (chapter-dominant form) |
| 30-torque:7 "follow the manuscript" understates the newer sweep | added "; the physical-validation paragraph of Section~\ref{sec:disc_pulley} additionally reports the newer nine-point sweep results, which postdate the manuscript draft" |
| INTRO-31 tracked edit (outcome-first Chapter-7 sentence) | advisor's wording applied verbatim at 10-introduction.tex:45 |

Escalated to Ben (NOT decided here): Gate-1 disposition is now moot (frames reworded),
but the flag is recorded; BG-24 "Delete?" ruling; SNS-14 depth-vs-reduce ruling;
first-person sweep scope (chs. 1–2 still "the author") + style-contract.md §2
supersession note; BG-03 and SNS-11 figure-level fixes.

## Final measurements (post-supervisor-fix)

| Chapter | Baseline narrative | Proposal narrative | Δ |
|---|---:|---:|---:|
| 10-introduction | 2,093 | 2,229 | +136 (repairs + INTRO-31) |
| 15-background | 5,217 | 5,392 | +175 (repairs; −1 figure) |
| 30-torque | 8,993 | 9,154 | +161 (verified sweep + reworded frames) |
| 50-database | 2,588 | 2,202 | −386 (−14.9%) |
| 60-sns | 10,851 | 7,910 | −2,941 (−27.1%) |
| Total (5 chapters) | 29,742 | 27,887 | −1,855 (−6.2%) |

PDF: 294 → 285 pages. Verification build of record: `tmp/build_live5.log`
(patches regenerated 2026-10-22 against the same fresh live baseline).


## NOTE (2026-10-22 round 2): measurements above are ROUND-1 numbers

Round 2 applied the principal's selection pass after this table was written.
Current numbers of record are in REPORT.md (SNS 6,794 / SADb 2,063 narrative;
280 pages). Round-2 additions since the round-2 tables: the SADb prune paragraph
now states the verified pre-prune base size of 1,393 records (957 post-prune +
436 removed; papers 561 + 382 = 943 audited; consistent with deletion_manifest.json
88+290+4 papers and 3+1+3+47 orphaned links) per Ben's recollection ("like 1300").
