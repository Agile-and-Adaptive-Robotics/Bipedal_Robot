# Ordered audit follow-up: Chapter 2, Ben's findings, and the other audit

30 October 2026. **Audit and proposed-change ledger only. No manuscript, figure, model input, measurement, optimization result, live comment, or tracked change was edited. No commit, push, or Overleaf upload was made.**

The three requested tasks were completed in order. Task 1 closed before Ben's findings were incorporated. Task 2 closed before the named other audit was opened. The comparison then expanded to include that file's later Round 3/4 additions, with both reviewed versions preserved. Subagents independently checked science, figures/layout, mechanics, the paper/supplement, and comment/reply context; the lead reconciled their findings.

## What this follow-up adds

| Requested task | Completed result | Navigation |
|---|---|---|
|1. Expand to Chapter 2,15-background |13 findings: six scientific/argument issues and seven layout/prose issues. All 290 source lines and all 32 printed pages 9-40 reviewed. | Chapter 2 audit and evidence sections below. |
|2. Incorporate Ben's findings | All 32 finding groups recorded, with 32a/32b separately tracked; all 42 saved advisor cards have current-state dispositions. One 159-word abstract proposal fits one page in the current dissertation formatting. | Ben's summary ledger, detailed paper/mechanics/SADb ledgers,42-card crosswalk, and abstract proposal. |
|3. Evaluate every item in the other audit | Initial report and Round 2 assessed; all 33 Round 3 B-items and all 39 Round 4 table entries assessed, plus additional claims/session notes. Accepted, qualified, rejected, and source-verification decisions are explicit. | Comparison sections for method/appendices, Chapters4-7, paper, Chapter 2, and the later rounds. |

Chapter 2's most important corrections concern the categorical exclusion of dielectric-elastomer legged robots, an electrically driven humanoid cited as a BPA example, the Geyer/Markin reflex-versus-CPG account, and the blanket healthy-human stimulation prohibition. The primary sources support narrower, accurate statements. The chapter already has substantial visual material; its main presentation needs are readable circuit/torque figures, better mechanical-to-neural float placement, and less repeated architecture/platform explanation.

Ben's findings expose specific misses beyond the preceding audit. Table C13 is an exact duplicate of the main biomimetic-fit table; C.9.1 repeats metric definitions. C.10 needs the master's-model provenance. Figure 4.10's northeast legend is still pending. Figure 4.14 is already correct, while the paper's same-named flexor asset differs and still has the wrong stacking/opacity/edges. The monthly watch is built and has saved output, but an unattended schedule is not established. The advisor's request for a relevant Table-tab screenshot in Chapter 5 remains unfulfilled.

The other audit contributes useful specific callout checks: six Appendix D figures lack explicit references, two Appendix E tables need purposeful references, and five supplementary displays lack direct main-text references. Its first-reference-order and rho-symbol checks also add value. Several proposed remedies need correction: the Chapter 5 future-work instruction was mapped to the paper; Appendix C table identities were wrong; current knee-route costs were attributed to Morrow; a universal figure-callout-first rule contradicts result-first guidance; some double-caption/panel-reversal claims are unsupported. Its initial15-display Frontiers blocker was later retracted, and that correction is recognized here.

The other audit's proposed abstract is206 whitespace-separated words, despite its 235-word label, and occupies **two pages** when compiled in the current class. Our 159-word proposal occupies **one page** without changing margins, font, spacing, or heading. Both are proposals, not active manuscript edits.

## Authority, locations, and verification limits

The current dissertation/paper ZIPs are the source authority. They remain unchanged. Isolated extracts are under `Dissertation/Archive/audit_20261030/dissertation` and `.../paper`; fresh builds are306 dissertation pages, 33 paper pages, and 10 supplement pages. The ZIP's bundled paper PDFs are stale 38/five-page artifacts; this is distinct from its live preview, which Ben reported recompiling at 33 pages.

All line anchors are one-based snapshot source lines. Chapter numbers and figure/table numbers follow the fresh build. A section range is a navigation aid, not a request to delete the whole interval. Artifact/package-level observations have no meaningful manuscript line; their source pointers or metadata are given instead.

| Short source key | Chapter / actual file |
|---|---|
|15 | Chapter 2,chapters/15-background.tex |
|30 | Chapter 4,chapters/30-torque.tex |
|50 | Chapter 5,chapters/50-database.tex |
|60 | Chapter 6,chapters/60-sns.tex |
|70 | Chapter 7,chapters/70-overall.tex |
|92/93/94/95/96 | Appendices A/B/C/D/E,chapters/92-AppendixA.tex through 96-AppendixE.tex |
|P/Main | Torque paper,Main.tex |
|S/Supplement | Torque paper,frontiers_SupplementaryMaterial.tex |

The named other audit was first frozen with SHA256 `859e9a7fc0738d92116287601c5157716ea3e8fc0e7cd7351ee484974185f772`. It changed during review; the later Round 3/4 version was separately frozen with SHA256 `0ec8ef1289fbe74504ca2bc8952ab20f66a4a3e31f0284d3114bae95414b9b13`. Earlier comparison sections identify the first version; later-round sections incorporate its added proposals and corrections. Neither original report was overwritten.

The saved 42-card review is from 29 October. It cannot certify every comment added live on 30 October. Cards4.1 and 5.25 retain uncertain live anchors, and 5.23's Table screenshot request remains unmet in the source. The UI accessibility state and screenshot identified different foreground windows during attempted live checking; no review action was taken from that unreliable state. No blanket comment resolution or tracked-change acceptance is claimed.

The FVU item remains handled by Ben's reply. The offset ruling, deferred classification work, excluded AI-tool assessment/comparison, rejected uncertainty exercise, waived individual pulley example, protected published Chapter 3, and frozen biomimetic results remain respected. Existing-record provenance checks are proposed; no fitting/optimization reruns are assigned.

## How to read the comparison

“Agree/accept” means a proposed change has merit. “Modify/qualify” keeps the useful concern while correcting its scope or remedy. “Reject” identifies an unsupported assertion or a conflict with Ben's ruling. “Verify” means an existing source/result/permission record must settle the point before wording is finalized. None means a manuscript change has been made.

The numerical manuscript scores are treated as subjective screening indices, as the other reviewer itself states. They are not competence grades or submission clearance. The assistant-failure rubric remains a separate feedback record: it penalizes missed instructions/replies, unobserved edits, unsupported claims, and unclear exposition. The comparison does not claim that manually assigned penalties are a trained predictive model.

The source-location index at the end covers all 37 original audit items. Chapter 2's13 items and every Ben/comparison item carry chapter/file/line evidence in their own sections. Historical findings remain available in the original23-page audit; this document records the additions and the explicit superseding decisions.


---

# Chapter 2 audit completed - Task 1

30 October 2026. Authority: the latest supplied dissertation ZIP, unchanged since the preceding audit, extracted under Archive/audit_20261030/dissertation. Source references are one-based lines of chapters/15-background.tex. The source was read in full, all 32 printed pages 9-40 (PDF 30-61) were visually reviewed, and targeted scientific claims were checked against the primary sources identified in CHAPTER2_SCIENCE_NOTES.md. No manuscript was changed.

This completes Task 1 before incorporation of Ben's new findings or comparison with the other audit. The thirteen findings below combine the independent science and layout reviews. They are new audit advice, not invented advisor directives. No FVU, AI-comparison, deferred-classification, or new-optimization task is introduced.

| ID | Chapter/section | Source lines | Finding and proposed response |
|---|---|---|---|
| C2-S1 | Ch. 2 §2.2.1 | 30 | DEA limitations did not keep them out of all legged robots: the primary DEAnsect example contradicts this. Scope the limitations to the human-scale load-bearing application. |
| C2-S2 | Ch. 2 §2.2.2; §2.1 | 61;8 | Asano's electrically actuated humanoid is cited as a BPA example despite the correct electrical description earlier. Remove the mismatched BPA citation, preserving the musculoskeletal example. |
| C2-S3 | Ch. 2 §2.9 | 244 | Geyer is not three positive feedback channels around every muscle; Markin's RG/PF model does not belong in the CPG-free class. Correct the pathway-specific description and taxonomy. |
| C2-S4 | Ch. 2 §2.10; §2.6 | 270;173 | The healthy-human prohibition groups destructive interventions with stimulation too broadly, despite the human-stimulation example already cited. Narrow the claim to destructive/highly invasive interventions. |
| C2-S5 | Ch. 2 §§2.5/2.9 -> §2.7 | 166,244 ->230-232 | Two references promise circuit/model content that the generic Computational Neuroscience paragraph does not provide. Point to the actual circuit sections/figures. |
| C2-S6 | Ch. 2 §2.2.3 | 65,95 | 'Even more accurate' compares different experimental domains, and 'none' overextends the named model survey. Retain reported fits, scope the comparison, and state the specific multi-length design gap. P3 conditional argument advice. |
| B2-L01 | Ch. 2 §§2.4-2.5 | 143-164 | Figures2.6/2.7 interrupt the first CPG paragraph across printed 22-25. Finish the mechanical float sequence before the neural section. |
| B2-L02 | Ch. 2 §2.4 | 143,154-158 | Figure2.7's six-panel comparison at 0.62 textwidth is too small for quantitative inspection. Enlarge or select representative panels while preserving master's attribution. |
| B2-L03 | Ch. 2 §2.6 | 175,190-203,207,211-226 | Figures2.9/2.10/2.12 contain small neural labels; height limits dominate the nominal width in 2.12. Improve final-size allocation or select fewer clearly readable examples. |
| B2-L04 | Ch. 2 §§2.5,2.6,2.8,2.9 | 164-166,175,207,236-248 | RG/PF explanations and platform inventories repeat; printed 31-36 are six consecutive text pages. Explain architecture once and use later paragraphs for specific differences and purpose. |
| B2-L05 | Ch. 2 §2.10 | 270-290, especially284/286 | Film-procedure chronology, the horse comparison, and ledger metaphor are peripheral to the engineering replacement argument. Compress while preserving the ethical stance and validation limits. |
| B2-L06 | Ch. 2 §2.2.1 | 34-36 | Duplicate Figure2.2 introductions and schematic-production narration. Keep one purpose-led callout. P3. |
| B2-L07 | Ch. 2 introduction; §2.8 | 4,240 | Repeated boundary self-narration and third-person references to Ben's program. Keep one concise roadmap and the sanctioned first-person voice. P3. |

The chapter already has 13 figures, one coefficient table, and a force equation. Its weakness is not a lack of illustrations. The most valuable changes are scientific precision, final-size readability, float placement, and removal of duplicate explanation. The primary mechanical-to-neural progression, human-versus-robot knee distinction, master's lineage, circuit source credits, and documentary-image disclosure are strengths to preserve.

Complete evidence, primary-source URLs, confidence, printed/PDF-page references, and retained strengths are in CHAPTER2_SCIENCE_NOTES.md and CHAPTER2_LAYOUT_NOTES.md. CHAPTER2_LEAD_NOTES.md records the lead's independent editorial read; its overlapping findings are not added again to the thirteen-item count.


---

# Chapter 2 scientific and argument audit notes

30 October 2026. **Task 1 only.** The source is the latest Overleaf snapshot `Archive/audit_20261030/dissertation/chapters/15-background.tex`. All references below identify Chapter 2, its compiled section number, and one-based source lines. No manuscript, bibliography, figure, or result was changed. Incorporation of Ben's new findings and comparison with the other audit were not started by this reviewer.

Read all 290 lines of the background source and checked the compiled contents/section mapping. This is a scientific/argument and targeted primary-source check, not visual certification of every page or a full literature systematic review. The parent handles layout and overall synthesis. The dissertation skill and relevant inherited AGENTS context were applied. No FVU, optimization, AI-tool-comparison, or deferred classification task is reopened.

## Ranked findings

### C2-S1 — P2: the categorical exclusion of dielectric elastomers from legged robots is false

**Location:** Chapter 2, Section 2.2.1 Overview, `15-background.tex:30`.

**Chapter wording:** kilovolt driving voltages and electrode/instability constraints “have kept them out of legged robots.”

**Evidence:** DEAnsect, published in *Science Robotics* in December 2019, has three legs actuated by low-voltage dielectric elastomers. The authors' EPFL laboratory describes those legs, autonomous operation, and the reduced-voltage electronics explicitly. This directly contradicts the broad exclusion, independently of whether such actuators are suitable for Ben's human-scale load-bearing leg.

**Minimal proposed change:** retain the engineering limitations but scope them to integration in the human-scale system studied here; do not claim that no legged DEA robots exist. One brief counterexample/citation is sufficient; no extended new technology review is needed.

**Primary source:** [EPFL's DEAnsect project](https://www.epfl.ch/labs/lmts/lmts-research/deansect/), linking Ji et al. (2019), DOI [10.1126/scirobotics.aaz6451](https://doi.org/10.1126/scirobotics.aaz6451). **Origin:** new audit finding. **Confidence:** high.

### C2-S2 — P2: an electrically actuated humanoid is cited as a BPA example

**Location:** Chapter 2, Section 2.2.2 Braided Pneumatic Actuators, `15-background.tex:61`; conflicting description in Section 2.1, `:8`.

**Chapter wording:** the BPA applications sentence cites `asano_design_2017` among humanoid examples. The opening robot paragraph correctly describes the same work's brushless motors, wire tendons, gearing, and tension measurement.

**Evidence:** the two passages assign different actuation technologies to the same cited example. The cited Asano work concerns the tendon-driven Kenshiro/Kengoro humanoids, not a demonstration that those humanoids use BPAs.

**Minimal proposed change:** remove that citation from the BPA-specific application claim, retaining it in the musculoskeletal/tendon-driven humanoid account. Do not replace it with an unverified citation. The other cited humanoid example can remain subject to its own source support.

**Source check:** [Asano et al., publisher DOI record](https://doi.org/10.1126/scirobotics.aaq0899), with bibliographic/abstract confirmation at [PubMed](https://pubmed.ncbi.nlm.nih.gov/33157878/). Publisher full-text access returned 403, so the strongest evidence available here is the chapter's own explicit electrical-actuation description and the identified primary work. **Origin:** new citation-consistency finding. **Confidence:** high internal mismatch; not a re-audit of every robot detail in the source.

### C2-S3 — P2: the reflex-only versus CPG taxonomy misdescribes the cited mechanisms and lineage

**Location:** Chapter 2, Section 2.9 Neuromechanical Simulation, `15-background.tex:244`.

**Chapter wording:** Geyer and Herr's model is described as using “positive force, length, and velocity feedback around each muscle”; the next sentence places Markin/Shevtsova cat models in the following generation of that reflex-based line.

**Evidence:** Geyer and Herr's actual Figure 1 and Sections II-A–C use different pathways for different muscles and gait phases: positive force/length feedback, negative force/length feedback, contact-triggered phase changes, and trunk proportional-derivative control. The model is correctly described as CPG-free, but not as giving every muscle the same three positive feedback channels. Markin et al. (2010), in contrast, explicitly contains a two-level CPG with RG/PF populations and afferent feedback. It should not be introduced as though it belongs to the CPG-free class just described or as an established causal descendant of Geyer.

**Minimal proposed change:** describe Geyer as muscle-/phase-specific reflex control without a CPG. Place the Markin/Shevtsova examples in the CPG-plus-feedback branch, or use them to explain overlap rather than a one-way lineage. No additional simulation or literature model replication is required.

**Primary sources:** [Geyer and Herr's published article PDF, hosted at CMU](https://www.cs.cmu.edu/~cga/tmp-public/song.pdf), Figure 1 and Sections II-A–C, DOI 10.1109/TNSRE.2010.2047592; [Markin et al. (2010), authors' hosted PDF](https://rybak-et-al.net/pdfs/Markin_et_al_2010.pdf), Figure 1, DOI 10.1111/j.1749-6632.2010.05435.x. **Origin:** new literature-description finding. **Confidence:** high for the pathway and CPG-class distinction; historical influence beyond these mechanisms was not reconstructed.

### C2-S4 — P2: the ethical argument groups destructive lesions with stimulation too categorically

**Location:** Chapter 2, Section 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case, `15-background.tex:270`; related human-stimulation example in Section 2.6, `:173`.

**Chapter wording:** the list includes muscle deletion, sensory-pathway blocking, RG lesioning, and tonic spinal stimulation, then states that these experiments cannot be performed on healthy humans.

**Evidence:** the destructive interventions and non-invasive stimulation are materially different ethical/methodological categories. The chapter itself already cites a human tonic-stimulation air-stepping experiment. Primary research also explicitly reports transcutaneous spinal stimulation in healthy adults. Those studies do not justify destructive lesion experiments, but they make the combined prohibition overbroad.

**Minimal proposed change:** limit the categorical claim to destructive or highly invasive interventions and explain that robotic models additionally offer repeatable, selective pathway interventions. Preserve the ethical case and Three Rs context; do not add a generic ethics-warning checklist or imply that robots replace biological validation.

**Primary source check:** [Shamantseva et al. (2024), healthy-adult transcutaneous spinal-stimulation experiment](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2024.1467182/full), DOI 10.3389/fnins.2024.1467182, Methods and introduction's non-invasive-study distinction. **Origin:** new scope/counterexample finding. **Confidence:** high; no legal or institutional-review ruling is inferred.

### C2-S5 — P2: two semantic cross-references point to a section that does not supply the promised content

**Location:** Chapter 2, Section 2.5 at `15-background.tex:166` and Section 2.9 at `:244`; target Section 2.7 Computational Neuroscience, `:230-232`.

**Evidence:** line 166 directs readers to Section `sec:bg_compneuro` for expression of the architectures as simulable circuit models; line 244 points there for cat hindlimb examples. Section 2.7 is a single generic paragraph defining computational neuroscience and contains no model-specific circuit account. The actual circuit descriptions and figures are in Sections 2.5–2.6. These references compile but do not help the reader find the promised evidence.

**Minimal proposed change:** point to the actual model-specific sections/figures or remove the redundant pointer. A substantive model discussion need not be duplicated in the short computational-neuroscience definition just to satisfy the reference.

**Origin:** new internal-reference finding. **Confidence:** high, source verified.

### C2-S6 — P3: cross-study accuracy language exceeds the comparison presented

**Location:** Chapter 2, Section 2.2.3 BPA Force–Length–Pressure Models, `15-background.tex:65,95`.

**Chapter wording:** the paragraph moves from an “impressive” R-squared for one actuator to an “even more accurate” model on different actuator diameters/resting lengths, then says none of the existing models accounts effectively for maximum-force differences with initial length.

**Assessment:** the quoted fit values can be retained, but no common test domain, normalization, or held-out comparison is supplied in this paragraph to support an overall model ranking across distinct actuator tests. The universal “none” also reaches beyond the named models. This audit has not established that the reported numerical fits or the dissertation's multi-length novelty are false.

**Minimal proposed change:** describe the literature as reporting high within-dataset fits or low reported errors, then state the actual unresolved design need: predicting force across resting lengths in the studied Festo actuator family. Scope the model limitation to the models evaluated here, rather than all existing models. Retain the specific numerical evidence and Chapter 3 comparison pointer; no new comparative experiment is requested.

**Origin:** new argument/register recommendation, not an advisor directive or independently established numerical error. **Confidence:** high that the paragraph supplies no common comparison; broader literature-exhaustiveness remains unverified.

## Strengths and boundaries

- Sections 2.3–2.4 clearly distinguish human Gait2392 translations from the crossed four-bar robot knee and state that different kinematics give different moment arms.
- The master's-work paragraph at line 123 identifies Ben's previous project directly, and the placement-optimization handoff is relevant.
- Section 2.3's signed fixed-axis moment arm is carefully defined; tendon-excursion geometry is not silently replaced with a simple fixed-axis expression for constrained/moving paths.
- Sections 2.5–2.6 include the relevant RG/PF, sensory, commissural, and biarticular ideas and actual circuit figures. The engineering motivation is present.
- The background/contribution distinction and Chapter 5/6 handoffs should be retained, though repeated statements of that boundary can be shortened by the parent flow audit.
- The source already distinguishes simulation, mechanism tests, and biological inference. These strengths should not be replaced by broader claims of universal physiological accuracy.

All six items are new audit findings/recommendations. No claim is made that the advisor specifically requested them, that his present live comments were fully recaptured, or that all references were independently validated. No item requires reopening Ben's settled FVU, AI-comparison, or optimization decisions.



---

# Chapter 2 prose, layout, and figure audit

30 October 2026. Task 1 only: expansion of the audit to Chapter 2. No manuscript edits, incorporation of the later user-findings list, or comparison with the other audit was performed by this reviewer.

## Authority and coverage

Read all 290 lines of `Archive/audit_20261030/dissertation/chapters/15-background.tex`. Inspected all 32 rendered Chapter 2 pages: printed pages 9-40, corresponding to PDF pages 30-61 of the fresh 306-page snapshot build. QA renders are in `Archive/audit_20261030/chapter2_layout_qa`. Source references below use this ZIP snapshot, so they remain stable even if the working manuscript changes later.

This is a prose, relevance, and visual audit. Literature-supported factual challenges belong in the separate scientific audit; no external factual error is inferred merely from a diagram's appearance. Recommendations are new audit advice unless explicitly identified otherwise. P2 means substantive reader/presentation repair; P3 means optional polish. No P1 scientific error was established by this layout review.

## Ranked findings

### B2-L01 — P2: mechanical figures interrupt the first CPG paragraph

**Location:** Chapter 2, Sections 2.4-2.5; `15-background.tex:143-164`, Figures 2.6/2.7; printed pages 21-25 / PDF pages 42-46.

**Evidence:** The first paragraph of Section 2.5 begins on printed p.22. Figure 2.6, which belongs to the preceding moment-arm discussion, then occupies the top of p.23; Figure 2.7 takes the entire next page. The CPG paragraph resumes beneath Figure 2.6 and continues on p.25 after another full-page interruption. The reader must carry an explanation of intrinsic rhythm generation across two unrelated mechanical figures.

**Minimal remedy:** finish the mechanical display sequence before Section 2.5, for example with a carefully placed float barrier or revised placement/size. Keep the CPG paragraph together as far as practical. This does not require adding content or deleting the mechanical evidence. **Confidence:** high, verified on all affected rendered pages.

### B2-L02 — P2: Figure 2.7's six-panel torque comparison is too small for its quantitative role

**Location:** Chapter 2, Section 2.4; `15-background.tex:143,154-158`; Figure 2.7, printed p.24 / PDF p.45.

**Evidence:** The complete six-panel figure is included at only `0.62\textwidth`. The surface plots' axis ticks, titles, and color-bar numbers are much smaller than the caption and surrounding body text. Large external margins remain on the float page. The paragraph asks readers to compare knee and hip torque surfaces, but the reduction makes the quantitative details difficult to inspect.

**Minimal remedy:** use more of the available page width and enlarge meaningful labels at final size. If the purpose is simply to establish a torque-profile mismatch, retain the clearest representative comparison rather than shrinking every legacy panel. Preserve the master's-work attribution and actual plotted data. **Confidence:** high visual finding; no unmeasured point-size value is asserted.

### B2-L03 — P2: several cited circuit schematics are reduced beyond comfortable label inspection

**Location:** Chapter 2, Section 2.6; `15-background.tex:175,190-203,207,211-226`; Figures 2.9, 2.10, and 2.12; printed pp.28-30 / PDF pp.49-51.

**Evidence:** Figure 2.10's bilateral network, seven-muscle hindlimb, feedback labels, and commissural cell names are all compressed into 62% of the text width. Figure 2.12 is a tall diagram constrained to 38% of text height: its nominal 88%-width setting does not make the actual image wide, because the height limit dominates. Figure 2.9 similarly presents many small interneuron/receptor labels at 62% width. The reader is asked to understand named neural populations and their routing, not merely recognize that a network exists.

**Minimal remedy:** give the most relevant foundational/implemented networks enough area to read; use a separate page or focused crop if needed. Keep the source citation and reproduction credits. Consider whether all five circuit reproductions need full display when two clearly selected examples plus the existing citations could establish the distinctions. This is a selection/scale recommendation, not a demand to remove cited neuroscience. **Confidence:** high final-page legibility issue.

### B2-L04 — P2: the two-level architecture is repeatedly reintroduced instead of progressively developed

**Location:** Chapter 2, Sections 2.5, 2.6, 2.8, and 2.9; `15-background.tex:164-166,175,207,236-248`; printed pp.22-27 and 31-34 / PDF pp.43-48 and 52-55.

**Evidence:** The RG/PF split, deletions, biarticular convergence, and afferent coupling are explained in 2.5; retold in the three-model circuit walkthrough at :175; developed again as the Rybak/Shevtsova lineage at :207; and summarized again through a platform/citation list in 2.8 and the three simulation traditions in 2.9. Each paragraph is locally intelligible, but their roles overlap. The printed pp.31-36 are six consecutive text-only pages despite the preceding concentration of circuit figures.

**Minimal remedy:** define the architecture once, use the circuit comparison to explain the specific differences that matter to this dissertation, and keep 2.8/2.9 focused on why SNS and the chosen software enable the subsequent experiments. Shorten duplicate model inventories and generic statements about comparing executable mechanisms. A comparison table is an option only if it replaces the repeated paragraphs; another figure is not the default cure. **Confidence:** high prose/flow finding.

### B2-L05 — P2: the ethics section has a long procedural film narrative that is peripheral to its engineering argument

**Location:** Chapter 2, Section 2.10; `15-background.tex:270-290`, especially :284 and :286; printed pp.36-40 / PDF pp.57-61.

**Evidence:** The section properly states replacement, repeatable interventions, and biological-validation limits. It then devotes a long paragraph to the sequence of a classroom film, including animal handling, anesthesia delivery, a comparison with the horse demonstration, carotid ligation, cannulation, trephination, and the transection level. The following paragraph uses phrases such as “honest accounting is murkier,” “cannot repay,” “knowledge purchased at that cost,” and “adding to the ledger.” This makes the passage longer and more rhetorical than the scientific replacement argument requires.

**Minimal remedy:** preserve the ethical stance, the fact that present knowledge relies on animal preparations, and the need for biologically validated robotic models. Compress the procedural chronology and horse comparison to the minimum needed to explain the preparation; cite the film/history rather than narrating it. Prefer a direct statement about reusing existing knowledge to the ledger metaphor. This is compression advice, not a request to soften Ben's position or hide the preparation. **Confidence:** high relevance/register finding.

### B2-L06 — P3: Figure 2.2 receives duplicate introduction and an editorial-style description

**Location:** Chapter 2, Section 2.2.1; `15-background.tex:34-36`; printed p.13 / PDF p.34.

**Evidence:** The adjacent paragraphs both begin “Figure ...” and say it compares/shows representative actuator mechanisms. The second repeats the panel descriptions already carried by the graphic/caption and closes with “The panels use a consistent schematic style,” which describes illustration production rather than the engineering comparison.

**Minimal remedy:** retain one purposeful callout explaining the mechanism/design contrast and remove the duplicate lead and schematic-style narration. No extra display is needed. **Confidence:** high exact-text duplication.

### B2-L07 — P3: chapter-boundary self-narration repeats and conflicts with the preferred first-person voice

**Location:** Chapter 2 introduction and Section 2.8; `15-background.tex:4,240`; printed pp.9 and 32 / PDF pp.30 and 53.

**Evidence:** The introduction explains that the author's program belongs in Chapters 5/6. Section 2.8 repeats the boundary explanation and ends “None of it is repeated as background here.” Both refer to “the author's” work/program despite Ben's standing preference for first person. The opening already supplies the needed navigation.

**Minimal remedy:** retain one concise chapter map, remove the later administrative boundary commentary, and use first person where retained text refers to Ben's own work. **Confidence:** high; the first-person direction is a standing author instruction, while the proposed cut is new editorial advice.

## Positive findings and limits

- The chapter provides a coherent mechanical-to-neuroscience progression and explicitly names the human reference and physical biomimetic knee.
- The master's project is identified at :123 and its comparison is tied to the later route-design work. That connection should survive compression.
- The introductory force-model equation and coefficient table are legible and fit their pages. They give a useful technical contrast with the later force characterization.
- Figures 2.1, 2.3, 2.4, and 2.5 have useful visual roles and are large enough for their qualitative purpose. The circuit reproduction captions carry source/figure credits.
- Figure 2.13 clearly discloses restoration/machine colorization and states that it is documentary context, not gait evidence. No hidden use of a generated image as measurement evidence is alleged.
- No visibly clipped or off-page Chapter 2 display was established. The issue is allocation and sequencing of space, not the mere existence of white space or a required figure-per-page quota.

No manuscript text, figure, model, or number was changed. These notes complete this reviewer's Chapter 2 layout/prose portion of Task 1.


---

# Ben's findings incorporated - Task 2 completed

Scope is **audit and proposed-change ledger only**, as Ben specified. No manuscript, figure, application, result, optimization, live comment, or tracked change was modified. Each finding below has a current-source response and verification boundary in the named detailed ledger. Source line numbers refer to the supplied ZIP snapshots; section numbers describe the current build and may change after eventual edits.

Detailed records: USER_FINDINGS_PAPER.md (**P**); USER_FINDINGS_DISS_MECHANICS.md (**M**); USER_FINDINGS_SADB_AND_COMMENT_CROSSWALK.md (**D**); ABSTRACT_PROPOSAL.md (**A**). The last record contains all 42 saved advisor cards, with Ben's replies/rulings and present-source dispositions. These are separate from Ben's new findings.

| ID | Chapter / source anchor | Disposition / proposed change | Detail |
|---|---|---|---|
| U01 | Paper §3.1 Main.tex:414,418; supplement:47-66 | Prefer removal/consolidation of Table S1; preserve any essential common-domain comparison. It is not a literal duplicate of main Table 1. Supersedes the earlier blanket Keep advice. | P |
| U02 | Paper §4.1 Main.tex:546; supplement:186-192 | Remove the slippage photograph from the journal supplement following Ben's unused-42-cm-test provenance; audit its callout and dependent claims. Supersedes earlier Keep-if-linked advice. | P |
| U03 | Ch. 4 Table 4.8,30-torque.tex:617-632; App. C:433-450 | Remove duplicate Table C13, retaining one authoritative main-text table. | M |
| U04 | App. C §C.9.1,94-AppendixC.tex:429-430; Ch. 4:316-336 | Remove repeated metric definitions and empty C.9 container after U03; no FVU revision. | M |
| U05 | App. C §C.10/Table C14,94-AppendixC.tex:452-495; Ch. 4:792 | Plain-language purpose and explicit master's-model attribution using the existing Bolen 2019 citation. | M |
| U06 | Paper §1 Main.tex:99-101; §4.4:637-656 | Replace anthropomorphic segue; include the existing improved 27-site model comparison with accurate scope. | P |
| U07 | Dissertation abstract,02-abstract.tex:1-7 | Proposed 159-word replacement compiled and visually checked at one page with current dissertation formatting. Not installed. | A |
| U08 | Paper §2.9 Main.tex:397-404; Ch. 4 figure:780-782 | Replace the old schematic with the updated crossed-link two-BPA mechanism, coordinated with U17/U26; move under U09/U14. | P |
| U09 | Paper §2.9:397-404; §§4.3.2/4.4:604-660 | Put the substantive pulley explanation after the meet/exceed hardware evidence. | P |
| U10 | Paper §§1/2.4/4.5,Main.tex:97,166,666 | Already satisfied: Morrow's cost is not reproduced. Preserve concise attribution. | P |
| U11 | Paper Figure 11,Main.tex:410,499-505 | Remove its hybrid series, legend entry, caption phrase, and linked discussion; retain unrelated pinned-knee evidence. | P |
| U12 | Paper §4.3.1 Main.tex:581-602 | Shorten repeated setup and numeric inventories; use a compact existing-result summary, with unique geometry in the supplement. | P |
| U13 | Paper §4.3.2 Main.tex:604-627 | Lead with all-nine/all-fifteen measured capacity outcomes, then margins, figure, and essential test details. | P |
| U14 | Paper §4.4 Main.tex:629-643 | Canonical home for the new pulley schematic and concise force-travel explanation. One move shared with U08/U09. | P |
| U15 | Paper Figure 13,Main.tex:606-623 | Use cam-assisted extensor route; match panel font/line scales; prediction before measured; scatter alpha 0.75, no edges. | P |
| U16 | Paper §4.4 Main.tex:631-643; supplement:255-296 | Recover the existing transmission-selection objective/ranking rule; give a concise description and exact supplementary expression/table as needed. Do not invent a method or rerun optimization. | P |
| U17 | Ch. 4 §4.3.5,30-torque.tex:763-782 | Use F_R for retained reaction labels; coordinate with deletion of mount-reaction treatment. | M |
| U18 | App. C Tables C5/C6,94-AppendixC.tex:105-130 | Convert dimension values/headings to meters; preserve pressure/force and nominal product names. | M |
| U19 | App. C §C.4,94-AppendixC.tex:205-214; Ch. 4:140-250 | Remove repeated compliance display; retain unique arrays, bounds, frames, and implementation mapping. | M |
| U20 | Ch. 4 Table 4.3,30-torque.tex:417-434 | Stack FLEX and EXT weights on separate display lines without changing values. | M |
| U21 | Ch. 4 §4.1.5,30-torque.tex:374,451-466 | Give constrained redesign a numbered subsection and shorten its tail; avoid a global style change to protected Chapter 3. | M |
| U22 | Ch. 4 Figure 4.10,30-torque.tex:594-599 | Still pending: legend is below panels. Move northeast inside B with consistent model-first order. | M |
| U23 | Ch. 4 Table 4.10,30-torque.tex:669-685 | Percent units/no plus sign; identify achieved versus required margins and retrieve the flexor worst-case angle from saved results. No guessed angle or normalization. | M |
| U24 | Ch. 4 §§4.3.3/4.3.4,30-torque.tex:715-733 | Explain manifold lower operating limit; distinguish Ben's approximately 200-kPa experience from the saved 150-kPa model assumption. | M |
| U25 | Ch. 4 Figure 4.14,30-torque.tex:742-747 | Verified already correct: predictions drawn before measured; alpha approximately 0.75 and no strokes. Preserve; legend order is separate. | M |
| U26 | Ch. 4 §4.3.5,30-torque.tex:753-788 | Remove mount-reaction derivation/example/discussion; use approximately two paragraphs and RoM. Supersedes the earlier favorable audit note. | M |
| U27 | Ch. 4 §4.3.6,30-torque.tex:792; App. C:454,458,500 | Briefly define the evaluated joint-coordinate grid/domain and coverage; replace unexplained graded terminology. | M |
| U28 | Ch. 4 §4.3.6:792; App. C Table C14:458-489 | State ideal diameter-squared maximum-force rationale for four 20-mm BPAs versus one 40-mm BPA, with matched-condition qualification. | M |
| U29 | Ch. 5 Figure 5.2,50-database.tex:39-44 | Crop the excessive bottom whitespace and inspect final-size readability. | D |
| U30 | Ch. 5 §5.3.2:61-67; Ch. 7 §7.2.3,70-overall.tex:59-61 | Move future improvements to Chapter 7, retaining the useful Chapter 5 segue and present coverage limits. | D |
| U31 | Ch. 5 Methods:14-16; Ch. 7:61; SADb_audit/watch | Describe the monthly watch as implemented; retain manual curation and distinguish the provided scheduling command from an installed schedule. | D |
| U32a | Ch. 5:8-10,24-44; App. E:16-21 | Bring a current relevant Table-tab example forward; preserve six bubble-map callouts and place the cropped pivot near its use. | D |
| U32b | All 42 saved cards,29 October ledger; Ch. 4/5 anchors in D | All saved cards dispositioned. Screenshot request5.23 remains unmet;4.1/5.25 anchors uncertain. New30 October live-comment completeness cannot be certified from saved evidence. | D |

## Reconciliation with the preceding audit

Ben's findings add specific corrections that the earlier audit missed, particularly the exact C13/C.9.1 duplications, C.10 master's provenance, paper carry-over changes, operating-pressure qualification, and the unfulfilled Table screenshot request. The earlier positive recommendation for detailed mounting reactions and retention of S1/S2 is superseded as stated above. Conversely, U10 and U25 are already satisfied and should not create unnecessary edits.

The 42-card crosswalk preserves Ben's no-change offset ruling, his completed FVU reply, deferred SADb classification work, and the distinction between a source change and a resolved live review thread. No AI-tool assessment/comparison task is added. The UI accessibility state and screenshot did not identify the same foreground window during the live check; no review action was taken on that unreliable state. This is why present-day live completeness remains unverified.

Task 2 is complete as an audit/ledger task: all 32 finding groups (with 32a/32b separately recorded), every saved advisor card, and the requested one-page abstract proposal have a response. No implementation is claimed. Task 3, comparison with the named other audit, begins only after this completion checkpoint.


---

# U07 - one-page dissertation abstract proposal

**Anchor:** dissertation front matter, chapters/02-abstract.tex:1-7. This is a proposal, not an applied manuscript edit. The active abstract is in the supplied ZIP; commented historical drafts below line7 are not active text.

The following 159-word proposal was compiled independently with the supplied main.tex preamble and psutex.cls at the existing font, spacing, margins, chapter heading, and roman page numbering. It occupies exactly one page with no compile error. Proof and render are in Archive/audit_20261030/abstract_proposal/abstract_proof.pdf and abstract_proof.png. No active source file, page geometry, or font was changed to make it fit.

## Proposed text

Biomimetic robots provide physical models for investigating neuromechanical control. This dissertation develops tools to design and simulate a human-proportioned lower-body robot actuated by braided pneumatic actuators (BPAs).

I characterized BPA force across resting lengths, contractions, and pressures, producing a resting-length-dependent maximum-force model and a three-coefficient nondimensional force surface. Knee experiments identified prediction errors from length offsets, effective compliance, and wrapping-induced loss of usable length. Incorporating these effects supported route redesign. Paired flexor and extensor BPAs met or exceeded selected human-muscle reference torques at all 15 flexor and nine extensor measurement angles, establishing discrete static capacity.

SADb Explorer organizes 1,393 literature and related records for investigating sensory feedback. I also developed neural-model tools and interfaces for musculoskeletal simulation in MuJoCo and Simulink/Simscape. Matched-input tests, stimulation, pathway removal, and independent replays support controller verification while separating neural effects from coupling errors. Supported evaluations demonstrate selected controller functions without establishing harness-free walking. These tools support first-pass mechanical design and repeatable neuromechanical experiments.

## Evidence retained

The force-model contributions are those already reported in Chapter 3. The 15 flexor/nine extensor measurements are in Chapter 4 §4.3.4 (30-torque.tex:729-745); the 1,393-record resource is Chapter 5 (50-database.tex:8/20); the tools and support-qualified evaluations are Chapter 6. The proposal avoids adding new performance estimates or treating discrete static tests as a continuous dynamic envelope.


---

# Ben's paper findings: audit and proposed-change ledger

30 October 2026. Task 2 only. No manuscript, figure, calculation, result, or optimization was changed. The external audit comparison has not been started.

## Evidence and status conventions

These responses use the latest supplied ZIP source, not older ProofFinal line numbers:

- **P** = `Dissertation/Archive/audit_20261030/paper/Main.tex`.
- **S** = `Dissertation/Archive/audit_20261030/paper/frontiers_SupplementaryMaterial.tex`.
- Paper page numbers are the fresh 33-page `paper_build/Main.pdf`; supplement pages are the fresh 10-page `paper_build/frontiers_SupplementaryMaterial.pdf`. The generated PDFs bundled inside the ZIP are stale and are not the reference here (prior T04).
- The prior audit is `Notes/publication_audit_20261030/TORQUE_AUDIT.md`; its T IDs below identify overlap, not advisor instructions. Ben's current findings take precedence over its proposals.
- **Proposed** means recorded for review, not implemented. **Already satisfied** means the current source meets that particular request. No item below is presented as an independently recovered advisor comment.

The actual Figure 11 and Figure 13 pages, supplement overview, old pulley asset, and newer dissertation pulley asset were visually inspected. Exact plotting properties require source/export verification; a thumbnail or legend order alone does not prove MATLAB object order or alpha.

## Item-by-item responses

### U01 — Supplementary Table S1 seems redundant

**Anchor:** Paper §3.1, P:414 and P:418, p17; Supplement §1.1, S:47–66, Table S1, p1. Related main Table 1 is P:451–495, p21.

**Current state:** The table supports the same pre-correction comparison already described in the text and Figure 8, so Ben's compression concern has merit. It is not a literal numerical duplicate of main Table 1: S1 uniquely supplies Hunt-model and hybrid fit statistics, and its two model predictions use their common defined domain. Main Table 1 instead compares the improved BPA baseline with the corrected model on the full tests. Those domain differences explain why simply treating the two tables' numbers as interchangeable would be wrong.

**Proposed change:** Prefer removing S1 from the submission supplement if its unique numerical breakdown does not support an essential claim. Retain the concise force-model comparison in §3.1, with Figure 8 as its evidence, and preserve the common-domain qualification wherever a quantitative comparison remains. If an indispensable numeric comparison is retained, preserve only that number and its domain in the text rather than automatically keeping the whole table. Remove or redirect both S1 references at P:414/418, and update the supplement inventory at P:708. Preserve the original table in the archive. This supersedes the earlier unqualified “Keep” disposition; it does not declare all S1 numbers duplicated.

**Verification boundary:** Check the final text has no unsupported quantitative model ranking, dangling table reference, or silent substitution of full-domain statistics for common-domain statistics. No recalculation, FVU revision, or refitting is needed.

**Prior audit:** Supplement relevance table's S1 row; T15 for the inaccurate supplement inventory. Updated disposition: conditional removal/consolidation in response to Ben, not mandatory retention.

### U02 — Supplementary Figure S2 is an unused 42-cm pinned-extensor test

**Anchor:** Supplement §1.3, S:186–192, Figure S2, p5; paper §4.1, P:546, p23; inventory P:708, p30.

**Current state:** Current Figure S2 is a two-view **slippage photograph**, not a torque-versus-angle plot. Its caption identifies only pinned-extensor slippage. Ben identifies its source as the unused 42-cm test; the current caption does not independently establish that resting length or inclusion status. His firsthand test provenance is the basis for the proposed removal. This is distinct from Supplementary **Table** S2 and Supplementary **Figure** S1.

**Proposed change:** Remove this photograph from the journal supplement as Ben proposes, together with its explicit callout and “actuator slippage” inventory entry. Review P:546's photograph-specific narrative and its ±20-mm displacement claim: retain a concise qualitative path-sensitivity statement if supported independently, but do not use this excluded test as validation evidence for the fitted/tested configurations. The photograph remains archived. A replacement image is unnecessary unless it supplies essential evidence from the reported study.

**Verification boundary:** Do not relabel the photo as an included test or imply that all observed slippage is disproved. Check every surviving numerical displacement statement has an appropriate record independent of this removed illustration. No experiment or fit rerun is proposed.

**Prior audit:** Supplement relevance table's Figure S2 row and T15. Ben's relevance/provenance finding supersedes the earlier “Keep if linked” advice.

### U06 — Introduction-to-Methods segue: “the force model must see beyond,” and the omitted 27-site result

**Anchor:** Paper §1, P:99–101, pp2–3; §2.1, P:108, p3. Existing 27-site results are §4.4, P:637–656, pp27–28; abstract P:81, p1; conclusion P:678–680, pp29–30.

**Current state:** P:101 contains the exact anthropomorphic “force model must see beyond” wording. Its contribution list omits the 27-site comparison and mentions hardware evidence only for the extensor, although Figure 13 includes both redesigned routes. The whole-leg results already exist later: human-torque coverage 37.28% to 65.19%, active-force coverage 92.44% to 98.81%, with 23 sites improved and four tied. These are saved model comparisons, not whole-leg physical validation.

**Proposed change:** Replace the phrase with a direct description such as “an empirical contraction correction representing wrapping-related loss of usable actuator length”; leave its mathematical definition to §2.7. Shorten the closing introduction around the study sequence: identified joint corrections, redesigned knee routes with physical tests, and the subsequent 27-site transmission/configuration comparison. Add a short statement of that third scope and, if useful, its existing torque-coverage result. Scope the comparison to inherited routing and its combined transmission/actuator-sizing/baseline-diameter changes; do not attribute the gain solely to placement or corrected force modeling. Let §2.1 describe the workflow without repeating the full contribution list. Coordinate abstract/conclusion framing with prior T12 rather than creating another lengthy list.

**Verification boundary:** Reuse the current recorded figures; no new 27-site optimization or coverage calculation. Preserve the distinction between discrete knee hardware tests and modeled whole-leg coverage, and between selected single-muscle references and summed joint-muscle capability.

**Prior audit:** T12 partially overlaps the contribution-framing omission; this specifically adds Ben's register/27-site finding. The prior audit's whole-leg scoping strength remains applicable.

### U08 — Replace the paper's old pulley schematic

**Anchor:** Paper §2.9, P:397–404, Figure 7, p16. Current include is `Figures/reverse_pulley_principle.pdf` at P:401. The newer ZIP dissertation asset is `dissertation/figs/30-torque/reverse_pulley_biomimetic.pdf`, included by dissertation Chapter 4 at `30-torque.tex:780`.

**Current state:** The paper still shows the old two-part schematic with a separate block of equations, a bottom disclaimer, and no anatomical knee. The newer dissertation asset shows two parallel braided BPAs, a crossed-link biomimetic knee, total force `F_BPA = 2F_20`, tendon force `F_t = F_BPA/n`, and tendon travel `n Δl_BPA`. It is a suitable starting asset, subject to Ben's new notation and mounting-reaction simplification requests.

**Proposed change:** Use the revised biomimetic two-BPA schematic in the paper, with the same accepted mechanism geometry and typography as the dissertation. Carry equations in the surrounding text; retain only natural force/travel labels in the schematic. Begin the caption “Example of …” and keep it short. Resolve Ben's `F_R` notation and removal of the `R_mount` calculation/discussion consistently before sharing the asset. Do not copy the old caption, equation panel, or bottom disclaimer into the replacement. Place the figure in §4.4 as U09/U14 require.

**Verification boundary:** Preserve the distinction between the number of parallel BPAs and the number of load-bearing cable parts. Check the graphical reeving matches `F_t = F_BPA/n` and `Δs_t = n Δl_BPA`, with `F_BPA = 2F_20` in this example. Do not restore the waived individual vastus/pulley optimization example. This is a replacement proposal only.

**Prior audit:** New explicit paper carry-over directive; prior audit's reaction-load discussion is superseded where it conflicts with Ben's requested simplification. Related T08 organization advice must be applied without overriding the location specified in U14.

### U09 — Put pulley discussion after the meet/exceed-human-torque discussion

**Anchor:** Present mechanism discussion/figure is paper §2.9, P:397–404, pp15–16. Human-torque redesign and hardware results are §4.3.1–4.3.2, P:594–627, pp24–26. Transmission/whole-leg comparison is §4.4, P:629–660, pp26–28.

**Current state:** The full pulley explanation and schematic precede Results. This gives transmission mechanics before the reader has reached the knee torque shortfall, redesign outcome, and measured capacity that motivate it.

**Proposed change:** Move the substantive pulley explanation and the replacement schematic to the beginning of §4.4, after §4.3.2 has established the measured meet/exceed results. Use the existing transition at P:627 as the bridge. In Methods retain only the minimum procedure statement/cross-reference required to identify the evaluated transmission configurations; avoid duplicate mechanism explanations.

**Verification boundary:** One pulley schematic and one full mechanism explanation. Rebuild figure numbering/callouts automatically and check first-reference order after the move. No physical-validation implication for the transmission survey.

**Prior audit:** T08 concerns organization but did not establish this specific placement. Ben's ordering governs. Same proposed action as U14; track both findings without duplicating the actual move.

### U10 — Morrow (2020) cost may be mentioned without reproducing it

**Anchor:** Paper §1, P:97, p2; §2.4, P:166, pp6–7; §4.5, P:666, pp28–29. The paper's own objectives are P:388–393 and P:583–592.

**Current state:** **Already satisfied in the current ZIP source.** Morrow is cited for earlier placement/routing work; its cost function is not reproduced. Equations 35–37 are the present route-redesign objectives, not a reproduction of Morrow's cost.

**Proposed change:** Retain the concise citation and comparison. Do not introduce Morrow's equation during the planned rewrite, and do not remove this study's required objective description under this finding.

**Verification boundary:** Verify attribution against the cited source if the comparison is substantively rewritten. This item requires no new objective calculation or manuscript edit as it stands.

**Prior audit:** No separate T ID. Keep distinct from T08's location issue and U16's missing transmission-selection objective.

### U11 — Remove the hybrid calculation from Figure 11, caption, and discussion

**Anchor:** Paper §3 opening, P:410, p16; §3.3, P:499–505, pp17/22; Figure 11 on p22. Included graphic is `Ext10mm_52cm_current_publication.pdf` at P:503.

**Current state:** The gold approximate-hybrid points and their legend entry remain in Figure 11. The caption names them, the §3 opening qualifies them, and §3.3 says measured torque fell below them. This paper carry-over gap is visually confirmed. The dissertation has already removed this specific series.

**Proposed change:** Remove the approximate hybrid series/legend entry from Figure 11 while preserving its apparatus panel, measured series, baseline, optimized curve, and numerical data. Remove its caption phrase and corresponding §3/§3.3 comparisons. Keep the important biomimetic routing limitation in direct language without a withdrawn diagnostic. Do not remove the pinned-knee hybrid method or its evidence elsewhere: Ben's instruction here is specific to Figure 11.

**Verification boundary:** Inspect the exported figure, caption, introductory Results paragraph, §3.3, and any alt text together; search both the figure label and “hybrid” for stale biomimetic-extensor references. Compare retained plot data before/after; no fit/statistic change.

**Prior audit:** T06 now becomes an explicit requested paper change, rather than a discretionary carry-over suggestion.

### U12 — Paper §4.3.1 is too long

**Anchor:** Paper §4.3.1, P:581–602, pp24–25; p25 is a dense text page. Existing common procedure is §2.9, P:384–397; route coordinates/release data are Supplement §§1.2/1.4/1.6, S:145–173, S:204–251, and S:298–310, pp3/4–7/9–10.

**Current state:** This subsection repeats the common objective/constraints and design setup, lists release angles and endpoint coordinates in prose, then gives long direct-route/cam/flexor narratives. The current cam/required-margin definitions also retain prior T01's unresolved consistency question; editing length must not silently choose new numbers.

**Proposed change:** Reduce the main subsection to a claim-first result paragraph and a short interpretation of the direct, cam-assisted, and paired-flexor differences, supported by a compact capacity table. Remove repeated setup already defined in §2.9; cross-reference the existing supplementary geometry/release inventories rather than repeating coordinate and release-angle strings. Place the full configuration-specific penalty inventory in one Methods/supplement location, with a brief in-text summary. Preserve the existing physical mechanism and paired-route distinction. Do not add a table while leaving all its values repeated in prose.

**Verification boundary:** Every moved or removed numeric inventory stays traceable in the supplement/archive. Resolve the written margin definitions from saved records before presenting a final comparison table; no optimization rerun or ad hoc recomputation of accepted outputs. Avoid silently deleting constraints needed to interpret the design.

**Prior audit:** T08 and T09 directly overlap. This strengthens the compression action while preserving Ben's §4.4 placement for the later pulley/whole-leg explanation.

### U13 — In §4.3.2, state important results before showing them

**Anchor:** Paper §4.3.2, P:604–627, p26; Figure 13, P:616–625, p27. Main outcome sentences are at P:610 and within P:612.

**Current state:** The subsection first lists the extensor test protocol and figure contents; its success statement comes in paragraph 2. The flexor paragraph starts with sequence/protocol and a rigid-versus-corrected prediction before arriving at the all-15-point hardware result. The main evidence is present but buried.

**Proposed change:** Open with the combined result: all nine extensor measurements and all fifteen flexor measurements meet or exceed their respective named human-muscle references at the tested angles. Follow with the measured margins and Figure 13 reference, then give only the test details needed to interpret them. Keep the discrete-point scope in one concise sentence. Relocate the installed-versus-model length inventory to a clear table/supplement reference after its definitions are reconciled, rather than allowing that inventory to obscure the result.

**Verification boundary:** Preserve all test counts, angles, pressures, and margins. Keep measured capacity separate from prediction agreement: the extensor's −86.5° point remains below its pressure-specific prediction, and the sparse tests do not establish a measured continuous full-pressure envelope.

**Prior audit:** T08 and T12 overlap; T10 remains a parameter-definition dependency. Claim-first presentation does not require moving Ben's specifically named dissertation section.

### U14 — The pulley schematic belongs in paper §4.4

**Anchor:** Paper §4.4, P:629–643, pp26–28; current schematic remains §2.9, P:399–404, p16.

**Current state:** §4.4 already contains the selected transmission comparisons and 27-site whole-leg survey. It is the natural point at which the force–travel mechanism explains the subsequent comparison. The existing figure location is inconsistent with Ben's requested sequence.

**Proposed change:** Make §4.4 the single canonical home of the replacement schematic and concise force–travel explanation, immediately following the knee capacity evidence and before the transmission comparison details. Share this one implementation action with U08/U09. Keep the master's-model attribution at P:631/639 visible.

**Verification boundary:** After eventual edits, verify only one active include/label and no forward reference stranded in the removed explanation. Check figure floats do not interrupt the knee results or appear before their intended introduction.

**Prior audit:** Explicit new placement directive related to T08; not a blanket endorsement of moving all cost formulas/results to an agent-selected structure.

### U15 — Figure 13 wording, matched plots, and panel-B scatter order/style

**Anchor:** Paper §4.3.2, P:606–623, pp26–27, especially graphic includes P:619/622 and caption P:623. Related “cam-assisted” wording is at P:600 and P:664. Figure 13 uses `Ext20_current_publication.pdf` and `Redesigned_flexor_pressure_publication.pdf`.

**Current state:** The caption says “Cam extensor.” The two included plots visibly use different final font/line scales; the lower plot has much larger labels. Follow-up inspection of the actual paper PDF confirms blue measurements are painted before the gold predictions; the predictions are opaque and have separate dark marker-edge strokes. The paper asset's SHA256 starts `79814B2A`, whereas the corrected dissertation asset starts `0721D2DD`. Thus the dissertation pass does not clear the paper. Detailed path, sequence, fill, alpha, stroke, and bounds evidence is saved in `B6_CURRENT_PDF_MARKER_EVIDENCE.json`. This follow-up closes the earlier export-property uncertainty without regenerating or recalculating either plot.

**Proposed change:** Use “cam-assisted extensor route” in caption/prose and a correspondingly concise plot label. Re-export both plots with common figure dimensions, usable axes area, font family/sizes, and line widths at their final inclusion size. For panel B draw the pressure-specific predictions before measured points, with measured on top. Set the scatter series to `MarkerFaceAlpha = 0.75` and `MarkerEdgeColor = 'none'`; set the legend order explicitly rather than relying on object creation. Preserve the human references, pressure-specific predictions, and the stated pressure-size encoding in panel A if it remains useful.

**Verification boundary:** Verify actual scatter properties/object order and the final exported page, not just legend text. Match final rendered typography, and compare every retained X/Y datum and curve before/after. Do not regenerate predictions or change model inputs to achieve visual consistency.

**Prior audit:** T07; Ben's specific order/alpha/no-edge requirements are additional to the previous readability finding.

### U16 — Give §4.4's cost function here or in the supplement with a table

**Anchor:** Paper §4.4, P:631–643, pp26–28; selected “saved transmission-ladder winner” at P:639. Coverage definitions are P:637/641. Supplement survey table is §1.5, S:255–296, Table S6, p8. The route costs at §2.9 P:388–393 and §4.3.1 P:583–592 are different objectives.

**Current state:** The manuscript defines reported torque/active-force coverage and describes re-solving lengths, but does not specify the objective/ranking rule that selected each transmission winner. Table S6 lists sites, primary/secondary joint coordinates, and actuator sizes; it is not an objective/weight/selection table. The knee route cost cannot automatically stand in for the transmission-ladder selection rule.

**Proposed change:** Recover the actual existing search/selection rule from the saved campaign code and records. State a brief account in §4.4, and place the exact objective/ranking expression with a compact table of terms, units/normalization, weights, candidate counts/ratios, constraints, and tie-breaks in the supplement if the full definition would interrupt the narrative. Include only quantities actually used. If any final choice was engineering judgment rather than a scalar optimizer, say how that choice was made; do not manufacture a formal objective or a new assessment procedure. Cross-reference that material from §4.4 and the survey table.

**Verification boundary:** This pass verifies the manuscript omission, not the saved selection algorithm. Read the existing code/records before drafting the mathematical rule. No rerun, reoptimization, new metric, new weights, or retroactive claim of a documented assessment method is authorized or needed.

**Prior audit:** New explicit reproducibility/selection-rule item. T08 covers only the different route-cost location; T09 covers compression. Do not merge these distinct costs merely because both occur near the redesign discussion.

## Proposed action groups and completion boundary

All twelve requested IDs above have an evidence-based response. They consolidate into: supplement relevance cleanup (U01/U02), contribution/segue rewrite (U06), one pulley replacement and move (U08/U09/U14), preservation of the already concise Morrow citation (U10), Figure 11 carry-over correction (U11), §4.3 compression and claim-first ordering (U12/U13), consistent Figure 13 export (U15), and an accurate account of the existing transmission selection rule (U16).

These are ledger decisions and verification criteria only. No manuscript implementation or acceptance of tracked changes occurred. The other audit file was not read for this task, preserving Ben's requested sequence.


---

# Ben's dissertation mechanics and Appendix C findings

30 October 2026. Task 2 only: audit and proposed-change ledger. No manuscript, figure, optimization, or numerical result was edited or recomputed. The source is the latest Overleaf ZIP, isolated in `Archive/audit_20261030/dissertation`. Source lines below are one-based and refer to that snapshot. `30-torque.tex` is Chapter 4; `94-AppendixC.tex` is Appendix C. Section, figure, and table numbers were checked against the fresh build's contents and labels.

This ledger incorporates the assigned findings from Ben's latest list. They are **Ben's directives/observations**, not newly attributed advisor comments. His instruction controls where the earlier audit recommended retaining something he now wants removed. The other audit comparison is a later task and is not performed here.

Current figure checks used the existing fresh-build page renders for printed pp.104,114,116 and the PDF skill's read-only inspection workflow. For Figure 4.14, the latest ZIP asset's vector drawing order, transparency, and absence of marker strokes were also inspected. Its SHA256 matches the repository asset: `0721D2DDFC6D354E1736FDF08245C50868C08AE2D45A55F3F090978DEE296AA5`.

## U03 — Table C13 repeats the main-chapter fit table

**Anchor/current state:** Appendix C, Section C.9, `94-AppendixC.tex:433–450`, Table C13 (`tab:TorqueFit`); Chapter 4, Section 4.2.4, `30-torque.tex:617–632`, Table 4.8 (`tab:BioTorqueGof`). The three configurations, both groups of metric columns, and all numerical entries are identical. Confirmed redundancy.

**Proposed change:** remove Table C13 and its introductory sentence. Use a short reference to Chapter 4 Table 4.8 only where the appendix needs to identify those statistics. Keep the complete pinned-knee tables C11/C12, which contain a different inventory.

**Verification needed:** check every `tab:TorqueFit` reference and replace it with the retained label as appropriate; rebuild labels, table list, and pagination. Compare the retained Table 4.8 against the removed duplicate to ensure no unique datum is lost. No metrics are changed.

**Prior audit overlap:** A 4's broader duplication/compression finding. This specific duplicate was not identified separately before Ben's observation. His instruction narrows the earlier general recommendation to retain authoritative tables: keep one authoritative copy.

## U04 — Section C.9.1 repeats the goodness-of-fit definitions

**Anchor/current state:** Appendix C, Section C.9.1, `94-AppendixC.tex:429–430`, repeats the three metric definitions already given in Chapter 4, Section 4.1.4, `30-torque.tex:316–336`, Equations `eq:gof_components`, `eq:gof_fvu`, and `eq:gof_max`. With U03, C.9 would contain no unique substantive material.

**Proposed change:** delete the duplicated subsection and, after removing Table C13, remove the empty C.9 container. Point any genuine metric-definition reference to the existing Chapter 4 definitions. Do not revise the metric, its interpretation, Ben's reply, or any optimization.

**Verification needed:** update references to `sec:app_torque_fit_details`; verify there are no dangling links or empty headings. Allow subsequent section numbers to update automatically.

**Prior audit overlap:** A 4 on repetitive model descriptions; new specific deletion directed by Ben. This is not a reopened FVU task.

## U05 — Section C.10 needs plain language and the master's-project provenance

**Anchor/current state:** Appendix C, Section C.10, `94-AppendixC.tex:452–458`, and Table C14 at `:456–495`. The section opens with “graded surface” and “multiarticular route” without identifying Ben's master's robot model. Its table caption also omits that origin. Chapter 4, Section 4.3.6, `30-torque.tex:792`, already identifies “my master's model” and cites `bolen_determination_2019`.

**Proposed change:** lead C.10 with the inherited robot-model provenance and the comparison's purpose. For example: “This comparison uses the 27 actuator sites of my master's robot-leg model [existing Bolen 2019 citation]. Table… lists the joint coordinates evaluated at each site and the actuator sizes used for the baseline comparison.” Replace “graded surface” with “evaluated joint-coordinate grid,” using U27's brief definition. Give Table C14 a caption that explicitly connects the inventory to the master's model. Retain the distinction between nominal 40-mm sites and four parallel 20-mm baseline actuators.

**Verification needed:** use the same master's citation and site inventory as Chapter 4; ensure the caption does not imply that the optimized transmissions were already implemented in the original master's hardware or converted neuromechanical model. Preserve the existing coordinates, dimensions, and coverage values.

**Prior audit overlap:** A 4's jargon/implementation-history compression and its inventory-retention assessment. The master's-provenance omission in this particular section is Ben's new finding. U27 and U28 share the explanatory paragraph and should be implemented together.

## U17 — Use F_R for the retained reaction-force notation

**Anchor/current state:** Chapter 4, Section 4.3.5, `30-torque.tex:763–774`; Figure 4.15 is included at `:780`, caption/label `:781–782`. The current graphic visibly labels `R_origin=F_BPA` and `R_mount`, while the text uses bold `R` for cable/mount reactions. The requested `F_R` notation is not present.

**Proposed change:** use `F_R` for retained reaction-force labels, with a location subscript only if needed, such as `F_{R,origin}`. Coordinate this with U26: remove the mount-reaction calculation and discussion, rather than renaming and retaining that derivation. The origin reaction can remain as `F_{R,origin}=F_BPA`. Remove a graphic label/arrow that exists only to support the deleted mount-reaction treatment, or use `F_R` if a descriptive reaction arrow is retained; do not add a replacement mount-force calculation.

**Verification needed:** inspect both the editable schematic and its exported PDF at final size; search retained Chapter 4 text for inconsistent standalone reaction `R` symbols. Preserve `F_BPA=2F_20`, the crossed knee linkage, force/travel transmission, and the “Example of…” caption opening.

**Prior audit overlap:** the earlier torque audit positively described the origin/mount distinction under “What already works.” That favorable assessment does not override Ben's new notation and deletion directives. U26 controls the scope.

## U18 — Report Tables C5 and C6 dimensions in meters

**Anchor/current state:** Appendix C, Section C.1, `94-AppendixC.tex:105–130`. Table C5 headers at `:111` and C6 headers at `:126` give diameter, resting length, and tendon length in mm. Tables C1/C4 already use meters for resting/tendon length. The appendix introduction at `:4` also makes the overly general statement that coordinates are reported in millimeters.

**Proposed change:** convert numerical linear-dimension columns in C5/C6 and their headings to m, preserving pressure in kPa and force in N. Explicit conversions: C5 resting lengths 415.0 mm become 0.4150m; tendons 11.0/14.0 mm become 0.0110/0.0140m; numeric diameters10/20 mm become 0.010/0.020m. C6 resting/tendon520.0/38.0 mm become 0.5200/0.0380m and diameter 10 mm becomes0.010m. Nominal “10 mm”/“20 mm” model names may remain nominal labels. Qualify the appendix-wide unit sentence so it reflects the actual tables rather than asserting a single unit for all quantities.

**Verification needed:** perform a cell-by-cell factor-of 1000 check for the converted dimension columns and verify the rendered headers and significant digits. Do not alter pressure, force, test identity, or effective-versus-measured length meanings. These are explicit display-unit conversions, not new experimental values.

**Prior audit overlap:** A 4's unit/configuration-specification concern at a broad level; the C5/C6 meter conversion is a new, specific Ben directive.

## U19 — Remove repeated compliance equations from C.4

**Anchor/current state:** Appendix C, Section C.4, `94-AppendixC.tex:205–214`, repeats projected/series compliance in Equation C5 and explicitly says it is “exactly as in” the chapter equations. Chapter 4 has the stiffness relation at `30-torque.tex:140–152`, equivalent stiffness near `:214–232`, and force balance at `:249–250`.

**Proposed change:** remove the repeated displayed derivation from C.4 and reference the existing Chapter 4 equations. Retain the unique per-configuration stiffness arrays, parameter bounds, units, and authoritative implementation mapping. Do not remove the appendix's unique geometry or contact-release calculation merely because this compliance equation is repeated.

**Verification needed:** redirect references to `eq:app_compliance`, confirm the chapter's retained equations define every quantity used in C.4, and keep the independent frame-convention reconciliation A 1 visible. No parameter or evaluator change is authorized by this deletion.

**Prior audit overlap:** A 4 and T13. A 4's general advice to retain useful equations is superseded for this specific duplicate by Ben's direction; T13's compression advice should result in one core presentation plus genuinely unique supporting details, not repeated copies in Chapter 4 and Appendix C.

## U20 — Put FLEX/EXT weights on separate displayed lines in Table 4.3

**Anchor/current state:** Chapter 4, Section 4.1.5, `30-torque.tex:417–434`, Table 4.3 (`tab:mesh_geometry_penalties`). Both strain rows at `:425–426` place `10^7 FLEX, 10^4 EXT` in one weight cell.

**Proposed change:** display two stacked lines in each affected weight cell: `FLEX:10^7` and `EXT:10^4`. A multiline cell is sufficient; duplicating the whole semantic row is unnecessary. Preserve each existing weight and its association with the configuration.

**Verification needed:** inspect the compiled table at final size to confirm both labels begin on their own lines, remain aligned with the correct penalty, and do not cause column overflow.

**Prior audit overlap:** no exact prior finding; new table-formatting directive from Ben. This is not a request to change or rerun the objective.

## U21 — Number the constrained-redesign heading and shorten its tail

**Anchor/current state:** Chapter 4, within Section 4.1.5, `30-torque.tex:374` is an ordinary `\subsubsection{Constrained Knee-Route Redesign}` and lacks a visible number in the current build/contents. It is not a manually missing number in the title string. The material continues through `:466`; the tail at `:451–466`, especially the single paragraph at `:464`, combines targets, clearance rules, solver budgets, cam motivation, margins, paired-route arrangements, and design trade-offs.

**Proposed change:** make the constrained redesign a clearly numbered Methods subsection, keeping the earlier mesh-search material in 4.1.5 and promoting this material to 4.1.6 with a dedicated label. Retarget relevant redesign links. This avoids a global numbering/style change to the protected published Chapter 3. Keep the cost definitions and tables; shorten the ending to the essential configuration distinctions and one Results handoff. Put unique budget/configuration inventory in an existing settings/parameter table or Appendix C rather than repeating it in prose. Identify every removed or relocated statement in the implementation ledger.

**Verification needed:** check ToC hierarchy, numbering, all mesh-versus-redesign references, and final rendered page flow. Preserve the two separately evaluated flexor routes, extensor force scaling,5% required allowance, and the actual solver settings wherever they are recorded. No new optimization or tolerance selection is requested.

**Prior audit overlap:** T09's dense-prose diagnosis applies to the paper; T13 addresses the chapter's repeated mechanics algebra. This specific heading/tail issue is Ben's new chapter finding. A 4 overlaps only if details are relocated to Appendix C; do not create a second duplicate there.

## U22 — Figure 4.10 legend belongs northeast inside panel B

**Anchor/current state:** Chapter 4, Section 4.2.3, `30-torque.tex:594–599`, asset `Ext10mm_52cm_current_publication.pdf`. Printed p.104 visibly has a shared horizontal legend below both panels, ordered Measured → Improved BPA → Optimized model. It is not northeast in panel B.

**Proposed change:** move the legend to panel B's northeast corner and order the entries consistently with the approved model-first torque figures: Improved BPA → Optimized model → Measured. Preserve the apparatus photograph, retained curves/points, and removal of the hybrid series.

**Verification needed:** compare the entry order with the approved current torque figures and inspect overlap/legibility at final inclusion size. Confirm data arrays are unchanged. Re-export the figure only; no model computation is needed.

**Prior audit overlap:** T06 previously treated the dissertation's outside-data legend as an improvement over the paper. Ben's new northeast directive supersedes that earlier favorable placement assessment. The earlier hybrid-removal conclusion remains valid.

## U23 — Table 4.10 needs consistent percentages, the flexor angle, and distinct achieved/required labels

**Anchor/current state:** Chapter 4, Section 4.3.3, `30-torque.tex:669–685`, Table 4.10. Rows `:680–681` mix extensor percentages with flexor fractions (`0.0501`, `0.05`), the extensor margin includes a plus sign, and the flexor worst-case angle at `:682` is blank. The required allowance is defined at `:383–390` and explained at `:464`.

**Proposed change:** express both columns in percentages and remove the explicit plus sign: the displayed fraction 0.0501 converts to 5.01%, and 0.05 converts to 5%. Use labels that distinguish the **achieved minimum predicted margin** from the **required design allowance**, naming their reference torque once the saved definition is checked. Retrieve the flexor worst-case angle from the existing adopted direct-route result; do not guess it. Keep predicted versus required torque clearly distinguished in the other rows.

**Verification needed:** check the existing saved margin definition, denominator, evaluated angle domain, and worst-angle field before naming the achieved margin “over human” or “over required.” The intended design requirement is 1.05 times the human target; that requirement alone does not establish how a saved margin field was normalized. Verify the flexor angle from the existing result/generator. Do not substitute a new calculation or rerun the search. Until verified, that angle remains an explicitly pending entry in the ledger.

**Prior audit overlap:** related to T01's margin-definition reconciliation, but T01 identifies a different issue in the cam Table 4.11. Do not collapse the direct-table presentation/missing-angle issue into the cam inconsistency or silently transfer a cam angle to the flexor.

## U24 — Pressure headroom also has a lower operating limit

**Anchor/current state:** Chapter 4, Sections 4.3.3–4.3.4, `30-torque.tex:715–718,729–733`. Table 4.11 lists269–507 kPa to meet the required envelope and a 150-kPa “regulated pressure floor.” The discussion does not explain Ben's observed manifold operating threshold of approximately 200 kPa or the resulting minimum controlled torque magnitude.

**Proposed change:** add a brief hardware-qualified interpretation: pressure headroom can reduce required supply pressure, while the manifold's pressure needed for reliable opening/closing can restrict the minimum controllable nonzero torque when the BPA is taut. Identify approximately 200 kPa as Ben's operating experience with the present manifold. Explicitly distinguish that observation from the 150-kPa assumption used in the saved design calculation.

**Verification needed:** retain the saved 150-kPa design/model constraint as a recorded assumption; verify any hardware-source reference already available for the reported 200-kPa observation. Do not change the calculation to 200 kPa, assert a universal BPA threshold, or produce a new torque envelope. If only Ben's experience supports it, say so directly rather than inventing a documented assessment method or test protocol.

**Prior audit overlap:** no exact prior audit finding. T01 remains a separate margin-definition check. This is Ben's supplied hardware observation and requested discussion, not a newly quantified result from this audit.

## U25 — Figure 4.14 prediction points before measured points

**Anchor/current state:** Chapter 4, Section 4.3.4, `30-torque.tex:742–747`, asset `Redesigned_flexor_pressure_publication.pdf`, printed p.114. **Already implemented for drawing order.** The current ZIP PDF paints the 15 gold prediction points at drawing sequence 72–86 and the 15 blue measured points at 87–101, so measurements draw last/top. Both scatter sets have fill opacity 0.7490196, the 8-bit PDF representation of the requested0.75, and no marker strokes. The retained local MATLAB editing record explicitly sets `MarkerFaceAlpha=.75`, `MarkerEdgeColor='none'`, and puts the measured group on top. The asset hashes match.

**Proposed change:** retain this verified implementation. No replot or data change is needed for this directive. The legend still lists Measured before Prediction; legend order is separate from draw order. Reorder those legend entries only if applying Ben's broader desired display-order consistency, and record it as a presentation change rather than claiming the draw-order request was missed.

**Verification needed:** after any later figure revision, preserve the measured-last stacking, marker settings, and data arrays. Existing evidence is the current asset plus `Figures/30-torque/Components/directives_20261029/edit_existing_plots.m` and `plot_verification.txt`; the PDF itself was checked here rather than relying solely on that older report.

**Prior audit overlap:** consistent with the earlier acknowledged successful physical-plot improvements and the 29 October directive record. This item should be marked verified complete, not added as a new defect.

## U26 — Delete R_mount calculation/discussion; use a short pulley explanation and RoM

**Anchor/current state:** Chapter 4, Section 4.3.5, `30-torque.tex:753–776`, contains the fixed-mount vector equation at `:763–768`, reaction-magnitude equation at `:769–773`, and numerical mount-reaction example at `:774`. The discussion continues at `:786–788`; “angle range” occurs at `:788`.

**Proposed change:** delete the `R_mount` calculation, vector/magnitude equations, and accompanying discussion as Ben directs. Do not replace them with an alternative reaction derivation. Reduce the section's main explanation to approximately two paragraphs around the schematic/ideal transmission equation: first define the two-BPA input, the reeving count, force loss, and tendon travel; second state what that trade-off contributes to RoM/capacity and refer detailed route searches to the existing appendix. Retain `F_t=F_BPA/n`, `Delta s_t=n Delta ell_BPA`, `F_BPA=2F_20`, and any retained origin-reaction label using U17's `F_R` notation. Replace “angle range” with “RoM.” Consolidate existing predicted pulley outcomes without silently deleting unique results; move supporting details to their existing appendix account when necessary.

**Verification needed:** confirm no mount-reaction equation, numeric example, or discussion remains in the section/caption after editing; update deleted equation references; inspect the matching schematic labels. Ensure the final section has the requested short structure and still distinguishes predicted pulley performance from physical verification.

**Prior audit overlap:** supersedes the earlier torque audit's favorable description of the detailed support-reaction discussion. The earlier instruction was to explain those reactions; Ben's latest deletion instruction controls now. No objection or re-derivation is proposed.

## U27 — Briefly define the evaluated or “graded” range

**Anchor/current state:** Chapter 4, Section 4.3.6, `30-torque.tex:792`, says each site is “graded” on a 5×5 grid but does not give a plain-language term for the evaluated range. Appendix C, Sections C.10/C.10.1, `94-AppendixC.tex:454,458,500`, uses “graded surface,” “graded coordinates,” and “graded domain.”

**Proposed change:** prefer “evaluated joint-coordinate grid” and define it once: the sampled values of the two joint coordinates over which each site's torque capacity is compared with its human-muscle reference. Explain briefly that reported coverage counts the included grid cells meeting the specified torque direction and magnitude; it is not a claim about every possible multi-joint pose or the whole physical joint RoM. Use a short cross-reference in Appendix C rather than repeat the definition.

**Verification needed:** match the existing grid/domain masks and coordinate bounds without inventing an assessment protocol or changing the denominator. If actual numerical range limits are added, retrieve them from the saved comparison, not assumed physiological limits. Preserve the distinction between torque coverage and positive-tension coverage.

**Prior audit overlap:** A 4's jargon compression and T09's broader prose-reduction recommendation; the explicit graded-range definition is a new Ben directive. Implement with U05 to avoid two new explanations.

## U28 — State the ideal force equivalence of four 20-mm BPAs and one 40-mm BPA

**Anchor/current state:** Chapter 4, Section 4.3.6, `30-torque.tex:792`; Appendix C, Table C14 caption at `94-AppendixC.tex:458` and rows `:465–489`. These identify the four 20-mm baseline substitution at eight nominal 40-mm sites but do not explain the theoretical force-equivalence rationale.

**Proposed change:** add one concise statement that, under matched ideal pressure, braid geometry, and contraction conditions, four parallel 20-mm BPAs have the same theoretical maximum force as one 40-mm BPA because the ideal force scales with diameter squared. Present this as the rationale for the substitution. Do not identify it as a measured equality of the actual actuators, an equality of their complete force–length maps, or an isolated explanation of the measured coverage gain.

**Verification needed:** ensure the ideal-scaling assumption and existing comparison definitions are explicit; preserve the actual actuator-size inventory and the existing statement that the comparison combines transmission, sizing, and diameter substitution. No new40-mm fitting or optimization is required. The protected published Chapter 3 is not edited.

**Prior audit overlap:** the earlier torque audit correctly preserved the mixed-comparison scope but did not state this theoretical substitution rationale. This is Ben's additional explanation, coordinated with U05/U27 rather than repeated in three places.

## Implementation boundary and completion check

All 15 assigned IDs are addressed: U03, U04, U05, U17, U18, U19, U20, U21, U22, U23, U24, U25, U26, U27, and U28, with U03/U04 deliberately separate despite sharing a deletion, and U17/U26 deliberately coordinated. U25 is verified already implemented. The remaining entries are proposed changes only, with the existing-result flexor worst-case angle and exact margin normalization explicitly left for provenance verification. No manuscript change is represented as completed, and no new calculation is used to fill a missing value.




---

# Ben's SADb findings and advisor-comment crosswalk

Audit and proposed-change ledger, 30 October 2026. **No manuscript, figure, application, or live review thread was changed.** This is Task 2, following completion of the Chapter 2 audit. The separate `AUDIT_ADVISOR_PAIN_20261030.md` comparison is not part of this file.

## Evidence and interpretation of status

Current manuscript evidence is the supplied ZIP extracted into `Archive/audit_20261030/dissertation`, rather than an older working manuscript. Source references below are one-based lines in its `chapters` files: `50` = Chapter 5, `50-database.tex`; `70` = Chapter 7, `70-overall.tex`; `30` = Chapter 4, `30-torque.tex`; and `96` = Appendix E, `96-AppendixE.tex`. Printed page references follow the audit build's table of contents, not a historical build's pagination.

The 42 comment cards are from the saved 29 October review ledger, `Notes/committee_revision_20261027/evidence/ch45_annotated_review_20261029/revised_review_ledger.json`. Ben's additions to the Word review were re-read through the companion `ben_plain_text_annotations.json`, including the confirmed paragraph anchors, the instruction to leave the offset and FVU items alone, the existing advisor deletions, and the deferred classification review. The original ledger preserves the comments and captured replies; the table below rechecks their conditions against the latest ZIP instead of inheriting an old “implemented” label.

**“Addressed in ZIP” describes the available source, not a resolved live Overleaf thread.** The saved cards date from 29 October. They cannot establish whether additional comments were added on 30 October, whether an anchor is still live, or whether a live thread is resolved. A standard source ZIP does not establish that review history. Those limits remain explicit in U32b.

FVU recalculation, optimization reruns, a new AI-tool assessment method or comparison audit, and the larger database classification exercise remain excluded under Ben's instructions. Comment 5.5 concerns stating the explorer's status as the primary deliverable; it is **not** a request for a documented assessment method.

## Ben's findings: current state and proposed changes

| ID | Finding and verified current state | Proposed disposition and evidence |
|---|---|---|
| U29 | **Figure 5.2 has excessive bottom whitespace. Confirmed.** The ZIP's `figs/50-database/sadb_pivot_app.png` is 1800 × 1100 px; the pivot table and controls occupy approximately its top 475 px, leaving most of the lower image blank. Figure 5.2 is on printed p.125. | **Propose cropping the asset**, with a small margin under the last table row, then sizing it for readable labels. Preserve the data and controls needed to understand the pivot, and inspect the final page for legibility. The image's obsolete header and technical rebuild text can be cleaned up in a replacement capture if one is made; cropping alone must not silently change the reported corpus values. Chapter 5:39–44, especially the include at 50:41. No crop performed in this audit. |
| U30 | **§5.3.2 contains future work, but its last paragraph supplies a useful segue. Confirmed.** Lines 63–65 discuss completing classifications, adding primary-source/pathway links, expanding vestibular/visual coverage, and requirements for new records. Line 67 connects evidence selection to implementation and testing in Chapter 6. | **Propose consolidating future improvements into Chapter 7 §7.2.3**, whose database future-work subsection already exists at 70:59–61. Remove the now-redundant Chapter 5 “Further Improvements” heading at 50:61 if its future-work paragraphs move. Retain the substantive segue at 50:67 as the close of the Chapter 5 Discussion. A short present-coverage limitation may remain where needed to interpret the Results; moving future work does not authorize hiding the existing coverage counts. Preserve the broad-to-specific uses discussion at 50:55–59. |
| U31 | **The monthly literature watch is built. Confirmed by implementation and a saved run.** Chapter 5 does not currently describe it, while Chapter 7 incorrectly says that an extension “is being implemented” at 70:61. `SADb_audit/watch/monthly_watch.py` contains the retrieval, deduplication, verification, and report workflow; `watch/reports/2026-10.md` and `2026-10_shortlist.json` are saved outputs. | **Propose describing the implemented watch briefly in Chapter 5 Methods**, after the curation workflow at 50:14–16, and changing 70:61 to describe the remaining extension/curation work rather than promising to build the existing script. Suggested content: “I implemented a monthly literature watch that queries OpenAlex and Crossref for papers matching selected authors and keywords, removes DOI and version duplicates, and produces a DOI-verified shortlist for manual curation.” Preserve human inclusion decisions. Do not claim an installed unattended schedule: `watch/README.md:110–114` explicitly says the scheduling command is provided but not installed. Do not claim complete OpenAlex coverage of the saved run: its run notes report an OpenAlex access limit and Crossref fallback. The built script is established without running a new search or inventing a validation campaign. |
| U32a | **Show app pictures earlier in Chapter 5. Supported.** Chapter 5 begins on printed p.119, while the large bubble-map figure appears on p.124 and the pivot on p.125. The Table screenshot remains in Appendix E, on p.283, despite advisor comment 5.23. | **Propose an early, current Table-tab example** next to the introductory app overview at 50:8–10 or the first Table workflow at 50:24. Capture a retained, relevant sensory-feedback record or filtered list with useful detail visible. Do not simply move the existing Appendix E image: `figs/96-AppendixE/sadb_explorer.png` visibly includes generic deep-learning and turbulence entries and an older header. Keep the existing bubble-map figure's six callouts and app-functionality explanation at 50:26–35. Keep the cropped pivot next to its Results explanation at 50:37–44. Avoid repeating an identical screenshot in both locations unless each use has a distinct purpose. |
| U32b | **Check that all advisor comments from today are dealt with, except Ben's exclusions/deferments. Partially verifiable, not certifiable from the saved evidence.** All **42 historical cards** are accounted for below. The latest ZIP still has the unaddressed Table screenshot/main-text request (5.23), and cards 4.1 and 5.25 retain anchor uncertainty. | **Ledger completion is achieved; blanket live completion is not claimed.** Retain the specific pending actions and uncertainties below. A fresh live 30 October comment/revision capture is needed to certify the word “all” for today's comments. Keep 4.3 unchanged, 4.5 handled by Ben's reply, and 5.17 deferred. Do not accept tracked changes or resolve threads merely because the old ledger said “implemented.” No live resolution or implementation was performed here. |

Implementation evidence for U31: `SADb_audit/watch/README.md:3–8` defines the monthly workflow and human review; lines 72–92 describe source queries, deduplication, Crossref verification, and output; `monthly_watch.py:394` begins the main workflow, and `monthly_watch.py:638` begins report rendering. The October report's method, saved shortlist, and run notes corroborate an actual run. The implementation uses the two named APIs; calling it the “monthly scraper” does not establish Google Scholar or Web of Science scraping.

## Chapter 4: all eight saved cards

The small quoted prompts identify the cards; full captured messages and replies remain in the original ledger. Each row states the current source condition and the allowed disposition.

| Card | Advisor issue and Ben's reply/ruling | Current ZIP evidence and disposition |
|---|---|---|
| 4.1 | “where is this?” Ben replied that he had uploaded the figure whose points should include pink. | **Ben replied; exact figure identity unverified.** The thread does not name the figure, so assigning a line or claiming the specific request complete would be speculative. Preserve the reply; recover the live anchor before certifying this item. No new upload proposed merely to satisfy an unidentified card. |
| 4.2 | Important result first, then show. Ben confirmed the pinned-extensor paragraph anchor. | **Addressed in ZIP.** Chapter 4:549 leads with reduced mean torque-prediction error in training and validation, then introduces the table/figure and interpretation. This meets the comment's narrative-order request. Further shortening elsewhere is a separate proposal. |
| 4.3 | Advisor suggested a zero offset. Ben asked what he meant and annotated “do not change.” | **Leave unchanged under Ben's ruling.** Chapter 4:549 retains the original offset observation; no new experimental-zero explanation or fit is requested. The advisor's unresolved interpretation does not authorize an additional analysis. |
| 4.4 | “in what equation?” Ben specified `eq:deltaL`. | **Addressed in ZIP.** Chapter 4:549 explicitly identifies the equation containing the wrapping term with `Equation~\eqref{eq:deltaL}`. |
| 4.5 | Advisor questioned FVU values above one. Ben replied and explicitly directed no further action. | **Handled by Ben; excluded.** FVU reporting remains in the torque results, including Chapter 4:592 and 647. Preserve his reply and existing metric/objective/results. No optimization, refitting, replacement metric, or additional FVU explanation is proposed. |
| 4.6 | “Unnecessary detail.” Ben said the advisor had already deleted the passages and that the earlier review missed the tracked deletions. | **Deleted passages absent in ZIP.** The biomimetic-extensor Results now proceed through the concise test/result paragraph at Chapter 4:592 and the figure at 594–602. The two previously recorded deleted explanatory passages are not present. This verifies the effective source outcome; it does not accept or reject a live tracked change. |
| 4.7 | Remove hybrid points from the biomimetic-extensor graph. Ben directed replotting/reassembly with apparatus (A), preserving readable legend placement. | **Removal addressed in dissertation ZIP.** Chapter 4:596 includes `Ext10mm_52cm_current_publication.pdf`; the surrounding result paragraph at 592 and caption at 597 do not describe the removed hybrid comparison. The existing dissertation figure has the retained apparatus and non-hybrid data; the figure audit supplies visual checking. Ben's newer request for a northeast legend is separately retained in the new-findings ledger, rather than claiming this old card proves that placement is already done. Paper Figure 11 is a separate deliverable. |
| 4.8 | Important result first, then show. Ben specified the biomimetic-flexor paragraph. | **Addressed in ZIP.** Chapter 4:606 leads with the error reduction in both tested configurations, then points to Figure `fig:FlxGrp` and the error table. The later shortfall description remains. |

## Chapter 5: all 34 saved cards

| Card | Advisor issue / Ben's clarification | Current ZIP evidence and disposition |
|---|---|---|
| 5.1 | Incomplete transition from existing reviews; Ben supplied the precise anchor. | **Addressed in ZIP, 50:6.** The transition names concrete tasks: comparing preparations, locating measurements, and following evidence to controller roles. |
| 5.2 | Better word than “following” groups of papers. | **Addressed in ZIP, 50:8.** The overview describes exploring citation links among papers. |
| 5.3 | “huh?” about offline views. | **Addressed in ZIP, 50:8.** Records, notes, and the citation graph are embedded, and the three offline tabs are named. |
| 5.4 | “huh?” about online search. | **Addressed in ZIP, 50:8 and 24.** The Search (online) tab is named, OpenAlex retrieval and the internet requirement are explained, and candidates are distinguished from curated additions. |
| 5.5 | Awkward wording; identify the explorer as the primary deliverable. | **Addressed in ZIP, 50:8–10.** The application is explicitly the primary deliverable, and Airtable's editable-source role is distinguished. The advisor did not request a documented assessment method. No assessment or AI-comparison task is added. |
| 5.6 | Were Zotero duplicates combined? Ben: yes. | **Addressed in ZIP, 50:14.** DOI/title matching and combining entries describing the same publication are stated. |
| 5.7 | What does “citation search” mean? Ben described references and citing papers. | **Addressed in ZIP, 50:14.** Backward reference tracing and forward citing-paper searches replace the unexplained term. |
| 5.8 | “huh?” about reviews guiding the search. Ben clarified that role. | **Addressed in ZIP, 50:14.** Reviews guide the search for primary studies. |
| 5.9 | Why import records and then remove them? Ben explained broad AARL/personal Zotero imports followed by trimming. | **Addressed in ZIP, 50:14–16.** Both collection sources and the initial broad import followed by scope-based narrowing are explicit. |
| 5.10 | Repeated duplicate-merging statement; Ben supplied the phrase to remove. | **Addressed in ZIP, 50:14–16.** Duplicate merging is described once in the collection paragraph, without repeating the old “Duplicate records and” wording in selection. |
| 5.11 | “Classification into what?” | **Addressed in ZIP, 50:16.** Experimental preparation, afferent class, feedback pathway, and mechanism are specified. |
| 5.12 | “I'm lost.” Ben identified tautological source-grounding wording. | **Addressed in ZIP, 50:16.** The tautology has been replaced by the actual fields recorded and the concrete selection workflow. |
| 5.13 | Advisor objected to unnecessary source-grounding detail and used an integrity accusation. Ben replied that the paragraph should use active voice and retain truthful scope/removal rules, including missing ResearchRabbit poster text. | **Ben's requested response is present, 50:16.** Active first-person verbs explain initial import, retention, scope exclusions, flagging, and removal of records whose source text could not be located. Do not respond to the accusation by inventing extra scientific procedures. A separate existing consistency check remains: Appendix E:12 says 38 retained records lack “grounding text”; clarify whether that means missing curation notes versus unavailable sources before claiming the entire cross-chapter account is unambiguous. |
| 5.14 | “What table view?” Ben: call it the Table tab. | **Addressed in ZIP, 50:8 and 24.** The four tab names are introduced, and the Table tab is named consistently. |
| 5.15 | Overview of views needed; Ben prefers “tabs.” | **Addressed in ZIP, 50:8.** Table, Pivot, Bubble map, and Search (online) are listed with their purposes. |
| 5.16 | “huh?” about a smaller maintenance subset; Ben also questioned its point. | **Addressed in ZIP, 50:20.** The unexplained aside is absent from Chapter 5. Appendix E:8–12 retains the actual separately scoped inventories rather than silently deleting their numbers. |
| 5.17 | Sparse classifications; advisor wanted the database called a start. Ben explicitly deferred the larger classification exercise. | **Deferred by Ben.** Current counts and coverage limitations remain visible at 50:22 and 63. No record sweep, reclassification, or extra completion demand is proposed. Moving the future-work paragraph under U30 must preserve an honest coverage caveat. |
| 5.18 | “which?” category. Ben specified Papers. | **Addressed in ZIP, 50:22.** The sentence explicitly begins “Within the Papers category.” |
| 5.19 | “huh?” about Louvain/topic clusters. Ben supplied that sentence as the anchor and wondered about the map's category legend. | **Wording clarified in ZIP, 50:22 and 26.** Citation-based grouping, 18 named topic clusters, and the layout use are explained. The saved evidence does not prove that an additional category-legend callout was requested or resolved. Preserve the existing six callouts; no seventh is invented. |
| 5.20 | Awkward “distributes the records” sentence; Ben supplied its anchor. | **Addressed in ZIP, 50:8.** The obsolete packaging phrase is absent; embedded records, notes, and graph support the concrete offline-use explanation. |
| 5.21 | Is the paragraph about all views or just Table? Use topic sentences. | **Addressed in ZIP, 50:24.** The paragraph opens with the Table tab and describes its filters and detail workflow. |
| 5.22 | Jumping map → table → map. | **Addressed in source order, 50:24–28 and 37.** Table precedes the map workflow, and Pivot has its own Results paragraph. U32a separately proposes placing the relevant screenshots closer to those descriptions; good source order alone does not certify final float placement. |
| 5.23 | Put Table figure in main text; irrelevant papers in a supposedly curated app were confusing. | **Not addressed in ZIP; proposal retained.** Chapter 5:24 still points to Appendix E's Table screenshot. Appendix E:16–21 includes `sadb_explorer.png`, whose visible archived/generic entries include deep-learning and turbulence papers. Capture a relevant retained-paper example with the archive state clear and introduce it early in Chapter 5 under U32a. Do not promote the old screenshot unchanged or reclassify the corpus under this task. |
| 5.24 | Toggle between layouts, or x/y axes? | **Addressed in ZIP, 50:26.** Two selectable layouts are named; the year–citation axes and spring-layout topic landscape are distinguished. |
| 5.25 | “huh?” with uncertain/dangling anchor. Ben suggested the topic-landscape sentence as a possible anchor. | **Possible wording addressed; live anchor unverified.** Chapter 5:26 explains the spring layout concretely, but the saved record cannot establish that this was the actual live anchor. Keep the item explicitly uncertain until the live thread is located. |
| 5.26 | “The user?” Ben supplied the subject correction. | **Addressed in ZIP, 50:26.** “The user controls” now introduces the controls. |
| 5.27 | What does “proposed” mean? | **Addressed in ZIP, 50:24 and 26.** The prose refers to notes on possible applications to a robot or simulation, rather than an undefined proposal. |
| 5.28 | Does two links mean two more or one more? Ben asked for app verification, perhaps a callout. | **Addressed in ZIP, 50:26 and 33.** The text names the 1/2/3-link selector, explains two links as one level beyond direct connections, and points to the existing neighborhood control callout (E). No new callout is required for the same control. |
| 5.29 | “huh?” about incoming citations. Ben asked for cites/is-cited-by clarity. | **Addressed in ZIP, 50:26–28.** Papers citing the focus and papers cited by it are explicitly distinguished; colors and connection-list navigation use those directions. |
| 5.30 | Cite both database presentations at the start. | **Addressed in ZIP, 50:4.** Both the 2023 and 2024 presentation references are present. |
| 5.31 | Begin broadly with database value, then model/simulation/robot design. | **Addressed in ZIP, 50:55–59.** General evidence navigation precedes neural-model design and the more specific sensory extensions. |
| 5.32 | Begin more positively with the working features. | **Addressed in current opening, 50:63; preserve the intent during U30.** The paragraph starts with searchable notes, comparison, and navigation. Current uses already receive substantive positive treatment at 50:55–59, so moving future work need not leave Chapter 5 defined by deficiencies. |
| 5.33 | Discussion reads like Results; explain why the resource matters using declarative uses. | **Addressed in ZIP, 50:55–59.** The text develops literature learning/evidence tracing and neural-model design, rather than listing controls again. |
| 5.34 | Another feature paragraph reads like Results. | **Addressed in ZIP, 50:37–46.** Pivot description, its figure, and pathway inventory are in Results, ahead of the Discussion at 50:52. |

## Completion boundary

The requested SADb findings and every saved Chapter 4/5 card have a current-state disposition and proposed action or explicit exclusion. The remaining concrete SADb proposals are an effective Figure 5.2 crop, an early relevant Table example, consolidating future work in Chapter 7 while preserving the Chapter 5 segue, and accurately describing the already implemented monthly watch. The Table screenshot request is still an actual manuscript miss; the two uncertain anchors and the unobserved 30 October live state remain verification limits. None is relabeled as completed implementation.


---

# Other-audit comparison: method, scoreboard, and appendices

Task 3 follows completion of Tasks 1 and 2. The reviewed document is the frozen copy Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030.md, SHA256: 859e9a7fc0738d92116287601c5157716ea3e8fc0e7cd7351ee484974185f772. The supporting rubric at tmp/audit_20261030/RUBRIC_ADVISOR_PAIN_20261030.md was also read. This evaluates proposals; it does not implement them or accept the other audit's scores as grades.

## Method and global assertions

| ID | Other audit's item | Evaluation and evidence |
|---|---|---|
| O-M01 | Latest ZIP authority;306-page dissertation/33-page paper | **Agree, with scope.** Independent fresh builds produce 306 pages, 33 paper pages, and 10 supplement pages. Source manifests identify the supplied archives. A source build says nothing by itself about newly added live comments. Current source/compiled output takes precedence over historical page counts. |
| O-M02 | Integrity, tracked-change, clarity, repetition, awkwardness, and float penalty classes | **Useful checklist; modify how applied.** These are sensible feedback categories when tied to an actual miss and Ben's reply. An advisor's insulting wording is evidence of a serious editorial failure to learn from, not evidence of scientific misconduct. Apply the largest relevant cost to a unique miss, avoiding multiple deductions for its repetitions. The existing REVIEW_RUBRIC.md evaluates assistant failures rather than grading Ben's chapters. |
| O-M03 | The “machine-learned part” | **Describe as feedback-derived rules.** No trained or validated predictive model is documented by the supplied report/rubric/scripts. Persist the specific corrected rules and counterexamples; do not present manually assigned weights as learned scientific measurements. |
| O-M04 |450-word wall threshold; default DELETE | **Modify.** Use density as a navigation trigger. A word count between display objects cannot establish needless prose, poor reasoning, or a need for a new figure. Inspect the argument and page layout; delete repetition, retain necessary assumptions/evidence, and move unique technical inventories where appropriate. Chapter 2's six consecutive text pages also show why a count threshold can miss reader burden (15-background.tex:236-248). |
| O-M05 | Rubric says Results paragraphs should open with a figure callout | **Reject as a universal rule.** Saved comments4.2/4.8 and Ben's U13 request say important result first, then show. Current Ch. 4:549/606 follows that instruction. Lead with the outcome and attach the relevant reference naturally; do not replace it with an empty “Figure… shows” opening. |
| O-M06 | Define every term at first use | **Modify.** Introduce unfamiliar abbreviations and explain local variables/interfaces. Consider the whole dissertation's earlier definitions and its specialist audience. A brief local reminder may help; defining every ordinary ME/neuroscience term again creates the repetition the rubric penalizes. Detailed term decisions appear in the Chapter 2 and Chapters4-7 responses. |
| O-M07 | Active voice; settings/bounds rather than workflow logs | **Agree with exceptions.** This matches Ben's style and generally improves Appendix C:151-216 and Ch. 7:8-16. Retain scientific chronology when it explains the tested configuration, selection criteria, or a failed comparison; do not erase Ben's requested Zotero import/selection account at Ch. 5:14-16. |
| O-M08 | Every unit below70 needs a pre-defense revision pass; all fixes bounded/under a day | **Do not use as a gate or schedule.** The score has no demonstrated calibration, normalization for unit length, or validated threshold. Use the concrete severity/provenance of each finding, Ben's instructions, and available time. Several proposed fixes need provenance or figure checks, so their duration is not established. |
| O-M09 | No displayed number may be touched; flag-only if it would | **Modify to Ben's actual rule.** No silent experimental/result changes. Ben has explicitly requested display-unit conversion and percentage formatting (U18/U23); those are authorized proposals, with cell-by-cell checks. A blanket numerical freeze would obstruct his instructions. No underlying measurements or optimization outputs are changed here. |
| O-C01 | Zero tracked-change artifacts everywhere | **Narrow the claim.** The audited TeX does not show visible revision markup. A normal Overleaf ZIP does not establish that live tracked changes are absent, reviewed, or accepted. Ben's plain-text Word edits and replies remain material evidence even without tracked markup. See the 42-card crosswalk and annotated-review JSON. |
| O-C02 | Every prose number reconciles; notation clean everywhere | **Not a general clearance.** Matching a table does not verify definitions or provenance. T01 finds the cam-margin definition inconsistency (Ch. 4:709-722/paper:600); T02 finds signed-versus-magnitude torque exposition (paper:176-194); T03/T10 find unlinked coordinate/length definitions. APP-A1 finds contradictory bracket frames (App. A:198/App. C:186,216,226-227). These require existing-record reconciliation. The narrower diameter/model vocabulary checks are useful but not a full equation audit. |
| O-C03 | Exactly one integrity-pattern survivor, all other units clean | **Reject the classification as a clearance or task.** The lexical class is not an assessment of scientific integrity. Its Ch. 5:50 item belongs to the AI-tool comparison that Ben explicitly excluded; it must not re-enter the task through a new label. Lack of this phrase elsewhere does not certify scientific quality. |
| O-C04 | All appendices are referenced from main chapters | **Agree at appendix level.** App. A-C support mechanics, App. D supports Ch. 6, and App. E supports Ch. 5. This does not imply that every figure is introduced or that every subsection needs a separate main-text citation. The specific missing callouts below are actionable. |
| O-C05 | Motor-schematic squares verified absent, open thread effectively clean | **Not sufficient to close a thread.** See Chapters4-7 response: inspecting five Ch. 7 assets cannot establish an unlocated live anchor's identity or resolution. No live review action was taken. |

## Scoreboard: response to every unit

The other audit explicitly calls these pattern-density indices rather than quality grades. Even on that interpretation, their unnormalized deductions and floor effects make the numerical ranking unsuitable for deciding scientific quality or clearance. The useful information is the underlying, verified finding.

| ID / unit / reported score | Merits and limitations; location |
|---|---|
| O-S01 Ch. 4,24 | Keep the targeted prose/notation checks; do not use the index as a priority rule. Existing correction/model/physical evidence is substantial. Prioritize margin definitions and Ben's concrete table/legend/pulley directives, not blanket jargon insertion. Ch. 4:665-788; U17-U28. |
| O-S02 Ch. 5,29 | A useful crop, screenshot, and organization pass remains. The alleged integrity line is excluded, and the import history is explicitly requested by Ben. Ch. 5:14-16,24,39-50,61-67; U29-U32. |
| O-S03 Ch. 6,36 | Repeated tools/capability framing warrants compression, but it supports a tools chapter and should not all disappear. Preserve functional tests and model distinctions. Ch. 6:4,29-39,63-67,137,234; SNS-S07. |
| O-S04 Ch. 7,27 | Long synthesis and detailed artwork/caption checks are worthwhile. Do not equate a901-word multi-paragraph section with a single paragraph or delete synthesis merely because chapters recur. Ch. 7:8-16,20-28,37-47. |
| O-S05 App. A,95 | Pointer-based compression is reasonable; the apparently high score overlooks the more serious frame contradiction already in APP-A1. App. A:198 versus App. C:186,216. |
| O-S06 App. B,83 | c1/V3 can be clarified, and repeated outcomes condensed. Preserve the unique conversion specification/calibration inventory. Correct the independently observed Figure B1 lettering mismatch, which this headline misses. App. B:162-201; APP-A2. |
| O-S07 App. C,65 | Strong agreement on reducing implementation history and duplicate Table C13. Ben's specific C.9.1/C.4 deletions, unit consistency, master's provenance, and graded-grid definition are more useful than a score. App. C:205-227,387-391,429-495; U03-U05/U18-U19/U27-U28. |
| O-S08 App. D,45 | The six missing figure callouts are real, useful new catches. The three allegedly orphaned subsections are relevant within the appendix and need purpose/connection checks, not automatic deletion. App. D locations below. |
| O-S09 App. E,95 | Two explicit table callouts would help. The score overlooks the irrelevant Table screenshot and source-versus-note wording tension identified by our audit/Ben's comments. App. E:12,16-21,25-27; DB-S1/DB-S2/U32a. |
| O-S10 Paper Main,3 | The paper needs a substantive presentation/reconciliation pass, but3/100 is not a meaningful journal-quality measure. Concrete sign/margin/parameter definitions, actual final-size fonts, and Ben's revisions outrank a raw wall count. Main.tex:176-194,581-643; T01-T15/U01-U16. |
| O-S11 Paper supplement,44 | Reconcile distinct coordinate/target/length definitions and reference placement. Distinguish true scientific ambiguities from source-only build hygiene and unsupported float-limit claims. Exact responses in the paper comparison. |
| O-S12 Ch. 2,0(raw-128) | The severe numerical floor hides which claims are scientifically wrong and which are optional local glosses. Our 13 findings distinguish factual/citation errors from readability. The other audit's three blanket CC-BY suggestions and missing-assets assertion need correction; see Chapter 2 comparison. |
| O-S13 Appendix means77/66; paper combined24 | Arithmetic summaries of these subjective, length-dependent scores add no independent evidence. Do not use them for scientific ranking or clearance. Preserve the item-level findings and source anchors instead. |

## Appendices: every named proposal and flag

| ID | Other audit item | Evaluation, source location, and proposed disposition |
|---|---|---|
| O-A01 | App. A/C compliance-array duplication; “pointer-ize” | **Agree selectively.** Keep one clear core mechanics presentation and the unique configuration specification. Ben's U19 explicitly removes the repeated C.4 display; do not delete essential Appendix A constitutive relations solely because arrays recur. App. A:189-210; App. C:205-227; APP-A4/T13. |
| O-A02 | App. B c1 and V3 unglossed | **Agree with a brief local explanation.** App. B:201 says both provide commissural pathways, so they are not wholly unexplained. Identify c1 as the implementation's population label and V3 by the modeled interleg pathway/function; do not invent a biological equivalence or gene-expression claim not supported by the source network. This also responds to the same cross-cutting jargon item. |
| O-A03 | App. B AnimatLab catalogue duplicates App. D | **Modify.** App. B:162-180 gives conversion rules and structural census; App. D:273-292 summarizes rhythm outcomes. These are complementary. Condense the duplicated outcome column/commentary or cross-reference it once, preserving the unique conversion parameters and failed-conversion evidence. APP-A5's removal of acceptance bookkeeping at App. B:194 remains appropriate. |
| O-A04 | App.C583-word wall/process narration at385-396 | **Agree on compression, not wholesale deletion.** The route-builder chronology, rejected-copy detail, and long crash inventory can move to repository notes. Retain supported versus incomplete scope, authoritative results, and the distinct main-chapter/cam-margin issue. App. C:387-391; APP-A4. Count is a heuristic, not a required word target. |
| O-A05 | App. C tab:TorqueFit is verbatim duplicate | **Agree, confirmed.** Table C13 at App. C:433-450 duplicates Ch. 4 Table 4.8:617-632. Remove the appendix duplicate and C.9.1's repeated definitions under U03/U04; keep the main table and pinned-test tables. No metric change. |
| O-A06a | App. D missing callout: phase-1 recordings | **Agree.** fig:animatlab_phase1_preliminary is labeled at95-AppendixD.tex:33 and never explicitly referenced. Attach its reference to the phase-1 explanation at:8 or the associated recording description. Distinguish the earlier instrumented network from the latched reorganization. |
| O-A06b | App. D missing callout: MuJoCo views | **Agree.** fig:mujoco_model_views at:47 has no explicit reference. Introduce the model/frame views in the model-construction explanation before the figure. Do not infer simulated performance from the render. |
| O-A06c | App. D missing callout: knee transient | **Agree.** fig:app_knee_transient at:70 lacks an explicit callout despite substantive results at:60-64. Cite it with the settled pulse/recovery result. No new transient analysis. |
| O-A06d | App. D missing callout: layer stimuli | **Agree.** fig:layer_stimulus_validation at:110 lacks a callout; the text at:79-83 cites only the table. Refer to the figure with the result it demonstrates, preserving non-physiological stimulus scope. |
| O-A06e | App. D missing callout: W2L perturbations | **Agree.** fig:app_w2l_perturbations at:320 lacks an explicit reference. Cite it in the pathway-deletion/stimulus explanation at:312-314. |
| O-A06f | App. D missing callout: force comparison | **Agree.** fig:force-compare at:453 lacks an explicit reference. Add the figure reference to:405-407, distinguishing normalized shape from the table's absolute force statistics. |
| O-A07 | App. D three orphaned sections | **Modify.** The likely intended locations are Archived AnimatLab Models(:115-119), Di Russo Replication(:323-327), and Cross-Platform Force Comparison(:402-407). Their absence of direct main-text section references does not make them irrelevant: they support the existing implementation/transfer account. Add a purposeful connection where needed; retain the failed replication scope, rather than deleting negative evidence. Multiple aliases on an already referenced heading also do not make separate orphan sections. |
| O-A08 | App. D former-startup narration at60 should become final state only | **Modify; preserve the essential comparison.** At:60 the original configuration did not settle;:62-64 identifies the modified six-parameter configuration and recovery result. Removing this distinction could attribute the valid transient result to the wrong configuration. Compress the old-startup aside, retain settled-versus-unsettled scope. This overlaps SNS-S03. |
| O-A09 | App. E tables lack explicit callouts; two-line fix | **Agree.** At 96-AppendixE.tex:25 the generic “tables” introduction should identify tab:sadb_afferents(:29) and tab:sadb_pathways(:46) explicitly. Cite them where their different purposes are explained; no extra table or duplicated inventory is required. |
| O-A10 | Ben flag7: App.A198 versus App.C227 frame line306/305 | **Important, but broader than a line-number typo.** App. C:186 says both pinned-flexor brackets use two rotations;:216/226-227 and App. A:198 say pitch-only. Reconcile saved identification provenance and present implementation, then fix stale line references. Do not pick a convention by majority wording or rerun fitting. This is APP-A1. |

## Suggested pass order and session notes

| ID | Other audit item | Evaluation |
|---|---|---|
| O-P01 | First priority: Ch. 5 integrity sentence | **Reject.** That excluded AI-comparison item must not become the first task through a new label. Use actual unresolved evidence/definition issues and Ben's explicit priorities. |
| O-P02 | Process narration, then~14 jargon glosses | **Modify.** Apply the source-specific accepted changes, retaining Ben's requested import history and valid configuration comparisons. Do not run a blanket gloss/delete transformation. |
| O-P03 | Walls, delete-first | **Modify.** Cut verified duplication and relocate unique inventories. Keep necessary reasoning and comparisons; no automatic word-count quota. |
| O-P04 | App. D callouts and App. C deduplication | **Agree.** These are concrete, low-risk proposals with identified destinations above and U03/U04/U19. |
| O-P05 | Ch. 7 captions require regeneration in MATLAB | **Verify the smallest remedy.** A confirmed wording/panel mismatch requires repair, but caption edits or removal of a baked heading may suffice. Do not assume MATLAB is the source of every schematic. See Chapters4-7 visual checks. |
| O-P06 | Paper pass may be post-defense; float limit and bib are quick | **Modify.** Ben controls timing. There is no established15-display Original Research cap here; bibliography placement is source hygiene in this build. Prioritize genuine definition, legibility, and current-artifact issues for submission. |
| O-N01 | Ch. 6 title synced from live;306-page clean rebuild | **Present title/build verified; historical action not independently certified.**60-sns.tex:1 reads Tools for Neuromechanical Simulation, and the fresh source build is306 pages. This audit made no synchronization or live edit. |
| O-N02 | Ben's figure drags all landed | **Confirm presence, not a universal history claim.** The newer crossed-knee asset is in the supplied dissertation ZIP and included at Ch. 4:780. Figure 4.14's drawing order/alpha is verified by U25. This does not prove every unknown drag or today's live state. |
| O-N03 | Retract the 38-page paper/live divergence flag | **Agree about the live preview; retain package hygiene.** Current source/live recompile is33 pages. The ZIP also contains an actual stale 38-page PDF and five-page supplement PDF, not only a stale log. T04 concerns accidental redistribution of those bundled artifacts, not a claim that live Overleaf currently compiles38 pages. |

An ancillary tmp/audit_20261030/ben_findings_incorporation.md session note maps Ben's “5.3.2” to the paper's hardware subsection. That mapping is incorrect: his finding follows Figure 5.2 and concerns dissertation Chapter 5's future-improvements discussion. U30 supplies the correct50-database.tex:61-67 to70-overall.tex:59-61 proposal. The named consolidated audit, its supporting rubric, and all its listed proposals are the scope of the item-by-item comparison; unrelated agent-session outputs are not silently treated as additional advisor directives.


---

# Response to the other audit: Chapters 4–7

Task 3, after the Chapter 2 audit and incorporation of Ben's findings. Audit/proposed-change ledger only; no manuscript or figure edits.

The reviewed document is the frozen `Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030.md`, not a subsequently expanded version of the original audit. Current source evidence is the dissertation ZIP snapshot in `Archive/audit_20261030/dissertation/chapters`. Below, **30** means Chapter 4, `30-torque.tex`; **50** means Chapter 5, `50-database.tex`; **60** means Chapter 6, `60-sns.tex`; and **70** means Chapter 7, `70-overall.tex`. All lines are one-based snapshot source lines. Figures refer to the assets actually included by that snapshot; five Chapter 7 assets were rendered afresh for this comparison.

**Verdicts:** accept = a supported proposed change; qualify = a useful concern requiring the narrower action stated; reject = unsupported or contrary to Ben's ruling; defer = explicitly left to Ben or the separately assigned response. None means a change was implemented or a live comment was resolved. The complete 42-card source/comment/reply crosswalk is in `USER_FINDINGS_SADB_AND_COMMENT_CROSSWALK.md`.

## Class 1: the proposed integrity-pattern priority

| Item | Other audit's proposal | Verdict and evidence |
|---|---|---|
| 1.1 | Treat 50:50, “I retained responsibility for inclusion, classification, and deletion decisions after checking the source,” as the single surviving integrity hazard; fix first. | **Reject as a task and penalty.** This is inside the AI-tool comparison that Ben expressly excluded after explaining that the ratings were his judgment of performance relative to his own. No new assessment method, comparison audit, or integrity-driven wording task is assigned. The actual advisor integrity accusation was comment **5.13**, whose confirmed anchor is the curation paragraph now at **50:16**, not the AI paragraph at 50:50. Ben's active-voice, scope-based retention/removal reply is present there. Conflating these two passages disregards the captured anchor and reply. Do not make this the first revision priority. |

## Class 2: all four process-narration claims

| Item | Site and proposed removal | Verdict and narrower disposition |
|---|---|---|
| 2.1 | 50:16, initial broad Zotero import followed by narrowing; replace with final criteria only. | **Reject blanket removal; qualify possible shortening.** Advisor 5.9 asked why records were imported then removed. Ben's reply explained bulk collection from AARL/personal Zotero and keyword searches followed by trimming. Advisor 5.13's later reply explicitly retained the scope/removal rules. **50:14–16 now answers those questions in active voice.** Keep the causal import → selection account and the archived/working-base distinction. Condensing duplicate words is permissible; removing the explanation would recreate the documented confusion. |
| 2.2 | 60:145, narrow versus broadened synaptic saturation ramps; state only final settings. | **Qualify.** The brief original-settings history can leave the main Results paragraph, but the changed ramp configuration is necessary to identify the sensory-current experiment and distinguish it from the earlier startup comparison. Retain the tested configuration and the fact that this is not a physiologically calibrated tendon-tap model. Detailed comparison belongs at **95-AppendixD.tex:60–64**, already pointed to by 60:145. Do not erase a consequential configuration change to make a result look cleaner. |
| 2.3 | Appendix D:60, former startup comparison; replace with final state. | **Qualify; appendix response owns final disposition.** This is the provenance that prevents the earlier nonsettled graded/spiking startup comparison from being mistaken for a replication of the later settled transient. The distinction matters at **60:145**, **D:60–64**, and **D:246–269**. A concise explanation is preferable to either an extended development story or silent conflation. See the separate Appendix D response and `SNS_DISCUSSION_AUDIT.md:S03`. |
| 2.4 | 30:342–348, placement-search lineage, called a spine-anatomy non sequitur; replace with final method. | **Qualify; compression is merited.** The lumbosacral pose grid, Cartesian perturbation details, and historical numeric weights interrupt a knee-design Methods section. Keep a short, cited explanation of how prior placement searches motivate the corrected knee-route method; move historical algorithm detail to the relevant appendix if retained. Do not erase Bolen's master's-project contribution or Morrow's lineage: Ben separately asked that Appendix C's routing narrative identify the master's project. **30:342–348** supplies the historical methods; **30:350 onward** develops the knee objective. A complete deletion of the lineage is not warranted. |

## Class 3: every listed terminology claim

These are local reading issues, not evidence of deficient scientific ability. A term defined in an earlier chapter is not automatically “never introduced.” The useful action is a brief reminder where the local context requires one.

| Item | Listed term and site | Verdict and evidence |
|---|---|---|
| 3.1 | **W2L**, 60:4 | **Accept moving the expansion to first local use.** “Walker Two-Layer (W2L)” is already defined at **60:12**. Move that expansion to line 4 rather than adding another explanatory paragraph. The historical model identity remains distinct from Gait2392 and the CAD-derived model. |
| 3.2 | **FVU**, 30:314 | **Reject further action under Ben's ruling.** The metric has equations in **30:320–326** and a textual expansion in the table caption at **30:486**. Its interpretation, selection, and recalculation remain closed. A terminology audit must not reopen advisor 4.5 through a different label. |
| 3.3 | **`l_620`**, 30:83 | **Qualify: a local reminder is useful; a dissertation-wide definition is not missing.** Chapter 3 defines it at **20-muscles.tex:25–28** as unloaded length at 620 kPa. A short parenthetical or cross-reference at 30:83 would aid chapter-level reading. Leave the protected published chapter unchanged. |
| 3.4 | **Hunt model**, 30:470; alleged uncited model | **Accept local identification/citation; qualify the absolute claim.** Chapter 4 uses the label at **30:470, 482, 486** without a local Hunt citation. The dissertation does cite `hunt_modeling_2017` in Chapter 3 at **20-muscles.tex:78**. Identify the comparison as Hunt's BPA force model and cite the existing appropriate record where first used; do not claim that Hunt is uncited everywhere. |
| 3.5 | **Mono/Bi-pulley algorithms**, 30:749 | **Accept plain-language replacement or definition.** The unexplained implementation label precedes the physical explanation at **30:753–761**. Name the tendon-routing variants in mechanical terms, or define exactly what the two algorithm names distinguish. Do not assume “mono/bi” means the number of BPAs or supporting cable parts; the source separately defines the latter as `n`. |
| 3.6 | **Insertion-side measure**, 30:749; apparent 85.56 versus 85.74 N·m discrepancy | **Qualify and retain as a definition/provenance check.** The source already attempts a gloss: torque applied at insertion after tackle force loss. At **30:786**, 85.74 N·m is explicitly knee-axis torque, with 0.803 N·m off-axis torque. These labels do not yet make the coordinate/projection relationship obvious. Specify which moment, point, axis, and optimized configuration each number represents, or omit the secondary summary value if redundant. The text alone does not establish a numerical contradiction, and changing either number is unsupported. No new optimization is required to trace the stored outputs. |
| 3.7 | **Production configuration**, 60:196 | **Accept a precise model/configuration name.** This is an implementation label in a scientific result. Use the identified tested configuration and a cross-reference to its appendix record, consistent with the controller families introduced at **60:85**. Do not casually substitute a neuron count without verifying that the named record is the one used for this result. |
| 3.8 | **Recovered**, 60:209 and caption 60:214 | **Accept clarification or deletion of this adjective.** The paragraph concerns a saved W2L controller replay with all effective gains delivered and fixed comparison conditions. “Recovered” does not identify what was recovered. “Saved W2L controller” plus the existing appendix reference is clearer, subject to preserving the specific replay identity. |
| 3.9 | **Isobaric**, 70:37 | **Optional small improvement, not a major defect.** “Constant-pressure” is a compact gloss useful alongside the proposed dynamic characterization modes. This is conventional mechanical-engineering terminology; it does not justify a large penalty or a new theory section. |
| 3.10 | **Pivot views**, 70:14 | **Qualify: naming consistency, not missing first use.** Chapter 5 already explains the **Pivot tab** at **50:8 and 37**. Use the same tab name in the Chapter 7 recap; do not explain the full interface again. |
| 3.11 | **c1**, Appendix B:201 | **Accept a precise legend-level identity if not already supplied; appendix response owns verification.** **93-AppendixB.tex:201** identifies c1 as commissural but does not say what the label denotes biologically or in the implemented circuit. Distinguish a model population label from an established physiological cell class. No neural redesign is implied. |
| 3.12 | **V3**, Appendix B:201 | **Accept a concise functional identity; appendix response owns final scope.** The same line places V3 in commissural pathways. A brief description tied to the actual modeled population/citation is more useful than acronym expansion alone. Do not assign a physiological role beyond the cited circuit. |
| 3.13 | **OpenAlex**, 50:8 | **Optional local gloss.** Its role is already concrete: the Search (online) tab retrieves candidate additions and requires the internet. Adding “bibliographic index” before the name could help a new reader. The advisor's offline/online confusion is already addressed at **50:8 and 24**; no new comparison or assessment task follows. |
| 3.14 | **Attached-PDF indicators**, 50:22 | **Accept a plain-language field definition.** Explain that 318 is the number of records flagged as having an attached PDF, preserving the count and avoiding an unsupported claim that all PDFs were available or read. This is a bibliographic metadata flag, not an assessment of source quality. |
| 3.15 | **`rho` collision**, paper supplement:242 | **Separate paper/supplement response.** This is outside Chapters 4–7 and must be evaluated against both definitions in that source. Retain it in the master response rather than dropping it from the other audit's list. It does not authorize changing dissertation notation globally. |

## Class 4: walls and the delete-first prescription

| Item | Other audit's prescription | Verdict and evidence |
|---|---|---|
| 4.0 | Default remedy for every >450-word run: delete. | **Reject as a universal rule; accept it as a useful screening prompt.** The reported “walls” aggregate several paragraphs between headings/floats, not necessarily one paragraph or one rendered page. Useful Methods detail is not dispensable merely because a counter crosses 450. Judge repeated claims, local reading load, actual float placement, and whether a table conveys the existing evidence better. Adding figures solely to break a count is equally unhelpful. |
| 4.1 | Chapter 7 §7.1, reported 901 words: split/trim the chapter recap. | **Accept targeted trimming.** **70:8–20** repeats four correction terms, prior placement history, SADb field counts, and tools before presenting Figure 7.1. Keep a synthesis of how the contributions connect and what they enable; cut detailed inventories already supplied in Chapters 4–6. Long single paragraphs at **70:8 and 12** mix result, method recap, and future speculation. The independent audit also identifies the sensory-extension paragraph at **70:75**, beyond this particular wall count. No extra figure is needed to cure duplicated prose. |
| 4.2 | Chapter 5 Results opener, 50:20–28, reported 506 words: delete first. | **Qualify; reorganize and condense.** Counts at **50:20–22** are needed to establish the delivered artifact and coverage. Workflow paragraphs at **50:24–28** are useful but can be shorter when a current Table example appears early. Consolidate repeated inventories with Appendix E and keep the headline denominators, archive distinction, and citation-direction explanation. This aligns with U29/U32a and the verified five-page delay before the first app visual; it does not call for deleting the application's actual contribution. |
| 4.3 | Chapter 6 §6.4.1, 60:71–93, reported 501 words: delete first. | **Qualify; move low-level conversion detail, preserve the implemented coupling.** **60:73–79** contains conversion-failure/repair detail already supported in Appendix D; a short methods summary and reference would suffice. **60:81** explains the BPA-force/sensory-feedback loop and its quasi-static scope, which readers need. **60:85** distinguishes mechanical models and controller families and would benefit from the compact configuration map proposed in the independent audit. The subsection ends at **60:87**; a quoted run through line 93 also reaches the next subsection and figure boundary. The density claim is useful, but its stated interval is not one unbroken subsection paragraph. |
| 4.4 | Appendix C:385–396, reported 583 words | **Separate appendix response, with supported overlap.** Ben already identified §C.9.1 as repeated Chapter 4 material and §C.10 as jargonistic/missing his master's project. Consolidate repeated material and preserve the distinct routing provenance; do not treat every appendix explanation as disposable. |
| 4.5 | Four paper walls | **Separate paper response.** The numerical lengths alone do not establish a delete quantity. The paper response must address each source interval and its claim-first/result/table opportunity. |
| 4.6 | Chapter 4 has zero walls, so compression held. | **Accept only as a result of that screening rule, not a clean prose bill.** Historical detail remains at **30:342–348**, and the reverse-pulley discussion at **30:749–786** can still be reduced under Ben's specific directive. A chapter can pass a 450-word screen and still need focused compression. |

## Class 5: float, callout, caption, and panel claims

| Item | Other audit's claim/proposal | Verdict and verified source |
|---|---|---|
| 5.1 | Six Appendix D figures have no callouts. | **Separate appendix response.** This requires label-by-label searches in narrative prose, distinguishing a caption's own label from a genuine callout. It cannot be certified by counting figure environments. Preserve the item for the root response. |
| 5.2 | 30:476 lacks the opening figure callout. | **Reject “missing” as an absolute claim; optional tightening only.** The preceding paragraph at **30:474** explicitly cites `fig:FlxPin_preopt` and explains the comparison. Line 476 then supplies the result and table reference. Adding the same figure reference to the claim-first sentence at 476 may improve locality if the two paragraphs merge, but no uncited figure is demonstrated. Follow advisor 4.2/4.8's claim-first principle without creating repetitive callouts. |
| 5.3 | Table 4.11 mixes +5.02% and 0.0501 fractions. | **Accept presentation repair; correct the table number.** The current ZIP/build labels `tab:redesigned_route_capacity` as **Table 4.10**, at **30:671–682**, printed p.109. The minimum and required margins mix percentages/fractions and the flexor worst-case angle is missing. Ben already instructed percentages, no plus sign, a sourced worst-case angle, and a clear minimum/required distinction. This is a formatting/definition problem, not a newly demonstrated result error. Preserve traceability when expressing 0.0501 as 5.01%. |
| 5.4 | `design_workflow.pdf` has a baked-in caption atop its LaTeX caption; regenerate. | **Reject for the current asset.** Fresh rendering and text inspection of the ZIP's `figs/70-overall/design_workflow.pdf` show diagram labels only, with no embedded caption paragraph. **70:24–25** contains the normal inclusion and LaTeX caption. The more substantive issue is separate: **70:20** describes two acceptance checks, but the graphic has only the initial “Meets targets?” branch and no pass/fail return branch after closed-loop verification (`SNS_DISCUSSION_AUDIT.md:S04`). Fix that mismatch if revising the figure; do not regenerate it for a nonexistent second caption. |
| 5.5 | `pipeline.pdf` has a baked caption; MATLAB regeneration required. | **Qualify.** The current asset at **70:90** includes an explanatory footer about computer-in-the-loop testing, while **70:91** supplies the figure caption. That footer can reasonably move into the caption/body to reduce duplicate explanation; it is not another numbered caption. A generator/vector edit is sufficient if chosen. The evidence does not establish that MATLAB is required, nor that the diagram must be reconstructed. |
| 5.6 | 70:47 reverses/contradicts the BPA panels' middle-rubber separation. | **Reject the asserted reversal; qualify the wording.** **70:45–47** includes A=`bpa_construction.pdf`, a layered FEA diagram explicitly labeling “Middle rubber,” then B=`bpa_exploded_cutaway.pdf`, the illustrative cutaway without that label. A and B correspond to their caption descriptions. The phrase “separated only in A” can be clearer about the specifically depicted layer, but the assets do not support simply reversing A/B. If edited, describe what each panel actually depicts; the layer's physical validity would require its construction source, not an inference from the word “exploded.” No numerical or model change is proposed. |

The five Chapter 7 asset checks used freshly generated PNGs in `tmp/audit_20261030/other_response_ch7_figures`: design workflow, hardware pipeline, layered construction, exploded cutaway, and sensory hierarchy. They establish the specific visual observations above. They do not establish the identity or resolution of an unlocated live advisor thread.

## Class 6: every repetition proposal

| Item | Repetition claim | Verdict and narrower action |
|---|---|---|
| 6.1a | Tools list occurrence 1, 60:4 | **Keep a concise overview.** This is a short motivation/contribution summary, not the detailed inventory. Do not remove the chapter's opening purpose merely because later Methods names the tools. |
| 6.1b | Tools list occurrence 2, 60:39 | **Accept compression.** This is the strongest duplication of **Table 6.1 at 60:43–61**. State the contribution/capability and point to the table instead of reciting its rows. |
| 6.1c | Tools list occurrence 3, Table 6.1 | **Keep the table.** It is the clearest compact model/tool/use mapping, introduced at **60:40**. Ben specifically restored it in prior work. It can replace prose inventories; it should not be deleted because those inventories exist. |
| 6.1d | Tools list occurrence 4, 60:63 | **Qualify; retain attribution.** This paragraph distinguishes third-party packages from Ben's model implementations, custom blocks, integrations, and procedures. That function differs from Table 6.1's uses. Keep it brief and cited rather than counting it as an identical fourth list. |
| 6.2a | Capability taxonomy occurrence 1, 60:67 | **Qualify.** A Methods orientation and Figure 6.2 workflow callout are useful. Shorten generic capability wording if it duplicates 60:137, while retaining the distinction between model construction, component/interface checks, and experiments. |
| 6.2b | Capability taxonomy occurrence 2, 60:137 | **Keep one clear definition of verification questions.** This introduces the Results structure and appendix scope. It is a better home than repeatedly warning in several later paragraphs that movement is not validation. |
| 6.2c | Capability taxonomy occurrence 3, 60:234 | **Accept deletion or integration into the Discussion.** The generic reminder repeats the distinctions after the individual experiments have already established them. Preserve experiment-specific assistance and comparison limits where they affect interpretation; remove repeated general warnings. |
| 6.3 | FVU values repeated at 30:606, 663 and the table | **Excluded as an FVU-directed action.** The numerical repetition exists, but Ben has closed the FVU item and this pass assigns no further FVU explanation, comparison, or calculation. General prose shortening can later cite the existing fit table without reciting every metric, but this audit does not reopen that closed item or change any results. |
| 6.4 | Paper abstract-to-introduction contribution restatement | **Separate paper response.** Some overlap is normal because the abstract is standalone; remove only redundant detailed claims or repeated explanation, not the introduction's contribution statement wholesale. |
| 6.5 | Appendix C `tab:TorqueFit` verbatim duplicate of Chapter 4 | **Separate appendix response; merits cross-reference-based consolidation.** Ben's redundant Table C13/§C.9.1 findings align with the concern. Verify actual values, row scope, and purpose before calling tables verbatim identical. Keep unique reproducibility constants and identify the master copy. |

## Ben-decision flags 1–4

| Flag | Other audit's proposed decision | Verdict |
|---|---|---|
| 1 | Harmonize percent/fraction columns in the capacity table; requires a Ben decision. | **Already directed by Ben; no renewed permission needed for a proposal.** The affected current table is **4.10**, **30:671–682**. Retain his percentages/no-plus-sign/missing-angle/definition requirements in the ledger. Converting units is a visible formatting change, not permission to invent a new margin or angle. Implementation remains outside this audit's scope. |
| 2 | Tape-measure resolution/uncertainty, 30:66; keep as a defense-Q&A flag. | **Respect the rejection; no task.** Ben rejected the added uncertainty exercise. The measurement method is stated at **30:66**. Do not keep resurfacing a declined study as a prerequisite under “defense Q&A.” No measurement campaign, bound, or speculative uncertainty value is proposed. |
| 3 | The 943-entry Papers denominator at 50:22 is redundant but a number of record. | **Accept denominator preservation and optional consolidation; no special decision needed.** **50:8** describes 1,393 total records, while **50:20–22** distinguishes the 943 Papers category and its field coverage. These denominators serve different purposes. Repeating “Within the Papers category” twice in line 22 can be tightened; retain its scope and every surviving reported count. This is not a classification-completion task. |
| 4 | Verify 70:47 against vector art; A/B may be reversed. | **Verification completed for the audit; reversal not supported.** See item 5.6. A is the labeled layered FEA panel and B the illustrative cutaway. Clarify the layer-specific wording if editing later; do not swap panels to satisfy an unproven diagnosis. |

## Motor-square clean claim and suggested priorities

The claim that the open motor-schematic square issue is “verified absent from Chapter 7's five figures” **does not certify that advisor item**. The current five assets render without an obvious missing-glyph square, and their intentional diagram boxes are not defects. However, the historical item concerns square positioning in a motor schematic; no reliable current live anchor tying it to these five assets was recovered. **Verdict: unverified thread, not verified clean or resolved.** Asset rendering and comment identity are different checks.

The other audit's strongest Chapter 4–7 recommendations are local naming clarity, consolidation of the historical placement block, concise results-first prose, consistent margin presentation, and removal of duplicated inventories. Its proposed first priority (the AI-comparison integrity sentence), blanket removal of the Zotero import/selection explanation, mandatory regeneration of two allegedly double-captioned figures, and presumed panel reversal are not supported under the current sources and Ben's replies. Its wall counts help locate dense passages but do not prescribe a fixed deletion quota or establish doctoral competence.

Suggested pass order for this scope is therefore: retain Ben's concrete figure/table directives and the genuinely unaddressed Table screenshot request; clarify model/configuration/torque-measure names; condense repeated inventories and Chapter 7 recaps; correct verified figure/prose mismatches. Keep the closed/deferred items closed. Appendix, paper, Chapter 2, and rubric/scoreboard responses are supplied separately in the master comparison, including the outside-scope terms and repetition items explicitly preserved above.


---

# Response to the other audit: torque paper and supplementary material

30 October 2026. Findings and proposed-change ledger only. No manuscript, figure, citation database, model input, result, or optimization was changed.

## Scope and source key

This responds to every paper/supplement claim or proposal in the frozen `Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030.md` (**O**). O line numbers below identify the assessed item. Broader dissertation/rubric/appendix findings are answered separately by the lead reviewer.

- **P**: `Archive/audit_20261030/paper/Main.tex`.
- **S**: `Archive/audit_20261030/paper/frontiers_SupplementaryMaterial.tex`.
- **B**: `Archive/audit_20261030/paper/BolenFrontiers22.bib`.
- **C**: `Archive/audit_20261030/dissertation/chapters/94-AppendixC.tex`.
- **Fresh builds**: `Archive/audit_20261030/paper_build/Main.pdf` (33 pages) and `frontiers_SupplementaryMaterial.pdf` (10 pages). Paper/supplement page numbers below refer to these outputs.
- **Our earlier findings**: T01–T15 in `TORQUE_AUDIT.md`; Ben's U IDs in `USER_FINDINGS_PAPER.md` and the combined user ledger.

Statuses: **Agree** = verified and useful; **Modify** = useful concern requiring a narrower interpretation or remedy; **Reject** = not supported as stated; **Verify** = a specific existing record must settle the question before an edit. These statuses are evaluations, not executed revisions.

The most consequential corrections are that the alleged general Frontiers 15-display limit is not established by current Original Research guidance, and the fresh builds do not reproduce the alleged duplicate-label warning. Citation order, missing specific supplement calls, repeated prose, and the rho symbol collision are real. Neither matching word counts nor a successful compile clears the paper for submission.

## Responses

### OP01 — Paper scores, combined score, and “one revision pass” estimate

**O:28–31, 11–13. Modify.** Anchors: P §§1/2.8/4.3.1/4.4 at P:90–103, 376–383, 592–603, 630–644; S:194 and 242; fresh builds pp1–3/14–15/24–28 and supplement pp4–7.

The low scores point toward genuine compression/organization work, but 3/100, 44/100, and the combined 24 are not independent measures of journal readiness, scientific correctness, or the author's ability. This review does not adopt or recompute those scores. The wall counts can be reproduced, whereas the adequacy of a single revision pass and the relative penalty calibration are not established by that reproduction. The six “submission blockers” also need the item-by-item corrections below. Preserve the actionable observations, not a grade-like interpretation of their sum. See T01–T15 and OP18 for additional unresolved issues.

### OP02 — Four main-paper walls; zero supplement walls; delete-first remedy

**O:28–29, 59–61, 100. Agree with the counts; modify the inference/remedy.** The supplied `tmp/audit_20261030/scripts/wall_measure_paper_20261030.py`, run read-only against the frozen sources, reproduces exactly 884, 529, 770, and 548 words. It counts prose between source-level section/display/float boundaries. It does not measure one paragraph, a rendered page, a universal journal limit, or scientific redundancy.

| Source and section | Reproduced count / fresh pages | Assessment and proposed response |
|---|---|---|
| P:90–103, §1 Introduction | 884 / pp1–3 | Several paragraphs form one uninterrupted introductory section. An introduction need not contain an illustrative figure. Shorten the repeated contribution enumeration and anthropomorphic wrapping phrase; add the omitted 27-site scope. Do not cut the question/literature gap merely to force fewer than 450 words. U06, T12. |
| P:376–383, §2.8 Identification Procedure | 529 / pp14–15 | The long stage/settings paragraph at P:382 is a stronger compression target than the count alone. Replace repeated settings/test-membership lists with one compact stage table or a concise supplementary specification. Retain training/validation distinctions and the shared stiffness constraint. T09; no FVU redefinition or new optimization. |
| P:592–603, §4.3.1 Route Redesign and Cost Functions | 770 / pp24–25 | Strong agreement: repeated setup, six release angles, coordinates, margins, and route mechanics obscure the outcome. Use a short result-first comparison and a compact capacity table; point to supplementary geometry. Resolve T01's definitions from saved records before finalizing the table. U12, T08/T09. |
| P:630–644, §4.4 Whole-Leg Capability Against Human Torque Surfaces | 548 / pp26–28 | Some length is necessary to define denominators, inherited routing, actuator substitutions, and the 27-site result. Introduce the relocated pulley schematic here and consolidate repeated setup/coverage explanations. Add the actual existing selection rule, rather than deleting the definitions that make coverage interpretable. U08/U09/U14/U16. |

**Supplement S:97–254, §§1.2/1.4, pp2–7:** no run exceeds the audit's 450-word threshold; the contact-release method still needs readable definitions and layout. The rho collision at S:242, crowded route labels, and the nearly empty p7 are unaffected by a “zero walls” count. Preserve the relevant equations instead of treating a zero count as proof that no editing is needed. See T07 and OP07.

### OP03 — “Zero tracked-change artifacts” and clean integrity-pattern classification

**O:35, 38–39. Agree only within the inspected static-source scope.** Anchors: P:62–747 and S:36–312, the complete document bodies; P:81/101/678–680 for the central claim framing.

The present TeX sources contain no ordinary `added/deleted/replaced`, `DIFadd/DIFdel`, or strikeout/highlight macros in their active content. No passage reviewed here warrants an allegation about Ben's scientific integrity. Static ZIP source cannot certify the separate live Overleaf review database or prove acceptance/resolution of comments and tracked changes. Do not turn the absence of such TeX macros into a statement that all live review state is clean. No tracked changes should be accepted as an audit “fix.”

### OP04 — “Every number reconciles” and “no notation violations”

**O:36–37. Modify; do not treat as numerical clearance.** Anchors: P:386 and 600 (§§2.9/4.3.1, pp15/25); P:176–194 (§2.5, p7); S:69–89 and 121–140 (Tables S2/S4, pp1–2); S:242 (§1.4, p6).

Prose-to-table transcription checks are useful but do not settle target/domain definitions. T01's two margin percentages remain inconsistent under the manuscript's stated ordinary margin interpretation; T02's nonnegative arm/force product is introduced as a signed torque; T03's two biomimetic-flexor coordinate inventories remain unmapped. The rho reuse below also narrows the broad “no notation violations” claim. These are definition/provenance questions, not permission to alter results. Keep Ben's existing number-change restrictions and verify saved records; no reoptimization or FVU recalculation follows.

### OP05 — Abstract-to-introduction repetition

**O:67–69. Agree with targeted compression.** Anchors: abstract P:81, p1; final introduction paragraph P:101, p3; Methods overview P:108, p3; conclusion P:678–680, pp29–30.

The near-parallel technical contribution lists repeat the same offset/stiffness/wrapping description. Normal overlap between abstract, introduction, and conclusion is appropriate; identical long enumerations are not needed. Keep the abstract outcome-focused, let the introduction end with the study question/scope and concise contribution, and let Methods explain the sequence. Incorporate both physical knee tests and the 27-site modeled comparison as Ben requests, with their distinct evidence scopes. Do not remove the central contribution to avoid all repetition. This supports T12/U06 and does not replace Ben's requested one-page dissertation abstract rewrite, which is a separate artifact.

### OP06 — Citation/jargon implications of the cross-cutting list

**O:53–58. Modify the paper application.** Anchors: P:97 (§1, p2), P:116–124 (§2.2, p3), P:356–366 (§2.8, pp13–14), S:242 (§1.4, p6).

The explicit paper/supplement defect in that list is rho, addressed next. The listed uncited Hunt-model complaint is a dissertation location; it should not be copied into the paper ledger: P:97 already cites `hunt_modeling_2017`. The paper defines `l_620` at P:116–120 and defines its fit metrics at P:356–366. Further plain-language compression can help, but a blanket missing-definition claim for those paper quantities would be wrong. Ben has explicitly excluded reopening FVU calculations; retain that ruling.

### OP07 — The supplement's rho collision

**O:58. Agree.** Anchor: S:242, §1.4 Contact and Release Angles of the Wrapped Route, p6; nearby geometry equations S:237–242.

The same paragraph uses `rho(theta)` as radius of curvature and then `rho = 3 mm` as a clearance relaxation. These are different quantities with the same symbol. Rename the relaxation to a distinct defined symbol, such as `delta_clear`, throughout its explanation/equations, retaining the existing 3-mm value and units. This is a notation edit, not a new clearance assumption or calculation. No separate prior T ID; it supplements the contact-release reproducibility review.

### OP08 — Ben-decision flag 5: vastus intermedius versus medialis

**O:80–81. Verify; a real labeling-purpose question, not a proven wrong muscle.** Anchors: S:69–89, especially S:82–83, Table S2, p1; P:386 (§2.9, p15), P:594/606/623 (§4.3, pp24–27).

Table S2 lists vastus intermedius origin/insertion among “different models”; the redesign and physical comparisons name vastus medialis. Both can be intentional if the table documents a historical baseline rather than the selected redesign target. Recover that table's role from existing model/figure records and either label the historical scope or replace it with the correct current-reference inventory. Do not substitute medialis for intermedius merely because the names differ. The broader conflict between old and current coordinate inventories remains T03. Reconcile purpose before adding missing Table S2 callouts (OP14).

### OP09 — Ben-decision flag 6: 25-mm route tendon versus 5-mm effective tendon

**O:82–83. Modify: the recorded conversion already exists, but the paper omits it.** Anchors: S:145–160, Table S5, p3; P:614 (§4.3.2, p26); C:151 (§C.1, after the physical-dimensions table).

Appendix C explicitly gives the existing relation: the 25-mm route parameter minus a 20-mm fitting adjustment yields the 5-mm effective prediction input. Main P:614 also distinguishes installed physical tendon lengths of 28/30 mm from prediction inputs. Thus these numbers need not be a scientific contradiction. Add the short physical/route/effective distinction and its existing conversion to the paper/supplement. Preserve all three categories and the measured dimensions. This narrows T10 to a disclosure/terminology repair; it does not require recalculating the plot or refitting.

### OP10 — Submission blocker 1: 17 displays exceed a Frontiers limit of 15

**O:88–89, 102–103. Reject as an established Original Research blocker.** Anchors: P:150–660 (13 figure and four table environments); template comment P:729; document/template identification P:9/14.

The count **13 + 4 = 17 is correct**. The claimed current general limit is not established. The official [Robotics and AI article-type guidance](https://www.frontiersin.org/journals/robotics-and-ai/for-authors/article-types) and [Robot Design guidance](https://www.frontiersin.org/journals/robotics-and-ai/sections/robot-design/for-authors/article-types) specify 12,000 words for Original Research without the alleged cap; the journal page's 15-display limit is under **Registered Reports**, with another 15 allowance under FAIR² Data. These are different article types.

P:729 is an inherited template comment, not sufficient authority for forcing two displays out of the paper. Do not merge/promote useful evidence solely to reach 15. Consolidate figures only for demonstrated relevance/legibility reasons, and apply the actual selected journal/section/type rules at submission. This confirms the primary-guidance check in our T07-related submission review. Ben's requested schematic replacement/move and figure corrections still stand independently of display count.

### OP11 — Submission blocker 2: bibliography sits in the middle of the supplement

**O:90. Modify: source/build hygiene is verified; a visible mid-document bibliography is not.** Anchors: S:194 between the §1.3 figures and §1.4 method; S:10–32 package/preamble; fresh `frontiers_SupplementaryMaterial.aux:23`; fresh `.blg:1–7`; fresh `.bbl` is empty.

The command is indeed mid-source. However, the supplement has no active citation commands and no bibliography style. Its BibTeX log reports both omissions and zero entries; the fresh 10-page PDF prints no reference list there. Moving the command to the end alone would not cure this unused BibTeX invocation. If the supplement continues to have no literature citations, remove the unused command and skip its BibTeX step. If references are added, supply the appropriate style and place the resulting list at the end. Keep cross-document figure/table references working. This is a small cleanup, not the claimed visible placement defect.

### OP12 — Submission blocker 3: multiply-defined slippage label from circular external documents

**O:91–92. Reject as a reproduced current-build defect; the requested verification is complete.** Anchors: P:41 and S:17 for reciprocal `externaldocument`; S:191 for the label; fresh `frontiers_SupplementaryMaterial.aux:33`; final `Main.log` and `frontiers_SupplementaryMaterial.log`.

The final fresh logs contain **no multiply-defined-label or undefined-reference warning**. The slippage label is defined once in the supplement aux and not locally redefined in Main's aux. Reciprocal external-document imports alone do not prove a collision. Retain the verified build chain/logs, and use prefixes only if a real namespace collision is later demonstrated. Ben's U02 photo-removal proposal would remove this label/callout as a relevance decision; do not present that removal as repair of a warning this build does not show. Other font/class/import warnings exist, so this finding is not a claim of zero warnings overall.

### OP13 — Submission blocker 4: Figure 11 precedes Figure 8 in first citations

**O:93. Agree.** Anchors: P:410 (Results opening, p16) first calls Figure 11B; P:414/416 first call Figures 8/9; P:442 first calls Figure 10. Main figure first-reference sequence is 1,2,3,4,5,6,7,11,8,9,10,12,13.

The [official figure/table guidance](https://www.frontiersin.org/guidelines/author-guidelines) requires figures and tables to be mentioned in numerical order. U11's removal of the biomimetic-extensor hybrid qualification at P:410 naturally removes this premature callout; verify the remaining sequence after that change and after the U08/U14 pulley move. This is not a reason to renumber figures manually or retain the unwanted hybrid series.

**Additional check, not credited to O:** Main Table 3 is first called at P:416 **before Table 1 later in that same paragraph** and before Table 2 at P:510. P:442 repeats that early Table 3 reference; its natural introduction is P:554. Defer/remove the unnecessary early table calls while retaining the stated parameter values, then verify table order 1–4. Also, Supplementary Figure S2 is first called at P:546 before Figure S1 at P:594; Ben's U02 removal resolves that sequence without adding a new callout.

### OP14 — Submission blocker 5: five supplementary displays lack direct main-text references

**O:94. Agree with the exact count; modify “uncited” and the blanket-blocker interpretation.** Anchors: below; inventory P:708, p30.

| Supplement item and exact source | Main-text callout | Specific internal supplement call | Proposed response |
|---|---|---|---|
| Table S2, `table:muscle_attachment`, S:68–92 (§1.1), p1 | None | S:98 | Reconcile historical/current coordinates and target identity first (T03/OP08); then cite the authoritative inventory from Methods. Do not add a callout that makes two unresolved coordinate sets appear equally authoritative. |
| Table S3, `tab:RoutePinnedBrackets`, S:100–116 (§1.2), p2 | None | None | Cite it where bracket reference frames/geometry are defined, around P:219/247, and give the reader an explicit reason to consult it. |
| Table S4, `tab:RouteBio10`, S:118–143 (§1.2), p2 | None | None | Cite it with the biomimetic route definition, around P:157/499/510 as appropriate after coordinate reconciliation. |
| Table S5, `tab:RouteRedesign20`, S:145–173 (§1.2), p3 | None | S:206/244 | Add one purposeful reference in the route-redesign Methods/results summary, together with the clarified route/effective tendon definitions (T10/OP09). |
| Figure S4, `fig:cam_route_end`, S:306–310 (§1.6), p10 | None | S:299 | Pair its reference with the existing Figure S3 call at P:600 if both parts of the cam sequence are retained. |

Thus only **Tables S3 and S4 have no specific callout anywhere**; the other three are internally referenced. A direct main-text cross-reference is useful when it connects an essential model definition to its supporting inventory, but adding five arbitrary references is not the goal. Use purposeful grouped calls and keep the supplement self-contained. Its published form is not typeset by the journal, making clear presentation the author's responsibility under the [supplement guidance](https://www.frontiersin.org/guidelines/author-guidelines).

All other supplementary displays have direct Main references: Table S1 at P:414/418; Figure S1 at P:594; Figure S2 at P:546; Table S6 at P:639; Figure S3 at P:600. Table S1's conditional consolidation/removal and Figure S2's proposed removal follow U01/U02. **Table S1 is not a literal duplicate**: it uniquely includes Hunt/hybrid common-domain comparisons. Figure S1 is a separate route-comparison graphic. Do not conflate their labels when repairing the callouts.

### OP15 — Submission blocker 6: duplicate, uncited Hunt bibliography record

**O:95. Agree with consolidation; modify “safe to drop” as a blind operation.** Anchors: B:369–417 (`hunt_modeling_2017`) and B:1020–1029 (`hunt_modeling_2017-1`); active citation P:97 (§1, p2).

The records identify the same work: identical DOI `10.1115/detc2017-67458`, matching normalized title, 2017 year, and corresponding authors. They are **not byte-for-byte identical records**: the URL, title capitalization/bracing, author initials, page representation, and other fields differ. Main/Supplement cite only the non-suffixed key; the suffixed key is not used by the fresh main bibliography. Retain the cited key and consolidate any useful verified metadata before dropping the duplicate from a submission-specific bibliography. This is low-priority database hygiene, not a demonstrated duplicate rendered reference or a scientific blocker. Preserve byte-safe bibliography handling. The materially more important bibliography defect is T05: B:1085–1090 still labels the now-published force-model paper unpublished.

### OP16 — Retraction of the 38-page/live-paper claim

**O:3–5, 117–118. Agree with the live-preview retraction; retain the distribution-package finding.** Anchors: fresh `Main.log:1477–1479` records 33 pages; fresh supplement `.log:1285–1287` records 10. Artifact-level checks: bundled `paper/Main.pdf` = 38 pages and bundled supplement PDF = five; fresh outputs = 33 and ten.

Ben's live-recompile report is compatible with the independently verified fresh local build. This review did not independently open the live Overleaf preview, so attribute the live fact to Ben rather than claiming direct live verification. The stale ZIP-generated PDFs/logs do not prove a stale live preview. Conversely, retracting that live-preview allegation does not make the bundled PDFs current: replace/exclude those stale generated files before sharing/submission (T04). Both points can be true without changing manuscript results.

### OP17 — “All of Ben's figure drags landed”

**O:115–116. Modify the scope to artifact presence.** Anchors: dissertation `30-torque.tex:780` includes the new biomimetic pulley asset; paper P:401/402 still includes/captions the old abstract schematic (Figure 7, p16); paper P:503/505 still includes the hybrid Figure 11 (p22); P:619/622/623 use the differently styled Figure 13 panels (p27).

Presence of the new asset in the dissertation ZIP supports the bounded statement that it arrived there; it does not establish that every requested paper update was implemented or visually verified. The paper's `Redesigned_flexor_pressure_publication.pdf` and the dissertation's same-named revised asset have different SHA-256 hashes (prefixes `79814B2A` versus `0721D2DD`). Consequently, satisfaction of the dissertation plot-order directive U25 does not clear paper U15. Actual paper scatter alpha/edge/order remain export-verification criteria, not properties inferred from the dissertation asset. U08/U11/U15 remain valid paper actions. Keep the existing crossed-knee/two-BPA asset as the replacement starting point, apply Ben's latest reaction-label/simplification requirements, and move it to §4.4 as U09/U14 require. No return to the waived individual vastus/pulley example.

### OP18 — Proposed paper pass order, “quick blockers,” and relation to our findings

**O:97–103. Modify.** Anchors: P:386/600, P:176–194, S:69–140, B:1085–1090, P:410/499–505, P:581–643, P:619–625, P:708.

Do not prioritize unnecessary display removal or a nonexistent printed supplement bibliography over established definition/provenance issues. The highest-value paper sequence is: reconcile written definitions/records without rerunning calculations (T01/T02/T03/T10); apply Ben's relevance/organization/figure directives (U01/U02/U06/U08–U16); correct the published bibliography record and inventory (T05/T15); verify final callout order, graphics, and build outputs. Authorship/title alignment remains a human confirmation item (T11), separate from scientific competence. Main compliance algebra can be compressed without discarding the requested mechanics (T13).

The other audit usefully adds the rho collision and precise callout/bibliography checks, and corroborates T09's dense prose. It does not replace T06's hybrid carry-over defect, T07's measured final-size legibility problem, T12's omission of flexor hardware evidence in the abstract, or Ben's matched Figure 13 styling/27-site contribution requirements. No decision to postpone authorized work until after the defense is inferred from its suggested schedule. This report completes comparison and ledger updates only.

## Checks completed and limits

Compared the frozen consolidated audit, both paper TeX files, relevant bibliography entries, fresh aux/log/BibTeX files, the named current figures, prior T findings, and Ben's U ledger. Reproduced the other audit's prose counts only; no scientific calculation, optimization, simulation, fit, or metric was recomputed. Label/callout evidence is saved in `Archive/audit_20261030/other_audit_paper_checks.json`. Official journal checks above were read on 30 October 2026. The actual final target journal/section/article type has not been selected by this review.


---

# Response to the other audit's Chapter 2 expansion

30 October 2026. Task 3 only. This responds to every distinct Chapter 2 proposal/content claim in Round 2, lines 123–148 of the frozen `Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030.md`. The root reviewer handles scores, penalty weights, and rubric/classification claims. No manuscript, figure, bibliography, or scientific result was changed, and no new scientific calculation was run.

Authority is the latest supplied dissertation ZIP, SHA256 `be280d1909cef092861de856e43923ac120e4da0de5c1e31ba67f6c890902d6a`, extracted under `Archive/audit_20261030/dissertation`. References below are one-based lines of `chapters/15-background.tex` unless another chapter/file is named. The prior independent audit is `CHAPTER2_AUDIT.md`: six science findings C2-S1–S6 and seven prose/layout findings B2-L01–L07. Every rendered Chapter 2 page was already reviewed in Task 1; targeted pages and current assets were reopened for this comparison.

Verdicts: **Agree** means adopt the useful proposal; **Modify** means retain its purpose but correct its scope/remedy; **Reject** means the factual assertion or proposed rule is unsupported by the current artifact; **Verify** means an actual source/permission/model-specific record must be established before the proposed wording can be used. These verdicts concern the proposals, not Ben's abilities or the other reviewer's motives.

## Terminology: each listed term

| ID | Other audit's claim/proposal | Current evidence and verdict | Proposed disposition |
|---|---|---|---|
| C2-O01 | BPA appears at Chapter 2 line 8 but is not defined until50, a first-use failure. | **Modify.** The chapter-local observation is true (§2.1:8; §2.2.2:50), but active Abstract `02-abstract.tex:7` and Introduction `10-introduction.tex:8` already define braided pneumatic actuators (BPAs). It is not a dissertation-wide undefined acronym or a 42-line barrier for a reader who has read the introduction. | An optional brief reminder at Chapter 2's first BPA use is reasonable if chapters should stand alone. Do not count it as an unexplained scientific concept or demand a full second actuator introduction. |
| C2-O02 | PAM at 123 is never introduced. | **Agree.** Section 2.4:123 uses “PAM-actuated”; the text/caption at 123/139 then uses PAM routes without introducing that abbreviation in the earlier chapters. BPA has been defined, but that does not automatically define a second acronym. | Write pneumatic artificial muscle (PAM) once if the broader term matters, or consistently use BPA for these actual braided actuators. Do not replace a broader technology class with a narrower one without matching the referenced robot. |
| C2-O03 | RoM is never introduced. | **Modify.** Section 2.4:123 spells out “range of motion” and uses RoM later in the same paragraph. Introduction:8,20 also uses the full phrase. Parentheses are missing, but the meaning is already locally supplied. | Make the first local phrase “range of motion (RoM)” if retaining the abbreviation. This is a small consistency fix, not evidence that the chapter leaves joint motion undefined. |
| C2-O04 | Ia afferents appear at 166 before definition at 173. | **Modify.** Section 2.5:166 previews Ia pathways; Section 2.6:173 explicitly identifies groupIa muscle-spindle signals and function. The section-level handoff is real, but the term is not never defined. | Preview the three afferent classes with a short forward pointer, or move their detailed inventory to the sensory section. Avoid defining the same classes in both paragraphs. |
| C2-O05 | Ib afferents appear at 166 before definition at 173. | **Modify.** Section 2.5:166 lists Ib; Section 2.6:173 identifies Golgi-tendon-organ force signals. This is the same ordering issue as O04, not a separate missing biological mechanism. | Resolve together with O04/O06 in one short transition or reordered sentence. Preserve the distinction between Ib force and spindle feedback. |
| C2-O06 | II afferents appear at 166 before definition at 173. | **Modify.** Section 2.5:166 lists II; Section 2.6:173 gives the length-related role. Ia/Ib/II are accepted afferent-class labels, not acronyms requiring invented word expansions. | Use the existing biological definitions at the reader's first substantive need. A single coordinated improvement resolves all three labels. |
| C2-O07 | “Fictive locomotion” is not glossed. | **Agree.** Section 2.5:164 uses this experimental term to support the RG/PF distinction but does not explain what is recorded when actual limb movement is absent. That distinction matters to a mechanical-engineering reader. | Add a short contextual gloss at 164, for example locomotor-like motor output recorded without actual stepping in an immobilized preparation, scoped to the cited experiment. The cited [Rybak primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890439/) describes its immobilized preparation and motor recordings. |
| C2-O08 | V3 is never glossed. | **Agree, with scope.** Section 2.6:175 names V3 in the interlimb model; caption 203 repeats it without explaining the population's role. V3 is a recognized neuronal-class label, not an acronym that needs a fabricated expansion. | Give one short functional description tied to the cited Shinohara model, or omit the class name from the prose if it adds nothing beyond the schematic. Verify that description against that model rather than treating all commissural populations as interchangeable. |
| C2-O09 | Renshaw cells are never glossed. | **Agree, with scope.** First prose use is Section 2.5:166; Section 2.6:175 and captions193/226 repeat the name. Their recurrent-feedback role matters more than the eponym itself. | Add one short functional gloss at the retained first substantive use, based on the cited circuit. No separate physiological review or new circuit figure is needed. |
| C2-O10 | Locomotor primitives are never glossed. | **Agree/Verify the intended model meaning.** Section 2.5:166 lists them with reflex populations but does not say what “primitive” represents in the cited Di Russo model. It is a specialist modeling concept, not an unexplained abbreviation. | Explain its model-specific meaning in a few words if it is important to the argument; otherwise remove it from this inventory. Do not add a generic definition that equates neural modules, muscle synergies, and the particular mathematical input without checking the cited paper. |

The valuable terminology advice is a small, coordinated improvement in teaching sequence. It should not produce ten independent paragraphs or treat every use of an established technical name as defective.

## Duplication, reference, caption, and prose proposals

| ID | Other audit's claim/proposal | Current evidence and verdict | Proposed disposition and relation to our audit |
|---|---|---|---|
| C2-O11 | Chapter 2:50's BPA paragraph is verbatim Introduction:8. | **Modify.** The force-to-weight/compliance sentence is substantially shared, but the paragraphs are not verbatim: Chapter 2 adds bladder/braid mechanics, the Chou–Hannaford geometric term, and its modeling handoff; Introduction:8 develops the force-characterization motivation. | Shorten identical promotional wording while preserving the background's mechanism explanation and the introduction's motivating gap. A chapter-local definition can remain useful. Do not delete the whole paragraph on the premise of verbatim duplication. This is additional cross-chapter detail beyond our B2-L04/L07. |
| C2-O12 | Chapter 2:61's three-differences list duplicates Introduction:8. | **Agree, with a bounded remedy.** Both enumerate resting-length maximum force, smaller contraction ratio, and nonlinear pressure/contraction dependence. Chapter 2 gives the detailed physiological contrast and references; Introduction supplies the project motivation. | Retain the detailed comparison in the background and compress the introductory list or vice versa according to the introduction's purpose. Do not lose the force-model motivation. Also repair the Asano BPA citation in both locations under C2-S2; copying the same citation does not validate it. |
| C2-O13 | Chapter 2:270's ethics enumeration duplicates Introduction:6. | **Modify.** The introduction gives two examples, muscle removal and sensory-pathway blocking; Chapter 2 includes those plus RG lesions and spinal stimulation. This is overlapping motivation, not a verbatim duplicated list. Brief reorientation across chapters can be appropriate. | Compress the overlap without stripping Chapter 2's ethical/model-validity argument. More importantly, correct the overbroad healthy-human stimulation prohibition identified by C2-S4. Deleting repeated examples alone would leave that substantive error. |
| C2-O14 | Section 2.9:244 sends readers to the wrong computational-neuroscience section; use the neural/sensory sections. | **Agree.** The target Section 2.7 at 230–232 is one generic paragraph without the cat-model account. The same mistaken promise also occurs at 166, which the other summary does not flag. | Correct both references to the actual Sections 2.5/2.6 or the relevant circuit figures. Do not add duplicate model descriptions to 2.7 just to satisfy the reference. Direct overlap with C2-S5. Separately correct the Geyer/Markin taxonomy at 244 under C2-S3. |
| C2-O19 | One 524-word wall at 284–290; delete about100 weak words, horse aside first. | **Modify.** Section 2.10:284–290 is indeed overlong and the horse comparison at 284 contributes little to the robotics argument. However, the passage consists of four paragraphs across printed pp.38–40, not one uninterrupted paragraph; the exact524 count depends on stripping LaTeX/citations and was not independently adopted as a quality criterion. | Agree with targeted cuts to procedural chronology, the horse aside, and the ledger metaphor. Preserve the scientific limits of replacement. Do not cut a fixed number of words to meet an arbitrary threshold or add figures simply to interrupt prose. Direct overlap with B2-L05. |
| C2-O20 | Duplicate callout sentences at 34/36. | **Agree.** Section 2.2.1:34 introduces Figure 2.2;36 reintroduces it, inventories its mechanisms, and describes the schematic style; caption 44 then supplies the panel account. | Keep one purpose-led callout and remove duplicate narration. Direct overlap with B2-L06. |
| C2-O21 | The moment-arm figure is double-captioned at 148. | **Reject for the current figure.** Figure 2.6 at 145–150 has one LaTeX caption. The current rendered printed p.23 shows short panel titles (“Force moment,” “Moment-arm projection”) and schematic equations, not a second figure caption baked into the image. | Do not delete useful panel titles or natural schematic equations on a false duplicate-caption diagnosis. If the other reviewer meant repetitive explanation, identify the specific text rather than classify all panel headings as captions. This was checked directly against `moment_arm_methods.pdf` and the fresh page render. |
| C2-O22 | Figure 2.6's caption sends Chapter 2 readers forward to a Chapter 4 equation number. | **Modify/Agree as usability advice.** At 148, `eq:project_ma` resolves to Equation (4.2), defined in Chapter 4 `30-torque.tex:49–51`. It is a valid cross-reference, not a compile error. Chapter 2 Section 2.3:105 already defines the signed fixed-axis moment arm. | Prefer the local explanation/Section 2.3 reference or make the forward link optional rather than necessary for understanding the caption. The schematic can remain self-contained. This is a useful additional navigation observation, separate from the rejected double-caption assertion. |

## Reproduced-figure claims: each source needs its actual reuse basis

The other audit correctly notices differing credit wording. It does **not** establish that all three works are CC BY 4.0 or that the same caption template is appropriate for them. A citation, an accessible article, and a reuse license are different records. This review does not declare any existing reproduction unlawful or assume that no prior permission exists.

| ID | Example and exact anchor | Primary-source evidence | Verdict and proposed disposition |
|---|---|---|---|
| C2-O15 | Kengoro Figure 2.1, Section 2.1:16–20, especially caption 19. | The caption cites Asano2017 without a reproduction/license line. The [primary article record in PubMed](https://pubmed.ncbi.nlm.nih.gov/33157878/) identifies authors/AAAS copyright with some rights reserved. Local `thesis.bib:220–228` explicitly records the Science Journals Default License. The publisher full-text request returned403; no CC BY 4.0 grant was established. | **Verify/Modify.** Confirm the actual image's source and existing permission/reuse basis, then add the correct credit. Do not insert “under CC BY 4.0” from the other audit's template. An authors/AAAS record plus the local default-license metadata is evidence against assuming a CC license, not a legal determination about this thesis's reuse. |
| C2-O16 | Rybak Figure 2.11, Section 2.6:211–217, caption 216; source is the 2006 article's Figure 1, not the 2015 review. | The [primary 2006 article](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890439/) gives authors/Physiological Society copyright. The caption already identifies the source and original figure. No CC BY 4.0 grant was established by this inspection. | **Verify/Modify.** Keep the accurate source/figure credit and establish the relevant reuse terms or existing permission before expanding it. Do not relabel the 2006 figure with a 2015 paper's license or add a CC BY declaration merely for stylistic uniformity. |
| C2-O17 | Shevtsova Figure 2.12, Section 2.6:222–227, caption 226, source Figure 5.8. | The [publisher's chapter page](https://link.springer.com/chapter/10.1007/978-1-4939-3267-2_5) explicitly identifies Springer Science+Business Media New York copyright and supplies a Rights and permissions link. It does not establish the proposed CC BY 4.0 caption. | **Verify/Modify.** Check the existing reuse/permission record and use its required credit. Do not substitute a generic CC BY line. The proposed citation-specific check is merited; the assumed remedy is not supported. |
| C2-O18 | These three lack the CC BY 4.0 pattern “the other five carry.” | The current Chapter 2 source contains **four** explicit “under CC BY 4.0” caption statements, at 119,184,193,203: Steele, Deng, Di Russo, and Shinohara. The other audit's count of five is not reproduced. This count verifies caption text, not the actual licenses of all four sources. | **Reject the count; retain source-specific verification.** Neither the presence nor absence of a formulaic CC line proves that a reproduction has the correct rights basis. No blanket license change follows from a caption-style mismatch. |

## Asset and rendering limitation claims

| ID | Other audit's claim/proposal | Current evidence and verdict | Proposed disposition |
|---|---|---|---|
| C2-O23 | `figs/15-background` is not in the ZIP extract. | **Reject for the current authoritative material.** The directory exists in the isolated extract with 20 files. The unchanged authoritative ZIP itself was opened read-only and contains the same 20 nonempty Chapter 2 figure entries. All 15 unique graphics referenced by `15-background.tex` exist there. Its SHA256 matches the source manifest. | Correct this claim. A previous partial extraction may explain the other reviewer's experience, but the assets are not missing from the latest supplied ZIP or this review's extraction. Do not delay on an unavailable-file premise or copy older figures over them. |
| C2-O24 | Caption-versus-panel rendering remains unchecked; a pass against live ProofFinal is pending. | **Modify.** That can be an honest coverage limit of the other audit, but it is not the current combined audit's limit. Task 1 reviewed all 32 Chapter 2 pages in the fresh snapshot build; this pass reopened printed pp.11,23,30 and checked the relevant schematic/credit assets. The latest ZIP, rather than a potentially different working ProofFinal tree, is the requested authority. | Use the completed current-source render evidence and retain only genuine checks that remain, such as exact reuse permissions. Findings B2-L01–L03 already identify float placement and final-size readability problems. Do not represent the available render review as still waiting for a live session. |
| C2-O25 | The deficit is concentrated in cheap mechanical fixes. | **Modify.** Most proposed glosses/cuts are small editorial changes, but the scientific/citation problems in C2-S1–S4 and actual reuse verification cannot be certified by replacing strings or counting words. | Treat the other audit as useful editorial triage, then resolve the specific scientific claims, citations, model taxonomy, and final-size displays. Do not infer scientific submission quality or Ben's competence from its mechanical-fix count. Scores and penalty calculations are assessed separately by the root reviewer. |

Read-only asset verification: 20 ZIP figure files, 15 unique included assets, zero missing included assets. Existing renders are in `Archive/audit_20261030/chapter2_layout_qa`; specifically `pdf32.png`, `pdf44.png`, and `pdf51.png` support the targeted checks above. The complete earlier render coverage and page mapping remain documented in `CHAPTER2_LAYOUT_NOTES.md`.

## Comparison against all 13 independent Chapter 2 findings

This comparison concerns what the supplied Round 2 summary actually reports; it does not presume access to a longer unprovided session-agent report.

| Our finding | Chapter/source evidence | Coverage by the other Chapter 2 summary |
|---|---|---|
| C2-S1 | Section 2.2.1:30, DEA exclusion from all legged robots | **Not reported.** The primary DEAnsect counterexample requires a scoped scientific correction, independent of jargon/caption work. |
| C2-S2 | Sections 2.1/2.2.2:8/61, electrically actuated Asano humanoid cited as BPA example | **Not reported.** Its introductory-duplication observation also encounters the same mistaken BPA citation in Introduction:8, but does not test it. |
| C2-S3 | Section 2.9:244, Geyer pathways and CPG-free versus RG/PF taxonomy | **Not reported.** Catching the wrong cross-reference does not correct the underlying mechanism/lineage description. |
| C2-S4 | Section 2.10:270 versus 2.6:173, blanket prohibition of healthy-human stimulation | **Not reported.** Ethics compression alone would leave this source-specific scope error. |
| C2-S5 | Sections 2.5/2.9:166/244 to 2.7:230–232 | **Partly covered.** It flags244 but misses the analogous promised circuit explanation at 166. |
| C2-S6 | Section 2.2.3:65/95, cross-study accuracy ranking and universal model-limitation language | **Not reported.** Conditional argument correction, not a demand for new model fits. |
| B2-L01 | Sections 2.4–2.5:143–164; printed pp.22–25 | **Not reported.** Mechanical floats interrupt the first neural paragraph across pages. |
| B2-L02 | Section 2.4:143,154–158; Figure 2.7 at 0.62 textwidth | **Not reported.** The six-panel quantitative torque figure is difficult to inspect at its current inclusion scale. |
| B2-L03 | Section 2.6:175,190–203,211–226; Figures 2.9/2.10/2.12 | **Not reported because of its stated render limit.** Several circuit labels need more page area or a more selective presentation. |
| B2-L04 | Sections 2.5/2.6/2.8/2.9:164–166,175,207,236–248 | **Not reported in this form.** Repeated RG/PF/model inventories and six consecutive text pages 31–36 are more consequential than first-use acronym count alone. Its cross-chapter introductory repetitions are useful additional observations, not a substitute for this internal flow review. |
| B2-L05 | Section 2.10:270–290, especially284/286 | **Covered substantially.** Both reviews recommend compressing procedure chronology, the horse aside, and rhetorical accounting. We preserve the ethical stance and validation limits. |
| B2-L06 | Section 2.2.1:34–36 | **Covered directly.** Duplicate figure introductions can be consolidated. |
| B2-L07 | Chapter opening:4 and Section 2.8:240 | **Not reported.** Repeated background/contribution boundary assurances and third-person self-reference remain candidates for concise revision. |

The other audit adds useful targeted terminology, cross-chapter repetition, and source-credit checks. Its strongest shared recommendations are the duplicate figure introduction, wrong neural cross-reference, and ethics compression. The current evidence does not support a blanket first-use indictment, automatic CC BY captions, a duplicated moment-arm caption, or a missing-figures explanation. Scientific precision and final-page readability still require the separate checks identified in our independent audit. No manuscript changes were made during this comparison.



---

# Round 3 response: mechanics and Appendix C

30 October 2026. This document evaluates only B1–B9, B15–B19, and B31 in the frozen `Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030_V2.md` (SHA256 `0EC8EF1289FBE74504CA2BC8952AB20F66A4A3E31F0284D3114BAE95414B9B13`). These are 15 separate responses. Round 4 is outside this reviewer's assignment. No manuscript, figure, or result was edited, and no model evaluation or scientific calculation was run.

Source authority is the unchanged latest Overleaf ZIP snapshot under `Archive/audit_20261030/dissertation` and `.../paper`. Line numbers are one-based. Dissertation file `chapters/30-torque.tex` is Chapter 4; `chapters/94-AppendixC.tex` is Appendix C. Printed table identities were independently checked against the fresh build's labels, rather than inferred from source position. The earlier controlling user ledger is `USER_FINDINGS_DISS_MECHANICS.md`.

## Per-item evaluation

| Other ID / U ledger | Verdict | Correct current anchor | Evaluation and proposed response |
|---|---|---|---|
| **B1 / U20** | **Modify** | Chapter 4, Section 4.1.5, `30-torque.tex:417–434`, Table 4.3; mixed weight cells 425–426. | The reported mixed FLEX/EXT cells are real. Use two displayed lines per affected weight cell, or an equivalent visual arrangement, preserving both weights. “Both weight tables” is overbroad: the specific defect is in Table 4.3; Table 4.1 at 360–372 is the earlier flexor-only mesh objective, and Table 4.2 at 401–415 already distinguishes its configuration-specific terms by rows. Do not change those tables merely to satisfy an inaccurately generalized target. |
| **B2 / U21** | **Modify** | Chapter 4, within Section 4.1.5, heading `30-torque.tex:374`; constrained-redesign content 376–466; dense ending especially 451–466 and 464. | The missing visible number and dense ending are confirmed. The exact class/package cause was not established by merely observing the unnumbered heading. Our concrete proposal is to promote this material to a numbered Methods subsection with a dedicated label, while leaving the earlier placement method in 4.1.5. Avoid a global numbering change to protected Chapter 3. The other audit's440–460 range includes a table and equations, not just excess prose; trim the actual tail. Its separate proposed lineage cut at 342–348 should not silently remove the earlier method's provenance or the mechanical reason for redesign. |
| **B3 / U22** | **Agree, with placement/order distinction** | Chapter 4, Section 4.2.3, `30-torque.tex:594–599`, Figure 4.10; included asset at 596. | The current rendered legend is below the panels, so northeast in panel B is still pending. Apply the requested legend placement and the approved entry order. Legend-entry order and artist drawing order are separate properties; changing one does not necessarily change the other. Preserve the apparatus, retained data, and hybrid removal. Re-exporting a figure does not require rerunning a model. |
| **B4 / U23** | **Agree/Verify saved definitions** | Chapter 4, Section 4.3.3, `30-torque.tex:669–685`, Table 4.10; margins 680–681 and blank flexor worst angle 682. Requirement definitions 383–390 and 464. | The presentation defect and missing angle are correctly identified. Use percentages consistently, remove the plus sign, and distinguish achieved minimum margin from the required design allowance. Verify the saved field's denominator/domain before calling a margin “over human” or “over required.” Retrieve the flexor angle from the existing adopted result; do not invent it or rerun optimization. This remains distinct from T01's cam-margin issue in Table 4.11. |
| **B5 / U24** | **Agree, add the existing-model distinction** | Chapter 4, Sections 4.3.3–4.3.4, `30-torque.tex:715–718,729–733`; Table 4.11 lists 269–507 kPa and the 150-kPa model floor. | Correctly treats approximately 200 kPa as Ben's present-manifold operating experience, not a newly established universal threshold. Explain that the manifold's opening/closing requirement can restrict minimum controlled nonzero torque when the actuator remains taut. Preserve and identify the 150-kPa saved design assumption rather than silently replacing it with 200 kPa. Attribution to his experience is sufficient if no separate measured record is being claimed; no new test protocol or calculation is requested. |
| **B6 / U25** | **Reject the shared-asset/current-status assertion; separate the two projects** | Dissertation Chapter 4, Section 4.3.4, `30-torque.tex:742–747`, Figure 4.14; paper `Main.tex:622`, Figure 13B. | The same basename refers to two different PDFs. The dissertation already has prediction-first/measured-last drawing order, both scatter sets at the requested 0.75 opacity, and no marker edges. The paper retains measured-first/prediction-last order, opaque gold prediction markers, and gold-marker edge strokes. Therefore U25 is already complete in the dissertation; the paper setting defect remains. A future common export can be deliberately copied into both projects after verifying data and styling, but the current files are not one shared artifact and no dissertation regeneration is needed to complete U25. See the exact read-only evidence below. |
| **B7 / U26** | **Agree on deletion; correct the section boundaries** | Chapter 4, Section 4.3.5, `30-torque.tex:751–788`; mount derivation 763–774; “angle range” at 788. | Ben's newer instruction supersedes the previous addition of the mount-reaction derivation. Delete that calculation/discussion; keep the force/travel relations and shorten the explanation to approximately two paragraphs. Use RoM. The other audit's749–798 range crosses subsection boundaries:749 is the preceding physical-section bridge;790 begins Whole-Leg Transmission Capability;792–798 contain the master's comparison and its table opening. Do not accidentally delete the whole-leg results while compressing the pulley subsection. Coordinate retained reaction labels with B31/U17. |
| **B8 / U27** | **Agree, with range/criterion distinction** | Chapter 4, Section 4.3.6, `30-torque.tex:792`; Appendix C Sections C.10/C.10.1, `94-AppendixC.tex:454,458,500`. | Correctly notes that the exact phrase “graded range” is absent and the sign/magnitude criterion exists at 792. The explanation is still dense. Define the evaluated joint-coordinate grid/range briefly, and keep its included-cell domain separate from the torque success criterion. Replace “graded surface/domain” with plainer wording where practical. Do not invent new ranges or change the coverage denominator. |
| **B9 / U28** | **Agree, retain the theoretical conditions** | Chapter 4, Section 4.3.6, `30-torque.tex:792`; Appendix C Table C14 caption/rows, `94-AppendixC.tex:458,465–489`. | This is the right location for the rationale behind the stated diameter substitution. Add one concise ideal force-scaling statement: four parallel 20-mm actuators equal one 40-mm actuator in theoretical maximum force under matched ideal pressure, braid geometry, and contraction conditions. Do not imply measured equality of the full actual force–length maps or that the aggregate coverage gain isolates diameter. No 40-mm model fit or search is requested. |
| **B15 / U03** | **Reject the table identity/remedy as written; agree with Ben's actual duplication finding** | Appendix C, Section C.9, `94-AppendixC.tex:433–450`, **Table C13 `tab:TorqueFit`**; Chapter 4 Section 4.2.4, `30-torque.tex:617–632`, **Table 4.8 `tab:BioTorqueGof`**. | C13 is the biomimetic goodness-of-fit table, with the same three configurations and numerical entries as Table 4.8. It is not the configuration-survey table at 456. Remove that duplicate and update its references, retaining the main table. The survey inventory at 456–495 is **C14**, not C13, and the fact that Chapter 4 points to it does not make its unique data redundant. Do not delete or move C14 under this mistaken C13 mapping. |
| **B16 / U04** | **Modify** | Appendix C, Section C.9.1, `94-AppendixC.tex:429–430`; Chapter 4 Section 4.1.4 definitions, `30-torque.tex:316–336`. | The metric-definition repetition is real, but Chapter 4 already defines the metrics explicitly in prose/equations before the Results captions. The audit's claimed first expansion at 486 is not the operative location. Delete the duplicated C.9.1 material; do not add a second new Chapter 4 definition or reopen FVU. With U03 removing C13, remove the empty C.9 container and redirect its genuine references to the existing Chapter 4 definitions. |
| **B17 / U05** | **Agree with the purpose; correct table anchors** | Appendix C, Section C.10, `94-AppendixC.tex:452–495`; **C14 at 456–495**, caption 458/label 459. Chapter 4 Section 4.3.6 master's credit at `30-torque.tex:792`. | The Appendix C survey account omits the master's-project attribution and uses avoidable “graded surface” jargon. Add the inherited master's-model provenance and plain-language purpose, preserving the model/version distinction. The other row's “Table C14 at 502” is incorrect:502 begins **C15**, the per-site coverage table with caption 504/label 505. Apply the attribution/caption change to C14 at 458, not by editing the wrong table. |
| **B18 / U18** | **Reject the identities/locations; keep the user's actual unit-conversion directive** | Appendix C, Section C.1: **C5** is flexor inputs at `94-AppendixC.tex:105–118` (header 111); **C6** is extensor inputs at 120–130 (header 126); **C7** is measured final dimensions at 133–147. C1 header 28 and C4 header 86 already give resting/tendon lengths in m. | The other audit is shifted by one table: its C5 is actually C6, its C6 is actually C7, and its “C4 at 105” is actually C5. Convert the requested C5/C6 numerical linear dimensions and headings from mm to m, preserving the physical values, pressure, force, and model identities. Do not change the measured dimensions table C7 as an accidental substitute. Exact display conversions and verification are already specified under U18. |
| **B19 / U19** | **Agree** | Appendix C, Section C.4, `94-AppendixC.tex:205–257`; duplicated display 207–214. Chapter 4 stiffness/compliance at `30-torque.tex:140–152,214–232,249–250`. | Replace the repeated projected/series compliance display with a pointer to the retained Chapter 4 relations. Keep the unique per-configuration stiffness arrays, bounds, units, and authoritative implementation mapping. This compression does not resolve or remove the separate A1 frame-convention contradiction. No mechanical parameter changes follow. |
| **B31 / U17** | **Agree on notation; coordinate scope with B7** | Chapter 4, Section 4.3.5, `30-torque.tex:763–774`; Figure 4.15 included at 780, caption/label 781–782; the current graphic labels origin/mount reactions with R. | Use `F_R` for retained reaction-force labels, with a location subscript if necessary. Implement together with B7's deletion so the removed mount-reaction calculation is not retained merely with renamed symbols. Preserve the two-BPA force notation, force/travel transmission, crossed linkage, and “Example of…” caption opening. The audit does not establish that MATLAB specifically is required; edit the existing schematic's appropriate source/generator rather than prescribing an unnecessary tool or model evaluation. |

## B6: compact current-artifact evidence

The read-only vector inspection is saved as `B6_CURRENT_PDF_MARKER_EVIDENCE.json` beside this report. It records the current asset paths, SHA256s, and colored vector paths with their sequence numbers, fills, opacities, stroke colors, and bounds. Legend markers are identifiable separately from the 15 plotted points in each series. No figure was regenerated.

| Property | Dissertation Figure 4.14 | Paper Figure 13B |
|---|---|---|
| ZIP-relative asset | `dissertation/figs/30-torque/Redesigned_flexor_pressure_publication.pdf` | `paper/Figures/Redesigned_flexor_pressure_publication.pdf` |
| SHA256 | `0721D2DDFC6D354E1736FDF08245C50868C08AE2D45A55F3F090978DEE296AA5` | `79814B2A964EA67B880911541DF16EEDCAF74EB5F8517AD509C7D8897D33CE1B` |
| 15 blue measured fills | Sequence 87–101 | Sequence 72–86 |
| 15 gold prediction fills | Sequence 72–86 | Sequence 87,89,…,115 |
| Measured fill opacity | 0.7490196, PDF's 8-bit representation of 0.75 | 0.7490196 for plotted points |
| Prediction fill opacity | 0.7490196 | 1.0 |
| Prediction marker edges | None | Separate dark-gold strokes, interleaved with prediction fills |
| Directive status | **Verified complete: predictions first; measured last; both sets at 0.75, no edges.** | **Verified incorrect current settings: predictions later, opacity 1, edges present.** |

This evidence supports closing dissertation U25 and replacing any pending paper verification status with a confirmed presentation defect. A shared filename and similar appearance are insufficient evidence of project synchronization.

## Documentation repair completed

Before this Round 3 comparison, the three accidentally spaced file/path names in `OTHER_AUDIT_CHAPTER2_RESPONSE.md` were repaired with exact literal replacements after checking their real names: `CHAPTER2_AUDIT.md`, `CHAPTER2_LAYOUT_NOTES.md`, and `Archive/audit_20261030/chapter2_layout_qa`. Technical identifiers were preserved verbatim while natural-language spacing was corrected. The previous Chapter 2 source/figure conclusions are unchanged.

All 15 assigned B rows are addressed. The highest-consequence corrections are the already-complete dissertation figure versus the different paper figure (B6), the C13/C14 mistaken identity (B15), and the C5/C6/C7 shifted identities (B18). Missing flexor-angle and margin-definition provenance remain checks against existing records; they were not filled by a new calculation. This file changes the audit/ledger only.




---

# Other audit, Round 3: SADb findings and session-note corrections

Audit/proposed-change ledger only. Reviewed the newly frozen `Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030_V2.md`, Round 3, rows B10–B14 and all four session notes. No manuscript, application, figure, or live review state was changed.

Locations below are one-based lines in the latest ZIP snapshot, `Archive/audit_20261030/dissertation/chapters`, unless explicitly prefixed **Paper**. `50` means Chapter 5, `50-database.tex`; `70` means Chapter 7, `70-overall.tex`; `30` means Chapter 4, `30-torque.tex`; `96` means Appendix E, `96-AppendixE.tex`. Paper lines refer to `Archive/audit_20261030/paper/Main.tex`. The earlier full evidence and all 42 historical comment-card dispositions remain in `USER_FINDINGS_SADB_AND_COMMENT_CROSSWALK.md`.

## B10–B14

| Row | Verdict | Current evidence and corrected proposed action |
|---|---|---|
| **B10 — crop Figure 5.2** | **Accept the crop; qualify the fixed-width prescription.** | **50:39–44**, image inclusion at **50:41**, `figs/50-database/sadb_pivot_app.png`, printed p.125. The 1800 × 1100 px capture has extensive blank canvas below the table. Crop to useful content with a small margin, then verify labels at final page size. Retaining `0.96\textwidth` is not itself a legibility check: the dense columns may require a deliberate subset or landscape presentation. Preserve the values and multi-tag counting explanation at **50:37**. |
| **B11 — §5.3.2 future improvements belong in Chapter 7** | **Reject the paper mapping; accept Ben's actual Chapter 5 restructuring.** | The current dissertation's **§5.3.2 is “Further Improvements,” 50:61–67**, printed p.127. Future curation/classification/sensory-coverage work is at **50:63–65**. The segue Ben praised is **50:67**, linking accessible physiological evidence to implementation/testing in Chapter 6. Consolidate future work into the existing database subsection **§7.2.3, 70:59–61**, while preserving the segue and necessary present-coverage caveat. Paper §4.3.2 at **Paper:604–627** is the separate physical-validation section; changing its “hedge tails” does not execute this Chapter 5 instruction. |
| **B12 — monthly scraper is built** | **Accept implementation status; reject unsupported “built and running” as a routine scheduling claim.** | Chapter 5 has no scraper account; its curation Methods are **50:12–16**. More directly, **70:61** incorrectly says the extension “is being implemented.” The implemented `SADb_audit/watch/monthly_watch.py:394` main workflow and saved `watch/reports/2026-10.md`/shortlist establish that it is built and has produced an output. `watch/README.md:72–92` specifies OpenAlex/Crossref queries, deduplication, DOI verification, and manual-curation output; **README:110–114 says scheduling is provided, not installed**. Propose a short completed-work description in Chapter 5 and correct Chapter 7's status. No extra description from Ben is needed to identify these verified functions. Do not claim unattended monthly execution or complete OpenAlex coverage of the saved run, whose run notes disclose a source-access limit. |
| **B13 — app pictures closer to the front** | **Accept the goal; qualify the chosen figure and promised page position.** | **50:8 already cites the bubble-map figure**, but its sideways float is at **50:30–35**, rendered p.124 after the chapter starts on p.119. Moving its source block may help, but cannot guarantee p.2 without checking float placement in a build. More importantly, advisor **5.23's unfulfilled request is a Table-tab figure in the main text**: **50:24** still points to Appendix E's example at **96:16–21**. That old image visibly includes irrelevant archived/generic records and should be replaced with a current retained, relevant example before promotion. An early Table example near **50:8–10 or 24** addresses both the front-loading goal and the specific advisor miss. Preserve the bubble-map's **six callouts A–F**, the explanation at **50:26–28**, and its useful interface organization. |
| **B14 — all today's advisor comments satisfied by Ben's list** | **Reject the completion claim; accept the list as directive/disposition evidence.** | The available captured ledger contains **42 cards from 29 October**, not a complete fresh 30 October live census. Current source still leaves advisor **5.23** unaddressed (**50:24; 96:16–21**); historical cards **4.1** and **5.25** have uncertain live anchors. A supplied source ZIP and Ben's findings list cannot prove all current live comments are resolved. Record every known item, preserve Ben's explicit exclusions/deferments, and state the missing live verification. Browser unavailability does not turn an unknown into “satisfied.” No blanket tracked-change acceptance or live-thread resolution is authorized by this audit. |

## All four Round 3 session notes

| Note | Verdict and evidence |
|---|---|
| **1. “Ben's 5.3.2 = paper Physical Validation”** | **Incorrect mapping.** The dissertation table of contents and **50:61–67** identify the actual §5.3.2 as “Further Improvements.” Its good segue and future-work paragraphs match Ben's wording directly. Paper **§4.3.2, Paper:604–627**, remains separately covered by his important-results-first instruction. Follow B11's corrected dissertation mapping. |
| **2. Paper Figure 11 = biomimetic extensor; Figure 13 = physical-validation pair** | **Accept the identities.** **Paper:503–505** includes the biomimetic-extensor graphic/caption; **Paper:616–625**, particularly includes at **619 and 622**, assembles `Ext20_current_publication.pdf` and `Redesigned_flexor_pressure_publication.pdf` into the physical-validation pair. The paper audit visually checked their rendered Figures 11 and 13. This accurate figure mapping does not support the incorrect section mapping in note 1. |
| **3. Remove `R_mount` under Ben's newer ruling despite earlier vector-equation directive** | **Accept the precedence.** Current dissertation **30:763–774** contains the mount-vector/magnitude calculations, and **30:776** discusses the related design load accounting. Ben now directs removal of that calculation/discussion and a shorter transmission account. Preserve the ideal transmission relations at **30:755–761**, including `F_t = F_BPA / n`, tendon travel, and the two-BPA input `2F_20`. The separate BPA-origin reaction relation at **30:774** is not itself a mount-reaction calculation; do not delete it indiscriminately with the surrounding paragraph. The new directive controls the proposed ledger; it is not an unresolved permission question. No deletion was implemented in this audit. |
| **4. Browser panel unavailable; treat Ben's list as sufficient to satisfy B14** | **Reject this inference.** A list can establish requested changes and exclusions, but cannot establish the current live comment population or anchor/resolution state. The source/capture mismatch and actual open source condition at **50:24; 96:16–21** remain. Keep the motor-square thread unverified unless its current anchor is recovered; do not infer its resolution from unrelated Chapter 7 graphics or assign ownership from an unsupported mapping note. |

These are corrections to the audit and proposed-change ledger. The later Round 4 AI-methodology withdrawal and its simultaneous retention of the excluded Chapter 5:50 integrity task are evaluated in the master response; this supplement does not reopen that task.


---

# Other audit, Round 3: paper/supplement items and abstract

30 October 2026. This evaluates only B20–B30 and B32–B33 (13 rows), plus the quoted abstract, in the frozen `Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030_V2.md`. No manuscripts, figures, model inputs, results, or scientific calculations were changed. Round 4 is evaluated separately by the lead reviewer.

**Source key:** V2 = that frozen comparison file; P = `Archive/audit_20261030/paper/Main.tex`; S = `.../paper/frontiers_SupplementaryMaterial.tex`; A = `.../dissertation/chapters/02-abstract.tex`. Paper/supplement pages refer to the fresh 33-page and 10-page builds. T IDs refer to `TORQUE_AUDIT.md`; U IDs to `USER_FINDINGS_LEDGER.md`/`USER_FINDINGS_PAPER.md`. All line numbers are one-based.

The revision incorporates most of Ben's requested paper changes accurately. Five corrections matter before it becomes a change plan: the present route objectives are not Morrow's cost; result-first does not mean figure-callout-first; an additional placement decision from Ben is unnecessary; Table S1 is not a literal duplicate; and the proposed abstract still occupies two pages in the actual dissertation formatting.

## Responses to all 13 assigned rows

| Other item / V2 location | U/T mapping and exact manuscript anchor | Verdict and evaluation |
|---|---|---|
| **B20**, V2:200 — segue/register and missing 27-site results | **U06; T12**. P §1:99–101, pp2–3; §2.1:108, p3; existing §4.4 results:637–656, pp27–28. | **Agree.** Replace the anthropomorphic “force model must see beyond” phrase with a direct description of the empirical wrapping correction. Include the existing 27-site transmission/configuration comparison in the contribution framing, alongside the corrected joint model and physical knee tests. Preserve its inherited-routing and combined transmission/sizing/baseline-diameter scope; do not imply that all 27 sites were physically validated or newly rerouted. |
| **B21**, V2:201/220–248 — one-page abstract proposal | **U07**. A:1–7 is the active abstract; A:9–15 contains commented older material. Abstract proof: `OTHER_ABSTRACT_PROOF.json` and `Archive/audit_20261030/other_abstract_proof/other_abstract_proof.pdf`. | **Modify content; reject the claimed word count and one-page delivery.** The active current abstract is 228 whitespace-separated words, not 606. The exact proposed quote is 206, not the claimed 235. An independent compile using the existing dissertation class/preamble, spacing, margins, font, heading, and page numbering returns code 0 but produces **two pages**. Thus it does not satisfy Ben's one-page requirement. Detailed content assessment follows below. |
| **B22**, V2:202 — old pulley schematic replacement | **U08**, coordinated with **U09/U14/U17/U26**. P §2.9:397–404, Figure 7, p16; replacement starting asset in dissertation `figs/30-torque/reverse_pulley_biomimetic.pdf`, active at Chapter 4 `30-torque.tex:780–782`. | **Agree with replacement; correct the stated defect.** The paper uses the older abstract schematic without the biomimetic knee. However, its actual current artwork does **not** carry `F_BPA(n−1)`; it shows the vector `R_c` relation. Do not repeat that false asset description as verified evidence. Carry the crossed-link, two-BPA schematic into the paper after applying Ben's current reaction-force notation and mounting-reaction simplification. Put equations in text except natural schematic labels, and use the short “Example of …” caption. This is a proposed figure change, not a reason to recalculate forces or revive the waived individual pulley example. |
| **B23**, V2:203 — pulley after meet/exceed section | **U09/U14**. Current P §2.9:397–404, pp15–16; physical evidence §4.3.2:604–627, p26; §4.4 starts at629, pp26–28. | **Agree.** Move the substantive mechanism explanation and replacement figure into §4.4 after the knee capacity evidence. The transmission comparison is already there at P:633 onward. Retain only the minimum Methods statement/cross-reference needed to describe what was evaluated; avoid two full explanations. One relocation serves B22/B23/B28 together. |
| **B24**, V2:204 — identify/remove Morrow's spelled-out cost | **U10 already satisfied; T08/T09**. Morrow is cited in P §1:97, §2.4:166, and §4.5:666. Present study's common route cost: P §2.9:388–393; configuration-specific `J_flx/J_ext`: §4.3.1:583–592 (actual align block584–591), pp24–25. | **Reject the attribution and resulting remedy.** The two expanded objectives are explicitly the present flexor/extensor route objectives under Eq.35. They are not identified as Morrow's master's cost. Current P already mentions Morrow without spelling out his objective, satisfying U10. Do not replace this study's method with a citation to a different formulation, or substitute Bolen 2019 for Morrow 2020 as interchangeable provenance. The present penalty inventory can be compressed/moved to Methods or supplement under U12/T08/T09, with its own accurate definitions retained. |
| **B25**, V2:205 — hybrid removal from Figure 11 | **U11; T06**. P Results opening:410, p16; §3.3:499–505, pp17/22; graphic include503; caption505. | **Agree.** The paper still contains the gold approximate-hybrid points, legend, caption phrase, and related narrative. Use the verified no-hybrid presentation while preserving the apparatus, remaining curves/points, and their data. Remove the premature hybrid qualification at410 as well as the §3.3 comparisons. Retain the pinned-knee hybrid method/evidence elsewhere. This also helps repair the first Figure11-before-Figure8 citation defect. |
| **B26**, V2:206 — §4.3.1 too long | **U12; T08/T09**, with **T01** definition dependency. P §4.3.1:581–602, pp24–25; the source-counted 770-word prose run is592–603. | **Agree with compression; reject B24 as its factual basis.** Cut repeated setup, replace numeric inventories with a compact capacity comparison, and refer to existing supplementary geometry/release details. Preserve the study's actual cost definition and paired-route distinction. Do not delete method/provenance clauses merely because they look administrative if they distinguish optimizer geometry, adjusted hardware inputs, or comparison targets. Resolve the written margin definitions from saved records before a final summary table; no optimization rerun. |
| **B27**, V2:207 — §4.3.2 important results first | **U13; T08/T12**. P §4.3.2:604–627, p26; main outcomes at610/612; Figure13 at616–625, p27. | **Modify the proposed order.** The identified problem is correct, but “figure callout first sentence, then meet/exceed headline” reverses Ben's requested priority. Open with the combined measured-capacity result—both redesigned routes met/exceeded their respective selected human-muscle targets at all nine extensor and fifteen flexor measurement angles—then cite Figure13 and give margins/essential conditions. Preserve the discrete-point scope and the distinction between human-target capacity and pressure-specific prediction agreement. A figure callout can be part of that result sentence; it should not replace the topic claim. |
| **B28**, V2:208 — pulley schematic in §4.4 | **U14/U08/U09**. P §4.4:629–643, pp26–28; old figure block399–404, p16. | **Agree.** §4.4 is the authorized canonical home, after the meet/exceed hardware results and before detailed whole-leg/transmission comparisons. Share the one move with B22/B23. Keep the master's-model attribution at P:631/639 and verify first-reference/float placement after eventual edits. |
| **B29**, V2:209 — matched Figure13 panels, wording, scatter style | **U15; T07**. P §4.3.2:606–623, pp26–27; asset includes619/622 and caption623. | **Agree, with two qualifications.** “Cam extensor” is verified in the caption at623; current prose already uses “cam-assisted extensor” at606/664/680, so there is no need to pretend that prose wording remains uniformly uncorrected. The panel font/line scales visibly differ. Export a matched pair and set panel-B prediction before measured, alpha 0.75, and no marker edges, preserving data. The paper's flexor PDF hash starts `79814B2A`; the dissertation's revised same-named asset starts `0721D2DD`. Therefore the dissertation's satisfied U25 does not clear paper U15, and “one regeneration fixes both” is true only if the correct export is actually installed and checked in both projects. No additional model evaluation is required. |
| **B30**, V2:210 — §4.4 cost or supplementary method/table | **U16**; distinguish **T08/T09**'s knee-route objective. P §4.4:631–643, pp26–28; “saved transmission-ladder winner” at639; supplementary site/size table S:255–296, p8; knee-route costs P:388–393/583–592. | **Agree with documenting the existing rule; remove the unnecessary Ben-choice gate.** Ben has already authorized either placement. A concise §4.4 explanation plus exact supplementary rule/table is a routine editorial choice, not missing permission. The actual selection expression must still be recovered from existing saved code/campaign records before being stated; a directive document alone is not proof of the implemented objective. Distinguish length optimization, feasible-candidate ranking, and any final engineering judgment. Do not invent a formal cost, weights, or assessment method and do not rerun the search. |
| **B32**, V2:217 — supplementary Table S1 redundant | **U01; T15**, prior supplement relevance advice superseded conditionally. S §1.1:47–66, p1; P §3.1:414/418, p17; main Table1 P:451–495, p21. | **Modify “confirmed duplication.”** The table overlaps the pre-correction discussion but uniquely reports Hunt/improved/hybrid statistics on a common model-defined domain. Main Table1's improved-versus-corrected full-test statistics are not interchangeable. Prefer removal/consolidation if the unique comparison is unnecessary; preserve any essential claim and its domain, and repair both callouts. Do not fold unmatched-domain numbers together or insist on keeping S1. Ben said “seems redundant,” which supports review and a removal proposal rather than a false claim of literal duplication. |
| **B33**, V2:218 — supplementary Figure S2, unused42-cm test | **U02; T15**. S §1.3:186–192, Figure S2, p5; P §4.1:546, p23; inventory708, p30. | **Agree with the removal proposal; qualify the claimed verification.** The current item is the slippage photograph. Ben supplies its unused42-cm-test provenance; the caption itself does not identify resting length. The claimed explicit paper-text exclusion was not found in P/S; the six-test inventory at P:382 does not by itself establish this photo's provenance. Remove the illustration, callout, and inventory entry as proposed, while checking whether the photograph-specific ±20-mm claim at546 has independent support. Do not erase a generally supported path-sensitivity explanation solely because the photo is removed. |

## The quoted abstract: merits, defects, and fit result

**Anchor:** V2:226–245, linked to B21 at201; current active abstract A:7. **Verdict: useful content draft, not a compliant one-page replacement.**

Its force-model → knee-model → tools sequence is coherent. It retains the three-coefficient force surface, the effective-length/compliance/wrapping contribution, both redesigned routes, selected human-muscle references, and the distinction between discrete isometric measurements and dynamic performance. These are appropriate claims already supported by the dissertation. It gives the database and neural tools more explicit roles than the current abstract's generic ending.

It can be tightened further: “long-running effort at Portland State University,” the extended system description, and repeated sequence scaffolding consume scarce space. First-person statements of Ben's contributions would also better match his dissertation style. Naming SADb Explorer and stating the concrete demonstrated neural-tool contribution would convey more than generic “supporting tools” language; qualify any broad “verified” wording to the selected implementations/tests. The existing 27-site comparison can be acknowledged if retained as a central contribution, but no new result needs to be added merely to make an abstract sound more comprehensive.

The fit evidence is decisive: `Notes/publication_audit_20261030/OTHER_ABSTRACT_PROOF.json` records **206 whitespace words, two pages, return code0**, with the exact quoted text and `source_unchanged=true`. The isolated proof PDF is `Archive/audit_20261030/other_abstract_proof/other_abstract_proof.pdf`; `build.txt` records the compile. The claimed 235-word count and one-page implication are therefore incorrect under the actual current formatting. The current A:7 has **228** whitespace words; 606 includes material outside the active paragraph, including commented old drafts, and must not be reported as the current active length.

Our separate `ABSTRACT_PROPOSAL.md`/`ABSTRACT_PROPOSAL_20261030.txt` contains the **159-word** U07 alternative, independently compiled at the same formatting and visually checked at one page. Both remain proposals. The evidence supports further shortening or adopting that already fit-tested alternative; it does not support shrinking the font/margins or claiming a two-page quote fulfills Ben's requirement.

## Boundary

All 13 assigned B rows have a response and U/T mapping above. B31, Round 4, and other dissertation/appendix rows are outside this subreview and are handled separately. These findings update the audit/proposed-change ledger; they do not authorize silent result changes or imply any live Overleaf comment was resolved.


---

# Response to the other audit's later Round 4

The named report changed during Task 3. Its later version was preserved as Archive/audit_20261030/OTHER_AUDIT_REVIEWED_20261030_V2.md (SHA256: 0ec8ef1289fbe74504ca2bc8952ab20f66a4a3e31f0284d3114bae95414b9b13). The earlier responses explicitly assess the first frozen version; this section evaluates the added Round 4, including its retractions. The added Round 3 B1-B33 is answered separately, one row per item.

## Every Round 4 response-table item

These rows respond to the other audit's evaluation of our earlier findings, rather than inventing a second set of scientific defects. The source anchors refer to the unchanged authoritative ZIP snapshots.

| Their item | Chapter / source anchor | Response to their proposed disposition |
|---|---|---|
| P1 A1 | App.A92:198; App.C94:186,216,226-227 | **Agree.** This is a substantive saved-fit/present-implementation contradiction, not merely line305 versus306. Reconcile existing provenance and align the documentation; no refit. |
| P1 T01 | Ch. 4:709-722; paper Main:386,600 | **Agree with the check; reject the guessed argmin explanation.** If required torque is exactly 1.05 times the same positive human target on the same domain, the relative margins have a pointwise affine relationship, so taking their minima preserves that relationship. Different argmin angles alone cannot explain the discrepancy. A genuinely different domain, target, or definition could; establish it from saved records rather than declaring it likely. No number is selected or changed here. |
| P1 T02 | Paper §2.5 Main:176-194; Ch. 4:78-82 | **Agree.** Use a magnitude relation or explicitly signed arm/axis projection. This fixes exposition without alleging incorrect plotted calculations. |
| P1 T03 | Supplement Tables S2/S4:69-89,121-140 | **Agree.** Identify frames, historical/current versions, and the plotted configuration. Keep separate from the vastus-name and tendon-parameter questions. |
| P1 SNS-S01 | App.D95:125,133,136,155,220; implementation kine_ref.py:271-275,320-331 | **Agree, bounded to historical scorer documentation.** Map each existing result family to its actual scorer; describe terms and failure sentinel from those records. Current code alone cannot establish the scorer used for an older result. No rerun, metric redesign, or AI-tool assessment method is proposed. |
| P1 T04 | Paper ZIP Main.pdf/supplement PDF versus current fresh outputs | **Agree.** Remove/exclude stale generated delivery artifacts. Keep separate from the correctly retracted allegation about the live preview. Artifact metadata has no manuscript line number. |
| SADb S1 | Ch. 5:24; App.E96:16-21 | **Agree.** A relevant current Table capture remains needed. It need not be a new task for Ben if later implementation can capture the existing app; do not simply move the criticized old image. |
| SADb S2 | Ch. 5:16/63; App.E96:12 | **Agree with the note/source distinction; reject keeping the separate excluded AI task.** Establish whether the 38 missing fields are curation notes or unavailable source text. This is a genuine corpus-description issue. It does not justify revising Ch. 5:50 under an integrity label after Ben excluded that comparison. |
| SADb S3 | Ch. 5 Figure 5.2:39-43 | **Agree.** Crop the blank bottom and judge final-size readability; a focused native table or column subset is an option if the screenshot remains unreadable. Do not silently change counts or data. |
| SADb S4 | Ch. 5 Figure 5.1/callouts:26-35 | **Agree, preserve structure.** Improve readable detail without discarding the six callouts or their app-functionality explanation. Local crops/clearer labels are alternatives to a wholesale redraw. |
| SADb S5 | Ch. 5:8,20,22,37,63; App. E:8-14 | **Agree selectively.** Consolidate repeated inventories while preserving collection-versus-Papers denominators, present coverage, and the early app explanation. This is not a mandate to eliminate every normal recap. |
| SADb S6 | Ch. 5:8-10,20-22; App. E:4-14 | **Agree as an availability proposal.** Identify the existing retrievable app/version/archive used for the dissertation. Do not invent a release, identifier, deposition, or mandatory new hosting service. |
| SNS-S02 | App. D Table D5:220-231 | **Agree.** Resolve the duty statistic's caption/zero-cycle meaning using existing saved definitions; do not recompute scores. |
| SNS-S03 | App. D:60-64 versus246-269 | **Agree with distinct-configuration wording.** Identify original startup versus modified settled-transient tests and their different purposes. Avoid calling the older comparison a matched settled replay. Retain the configuration distinction rather than deleting all development context. |
| SNS-S04 | Ch. 7:20-28,Figure 7.1 | **Agree.** Align the claimed second acceptance check with the diagram, using a valid return route or narrower prose. This is the actual diagram mismatch, unlike the unconfirmed double-caption claim. |
| SNS-S05 | Ch. 6 Figures6.3/6.4:102-104,149-154 | **Agree with the visual remedy; do not assume its tool.** Enlarge the useful plot area and final-size labels, or re-export from its actual generator. A figure being small does not prove MATLAB regeneration is required. |
| SNS-S06 | Ch. 7:8,12,65,75 | **Agree with targeted compression.** Remove repeated inventories and split overloaded arguments. Do not adopt a universal delete-first word quota or remove useful synthesis. |
| SNS-S07 | Ch. 6:67,137,234-244; Ch. 7:16,65 | **Agree.** Explain verification questions once and keep experiment-specific assistance/limits where they affect interpretation. Remove generic repeated warnings. |
| SNS-S08 | Ch. 6:85; App. D:298-360 | **Agree, with placement flexibility.** One compact model/controller/support/result map would help. Prefer revising/placing the most useful existing inventory; a new appendix table is an option, not a requirement to add another duplicate. |
| SNS-S09 | Ch. 7:8,53 | **Agree.** Scope accuracy to first-pass static sizing supported by the reported comparisons; retain dynamic/controller integration as unfinished work. |
| SNS-S10 | Ch. 6 Figure 6.8:209-216 | **Agree.** Expose the saved phase statistic/criterion band and use comparable axes. Preserve actual values and distinguish illustrative traces from the aggregate statistic. |
| SNS-S11 | App. D:231,285,302,358 | **Agree.** Replace debugging aliases/adjectives with scientific model/status descriptions, retaining failed/assisted results and their evidence. |
| T05 | Paper bibliography:1085-1090 | **Agree.** Correct the published Actuators record from the primary DOI record, with byte-safe handling. This is more consequential than dropping an uncited duplicate Hunt entry. |
| T06 | Paper Figure11 Main:410,499-505; Ch. 4:592-598 | **Agree for the paper; distinguish work already done.** Dissertation hybrid removal is already complete. Carry it into paper graphics/caption/discussion; separately apply Ben's newest northeast legend directive. Do not claim both manuscripts still need hybrid removal. |
| T07 | Main Figures9/10/13; supplement Figure S1 | **Agree.** Fix actual final-size font/line scales. Do not infer compliance from the generating figure's nominal font size or from similarly named dissertation files. |
| T08 | Paper Main:581-626 | **Partly agree with their caution.** Cost placement/compression is useful. Moving the paper's hardware evidence into Results remains an optional journal-structure proposal; the earlier audit already explicitly excluded silently moving dissertation §4.3.4. Keep Ben's requested current locations and ordering as the ledger's implementation baseline. No new approval is needed merely to record alternatives. |
| T09 | Paper Main:382,596-602 | **Agree.** Replace repeated inventories with compact settings/capacity tables and concise interpretation. Remove duplicated prose after consolidation, preserving unique constraints and outputs. |
| T10 | Paper Main:614; supplement:145-160; App. C:151 | **Agree, now narrower.** The dissertation already states the 25-mm route parameter minus 20-mm fitting adjustment equals 5-mm effective input. Carry this definition into the paper/supplement, distinguishing physical 28/30-mm cables. This is disclosure, not an unresolved need to refit. |
| T11 | Ch. 4 publication note:5; paper author:48/contributions:696 | **Agree.** Human authors confirm the list/title; an agent does not add or remove coauthors. Record the mismatch without alleging misconduct. |
| T12 | Paper abstract:81; hardware result:612 | **Agree.** Include the flexor as well as extensor measured capacity evidence, scoped to selected references and discrete tests. “Extensor only” is an omission in contribution framing, not evidence of inflated measurements. |
| T13 | Ch. 4:140-265; paper:245-318 | **Agree moderately.** Retain the requested mechanics and core assumptions, condense equivalent substitutions, and use the technical appendix for genuinely useful derivation. Coordinate with Ben's removal of duplicate C.4 equations. |
| T14 | Ch. 4 Figure 4.12:697-703 | **Agree.** Restore frame/axis names and units at actual inclusion size and separate crowded contact labels. Use the existing accurate geometry rather than a speculative replacement. |
| T15 | Paper Supplemental Data:707-708 | **Agree; distinguish two different defects.** Correct the inaccurate inventory, then give purposeful references to retained supporting items. Missing item calls are not the same as promising the wrong fit table. Ben's U01/U02 override the earlier S1/S2 retention advice; S1 contains unique common-domain statistics, so preserve only essential evidence when consolidating/removing it. |
| A2 | App. B Figure B1:190 | **Agree.** Use the actual three row/condition letters A/B/C; do not assign nonexistent individual subplot letters. |
| A3 | App. C:322-333,397-408 | **Agree with crop/selection/readability advice.** Do not add an unsupported letterless-panel diagnosis. Check final-size labels and whether C2 contributes unique geometry before removing it. |
| A4 | App. C:11-20,99-105,151-255,387-391 | **Agree.** One authoritative configuration/source map plus unique equations/inventories is preferable to scattered implementation archaeology. Respect Ben's specific removals and master's provenance requirement. |
| A5 | App. B:148,194 | **Agree as optional presentation.** Retain headless-test scope, identity-test outcome, and failed dynamics. Repository notes can hold GUI task status and full identity strings. Do not silently change numerical values/precision. |
| Rubric | REVIEW_RUBRIC.md; other audit's supporting rubric | **Agree that they assess different objects.** Keep the concrete feedback checklist, but do not use the manuscript index/threshold to grade competence or certify readiness. Our rubric penalizes actual assistant failures and respects Ben's exclusions. |
| Frontiers check | Paper Main:729 template comment; official article-type guidance | **Agree with the retraction.** The initial15-display Original Research blocker is withdrawn in the newer audit. Retain article-type-specific word/display requirements, ordered calls, and axes/units checks. Do not describe the newer audit as still insisting on that cap. |

## Round 4's additional claims, corrections, and merged queue

| ID | Added item | Evaluation |
|---|---|---|
| R4-X01 | Our earlier audit lacked Chapter 2/terminology/source-markup coverage | **Historically fair for its bounded scope; now partly superseded.** That report explicitly audited Chapters4-7 and appendices/paper. Task 1 of this ordered follow-up completed the entire 290-line/32-page Chapter 2 review and 13 findings. Targeted terminology is now evaluated item by item. Static markup scanning is useful but cannot certify live review history. |
| R4-X02 | Correction1, retract float cap | **Accepted and recorded above.** The older claim remains visible in the frozen first version; the final combined verdict incorporates its retraction. |
| R4-X03 | Correction2, add five provenance classes | **Agree.** Preserve the specific source/result-family checks without new fitting, optimization, or fabricated definitions. They complement editorial screening. |
| R4-X04 | Correction3, withdraw AI-tool rating methodology | **Agree; apply consistently.** The same Round 4 retains Ch. 5:50 in “Where my audit stands” and the merged integrity-priority list. That would reintroduce the excluded comparison through an integrity wording task. Remove it from this task's queue as well; no substitute assessment or attribution exercise. |
| R4-X05 | Correction4, S2 photo may be Ben's excluded42-cm test | **Agree with Ben's removal ruling.** Current Figure S2 is the slippage photograph, identified by its asset/label; the resting-length provenance comes from Ben. Do not claim the generic caption independently establishes that test's exclusion. |
| R4-X06 | Distinctive contributions: first-use sweep, markup scan, process classes, licensing, wall counts | **Useful with qualifications already supplied.** Local terminology and callout checks add value. Avoid false missing-definition claims, automatic CC-BY captions, erasing requested selection history, or treating counts as competence measurements. Chapter 2's substantive scientific/citation/layout findings extend beyond those mechanical screens. |
| R4-X07 | Merged provenance-first queue | **Agree as a priority option, not a new gate.** Frames, margins, sign definitions, scorer provenance, citations, and authorship deserve clear dispositions. Record-level checks can proceed independently of straightforward authorized presentation proposals. No optimization rerun is warranted. |
| R4-X08 | Merge every B1-B33 directive and original P1-P3 without filtering | **Reject the unfiltered union.** Round 3 misidentifies Table C13, Tables C5/C6, §5.3.2, existing Morrow costs, and current Figure 4.14 state. The Round 3 responses correct each. Ben's actual intent outranks those mistaken mappings. |
| R4-X09 | Consolidated MATLAB regeneration queue | **Modify to verified need and actual generator.** Dissertation Figure 4.14 is already correct; the paper asset differs. Chapter 7's design workflow has no second caption. Target readability/caption/data-order changes in the affected assets, preserving already-correct output. A crop, caption edit, or vector/generator change may suffice; not every schematic is a MATLAB figure. |
| R4-X10 | No number changes anywhere | **Preserve result integrity, honor explicit formatting directions.** Meter conversion and percent display changes are Ben's requested proposals. Unknown margins/angles/frames must come from existing records, not guesses. No manuscripts were edited in this pass. |

No item in either version creates permission to revive FVU, the rejected uncertainty exercise, deferred classification work, the excluded AI comparison, or the waived individual pulley/vastus example. Both frozen versions and the source-specific responses are preserved so later readers can distinguish a corrected claim from an unreviewed one.


---

# Location index for the original audit and Chapter 2 expansion

Line numbers are one-based in the authoritative supplied ZIP snapshots under Archive/audit_20261030/dissertation/chapters and Archive/audit_20261030/paper. A line range is a navigation aid, not a claim that every sentence in the range needs revision. Artifact-level issues have no meaningful manuscript line; their relevant source pointers are supplied instead.

| Item | Chapter / section or artifact | File and source lines |
|---|---|---|
| T01 | Ch. 4 §4.3.3, Table 4.11; paper §§2.9/4.3.1 | 30-torque.tex:709-722; Main.tex:386,600 |
| T02 | Paper Methods, hybrid torque calculation | Main.tex:176-194; compare Ch. 4 30-torque.tex:78-82 |
| T03 | Paper supplement Tables S2/S4 | frontiers_SupplementaryMaterial.tex:69-89,121-138 |
| T04 | Paper ZIP's stale uploaded PDFs | Main.pdf/frontiers_SupplementaryMaterial.pdf metadata/content; current source Main.tex:655 and Figure12 at516-518; no PDF source line |
| T05 | Paper References, published force-model citation | BolenFrontiers22.bib:1085-1090; comparison thesis.bib:1086-1095 |
| T06 | Paper §3.3 and Figure11 | Main.tex:410,499-505; Ch. 4 Figure 4.10 at30-torque.tex:592-598 |
| T07 | Paper Figures9/10/13; supplementary Figure S1 | Main.tex:431,448,619-625; frontiers_SupplementaryMaterial.tex:179-183 |
| T08 | Paper §§4.3.1-4.3.2 | Main.tex:581-626 |
| T09 | Paper identification and redesigned route discussion | Main.tex:382,596-602 |
| T10 | Paper §4.3.2 and supplementary Table S5 | Main.tex:614; frontiers_SupplementaryMaterial.tex:149-157 |
| T11 | Ch. 4 publication note; paper author metadata | 30-torque.tex:5; Main.tex:48,696 |
| T12 | Paper Abstract and §4.3.2 | Main.tex:81,612 |
| T13 | Ch. 4 §4.1.4; paper compliance methods | 30-torque.tex:140-265; Main.tex:245-318 |
| T14 | Ch. 4 §4.3.3, Figure 4.12 | 30-torque.tex:697-703 |
| T15 | Paper Supplemental Data inventory | Main.tex:707-708; frontiers_SupplementaryMaterial.tex:47-66, Table S1 |
| DB-S1 | Ch. 5 §5.2; Appendix E FigureE1 | 50-database.tex:24;96-AppendixE.tex18-20 |
| DB-S2 | Ch. 5 methods/improvements; Appendix E coverage | 50-database.tex:16,63;96-AppendixE.tex12 |
| DB-S3 | Ch. 5 Figure 5.2 | 50-database.tex:41-43 |
| DB-S4 | Ch. 5 Figure 5.1 and callout walkthrough | 50-database.tex:26-35 |
| DB-S5 | Ch. 5 and Appendix E repeated coverage | 50-database.tex:8,20,22,37,63;96-AppendixE.tex8-14 |
| DB-S6 | Ch. 5 resource availability | 50-database.tex:8-10,20-22;96-AppendixE.tex4-14 |
| SNS-S01 | Ch. 6 scoring; Appendix D scoring definition | 95-AppendixD.tex125,133,136,155,220; current implementation comparison kine_ref.py271-275,320-331 (historical-version mapping required) |
| SNS-S02 | Appendix D Table D5 | 95-AppendixD.tex220-231 |
| SNS-S03 | Appendix D knee startup versus transient | 95-AppendixD.tex60-64,246-269 |
| SNS-S04 | Ch. 7 first-pass workflow, Figure 7.1 | 70-overall.tex:20 and figure block22-28 |
| SNS-S05 | Ch. 6 Figures6.3/6.4 | 60-sns.tex:102-104,149-154 |
| SNS-S06 | Ch. 7 long synthesis/extension paragraphs | 70-overall.tex:8,12,65,75 |
| SNS-S07 | Ch. 6/7 repeated verification framing | 60-sns.tex:4,14,29-39,67,137-139,234-244;70-overall.tex:16,65 |
| SNS-S08 | Ch. 6 configurations; Appendix D model/support map | 60-sns.tex:85;95-AppendixD.tex298-360 |
| SNS-S09 | Ch. 7 first-pass prediction claim | 70-overall.tex:8,53 |
| SNS-S10 | Ch. 6 Figure 6.8 feedback/phase comparison | 60-sns.tex:209,211-216 |
| SNS-S11 | Appendix D bookkeeping/table wording | 95-AppendixD.tex231,285,302,358 |
| APP-A1 | Appendices A/C contradictory bracket frames | 92-AppendixA.tex198;94-AppendixC.tex186,216,226-227 |
| APP-A2 | Appendix B Figure B1 letters/caption | 93-AppendixB.tex190 |
| APP-A3 | Appendix C FiguresC1-C4 | 94-AppendixC.tex322-333,397-408 |
| APP-A4 | Appendix C repeated inventories/history | 94-AppendixC.tex11-20,99-105,151,155,186,214-255,387-391 |
| APP-A5 | Appendix B administrative/precision narration | 93-AppendixB.tex148,194 |

Chapter 2's thirteen items have complete chapter/section/line locations in CHAPTER2_AUDIT.md. Ben's findings U01-U32 and the comparison audit's items are indexed separately in their stage-specific ledgers to preserve their provenance and avoid confusing them with new reviewer recommendations.
