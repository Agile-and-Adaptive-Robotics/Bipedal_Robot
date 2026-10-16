# migration_map.md — Phase-A restructure of ProofFinal into topic chapters

**Run date:** 2026-10-15. **Machine:** EB475WS4. **Actor:** Restructurer subagent (workflow).
**Directive:** advisor-approved architecture from the campaign ledger
(`comment_ledger.json` in this folder; the migration map = each entry's `oldSection → newHome`
fields) plus the approved chapter plan: `20-muscles, 30-torque, 40-ankle, 50-database
(before 60-sns), 60-sns, 70-overall, 80-conclusion`; `10-introduction` and `15-background`
stay; the Prior AnimatLab Walker section moves into `15-background`; appendices 92/93/94
unchanged.

**Method:** sections moved VERBATIM as exact line-range slices of the donor files (script
`tmp/disscampaign_20261015/restructure_phase_a.py`, run log `restructure_log.txt` in this
folder). Zero prose rewritten. The only composed text is (a) the deduplicated publication
notes (one per topic chapter, as directed) and (b) skeleton scaffolding (chapter headers,
`\graphicspath` lines, 50-database section headings, `%` provenance comments).

---

## 1. Section map — every old section → its new file

Provenance notation: `donor lines a–b` are 1-based inclusive line ranges of the donor file
**as it was before this restructure** (recoverable from git history if ever needed).

### From `chapters/20-methods.tex` (deleted)

| Old section (label) | Donor lines | New home |
|---|---|---|
| `\chapter{Materials and Methods}` + `\graphicspath` + publication note | 1–5 | header dropped; note deduplicated (see §3) |
| Overview (`sec:methods_overview`) | 7–11 | **70-overall** |
| Actuator Force Model (`sec:force_model`) | 13–68 | **20-muscles** |
| BPA Isometric Force Characterization Experiment (`sec:force_experiment`) | 70–87 | **20-muscles** |
| Robot Leg Model and Knee Test Stands (`sec:robot_model`) | 91–125 | **30-torque** |
| Knee Torque Data Collection and Equipment (`sec:torque_equipment`) | 128–137 | **30-torque** |
| Hybrid Torque Calculation Method (`sec:hybrid_method`) | 140–155 | **30-torque** |
| Improved Torque Model (`sec:improved_model`) | 159–367 | **30-torque** |
| Optimization Algorithms and Cost Functions (`sec:optimization`) | 369–385 | **30-torque** |
| Balance Platform and Inverted-Pendulum Testbeds (`sec:balance_testbeds`) | 387–404 | **40-ankle** |
| Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker (`sec:prior_walker`) | 406–412 | **15-background** (inserted at new line 152, between the locomotor-circuit figure of §`sec:bg_neural` and §`sec:bg_sns`) |
| Preliminary Neuromechanical Simulation Workflows (`sec:preliminary_sim_methods`) incl. subsections AnimatLab Body and Controller Organization, MuJoCo Conversion and Synthetic BPA Knee, Simulink SNS Library and Reflex Topology, Cross-Platform Verification (`sec:crossengine_methods`), Spiking-Neuron Mirrors (`sec:spiking_methods`) | 414–519 | **60-sns** |
| `%%%%%%%%%%` divider comment | 89 | dropped (cosmetic comment line only) |

### From `chapters/30-results.tex` (deleted)

| Old section (label) | Donor lines | New home |
|---|---|---|
| `\chapter{Results}` + `\graphicspath` + publication note | 1–5 | header dropped; note deduplicated (see §3) |
| Completed Foundation Studies (`sec:foundation`) | 7–9 | **70-overall** (rationale: it summarizes four foundation studies spanning muscles/torque/ankle/database topics — a program-level overview; see §7 flags) |
| BPA Force Characterization (`sec:results_force`) incl. subsections | 11–230 | **20-muscles** |
| `%%%%%%%%%%` divider comment | 231 | dropped (cosmetic comment line only) |
| Isometric Knee Torque (`sec:results_torque`) incl. subsections | 233–330 | **30-torque** |
| Follow-Up Identification and Route Redesign (`sec:ongoing`) incl. Redesigned-Route Results (`sec:route_redesign_results`) | 332–376 | **30-torque** |
| Balance Platform and Inverted-Pendulum Demonstrations (`sec:results_balance`) | 378–384 | **40-ankle** |
| Preliminary Neuromechanical Simulation Recordings (`sec:preliminary_sim_results`) incl. Cross-Engine Verification (`sec:crossengine_results`), Headless Validation (`sec:animatlab_validation`), BPA-Actuated Walker (`sec:bpa_walker`), Gait-Validation Sweep (`sec:gait_library_results`), Walker-Variant Comparison (`sec:walker_variant_results`) | 386–473 | **60-sns** |
| Spiking-Neuron Mirror Results (`sec:spiking_results`) incl. subsections (with both ZCODE fence comments) | 475–598 | **60-sns** |

### From `chapters/40-discussion.tex` (deleted)

| Old section (label) | Donor lines | New home |
|---|---|---|
| `\chapter{Discussion}` + `\graphicspath` + publication note | 1–5 | header dropped; note deduplicated (see §3) |
| Effects of Resting Length (`sec:disc_effects`) | 8–19 | **20-muscles** |
| Nondimensionalized Force Prediction (no label) | 22–30 | **20-muscles** |
| Comparison with Other Published BPA Models (no label) | 32–93 | **20-muscles** |
| Implications for Real-Time Control and Simulation (no label) | 96–101 | **20-muscles** |
| Biomimetic Artificial Muscle Analogy (`sec:disc_bio`) | 104–112 | **20-muscles** |
| Sources of Error in Isometric Torque Generation (`sec:disc_error`) | 117–135 | **30-torque** |
| Synthesis: From Actuator Characterization to Neuromechanical Control (`sec:synthesis`) | 138–144 | **70-overall** |
| What Spiking Changes and Does Not (`sec:disc_spiking`) | 146–160 | **60-sns** |
| What Mono- and Bi-Pulley Routing Buys (`sec:disc_pulley`) | 162–169 | **30-torque** |

### From `chapters/50-futurework.tex` (deleted)

| Old content (label) | Donor lines | New home |
|---|---|---|
| `\chapter{Future Work}` + `\graphicspath` | 1–2 | header dropped; `ch:futurework` re-attached to 70-overall (see §4) |
| Chapter-opening paragraph (controller/robot framing) | 4–58 | **70-overall** |
| A Neuromechanical Simulation Pipeline (`sec:fw_pipeline`) | | **70-overall** |
| A Two-Layer CPG with Sensory Feedback (`sec:fw_cpg`) | | **70-overall** |
| Extending the Sensorium (`sec:fw_senses`) | | **70-overall** |
| Physical Realization (`sec:fw_hardware`) | | **70-overall** |
| Toward Dynamic Whole-Body Maneuvers (`sec:fw_kickflip`) | | **70-overall** |

(Per ledger entry T11's `newHome: overall`.)

### From `chapters/60-conclusion.tex` (deleted)

| Old content | Donor lines | New home |
|---|---|---|
| `\chapter{Conclusion}` + `\graphicspath` | 1–2 | header rebuilt in **80-conclusion** (same title, same `ch:conclusion` label) |
| All four closing paragraphs | 4–11 | **80-conclusion** (verbatim) |

### New files with no verbatim donor section

- **50-database.tex** — SKELETON ONLY: `\chapter{The Sensory Afferent Database}\label{ch:database}` + `\graphicspath{{../Figures/50-database/}}` + four section headings (*A Structured, Model-Oriented Afferent Resource; Audited Corpus and Curation Coverage; AI-Assisted Literature Triage; Toward Robot and Simulation Use*) with `% Phase B: content.` placeholders. See §7 flag 1.
- **70-overall.tex** — skeleton + the moved sections listed above (methods Overview, foundation studies, synthesis, entire former Future Work chapter), under `\chapter{Overall System Integration and Future Work}\label{ch:overall}\label{ch:futurework}`. No composed prose.

---

## 2. Figure moves (MOVE, not copy; subfolder structure preserved)

**126 figure files moved in total** (20-muscles 18, 30-torque 22, 40-ankle 1, 60-sns 85 —
the 60-sns count includes the two wholesale subfolders). Full per-file list also in
`restructure_log.txt`. All moves are within
`Documentation/Reports and Papers/Dissertation/Figures/`. Nothing about any figure file was
modified — paths only.

**→ Figures/20-muscles/** (15 files + Deprecated)
- from `30-results/`: `02a_MaxForce10.eps` + `02a_MaxForce10-eps-converted-to.pdf`, `02b_MaxForce20.eps` + `…-eps-converted-to.pdf`, `03a_MaxForce10_all.eps` + `…pdf`, `03b_kmax10.eps` + `…pdf`, `04a_FStar10.eps` + `…pdf`, `04b_FStar20.eps` + `…pdf` (figs `fig:Fmax_fit`, `fig:F-P_contraction`, `fig:GenForceFit`)
- from `40-discussion/`: `05_Comparison.eps` + `05_Comparison-eps-converted-to.pdf` (fig `fig:ModelCompare`)
- from `20-methods/`: `01_testJigs1.pdf` (fig `fig:testJig` — the single-actuator force test stand)
- from `30-results/Deprecated/`: `02_Fmax_fit.eps`, `03_F-P_contraction.eps`, `04_GeneralizedForceFit.eps` → `20-muscles/Deprecated/`

**→ Figures/30-torque/** (19 files + Deprecated)
- from `20-methods/`: `testJigs_reviewed.pdf` (fig `fig:testJigs`), `KneeICR.eps` + `KneeICR-eps-converted-to.pdf` (fig `fig:ICR`), `bktFrame.pdf` (fig `fig:bktFrame`), `xiFrameTransform.pdf` + `.png` (fig `fig:xiFrameTransform`), `xiProjection.pdf` + `.png` + `.svg` + `.pptx` (fig `fig:xiProjection`)
- from `30-results/`: `FlxPin_group.eps` + `…pdf`, `ExtPin_group.eps` + `…pdf`, `Ext10mm_52cm.eps` + `…pdf`, `FlxGrp.eps` + `…pdf`, `KneeExt_10mm_passthru.png` (figs `fig:FlxPin_group`, `fig:ExtPin_group`, `fig:Ext10mm_52cm`, `fig:Ext_passthru`, `fig:FlxGrp`)
- Deprecated: `testJigs.pdf`, `testJigs_alternate1.pdf` (from `20-methods/Deprecated/`); `paired_coefficient_tables_2026-09-29.jpg` (from `30-results/Deprecated/`)

**→ Figures/40-ankle/**
- from `20-methods/`: `inverted_pendulum_photo.jpg` (fig `fig:balance_inverted_pendulum`)

**→ Figures/50-database/** — created, empty (skeleton chapter; no figures exist for it yet)

**→ Figures/60-sns/** (22 files + 3 subfolders)
- from `20-methods/`: `animatlab_body_gui.png`, `animatlab_network_gui.png`, `mujoco_converted_robot.png`, `mujoco_synthetic_knee.png`, `simulink_knee_reflex_arranged.png` + `.pdf`, `simulink_neuron_detail.png` + `.pdf`, `simulink_sns_library_grid.png` + `.pdf`, `mujoco_render_provenance.json` (provenance sidecar of the MuJoCo renders)
- from `20-methods/`: whole subfolder **`AnimatLab_reference/`** (9 files: `OurWalker.png`, `2layerImprove.png`, `RGnetwork.PNG`, `PFnetwork.PNG`, `sensoryMotor_Network.PNG`, `CPG_Output_MN.PNG`, `footMechano.png`, `fullNetwork.PNG`, `jointAngles.PNG`)
- from `30-results/`: `animatlab_phase1_joint_angles.pdf`, `sns_w2l_traces.png`, `spiking_mirror_rhythm.pdf`/`.png`/`_alt.txt`, `spiking_mirror_ground_gait.pdf`/`.png`/`_alt.txt`, `simulink_knee_spiking_twin.pdf`/`.png`/`_alt.txt`
- from `30-results/`: whole subfolder **`CPG_airstepping_figs/`** (49 files — unreferenced by any `\includegraphics` in ProofFinal; moved to keep it with its topic; nothing references its old path inside ProofFinal)
- Deprecated: `simulink_knee_reflex.pdf`/`.png`, `simulink_sns_library.pdf`/`.png` (from `20-methods/Deprecated/`); `animatlab_phase1_preliminary.pdf` (from `30-results/Deprecated/`)

**→ Figures/70-overall/**, **→ Figures/80-conclusion/** — created, empty (no figures in their moved content)

**Removed after emptying:** `Figures/20-methods/` (incl. its empty `Components/`), `Figures/30-results/` (incl. emptied `Deprecated/`), `Figures/40-discussion/` (incl. empty `Components/`), `Figures/50-futurework/`, `Figures/60-conclusion/`.
**Untouched:** `Figures/10-introduction/`, `Figures/15-background/` (incl. its `Components/` and `Deprecated/`), `Figures/92-AppendixA/`, `Figures/93-AppendixB/`, `Figures/94-AppendixC/`.

---

## 3. Publication-note deduplication (one per topic chapter)

The three near-identical notes (20-methods, 30-results, 40-discussion) were replaced by two:

- **20-muscles.tex** keeps the **bolen_2026** note: claims that `sec:force_model–sec:force_experiment`, `sec:results_force`, and `sec:disc_effects–sec:disc_bio` (with figures/tables) appeared in the *Actuators* 2026 paper; adds one clause pointing the knee-torque manuscript at Chapter `ch:torque` (needed so the verbatim "For both works, B.B. performed …" contribution sentence stays truthful); contribution sentences verbatim from the old note.
- **30-torque.tex** keeps the **bolen_isometric_2025** note: claims `sec:robot_model–sec:improved_model`, `sec:results_torque`, and `sec:disc_error` as the manuscript's subject (merging the three old notes' claims about that manuscript), plus "Author contributions for both works are stated in Chapter~\ref{ch:muscles}." (mirroring the old 30-results phrasing).
- The 30-results and 40-discussion notes were dropped as duplicates.

---

## 4. Label attachments (every old `\ref` keeps resolving)

New chapter labels: `ch:muscles`, `ch:torque`, `ch:ankle`, `ch:database`, `ch:sns`, `ch:overall`.
Old labels re-attached to moved content (as extra `\label`s on the new `\chapter` lines):

| Old label | Attached to | Why |
|---|---|---|
| `ch:methods` | 20-muscles (Ch. 3) | the old methods chapter opened with the actuator force model; keeps `Chapters~\ref{ch:methods}--\ref{ch:results}` ranges monotone (renders "3–4") |
| `ch:results` | 30-torque (Ch. 4) | largest results block; keeps ranges monotone |
| `ch:discussion` | 20-muscles | unreferenced anywhere; 5 of 9 old discussion sections moved here |
| `ch:futurework` | 70-overall (Ch. 8) | the former Future Work chapter moved here |
| `ch:conclusion` | 80-conclusion (Ch. 9) | same chapter, same title |

All 190 labels across the chapter set are unique (verified); zero unresolved `\ref`/`\eqref` targets in all chapters + main.tex (verified — see §6).

---

## 5. Advisor tracked changes T02–T07 (disposition)

All six were **already applied** to ProofFinal by the 2026-10-15 Overleaf sync
(`source_sync.md` §"Ports applied"), which pre-dates this restructure. The restructure moved
the containing text verbatim, so each TC is present in its new home — each re-verified by
exact-string search after the move:

| TC | Text (verbatim) | Now in |
|---|---|---|
| T02 | "first **improvement** of the torque tools" | 15-background (unchanged section §`sec:bg_bpa`) |
| T03+T04 | "**I** used these tools to convert" | 15-background (§`sec:bg_tools`) |
| T05 | "animals **also** carry real moral cost" | 15-background (§`sec:bg_ethics`) |
| T06 | "encountered by a humanoid" (deleted " walking") | **40-ankle** (§`sec:balance_testbeds` moved) |
| T07 | "simplifying the system" (was "increasing foot compliance and conserving energy") | **15-background** (§`sec:prior_walker` moved) |

---

## 6. Checks run (this session, exact commands)

1. **Conservation (verbatim-move proof)** — script `restructure_phase_a.py`, run
   2026-10-15: non-blank line multiset of the five donors (951 content lines, minus the
   deliberately dropped headers/notes/dividers) vs. the new files. Result: 0 donor lines
   missing (the 4 lines initially reported missing were the prior-walker lines landing in
   15-background, confirmed present at 15-background.tex:152–160 by `verify_phase_a.py`);
   38 added lines, all accounted for as chapter headers, `\graphicspath` lines, the two
   publication notes, and skeleton `%` comments/headings.
2. **Label uniqueness + `\ref` resolution** — `verify_phase_a.py`: 190 labels, 0 duplicates;
   0 unresolved `\ref`/`\eqref`/`\pageref` targets across all 18 chapter files + main.tex.
3. **Compile** — from `ProofFinal`:
   `pdflatex -output-directory=build_wf -interaction=nonstopmode main.tex` (exit 0),
   then `bibtex` from `build_wf` with absolute `BIBINPUTS`/`BSTINPUTS` pointing at
   ProofFinal (a first attempt with `BIBINPUTS=..\` silently opened the wrong thesis.bib
   and produced 130 missing-entry warnings — the absolute-path run produced a 56,953-byte
   `main.bbl` with **0** missing entries), then three more pdflatex passes (exit 0 each).
   Final `build_wf/main.log`: **0** undefined references, **0** undefined citations,
   **0** errors, no rerun warnings; `Output written` with **main.pdf = 10,794,410 bytes**.
   Built TOC: 1 Introduction, 2 Background, 3 Artificial Muscles, 4 Isometric Knee Torque,
   5 Inverted-Pendulum Balance Platform, 6 The Sensory Afferent Database, 7 Synthetic
   Nervous System Simulation, 8 Overall System Integration and Future Work, 9 Conclusion,
   Bibliography, Appendices — 50-database before 60-sns and 80-conclusion after 70-overall,
   as approved.
4. **Build-artifact hygiene** — stale `build_wf` aux/toc/lof/lot/out/log/pdf (including
   `.aux` files of the deleted chapters) were deleted before the fresh build so the
   reference check could not be satisfied by stale data.

Not run: no git commands (per campaign rules); no Overleaf-side changes; `upload/` does not
exist under `Dissertation/` (confirmed 2026-10-15 by `source_sync.md`, unchanged since).

---

## 7. Flags for Ben / Phase B

1. **50-database received no verbatim section.** The only SADb prose inside the deleted
   files is the second study of `sec:foundation` (now in 70-overall) — a single paragraph
   spanning four studies that cannot be split without breaking the verbatim rule. The
   substantive SADb description lives in `15-background` §`sec:bg_sns` (that chapter stays
   put this phase). Phase B: expand 50-database from that material; its four headings are
   placeholders, rename freely.
2. **Two references resolve but are one chapter off semantically** (deliberate trade-off
   to keep `Chapters~\ref{ch:methods}--\ref{ch:results}` ranges monotone as "3–4"):
   `15-background.tex:134` ("Chapter~\ref{ch:methods} discusses the algorithms used to
   optimize…") now renders Chapter 3 but the optimization section lives in Chapter 4;
   `15-background.tex:86` ("Chapter~\ref{ch:results} compares these models…") renders
   Chapter 4 but the model-comparison section lives in Chapter 3. Phase B: re-point to
   `ch:torque`/`ch:muscles` respectively.
3. **10-introduction §"Dissertation Organization" (lines 36–50) is now stale**: it
   hard-codes the old monograph structure (Chapters 2–7 as Background/Methods/Results/
   Discussion/Future Work/Conclusion). All its section `\ref`s still resolve; the chapter
   numbers and descriptions do not match the new topic structure. Phase B rewrite.
4. `ch:discussion` is attached to 20-muscles only to preserve it; nothing references it.
   Drop it during Phase B if unwanted.
5. The moved 50-futurework chapter-opening paragraph still says "This chapter lays out the
   planned work in four stages" — accurate inside Chapter 8, but Phase B may want to
   reframe it around the integrated chapter.
6. Two cosmetic `%%%%%%%%%%` divider comment lines (old 20-methods:89, 30-results:231)
   were dropped; every other comment (including all ZCODE campaign fences and FLAG-BEN
   notes) moved verbatim with its section.
7. `Figures/60-sns/CPG_airstepping_figs/` (49 files) is not referenced by any
   `\includegraphics` in ProofFinal — it moved wholesale to keep it with its topic. If any
   external draft (e.g., `Dissertation/CPG_spinal_section_draft.tex`) references its old
   `30-results/` path, that draft will need its path updated.
8. Build note for reproduction: bibtex must run from `build_wf` with
   `BIBINPUTS`/`BSTINPUTS` pointing at ProofFinal (relative `..\` failed silently this
   session; absolute paths worked).

---

## 8. Files touched

- **Created:** `ProofFinal/chapters/{20-muscles,30-torque,40-ankle,50-database,60-sns,70-overall,80-conclusion}.tex`;
  `Figures/{20-muscles,30-torque,40-ankle,60-sns}/` (+ their `Deprecated/` subfolders),
  empty `Figures/{50-database,70-overall,80-conclusion}/`; this log; `restructure_log.txt`.
- **Modified:** `ProofFinal/chapters/15-background.tex` (prior-walker insertion only),
  `ProofFinal/main.tex` (include order + structure comment only).
- **Deleted:** `ProofFinal/chapters/{20-methods,30-results,40-discussion,50-futurework,60-conclusion}.tex`;
  emptied `Figures/{20-methods,30-results,40-discussion,50-futurework,60-conclusion}/`;
  stale `build_wf/` artifacts (regenerated).
- **Not touched:** `upload/` (absent), all code/data/result files, the three appendices,
  10-introduction, front matter, `thesis.bib`, `psutex.cls`, `15-background` figures.

---

## 9. Audit round 1 addendum (2026-10-15, same session)

Audit verdicts: build errors **[]**, ref/figure problems **[]**; the re-anchor list held
one genuine document defect and a set of degenerate-anchor artifacts. Actions:

1. **T14 APPLIED (the one genuine defect).** "Regenerating the results files updates
   this table." existed twice in `94-AppendixC.tex` — inside the `\caption` of
   Table~`tab:app_xi_values` (old line 205) and at the end of the adopted-values
   paragraph (line 223). The ledger's T14 `oldSection` ("Adopted Correction-Term Values
   **(table captions)**") identifies the **caption** instance as the tracked deletion's
   anchor; source_sync flag 2 had assumed the paragraph copy was the target only because
   it did not know the caption copy existed (the audit corrected this). Fix: the caption
   instance was deleted verbatim per the ledger TC (caption now ends "…pick 32 at lines
   164--171.}"); the paragraph instance at line 223 is targeted by no ledger entry and
   was left verbatim — **flagged for Ben** (see below). Post-fix grep: exactly 1
   occurrence document-wide; the T14 anchor is now unique.
2. **Degenerate anchors (C01/C02/C04/C05/C08/C09/C18, T02/T04/T05/T06 quotes)** — no
   document defect; the one-word quotes are review-tool truncations/common words. The
   intended anchor paragraphs were re-verified present this round by exact-string grep:
   C01 anchor ("Artificial muscles are a key enabling technology…") at
   `10-introduction.tex:8`; C02 anchor ("Second, actuator force alone does not
   determine…") at `10-introduction.tex:22`; post-TC strings for T02 (15-background),
   T03+T04 ("I used these tools to convert", 15-background), T05 (15-background), T06
   (40-ankle.tex:8), T07 (15-background, prior-walker section).
3. **Blob entries (T01, T03, T07–T13, T15)** — audit-verified expected absences/
   concatenation artifacts; deleted fragments confirmed 0 hits, kept fragments confirmed
   present in their new homes. No action required or taken. One observation stood:
   T12's "predicts with the wrong trend" survives only inside a **commented** line
   (`02-abstract.tex:11`, begins `%First, …`) — not live prose; left untouched (removing
   comment text has no ledger basis).
4. **Rebuild after the T14 edit** — `pdflatex -output-directory=build_wf` ×2 (exit 0
   both): 0 undefined references/citations, 0 errors, no rerun warnings;
   `main.pdf` = 10,794,353 bytes.

**Round-1 flags for Ben (new):**
- The paragraph copy of the T14 sentence at `94-AppendixC.tex:223` ("… Regenerating the
  results files updates this table.", end of the adopted-values paragraph) remains: no
  ledger TC targets it (T14's oldSection is table captions). Delete it in the Overleaf
  review UI or a later phase if you want the maintenance note gone entirely.
- `%`-commented legacy text in `02-abstract.tex:11` still carries T12-deleted phrasing;
  harmless to the build, left verbatim.

---

## 10. Audit round 2 addendum (2026-10-15, same session)

Audit round 2: build errors **[]**, ref/figure problems **[]**, T14 fix confirmed
**FIXED**, and the 21 listed re-anchor items self-annotated as degenerate-quote or
expected-absence artifacts ("NOT document defects"). Independent verification this round
(script `tmp/disscampaign_20261015/audit_round2_verify.py`, run over all 18 chapter files
+ main.tex): every audit count reproduced exactly — multi anchors C01=11, C02=8, C04=182,
C05=2, C08=9, C09=4, C18=36, T02=9, T04=2, T05=31, T06=28 (all legitimate extra hits of
ordinary words/substrings, e.g. "Artificial" ×4 inside `\chapter{Artificial Muscles}`,
"begin" dominated by `\begin`); expected-absence anchors all 0 hits (T03, T07, T01-blob
rewritten lines, T08–T13/T15 deleted fragments); kept-content anchors all exactly-once in
their new homes (incl. `{bktFrame.pdf}` at 30-torque.tex:119 and the T11 fragments in
70-overall.tex:46/73); T14 sentence exactly 1 (94-AppendixC.tex:223); all six TC
post-states exactly-once at their recorded lines. **No fixes applicable** under
moves-and-path-fixes-only: making degenerate anchors unique would require deleting
legitimate prose; making expected absences present would revert applied TCs. Zero files
changed this round; one confirmation `pdflatex` pass (exit 0): 0 undefined
references/citations, 0 errors, no rerun warnings, PDF written. Standing flags from §9
remain open for Ben.

## 11. Audit round 3 addendum (2026-10-15, same session)

Round 3 re-ran all counts ("identical to round 2" per the audit header, which now
adjudicates **all 21** listed items as artifacts, NOT document defects). My independent
sweep (`audit_round2_verify.py`, re-run this round) reproduced every count exactly as in
§10 — multi anchors C01=11, C02=8, C04=182, C05=2, C08=9, C09=4, C18=36, T02=9, T04=2,
T05=31, T06=28; all expected-absence anchors 0; all kept-content anchors exactly-once at
the same lines; T14 = exactly 1 (94-AppendixC.tex:223). **No fixes applicable; zero
source files changed this round.** Confirmation `pdflatex` pass: exit 0, 0 undefined
references/citations, 0 errors, PDF written. The audit has reached a fixed point: the
remaining re-anchor entries are properties of the ledger's one-word/blob quote quality
(review-tool truncation), not of the document, and no move/path edit can or should
change them. Open items remain the two Ben-ruling flags of §9.
