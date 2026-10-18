# CONTENT_SPEC — defense deck update to the 2026-10-17 topic-chapter dissertation

Authoritative inputs, in priority order:
1. The dissertation chapters: `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\*.tex`
   (10-introduction, 15-background, 20-muscles, 30-torque, 50-database, 60-sns, 70-overall,
   80-conclusion; appendices 92/93/94). **Every number on a slide must come from a chapter, with a
   file:line anchor recorded in the evidence manifest.** Where the deck shows a result that is NOT
   in the dissertation (post-freeze work), label the slide "since submission" and cite the artifact
   (mat/log/report path) — allowed only in clearly-marked slides.
2. `Documentation\Project_Digests\digest_2026-10-15.md` (cross-check only; where it disagrees with
   a chapter, the chapter wins — e.g. the cam-redesign margins in 30-torque.tex:512 are
   +30.83%/+39.77% on the FINAL geometry, superseding the digest's +31.30/+42.34 winner-only
   numbers).
3. `spinal\reports_20261016\` (easteregg2, uncommitted) for the since-submission wean campaign.

## Hard deck rules (Ben's directives)

- Assertion/evidence style: the assertion title carries the claim; a small `cite=` source line
  carries the evidence (number + source). Every quantitative claim needs both.
- SLIDES STAY SPARSE: near-zero prose, one message per slide, figures carry the weight. If a slide
  needs a paragraph, it is two slides. Prefer figures/GIFs/diagrams over text cards.
- Speaker notes on every slide: the evidence detail, provenance, and one anticipated committee
  question with its answer.
- Chapter kickers replace the old "Aim 1/2/3" kickers — the dissertation now frames FOUR research
  objectives (10-introduction.tex, sec:objectives):
  1. Characterize isometric BPA force (Ch. muscles)  2. Validate/improve knee torque prediction
  (Ch. torque)  3. Consolidate the afferent literature into the SADb (Ch. database)
  4. Develop tools + framework for neuromechanical control (Ch. sns + overall).
  The introduction names FOUR gaps (read §Motivation/§Research Objectives for their exact
  wording); the old deck's three-gap slide must be updated to the dissertation's four.
- Standing prose rules: Oxford comma; "model" never "plant" for the robot/bodies; American
  spelling; diameter as engineering symbol where it appears (⌀ is fine on slides); no emoji.
- Palette/typography: keep deck_lib constants (BROWN assertion titles, Calibri, 13.33×7.5). Reuse
  `content_slide`, `pic`, `node`, `stat`, `bar_chart` helpers. Keep `pageno` on every content slide.
- Keep each slide's `build()` numbering consistent: `build(prs, start=N)` receives its first
  slide number from build_deck.py (which chains a→b→c→d), and every internal call passes
  `start + k`. Never hardcode absolute slide numbers — batches are edited in parallel.
  build_deck.py prints the final count.
- DO NOT modify anything outside `Defense\deck_builder\` and the deck output. Dissertation sources
  are read-only.

## New deck outline (target ~52 slides)

### slides_a.py — open + Ch. 3 muscles (published: Bolen et al. 2026, Actuators)
- s01 title: update date line to "November 2026" (defense 2026-11-06).
- s02 roadmap: FOUR objective cards (force → torque → database → SNS control), evidence tags per
  card from the chapters (e.g. 537 tests; meets/exceeds human torque; 943-paper corpus; verified
  toolchain). Mirror the intro's wording.
- s03-s06 project/motivation/BPA/gaps: keep, light edits; gaps slide → four gaps per the intro.
- s07-s14 (Ch. 3 muscles): keep the science and numbers (verified in the revision campaign);
  update kickers to "Ch. 3 · muscles" style; verify each number still matches 20-muscles.tex
  (RMSE/GoF tables, 537 = 321+216, 42.1% Festo error example, arctan coefficients).

### slides_b.py — Ch. 4 torque (manuscript bolen_isometric_2025 + redesign)
- s15-s24: keep (stands, hybrid, χ corrections, identification, flexor, extensor, values of
  record, biomimetic flexor/extensor). Update kickers + cites to Ch. 4 / manuscript; verify
  numbers against 30-torque.tex.
- s25 REWRITE as the redesign ladder (30-torque.tex ~lines 315-545):
  * Direct route: flexor +5.01% at −107°; extensor +5.02% worst margin (superseded) — line ~504.
  * Cam redesign (current design of record): min margin +30.83% vs required envelope / +39.77% vs
    plain human; pressures 269–507 kPa (VTUG floor 150 kPa); bend radius 96.9 mm; condyle
    clearance 43.1 mm at 10° extension; moment arm 59.4–102.8 mm, mean 72.4 mm — line ~512.
  * Mono-pulley: peak flexion torque 81.23 → 85.74 N·m (+4.51 N·m, 5.6%); off-axis ≤0.803 N·m
    (<1% of peak); human flexor scale 23.16 N·m (factors 3.51/3.70) — line ~542.
  * Bi-pulley screening: all five sites nonpositive/infeasible; clearances 4.26–19.14 mm vs 25 mm
    required — geometry binds, not torque.
  * Hardware: "hardware torque tests of the updated compliant system have begun with the extensor
    route" — if ExtTest_1_01 numbers (1.1630/7.8494/23.4367 N·m) are in the chapter cite it there,
    else a since-submission tag citing `Testing_Data\2026_06_Festo\Ext_20mm\ExtTest_1_01.mat`.
  Consider 2 slides (ladder overview + cam detail with route/winner figure if a legible one
  exists in Figures\30-torque or Mesh_Optimization Results; do NOT use print-FAILED figures —
  see BEN_TODO #7; the GIFs Knee_{Flexor,Extensor}_20mm.gif may be embedded; PowerPoint animates
  them live, the COM render shows frame 1 — only use if frame 1 is acceptable, else skip).

### slides_c.py — Ch. 5 database + Ch. 6 SNS + Ch. 7/8 close  (biggest rewrite)
Ch. 5 (2-3 slides):
- Structured afferent resource: audited corpus 943 papers (382 archived → 561 retained);
  reconciliation story (Zotero 538 / group 911 / Airtable 99 → 943); citation graph 730 papers,
  11,130 directed pairs, 18 Louvain clusters; working base 957 records (561/39/244/109/4).
  Retained-paper coverage: 523 notes, 426 animal-prep, 118 pathway-linked, 86 afferent-classified,
  258 translation notes. Figure: `Figures\50-database\sadb_explorer.png`.
- AI-assisted triage + toward robot/simulation use (numbers from the chapter).
Ch. 6 (6-8 slides; the old Aim-3 slides are raw material — s26 pipeline, s27 knee loop, s29 walker
lesson, s30 sensory, s31 verification are reshaped into the chapter's actual arc):
- Prior work: AnimatLab bipedal walker (figs in Figures\60-sns: animatlab_body_gui.png,
  animatlab_network_gui.png, animatlab_phase1_joint_angles.pdf).
- Toolchain: OpenSim→MuJoCo conversion + synthetic BPA knee + Simulink SNS library
  (mujoco_converted_robot.png, mujoco_synthetic_knee.png, simulink_knee_reflex_arranged.png).
- Contact-driven rhythm core: cross-engine verification (subsection exists twice — engine A/B),
  headless validation of the AnimatLab reference models.
- BPA-actuated walker in MuJoCo + gait-validation sweep vs the locomotion reference library.
  Assets: `Defense\slideshow\` s3k_mutfix_joints.png, s3k_mutfix_skeletal.gif, syn6_*, w2lvar_*.
- Stimulus and deletion experiments on the walker networks (the chapter's distinctive content —
  extract its actual experiments/numbers).
- Spiking-neuron mirrors: results + "what spiking changes and does not" (numbers from the chapter;
  slideshow assets bpa_spiking_* vs bpa_nonspiking_*, spiking_mirror_*,
  simulink_knee_spiking_twin.png; connectome render `bilateralrg_sns_diagram.png` fits here or
  the rhythm-core slide — user QC directive: circuit renderings get a correctness check).
Ch. 7/8 close:
- Synthesis chain slide (old s33, now 5 stages incl. database; cite Ch. 7).
- Honest boundaries (old s34; content now from Ch. 7 + chapter-level limitations; keep the
  bio-extensor FVU honesty card).
- Future work (old s35; now from 70-overall §Future Work subsections: closing the loop in
  simulation; sensorium, physical realization, whole-body maneuvers).
- Acknowledgments + thanks (keep; date/emails unchanged).

### slides_d.py — appendix (keep A-J shape)
- Update every chapter/appendix cross-ref to the new numbering (muscles = Ch. 3, torque = Ch. 4;
  appendices now 92/93/94). a10 SNS implementation: refresh with the Simulink spiking twin /
  library figures. Optional new appendix slide: walker-comparison table with per-model numbers
  ONLY if anchored in 60-sns.tex.

## Asset wiring

Already in deck_assets/figs: walker_arch, walker_body, rg/pf/sensorimotor/full network renders,
sns_library, sns_reflex, mujoco_robot, mujoco_knee, all Aim-1/2 figure PNGs.
New raster assets to copy into deck_assets/figs (or reference by absolute path):
`Defense\slideshow\*.png|*.gif` (16 files), `Figures\50-database\sadb_explorer.png`,
`Figures\60-sns\*.png`. PDFs in Figures\ need convert_figs.py (repointed; Ghostscript at
D:\Anaconda\envs\gs). GIFs: embed only if frame-1 acceptable (COM render check).

## Evidence manifests (required deliverable per batch)

Write `deck_builder\evidence_<a|b|c|d>.md`: one row per slide — slide #, assertion, every number
on the slide, and its source as `chapters/<file>.tex:<line>` (or artifact path for
since-submission items). The QC gate verifies a sample of these rows against the sources.

## Build + render (verified working)

- Build: `D:/Anaconda/python.exe build_deck.py` (base Anaconda has python-pptx 1.0.2 + PIL).
- Render: `D:/Anaconda/python.exe render_deck.py` (PowerPoint COM → deck_render/slide_NN.png,
  1600×900). GIFs render as frame 1.
- Output deck: `Defense\Bolen_Dissertation_Defense.pptx` — the up-to-date deck always lives at
  this plain name (Ben's ruling 2026-10-17); the pre-update deck stays archived at
  `Bolen_Dissertation_Defense_20260910.pptx`.
