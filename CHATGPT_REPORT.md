# ChatGPT report for ZCode — 2026-09-08

> Latest handoff: the **2026-09-09 dissertation neuromechanical figures** session
> is appended below. The compilation status in this older section is historical
> and does not verify the new local edits.

## Final handoff status

Completed: manuscript-status corrections, publication metadata consistency, 1-inch default left margin, explicit bar unit definition, and missing Festo institution fields. Seven Overleaf files were edited: `chapters/20-methods.tex`, `chapters/30-results.tex`, `chapters/40-discussion.tex`, `thesis.bib`, `BolenFrontiers22.bib`, `psutex.cls`, and `main.tex`.

The last observed Overleaf compilation had **zero errors and one warning**, the subfloat counter advisory. Both subfig and subfloat remain installed, as Ben wishes to retain multiple images within a figure. Detailed layout/typesetting messages still require a separate review. The last page count observed after the margin change was 105; no new page-count check was made after the unit/institution cleanup.

**Commit scope:** this session saved `D:/GitHub/Bipedal_Robot/CHATGPT_REPORT.md` locally but did not download or synchronize the edited Overleaf source files into the repository. Committing this report alone does not commit those LaTeX changes. ZCode should obtain the current Overleaf version before integrating dissertation sources. No commit or push was performed by ChatGPT.

**Unfinished:** wide-table and oversized-figure layout work was discussed but not implemented. Ben asked whether wide tables should use rotated pages with minimal interruption to the text; that question remains to be addressed. No tables or figures were rotated, resized, or split in this session. A complete PSU formatting audit and verification of published author contributions were not performed.

The sections below preserve the specific changes and verification history. Earlier warning counts are intermediate observations, superseded by the final count above.

## Work completed

Ben requested correction of false claims that the knee-torque manuscript had been submitted to Frontiers. Ben explicitly confirmed that it remains in preparation and has not been submitted; separately, he confirmed that the actuator force-characterization paper is published.

Edits were made directly in the Overleaf project **Bolen-Dissertation**, https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb . Repository dissertation source files were not synchronized or changed. The standalone journal-paper Overleaf project was not edited.

## Files created or modified

Overleaf paths:

- `chapters/20-methods.tex`: changed the knee-torque publication note from “appeared in” and “submitted to Frontiers in Robotics and AI” to an unpublished manuscript in preparation, explicitly stating it has not been submitted. Preserved the existing author list and contribution statement.
- `chapters/30-results.tex`: equivalent correction in the Results note.
- `chapters/40-discussion.tex`: equivalent correction in the Discussion note.
- In those three notes, replaced `\singlespacing\small\itshape` with `\doublespacing\normalsize\upshape`; main.tex specifies 12pt and the class loads setspace. No heading pages were added for the unpublished manuscript.
- `thesis.bib`: retained key `bolen_isometric_2025`, title, authors, year, and `@unpublished` type. Changed note from `Manuscript submitted for publication` to `Unpublished manuscript in preparation; not submitted for publication`.
- `BolenFrontiers22.bib`: after Ben confirmed publication, replaced the stale `@unpublished{bolen_2026,...}` entry with the existing published `@article` entry from thesis.bib: Actuators 15(5), article 230 (2026), DOI 10.3390/act15050230. This is the actuator paper, not the knee-torque manuscript. Both bibliography replacements were checked by exact full-editor text comparison to their intended content.

Local files:

- `C:/Users/Ben Bolen/.codex/.chatgpt-projects/g-p-681bb99d9c088191ad86129c8548cbd9/PROJECT_INSTRUCTIONS.md`: created earlier in this task as reusable custom-instruction text. It is not an installed project-settings change.
- `D:/GitHub/Bipedal_Robot/CHATGPT_REPORT.md`: this final report; a staging copy also exists in the local ChatGPT project workspace.

## Model changes and runs

No model constants, transforms, stiffness ordering, bounds, flags, scientific results, or MATLAB files changed. No optimizations ran. No commits or pushes.

Overleaf automatically recompiled during edits. After the initial publication-note correction, the PDF had 113 pages (initial preview: 112), zero compilation errors, six warnings, and 23 informational/typesetting messages. Subsequent results are recorded below and in the final status above. Per-compile wall times were not measured. Corrected note content was verified in the editor and compiled preview; all three chapter notes were visually inspected.

## Formatting guidance and limits

Read PSU's live requirements:

- https://www.pdx.edu/gradschool/etd-formatting-requirements
- https://www.pdx.edu/gradschool/chapter-heading-page

PSU recommends, but does not require, a Chapter Heading Page for unpublished co-authored material. Ben specifically emphasized this distinction. Retaining corrected status text in the existing notes does not mean a heading page is mandatory for the knee-torque manuscript.

This was a targeted status correction, not a full ETD compliance audit. The published Actuators material has separate chapter-heading/citation/contribution requirements, and the existing combined notes should not be represented as having passed a complete audit of those requirements. Existing contribution wording was not independently verified against the published paper. DOI and publisher retrieval attempts returned HTTP 429; published bibliographic metadata came from the existing thesis.bib entry, with publication status confirmed by Ben.

At the initial publication-note verification, diagnostics included a Methods float too large by 16.17862pt, Results overfull vertical boxes, other overfull/underfull boxes, subfloat and deprecated-unit warnings, and missing institution fields for Festo bibliography entries. The bar and institution warnings were later fixed. No diagnostics were reported at the modified note lines. The six-warning count was already present on initial inspection; the individual pre-edit diagnostics were not captured, so do not claim all typesetting messages were proven pre-existing.

## Follow-up for ZCode

Pull or export the current Overleaf sources through Ben's established workflow before making local dissertation changes; this session did not synchronize the repository mirror. Preserve the corrected manuscript status. Review broader ETD formatting separately if requested. The knee manuscript's existing 2025 citation year was retained as draft metadata, not a publication date.

## Subsequent requested margin change — 2026-09-08

Ben authorized changing the class left margin to 1 inch under the updated PSU minimum.

- Modified Overleaf `psutex.cls` only. Added a `1in` option with `@lmarwidth=0in`, `@smarwidth=0in`, and `@texwidth=6.5in` (LaTeX's default origin adds the physical first inch).
- Changed default options from `12pt,1.7in,double,submission` to `12pt,1in,double,submission`. The old effective left margin was 1.7 inches, not 1.5; text width increased from 5.8 to 6.5 inches, retaining a 1-inch right margin on US Letter. Explicit legacy `1.5in` and `1.7in` options remain available. Updated corresponding class comments.
- No vertical layout, font size, scientific content, or MATLAB changes. Verified the entire saved class source exactly equals the intended edit.
- Automatic compilation completed with zero errors, six warnings, and 12 informational/typesetting messages; document reflowed from 113 to 105 pages. Visually inspected a body page with symmetric 1-inch left/right margins.
- Reflow changes float layout: the existing Methods oversized-float warning at line 101 now reports 54.8025pt instead of 16.17862pt. Other overflow messages remain, so the document still needs a broader float/layout review before final submission. No claim of complete ETD compliance.
- The Results note was also visually verified after the previous report was first saved.

## Unit and institution cleanup — 2026-09-08

- Ben requested an explicit bar definition and wished to retain multi-image figure support.
- Overleaf `main.tex`: inserted `\DeclareSIUnit\bar{bar}` immediately after loading siunitx. Left subfig/subfloat packages and counters unchanged.
- Overleaf `thesis.bib`: added `institution = {Festo AG \& Co. KG}` to `festo_2013`, `festo_2018`, and `festo_2026`, matching the corporate author already recorded in each entry. No citation dates, authors, URLs, or titles changed.
- Read fresh editor contents before both edits and verified the entire resulting files against intended replacements.
- Latest compilation UI shows zero errors and one warning, the subfloat counter advisory. Deprecated bar and missing-institution warnings are gone. No figure-size/layout edits were performed in this cleanup; do not attribute other diagnostic changes to it. Typesetting messages were not audited in this pass.

## Handoff — dissertation neuromechanical figures, 2026-09-09

Ben requested this handoff and ended the session. Work was performed locally on
DESKTOP-5Q16KE9 under `C:/Users/Ben/Documents/GitHub/Bipedal_Robot`.

### Completed scope

Added a preliminary simulation Methods section covering the assembled
OpenSim-to-MuJoCo conversion workflow, Python quasi-static BPA coupling,
reciprocal-inhibition SNS on the separate synthetic antagonist knee,
seven-block Simulink SNS library/reduced-order reflex topology, and the saved
10-second AnimatLab phase-1 recording. Following Ben's feedback, moved this
section from the overview to the end of Methods, immediately after the historical
AnimatLab walker section. The actuator-force section remains first in the
technical Methods sequence.

Expanded the visual material to show the AnimatLab body, hierarchy and circuit
schematics, actual AnimatLab GUI views, native MuJoCo model renders, and native
Simulink implementation diagrams. Replaced the included neural/contact diagnostic
Results plot with bilateral hip, knee, and ankle angles. The diagnostic remains
available as a supporting asset. All new results remain explicitly preliminary.

### Files changed or created

All paths below are relative to `Documentation/Reports and Papers/Dissertation/`
unless otherwise stated.

- `ProofFinal/chapters/20-methods.tex`: new section
  `sec:preliminary_sim_methods` with AnimatLab, MuJoCo, and Simulink subsections;
  six Methods figure environments, including a sideways full reflex diagram.
- `ProofFinal/chapters/30-results.tex`: section `sec:preliminary_sim_results` and
  figure `fig:animatlab_phase1_preliminary`, now using
  `figs/Preliminary/animatlab_phase1_joint_angles.pdf`.
- `ProofFinal/figs/Preliminary/AnimatLab_reference/`: unchanged copies of all nine
  images Ben supplied. Included architecture assets: `OurWalker.png`,
  `2layerImprove.png`, `RGnetwork.PNG`, `PFnetwork.PNG`, and
  `sensoryMotor_Network.PNG`. `jointAngles.PNG`, `CPG_Output_MN.PNG`,
  `footMechano.png`, and `fullNetwork.PNG` remain reference material; the older
  supplied plots were not relabeled as new run output. Exact historical image
  dates were not established.
- `ProofFinal/figs/Preliminary/animatlab_body_gui.png` and
  `animatlab_network_gui.png`: real captures of the saved
  `Neuromechanical_Models/Biped_2xCPG_wSubs/Biped_2xCPG_wSubs.aproj` in AnimatLab Pro.
- `ProofFinal/figs/Preliminary/mujoco_converted_robot.png` and
  `mujoco_synthetic_knee.png`: native MuJoCo renders, with settings and source
  hashes in `mujoco_render_provenance.json` beside them.
- `ProofFinal/figs/Preliminary/simulink_sns_library_grid.png`,
  `simulink_neuron_detail.png`, and `simulink_knee_reflex_arranged.png`: included
  high-resolution native Simulink exports. Vector PDF equivalents and the
  original unarranged library/reflex exports are retained beside them. PNGs were
  chosen for inclusion because PDF printouts carry excess page margins and
  unembedded-font dependencies.
- `ProofFinal/figs/Preliminary/animatlab_phase1_joint_angles.pdf`: included Results
  figure. `animatlab_phase1_preliminary.pdf`: original supporting diagnostic.
- `Notes/build_preliminary_animatlab_figure.py` and
  `preliminary_animatlab_verification.json`: reproduce and audit both phase-1 plots.
- `Notes/export_mujoco_models.py`: static native render reproduction; optional
  `--viewer` opens the native viewer without stepping.
- `Notes/export_preliminary_simulink_figures.m`, `simulink_figure_export.log`,
  `simulink_figure_provenance.md`, `simulink_model_hashes_before.json`, and
  `simulink_model_hash_verification.json`: native export workflow and verification.
- `Notes/build_neuromechanical_review.py` and `neuromechanical_figure_review.pdf`:
  seven-page figure-review packet, **not a compiled dissertation**.
- `Notes/neuromechanical_figure_integration.md`: integration/provenance detail.
- Repo-root `CHATGPT_HANDOFF.md` and `CHATGPT_REPORT.md`: this handoff, preserving
  earlier session material.

### Verification and model state

No Xi values, bracket points, transforms, stiffness ordering, bounds, or solver
flags were changed. No conversion, optimization, or new dynamic simulations ran.
Total wall-clock time was not recorded; do not infer runtime performance from
this figure-preparation session.

AnimatLab plotting used the existing `DataTool_7.txt` and `DataTool_8.txt` in
`Neuromechanical_Models/Biped_2xCPG_wSubs/`, associated with
`walk new new tester added 2 axis_phase1.asim`. Both exports contain 50,010 rows,
with aligned, uniformly increasing 0.0002-second timestamps and finite values.
The configured chart interval is 0 through 10 seconds: retained 50,001 samples
and excluded nine trailing all-signal-zero padding rows after 10 seconds.
Within-window zero-voltage samples were retained. Joint channels are explicitly
`JointRotationDeg`; no sign reversal, smoothing, or resampling was applied.
Contact-labelled channels are sensory-neuron membrane voltages, not measured
contact forces or Boolean contact states. The model includes an external
horizontal force during the first second. Body-position channels were recorded
but not plotted because their unit convention was not resolved.

MuJoCo 2.3.7 renders load the existing converted robot at keyframe 0 and the
separate demo knee at an illustrative static -35-degree angle. Camera, lighting,
color, visual tendon width, and knee transparency were changed in memory only.
The source geometry and physics files were not edited; no simulation steps ran.

MATLAB R2025b exported the actual saved Simulink blocks. Library grid positioning
and full-model auto-layout were display-only changes in memory. Every input
port's source block/port was compared before and after arrangement. Both saved
`.slx` files remained byte-identical according to SHA-256 checks. The neuron
interior was re-exported after opening it to populate its native primitive icons.

Visually inspected all seven review pages and the individual figures; enlarged
the Simulink contents by using the tightly bounded PNG exports. Verified all
13 included image paths, nine unique new labels, and their references.
Targeted `git diff --check` passed. Final LaTeX pagination is not verified.

### Claim limits — preserve these

Do not claim SNS control of the converted Gait2392 robot, validated Xi-corrected
MuJoCo torque, a completed CAD-coupled Simscape Multibody plant, stable AnimatLab
walking/effective ground-contact gating, or quantitatively valid Simulink knee
regulation from the current saved run. The diagrams document implementation;
the phase-1 plot documents recorded joint motion.

### Remaining work and integration

1. Compare the edited local chapters with the current Overleaf project
   `Bolen-Dissertation` before copying, to preserve changes made in other sessions.
   Transfer the new sections and referenced `figs/Preliminary/` assets with their
   relative paths. Compile in Overleaf and inspect float placement, caption fit,
   and the sideways topology page. No new packages were introduced.
2. The current review PDF is a separate preview. The full dissertation PDF,
   `Bolen_Dissertation.zip`, and online project were not updated in this session.
3. A genuine MuJoCo viewer screenshot was attempted. The viewer opened, but the
   installed desktop capture tool hung. No GUI PNG was saved; the included
   images are native renders. The viewer subsequently exited cleanly without
   stepping. No active viewer remains from that attempt.
4. Preserve existing/unrelated dirty files. Other sessions already had MATLAB
   source/model/results, AGENTS, dissertation ZIP, dedication, future-work, and
   appendix changes. This session did not produce those scientific changes.
   Python imports also changed two tracked MuJoCo `__pycache__` files; these are
   runtime byproducts, not intended dissertation changes.

No git commit, amend, merge, or push was performed. Ben handles commit/push in
GitHub Desktop. This session is finished at Ben's instruction; the Overleaf
integration above is a handoff, not an active continuation request.

## 2026-09-17 — MuJoCo/SNS model, plot, and dissertation-draft verification

### Files created or modified

- `Code/MuJoCo_SNS/spinal/runner.py`: fixed the built-in summary/PNG to
  use named neural channels, restrict rhythm metrics to the true walk
  window, and convert joint radians to degrees before plotting.
- `Code/MuJoCo_SNS/spinal/draw_circuit.py`: failed curriculum sentinel
  results are no longer labeled as winners; current stage-1 rhythm gains
  are combined with explicitly representative feedback gains; the
  compiled-edge contract now rejects both missing and extra edge groups;
  the V3 target/label was corrected from contralateral RG-E to InE.
- `Code/MuJoCo_SNS/spinal/_render_panels.py`: removed the retired
  DRIVE->PF edge and stopped applying the failed stage-3 `score=-100`
  parameter file as a tuned winner.
- `Code/MuJoCo_SNS/spinal/spinal_layers.py` and `_panels_check.py`:
  removed the retired DRIVE->PF renderer input port and added a regression
  assertion that it remains absent.
- `Code/MuJoCo_SNS/spinal/_figure_hindlimb_style.py`: resolves RG/PF
  channels by `neuro_names` rather than hard-coded column indices.
- `Code/MuJoCo_SNS/spinal/_live_plant_audit.py`: new reusable check of
  the plant after `runner.apply_harness`, including the eight Fmax repairs
  and a stock-OpenSim tolerance summary.
- `Code/MuJoCo_SNS/spinal/DESIGN.md`: added the verified 2026-09-17
  status at the top.
- `Documentation/Reports and Papers/Dissertation/CPG_spinal_section_draft.tex`:
  corrected curriculum status, air-run metrics, and Fmax fidelity claims.
- Regenerated `Code/MuJoCo_SNS/spinal/figures/circuit_dengstyle.{pdf,svg,png}`,
  `sns_diagram_panels.png`, `sns_layer_{rg,pf,motor}.png`, and
  `hindlimb_style_nap_air.{pdf,png}`; synchronized the dissertation copies
  under `Documentation/Reports and Papers/Dissertation/CPG_airstepping_figs/`.
- `Code/MuJoCo_SNS/spinal/spinal_run.npz/.png` now contain the successful
  NaP air smoke run made during this verification.

### Model changes

No physical model constants, muscle routes, joint ranges, bracket points,
transform conventions, Xi values, optimizer bounds, or solver settings were
changed. Changes are audit, reporting, and figure-generation corrections only.

### Runs performed

All Python commands used
`C:\Users\Ben Bolen\.conda\envs\myo\python.exe` from
`Code/MuJoCo_SNS/spinal`.

- `_fix_check.py`: PASS (heel edges single, Ia->IaIN present, central
  pathways present).
- `_panels_check.py`: PASS (NaP RG, laminated RG/PF, IaIN and mutual RC
  topology).
- `audit_signs.py`: PASS for every listed muscle/joint anchor on both sides.
- `_live_plant_audit.py`: PASS for all eight live Fmax repairs; tolerance
  exceptions reported.
- `draw_circuit.py --which deng --source best --fmt pdf,svg,png`: PASS;
  compiled missing-edge set and extra-edge set both empty.
- `_render_panels.py`: PASS with the failed stage-3 sentinel ignored.
- `_smoke_nap.py`: PASS, approximately 85 s wall time; full patched
  MuJoCo--SNS loop at 2 ms, finite output, five RG rises, no crash.
- `_figure_hindlimb_style.py nap_air_walk.npz nap_air`: PASS.
- `plot_run.py spinal_run.npz --side r`: PASS and located
  `subject01_walk1_ik.mot`; temporary plot-suite PNGs were removed after
  verification.
- Python syntax compilation passed for all edited/new Python files.

### Numbers and conclusions

- Environment: Python 3.10.21, NumPy 1.22.4, SciPy 1.9.3, MuJoCo 2.3.7,
  SNS-Toolbox 1.5.2.
- Smoke run: 5 RG bursts per leg in the 5--15 s walk window, period
  2.11 s (0.47 Hz), E-duty 0.67, right knee -106.3..+11.6 deg, right hip
  -20.1..+50.9 deg; finite throughout.
- Dissertation air NPZ (`nap_air_walk.npz`): four complete cycle
  intervals/five onsets, mean period 2.108 s, cycle duty 0.691, mean-cycle
  knee -106.1..+10.0 deg, hip -19.9..+50.0 deg.
- Curriculum stage 1 remains valid (score 137.254, trial 78). Stage 2 and
  stage 3 files both contain the `-100` no-countable-cycles sentinel and
  are diagnostic records, not converged winners.
- Live Fmax repairs are exact: ercspn 2500 N, intobl/extobl 900 N,
  ext_hal 162 N bilaterally. Of 86 non-pruned active actuators, 78 are
  within 15% of stock; eight exceed 15%.
- The circuit diagram is now structurally tied to the compiled network in
  both directions. The previous figure's V3->RG-E target was inaccurate;
  the implementation and corrected figure use V3->contralateral InE.

### Unfinished business

- Re-run the corrected stage-2 and stage-3 curriculum campaigns after
  deciding whether to retain/purge the all-sentinel studies. No optimizer
  or long curriculum run was launched in this session.
- `opensim_overlay_gait_cycles.png` and a new supported-ground figure must
  wait for a valid NaP stage-3 winner; the draft now says so explicitly.
- Compile the dissertation after integration and visually inspect final
  float placement. The figure generators and PDFs were validated, but a
  full LaTeX build was outside this request.
- Tracked and untracked `__pycache__` byproducts created by the verification
  runs were restored/removed and are not part of this handoff.

No commit or push was performed.

## 2026-09-17 — Follow-up visual review and literature comparison

Ben asked to continue the plot review and supplied the historical AnimatLab
figures plus recent Zotero PDFs. Work was stopped here at Ben's warning to
preserve credits; this section records the exact state.

### Additional figure changes completed

- `_figure_hindlimb_style.py`: the air panel now ends at 15 s instead of
  20 s (`6--15 s` shown). The previous final five seconds were the commanded
  DRIVE ramp-down/holding phase, but without a phase annotation they looked
  like an oscillator failure. The corrected figure shows only the settled
  commanded walk interval. `hindlimb_style_nap_air.{png,pdf}` was regenerated
  and the source/dissertation copies are SHA-256 identical.
- `_render_panels.py`: replaced the unreadable single horizontal strip with a
  two-row dissertation-scale composition: A/B = RG and PF side-by-side;
  C = the wider right-knee motor/reflex circuit. Added readable subpanel
  labels, a wrapped provenance note, and `sns_diagram_panels.pdf` output.
- `_render_panels.py`: removed duplicate visual input ports. The RG panel now
  shows the actual compiled DRIVE and POSTURE populations rather than those
  populations plus the reusable subnetwork's direct composition ports. The PF
  panel keeps one RG-E and one RG-F input and labels their live `rg_to_pf`
  gain. The motor panel distinguishes the tonic POSTURE population from the
  per-muscle `POST_i` inputs.
- `_panels_check.py` was rerun after these changes and all NaP/laminated RG,
  laminated PF, IaIN, mutual-Renshaw, and Ib assertions still pass.

The newly recomposed `Code/MuJoCo_SNS/spinal/figures/sns_diagram_panels.*`
files and their dissertation copies are current and SHA-256 identical. A
transient Windows preview lock initially blocked replacement; after it cleared,
all final source/dissertation figure pairs synchronized successfully.

### Visual findings

- The corrected hindlimb figure visibly shows four complete settled cycles
  (five onsets) from 6--15 s. RG/PF switching, knee-extensor release, deep knee
  flexion, and hip reversal occur in the expected sequence. The flat post-walk
  phase is no longer mixed into the gait figure.
- `circuit_dengstyle` is structurally valid and appropriate as a detailed,
  zoomable audit/reference schematic, but its parameter annotations are too
  small for it to be the only conceptual figure at ordinary print scale.
- The recomposed toolbox A/B/C figure is now readable at dissertation scale
  and is the better explanatory bridge. It remains a representative right-side
  RG/PF plus right-knee column, not a claim that every one of the 92 muscle
  columns is drawn.
- The older `fullNetwork.PNG` has the same reduction-density problem as the
  detailed current schematic. The older `RGnetwork.PNG`, `PFnetwork.PNG`, and
  `sensoryMotor_Network.PNG` succeed because they use symmetry, hierarchy, and
  only one mechanism per panel. The stacked historical activity plots support
  the organization of the new hindlimb trace figure.

### Reference-PDF review completed

Read-only visual review used local Poppler thumbnails and high-resolution page
renders. No Zotero files were modified.

- Shinohara et al. 2025 Fig. 1 is the closest conceptual precedent: bilateral
  RG -> PF -> MN -> musculoskeletal layers, centered C1/V3 commissurals, and
  sensory paths returning upward.
- Shevtsova et al. 2026 Fig. 2 and Rybak et al. 2024/2025 show that a dense
  full-page circuit can work when it is bilaterally symmetric, edge text is
  minimized, and connection type is encoded by arrowheads/circles and color.
- Deng et al. 2019 Fig. 2 supports the two-level RG/PF organization plus a
  separate detailed joint-controller inset.
- Jankowska's Ib/II schematic supports grouping sensory populations by spinal
  layer and explicitly separating ipsilateral and contralateral outputs.
- Rahmati and Klishko were reviewed mainly for kinematic/synergy plot style;
  they do not supply a better whole-network schematic template.

The modern references consistently map extensor populations to blue and
flexor populations to red/orange. `_figure_hindlimb_style.py`,
`spinal_layers.py`, and `draw_circuit.py` now use that convention consistently
with color-blind-safe blue and vermillion/orange. This is presentation-only;
labels, channels, weights, and topology did not change.

Di Russo et al. 2023 was reviewed after the credit reset: Fig. 2 separates the
system-level controller, Fig. 3 isolates phase primitives, Fig. 4 presents one
reflex motif at a time, and Fig. 5 draws one antagonist pair. This strongly
supports retaining `circuit_dengstyle` as the audit figure while using the
toolbox A/B/C panels as the readable conceptual explanation.

### Cleanup and final validation

- Temporary PDF thumbnails/contact sheets and `_preview_*`/`_final_*` visual
  QA images were removed. Zotero originals were read-only and unchanged.
- Added a render-beside-and-atomically-replace helper to
  `_figure_hindlimb_style.py` so transient Windows preview locks do not corrupt
  or partially overwrite the final PNG/PDF.
- Final source/dissertation copies of `circuit_dengstyle.{png,pdf,svg}`,
  `sns_diagram_panels.{png,pdf}`, and
  `hindlimb_style_nap_air.{png,pdf}` are SHA-256 identical.
- PNG decode/dimensions: circuit 5482x5559, panels 2828x2628, hindlimb
  1800x1639. All three PDFs are valid single-page files. Python syntax checks
  pass for all four affected generators; `_panels_check.py` still passes.

No model constants, network weights, optimizer studies, git commits, or pushes
were changed/performed during this follow-up.

### Draft figure order

`Documentation/Reports and Papers/Dissertation/CPG_spinal_section_draft.tex`
now introduces the architecture with `circuit_literature.pdf`, a bilateral
RG--PF--MN--muscle hierarchy plus one complete right-knee antagonist/reflex
motif. The toolbox-rendered `sns_diagram_panels.pdf` remains a software-native
artifact but is no longer the primary explanatory figure. The dense
`circuit_dengstyle.pdf` follows as a dedicated-page compiled-edge audit rather
than serving as the reader's first and only circuit explanation. Cross-references
and captions distinguish what is omitted for legibility from what is verified
against the compiled network.

### User-identified connection defects and correction

Ben's annotated review correctly identified real drawing defects in
`circuit_dengstyle`: both MN-to-activation lines began at hard-coded points
outside the MN glyphs; the flexor activation map had no output to the flexor
muscle; and the flexor muscle had no L/Ldot/F encoder paths back to its
Ia/II/Ib afferents. Additional raw-coordinate starts were found in the
heel/toe-to-InE and PF-F1-to-IaIN paths. The old compiled-edge assertion did
not catch these because it covered only SNS synapses and collapsed them into
population classes; MuJoCo plant/conversion paths are not in
`net.net.connections`.

Corrections:

- all identified raw-coordinate starts now originate at real glyphs or at
  explicitly connected routing buses;
- both knee antagonists now have complete MN--activation--muscle and
  muscle--Ia/II/Ib paths;
- partial hip MN glyphs were removed rather than left without plant/reflex
  connections;
- a separate plant-interface contract now requires all displayed activation
  and encoder paths, alongside the existing compiled-SNS edge-class contract;
- `draw_literature_circuit.py` generates the new literature-facing bilateral
  overview and complete knee motif. Every solid neural edge it displays is
  asserted against a freshly compiled representative network.

A second pixel-level review of the new figure caught three initially omitted
implemented relationships before promotion: RC->IaIN recurrent disinhibition,
RG-E stance gating of IB-EXC, and the overview's sensory projections to RG as
well as PF. All are now drawn and included in the neural-edge contract.

Both corrected figures were regenerated in the `myo` environment, inspected
pixel-by-pixel through review-scale renders, and synchronized to the
dissertation folder. `_panels_check.py` passes, both PDFs are valid one-page
files, and all PNG/PDF/SVG source/dissertation pairs are SHA-256 identical.

## 2026-09-17 — KINH and the four-PF/synergy question

### What KINH is

`KINH` is a project-local name for a **conditional swing-phase extensor
inhibitory interneuron**, not a recognized anatomical interneuron class. When
enabled, PF-F1 excites KINH with conductance 1.5; KINH then inhibits every
primary knee-extensor MN pool on that side through
`G["f1_kneext_inh"]`. The same cell can optionally inhibit ankle
plantarflexor pools through `f1_anklepf_inh`. It was introduced as an
engineering response to extension-dominant swing: suppress quadriceps during
the F1 window so the knee can flex. It is functionally compatible with a
flexor-phase inhibitory pathway but must not be presented as a specifically
identified V-class or literature-named population.

The topology is conditional: `params.py` defaults `f1_kneext_inh` to 0, in
which case KINH is absent. The figure audit deliberately builds a
representative network with gain 0.59 so the available pathway is visible.
The diagrams now label KINH as conditional.

### Does four PF populations equal four validated muscle synergies?

**No—not in the current implementation.** Ben's intended architecture is one
PF population per muscle synergy. In that architecture, each PF activity is
the synergy's temporal recruitment coefficient and its PF-to-MN weight vector
is the synergy's spatial muscle weighting. The current E1/E2/F1/F2 system was
constructed in the opposite order: four phase labels were chosen a priori,
their waveforms were generated from two RG drives, and `fit_pf.py` then used
those four fixed columns in NNLS to approximate functional-group activation
profiles. The four PF populations were not obtained from four NMF synergy
components, and their PF-to-MN weights are not the NMF component weights.

The count **four is plausible as a minimum**, but is not uniquely established:

- Literature commonly reports four or five human walking modules/temporal
  components, depending on muscles, preprocessing, task, and factorization.
  Four therefore sounds reasonable, but literature prevalence is not a
  validation of these four particular cells.
- The available SO backsolve contains only 126 frames from 0.49--2.49 s
  (about 1.6 gait cycles). A new per-leg NMF audit gives VAF 0.919 right and
  0.925 left at four components; five gives 0.949 and 0.952. Thus four clears
  a conventional 90% VAF threshold, but there is no sharp elbow at four and
  information criteria do not select four consistently.
- The older combined-leg NMF was not a valid PF-count test because
  left-versus-right phase consumed components. Its rerun gives VAF
  0.819/0.886/0.921/0.948 for 3/4/5/6 components and a continuously improving
  BIC-like score through eight.
- Most importantly, `_pf_basis_audit.py` shows that the four actual PF
  waveforms are effectively only **two temporal degrees of freedom**:
  corr(E1,E2)=0.998 and corr(F1,F2)=0.990; singular-value energy fractions are
  0.719, 0.280, 0.001, and approximately 0; the four-column condition number
  is 55.68. E1/E2 share the same peak and half-maximum duty, as do F1/F2.
  The matrix is numerically rank four but practically rank two, so individual
  E1-versus-E2 and F1-versus-F2 weights are poorly identifiable.

### Verdict and required correction

Graphically, four boxes do not justify four synergies. The current E1/E2/F1/F2
cells should be called **four PF phase channels/windows**, not four validated
motor primitives or muscle synergies. Their effective two-dimensional
temporal basis diagnoses redundancy in the present implementation; it does
not establish that the correct biological/controller architecture contains
only two PF layers.

**Correction after Ben's PF-layer clarification:** the earlier
“one-PF-per-synergy” requirement was too literal. Synergy count should be
estimated from multiple aligned cycles/trials using held-out reconstruction
and component stability. PF architecture must be evaluated separately:
joint/functional PF layers may each contain extensor and flexor half-centers,
several PF layers may converge on one MN pool, and their combinations may
generate the observed synergies. The reconstruction target is therefore the
set of synergy and joint-angle patterns—including biarticular contributions
and both knee-flexion episodes—not an identity map from each NMF component to
one PF cell. The later 2026-09-17 section gives the corrected architecture and
Ben's full S1--S6 interpretation.

Reproducible diagnostics added: `_pf_basis_audit.py` and
`_synergy_count_audit.py`.

## 2026-09-17 — Directional commissural correction and FSA backsolve

### V3/C1 topology and the literature diagram

Ben's correction is anatomically and graphically important. The bilateral
pathways must not be drawn as one shared V3 cell or one shared C1 cell. Each
direction uses its own commissural interneuron:

- left RG-E half-center -> excitatory left-to-right V3 -> excitatory right
  RG-E interneuron (`InE_r`), with a separate mirrored right-to-left V3;
- left RG-F half-center -> excitatory left-to-right C1 -> inhibitory right
  RG-F half-center, with a separate mirrored right-to-left C1.

The compiled SNS network already represented these as four source-indexed
cells (`CIN_E_r/l` and `CIN_F_r/l`); the error was in
`draw_literature_circuit.py`, which visually collapsed them into one shared
V3 and one shared C1 population. The figure now draws all four directional
cells and asserts each displayed source, destination, and sign against a
freshly compiled network. `_fix_check.py` passes all eight legs of the four
two-synapse pathways. The representative figure build now sets
`rg_weak_exc=0`, and the audit verifies that no direct ipsilateral
RG-E<->RG-F excitatory connection is present. The optional code parameter is
retained at its default zero for controlled experiments; the change does not
claim that a tuned forward gait has already been re-optimized without it.

### Activation-source provenance correction

`_activation_provenance_audit.py` demonstrates that
`ResultsBSolve/zz_bsolve_StaticOptimization_activation.sto` is **not** the
current converted-MuJoCo ridge/NNLS array in `bsolve_out.npz['acts']`: across
the 121 common frames the RMSE is 0.3924, maximum absolute difference is
0.9900, and correlation is -0.0119. The current NPZ contains no parsed
`so_names`, so the STO cannot be positively reclassified from the surviving
artifacts, but it can be ruled out as the current `acts` array. Consequently,
the preceding four/five-component NMF numbers describe the ambiguous STO
source, not the converted-model activation target requested for the circuit
backsolve.

The FSA analysis therefore uses the unambiguous
`bsolve_out.npz['acts']`, produced by `bsolve_ik.py`'s converted-model,
ridge-regularized bounded least-squares solve along
`subject01_walk1_ik.mot` with measured GRF/CoP from
`subject01_walk1_grf.mot`.

### OpenSim `normal.mot` experiment and why it is quarantined

OpenSim 4.6 was run on `gait2392_simbody.osim` with the Tutorial-1
`normal.mot`. The exact model could not initialize because `lat_gas_r` failed
muscle equilibrium. An analysis-only OpenSim-updated derivative,
`gait2392_simbody_normal_so.osim`, disables only that actuator; the source
model is unchanged. Static optimization remained infeasible until six
explicit pelvis residual actuators were added, after which all 51 frames
converged with constraint violations approximately 1e-12--1e-8.

This is not a physiological walking-activation solution: `normal.mot` has no
measured ground reactions, and the pelvis residuals supply a mean vertical
force of 739.4 N (peak 782.0 N), essentially supporting body weight. Other
residual RMS values are 53.4 N fore-aft, 31.3 N transverse, and
10.3/8.0/40.3 N-m for the three pelvis moments. The resulting activations
(0.010--0.488) are retained in `ResultsNormalSO` as a reproducible
kinematics-only/residual-supported experiment and are **not** used to infer
PF count or FSA synapses. `_normal_so_audit.py` writes the numerical audit.

### Function Subnetwork Approach backsolve

`fsa_backsolve.py` maps each target activation to a nonspiking MN voltage by
`V_MN = 5 mV * activation` and uses the conductance-based leaky-integrator
equation

`C_m dV/dt = G_m(E_r-V) + sum[g_i u_i(E_i-V)] + I_bias`.

Analytical PF-to-MN conductances use the Szczecinski, Hunt, and Quinn (2017)
FSA transmission relation `g = k R G_m / (Delta E_s - k R)`. A separate
dynamic nonnegative-conductance fit tests the same signals in the implemented
leaky membranes. The inverse PF excitation/inhibition traces reconstruct the
six target PF waveforms to numerical precision, but they are requirements on
an upstream RG/PF realization, not proof that the present two-half-center RG
already generates those waveforms.

For the converted-MuJoCo activation target, four components fail the held-out
test: interleaved-frame centered VAF is 0.852 right and 0.823 left. Six is the
smallest **shared bilateral PF count** that clears 90% held-out VAF:

- right: five is the in-sample minimum; at six, NMF centered VAF is 0.926
  and held-out VAF is 0.924;
- left: six is the minimum; at six, NMF centered VAF is 0.913 and held-out
  VAF is 0.912.

Across ten seeds and both alternating phase-frame splits, the six-component
held-out standard deviation is about 0.003 per side. The left spatial
components are fully stable under optimal matching. The right mean matched
cosine is 0.966, but one matched sixth component has a minimum cosine of
0.117. Together with the short 2.0-s record (about 1.6 gait cycles), this
means **six is the justified engineering target for this dataset, not yet a
claim of six universal biological primitives**. More cycles/trials are
needed before hard-coding the architecture.

The dynamic excitatory PF-to-MN FSA fit achieves centered VAF 0.863 right and
0.835 left (activation RMSE 0.125 and 0.133). Approximately 16.1% and 16.8%
of right/left MN samples require net negative current. Thus an
excitatory-only PF projection plus leak cannot reproduce every falling phase;
phase-specific inhibitory pathways (or another explicit negative-current
mechanism) are needed for higher-fidelity closure.

Reproducible outputs are under `Code/MuJoCo_SNS/spinal/fsa_results/`, including
rank selection, six PF time courses and muscle weights, MN reconstructions,
upstream current requirements, JSON/NPZ data, and both PNG and vector PDF
figures. Robustness results are produced by `_fsa_rank_robustness.py`.

## 2026-09-17 — Phase-normalized synergies, joint kinematics, and ankle capacity

### ZCode PF phase channels versus the six extracted synergies

The existing ZCode PF1--PF4 construction and the six NMF components are
different levels of description. E1/E2/F1/F2 were designed as early/late
stance and early/late swing **pattern-formation channels**. S1--S6 are
empirical muscle-recruitment synergies: each has a temporal coefficient and a
spatial muscle-weight vector. A synergy is not itself a PF neuron, half-center,
or layer, and the number of extracted synergies need not equal the number of
PF layers.

A useful PF layer can contain an extensor and a flexor half-center. Those two
outputs may correspond to, or combine into, two observed synergies. Conversely,
one observed synergy may be produced by coordinated outputs from several
joint-specific PF layers. Motoneuron pools can therefore receive convergent
signals from multiple PF layers. This is especially appropriate for
biarticular muscles: their MN pools may combine hip+knee PF drive or
knee+ankle PF drive instead of being assigned exclusively to one joint.
Mathematically the useful model is
`V_MN,m(t) = f(sum_j g[j,m] P_j(t) + sensory + descending bias)`; the NMF
synergies can emerge from correlated PF outputs and the projection matrix
`g`, rather than requiring the identity “one synergy = one PF cell.”

The current `fsa_backsolve.py` lets all six extracted temporal coefficients
project to every MN through fitted conductances. That is a useful unconstrained
reconstruction benchmark, but it does not yet impose the desired anatomical
organization into hip, knee, ankle, extensor, and flexor PF half-centers.
Figures therefore label the factors only S1--S6, not “PF candidates.”

Ben's S1--S6 visual/mechanical interpretation, retained as the working
hypothesis, is:

- **S1:** hip-extensor muscles, including biarticular knee flexors. Because
  these muscles also create abduction, adductor recruitment counteracts that
  action. Trunk muscles provide what is effectively a positive OpenSim-X
  torque to maintain balance while the ipsilateral leg is in stance and
  extending and the contralateral leg is in swing.
- **S2:** slight hip flexion, knee extension, talocrural dorsiflexion, MTP
  extension, and subtalar muscles counteracting inversion/eversion. Ben's
  initial timing hypothesis was the second half of swing; the present
  coefficient's largest peak is instead near 40% on the stance-rescaled axis
  (late stance), so its phase interpretation remains open.
- **S3:** hip extension and flexion, with gluteal recruitment probably
  maintaining sagittal alignment; knee flexion; ankle plantarflexion;
  possibly MTP flexion; and internal trunk rotation.
- **S4:** hip flexion, sagittal-plane hip stabilization, slight knee flexion,
  and internal trunk rotation.
- **S5:** trunk/pelvis flexion, ankle dorsiflexion, and hip abduction.
- **S6:** hip extension and knee flexion.

The stance-rescaled timing also shows why these should not be mistaken for
four chronological ZCode windows. Mean peak phases are approximately S1
95--97%, S2 40%, S3 10%, and S4 81--84%. S5 peaks at 26% right versus 59%
left, while S6 peaks at 4% right versus 23% left. The first four are broadly
bilaterally interpretable; S5/S6 carry residual/asymmetric structure.

`fsa_pf_synergies_r/l` now plots each component on 0--100% gait phase with
heel strike at 0, measured toe-off rescaled to 50, and the next heel strike
at 100. Measured duty was retained in metadata (0.622 right, 0.615 left).
Only one complete stride per side lies inside the activation recording:
right 0.630--1.863 s and left 1.257--2.470 s. Thus this is correctly
cycle-normalized but not yet a multi-cycle statistical mean; no artificial
variability band is shown.

The weight assignment is not final. NMF admits component rotations and can
split one mechanical pattern or merge two correlated patterns. Visual review
should explicitly test whether part of S3 belongs in S6 by plotting their
weighted muscle contributions against joint kinematics and by separating
monoarticular from biarticular muscles. A PF-informed constrained refit can
then compare that interpretation quantitatively instead of manually moving
weights after the fact.

The “double-knee” pattern is a key discriminator. Human walking contains an
early-stance/load-acceptance knee bend and a later swing-phase knee bend.
Plotting both bends together with ankle angle and the hip-, knee-, and
ankle-PF half-center outputs can reveal whether a knee-flexor synergy is being
generated twice by different combinations of joint PF layers. Gastrocnemius,
hamstrings, and other biarticular paths make convergence essential: the same
MN pool can legitimately receive signals from more than one PF layer, and a
single extracted synergy may combine those signals.

The architecture decision is therefore not “four PFs or six PFs.” It is:
(1) how many joint/functional PF layers are needed; (2) which extensor/flexor
half-centers each layer contains; (3) which MN pools receive convergent
outputs from multiple layers; and (4) whether that structured circuit
reconstructs S1--S6 and both knee-flexion episodes without arbitrary
cross-loadings.

### Separate OpenSim joint-angle phase figure

`plot_gait_joint_angles.py` creates `gait_joint_angles_phase.{png,pdf}`.
It uses the same measured-event normalization and preserves OpenSim
coordinate signs rather than remapping them to MuJoCo axes. The left column
contains lumbar extension, hip flexion, knee, talocrural ankle, and MTP.
The right column contains lumbar bending, hip adduction, a deliberately blank
knee panel, subtalar angle, and a deliberately blank MTP panel. Both legs are
shown relative to their own stance/swing cycles. The source IK has no MTP or
subtalar excursion (both are effectively zero); those flat values are data,
not plotting failures.

### OpenSim versus MuJoCo right dorsiflexor capacity

`compare_ankle_df.py` reads Ben's extensionless OpenSim force and torque
tables as whitespace text and evaluates the patched converted MuJoCo model at
the same ankle angles, activation state 1, and zero velocity. Although the
tables say `inDegrees=no`, their -90 to +90 coordinate sweep is necessarily
interpreted as degrees. The measured IK gait range is -8.84 to +16.02 deg.

Within that gait range, the converted MuJoCo dorsiflexors are **not
underpowered at activation 1**:

| muscle | force capacity MJ/OS | torque capacity MJ/OS |
|---|---:|---:|
| ext_dig_r | 1.007 | 1.021 |
| ext_hal_r | 0.966 | 0.974 |
| per_tert_r | 0.945 | 0.973 |
| tib_ant_r | 1.033 | 1.041 |

The summed dorsiflexor torque ratio is 1.023 in the gait range. Across the
entire -90..+90 sweep it is 1.008, with torque-curve correlation 0.999 and
100% sign agreement away from zero. Therefore a dorsiflexion deficit in the
walking simulation should first be sought in recruitment magnitude/timing,
antagonist co-contraction, or force-velocity dynamics rather than Fmax or
static ankle moment-arm capacity.

The comparison uses the same equality-aware central-difference tendon moment
arms as `bsolve_ik.py`. MuJoCo's raw `actuator_moment` is exactly zero for
these converted ankle paths even though their tendon lengths change with
ankle angle, so raw actuator moments cannot be used for this audit.
Reproducible JSON/NPZ data, a markdown table, and PNG/PDF force/torque plots
are under `Code/MuJoCo_SNS/spinal/ankle_df_results/`.

## 2026-09-17 — Bilateral architecture visual-audit stop point

Ben's annotated screenshot identifies four unresolved presentation defects in
panel A of `circuit_literature`. No diagram code or generated figure was
changed in this stop-point pass, so there is no partially completed repair.

### 1. Mirror the right-side RG/PF/motor columns

Both limbs currently use
`ext_x = cx - 0.78; flx_x = cx + 0.78`. This correctly puts the left
RG-F on the medial/inside edge, but puts the right RG-F on the lateral/outside
edge. The right limb should be mirrored so both RG-F half-centers face the
midline. The complete right functional columns—not only the RG circles—must
be mirrored together: RG-F, PF-F, MN-F, and flexor muscle on the inside;
RG-E, PF-E, MN-E, and extensor muscle on the outside. The associated InE/InF
positions and commissural routes must then be rerouted rather than allowed to
cross through the group.

### 2. Pastel-yellow: PF-box connections appear unattached

The intended code paths are RG-E->PF-E, RG-F->PF-F, PF-E/F->MN-E/F, and the
cross-inhibitory PF-IN paths. The generic glyph-to-glyph edge helper does not
give the rectangular PF boxes explicit top/bottom/inner ports. Consequently,
some arrow stems and terminal glyphs stop beside a box or appear to pass
behind it. This is a geometry failure even when the collapsed semantic edge
class passes the compiled-network contract.

Repair plan: give every PF rectangle explicit named anchors (RG input at
top-center, MN output at bottom-center, PF-IN excitation at bottom-inner, and
cross-inhibition at the opposite inner edge). Route each connection between
those anchors and clip it exactly at the rectangle boundary.

### 3. Grass-green: PF interneurons appear self-exciting

The intended edges in code are PF-E->IN-E excitation and PF-F->IN-F
excitation, followed by IN-E -| PF-F and IN-F -| PF-E. The close geometry,
curvature, and excitatory triangle placement make the first two lines look as
though they originate on the interneuron itself. This must be treated as a
failed drawing even though the source/target tuple is correct in the Python
edge registry.

Repair plan: draw the PF-HC-to-PF-IN excitatory stems from the PF box's
bottom-inner port to the interneuron perimeter, place the excitatory triangle
immediately before the interneuron target, and route the inhibitory return
paths on a visibly separate curve. Add instance-level assertions for all four
per-side paths and explicitly forbid an excitatory PF-IN self-edge.

### 4. Indigo: feedback lines appear to originate from nowhere

The dashed aggregate sensory-feedback edges are currently made by direct
`dc.syn(sens, ...)` calls from the shared “Ia / II / Ib + foot contact”
rectangle. They bypass the `neural()` registry used by the compiled-edge
contract, and their overlapping fan-out makes their lower endpoints look
like unattached lines between the muscles. This explains how they survived
the previous contract checks.

Repair plan: connect the sensory rectangle to one explicit, labeled fan-out
bus/port, then branch that bus to PF-E, PF-F, RG-E, and RG-F. Register these
four aggregate visual paths in a separate contract (they intentionally
represent several compiled sensory populations) and ensure the two
muscle-to-sensory encoder paths terminate visibly on the sensory box.

### Safe implementation and verification sequence for the next session

1. Add rectangle-port and feedback-bus primitives without regenerating the
   promoted figure.
2. Mirror the entire right E/F column and reroute its local and commissural
   paths.
3. Replace all PF and aggregate-feedback panel-A edges with named-port routes.
4. Add instance-level contracts containing side, exact source glyph, exact
   target glyph, and sign; retain the existing compiled semantic contract.
5. Add geometry assertions that every path begins on its registered source
   boundary, ends on its target boundary, has nonzero length, and cannot be
   interpreted as a self-edge.
6. Render a temporary review PNG/SVG and inspect all four annotated regions
   before replacing any source or dissertation artifact.
7. Run `_fix_check.py` and the figure contracts, then regenerate PNG/PDF/SVG,
   synchronize the dissertation copies, and verify byte-identical hashes.

Until this sequence is completed, the current `circuit_literature` should
not be described as visually final even though its compiled neural edge-class
contract passes.
