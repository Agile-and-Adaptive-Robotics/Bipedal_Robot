# ZCode task: show the geometry and configurations behind Table 4.12

Create reviewable figures and a per-route configuration/geometry inventory for the **actual saved configurations used in dissertation Table 4.12**. Read existing sources and saved outputs; do not optimize, refit, simulate new outcomes, modify manuscripts, or overwrite results. Check `hostname` first and follow that machine's `AGENTS.md` paths. The paths below use the desktop checkout `D:\Github\Bipedal_Robot`; map the repository root to the machine you actually use.

## Identify the correct table and comparison

Working source:

`D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`

The current 297-page dissertation build is:

`D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\second_pdf_comments_20261030\build_dissertation_r2\main.pdf`

Table **4.12**, label **`tab:CapabilitySurvey`**, is on printed page **115**, PDF page **136**. Its caption is “Mean torque and active-force coverage across 27 actuator sites at 620 kPa.” Confirm against `build_dissertation_r2\chapters\30-torque.aux`. Do not confuse this table with **Figure 4.12**, which shows the cam-extensor route.

The table reports baseline/transmission means of **37.28/65.19% human-torque coverage** and **92.44/98.81% active-force coverage**. These numbers are frozen. The active-force definition is positive insertion tension with slack/infeasibility excluded, not the older reconstructed-strain criterion. Do not recalculate or substitute another metric.

**Critical scope:** this 27-site comparison retains inherited master's-era attachment routing and re-solves actuator/tendon lengths while selecting transmission configurations. It does not, by itself, establish 27 newly optimized origins, insertions, or via-point sets. Extract the actual saved geometry. If an origin/insertion/via point is unchanged, label it inherited/unchanged; do not present initial master's coordinates as newly optimized coordinates. Do not mix the separate knee-route redesign or incomplete Routing27 campaigns into the configurations producing Table 4.12.

## Sources of record

Read these existing files; do not execute the historical optimizer or replay/assembly scripts:

1. **Current per-site comparison and exact saved-file roster:**
   `D:\Github\Bipedal_Robot\tmp\torque_comparison_20261028\goal_table_torque_active_corrected_20261028.csv`
   Its 27 rows give `site`, baseline/updated diameter, `nBPA`, `G`, and `original_mat`/`updated_mat`. All 27 updated MAT paths were verified present locally. This corrected CSV is the metric authority.
2. **Canonical 54-design records and provenance:**
   `D:\Github\Bipedal_Robot\tmp\torque_comparison_20261028\canonical_replay.json`
   `D:\Github\Bipedal_Robot\tmp\torque_comparison_20261028\canonical_replay.csv`
   `D:\Github\Bipedal_Robot\tmp\torque_comparison_20261028\verification.json`
3. **Saved baseline and selected transmission designs:**
   `D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\BiPulley_surveys\`
   Follow the exact `original_mat` and `updated_mat` paths in the corrected CSV. Do not select a different or newer-looking MAT by filename. The selected updated designs use names such as `BiPulley_survey2_ercspn_r_nB1_G1_20261023_051731.mat`.
4. **Original selection roster, for provenance rather than current active-force percentages:**
   `D:\Github\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022\BiPulley_survey2_perMuscle_20261023.csv`
5. **Existing replay/assembly definitions, inspect only:**
   `D:\Github\Bipedal_Robot\tmp\paper_bridge\figscripts\regrade_capability_active_20261028.m`
   `D:\Github\Bipedal_Robot\tmp\paper_bridge\figscripts\assemble_capability_active_20261028.py`
6. **Geometry/kinematics and transmission implementation:**
   `D:\Github\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022\biPulleySpecsFromMasters_20261023.m`
   `D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\biPulleySpecsFromMasters.m`
   `D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\BiPam_pulley2_20261022.m`
7. **Original route ingredients, only when tracing a saved point to its named source:**
   `D:\Github\Bipedal_Robot\tmp\old_code_copy_20261022\AttachPoints_Robot.mat`
   `D:\Github\Bipedal_Robot\tmp\old_code_copy_20261022\PAM_Prop.mat`
   `D:\Github\Bipedal_Robot\tmp\old_code_copy_20261022\MovePoints.mat`
8. **Authoritative context/rulings:**
   `D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\PULLEY_SWEEP_DIRECTIVE_20261021.md`
   `D:\Github\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022\FIX_SPEC_20261023.md`
   `D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\committee_revision_20261027\WHOLELEG_TRANSITION_20261028.md`
   `D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\committee_revision_20261027\PIPELINE_REVISION_20261028.md`
   `D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\94-AppendixC.tex`, particularly `tab:survey_dofs_sizes` and `tab:capability_by_site`.

The old sweep directives describe how the saved studies were run. They are **not permission to run another sweep**. The individual vastus-medialis placement example was waived by Ben; do not reopen it.

## Extract the actual geometry and configuration

The inspected selected MAT contains top-level `spec`, `points`, `T`, `xBest`, `g`, and `surveyMeta`. `spec` includes `points`, `cross`, `frames`, `diameter`, `pulley`, `bundleOffsets`, `kin`, `actuatorId`, `masterName`, and notes. `surveyMeta` includes `nBPA`, `G`, diameter, site, and grading coordinates. The replay accesses the saved object's `Location`, `Cross`, `Diameter`, `TransformationMat`, `RestingL`, `TendonL`, and related fields through `S.g.obj`. Read these existing values; do not reconstruct a “best” design by rerunning its class or solver. If a value is absent or ambiguous, report that field as unavailable with its source issue rather than inventing it.

Produce two linked inventories:

- **One configuration row for each of the 27 sites:** site ID/name, nominal/reference and actual saved BPA diameter in mm, number of parallel BPAs `N_BPA`, tackle line parts `G`, saved resting length and tendon length in m, evaluated joint coordinates, source MAT and field names, and whether the attachment/via geometry differs from the baseline. Include baseline columns or a concise baseline-to-selected comparison so diameter substitution is visible.
- **One ordered point row per physical route point:** site ID, point order/name, role(origin/via/insertion), parent body and coordinate frame, x/y/z in m, fixed versus moving/conditional status, function/coordinate dependence when applicable, and exact source MAT/field or named ingredient. Retain raw precision in CSV. A readable display table may round for presentation with its rounding stated. Preserve duplicate bookkeeping rows in the raw export, but distinguish them from additional physical guides in the drawing.

Check diameter and count separately. The baseline substitutes four parallel 20 mm BPAs at eight nominal 40 mm sites; selected configurations retain their inherited diameter specification. Do not draw every route as 20 mm or assume one 40 mm actuator means four physical 20 mm actuators on the selected side.

Use `G` to mean **supporting tackle line parts/travel multiplier**, not pulley-wheel count. Report `spec.pulley.tackleLineParts` for both modeled crossings where present, and the primary `surveyMeta.G`. A simple redirect can contain a sheave with unity transmission. State physical sheave count only if an actual saved routing/schematic establishes it; otherwise say it is not specified by the abstract transmission record. Do not infer a realizable packaging layout solely from `G`.

## Frame and moving-point constraints

Use the documented anatomical convention **x anterior, y superior, z right**, with knee flexion negative. The saved route points are body-local and must use their recorded parent frames, pivots, joint-axis signs, crossing indices, row folds, and transforms. Never rotate arbitrary CAD placement and call it the anatomical frame. State the displayed pose and all relevant joint coordinates. Prefer stored `T` grid poses; transforming saved geometry for a drawing is allowed, but computing new torque results is not.

Preserve OpenSim moving/conditional insertion and via-point functions, including their knee-angle dependence, when depicting the human reference or any design that actually retains those functions. Ben's origin-only ruling does not permit moving or freezing an insertion to make a figure look convenient.

There is an existing source distinction to disclose: the master's builder used for this survey explicitly substitutes fixed tibia point `t3` for the rectus-femoris/vastus-intermedius patella function points, a rigid-patella approximation. Do not silently replace that frozen survey geometry with a different OpenSim function and still claim it produced Table 4.12. Show the actual saved survey path, clearly identify its approximation, and keep any original OpenSim reference path/function separately labeled. If this prevents the requested figure from honestly being called “updated attachment locations,” say so and title it “Saved transmission configurations on inherited routing.”

If a saved parallel group has only a centerline and configuration count, do not invent individual physical origin offsets, mounts, or spacing. Use saved `bundleOffsets` when supplied; otherwise show a labeled schematic bundle and report that physical placement is unspecified. Avoid unsupported claims of improved packaging or clearance.

## Visual deliverables

Create an actual publication-ready visual set: a readable whole-leg/torso overview with indexed routes, complementary views or route-group closeups exposing origins/insertions/vias, and a compact configuration table. Match figure indices to the 27-site inventory. Avoid 27 tiny illegible panels or a tangle of unlabeled paths. Use registered geometry only; otherwise draw the correct frame/link skeleton rather than guessing bone/CAD registration.

Use clear markers to distinguish origin, insertion, and via points, and labels distinguishing saved selected configuration from inherited geometry. Use Arial at least10 pt at final size, consistent line widths, and CVD-safe colors plus distinct markers/line styles. Put explanatory equations in proposed manuscript/caption text, not pasted into the visual. Do not add decorative or fabricated actuator/linkage details.

Save new outputs in a dedicated review folder, for example:

`D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\table412_geometry_20261030\`

Deliver:

1. Vector PDF visual(s), PNG preview(s), and editable vector/native source plus reproducible render script. Put ingredients in `Components` within that review folder.
2. `route_configurations.csv` and `route_points.csv`, plus a readable PDF/LaTeX table with units and repeated headers if long.
3. A provenance manifest with input hashes, source MAT/field per row, frame/pose definitions, saved approximations, and a 27/27 completeness check.
4. A concise factual report identifying which locations changed and which were inherited, all unavailable fields, and the exact diameter/parallel/tackle combinations actually used.
5. Ready-to-review caption and short LaTeX callout using unique labels such as `fig:wholeleg_saved_routes` and `tab:wholeleg_saved_configurations`. Do not insert them into either working manuscript.

Acceptance checks: all 27 Table 4.12 sites represented; exact selected MAT provenance; no invented configuration or attachment update; point/frame/function integrity; original/selected diameter substitution explicit; tackle line parts distinguished from sheave count; existing four aggregate percentages unchanged; no optimization/refit/simulation; no manuscript, bibliography, protected Chapter 3, or existing figure/result overwritten; no commits or pushes.
