# v20 directive and change record

Finalized 2 November 2026 from the authoring sources, final package/mapping, advisor annotations, four independent perspective reports, combined static review, final validation reports, and completed bounded media checks. The final PPTX/PDF is promoted as the current v20 delivery. Cleanup is complete; spoken timing remains unrehearsed.

**Status definitions:** **Static review complete** means the feature is covered by the combined full-deck review and subsequent changed-slide checks described below; the final rerender matches those images exactly. **Package verified** means the final-file validator checked the stated property. Separate media reports establish the stated playback results. These checks do not establish spoken timing or constitute a new numerical recomputation. All positions use the final 127-slide order; authoring source IDs remain stable. Source line references identify the implementation snapshot reviewed.

## Main narrative and presentation requirements

| Directive | Implementation and evidence | Status / final check |
| --- | --- | --- |
| Establish the case before the research goals, with about ten setup slides | Final 2–10 cover terrain, actuator rationale, mechanics, feedback, anatomy, range/torque mismatch, reciprocal research, and prosthetics; roadmap at 11. `pass2_content.mjs:14–40`; opening evidence map. | Static review complete; final 3 and 9 rechecked after pass 3 |
| Make biology ↔ robot research explicit | Final 9 has biological evidence → robot principles and robot experiments → biological predictions, with a scoped hypothetical intervention in notes. `pass2_content.mjs:30–34`; `pass3_content.mjs:12–13`. | Static review complete; final 9 rechecked after pass 3 |
| Explain why BPAs are used | Final 3 states lightweight, compliant actuation and uses an actual BPA and Puppy robot; first use expands the acronym. `pass2_content.mjs:14–16`; `pass3_content.mjs:11`. | Static review complete; final 3 rechecked after pass 3 |
| Link limited contraction and force to range of motion and torque | Final 7 shows the historical range comparison; 8 historical torque mismatch; 4 explains the moment arm. Notes identify model evidence and interference as another range constraint. `pass2_content.mjs:17–29`. | Static review complete; source/model conditions retained |
| Use an explicit feedback loop without a generic negative summing sign | Final 5 is an editable controller/body/sensory loop, with descending drive and returned length/force/contact. `revise.mjs:42–49`. | Static review complete |
| Use relevant authentic imagery, including a prosthesis example | Original anatomical/robot assets replace the rejected opening pairing; MIT photograph credited in final 10. `pass2_content.mjs:20–37`. | Static review complete; source references retained |
| Keep SADb out of the main talk | Three-gap roadmap and three-contribution conclusion; latest conclusion notes do not mention SADb. `pass2_content.mjs:38–40`; `pass3_content.mjs:41–55`. SADb source 117–123 maps to final 121–127. | Package verified: no SADb in main text/notes; 75–127 hidden |
| Keep approximately 70 main slides and ≤45-minute talk | Source 14 and 50 omitted; order is 74 main + 53 questions. `revise.mjs:146`; `package_current.py:14`; `build/slide_order.json`. | Package verified: 74 main; 37-minute baseline, 39-minute target, hard 45-minute limit; unrehearsed |
| Emphasize published force and experimentally grounded torque work | Final 12–53 is the force/torque sequence, 42 of 74 main slides including dividers. Neural section begins at 54. | Package verified; static review complete |
| Preserve the supplied C3NS/RoseGreen design | Packager retains v19 masters/layouts; supervisor saw small branding and honeycomb in the prior opening render. `package_current.py:90–94`. | Package verified: original title, masters, layouts, themes, and template media unchanged; static review complete |
| Use chapter-title dividers | BPA, knee torque, and neuromechanical model/tool dividers remain at final 12, 27, and 54. | Static review complete |

## Force and knee mechanics

| Directive | Implementation and evidence | Status / final check |
| --- | --- | --- |
| Existing model deficiencies before the proposed model | Manufacturer discrepancy at final 15; existing alternatives grouped as Others at 16, without the proposed trace. `revise.mjs:69–70`. | Static review complete |
| Explain the experiment before its results | Test apparatus and current length at final 17; maximum-force protocol at 18 holds `lm=l0`, pressurizes to 620 kPa, and measures F620. `pass2_content.mjs:47–49`. | Static review complete |
| Show maximum-force reasoning before diameter results | Protocol 18 precedes 10 mm result 19, 20 mm result 20, and pressure extension 21. `revise.mjs:72–81`; order JSON. | Static review complete; final order verified |
| Define clamp-to-clamp initial/current/free-contraction length | Final 13 retains source dimensions and explicitly labels free contraction at 620 kPa. `pass2_content.mjs:42`. | Static review complete |
| Preserve individual free-contraction measurement requirement | Corrected notes on final 22 reject a reliable resting-length prediction from the scattered fit. `pass2_content.mjs:106`. | Notes retained in final package; static review complete |
| Do not put force-data sample counts on main slides | 537-point/80%-20% record is in notes for source 24/final 23 only, with the non-held-out-actuator limitation. `pass3_content.mjs:14–15`. | Main count omitted; source note retained; final 23 rechecked after pass 3 |
| Close the 20 mm, 450 mm deficiency with corresponding evidence | Final 26 uses the published 20 mm 300/450 mm pair. `pass2_content.mjs:44–45`. | Static review complete; published figure retained |
| Use authentic moment-arm and bracket-frame methods | Final 29 retains both source torque constructions; 35 bracket CAD; 37 two-rotation figure. `revise.mjs:82–83,114–117,128–132`. | Static review complete |
| Correct signed χ0 meaning without invented photo dimensions | Final 34/36 use the signed no-tendon equation and describe the unloaded length discrepancy; tendon subtraction is in notes. `revise.mjs:84–94`; `pass2_content.mjs:50`. | Static review complete; no original result recomputation |
| Explain the inner δ solve separately from parameter fitting | Final 38 uses the original projected-displacement detail and reference-scaled force balance; δ and dimensionless response definitions added. `pass2_content.mjs:71–75`; `pass3_content.mjs:16–19`; `make_equations.py:10`. | Static review complete; final 38 rechecked by parent and supervisor |
| Identify χ before showing corrected curves | Final 39 has flexor/extensor fit and validation stages, fixed shared χ1/χ2, residual objective, and selection caveat. `pass2_content.mjs:76–80`; `pass3_content.mjs:20–21`. | Static review complete; final 39 rechecked by parent and supervisor |
| Retain all three requested goodness-of-fit measures | Final 41 is an editable source-value table with RMSE, FVU, maximum error, separate training/validation rows, and means. `revise.mjs:95–113`. | Static review complete; source-value table retained; not a new arithmetic recomputation |
| Match apparatus photographs to the actual experiment | Requested pinned-extensor photo on source 40/final 42; original biomimetic apparatus retained on source 42/final 44. `pass2_content.mjs:51–55`; `pass3_content.mjs:28–31`. | Static review complete; final 42 and 44 repairs rechecked |
| Separate correction-factor optimization from route redesign | Source 126/final 45 legend says rigid-body/corrected; same repair added to source 42/final 44. `pass2_content.mjs:63–65`; `pass3_content.mjs:33`. | Static review complete; final 44 rechecked by parent and supervisor |
| Preserve physical shortfall → design change → GIF → result causality | Final 45–51: single-BPA shortfall, corrected route design, flexor GIF/results, direct extensor GIF and hardware overload, cam mechanism/results. `pass2_content.mjs:57–70`; order JSON. | Final sequence and static size verified; native PowerPoint-export movement verified for route slides 46, 47, and 49 |
| Use actual route GIFs at useful size | Existing full-cycle GIFs are embedded with native crop; aspect ratio is preserved. Notes identify illustrations and the flexor provenance limit. `pass2_content.mjs:57–62`; `revise.mjs:17–23`. | Static size/aspect ratio reviewed; movement verified at six times across a source cycle on 46, 47, and 49; see bounded GIF report |
| Show direction as well as torque magnitude in the design requirement | Final 46 says required direction and “across the joint range”; 1.05 magnitude target retained. `pass3_content.mjs:22`; `make_equations.py:11`. | Static review complete; final 46 rechecked by parent and supervisor |
| Explain cam contacts and geometric constraints | Final 50 adds BPA-path, contact-point, and clearance/bend-constraint labels. `pass3_content.mjs:25–27`. | Static review complete; final 50 rechecked after pass 3 |
| Scope whole-leg comparison and avoid rejected terminology | Final 52 says 27-site model, 620 kPa, actuator sizing/pulley ratios, and angular-domain coverage. `pass2_content.mjs:103–105`. | Static review complete; final rejected-wording scan clear |
| Keep physical claims at measured angles | Flexor and cam-extensor plots retain the selected targets and measured-angle claims; no dynamic performance claim added. Independent source reviews and unchanged source figures govern. | Source review and static check complete; native chart/workbook bytes unchanged |

## Neural models and question slides

| Directive | Implementation and evidence | Status / final check |
| --- | --- | --- |
| Follow Table 6.1's model/tool flow | Final 55–68 proceeds through walker/controller, Gait2392, blocks, CAD, bridge, training, intervention, and acceleration. Current mapping; submitted `60-sns.tex:84–91`. | Static review complete; final order verified |
| Explain newer bilateral rhythm/contact feedback without clipped fragments | Final 57 uses actual partial subnetwork crops with explicit retained type-II feedback; pass 3 restores the top source heading. `pass2_content.mjs:82–86`; `pass3_content.mjs:35–37`. | Static review complete; final 57 rechecked by parent and Reddy-perspective reviewer |
| Do not conflate W2L108 and recovered-W2L results | Final 59 is the newer supported trial; final 60 explicitly names recovered W2L, fixed controller/drive/support, analysis window, and alternation reference. `pass2_content.mjs:87–92`. | Static review complete |
| State what the contact-phase statistic summarizes | Final 60 identifies 4.5–30 s median, with 12–15 s plotted excerpt retained. Same source as above. | Static review complete |
| Identify reduced-circuit perturbation and exact protocol | Final 67 names reduced Deng circuit, +50 nA, and 3.876 s; notes say fourteen neurons and neural-response test. `pass2_content.mjs:108–109`. | Static review complete |
| Correct duty-side labels without changing results | Final 81 caption becomes 0.71 right/0.65 left/reference 0.61, checked against the existing JSON by neural reviewer. `pass2_content.mjs:93–94`. | Static review complete; label correction only |
| Explain own-peak normalization | Final 82 note/label says each curve is divided by its own peak; label moved below the plots. `pass2_content.mjs:95–97`; `pass3_content.mjs:38`. | Static review complete; final 82 repair rechecked by parent, supervisor, and Reddy-perspective reviewer |
| Distinguish imported rat diagram from reduced-circuit response | Final 107 labels the model organization and separate reduced-circuit test explicitly. `pass2_content.mjs:100–102`. | Static review complete |
| Preserve SADb question content and authentic app views | `sadb_questions_fix.mjs` edits sources 117,118,120,121 with authentic cropped screenshots and separate literature finding/proposed test; 119 and 122 retained. | Static review complete; all seven questions hidden at final 121–127 |
| Preserve native SADb offline walkthrough | Source 123/final 127 is deliberately not marked changed or transplanted from the imported artifact representation. `sadb_questions_fix.mjs` final comment block; `package_current.py:38–55`. | Package and playback verified: native mapping/bytes preserved; click-to-play clip advances with its original no-loop setting |
| Preserve original charts, embedded workbooks, movies, and template | Packager reuses native parts and checks equality to v19. `package_current.py:90–94`. | Package verified: 7 chart and 7 workbook parts unchanged; all original MP4 assets retained |
| Freeze manuscript, model, and numerical result sources | Presentation source and documentation only. No manuscript, protected Chapter 20, database, or simulation result edits requested. | Submitted PDF hash unchanged; no scientific source edits from this documentation task |

## Recorded changes and exclusions

- The current source order adds six mechanics/route slides (124–129) and omits two redundant source positions (14 and 50), yielding 127 total.
- Numerical findings are retained. Corrections concern scope, statistical interpretation, and labels; the duty-side swap is explicitly recorded above.
- Existing actual images/GIFs are reused. Crops and small outside-geometry masks change the presentation window, not underlying experimental data. `build/gif_mask_check.json` records the narrow 199-frame flexor mask check, not complete animation verification.
- All 127 final-output images exactly match the inspected pass-3 images. Historical reports retain their snapshot scope; the full-deck coverage and changed-slide rechecks are recorded below.
- The parent promoted the current-delivery pointer to v20 after media checks and completed the documented cleanup. No commit, push, or Overleaf upload was performed by this handoff task.

## Final verification evidence and remaining work

The final package, static review, and bounded media checks are complete. The parent promoted v20 as current and completed cleanup. Spoken timing remains unrehearsed.

| Evidence | Verified result and boundary |
| --- | --- |
| `build/delivery_checks.json` | Final PPTX **90,314,963 bytes**, SHA-256 **`69ea30c4fe7ab533f1c8a74b258cf431ad3025c4c4b8a946215916393aa9d272`**; all package checks pass, with no listed failures. |
| Final order and PDF | **127 total, 74 main, 53 hidden**; hidden positions 75–127 exactly; PDF 127 pages. |
| Final text and notes | No main-talk SADb or rejected wording found; no missing slide relationships. |
| Native quantitative material | All **7 charts and 7 embedded workbooks** are byte-identical to the v19 baseline; owners 15, 25, 30 (two), 52, 101, and 102. This is conservation, not a new experimental-value audit. |
| Original appearance and movies | Title, masters, layouts, themes, and template media unchanged. Four original MP4 assets retained, with three active clips on slides 56 (two) and 127 (one); one was already inactive. Seven GIFs retain multi-frame content. |
| `build/validation_delivery.json` | Package integrity, font policy, native chart/workbook structure, and first-party import pass. Zero blocking findings; three nonblocking layout heuristics concern slides 36, 41, and 60, all covered by rendered inspection. Table arithmetic skipped cases without unambiguous totals. |
| Pass-2 full-deck static review | Final 1–12 supervisor; 13–53 Turcic-perspective reviewer; 54–105 Reddy-perspective reviewer; 106–127 parent editor. |
| Pass-3 changed-image inspection | Parent re-opened all 15 changed positions: **3, 9, 23, 38, 39, 42, 44, 46, 47, 50, 57, 74, 82, 105, 127**. Supervisor independently rechecked **38, 39, 42, 44, 46, 74, 82**; Reddy-perspective reviewer independently rechecked **57, 74, 82, 105**. |
| `build/final_render_comparison.json` | All **127 final-output images** exactly match the inspected pass-3 images; zero differences. Static coverage is combined, not an assertion that each reviewer inspected every slide. |
| Submitted dissertation | PDF unchanged, SHA-256 **`447c2cc5a9dd219fa0857e298727fb68d4797af8f406251865eb22f4b126361f`**. No manuscript/results edits were part of this documentation task. |
| `build/native_gif_playback.json` | Native PowerPoint video export shows movement at six times across a source GIF cycle for final slides **46, 47, and 49**. The parent checked two exported frames per slide for clarity. This bounded result does not certify all seven GIFs. |
| `build/native_video_playback.json` | Both clips on **56 advance and restart their loops**. The offline walkthrough on **127 advances with its original click-to-play, no-loop setting**. All three active native MP4 clips have completed playback checks. |
| Current pointer and cleanup | `Defense/current_deck.json` points to the verified v20 PPTX/PDF. `build/cleanup.json` records deletion of the v18 final PPTX, v20 draft/candidate packages, and stale slide-128 render, totaling 307,614,121 bytes. Final PPTX revisions **19 and 20** and the original v18 source assets are retained. |
| Remaining pacing check | **37-minute baseline, 39-minute target, hard 45-minute limit; unrehearsed**. The target is not a measured duration. |

The review reports are independent analytical committee perspectives, not statements by the actual committee. The supervisor's and reviewers' iteration-specific findings remain in their original reports, including corrections that prompted pass 3. No blanket independent-review approval is inferred from those historical reports.

`HANDOFF.md` contains the explicit ZCode source/builder/assets tree. Original v18 submitted source, figure assets, and `transplant_v18.py` remain retained and required after deletion of its standalone final PPTX. The v19 baseline PPTX/assets and v20 runtime source dependencies are retained. The parent completed promotion and cleanup; this documentation task changed only these two records. Only spoken timing remains unmeasured. A future new physical photograph would be an optional user-directed revision, not an outstanding delivery requirement.
