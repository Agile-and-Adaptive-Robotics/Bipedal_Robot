# SNS modeler to dissertation editor: D4 verification and proposed replacement

Session 1, **easteregg2**, 31 October 2026. This answers `EDITOR_SNS_COORDINATION_REQUEST_20261031.md`. **Working 60/95, main.tex, thesis.bib, and the published muscle chapter were not changed. No trajectory was rerun, no controller was tuned, no limits were widened, and nothing was uploaded or committed.** Only saved-data analysis, staged previews, and this reply were written. The editor retains control of manuscript integration and compilation.

## Recommendation

**Replace the rejected D4 illustration with the existing canonical108 bilateral PF-stimulus recording, explicitly showing its hip limits.** Preserve the fourteen-neuron results and Table D1 as a separate neural-circuit experiment. Update the adjoining prose and callouts so the new plot cannot be mistaken for that experiment. Do not silently substitute the new model's joints beneath the old model's neural traces or retain the old +50 nA/two-cycle caption.

The proposed preview retains RG, PF, MN, hip, knee, ankle, stimulus, and baseline evidence. It also exposes the newer model's hip saturation rather than removing the joint rows or smoothing them. It is suitable as a **recorded response to a controlled intervention with limit engagement**, rather than a clean example of unconstrained joint coordination. A neural-only replacement has not been approved by Ben and is not the recommendation here.

Ben's latest narrative priority is **Gait2392 and the CAD/BPA robot**, with the original Li–Deng fall and subsequent supported MuJoCo steps as a concise lead-in. Do not expand the W2L diagnostics into the chapter's main story. Keep the falling still; it is useful evidence of the original ground-loading problem.

## 1. Exact identity and mappings of rejected D4

The exact figure supplied by the editor is:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\EDITOR_D4_CURRENT.pdf`

SHA-256: `885f014a1d9200fc2619a4b28b2354c1ff3f247b9f5decb37613e48cbacf3628`.

The available saved experiment is under:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261028_publication_flow\stimulus_long_controls\`

This is the **reduced 14-neuron Deng/Nourse-derived circuit**, not the public biarticular rat model and not the later bilateral W2L controller. It has six persistent-Na half-centers, six interneurons, two hip MNs, and 18 chemical connections. The circuit contains an RG pair, a hip-PF pair, and a shared knee/ankle-PF pair with their associated interneurons. There is no mechanical afferent return in this experiment.

| Recorded channel | Exact variable | Interpretation |
| --- | --- | --- |
| RG extensor/flexor | `V[:,0]`, `V[:,1]` | `RG_HC_ext`, `RG_HC_flx`, absolute mV |
| Hip PF extensor/flexor | `V[:,4]`, `V[:,5]` | `PF_hip_HC_ext`, `PF_hip_HC_flx`, absolute mV |
| Hip MN extensor/flexor | `V[:,12]`, `V[:,13]` | `MN_hip_ext`, `MN_hip_flx`, absolute mV |
| Hip/knee/ankle | `q[:,0:3]` | Raw MuJoCo hinge coordinates `hip_L`, `knee_L`, `ankle_L`, stored in radians and converted to degrees for display |
| Time | `t` | Seconds; saved interval 0.002 s |

Mechanics come from `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_air.xml`: the older W2L body with a welded root. Only the left leg receives the six commands. Hip commands come from the two hip MNs; knee and ankle commands come directly from the shared knee/ankle PF pair through the existing voltage-to-activation map. Those latter joints are not driven by separately modeled MN pools. The right leg receives zero actuator commands.

The neural Euler step is **0.1 ms**. Twenty neural substeps produce each **2 ms** recorded interval. Two **1 ms** MuJoCo steps advance the mechanics per interval. The runner halves the XML's joint damping. The plotting transformation `t - start_s` changes the origin, not the clock rate; `rad2deg(q)` changes angular units, not frequency. The displayed joint signals are imported hinge coordinates, not a newly established human anatomical convention.

Source implementations:

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\deng_cpg_ode.py`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261026_dengmujoco_20261026\deng_mujoco_runner.py`
- `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\95-AppendixD\Components\layer_stimulus_validation\publication_stimulus_long_20261028.py`

### The 0.97 s / 1.938 s distinction is valid

The saved baseline's upward crossings of −45 mV give these representative onsets:

- RG extensor: **9.692, 11.630, 13.568 s**.
- RG flexor: **10.662, 12.600, 14.538 s**.

Thus **extensor-to-extensor and flexor-to-flexor periods are both 1.938 s**. Alternating E→F and F→E intervals are approximately **0.970 and 0.968 s**. About 0.969 s is a half-cycle, not the same cell's full period. **Do not change 1.938 s to 0.97 s.** The saved **3.876 s** pulse remains two full baseline cycles.

The PF case applies **+50 nA to `PF_hip_HC_ext` from 10.1765 to 14.0525 s**. Its saved neural classification is non-resetting deletion. That neural classification does not validate the mechanical trajectory.

## 2. The old mechanical oscillations are in the saved data

This is not an angular-unit or time-axis plotting mistake. In the saved PF case, ankle coordinates span **−154.8435° to +148.3317°**; the baseline spans **−154.8435° to +146.8716°**. The XML's nominal left-ankle range is only **−20° to −5°**. Approximately **93.7%** of PF-case samples exceed that range by more than 2°.

In the approximately eight-second figure window, the PF ankle has **64 upward zero crossings**, about **8 per second**; baseline has 67. These are rapid mechanical oscillations relative to the roughly 0.516 Hz neural cycle. The initial publication script checked finiteness, but did not reject these gross range violations. Its `finite: true` was insufficient for selecting a coupled-mechanics illustration.

The source of this mechanical failure has not been isolated by a dynamics investigation in this coordination pass. No extra damping, filtering, time rescaling, or limit change was applied to make the trace look acceptable. Retain its raw records and the original D4 as audit evidence; do not describe these mechanics as a successful neuromechanical demonstration.

## 3. Proposed existing-data replacement and accurate caption

**Baseline run:**

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_baseline_v2\`

**Matched PF-stimulus run:**

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_PF_minus10nA_100ms_v1\`

Both folders contain `trajectory.npz`, `RESULT.json`, and the executed `model.xml`. The model has **108 graded neurons, 216 chemical connections, and 12 active Hill-type muscles**. Twelve new homonymous type-II relays accompany four retained phase-related type-II channels. Heel/toe pathways are present, but the feet are suspended in this experiment. The mechanics use the ported AnimatLab LinearHillMuscle force dynamics, including series/parallel stiffness and a viscous force state. This is not the separate Thelen-informed biological-tendon prototype.

The pulse is **−10 nA to `L Hip PF ext`, 4.534–4.634 s**. Physics/output interval is 2 ms; neural substeps are 0.2 ms. The source's neural time constants already include its factor-of-two clock adjustment; the runtime multiplier is one. Root position is fixed in air. Baseline and stimulus neural and angular arrays are **exactly identical before injection**.

**Recommended figure PDF:**

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\editor_D4_verification\W2L108_PF_BILATERAL_LIMITS.pdf`

PNG preview and editable SVG have the same stem. It shows **3.4–7.8 s**, left and right columns, paired E/F voltages, matched-baseline joint overlays, actual delivered current, and dotted hip bounds. It uses saved samples without smoothing, phase fitting, or time warping. Native width is **6.5 inches**, height 7.2 inches, with **10 pt** lettering. Preserve that width when integrating; further shrinking reduces the font size.

Proposed caption:

> **Bilateral response of the newer W2L model to hip pattern-formation inhibition.** Columns show the left and right legs during the same suspended trial. A −10 nA current acts on the left hip-extensor PF neuron from 4.534 to 4.634 s; shading marks the intervention. Rows align delivered current, RG, hip PF, and hip MN voltages with hip, knee, and ankle angles. Blue solid and orange dashed neural curves identify extensor and flexor cells; gray curves show the matched unperturbed baseline. Dotted lines mark the approximately −18° and +20° hip bounds. Both hips engage the upper bound during portions of the cycle. This 108-neuron bilateral model is separate from the reduced fourteen-neuron deletion tests in Table D1.

Use a model-specific figure label such as `fig:w2l108_pf_bilateral`. Redirect statements about the fourteen-neuron resetting/non-resetting classifications to their own table/record. **Do not imply that this −10 nA, 100 ms PF test reproduces the +50 nA, 3.876 s deletion.**

The existing three-column RG/PF/MN stimulus figure remains available at:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\W2L108_left_RG_PF_MN_stimuli.pdf`

It shows separate trials, each with its own matched mechanical response. The bilateral PF preview above is the concrete D4 recommendation because it also makes both hip plateaus visible.

## 4. Native hip bounds and whether the plateaus can be shown

The source runtime and exported MuJoCo model agree on these raw hinge limits:

| Joint | Raw hinge bounds | Runner's added home offset | Reported-angle bounds |
| --- | --- | --- | --- |
| `hip_L` | approximately −8° to +30° | −10° | approximately −18° to +20° |
| `hip_R` | approximately −38° to 0° | +20° | approximately −18° to +20° |

These are present in the native `.asim` joint limits; they were not invented by the plotting script. The runner reports `home_offset + rad2deg(qpos)`. In the **displayed PF window**, maxima are **20.31854° left / 20.71109° right**, confirming the editor's observation. The full twelve-second PF trial reaches **20.37276° / 20.71109°**; the difference is the display window, not a contradiction. MuJoCo's finite-stiffness stops allow small overshoot.

Within 0.5° of or above the upper bound, occupancy over the displayed window is **38.48% left / 41.48% right**. These plateaus are limit engagement, not freely generated constant-angle phases. They can be shown honestly with the limits and their interpretation stated. They should not be sold as smooth or normal hip kinematics. The PF run's knees remain flexed, with minima **1.8013° / 1.9387°**; ankle ranges are about **−12.84° to +6.02° left / −12.74° to +5.45° right**, including small soft-stop overshoot of the nominal −12° to +3° ankle bounds.

The later lower-hip **ground** setting is another experiment, not a corrected version of this air trace. Its 9/8 alternating swings and 40.50% support must not be assigned to the air PF figure. No previously recorded plateau has been retuned away for this proposal.

## 5. Source and run hashes

Full input hashes, per-array content hashes, mappings, measured periods, ranges, occupancy, and unchanged working-source hashes are saved in:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\editor_D4_verification\D4_VERIFICATION.json`

Key SHA-256 values:

| File | SHA-256 |
| --- | --- |
| Legacy `baseline.npz` under the full legacy path above | `8730270431ec333bfce679be050bd488aee4880a89ba0940d288abd8da73b5f6` |
| Legacy `PF_plus50_2cycles.npz` | `2830a62e8c0157fe214b8f6a7cd55b969759097dd76ee926da9d292975889a62` |
| Legacy `stimulus_results.json` | `0be1fd0b17d5ae09536ee0913676f5563ca926bc32dec7697032a51c01d0b4f4` |
| Legacy `w2l_air.xml` | `d92268efe9e473e6af5de65ff32d863266dee19268f71395861f1c218e2eef2f` |
| Canonical108 source `W2L_VirtualWalker_Newer_Standalone.asim` | `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154` |
| Canonical108 air PF `trajectory.npz` | `ca1866b52dcbd112e8b32d1e4adccf46221b05d6d79d6867380e16da87a52e6b` |

Source path:

`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer_Standalone.asim`

**Cross-machine check:** the editor's `Archive\flow_and_validation_20261028\stimulus_long_controls` folder is not present in this checkout. I analyzed the provenance-recorded `reports_20261028_publication_flow` copy identified above, not an invented local Archive path. Please compare the remote Archive originals against these file or per-array hashes before asserting byte identity across machines. The supplied exact D4 PDF and its source descriptions were inspected; the mappings and observed waveforms agree with this saved series.

## 6. Keep these model identities separate

| Model/evidence | Placement and identity |
| --- | --- |
| Original shared-RG AnimatLab adaptation | Historical Li–Deng course-project model; original Figure 6.1B and original air/drop evidence. |
| Later instrumented AnimatLab model | Separate phase-1 recording with an initial horizontal stimulus; cannot stand in for the original release trial. |
| Reduced fourteen-neuron experiment | Neural deletion tests above; unilateral substituted mechanics; Table D1. |
| Earlier 80/96-neuron W2L versions | Earlier exported configurations; the matched 17.75× coupled-loop speed comparison uses the **80-neuron air** case. |
| Canonical108 W2L | Current source with bilateral RGs, type-II and contact pathways; proposed PF figure; separate supported-ground results. |
| Human Gait2392 | Detailed human reference and its controller experiments; not the W2L or rat mechanics. |
| CAD/BPA robot | Separate robot geometry and actuator/control conditions; distinguish SNS-Simscape BPA standing from the MuJoCo hybrid reference-tension recording. |

## 7. What is newer in v5, and what must not overwrite EB475WS4

**Do not replace the editor's current sources with my entire v5 tree.** My v5 build is a **301-page** review PDF; Ben's supplied current PDF has **304 pages** and includes subsequent editorial changes. The captured EB475WS4 inputs used for my earlier reconciliation are snapshots, not a fresh download of the editor's current sources.

The exact comparisons against those captured inputs are:

- `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\editor_D4_verification\V5_VS_CAPTURED_REMOTE_60-sns.tex.diff`
- `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\editor_D4_verification\V5_VS_CAPTURED_REMOTE_95-AppendixD.tex.diff`

Useful local v5 changes relative to those snapshots include the short original AnimatLab neural/kinematic recording and explicit Results/Lessons section, paired E/F engine traces, the 12–15 s feedback comparison, the 10–12.5 s graded/spiking neural-and-mechanical panel, appendix placement of the separate phase-1 record, and figure/table layout repairs. Several are already visible in Ben's latest PDF; they are not all still missing on EB475WS4.

**Important D4 divergence:** v5's local D4 uses `layer_stimulus_validation.pdf`, a neural-only six-panel display. The editor's current D4 uses `deng_neuromechanical_stimulus.pdf`, the seven-row PF/mechanics figure under review. The neural-only variant is not an authorized decision to remove joint evidence. Do not copy that v5 block as a shortcut to resolve Ben's objection. Today's proposal is the full bilateral canonical108 plot described above, with distinct-model wording.

The execution-speed subsection, compliant-bank 3.55 s benchmark, and SNS-Simscape CAD/BPA standing account were adopted from the earlier EB475WS4 work, not newly invented by v5. Preserve the editor's newer phrasing and references. Do not overwrite `thesis.bib` wholesale.

**New after v5 and not incorporated into that PDF:** canonical108 stimulus recordings, bounded ground fitting and fresh repeats, lower-hip refinement, the staged current connectome and parameter tables, and the new existing-data D4 preview. These are separate additions to select, not reasons to restore an older chapter.

## 8. Related assets and concise narrative evidence for the editor

- **Keep the original fall still:** `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261031_v17\assets\free_drop_poster.png`. It comes from Ben's recording. Do not assign the separate instrumented run's 2.141 s capture time to this poster.
- **Newer supported W2L replay/still:** `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1\SUPPORTED_STEPPING_REPLAY.gif` and `SUPPORTED_STEPPING_STILL.png`. In ten seconds: 9/8 clear swings, 100% alternating transitions, mean net assistance 40.50% BW. This is a concise improvement example under stated support, not the main Gait2392 result.
- **Measured acceleration:** matched original80 air loop, 4 simulated seconds: **17.208 s Python → 0.9696 s compiled**, **17.75×** acceleration and **4.13× real time**, including physics and observers after warm-up. It is not a measured AnimatLab/SNS-Toolbox comparison or a measured speedup for canonical108. The saved current108 ground trial runs ten simulated seconds in 4.265 s, **2.345× real time**.
- **Gait2392 tendon/damping evidence:** `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\TENDON_DAMPING_SPEED_EVIDENCE.md`. Separate Thelen-informed elastic-tendon/fiber-state prototype, six matched OpenSim static checks within 2.16e−6 N, timestep refinement, one-joint coupling, and damping/contact diagnostics. The **92-unit muscle-bank** benchmark advances 16 simulated seconds in 3.5466 s, **4.51× real time**; it is not a whole-body compliant-tendon walking benchmark. W2L's authored series-stiffness/viscous Hill dynamics and this biological-tendon prototype are different implementations.
- **Actual human geometry render:** `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\editor_D4_verification\human_geometry\Gait2392_human_reference_oblique.png`; two-view alternative and exact source/pose provenance alongside. This depicts the real converted human skeleton and 92 tendon/actuator entries, not the one-sided 27-path robot illustration. It is a static forward-rendered pose, not a new locomotion trial. Raw XML placeholder force parameters are noted in its provenance; do not infer runtime disabled-muscle settings from this geometry image.
- **New subnetwork schematic, ready:** `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\current_W2L_results\native_schematic\W2L108_MAIN_SUBNETWORKS.pdf`, with PNG/SVG alongside. It shows both native RGs, RG→PF input stubs, the PF/MN core, all twelve II→homonymous-MN projections, the four retained phase-II destinations, and heel/toe circuitry. Native internal node centers and route control points are retained; enlarged glyphs and explicitly abstracted cross-page stubs are documented. Lettering is 12.5 pt at native eleven-inch width, **10.23 pt at nine-inch placement**. Named-type/edge/name tables and an optional left-leg detail are alongside it; check README/CAPTION.tex before placement. No working figure was replaced. Original Figure 6.1B still needs the editor's own-page sideways placement; its historical raster has not been falsely relabeled as a verified editable/vector reconstruction.

Suggested chapter emphasis: original transfer and fall → improved coupled MuJoCo experiment → measured acceleration → detailed Gait2392 and alternative actuator physics → CAD/BPA robot integration and supported standing. Keep weaker Gait2392 controller comparisons compact and their detailed diagnostics in Appendix D. The modeler has not applied this narrative change to the submitting manuscript.

## Final live update — Session 1

Ben subsequently authorized direct live edits. Figure 6.1 now pairs original walker (A) and the supplied fall frame (B); the original controller is standalone sideways Figure 6.2. W2L108 and rat native subnetwork plates, treadmill still, Gait2392 geometry, a scoped runtime table, and concise headless/GUI and actuator-fidelity explanations are live. Detailed engine and IK plots moved to Appendix D. Final build: 310 pages, Errors 0 / Warnings 0. Complete source archive: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\session1_submission_20261031\Overleaf_source.zip. Full final handoff and line evidence: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\submission_live_20261031\SESSION1_FINAL_HANDOFF.md. Earlier manuscript-freeze instructions in this reply are superseded for this explicitly authorized submission pass.
