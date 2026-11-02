# FSA and scaling design record

Session 1. Machine: **easteregg2**. Repository: `D:\GitHub\Bipedal_Robot`.

This record is separate from CHATGPT_REPORT.md and the historical training records. No manuscript or Overleaf edits belong to this campaign. ZCode owns the eight-synergy investigation; the existing K6 basis is preserved.

## Required sequence

Ben directed Gait2392 design first, followed by the other active families. Appropriate functional design, signal and physical-time scaling, individual pathways, and control coverage must be checked across **all** families before the next Hunt training or curriculum calculations. Design calculations and isolated electrical tests are permitted. No family receives a training release merely because another family passes.

`CAMPAIGN_GATE.json` remains closed. `design_gate.py` requires reviewed per-family evidence and current hashes for source, applied profile, design, and pathway results. Six prerequisite tests pass: reviewed fixtures release, boolean-only release fails, changed source invalidates evidence, missing control coverage fails, a second pending family prevents release, and unknown operations are rejected. The test fixtures are deleted from their owned temporary directories.

The existing historical runners are preserved. New campaign launchers must call `require_calculation_release()` before launching numerical training. This helper is not a claim that every old standalone script has been retrofitted with a global lock.

## Completed handoffs

- `DEFENSE_DESIGN_DECISIONS.md`: committee questions, defensible design choices, and correction record.
- `AUDIT_ANCHORS.md`: named source/reference and diagnostic conditions, full manifests, and frozen checkpoints. Prepared in response to Ben's ambiguity question; neither modification time nor the latest experimental folder defines an accepted condition.
- `LAPTOP_CAD_SIMSCAPE_REPAIR_HANDOFF.md`: exact machine/file map, mechanical defects, and copyable laptop prompt.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\FUNCTIONAL_SUBNETWORK_OTHER_MACHINE_PROMPT.md`: actual-model graph/pathway audit and editable functional-subnetwork drawing prompt.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\CONNECTOME_TEMPLATE_CATALOG.md`: old templates hidden behind Legacy, preserved data, and current runtime graphs marked import-required. No user localStorage or live drawing was changed.

## Design evidence

### Dimension and torque coverage

`MODEL_DIMENSION_AND_CONTROL.md` explains six pelvis position/orientation coordinates plus 17 internal coordinates (seven per leg and three lumbar). It separates low-dimensional rhythmic coordination from the additional feedback, passive dynamics, ground contact, and apparatus stabilization.

The static `syn6_torque_span.py` diagnostic used four existing walking postures, actual native muscle active-force derivatives, and equality-aware anatomical moment arms. At each posture:

| Activation-to-internal-torque map | Numerical rank |
| --- | ---: |
| Full 92 muscles | 17 |
| 86 muscles represented in the current NMF basis | 17 |
| Fixed bilateral K6 activation basis, 12 primitive inputs | 12 |

Moment-arm finite-difference discrepancy at the checked posture was `1.8189912900576366e-11 m`. Source/model/basis/reference hashes were unchanged. The record is `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\SYN6_TORQUE_SPAN.json`.

These are static activation-map results, not a measured rank of the entire nonlinear SNS or a walking-stability result. Reflex/posture channels are additional inputs. Input bounds, unilateral muscle forces, torque margins, contact modes, and temporal dynamics require their own checks.

### Graded Gait2392 design

Optional profiles and default-off entry points are being checked for s3k, w2lvar, and syn6. The profile targets explicitly declared electrical functions; it preserves NaP rhythm dynamics, PF timing, MN/actuator dynamics, and existing default numerical behavior. A unit-current calibration amplitude is not called a universal physical sensor maximum. Full incoming-pool shunting and final compiled/mirrored/accelerated consumption must be checked separately from isolated endpoint algebra.

Source capacitance numbers use a seconds clock. Numerical C=0.03 with G=1 has a 30 ms leak time constant and corresponds to 30 nF if the physical conductance is 1 microSiemens. It must not become 0.03 nF in the Simulink physical-unit convention.

### Hunt/PT and four-cell runtime correction

`HUNT_GRADED_DESIGN.md` documents a runtime omission: 92 allocated homonymous II-to-MN edges were zero at every tested command because the installed Hunt reflex set omitted IIX. The new opt-in `hunt_graded_profile.py` initializes the four local FSA functions and updates the actual gain buffers used on every command update. Original saved conditions are preserved. Electrical tests on the 1033-neuron Hunt/PT and 1027-neuron four-cell conditions pass: 92 restored II paths, 356 tied local edge pairs, command-policy persistence, exact default no-op, and identical R/L soleus II responses of +0.5 mV. Records are under `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_local_reflex`.

### Spiking functional design

`table2_glif.py` implements the appropriate seven-step pure-GLIF derivation, explicit physical units, and actual reset-and-decay synaptic conductances. NaP/Ca bursting states are explicitly outside that pure-GLIF derivation. They are retained or require their own ionic-spiking design.

The literal published conductance mapping and the selected spiking receiving-cell model produced an approximately 11–12% transmission excess in isolated tests. This negative result is retained; it is not repaired by an undocumented gain multiplier. `conductance_matched_design.py` derives a separate full-conductance firing-rate equation and solves it numerically by monotone scalar bisection. This is a design root solve, not Hunt training or an empirical controller sweep. It explicitly distinguishes rate matching in a spiking receiver from voltage matching in a graded receiver.

The corrected isolated transmission tests reported approximately unit gain across 25%, 50%, and 100% input and matched scalar/accelerated states within floating-point roundoff. Full results and applicability limits are in the agent's compact external design records. These module results do not themselves qualify a complete controller family or release the global gate.

### Non-Gait controller contracts

`NON_GAIT_DESIGN_INDEX.json` and `non_gait_ledgers` cover actual W2L, rat, CAD, SNS-Simscape, and reduced Deng sources. They distinguish verified electrical adapters from uncompleted family design/scaling/pathway work.

Clock factors must be interpreted from their actual implementations. In the current W2L spiking runtime, time and interval are divided by `clock_scale`; a factor 0.2 therefore gives five times the neural rate in physical seconds and one-fifth the neural time constants. It is not a cadence-only factor. Other controllers multiply time constants directly. A canonical profile whose source parameters already incorporate scaling must not be scaled twice. Physiological cadence, sensory encoding, spike carriers, synaptic/readout filters, and actuator dynamics need separately stated physical-time contracts.

## Current release status

The all-family gate is **closed**. Module-level checks are progressing; the full per-family prerequisite set has not yet passed. There have been no new Hunt parameter-training jobs, curriculum trials, whole-body optimization, source-model overwrite, or NMF refitting in this design phase.

### Composition and mechanics review checkpoint

The optional syn6 canonical-right policy now mirrors actual right-side parameters and topology, with initial/dynamic states preserved. It avoids averaging nonhomologous left/right K6 components. The K6 basis file is unchanged. Static parameter equivariance is within floating-point roundoff, and the final dense/accelerated neural trajectories agree exactly. See CANONICAL_RIGHT_DESIGN.md and the external fsa_design\canonical_right_syn6 checkpoint.

CENTRAL_AFFERENT_BUDGET.md describes the existing w2lvar fan-in normalization and the new per-type 0.5 nA current budgets. The source already divided central gains by muscle-group count; edge count alone was not evidence of unnormalized amplification. The new profile covers 520 existing edges in 24 groups and repeats exactly in fresh processes. The saved s3k/syn6 circuits have no qualifying active direct central pathways; their inventory checks are not central transmission passes.

IBEXC_COINCIDENCE_DESIGN.md separates the original additive branches from an optional thresholded phase/load coincidence circuit. All 30 canonical gates pass the specified truth table and release checks. Ib normalization uses the actual canonical stance current of 0.225 nA at F/F0=1. Hunt has a separate wrapper-aware adapter; historical masked outputs remain masked and are not silently reactivated. Combined Hunt command-policy tests are being completed with the local-II and timing modules.

Timing re-encoding now distinguishes function-equivalent normalization, saturation-preserving normalization, and a separately declared canonical-right 0.5 nA-per-RG-pool candidate. The older four-cell source contains substantial left/right timing asymmetry and saturates its original 5 mV source window. Preserving its old function is not a symmetry repair. New candidate timing profiles use the right-side pattern and mirror it; remaining controller buffers and shared midline pathways require their own audit. Exact settings and isolated tests are under D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_timing.

The physical-domain spiking probes preserve a useful negative result: a 16-cell, common-5-nA representation can pass DC and phase checks while giving approximately 45.5% error relative to the small normal-Ib receiver-voltage modulation. That condition is not accepted as a faithful dynamic replacement. A separately named per-afferent-domain population design is being checked against an explicit modulation-relative error criterion; no body trials or learned gain fitting are involved.

The mechanical audit distinguishes 75.1646 kg anatomical mass, 3.2 kg moving auxiliary path-point bodies, and 0.1 kg fixed ground. Both real tibiae are included in anatomical mass. The same saved central-v4 support force is 0.5944% of moving-model BW after excluding the fixed ground, rather than the historical 0.5937%; raw trajectories, files, and acceptance decisions are unchanged. The separate correction is D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\GAIT_V4_SUPPORT_MASS_CORRECTION.json.

DAMPING_REPAIR_USAGE.md documents a default-off repair for inherited joint tags missed by the existing textual patch. It changes 13 compiled internal damping coefficients to their already declared runtime values; other compiled model fields, mass matrices, and stiffness matrices are preserved. All apparatus generalized forces and power are included. Frozen-posture poles and accelerating-state Jacobians are scoped diagnostics, not gait-stability proofs. Evidence is under external fsa_design\damping_repair_v1.

The torque-rank result above uses raw, unpatched CVT3 mechanics. Runtime force and hip-axis patches were not applied; the record must not be presented as the final patched controller's rank.
### Ben's approval and physics clarification, 2026-11-01

Priority is now controllers documented in the submitted dissertation, followed by other active variants. Finish FSA design, scaling, and pathway verification; do not begin Hunt/Szczecinski training or curriculum until Ben explicitly says OK. Completing design checks does not itself authorize training. The release helper now checks this separate user approval.

Compliant tendons, joint damping changes, and muscle damping changes are separate optional experiments. Each requires Ben's separate approval; none is enabled by the neural design work. Preserve existing reference-model mechanics. The worker historically named tendon_damping_physics is currently assigned only Gait neural pathway composition. No new tendon, damping, or body simulations occur in that assignment. Earlier optional mechanics records remain historical evidence. New campaign launchers must declare any requested physics_changes to the release helper.
### Submission mapping and rank-eight information

The submitted-source snapshot is Dissertation\Archive\session1_submission_20261031\source, not the current research-condition list. DISSERTATION_MODEL_MAP.md records this correction. Chapter 6 explicitly describes a 120-neuron public rat precursor; Rat134 is a later curated variant. The original public rat120 and reduced Rybak circuit are now added to the design registry and prioritized with the other submitted controllers. Later research conditions do not replace the submission's results. POSTURAL_COMMAND_DESIGN.md explains the existing POSTURE and PT routes and the incompletely implemented joint-level descending posture organization.

Ben supplied ZCode's Code\MuJoCo_SNS\spinal\fsa_results\rank8\fsa_rank8_report.md. This is a matched inverse-reconstruction comparison: the K6 factors/targets/time/muscle ordering match production synergy_basis.npz, and K8 uses the same targets/time/muscle ordering. Dynamic reconstruction centered VAF improves from 0.8626 to 0.8961 right and 0.8350 to 0.8747 left. These are reconstruction scores, not ground-stepping results. Only one complete cycle per side is available. Required negative MN current remains 16.1%/16.8% of frames, and forward RG-to-PF waveform generation remains unverified. K8 is preserved as a separate analysis; no replacement of the production basis or refit was made.
### Dissertation-first pass/fix review and exact remote drawings

DISSERTATION_FIRST_PASS_FIX.md separates passed numerical/electrical contracts from still-failing or incomplete whole-circuit requirements. Root independently verified all77 hashed files across15 review rows without mismatch. Completed standing-caller supplements now account for all376 s3k and382 W2L/Syn6 ports using source formulas and validated static calibration; W2L contact inputs span0..0.75nA, while selected Syn6 contacts are zero. This does not qualify future changed contact gains or training protocols automatically.

Ben identified the actual route-builder ledgers on EB475WS4. Root checked remote hostname and all3SHA256 values through SSH; every ledger matches easteregg2 byte-exact. ROUTE_BUILDER_EXPORT_CHECK.md records names, counts, hashes, settings, and actual edge counts. They are standing exports, not deafferented air-run copies. Both s3k/Syn6 exports have0direct Ia/II/Ib-to-central edges; W2L has152/152/216. Source s3k supports conditional central routes with gains zero; Syn6 lacks general Ia/II central implementation and has an optional disabled Ib block.

Syn6 hardwires ipsilateral heel excitation of both RG interneurons; their outputs inhibit opposite RG half-centers. The chosen standing caller silences those contact inputs, so it did not test the intended heel-reset function. Earlier operator/unit/cache tests were insufficient to establish correct end-to-end neural routing. A separate functional topology audit and explicitly named default-off corrections are now required. Training, all body integration, compliant-tendon changes, and damping changes remain prohibited until Ben's OK. The existing physics reference is unchanged.
### Historical-trial correction and defense-deck handoff

DISSERTATION_TRIAL_PROVENANCE.md corrects the earlier scope: later standing exports do not establish the topology used for dissertation full-s3k trials. The exact pruning results match Appendix D's rounded scores; full-s3k loads trial34 through stage5, with nonzero centralIa/II/Ib, full_rules1, IaINenabled,and Renshaw0.5. Source-only reconstruction with current code gives878neurons; this is explicitly not historical byte-identical replay. HistoricalSyn6 has zero centralIa/II/Ib but nonzeroheel/toe gains. The pass/fix report now leads with this distinction.

Requested defense-deck summary: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\SESSION1_NEURAL_CONTROLLER_HANDOFF.md. All30 absolute paths in the document were checked and exist. It links the actual trial provenance, design reviews, media manifests, measured support, committee answers, and records. It documents zero-conductance semantics and the current sparse backend's construction-time nonzero filtering/cachedconductance copy; stale cache is a possible explanation, not a proven diagnosis of Ben's prior observation.

The exact fourth EB475WS4 W2L108 export was copied read-only to an owned audit directory. All108neuron fields,216connections,and47synapse types match the canonical source; completeXML differs only inProjectPath. No remote source or drawing was changed. No training, body run, tendon change, or damping change occurred for this handoff.