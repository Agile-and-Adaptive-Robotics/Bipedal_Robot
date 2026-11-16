# Independent air and W2L108 evidence audit

Canonical Chapter 6 and Appendix D were read before analysis. This audit reads saved recordings and models only. No simulation, training, controller edit, or source-recording change was performed. Reproduce with `D:\Anaconda\envs\myo\python.exe AIR_audit.py`; results, exact file hashes, arrays, and runtime are in `AIR_RESULTS.json`. The successful read-only calculation took approximately 1.6 s.

## Air-stage identity and command contrast

Evidence root: `D:\Users\Ben Bolen\ControllerTraining_20261104\hunt_pose_air_validation_v3`. The two stage directories contain `body.json`, `SOURCE.json`, `body.npz`, `body.assistance.npz`, and `capture\{aligned_trace.npz,ALIGNED_TRACE_METADATA.json,model.mjb}`. The stage-2 compiled-model SHA256 is `a9a82bdfa201510d9976f9f1cc4a55ec701d1756707d92fe376cadce591ea1fc`. The aligned recording has 7,700 samples at 2 ms; accepted state/readout times are 0.602–16.000 s, input times 0.600–15.998 s. The first 300 held-pose preparation intervals are excluded, rather than invented as motion.

`ordered_output_fits_v5\hunt\C\MANIFEST.json` binds the same per-muscle fit to normalized command_F=command_E=0.15. The retained, later body-probe source and runtime declare the following stage definitions. Numeric F/E input samples were not separately logged in the adopted air captures, so this is a reconstruction from manifests and compatible adapter code rather than independent recovery of consumed command samples:

| Condition | Stage 1 | Stage 2 |
|---|---|---|
| Sustained normalized F/E commands, both sides | 0 | 0.15 |
| Sensory return | Disabled; sensory ports zeroed and central afferent gains zeroed | Live Ia/II/Ib and declared central sensory routes enabled |
| Contact-event return | Disabled | Disabled |
| Common onset seed | +0.4 to E_r and F_l for the first 50 neural calls (100 ms) | Same |
| Direct posture, balance, and function channels | Clamped off | Clamped off |
| Legacy shared DRIVE/HUNT_INIT/ALPHA/GAMMA ports | Disabled by declared runtime | Same |
| PF-to-MN expression | 0.1 | 0.1 |

Thus unchanged fitted weights do not isolate sensory feedback as the cause of the 2.466→1.300 s period change. The retained stage definitions change both command and sensory conditions. **Source provenance limit:** the exact original `body_probe.py` hash in SOURCE.json (`d73dee...`) was not found among the four later snapshots; the reviewed `source_snapshot_final` version hashes `acc6c287...`. A subsequent read-only check of all three Git-history blobs also found no match; AIR_source_check.py and AIR_SOURCE_CHECK.json preserve this bounded search. The numerical evaluator `existing_dependencies\hunt_training__physical_eval.py` does match the exact captured hash `ecda3d...`, but delegates neural construction to the unavailable exact-hash launch adapter. The saved body configurations do not themselves specify numeric F/E command values. Consequently the exact consumed commands are not independently proven by that evaluator/configuration pair alone. Manuscript wording should attribute the 0/0.15 contrast to the retained adapter and calibration manifest, and explicitly preserve this source-binding limitation, rather than imply complete launch-source recovery or measured input-command parity.

Principal-cycle scoring uses `5 <= t_state < 14 s`, with knee signs inverted solely for flexion-peak detection. Retain 4R/4L and 7R/7L, periods 2.466/1.300 s, and antiphase errors 0.0336667/0.0932768 cycles from COMPLETE.json. The source explicitly preserves the failed legacy stage-1 assessment (secondary peaks counted as cycles). Raw knee maxima are independently recomputed: stage 1 R −1.299960°, L −1.325969°; stage 2 R −1.717975°, L −1.800829°. Every recorded foot-force component is exactly zero in both aligned captures.

**Cycle reexecution:** the myo environment initially failed because SciPy's initialization calls `os.chdir('D:\Anaconda\envs\myo\lib\site-packages\scipy\.libs')`, raising `PermissionError: [WinError 5] Access is denied`. This bounded failure remains recorded in AIR_RESULTS.json. The archived detector was then successfully executed with available `D:\Anaconda\python.exe` (NumPy 1.26.4, SciPy 1.17.1), using AIR_cycles.py. Its antiphase/double-hump, synchronous, frozen, and fast controls pass. Both raw captures reproduce the adopted counts, period, phase error, and period CV exactly: scalar differences are 0.0. AIR_CYCLE_RESULTS.json preserves the successful independent calculation and environment. No dynamics was integrated in either calculation.

## Body-weight denominator and clock distinction

The compiled model confirms a 0.100000 kg fixed ground body. The legacy evaluator sums every model body's mass, including that fixed body. The aligned recorder uses the pelvis subtree. Both include the model's moving path-point bodies; neither denominator is a newly estimated biological body mass.

| Statistic | Mass / denominator | Window and observed value |
|---|---|---|
| Adopted stage-2 mean suspension support | Total model mass 78.464600 kg; 769.737726 N | Legacy pre-interval `t>=5`, samples 5.000–15.998 s; mean 768.661982 N; **99.860245%** |
| Same force divided by pelvis weight (diagnostic only) | Pelvis subtree 78.364600 kg; 768.756726 N | Same numerator/window; **99.987676%** |
| Aligned superior harness-force mean | Same pelvis denominator | Accepted post-interval `t_state>=5`, 5.000–16.000 s; 768.660991 N; **99.987547%** |

The small force difference between rows 1 and 3 is due to time masks/sample phases, in addition to the distinct denominator. Do not silently replace 99.860% with the metadata-normalized value or normalize plotted aligned forces using the legacy statistic. Plotting new air assistance in newtons is preferable. The existing phrase "from 5 s onward" is correct for the adopted support statistic; `5<=t<13` in the evaluator applies separate foot/horizontal/lateral diagnostics, not mean suspension support.

## Phase-specific current pulses

Evidence: `hunt_layer_stimulus_2nA_v2\Baseline`, its RG/PF/MN/PF_inhibition conditions, and `hunt_active_phase_stimulus_v3\{RG,PF,MN,PF_inhibition}`. Actual logged pulse arrays independently confirm 50 samples at 2 ms per condition. Positive pulses are +2 nA; the inhibitory PF control is −2 nA. These are external bias-current additions, not voltage clamps.

| Condition | Right-side target | Current interval | Maximum-difference scoring window |
|---|---|---|---|
| Earlier, PF silent | RG_F_r, PF_HIP-F_r, or MN_iliacus_r, separately | [7.100,7.200) s | [7.100,10.500) s |
| Later, PF active | Same three targets, separately | [7.650,7.750) s | [7.650,11.050) s |

The maxima use the neural readout clock; the interval marks use `t_neural_input`. Reported differences are between separately stimulated trajectories and the no-stimulus baseline. Full pre-onset saved neural traces agree exactly for the three compared channels.

The outgoing PF ramp begins at 0 mV and saturates at 5 mV (rest-relative voltage): `clip(V/5,0,1)`. During the silent-PF excitatory pulse the maximum PF voltage is −3.057632 mV. The pulse therefore changes PF voltage (maximum difference 0.550443 mV) while transmitted activity stays zero and iliacus MN difference is only 5.196e−8 mV. At the active phase, the PF pulse produces maximum PF/MN differences 1.141985/0.161378 mV. The direct MN pulse produces 1.790992 mV MN difference at 7.750 s, followed by the 6.038807 mV RG-F difference at 9.146 s. These are response-window maxima, not direct synaptic gains or instantaneous transmission ratios. The same active-phase RG pulse gives only 0.333859 mV RG difference and approximately 1.34e−6 mV iliacus MN difference; its effect cannot be generalized from the silent-phase RG condition, whose RG difference reaches 10.095661 mV.

The existing two Appendix figure strips come from `Defense\Session1_Handoff_20261104\experiments_20261104\BiW2LvarX_active_phase_stimulus.pdf`, as explicitly bound by `Notes\ch6_revision_20261116\split_biw2l_stimulus.py`. PDF text confirms the 7.5–11.0 s ticks and right hip-flexor PF/iliacus readouts. These are the **7.65 s active-phase** experiment, not the 7.1 s silent-phase controls. Captions should state this distinction; the silent-phase control is additional textual evidence.

Recommended interpretation change: bind the chapter's small PF response explicitly to the right hip-flexor PF at 7.65 s, distinguish the silent 7.1 s phase, and state the 3.4 s response window in Appendix D. The existing 1.14, 0.16, 1.79, and 6.04 mV rounded values are retained, rather than altered.

## Fixed-weight RG-current noise

`hunt_rg_noise_v1\COMPLETE.json` retains six independent noisy cases (seeds 41001–41003 at stationary SD 0.01 and 0.025 nA), one zero control, and one repeated noisy seed. All eight stored assessments pass; six are noisy cases rather than eight independent noisy trials. The selected repeated seed's raw `qpos`, `neural_voltage`, and actual `stimulus_nA` arrays were independently compared in this audit and all maximum absolute differences are exactly zero. The 20 ms correlation time and four RG inputs remain correctly scoped to suspended motion; no ground-balance or parameter-basin claim follows.

## W2L108 side binding, assistance, and independent replay

Evidence root: `D:\Users\Ben Bolen\NeuroCampaign_20261031\LiDeng_W2L_20261104`. The selected recording is `body_trials\W2L_VirtualWalker_Newer\rustuna108_trial_0010_1793854859534189500\trajectory.npz`; independent replay is `fresh_repeat\W2L_VirtualWalker_Newer\repeat\trajectory.npz`. All **32** numeric arrays were compared again here and agree exactly. This is reanalysis of saved repeat files, not a newly simulated repeat. The fresh model XML hashes `296422d21116b22b88cd167c9ad2b0e13512aa2fda8c45ac2d9b3dc013dd5e74`.

The recording orders joints **hip_L,knee_L,ankle_L,hip_R,knee_R,ankle_R**, so minima are **left 1.946153° and right 2.086698° anatomical flexion**, over the post-startup t>=1 s assessment. Explicit side labels are needed wherever the pair is quoted. These positive anatomical flexion values are distinct from the converted human model's native negative-flexion coordinates.

The mass is 41.919002 kg and BW=411.22540962 N. Upward assistance independently recomputes to **29.5187688%** over t>=1 s (1.000–9.998 s); separately clamp rail, guide, and harness superior forces to nonnegative before summing. Whole-trial signed net assistance is **35.0685746%** over 0–9.998 s, with no positive clipping. Guide upward and downward forces can cancel only in the latter statistic. Continuous harness force is zero, but rail and guide forces remain.

Retain eight left/seven right qualifying swings and 14 alternating transitions. The retained gate requires whole-foot clearance >1 mm for at least 80 ms with peak >=10 mm, ipsilateral load <0.025BW, and contralateral load >0.05BW for at least half the interval; at least three intervals per side and >=80% alternation. This is engineering supported-stepping acceptance, not physiological validation or unsupported walking. The source contains 24 actuator entries for twelve antagonist muscle routes; do not infer the actuation count solely from the phrase twelve-muscle body.

The later `W2L_short_slow_shuffle_20261104\EXPERIMENT_RECORD.md` and associated FRESH_VERIFICATION retain a separate 16 s candidate, seven swings per side, and 10.105% upward assistance over 2–16 s, but failed short excursion and sustained double-support criteria. Preserve this failure; do not promote it as a more favorable replacement of the earlier accepted supported-stepping condition.

## Numerical/interpretation ledger for proposed manuscript changes

| Old | New / retained | Reason |
|---|---|---|
| Air stages' commands unspecified; differences described generally | Explicit stage 1 F/E=0 versus stage 2 F/E=0.15 plus stated sensory contrast | Same weights do not isolate feedback causality |
| 99.860% with unspecified BW denominator | Retain 99.860%; disclose 78.4646 kg/769.737726 N legacy denominator, fixed-ground inclusion, and 5–16 s window | Avoid silent renormalization when aligned metadata uses 78.3646 kg |
| Pulse target/phase/onset/window partly unstated | Right-side named cells; 7.1 versus 7.65 s onset; 100 ms pulse; 3.4 s response window; 0/5 mV PF ramp | Maxima and silence depend on operating phase and threshold |
| W2L108 1.946/2.087° with no side binding | Left 1.946°, right 2.087° | Raw joint ordering and per-leg assessment are explicit |
| W2L108 assistance windows described generally | Retain 29.519% upward after1 s and35.069% whole-trial signed net; define clipping and mass | Different metric and window, not contradictory estimates |

No model points, transforms, stiffnesses, bounds, flags, actuator laws, or numerical recordings were changed.
