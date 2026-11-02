# Di Russo paper and Gait2392 example: verified comparison

Machine: **easteregg2**. Review date: 2 November 2026. This is an audit of the
existing example and a narrowly defined, opt-in input/event repair. It is not a
claim to have reproduced the published optimization or walking results.

Primary source: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\_dirusso_2023.pdf`.
Page references below give the **PDF page, then the printed journal page**.
Source files are in
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\`.
The three audited implementation files and `dirusso_fast.py` remain unchanged.

## Main finding

The example implements the recognizable architecture of a bilateral phase
oscillator, five signed motor primitives, local reflex interneurons, and grouped
motor outputs. Its executable input path contains three concrete faults:
length and velocity arguments are exchanged; contact tangent force is read as
normal force; and an actual heel-contact onset is replaced by a delayed, held
reset condition. Those faults are separate from the substantial, explicitly
unresolved differences in normalization, balance, delays, parameter bounds,
muscle mechanics, and optimization objective.

`paper_checked_variant.py` corrects only these three input/event faults in the
named `input_event_v1` condition. `legacy_inputs` uses the original input/event
semantics. Both use the same original neural core, gains, signed `W_stage1`,
muscle model, sagittal constraints, initial state, and prescribed support. Both
apply the existing Hunt protective knee bumper, so they are **protected
reconstructions**, rather than exact replays of an unlimited-knee historical run.

## Paper-to-code matrix

| Element | Primary source | Actual implementation and verdict |
|---|---|---|
| Neural state and units | PDF 4/printed 3, Eqs. 1–2: first-order state with 10 ms time constant and weighted inputs | `controller.py:33–34,109–195` and `controller_numba.py:35–79` use 10 ms states and 2 ms Euler updates. Activities and weights are normalized numerical quantities, not SNS conductances in µS or voltages in mV. Applying a conductance-based GLIF formula directly to these weights would require a separate port. |
| Bilateral rhythm | PDF 6/5, Eq. 6 | `controller_numba.py:23–28` updates both phases from the previous phase pair. The division by `2π` correctly converts the paper's angular frequency to the code's cycles-per-second phase derivative. The coupling also contains the corresponding `2π` phase conversion. |
| Contact reset | PDF 6/5, Eq. 7: phase reset at heel strike | The controller's boolean reset is appropriate **when it receives an event**. `dirusso_v3.py:138–169` and `dirusso_fast.py:215–249` calculate a true rising edge, then fail to pass that edge. They suppress reset on the strike and subsequently reset repeatedly while a long contact remains loaded. This is a caller fault, not an oscillator equation fault. |
| Five primitives | PDF 6/5, Eq. 8 and Fig. 3 | `controller.py:62–65` and Numba lines 61–63 use five raised cosines, centres `.1,.3,.5,.7,.9`, half-width `.2`. They additionally wrap phase distance across the cycle boundary. The printed equation/figure show truncated supports; the code therefore adds boundary activity for the first/last primitive. Author implementation confirmation is needed before treating periodic extension as an error. It is retained in the matched repair. |
| Signed output weights | PDF 6/5, Eq. 9, weights constrained to `[-1,1]` | The example has an independent `18×5` array, rather than one shared nine-muscle array. Actual `dirusso_W_stage1.npy` ranges from **−3.1702765530 to +3.4043088587**, with 51 negative entries and maximum left/right difference **1.1026286136**. Thus it is neither paper-bounded nor bilaterally tied. The repair deliberately preserves it, rather than silently clipping or retraining it. |
| Stage-1 generator | Published imitation stage versus local fitting recipe | `dirusso_run.py:67–89` builds reference index 0 from **right** activations and index 1 from **left**, then fits index 0 against **left** heel phases and index 1 against **right**. It uses unconstrained `np.linalg.lstsq`, and each cycle overwrites the preceding cycle's weights. These are concrete generator defects for an ipsilateral, pooled-cycle imitation target. The stored array has no generator manifest proving that this exact revision produced it; that provenance limit is retained. The generator was not rerun. |
| Ia receptor | PDF 7/6, Eq. 11: `.325 sqrt(max(ṽ,0))` | The receptor function is correct: `controller.py:105–107` expects `(length, velocity, force)`. Both physical callers pass `(velocity, length, force)`. The repair changes the argument order. |
| II and Ib receptors | PDF 7/6, Eq. 11: normalized length and force | II consequently receives velocity in the old physical caller. Ib receives the force argument, but group normalization differs from the paper's muscle-normalized signals: the code divides the **mean** component force by the **sum** of component maximum forces. This reduces multi-component group inputs by the number of components. The matched repair retains this separate adaptation. |
| Length and velocity normalization | PDF 7/6, receptor normalization discussion | The example uses the first component's actuator length-range span for the whole group. `V/span` has units `s⁻¹` before numerical clipping; equivalence to the paper's normalized fibre velocity is not established. No optimal-fibre/Vmax domain calibration is demonstrated. |
| Afferent delays | PDF 7/6, Eq. 10 and muscle-proximity discussion | No explicit transport-delay buffer exists in either core. Old sensory/interneuron states introduce discrete pipeline delays, but these are not individually parameterized receptor-to-cord delays. Header descriptions of one-step delay do not establish paper delay equivalence. |
| Ia pathways | PDF 6–9/5–8, reflex text, Figs. 4–5, Table 1 | Own-muscle Ia excitation reaches MN; Ia inhibitory interneurons inhibit listed antagonistic MNs and inhibit one another. RCs inhibit the corresponding Ia interneuron. Paths are ipsilateral. Antagonist lists follow Table 1, which is not perfectly reciprocal: ILPSO lists HAMS, whereas HAMS lists VAS/RF. Do not describe every enumerated pair as symmetric without resolving that convention. |
| II pathways | Same reflex source | Own-muscle II excitation reaches MN through `IIE`; inhibitory `III` pathways reach listed antagonists. There is no explicit locomotor-phase switch for these reflexes, consistent with the paper's description. |
| Ib pathways | Same reflex source | Own-muscle Ib autogenic inhibition, mutual antagonist Ib-interneuron inhibition, and positive-force excitation for the code's extensor classification are present. HAMS/RF/GAS have biarticular roles; the extensor classification is an explicit implementation choice. No generic stance-gating claim should be added to these equations. |
| Renshaw pathways | PDF 6/5 prose, PDF 8/7 Fig. 4e, PDF 9/8 Table 1 caption | Own MN excites RC; RC inhibits its own MN and Ia interneuron. The code uses inhibitory antagonist RC coupling. **The paper is internally inconsistent:** the prose describes reciprocal RC inhibition, while figure/table captions describe reciprocal excitation. The code follows the prose. This conflict is not resolved by silently flipping a sign. |
| Negative internal states | PDF 4/3, Eq. 1 and output clipping description | Numba line 78 lower-clamps RC/IBE/IIE states to zero and safety-clamps other states at ±50; the paper's internal integrator is not generally constrained to nonnegative values. Clipped output is distinct from clipped internal state. This can change inhibitory memory and transients. |
| Cutaneous excitation | Paper Eq. 11 includes a cutaneous signal; contact resets the CPG | `controller.py:175` broadcasts `cut_exc[m]×cut_aff×cut` directly into all nine MN groups. The audited paper's five-reflex discussion and diagrams do not establish that all-MN broadcast as a published rule. It is retained and marked an unresolved additional branch. |
| Balance | PDF 5/4, Eqs. 3–5 | The paper has delayed angle/velocity feedback, including **5 ms hip delay**, with muscle-specific parameters for ILPSO/GMAX/HAMS. The code has immediate finite-difference pitch feedback, default shared `kp=8,kd=1.5,target=0`, an extra `.1` multiplier, and a contact-derived stance threshold. Feedback is added after MN integration, consistent with the placement in Eq. 4, but the delay, sign/frame, gains, and stance rule are not a verified paper implementation. |
| Activation offsets and fitted dimensions | PDF 14/13, `M.NA.w0`; PDF 6/5 parameter discussion | There are no per-neuron activation offsets. Seventeen nine-element gain arrays plus `ff` and `cut_aff` provide 155 numerical gain entries, with entries reused across several antagonist edges. They are not the paper's complete independently fitted parameter set. The paper's stated 48 CPG parameters and the simple 45-weight-plus-frequency count are not reconciled here by inventing missing parameters. |
| Musculoskeletal model | PDF 4/3, model description | Published model: planar nine anatomical DOFs, 18 muscle units, 75.16 kg, 1.8 m. Adaptation: converted Gait2392 with 34 actual actuators, controlled as nine groups per leg; 58 other actuator active/passive force scales are zeroed. MTP remains free, leaving 11 anatomical DOFs after the 12 named sagittal locks. Body mass is retained. This is a grouped adaptation, not the exact published mechanical model. |
| Muscle dynamics | Paper's Hill-type/Thelen model | Default `plant.py:46–89` preserves MuJoCo muscle dynamics, not a verified Thelen compliant-tendon implementation. `linear_muscles=True` is an experimental alternative and is **not used**. The old comments calling linear force the author's assumption are unsupported by the paper. |
| Objective | PDF 10/9, Eqs. 13–14 | `dirusso_v3.py:220–229`: imitation and joint-limit costs are zero; effort is based on the last activation sample, not integrated metabolic energy per mass/distance; the head term is final pitch beyond 15°, not head acceleration; published weights are not reproduced. Its objective must not be labelled an implementation of Eq. 13. |
| Training and support | PDF 9–10/8–9, three optimization stages | The example's 90/70/50/30/0% support ladder is an additional engineering curriculum. It does not reproduce the paper's imitation, locomotion, and reflex co-optimization stages. Prescribed vertical support and exact sagittal equality constraints remain external assistance. |
| Reference initialization | Local adaptation prerequisite | `FastEnv.set_reference_ic()` copies selected anatomical reference coordinates but leaves 82 dependent moving-path/joint coordinates at keyframe values. Initial native polynomial position mismatch reaches **29.85 mm** and velocity mismatch **.547 m/s**. A new, explicit `equality_consistent` initializer analytically projects dependent qpos/qvel while retaining every driver qpos/qvel. Native equality residual then falls to `2.13e-17`. A separate `equality_static` condition additionally declares all initial velocities zero. Neither changes a muscle law, tendon, contact parameter, or the source weights. |

## Which heel policy belongs to which architecture?

For this **abstract phase oscillator**, an actual heel-contact onset may apply the
paper's hard phase reset once. Continuous force can independently supply an Ib or
cutaneous reflex. A force level should not repeatedly hard-reset phase during the
whole stance period.

For a **persistent-Na half-centre**, there is no abstract `φ` state to overwrite.
Contact must instead enter a declared neural pathway, with its sign, duration,
interneuron target, and phase response verified. A brief contact-triggered input
is one possible engineered mechanism; it is not Eq. 7 transplanted literally.
This distinction does not imply that all continuous load feedback is erroneous.

The dissertation's **Syn6** is another architecture: a six-component NMF synergy
basis with RG/PF and reflex additions. It is not the paper's five signed analytic
raised cosines. The K6 basis is owned by ZCode and was not modified in this audit.

## Isolated proofs and runnable comparison

Entry point:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\paper_checked_variant.py`

External record:

`D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\literature\dirusso_input_event_v1\`

The 2,000-step Python/Numba neural parity proof has maximum output difference
`3.3306690738754696e-16`. In a two-bout contact fixture, the old caller resets
349 times; the repaired event caller resets exactly twice, at the two onsets.
For length `.8`, shortening velocity `−.5`, and normalized force `.3`, the old
caller returns Ia `.2906888371` and II `−.5`; the repaired caller returns Ia `0`
and II `.8`. A synthetic contact-frame `[100,−4,−2] N` has 100 N normal force,
whereas the old tangent-index measurement reports zero after rectification.

Warm isolated neural timing was 5.272 s Python versus .006746 s Numba over 2,000
steps, approximately 781.5×. **This is a neural-kernel measurement, not full-body
acceleration.** Both full-body pilots use the same accelerated core; their timing
is reported separately in their result files, including setup and neural warmup.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' -X utf8 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\paper_checked_variant.py' --mode electrical --out 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\literature\dirusso_input_event_v1\electrical'

& 'D:\Anaconda\envs\myo\python.exe' -X utf8 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\paper_checked_variant.py' --mode pilot --variant both --initializer equality_consistent --duration 4 --support 0.7 --out 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\literature\dirusso_input_event_v1\new_consistent_repeat'
```

The first pilot pair is retained. Its independent ground-load logger incorrectly
required the ground geometry to belong to worldbody 0; the actual plane belongs
to the named `ground` body. Therefore its reported zero foot-ground loads are
**invalid telemetry**, not evidence of airborne feet. This measurement defect
does not alter either controller or the body dynamics. A logger-only correction
and matched repeat are required before interpreting load/unloading. Initial
protected knees are −2.136° left and −19.938° right. Both first pilots remain
within the protected range, but flex to approximately 118°: completion is not a
walking pass. No optimizer or training was run.

## Completed matched initial-condition diagnosis

`DIAGNOSIS_CHECKPOINT.json` at the external record root is the compact index.
All six valid-load comparison cases have **byte-identical compiled model
hashes**, the same 70% prescribed upward support, the same W, and the same
protected knee parameters. Source snapshots, including the exact older runner
versions, are in the single external `source` directory. The measurement-only
repeat reproduces the original qpos and MN arrays exactly, maximum difference 0.

| Initial condition | Inputs | COM below 80% of initial height (s) | Mean L/R foot upward load (N) | Peak torso-ground normal load (kN) | Physical verdict |
|---|---|---:|---:|---:|---|
| Original anatomical and inconsistent dependent state | Legacy | .386 | 19.95 / 26.62 | 11.83 | Fall |
| Same | Corrected | .392 | 36.04 / 25.79 | 14.39 | Fall |
| Equality-consistent qpos/qvel, source driver velocities | Legacy | .388 | 30.39 / 19.77 | 11.12 | Fall |
| Same | Corrected | .394 | 39.45 / 19.99 | 14.61 | Fall |
| Equality-consistent qpos, explicitly zero velocity | Legacy | .382 | 21.83 / 22.75 | 11.26 | Fall |
| Same | Corrected | .388 | 7.12 / 16.16 | 13.29 | Fall |

Prescribed upward support is **516.16 N**, against body weight **737.36 N**.
No downward support force is prescribed. Sagittal equality reactions remain
external constraints; their generalized-force channel is saved. The only
additional command beyond the 34 grouped muscle activations is this vertical
support force. The inherited neural balance term and direct cutaneous-to-MN
branch remain present in both cases. There is no online force allocator or new
joint-torque command.

All six cases retain negative knee angles, with maximum extension no closer to
neutral than about **−1.60°**. Flexion reaches about **−118.7°**. Corrected-contact
edge counts during a fall are **not step counts**. A successful four-second
integration therefore does not pass standing, unloading, or walking criteria.

The raw reference starts with knee velocities **+35.3°/s left and +365.5°/s
right**. Removing those velocities does not remove collapse, so neither the
dependent-state mismatch nor initial momentum alone explains the failure.
The repaired input path also does not repair the motor primitive itself:
at .1 s in the original-position comparison, left HAMS/BFSH feedforward drives
are approximately **2.10/1.81**, while VAS drive is negative and its MN output is
zero. This is evidence of an inappropriate stance-recruitment pattern, not proof
that one particular afferent gain causes the entire collapse. The grouped-force
normalization separately attenuates two- and three-component group Ib signals.

Raw joint-coordinate actuator moments must not be mistaken for anatomical
muscle torques: the 105-coordinate converted model contains dependent moving
path points. An anatomical torque audit must include the native equality
Jacobian. Saved direct-coordinate group moments are diagnostic channels, not a
claim that RF has reversed anatomical function.

The legacy integration loop also exposes cached geometry and force fields from
the preceding force evaluation alongside the newly integrated qpos/qvel state.
After the first sample this produces approximately one 2 ms step of alignment
difference. This inherited caller timing is preserved and documented in the
latest configs; it is not a verified receptor-specific physiological delay.

With the full diagnostic logger, the final static-start cases required **4.93 s
and 4.48 s** of integrated-loop wall time for four simulated seconds. Their total
case times, including setup, neural warmup, and compressed output, were **5.83 s
and 4.85 s**. These cases therefore run near real time, despite the much faster
isolated neural kernel. No full-body Python/Numba speed-up ratio was measured.

### Specific next corrective study, not executed

1. Audit the fixed-kinematics expected inputs using a coherent group-force
   definition `sum(abs(F_i))/sum(Fmax_i)`, compared with the inherited
   `mean(abs(F_i))/sum(Fmax_i)`. Monoarticular single-component SOL/TA/BFSH/RF
   stay unchanged; IP/GAS scale by two and GMAX/HAMS/VAS by three before clipping.
   Length/velocity normalization requires its own native-muscle operating-range
   derivation; it must not be presented as paper fibre normalization by naming.
2. Correct the imitation generator's side binding and pool registered complete
   cycles. Only after reference and muscle-force feasibility are established,
   fit bounded `[-1,1]`, functionally mirrored primitive weights against those
   targets. This would be new training, not a relabelled historical checkpoint.
3. Verify the resulting neural waveform against the expected-input fixture
   before another body trial. Do not fit against heel chatter after torso impact,
   and do not use the existing surrogate objective as if it were paper Eq. 13.

At completion of that diagnosis, no normalization change or corrective fitting
had been performed. The separately authorized small fit below supersedes only
the fitting status; the six physical failures and their original weights remain
unchanged.

## Authorized small, bounded corrective fit

Entry point:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\bounded_primitive_fit.py`.
It does not launch physics. The existing source generator and source W were not
edited. Results are in the external `bounded_right_v1` directory.

The actual back-solved activation dataset has **121 frames from .5 to 2.5 s**.
Only one complete right cycle is available: **.6183–1.8533 s, 74 frames**. The
complete left cycle is **1.2467–2.4600 s, 73 frames**. Consequently, no pooled
multi-cycle or independent-right-cycle validation is claimed.

The fit uses **59 right-cycle frames**, reserving every fifth frame (**15**) for
an interpolation check. The entire left cycle is held out. Each reference
activation channel is selected explicitly by its matching `_r` or `_l` suffix;
the old source's right/left indexing reversal is removed. A bounded linear solve
fits 45 canonical-right coefficients in `[-1,1]`; the same nine-muscle matrix is
copied exactly to the left. The existing periodic cosine basis is retained.
Targets are model-back-solved activations, not measured EMG or the author's
published training data.

Actual fitted range is **−.1147836628 to 1**, and left/right coefficient difference
is exactly **0**. The fitted file is `W_mirrored.npy`, SHA-256
`cb4e74b2b07a429b17b5d3765a3112fa8012c9d075066df40363bd7ba73fd9d6`.

The fixed reconstruction checks were per-muscle RMSE ≤ .20 and correlation ≥ .50
when target standard deviation is at least .05. **Both validation sets fail**:

| Check | Main remaining discrepancies |
|---|---|
| Right interpolation holdout | ILPSO RMSE .255; BFSH .233; TA .221; GMAX correlation .422 |
| Entire held-out left cycle | ILPSO RMSE .309; HAMS .233; BFSH .240; RF .219; TA .239; GMAX correlation .306 |

This is nevertheless a concrete reconstruction improvement: descriptive
all-frame, output-clipped activation RMSE falls from **.458 to .177 on right** and
**.472 to .209 on left**. These descriptive comparisons do not replace the failed
preregistered held-out checks or establish physical improvement. The downstream
reflex network and activation dynamics were not fitted by this linear solve.

`--weights ABSOLUTE_PATH` is consumed by the actual `CompiledDiRusso` object. The
isolated 500-step consumption proof verifies exact array identity and a maximum
MN-output difference of **.6033** against otherwise identical zero weights. This
is a neural operator check, not a physical walking result.

A separately named `--phase-init reference_antiphase` condition registers right
phase to the reference start `.5 s` and right heel onset `.6183 s`, using the
measured right period `1.235 s`. The initial phases in cycles are **left
.4042105263, right .9042105263**; left is right plus half a cycle modulo one.
This extrapolates the measured periodic phase backwards to the reference start;
it is declared explicitly, and no independent left phase is fitted. Source
`omega=5.11 rad/s` remains unchanged. `source_default` preserves `[0,.5]`.
Only static preflight and electrical consumption checks have been performed for
the new fitted/phase-registered configuration at that stage. The separately
authorized final diagnostic physical comparison is reported below; it establishes
no body success.

The full restraint ledger is `MECHANICAL_RESTRAINTS.json` at the external root.
The model has 105 generalized coordinates, 82 dependent-joint equalities, and
12 additional sagittal/rigid-trunk equality restrictions: `pelvis_list`,
`pelvis_rotation`, `pelvis_tz`, both `hip_adduction`, both `hip_rotation`, both
`subtalar_angle`, and all three lumbar coordinates. Thus the sagittal pelvis
translations/pitch and hip/knee/ankle motions remain free, with **both MTP joints
also free**, giving 11 anatomical DOFs rather than the paper's nine. The
restrictions are solver equalities, not an online trajectory clamp. The fixed
upward pelvis support, preserved neural balance branch, and identical knee
bumpers are additional declared conditions.

```powershell
# Offline fit only; writes a new bounded matrix and its held-out results.
& 'D:\Anaconda\envs\myo\python.exe' -X utf8 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\bounded_primitive_fit.py' --execute --out 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\literature\dirusso_input_event_v1\new_bounded_fit'

# Electrical consumption proof only; no physical integration.
& 'D:\Anaconda\envs\myo\python.exe' -X utf8 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\paper_checked_variant.py' --mode weights_probe --phase-init reference_antiphase --weights 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\literature\dirusso_input_event_v1\bounded_right_v1\W_mirrored.npy' --out 'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\literature\dirusso_input_event_v1\new_registered_weight_probe'
```

## Final registered-phase physical diagnostic

Two additional four-second cases were expressly authorized as diagnostics after
the failed reconstruction checks. Both use corrected inputs/events,
equality-consistent zero velocity, reference-registered phases, identical 70%
upward support, and identical compiled mechanics. Only old versus fitted W
differs. Initial qpos/qvel arrays and model hashes match exactly.

| Measure | Old source W | Bounded mirrored W |
|---|---:|---:|
| COM below 80% of initial height | .428 s | .476 s |
| Mean left/right foot upward load | 36.23 / 43.55 N | 57.97 / 63.93 N |
| Peak torso-ground normal load | 13.30 kN | 11.99 kN |
| Maximum left/right knee flexion magnitude | 118.31 / 118.19° | 81.51 / 118.03° |
| Physical verdict | Fall | Fall |

Prescribed upward support remains **516.16 N**, with zero prescribed downward
support. No extra applied torque, force allocator, or trajectory clamp is added.
The 48 ms difference in the fall criterion does not establish walking improvement.
The fitted case remains unqualified, and no fitting was performed on these falls.

The compact, visually checked comparison is
`diagnostic_registered_weights\REGISTERED_WEIGHT_DIAGNOSTIC.pdf` (also PNG/SVG).
It shows actual oscillator phase, the five analytic primitives, selected actual
bilateral MN outputs, hip/knee/ankle coordinates, foot-ground loads, and COM
height. It does not invent RG/PF membrane voltages. Grey post-fall shading is
explicitly labelled. All paired axes use common scales.

Two additional physical limitations are now quantified:

- The inherited source hips and ankles are **not limited** in the compiled
  model. Old left ankle reaches −214.88°; the fitted left hip reaches −130.00°.
  Knee protection works in both cases, but it does not make other joints valid.
  `ROM_DIAGNOSTIC.json` gives every actual range, nominal range, and stop flag.
- The shared reference pose starts with **zero ground contacts**. Actual contact
  mesh minima are left toe **81.63 mm**, left calcn **93.89 mm**, right calcn
  **114.58 mm**, and right toe **130.12 mm** above the plane. These are therefore
  partially weight-supported **drop/landing diagnostics**, not initially
  grounded standing. `INITIAL_FOOT_CLEARANCE.json` is a static-FK proof. No
  height correction was silently applied to the completed comparisons.

Consequently, a later grounded comparison would require an explicitly declared
initial-height alignment before treating failure as a pure neural-control result.
No normalization changes were made in these cases. The expressly authorized
grounded follow-up below changes only the initial height and adds rejection
checks, while retaining the failed drop cases.

## Grounded follow-up and closure

The shared `pelvis_ty` coordinate was lowered **82.631719 mm**, placing the left
toe mesh minimum at **−1.000 mm**. Actual initial plane contact was verified;
right calcn clearance was **31.9487 mm**. Only the root coordinate changes.
Prescribed support is present during the initial contact evaluation as well as
throughout the run. Both compared cases have identical model hashes, qpos/qvel,
phases, support, and original muscle/contact/tendon/damping parameters.

A rejection guard checks 18 non-knee rotational coordinates against the actual
input OpenSim model ranges, with .05° numerical tolerance. These stored ranges
are broad (hips ±120°, ankles ±90°), so they are engineering rejection bounds,
not a healthy-walking envelope. No joint limit is added and no state is clipped.
The existing, stricter protected-knee guard remains active.

| Grounded diagnostic | Stop | Joint and angle | Knee maximum L/R | Torso impact before stop |
|---|---:|---|---:|---|
| Old source W | .138 s | Right ankle +90.86° | −2.14 / −16.30° | None recorded |
| Bounded mirrored W | .630 s | Pelvis tilt +90.24° | −1.99 / −11.76° | None recorded |

Both runs fail before a complete cycle; neither crosses the COM-80% threshold
before rejection. Thus COM height alone would miss the early orientation/ROM
failure. The newer weights avoid the immediate ankle excursion but do not
establish balance or alternating ground steps. The requested four seconds are
an upper bound; the guard correctly stops the trajectories early.

Saved-state inspection confirms that the fitted case pushes off and becomes
airborne: at .3 s both foot-ground forces are zero and minimum foot mesh gaps
are **41.24 mm left and 32.82 mm right**. At .1 s, its inherited cutaneous
broadcast contributes **.541 normalized input** to every left MN, while the
effective balance addition subtracts **.259** from each left ILPSO/GMAX/HAMS
output. These are measured pathway contributions, not a completed causal
ablation or proof of a single cause. They identify the unverified broadcast
branch and common-polarity balance defaults as specific mechanisms to test
before another blind strength search.

The final, visually inspected figure is
`grounded_weights_diagnostic\GROUNDED_WEIGHT_DIAGNOSTIC.pdf` (also PNG/SVG).
Hatching starts at actual rejection time; no data are padded after stopping.
`GROUNDED_WEIGHT_COMPARISON.json` records identical physical/initial states and
each stop. `SAVED_STATE_MECHANISM.json` contains the clearance/current audit.
The complete current index is `LITERATURE_CURRENT_CHECKPOINT.json` at the
external record root. No additional fit or body run follows this closure.

MuJoCo force-frame evidence: the official simulation documentation states that
`mj_contactForce` reports three force and three torque components in the contact
frame, whose first axis is normal. A world-frame load requires the contact-frame
rotation and the appropriate force-on-body sign:
https://mujoco.readthedocs.io/en/latest/programming/simulation.html .
