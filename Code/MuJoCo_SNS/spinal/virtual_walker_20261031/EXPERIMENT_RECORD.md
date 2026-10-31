# W2L virtual-walker experiment record

## Latest actual108 checkpoint and stimulus results

The current canonical NEWER W2L virtual walker has **108 neurons/216 chemical edges**: the 96-neuron bilateral-RG/contact/phase-II circuit plus **12 independent homonymous excitatory II relays**, covering flexors and extensors at hip, knee, and ankle on both sides. All four older phase-II pathways remain separate. The twelve new physical inputs are matched to actual actuator route geometry and motor-neuron bindings; all twelve carry measured nonzero II→MN currents (means approximately 0.25–0.32 nA, maxima 0.34–0.42 nA) at 0.005 µS. This is W2L, **not the separate 120-neuron rat controller**.

Both actual runtime filenames now have SHA-256 `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`. The first attempt labelled `canonical108_air_baseline_v1` accidentally read a stale96 standalone because an explicit exporter output list omitted it; that result is marked invalid. The corrected `canonical108_air_baseline_v2` and all following stimulus trials enforce the actual108/216 census. The GUI is styled separately without changing these numerical runtime files. Actual source membrane/gating time constants already include clock2; runtime multiplier is **1**, not2.

The actual108 afferented fixed-root air baseline runs for12 s, with period approximately1.072 s and bilateral knee excursions over22°. Three stable cycles complete at4.510 s: RG-period CV0.088%, knee-amplitude CV2.6%/4.9%. The predeclared shared active-phase stimulus onset is4.534 s. Each of three matched runs applies **−10 nA for100 ms** to respectively L RG ext, L Hip PF ext, or L Hip MN ext. Actual substep-current integration is −1.000 nA·s, and pre-onset neural/physical trajectories are bit-identical in every case. All complete without numerical failure or knee hyperextension.

| Stimulated layer | Measured later L RG-ext peak shift | Main observed response |
|---|---:|---|
| RG | Approximately +14 ms | RG interruption/rebound and downstream changes; maximum left hip change22.2° |
| PF | Approximately −2 ms | Strong local PF/MN interruption; maximum left hip change35.8° while RG timing remains close to baseline |
| MN | 0 ms at2 ms resolution | Local attenuation; RG/PF peak times unchanged, maximum left hip/knee changes1.15°/2.45° |

These are single matched perturbations, not statistical robustness claims. A missing local voltage maximum is not automatically a deleted whole burst. The early RG rebound must not be treated as a full-cycle phase reset by matching it to the wrong later baseline peak. Three already-started optional positive-current probes also completed (PF+2 nA short/long and bilateral RG+2 nA); their data remain separate, with no added gain sweep.

Actual data and short figures: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_air_baseline_v2\` and sibling tags `canonical108_air_RG_minus10nA_100ms_v1`, `canonical108_air_PF_minus10nA_100ms_v1`, and `canonical108_air_MN_minus10nA_100ms_v1`. Shared summary: `runs\CANONICAL108_STIMULUS_TRIPLET.json`; actual twelve-channel check: baseline `HOMONYMOUS_II_CURRENT_VERIFICATION.json`. These new108 tests are air stepping; prior96 supported-ground percentages below do **not** constitute a108 ground pass. No nondefault segment-mass experiment was run. Current work now transfers the distinct genuine120-neuron rat controller to only approved named Gait2392 muscle channels.

**Selected paper model: canonical NEWER only**, `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`, with its synchronized standalone in the same folder. Its intended architecture is 96 neurons with bilateral RG modules, four separate heel/toe contact channels, and four actual type-II physical adapters: hip and knee channels on both sides. There are no ankle or extensor-labelled II channels. Older80 results below are historical comparisons, not the selected paper model.

The bounded NEWER96/204 diagnostic previously passed air and supported-ground tests with tied 0/0 RG-extensor tonic currents and V3 outputs zero (6/8 clear swings, 84.6% alternation, 35.2% BW assistance). A single **positive V3 scale 0.01** candidate retains those pathways and also passes: deafferented air 7/7 full knee cycles, phase 0.487; afferented air 6/7, phase 0.509; ground 8/8 clear swings with 100% alternation, **73.53% BW assistance**. Effective positive V3 output is 0.02749 µS, not zero. Data: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\POSITIVE_V3_0P01_V1.json`.

These results use the preserved GUI-authoritative pre-template-numeric-update snapshot, now an immutable historical input under `runs\source_inputs\`, not an active family copy. The canonical NEWER is currently being updated to Ben's exact reference names and numerical synapse fields; **the prior passing numbers cannot be assigned to that changed canonical export without a fresh check**. No completed five-stage/Hunt curriculum, independent balance, or overground walking is claimed. Six redundant generated family folders have been archived for cleanup; actual inputs, hashes, records, and raw runs remain preserved.

Detailed canonical pathway audit: `runs\canonical_newer_audit\AUDIT.json`. The snapshot audited four real II adapters in GUI and runtime, rather than inferring them from names. Effective nonspiking gains come from the referenced synapse type, not the separate GUI link-local conductance. Both hip II-labelled receptor routes have extensor-sign geometry at the assembled pose; this is recorded as a source mapping issue/possible heteronymous feedback, not silently rewired.

Machine: **easteregg2**. Session 1, 2026-10-31. This record is separate from the ZCode handoff and the Gait2392 training record.

## Scope and controller identity

Both requested W2L standalone models have runnable direct MuJoCo imports. The selected original uses its actual **80 nonspiking neurons and 150 chemical connections**, including a **single left-labeled RG half-center pair with its two RG interneurons**, bilateral hip/knee PF layers, twelve motor neurons, Renshaw neurons, Ia/Ib pathways, and shared knee–ankle coordination. It is not the generic 95-neuron W2L network, the Gait2392 biw2lvar-Hunt-PT controller, or a biarticular model. Twelve physical actuators and twelve source stretch-receptor elements retain their native sigmoid and Hill dynamics. No allocator supplies reference tensions.

The source-faithful equations are implemented in sparse, preallocated Numba loops. Neuron voltage/capacitance/current units, exponential gating, receptor outputs, saturation, and internal Hill tension dynamics were checked against public AnimatLab source. A shared time scale of 2 changes membrane and intrinsic gating time constants equally; it leaves neural/physics timesteps and physical muscle parameters unchanged.

Source original: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\W2L_VirtualWalker_Standalone.asim`, SHA-256 `5beb078060b5fd3193f315924b32ef6ce5860fc4678687848729e3d80949d535`.

## Structural checks before calibration

- Body positions and routed lengths reproduce the explicitly reflected native initial frame to numerical precision. The anatomical target is x anterior, y superior, z right. Gravity is −9.81 m/s² along Y; source gravity was −9.8.
- Moving-pose native marker checks gave maximum RMS errors of 0.659 mm for the original and 0.818 mm for the newer model; incorrect hinge signs produced 88–115 mm errors. The residuals are consistent with native finite joint compliance, not exact cross-engine dynamics parity.
- A separate unchanged-source native 2 s belt probe measured world-X displacement −0.12743 m. After 0.6 s, native belt velocity averaged −0.08963 m/s while its prismatic velocity averaged +0.08954 m/s under command +0.1 m/s. World displacement equals negative joint displacement within 0.105 mm, confirming the port's motor/slider sign independently of the hinge convention. Native motor tracking differs from the port's finite velocity servo. Native versus ported root progression also differs, so matching initial poses and route geometry is not a claim of complete cross-engine motion parity. Evidence: `Neuromechanical_Models\W2L_VirtualWalker_20261031\build\native_belt_sign\BELT_SIGN.json`.
- All twelve muscle paths are uniarticular. Right-hip physical names are reversed, but source neural output adapters already compensate. No right-hip signal swap was made. Ankle naming alone was not treated as proof of wrong phase.
- A labeled bilateral ankle PF phase-swap variant and a shared knee-flexor gain increase were tested and rejected. Broad parameter sweeps were not used.
- The source support consists of passive root-mounted wings on horizontal rails, not a hidden root weld. The imported model retains the rails and treadmill. Physics uses MuJoCo 2.3.7, implicitfast, 2 ms; neural updates are 0.2 ms. Joint stops prevent anatomical knee hyperextension, and a runtime guard rejects any negative knee flexion.

## Timing-first checkpoint and curriculum limits

`ORIGINAL_TIMING2_SAGITTAL.json` saves the selected original-controller configuration. With explicit continuous sagittal guidance and rails lowered 20 mm, its fresh 10 s repeat produced 10 left and 10 right unloaded, clear swing intervals with 19/19 alternating transitions. Minimum knee flexion was 2.307° left and 2.016° right. Peak whole-foot clearances reached 17.70 mm left and 13.20 mm right. Rails carried 85.70% of body weight. The root moved −0.157 m in world X while the belt moved −0.937 m: +0.780 m relative to the belt. This is supported treadmill stepping, not forward overground walking or independent balance.

The saved repeat was bit-identical in a fresh process. A 1 ms physics comparison retained the supported-stepping gate but not identical contact waveforms; this is qualitative gate persistence, not a converged contact-force validation.

The same parameters passed both deafferented and afferented air checks at 10/10 full knee excursions, period about 0.882 s, and contralateral phase about 0.48 cycles. A naive local-maximum counter counted a small extension recoil as a second left cycle in the deafferented case. Raw peak counts remain recorded; the corrected measure applies the same predeclared full-excursion hysteresis to both legs and all variants: ≥10° total range, arm below 25% excursion, count above 65%, then require return low before another count. It does not replace the ground contact/clearance gate.

**A complete Hunt four-part training sequence and corrected five-stage curriculum have not been completed.** These results establish a structurally audited port, timing calibration, air checks, and supported-ground diagnostics. Independent standing and harness-free walking are not established. Raw ground forces have high-frequency contact chatter and must not be presented as validated physiological GRFs.

## Physical contact guides allowing bounded wobble

The optional engineering apparatus is a 2 kg passive carriage on a finite ±1.9 m fore–aft slider. Eight collision boxes form lateral walls, front/rear wing posts, and overhead caps. There is no root weld, root-following geometry, or continuous restoring controller in guide-only conditions. The source horizontal rails keep their original length. Both 10 mm and 30 mm neutral lateral/fore–aft clearances were tested. Vertical cap clearance starts at the same value but increases as the wings settle; its measured supported bounding-box gap is reported separately.

`verify_guides.py` confirmed zero guide contacts at neutral, half-gap lateral displacement, and half-gap yaw displacement, then actual contact and correctly signed restoring reactions beyond the gaps. Geometry, contact assumptions, carriage inertia, and exact model XML are retained per run. A 5 N·m world-Y torque for 100 ms is triggered at the first left RG-flexor −59 mV upward crossing after 2 s; it is an explicit perturbation, not a balance controller.

All cases use the original 80/150 graph and physical muscle parameters, 20 mm lowered rails, and no continuous harness unless stated. Analysis starts at 1 s. The primary swing gate requires whole-foot clearance >1 mm for ≥80 ms, peak ≥10 mm, foot load <2.5% body weight, and opposite-foot load >5% body weight for at least half the interval. At least three swings per leg and 80% alternating transitions are required. This is a supported-stepping gate only.

| Shared neural time scale | Guide slack | Clear swings L/R | Alternating fraction | Double support | Net vertical assistance | Relative-belt stride L/R |
|---|---:|---:|---:|---:|---:|---:|
| 0.5, faster | 10 mm | 0/0 | 0% | 22.6% | 81.2% BW | No qualified strides |
| 0.5, faster | 30 mm | 0/0 | 0% | 22.9% | 80.6% BW | No qualified strides |
| 1 | 10 mm | 20/20 | 100% | 3.3% | 87.5% BW | −2.1/+2.7 mm |
| 1 | 30 mm | 20/20 | 100% | 4.6% | 87.1% BW | −12.9/−8.8 mm |
| 2 | 10 mm | 10/10 | 100% | 25.1% | 85.6% BW | 56.3/56.4 mm |
| 2 | 30 mm | 10/10 | 100% | 26.3% | 85.5% BW | 52.5/52.6 mm |
| 4 | 10 mm | 5/6 | 100% | 28.0% | 84.7% BW | 146.9/137.1 mm |
| 4 | 30 mm | 5/6 | 100% | 28.8% | 84.6% BW | 143.6/134.2 mm |
| 6, slower | 10 mm | 4/4 | 85.7% | 26.1% | 86.0% BW | 138.8/190.8 mm |
| 6, slower | 30 mm | 4/6 | 66.7%, fails | 27.4% | 85.8% BW | 138.2/117.1 mm |

The fastest case fails unloading/clearance. Scale 1 can clear and alternate yet has negligible or negative steady relative-belt stride; therefore clearance alone is insufficient evidence of useful forward gait. Slower is not universally better: scale 6 with wider guides fails the alternating criterion and has less consistent bilateral stride. Scale 2 gives matched bilateral strides and is retained as the conservative checkpoint; scale 4 is a promising longer-stride condition requiring additional repeats and basin checks.

At scale 2, 10 mm and 30 mm guides permitted approximately 2.70° and 4.31° yaw range, with guide contacts on approximately 16% and 7% of samples. The matched rail-only perturbation control reached about 6.03° yaw range and still passed the clearance gate; the continuous sagittal-harness comparator reached 1.85°. Thus these data support an effect on rotation and mechanical interaction, not a claim that the guide rescued an otherwise inevitable fall. A no-perturbation 10 mm guide control also passed 10/10 swings, with about 2.68° yaw range. Detailed contact/torque windows distinguish ordinary gait contacts from perturbation response.

Double support is an **offline** estimate from a causal 10 ms force average, contact-on threshold 2.5% BW held 40 ms, contact-off threshold 1.25% BW held 20 ms. Raw simultaneous-contact fractions and all forces are also retained. These filters do not feed the neural controller. Repeated within-stance recontacts are listed separately from qualified strikes after whole-foot swings. Stride is the displacement between successive same-foot qualified strike centers minus actual belt translation; signed negative stride is preserved. Center-of-pressure strike locations are also retained.

Logged signals include pelvis quaternion and declared Euler decomposition, world angular velocity, pelvis-only rotational kinetic energy using its actual inertia, whole-robot center of mass and velocity, signed guide/rail forces and torques, foot/toe forces, contact centers, foot clearance, carriage motion, and all neural/actuator states. Net vertical assistance includes signed rail, guide, and harness forces; it is not inferred from the harness command.

## Matched speed and numerical checks

With the current extended observers enabled, the same 4 s air simulation took 17.208 s in the scalar Python reference and 0.970 s in Numba: **17.75× acceleration and 4.13 simulated seconds per wall second**. Source/model loading, JIT warm-up, saving, and rendering are excluded. Both execute the same original 80-neuron/150-edge equations, 24 Hill/receptor elements, MuJoCo physics, and observers; this is not a comparison against AnimatLab or the SNS-Toolbox dense backend. Maximum differences were 3.27×10⁻¹³ mV in neural voltage, 7.94×10⁻¹¹ N in force, and 7.53×10⁻¹³° in angles. Supported contact runs have additional collision/reaction work and are timed separately in each RESULT.json. An earlier, lighter observer set was faster; its benchmark is archived separately rather than attributed to the current executable.

## Newer source remains a distinct unsuccessful branch

The supplied newer native file has 96 neurons/208 edges; its editable project has 96/204 plus different pathways, synapse types, and a left RG-extensor 2 nA tonic drive absent from the old native export. The GUI-authoritative source is the appropriate editable branch, not a silent native-208 substitute.

Installed AnimatLab contact telemetry proved that the supplied native polynomial C=20 actually injects 20 A at contact count 1. GUI prefixes were checked, not ignored: the project's polynomial scalar is Value=1, Scale=None, Actual=1, whereas the screenshot sigmoid with a nano-prefixed amplitude is a different gain. Separate GUI-authoritative 1 nA and 20 nA engineering companions preserve their explicit magnitude/offset choices. Neither optimum nor intended magnitude is established. Both corrected native runs fall despite finite voltages; the fall occurs inside the longitudinal rail span and cannot be explained solely by running off the track. Mirrored tonic-drive probes in MuJoCo have not produced a passing oscillatory supported-ground checkpoint for this branch.

## Provenance, artifacts, and next work

Raw data, actual model XML, hashes, and complete settings are under `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\`, which is gitignored. Principal summaries are `GUIDE_STUDY_V1.json`, `GUIDE_EXTREMES_V1.json`, `GUIDE_DEADBAND_VERIFICATION.json`, and `BENCHMARK_FINAL.json`. The original timing checkpoint is versioned beside this record. Every guide run has `PARAMETERS.json`, `RESULT.json`, `ASSESSMENT.json`, `GUIDE_ASSESSMENT.json`, and `trajectory.npz`; full rerun commands are in README.md. CLI flags explicitly override configuration-file values.

The previously rendered checkpoint GIF and still are at `runs\W2L_VirtualWalker\timing2_fresh_ground\SUPPORTED_STEPPING_REPLAY.gif` and `SUPPORTED_STEPPING_STILL.png`; their label is the continuous-guidance condition, not the new physical-guide condition. Short neural/mechanical figures are in that folder and `timing2_air_deafferented`. They show a few cycles rather than long compressed traces.

The public `vipamol/Rat-SNS-design-and-tuning/Biarticular_net` repository contains an actual 16-muscle/120-neuron biarticular precursor, including RF, BFP, and gastrocnemius. It is not authenticated as the exact 2023 release: its dates and multiple published conductances/integration coefficients differ. The separate BPA/literature agent cached original files and documented the differences. This is a useful audited starting point for a later biarticular branch, not evidence that a Deng 2023 replication has already run.

### Longer cadence comparison, completed

The 10 s slow-case failure contained relatively few intrinsic cycles, so the same fixed parameters were repeated for 20 s at 30 mm clearance. These three runs explicitly extend the treadmill motor end from its source default 10 s to 20 s; this engineering extension is saved as `belt_motor_end_s`, not silently applied to source baselines. The finite belt, original rails, and ±1.9 m carriage remain unchanged.

| Time scale | Actual RG peaks / median period | Clear swings L/R | Alternating fraction | Double support | Mean pelvis rotational energy | Qualified relative-belt stride L/R |
|---|---:|---:|---:|---:|---:|---:|
| 2 | 22 / 0.882 s | 21/22 | 100% | 26.77% | 0.09414 J | 51.0/51.1 mm |
| 4 | 10 / 1.760 s | 11/11 | 100% | 28.18% | 0.07985 J | 146.1/140.7 mm |
| 6 | 7 / 2.640 s | 8/11 | 77.78%, fails | 27.46% | 0.07218 J | 186.0/146.0 mm |

All three remained finite and nonhyperextending. Across them the conservative whole-foot-to-belt-edge margin stayed ≥0.403 m, wing-to-rail longitudinal margin ≥1.160 m, and carriage-limit margin ≥1.118 m. Thus the remaining slow-case coordination failure is not an arena-end effect. It is a declared-gate failure over seven observed RG peaks, not a statistical claim of a robust bifurcation. The extra right swing events also show that mechanical swing counts need not equal intrinsic RG peak counts.

In this comparison mean pelvis rotational energy decreases with slower timing, but double support is nonmonotonic and coordination does not improve monotonically. Energy 95th percentiles are 0.5066, 0.4420, and 0.4285 J, respectively; post-startup peaks are 1.3545, 0.9235, and 0.8673 J. These are pelvis-only values using measured angular velocity and the actual pelvis inertia, not total-body mechanical or metabolic energy. The 20 s simulations took approximately 9.50, 9.62, and 9.49 s of warmed wall time with all observers enabled. Data are in `GUIDE_LONG_COMPARISON_V1.json`.

### Newer GUI-authoritative paired guide test, completed

The guide grid uses the **original 80/150 circuit** and does not explain or repair the **newer GUI-authoritative 96/204** branch automatically. Four 10 s paired tests preserved that branch's supplied left/right tonic 2/0 nA, shared neural time scale 2, lowered rails, and existing 1 nA or 20 nA candidate contact gains. No weights were fitted, and no perturbation was applied.

| Newer candidate | Support | Minimum root Y | Yaw range | Roll range | Qualified swings L/R | Alternating fraction |
|---|---|---:|---:|---:|---:|---:|
| GUI 1 nA | Original rails | 0.94055 m | 31.85° | 11.44° | 4/3 | 50% |
| GUI 1 nA | Rails + 30 mm guides | 0.94050 m | 8.79° | 2.92° | 4/5 | 50% |
| 20 nA engineering | Original rails | 0.94073 m | 45.50° | 13.41° | 3/5 | 71.4% |
| 20 nA engineering | Rails + 30 mm guides | 0.94075 m | 8.79° | 3.72° | 3/3 | 80%, only four alternating transitions |

All four retain body height for this 10 s calibrated condition, including the rail-only controls. This is not a rescue of the earlier source-condition native fall: cadence and rail height differ from that baseline. Guides substantially reduce rotation, but **all four fail the full supported-stepping gate**. All four RG neurons have zero >1 mV prominent peaks after 1 s; their voltage ranges remain below 0.41 mV. Their observed leg clearances therefore do not establish functioning rhythmic RG coordination. In particular, a physical apparatus holding body height does not repair a latched neural oscillator. The source tonic asymmetry remains explicitly preserved rather than hidden by fitting. Full results: `NEWER_GUIDE_PAIRS_V1.json`.

A bounded isolated-core tonic-current A/B found that **both source 2/0 nA and tied 0/0 nA remain latched without any physical afferents**. Thus level-contact feedback is not required for this latch, and removing tonic asymmetry alone is insufficient. The loaded GUI-authoritative source already contains two phase-init pulses: **10 nA for 10 ms to L RG ext and R RG flx**, both starting at 0 s. They are retained exactly; no duplicate pulse is added. The 1 nA and 20 nA candidate intrinsic neural cores were checked for exact identity before sharing this deafferented test. Final per-RG leak, calcium, tonic, afferent, and individual synaptic driving currents are recorded.

At tied 0/0, the steady flexor balance is approximately −3.811 nA from the extensor interneuron, +2.086 nA calcium, +1.130 nA leak, and +0.595 nA direct extensor excitation, leaving essentially zero net current. The original versus newer left RG intrinsic parameters and initialization are identical; common RG synapses differ only by floating-point serialization. The newer adds V3 excitation of the contralateral extensor interneuron and c1 inhibition of the contralateral flexor. c1 input is negligible at the latched state.

With tonic currents **still tied at 0/0**, zeroing only the two mirrored V3→contralateral RG-extensor-IN output conductances restores all four isolated RG oscillations: seven prominent peaks each, approximately 0.980 s period over 8 s. The original80 positive control gives eight peaks per RG cell, approximately 0.882 s. This identifies an effect of the V3 relay **conditional on tied 0/0 tonic currents**. The additional isolated source 2/0 + V3-zero control fails the four-cell gate: left RG-extensor range remains only 0.217 mV with no prominent peaks, although the other three RG traces vary rhythmically. Thus V3-zero alone does not repair the supplied source condition. No source file is overwritten. Zero output conductance disconnects those two pathways mathematically; the V3 neurons and c1 connections remain present. Data: `runs\newer_rg_isolation_v1\`.

### Diagnostic afferent return, completed

The explicitly named **newer96 / GUI1nA / tonic0/0 / V3-output0** diagnostic retains c1 and the original startup pulses. The same neural parameters pass both bounded air tests and the supported-ground clearance criterion:

| Condition | Observed motion | Contralateral phase | Result |
|---|---|---:|---|
| Fixed-root air, deafferented, 8 s | 7/8 full knee excursions | 0.488 knee; 0.502 RG | Air gate passes |
| Fixed-root air, afferents returned, 8 s | 6/7 full knee excursions | 0.496 knee; 0.503 RG | Air gate passes |
| Rails + 30 mm physical guides, 10 s | 6/8 unloaded-clear swings; 11/13 alternating transitions | 0.503 RG | Supported-stepping gate passes at 84.6% alternation |

The RG median period is approximately 0.980 s deafferented, 1.080 s with air afferents, and 1.094 s on the ground. Ground minimum knee flexion is 1.124° left and 1.604° right; peak whole-foot clearance reaches 32.16 and 46.48 mm. The first qualified left swing occurs at 4.192 s, so the initial loading transient is not bilaterally successful throughout. The rendered 5–8.2 s window is explicitly a later settled segment, not the entire-trial behavior.

Whole-trial net vertical assistance is **35.176% BW**, almost entirely upward rail support (35.180% BW). Guide vertical force averages −0.00385% BW; upward and downward net-channel means are 0.00590% and 0.00975% BW. Continuous harness force is zero. A replay of contact geometry at all 5,000 saved configurations found **zero cap contacts**. Thus this reduction is not a large downward cap load canceling upward support. After 1 s, net assistance averages 32.617% BW. The rail reaction nevertheless contains sharp peaks up to 6.874 kN: these hard-contact waveforms remain uncalibrated and are not validated hardware loads. Full signed components and scope are in `SUPPORT_BREAKDOWN.json`.

The root still moves backward in world X with the belt. Relative-belt progress is +0.117 m over 10 s and +0.051 m over the 1–10 s evaluation window. This is actual supported alternating stepping with unloading and clearance, not overground walking or a completed weaning curriculum. The lower assistance belongs to this different, explicitly altered newer controller; it is not a same-controller weaning result from the original80 checkpoint.

The default source `v3_output_scale=1` is retained. A fresh default-path regression reproduced earlier source GUI trajectories bit-identically for q, neural voltages, forces, root position, and foot forces. No source graph or native project is overwritten. Canonical diagnostic neural configuration: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\NEWER_GUI1NA_TIED0_V3ZERO_DIAGNOSTIC.json`. Full bounded results: `NEWER_V3_RETURN_V1.json`; individual runs contain `V3_DIAGNOSTIC.json` and saved settings.

The separate native AnimatLab diagnostic companion is now independently verified at `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_NewerV3Guided_20261031\W2L_VirtualWalker_NewerV3Guided_Standalone.asim`, with a matching editable project. It retains 96/204 entries, including the two zero-output edges, c1, and GUI1nA contacts. Its 10 s native run exited normally, retained root Y ≥0.942147 m, had minimum knee flexion 2.278°/2.365°, and produced bilateral RG period about 1.074 s with phase about 0.505. Native foot unloading/force gates were **not evaluated**, so it is not assigned the MuJoCo ground-passing claim. See its `V3_NATIVE_VERIFICATION.json` and the native agent's separate record.

The diagnostic GIF, still, short neural/mechanical figure, actual model, and reactions are at `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer_GUI1nA\newer_v3_diag_v1_ground_guides\`. The render labels tonic0, V3-output0, 30 mm guide slack, and measured assistance. No background numerical job remains running from this package at handoff.

The versioned physical-guide configuration was separately rerun through the actual CLI in a fresh process. An explicit 5 N·m perturbation flag correctly overrode the undisturbed JSON setting, and q, angles, neural voltages, forces, root/COM positions, guide reactions, and foot forces were bit-identical to the earlier case. Mean external vertical force was 1.000412 BW versus 1.000536 BW expected from gravity plus COM velocity change: a mean closure discrepancy of 0.000124 BW. This checks reaction signs and bookkeeping, not pointwise contact-force convergence. Evidence: `FRESH_GUIDE_CLI_PARITY.json`.

Next: quantify a basin for the retained physical-guide checkpoint, then lower rail support in small explicit steps without changing the controller invisibly. Preserve failed cases. Complete the standalone standing/Hunt gates before presenting a full curriculum pass; retain nonhyperextending knees. Native physical-guide companions are completed; their separate record is `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided_20261031\EXPERIMENT_RECORD.md`, and the newer diagnostic has its independent verification under `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_NewerV3Guided_20261031\`.

## Session 1: Rustuna driver capability and current scope
Rustuna 0.2.1 is isolated in ignored tmp/rustuna_runtime_py311 under existing Python 3.11.5; myo Python 3.10.21 and its scientific packages are unchanged. The cheap 256-trial quadratic measured 0.0559 s versus Optuna 5.0.0's 1.8263 s, about 32.7x optimizer bookkeeping only, with different suggestion histories and solution quality. Same-backend seeded histories repeat exactly. The persistent myo protocol passed 100 exact toy requests in one worker, averaging 0.128 ms per round trip. A 64-trial ask/tell history resumed exactly after a trial-32 JSON checkpoint, including an enqueued zero-valued candidate; corrupted replay was rejected. A false baseline gate prevents any optimization.
Current user direction pauses rat/Gait transfer and synthetic mass changes. The driver is ready for the actual newer108 W2L contact-feedback ablations and bounded training owned by the physics agent, with 64 trials/900 s as the first-stage budget. The adapter must validate consumed symmetric contact gains, score physically unloaded/clear alternating ground steps, reject drag/nonfinite/hyperextension/apparatus-only support, and preserve failed cases. Best candidate is not a pass until fresh-process confirmation. No model optimization was launched in this capability check. Full API, benchmark caveats, protocol, and exact paths are in RUSTUNA_CAPABILITY_RECORD.md; raw toy histories and checks are under ignored runs/rustuna_capability and runs/rustuna_protocol_probe. Existing Optuna SQLite files were not opened or migrated.

