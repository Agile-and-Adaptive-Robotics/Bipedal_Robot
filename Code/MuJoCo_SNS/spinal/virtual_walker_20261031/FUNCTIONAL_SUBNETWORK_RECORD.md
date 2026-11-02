# Editable W2L subnetwork verification

Machine: **easteregg2**. One isolated review copy; canonical projects, saved controller results, the user-selected publication extraction, and live applications were not changed.

## Concrete editable deliverable

Open `D:\GitHub\Bipedal_Robot\tmp\w2l_functional_subnetworks_session1\editable\W2L_VirtualWalker_Newer.aproj` in AnimatLab. Its matching native files are `W2L_VirtualWalker_Newer.asim` and `W2L_VirtualWalker_Newer_Standalone.asim` in that same folder. They describe the **original canonical 108-neuron configuration**, not the subsequently fitted gentle-hip MuJoCo runtime recipe.

The existing RG, hip PF/MN, shared knee/ankle PF, knee/ankle MN/reflex, contact, and commissural subnetworks remain. Six new native subnetwork pages organize the existing homonymous II channels: `L/R Hip homonymous II`, `L/R Knee homonymous II`, and `L/R Ankle homonymous II`. Each is opened from its green subsystem node on the corresponding MN page. Each page contains two real signal chains:

`existing physical stretch receptor → existing current-conversion adapter → existing II neuron → OffPage reference to the existing homonymous MN`.

This is an editable hierarchy change, not merely a publication schematic. No neurons or biological connections were added. All 12 II neuron IDs, adapter IDs, chemical-edge IDs, physical source IDs, target MN IDs, neural intrinsic parameters, gain coefficients, and named synapse definitions are retained. New IDs identify only visual physical-source references, MN OffPage references, and subsystem containers. The four existing phase-related II channels remain separate and unchanged.

## Actual checks and results

- **108 neurons, 216 chemical connections, 47 native synapse-type definitions, and 56 adapters** match the source by ID and all serialized parameters. Physical body payloads also match.
- All **24 saved diagram pages** parse independently. Their node/link counts and interleaved AddFlow Org/Dst handles agree with their actual logical endpoints.
- The **92 local edges using 19 historical named types** retain exact reference names, equilibrium potentials, maximum conductances, presynaptic thresholds, and saturation levels from the historical uniarticular LH template at git `23b0e785`.
- A fresh native baseline and the reorganized native model both completed **2.02 seconds**, with 1,001 finite recorded samples through 2 seconds. All **1,001 × 146 recorded values were bit-identical**, with maximum absolute difference **0.0**. The neural kernel, feedback, muscle mechanics, and support were not changed for these tests. Source native physics is 1 ms and neural integration 0.2 ms; this is not a new cross-engine convergence result.
- Baseline root y stayed between **0.960699 and 0.982885 m**. Both RG-flexor traces varied within approximately −63 to −52 mV. These short runs verify unchanged execution after reorganization; they are not a new curriculum or ground-stepping acceptance test.
- Source connection evaluation order was retained in the standalone files so that diagram reparenting did not introduce a floating-point summation-order change.

## Pre-existing editable-project defects repaired in the copy

The LH and RH overview pages each contained duplicate logical definitions of the RG-extensor and knee-PF-extensor OffPage objects: **four duplicate proxy definitions** in total. Each duplicate pair had only one actual drawn glyph. Their link caches were merged, the redundant definitions removed, and the header counts corrected; all original unique object IDs and effective neural endpoints remain. Two duplicate RG proxies also carried unsupported zero-valued `TonicStimulus` metadata. Those fields belonged to OffPage references, not neurons; no actual neuron tonic current was changed.

The source physical-receptor proxy names sometimes described legacy Ia labels even when their body reference supplied the new II channel. New proxy labels therefore use the **actual referenced physical body name** and `II length signal`, rather than copying a misleading old label. Physical binding and gain coefficients remain unchanged. The display-label and duplicate-object ledgers identify every such repair.

Receptor-name suffixes also do not reliably identify anatomical action. For example, the left hip extensor channel uses the physical object named `hip_L_flx_rec`. The original channel-construction ledger matched this receptor to the `hip_L_ext` actuator by its full attachment-route signature, with zero route mismatch. The saved-page preview deliberately exposes the exact physical name rather than silently relabeling that body or treating its suffix as evidence of a crossed connection. The mapping remains recorded in `Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\canonical_12II\BUILD_LEDGER.json`.

## Saved-page previews and verification files

Workspace: `D:\GitHub\Bipedal_Robot\tmp\w2l_functional_subnetworks_session1`.

- `saved_page_previews\L_Hip_homonymous_II.pdf` / `.svg` / `.png`: an actual saved functional subnetwork page; five other joint/side pages accompany it.
- `FUNCTIONAL_SUBNETWORK_VERIFICATION.json`: ID/parameter comparison, all drawing checks, native run results, and exact trace comparison.
- `HISTORICAL_NAMED_TYPE_CHECK.json`: per-edge comparison of all 19 reference type definitions.
- `SUBNETWORK_REORGANIZATION.json`: every moved channel and its physical source, adapter, neuron, chemical edge, and MN reference IDs.
- `DUPLICATE_PROXY_REPAIR.json` and `PHYSICAL_PROXY_DISPLAY_LABELS.json`: explicit repair ledgers.
- `SOURCE_MANIFEST.json`: original source hashes and isolated-copy locations.
- `baseline_native\` and `reorganized_native\`: owned-run logs, raw traces, and recorded process IDs. Both simulator processes have exited.

The previews render stored positions and route points rather than substituting a generic neural layout. The muscle-spindle image is the actual bundled `AnimatGUI.MuscleSpindle.gif`, extracted from the installed AnimatGUI assembly; it is a raster icon within the vector figure. These are reconstructed saved-page views, not screenshots of a verified GUI opening.

The isolated copy and raw traces are ignored by Git. No additional active model family was created. The two user-supplied model folders remain intact. This review copy has not been promoted over the canonical project, whose duplicate-proxy defects therefore remain pending that explicit integration.

## Scope and next bounded check

Interactive GUI opening has not been exercised because the user's live window must remain undisturbed. XML, graph, drawing-handle, native export, and native execution checks have passed. The model uses **graded NonSpiking neurons**. A spiking implementation would require explicit threshold/reset, spike-event synapses, and sensory/MN-interface verification; changing a class name alone would not make a valid spiking counterpart.

## Public rat precursor: editable/runtime reconciliation and measured source problem

Originals remain unchanged in `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat`. The public repository is `https://github.com/vipamol/Rat-SNS-design-and-tuning/tree/main/Biarticular_net`; the archived project is a public precursor, not a verified reproduction of the final Deng 2023 model.

The supplied editable project contains **132 enabled graded neurons and 252 chemical connections**. The supplied native standalone contains **120 neurons and 238 connections**. All 12 missing neurons and 14 missing edges belong to the additional subnetwork named `1`. Its physical references share the actual left biceps-femoris-anterior/posterior muscles and receptors, so the additional page cannot be dismissed as disconnected drawing metadata. The original 120-node primary circuits agree by effective endpoint and type; all 21 saved pages pass XML and drawing-endpoint checks.

An isolated GUI-authoritative export now exists at:

`D:\GitHub\Bipedal_Robot\tmp\w2l_functional_subnetworks_session1\rat_verification\editable\Biarticular_net.aproj`

Matching `Biarticular_net.asim` and `Biarticular_net_Standalone.asim` sit alongside it. These include the actual 132-node/252-edge editable graph; no neuron, synapse, or gain was invented or fitted. `SUPPLIED_120_REFERENCE.asim` preserves the supplied 120-node native model separately within the verification workspace. Original source files and their physical definitions remain untouched.

Both native versions complete 2.02 s and produce 1,001 finite recorded samples, but **both fail the broad voltage-range diagnostic** of −500 to +200 mV. The supplied 120-node source reaches 304,866,000 V in `RH_SO- Ext Ib`; the full 132-node export reaches 337,349,000 V in `RH_IP-    FLX Ib`. An ordinary-looking RG trace and a successful simulator exit therefore do not establish valid sensory feedback. The issue already exists in the supplied standalone and was not introduced solely by exporting the additional page.

### Direct physical-input/current measurement

A separate 0.52 s run of the unchanged supplied 120 model records physical soleus-Ib output, target-neuron `AdapterCurrent`, and target `MembraneVoltage` through 0.5 s. It changes only chart instrumentation.

| Item | Left soleus-Ib path | Right soleus-Ib path |
|---|---|---|
| Adapter ID | `05f581fc-44c1-417e-b9b5-193bd1b422e4` | `376084f4-da7a-4758-be65-b4f7564a2a7f` |
| Native polynomial | C=7.9e−12, D=0 | C=1, D=0 |
| Limits | Disabled | Disabled |
| Physical Ib range, source units |0 |0.0177879–572.383 |
| Actual adapter-current range, A |0 |0–572.321 |
| Neuron-voltage range, V |−0.06 |−0.06–237,162,000 |

The right GUI coefficient is explicitly `<C Value="1" Scale="None" Actual="1"/>`, with `TargetDataTypeID=ExternalCurrent`; the native coefficient is also 1. This specific measured path has no effective nanoampere conversion. The displayed `Scale` and `Actual` values, rather than just a GUI number, were checked. The left physical signal happens to be zero in this short probe, so its gain is a source contrast, not a validated left-side physiological calibration. No automatic 1e−9 correction or mirrored gain replacement has been applied.

Artifacts under `D:\GitHub\Bipedal_Robot\tmp\w2l_functional_subnetworks_session1\rat_verification`:

- `RAT_EDITABLE_AUDIT.json`, `GUI_ONLY_NEURONS.json`, `GUI_ONLY_CONNECTIONS.json`, and `PROTOTYPE_PHYSICAL_BINDINGS.json`: exact source differences and physical bindings.
- `FULL_GUI_EXPORT.json`: GUI-to-native synchronization ledger.
- `SOURCE_VOLTAGE_BOUNDS_COMPARISON.json`: both complete voltage-range results.
- `ADAPTER_UNITS_MEASUREMENT.json` and `adapter_measurement\MEASUREMENTS.npz`: exact coefficients, channel bindings, and measured current/input/voltage samples.
- `native_smoke\` and `full_GUI_native\`: separate source 120 and exported 132 traces and logs.

All owned native runs have finished. The rat export is an editable diagnostic deliverable, not a newly validated locomotor controller. A subsequent correction should audit all receptor-to-current gains, preserve an explicit before/after ledger, and verify actual input/output bounds before training; the RH soleus measurement alone does not identify calibrated gains for every channel.

## Durable W2L checkpoint

The verified W2L editable checkpoint is preserved outside the disposable workspace:

`D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\W2L108_verified\editable\W2L_VirtualWalker_Newer.aproj`

Its matching native files, forms, bitmap dependencies, six saved-page vector previews, verification ledgers, and `MANIFEST.json` accompany it. The manifest hashes 58 curated files. Only embedded absolute workspace paths changed during relocation; native XML equivalence with the verified copy was checked after normalizing `ProjectPath`. The earlier bit-exact paired native test remains the execution verification; no additional physical or neural change was made for archival.

## Rat: complete adapter audit and experimental bounded conversions

All **44 physical-to-neural adapters** were instrumented in the full 132-node model. **19 right-side adapters and six prototype-page adapters** retain identity polynomials with C=1, `Scale=None`, and `UseLimits=False`. This affects Ia as well as Ib and II. The remaining 19 left-side adapters contain nontrivial calibrated slopes and offsets. Six prototype nodes labeled Ia/Ib/II actually all read the physical `Ia` getter; their actual channel bindings were retained rather than silently inferred from names.

The installed engine's matching public SDK source establishes:

- `LinearHillStretchReceptor.cpp:172–173`: Ia = discharge constant × series-element length, and II = discharge constant × parallel-element length.
- `LinearHillMuscle.cpp:498`: Ib = discharge constant × current tension.
- `Neuron.cpp:1287–1288`: incoming adapter current is converted from A to internal nA, and `AdapterCurrent` retains the A-valued measurement.
- `Neuron.cpp:1106–1107` and `1235–1241`: applied nA current is divided by the native `Size` before entering the relative membrane-voltage update in mV.

The local source copies are under `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs`. The [official SDK implementation](https://www.animatlab.com/SDK/_linear_hill_stretch_receptor_8cpp_source.html) also documents the stretch-receptor calculation. These are model-specific rate signals, not force values that can be passed straight to a current input. The supplied 2023 manuscript describes the feedback dependencies and manual feedback tuning, but does not supply a verified numerical conversion for these unfinished identity maps; no literature calibration was claimed.

Three explicitly experimental 132-node/252-edge editable/native pairs were tested, all below the existing `rat_verification` workspace. Original source projects, source neuron properties, synapse types, chemical connections, and physical definitions remain unchanged. These diagnostic folders are not additional active model families.

1. **`experimental_bounded_maps`**: physical-domain normalization for all 44 adapters. Ib limits derive from the source maximum tension and discharge multiplier. Ia limits derive from source series-element resting length and the maximum elastic extension `MaximumTension/Kse`. II limits use the source minimum parallel-element length and the length-tension envelope. Their current endpoints derive from each target's actual outgoing synapse thresholds/saturation and native current/voltage equation. Explicit input/output limits prevent values outside these engineering envelopes from producing unbounded current.
2. **`experimental_source_LH_bounded`**: preserves all 19 original LH polynomial slopes/offsets, adds limits obtained by inverting those polynomials at the target synaptic current endpoints, and retains the first candidate's corrections to the 25 identity maps. Thus it tests the effect of replacing existing LH calibration separately from correcting identity conversions.
3. **`experimental_prototype_drive_ablation`**: retains the second candidate and disables only the two additional prototype physical motor adapters listed below. All 132 neurons and 252 chemical connections remain present and enabled as before.

Every `MAP_LEDGER.json` records the actual physical channel, coefficients before/after, physical-domain derivation, native size/rest potential, outgoing type thresholds/saturation, chosen current range, and export ledger. All adapter `Enabled` flags and numerical coefficients match the GUI and both identical native exports. A reusable isolated-export helper omission was also corrected: it previously retained old native adapter `Enabled` flags instead of propagating the GUI property. Existing source-enabled adapters are unaffected by that correction.

### Additional prototype motor drivers overwrite existing muscle control

`MuscleBase.cpp:233–238` assigns `m_fltVm=fltInput`; multiple motor adapters do **not** sum into the muscle. The additional enabled prototype adapters address the same actual left BFA/BFP muscles as the primary controller, creating competing assignments whose execution order matters:

| Prototype driver disabled in diagnostic only | Adapter ID | Existing physical target |
|---|---|---|
| HE MN | `29d70e3d-513c-476a-bed6-659ebdb677c4` | Left BFA, `08b9c1bc-15cc-4dd4-95cc-4738f99a49d2` |
| Biar1 MN | `752936d8-3a5c-4a01-ab0a-587fe8603d22` | Left BFP, `35e368dc-4be8-456d-8baa-ea7f5a5c9e51` |

The supplied 120-node standalone produces three left RG-extensor peaks from 0.5–2 s, with median period 0.481 s. The uncorrected GUI-authoritative 132-node export already loses those cycles. Both bounded-map candidates remain latched while the extra prototype drivers are enabled. Disabling only those two adapters in the otherwise identical second candidate restores the left layer rhythm, providing a direct intervention test of the competing-driver mechanism.

### Native results and remaining limitation

Each experimental case completed **5.02 s**, recording **2,501 finite samples** through 5 s. All 132 voltages remained between **−100 and −30.025 mV**, instead of the original extreme values. The first two cases fail the layer-rhythm screen. With the two prototype drivers disabled:

- Left RG-extensor/flexor show eight peaks after 1 s, periods approximately **0.484/0.486 s**, and voltage ranges **4.377/4.104 mV**. Hip PF and MN layers also cycle.
- Right RG-extensor/flexor show eight smaller peaks using a stated **0.2 mV prominence** criterion, periods approximately **0.486/0.484 s**, and ranges **1.826/1.662 mV**. Right hip PF modulation is weaker still, approximately **0.68–0.81 mV**, while its MN outputs cycle strongly.
- A separate stricter 2 mV burst-prominence screen therefore reports zero right RG peaks. That means the screen fails; it does **not** mean the right traces are flat or contain no rhythm. Existing cross-side pathways and strong MN activity prevent interpreting periodic MN output alone as an independently functioning right RG/PF controller.

The source pelvis remains fixed. These tests establish bounded sensory conversion and identify the competing physical drivers; they do not establish a calibrated final rat controller or a locomotion acceptance result. No experimental candidate was promoted over the original source.

Additional files in the verification workspace:

- `ALL_SENSORY_ADAPTER_AUDIT.json`: all 44 original conversions, GUI scales, actual source ranges, current ranges, and neuron ranges.
- `BOUNDED_MAP_COMPARISON.json`: three candidate summaries and GUI/native parity checks.
- `experimental_prototype_drive_ablation\WEAK_RHYTHM_ANALYSIS.json`: explicit small-amplitude cadence measurements.
- `RAT_BOUNDED_MAPS_LAYER_COMPARISON.pdf`, `.svg`, and `.png`: inspected bilateral RG/hip-PF/hip-MN comparison over the same 0–2 s interval, with original 120 and GUI 132 baselines overlaid. This is a diagnostic figure, not an active dissertation asset.

All owned simulator processes have exited. Further work should audit the right RG/PF drive and source cross-side bypasses with the bounded conversion held fixed, before fitting or promoting a bilateral controller.

## Subsequent bilateral diagnostic

The next bounded pass is recorded separately in `RAT_BILATERAL_DIAGNOSTIC_RECORD.md` alongside this file. It preserves the rat editable checkpoint externally, verifies matched recorded-current neural replay, restores full-size right RG/hip-PF modulation by explicit operating-range registration, and identifies missing executable RG→PF connections in this public precursor. Right distal PF and mechanical motion do not yet pass. The record includes an inspected aligned neural/joint figure, a reviewable 12-edge reconstruction proposal, and a scoped 255.9× compiled-versus-normal neural calculation benchmark. No proposed connections were applied, and no complete coupled MuJoCo result is claimed.
