# Four-command and central-feedback preparation — Session 1

Machine: easteregg2. 4 November 2026. Ben authorized the new four-command interface, missing central sensory pathways, and contact-fanout checks before retraining. MATLAB parallel work takes priority: **zero neural steps, zero physics steps, and zero optimizer trials** were executed. Static network construction/compilation and graph fault-injection checks were performed sequentially. Rustuna was imported, but no study was created. No tendon, damping, CAD, or mechanical-model changes.

## Current files

Code factory: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\profile.py
Static preparation: same folder prepare.py.
Deferred bounded optimizer: same folder train.py.
Persistent scientific adapter: same folder replay_adapter.py.
Physical model whose actuator roster is used: D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml. Preparation reads XML names; it does not load/integrate the body. Physical evaluations must use the repaired, frame-adapted, knee-protected condition, not run this raw converted XML directly.

All experiment records: D:\Users\Ben Bolen\ControllerPreparation_20261104. Subfolders hunt, s3k, w2lvar, syn6 each contain PROFILE.json, STATIC_CHECK.json, WIRING_CHANGES.json, SOURCE_HASHES.json, and BUILD_SETTINGS.json. hunt additionally has HUNT_OPTIONS.json, including the saved symmetric motor-weight source. These are revised neural candidates, not standalone new mechanical MJCF files or trained winners.

## Static results

| Family | Neurons | Inputs | Actual chemical slots | Ia/II to PF (each) | Ib to PF | Result |
|---|---:|---:|---:|---:|---:|---|
| hunt | 1029 | 424 | 9424 | 104 | 14 | Static pass |
| s3k | 420 | 386 | 1604 | 162 | 26 | Static pass |
| w2lvar | 892 | 386 | 3296 | 106 | 14 | Static pass |
| syn6 | 798 | 386 | 3798 | 318 | 18 | Static pass |

Every family has 18 Ia->RG routes, 18 II->RG routes, and 14 Ib->RG routes. Every explicitly selected load source (bilateral soleus, three vasti, and three gluteus-maximus parts) has Ib access to RG and its own extensor PF channel(s). All added parameters are mirrored. All new edges match actual compiled slots. No duplicate endpoints or same-sign contact fanout remain. Injected excitatory and inhibitory errors to both antagonistic RG or PF channels are detected in all four families. These checks do not establish cadence, stability, body movement, or a physical curriculum pass.

The stronger per-source gate caught a soleus-to-ankle-PF omission in the first Hunt draft. A name-prefix test had classified PF_FUNC as PF_F. Corrected phase recognition and the per-source Ib check are retained as regression gates. This was a preparation error; no trial used the incomplete version.

## What changed

Four channels: RG_F_drive_l, RG_E_drive_l, RG_F_drive_r, RG_E_drive_r. Commands are normalized 0..1, encoded as 0..5 nA into 50 ms command-relay neurons with 0..5 mV output range. The initial conductance budgets use a declared 0 mV postsynaptic reference and 1.2 nA excitation budget. Flexor command also reaches the opposite flexor at one-third ipsilateral weight and the ipsilateral flexor inhibitory relay; available V3-E/F targets receive scaled input. This is a human-controller adaptation, not full reproduction of the cat network or literal transfer of its weights.

Shared DRIVE/legacy alpha/gamma central and MN projections are removed in the new candidate, and obsolete input currents are cleared. Hunt's added HUNT_INIT tonic inputs are cleared. Historical builders, checkpoints, and controller.py source remain untouched. Original controller.py SHA256: 00112d00d7e65c700c1c89e651f57354588d88335386209e5803b67f7154b57c.

Ia/II feed PF channels that actually recruit the corresponding muscle. Explicit monoarticular timing-source sets project to appropriate RG half-centers, avoiding a single global flexor/extensor label for biarticular muscles. Ib reinforcement is restricted to selected stance extensors, and is gated by measured LOAD input in a mechanical run. Each afferent-type/target budget is normalized for fan-in and mirrored. Runtime sensory gain is divided by 1+k*alpha_side using the filtered flexor-command signal; applying this to PF as well as RG is a declared project extension. Existing local MN/reflex circuits are retained.

Heel contact becomes a rearmed loading event with hysteresis and decay, not a sustained reset signal. It excites ipsilateral InE and contralateral RG-F; at PF it excites only the extensor-side inhibitory relay, suppressing flexor PF through lamination. Toe/whole-foot load remains distinct from the heel event; sustained toe contact does not supply a second central phase reset in this profile. Contact encoder scaling is explicitly initialized. The graph audit stops at central entry, so it does not mistake normal reciprocal inhibition inside a half-center for a second independent contact input.

Syn6's original heel input excited both ipsilateral RG interneurons. The before/after ledger retains that flag and the removed endpoints. Its original left/right synergy membership also differed. The new profile mirrors the **right PF dynamics, incoming circuit, PF->MN readout, and PF->IaIN gates** to the left, preserving the archived synergy_basis.npz and original left basis. This is an explicit new symmetric rank-six candidate, not a claim that the old NMF decomposition was symmetric or a change to ZCode's rank-eight work.

## Deliberate scope limits and remaining gates

This is the graded branch. Spiking re-derivation and new-network dynamical equivalence are not claimed. Initial gains/time scales need fitting. Generic robustness, source-muscle fidelity, and full 23-DoF stability are not proved by central coverage.

Phase feedback follows existing readouts: s3k's twelve primarily adductor muscles have no phase-PF readout; W2L/Hunt's six trunk actuators use separate engineering/postural inputs; syn6 retains the six source-pruned muscles without synergy drive. These are listed in each WIRING_CHANGES.json, not silently assigned to arbitrary gait phases. The requested **postural premotor replacement remains pending**; it is a gate before new stage-3/4 claims. When composing that full controller, repeat the earlier physical gates rather than inheriting a different topology's passes.

## Training prepared, not executed

The launcher uses an isolated Rustuna 0.2.1 Python-3.11 runtime copied from disposable repo tmp to D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311 (4.8 MB). Scientific work remains in D:\Anaconda\envs\myo\python.exe through the existing persistent JSON-lines worker. No environment was upgraded or installed.

On easteregg2, this command prints a plan and starts no worker or optimization:

    & 'D:\Anaconda\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\train.py' --family hunt --prepared 'D:\Users\Ben Bolen\ControllerPreparation_20261104\hunt' --phase A

After Ben releases the computer, add --execute --computer-free. Substitute s3k, w2lvar, or syn6 in both family and folder. Default block: 16 trials, 900 s maximum; run one initially, measure actual timings, and do not start parallel heavy work blindly. No background watcher or automatic restart was installed. The first real preflight performs the first neural replay of the new topology; **that runtime check is still pending**.

A fits shared intrinsic tau/ionic gain for twice the target period, with only a brief startup pulse. B freezes that result and fits shared RG sensory gains. C freezes timing and performs a bounded six-sagittal-group feedforward-peak recalibration. D freezes feedforward and refines local Ia/II/Ib gains against desired command waveforms. This initial small C block is not an unrestricted refit of every PF coefficient or every 3-D motor primitive. Later blocks require an accepted FINAL.json from earlier blocks. Period, cycle detection, left/right alternation, and PF activity gates prevent a frozen result from passing on low command error alone.

The adapter replays the existing selected IK-derived length/velocity/force targets. It integrates the actual compiled neural network, but imposes expected afferents rather than integrating mechanics. Its MN/activation and kinematic-reference information is not new simulated joint motion, EMG, physical foot loading, or walking. New local results are stored under <family>\training\<A|B|C|D>, with compact trial rows, best parameters, manifests, source hashes, and explicit failure diagnostics. Existing folders are preserved; this first launcher does not automatically resume an interrupted study. Fresh-process finalists, full motor fitting where needed, and the corrected body curriculum remain subsequent gates.

## Dissertation revision

Authorized separate source workspace: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Final\main.tex. Chapters 6/7 and Appendix D have the first revision; ProofFinal and live Overleaf are unchanged. Final\BUILD_REVIEW.md describes coverage and conservation. Built PDF: D:\Users\Ben Bolen\ControllerPreparation_20261104\final_tex_build_v1\main.pdf (310 pages).

Important correction: the historical human drive-weaning script explicitly kept shared DRIVE=2.0773198435110425 while scaling only supplementary left/right channels to zero. Final Chapter 7 and Appendix D now describe that scope. This finding corrects the prior 'zero external neural drive' interpretation, not the archived result numbers.
