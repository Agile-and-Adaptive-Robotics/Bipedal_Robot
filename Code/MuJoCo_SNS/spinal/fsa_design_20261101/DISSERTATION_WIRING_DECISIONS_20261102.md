# Dissertation controller wiring: selected decisions

Session 1, **easteregg2**. This is a paper-to-code design disposition, not a manuscript revision or a new body result. It prioritizes the dissertation's actual recipes, distinguishes later standing snapshots, and identifies which new modules have really been composed. Bounded body pilots authorized separately by Ben and the campaign owner do not inherit results from these electrical checks.

## 1. Which historical graph are we describing?

The full dissertation s3k recipe **does contain central muscle feedback**. A later 410-cell standing export cannot stand in for it. Counts below are current-code reconstructions of saved recipes, not proof of the exact historical executable: the original trial logs lack source/model hashes.

| Reconstructed dissertation condition | Cells / nonzero compiled pairs | Direct central Ia / II / Ib pairs | Local relay populations and disposition |
| --- | --- | --- | --- |
| s3k full | 878 / 8,340 | 152 / 368 / 216 | 92 each IaIN, RC, IIX, IIIN, IBIN, and 10 IBEXC. Retain this recipe as the historical full-controller reference. |
| s3k combined prune and two pruned retunes | 868 / 8,084 | 152 / 368 / 0 | Local Ib/reflex circuits remain. This is removal of direct central Ib, not removal of every Ib signal. |
| w2lvar full | 888 / 7,868 | 152 / 152 / 216 | The same five 92-cell relay classes and 10 IBEXC; central feedback is already present. |
| syn6 full | 794 / 3,037 | 0 / 0 / 0 | IaIN, IIX, IIIN, IBIN, and IBEXC remain; no RC population. Contact input is enabled in this recipe. |

Source chain: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\prune_matrix.py`, its actual family curriculum setter, saved parameters, then the final compiled matrices. Full s3k loads `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260923\s3k_trial34_full_params.json`; full Syn6 loads `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\curriculum_syn6_stage5.json`. Historical reconstruction and exact evidence paths are in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\GAIT_DISSERTATION_RECIPE_REVIEW.md` and `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DISSERTATION_TRIAL_PROVENANCE.md`.

Selected correction: call an omission **historical** only when the actual recipe lacks the implemented pathway. The raw/default and later standing s3k graphs have no direct central Ia/II/Ib; that statement is false for full dissertation s3k. Syn6's direct central Ia/II/Ib omission is real in its full recipe, while its local reflexes remain. Its saved heel/toe gains are 0.12312339160142002 / 0.10178084469071401, not the zero-contact standing caller.

## 2. Heel input: first-hop excitation is not the final RG effect

In these graded Gait builders, **InE inhibits RG-F and InF inhibits RG-E**. The compiled current is

`I = gmax × clip((Vpre − elo)/(ehi − elo), 0, 1) × (Erev − Vpost)`.

At the stated operating voltages, +8 mV reversal is excitatory and −5 mV inhibitory; the current's driving force must still be checked. The raw-template RG inhibitory reference is 4 µS; actual recipe multipliers change it (for example, full dissertation Syn6 has 4.594286048130957 µS). Use final fields for the selected condition. A positive input into an inhibitory relay does not excite its final RG target.

The tags describe reversal relative to rest, not a voltage-independent current sign. A rest-relative +8 mV synapse supplies **negative** current when a NaP half-center is above +8 mV, and zero at +8 mV. The two-hop table describes the declared operating regime, not a complete phase-response curve. Likewise, the 0.5 nA normalized budget is evaluated at **2.5 mV**, not delivered as a fixed current at every voltage. Preserve actual reversals and examine recorded driving forces; do not alter reversal potentials merely to force the desired sign.

| Existing policy | Actual first hops | Effective two-hop direction |
| --- | --- | --- |
| Full-rules s3k | Heel excites ipsilateral InE and inhibits ipsilateral InF | Suppresses ipsilateral RG-F; can disinhibit ipsilateral RG-E. |
| w2lvar | Heel excites ipsilateral InE and **contralateral InF** | Suppresses ipsilateral RG-F and contralateral RG-E. The source comment claiming partner-flexor suppression is inconsistent with these endpoints. |
| syn6 | Heel excites ipsilateral InE and ipsilateral InF | Both ipsilateral RG half-centers receive additional inhibitory influence. This is a source policy to test, not a selective-extension mechanism. Its separate PF relay branch also remains. |

For full dissertation s3k, the final heel→InE conductance is **0.24354654619098992 µS**, whereas heel→InF inhibition is **0.37530421840765804 µS**. The later PF-building block writes heel→InE again using `ib_e_central`; pinned SNS assigns one final matrix entry per source/target pair. Therefore the configured `heel_rge` is not the final conductance of both paths. Likewise, an optional excitatory heel→InF definition cannot be described as safely coexisting with an inhibitory definition on the same compiler entry. Inspect the final matrix, not just source declarations.

Selected decision: retain historical graphs; do not silently choose an alternative sign or add direct contralateral RG-F excitation on top. The optional A/B policies remain distinct hypotheses in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\HEEL_STANCE_POLICY_REVIEW.md`. Their isolated sign proofs do not show a full recurrent controller's phase response.

**Ben's subsequent explicit Syn6 decision:** the repaired candidate retains **HEEL_s→InE_s excitation only**, removing HEEL_s→InF_s. It does not add inhibition onto InF or select either earlier A/B policy. Historical source and result files remain unchanged. The named opt-in implementation is `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\syn6_single_heel_in.py`. It zeroes exactly two compiled conductances, preserves all other fields and heel→PF routes, and rejects an existing direct/crossed reset rather than stacking policies. Apply before accelerator construction and reverify any later full-graph design composition. This decision is separate from Hunt's masked inherited IN branches and retained direct crossed event reset.

**Held contact versus event:** the generic runner writes continuous load-dependent contact current and can add onset transients; the presence of an onset branch does not itself remove the tonic component. The selected current Hunt1033 caller is different: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py` replaces its heel inputs with a hindfoot landing-event encoder when configured. Its packaged condition uses cross-contact conductance 6 µS and event amplitude 2 nA, armed after whole-foot unloading. It has hysteresis, a capture window, and finite decay, rather than a reset sustained throughout stance. The new integrated profile preserves this actual direct heel→contralateral RG-F route and **does not stack heel-policy B** over it. Raw heel and toe were historically duplicated contact proxies; regional contact and independent once-counted force measurements must be distinguished from that old logger.

## 3. Ia, II, and Ib: select functions by architecture

The Shinohara **2025 preprint**, not a 2019 paper, uses positive flexor length/velocity feedback and extensor force feedback into specified central populations (equations 9–11). Its muscle signals, nonlinear velocity encoding, and dimensionless coefficients multiplying a 10 nS synaptic scale are not our numerical SNS gains. [Primary preprint](https://doi.org/10.1101/2025.11.11.687930). Cached text: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\shinohara_2025_biophiv687930_fulltext.txt`.

- **w2lvar already implements** flexor-family Ia/II→same-side RG-F, InF, hip-F PF, and knee-F PF; extensor-family Ib→RG-E, InE, hip-E PF, and knee-E PF. Conductance is divided by source-class membership: 0.53/19 and 0.5/27 in the reviewed template. There is no basis for calling these hundreds of individually full-strength copies. This block is not a separate ankle, lumbar, frontal, or transverse PF controller.
- **s3k has optional** F Ia/II→F and E Ib/II→E routes. They are active in the dissertation full recipe, independent of whether some later caller sets them to zero. Local IaIN/RC enables are also independent of `full_rules`; one flag does not certify complete lamination.
- **syn6 has optional extensor Ib central/LBIN routes**, but not the stock central flexor Ia/II or extensor II block. An unused similarly named gain is not an implemented route. Absence alone does not make it a failed Di Russo reproduction, because Di Russo uses a different rhythm architecture.
- **Local autogenic Ib inhibition and positive central/load feedback are different paths.** An extensor force signal can inhibit its own MN locally while supporting an extensor-phase central population. A phase/load coincidence gate does not close bypass routes; w2lvar retains 216 direct central Ib edges outside IBEXC.

Selected new rule: the optional monoarticular overlay covers only 22 reviewed muscles, derives actual supported PF targets from positive source PF→MN projections, and does not project indiscriminately to all PFs. Flexor Ia/II and extensor Ib mechanisms have a declared source precedent; added extensor II is an explicit engineering extension. Biarticular and non-sagittal functions need their own signed mechanical maps. This overlay's 0.5 nA per-target/type current budget at 2.5 mV is an engineering reference, not measured physiology. It remains a separately verified proposal, not an installed repair of the historical trials.

## 4. Syn6 is not Di Russo's five-primitive generator

| Function | Actual production Syn6 K6 | Di Russo et al. (2023) |
| --- | --- | --- |
| Rhythm | Persistent-Na RG half-centers with inhibitory relays | Two coupled phase oscillators with heel-strike resets |
| Primitive generation | Six NMF-derived channels, E/F classification, phase-ranked PF time constants, and inhibitory lamination | Five equally spaced analytic raised-cosine primitives |
| Motor map | Nonnegative spatial basis transformed by Szczecinski 2017 Eq.18; other inhibitory/reflex circuits are separate | Signed positive/negative primitive→MN weights plus local reflex/Renshaw circuits |
| Fitting evidence | One computed inverse-activation cycle per side; production K6 basis and separate dynamic backsolve | Staged imitation, physical optimization with effort/stability terms, then removal of imitation |

[Di Russo primary paper, §§2.2.1–2.3 and equations 6–14](https://doi.org/10.1088/1741-2552/acfdcc); [author-hosted full text](https://infoscience.epfl.ch/server/api/core/bitstreams/4e04db10-e21e-4aea-b772-4a50bf1585aa/content). [Szczecinski 2017, Eq.18](https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2017.00037/full).

Selected wording: **Syn6 is an engineered NMF/half-center hybrid**, not a literal Di Russo circuit. An Eq.18 equilibrium transform is not proof of complete dynamic waveform fitting, full Szczecinski 2020 Table2 design, or biological identification of six spinal synergies. K8 improves inverse reconstruction on the same dataset but has not replaced production K6 or shown matched forward-controller/body improvement. See `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\SYN6_K6_K8_DESIGN_REVIEW.md`. Do not interpret independently learned left/right component indices as homologues; the right-canonical mirror is an explicitly new design policy.

## 5. Which new designs are actually composed?

| Design | Verified composition | What remains separate |
| --- | --- | --- |
| Graded s3k / w2lvar / Syn6 electrical candidates | Graded terminal maps, qualifying **existing** central budgets, eligible IBEXC coincidence, final sparse caches, and dense/sparse neural checks are composed. | Saved standing masks are retained. Active expression fixtures are different conditions. No historical body result transfers. |
| Six dissertation recipe clones | Terminal maps, input-domain bounds, and final cache identity are statically checked for the reconstructed full/pruned recipes. | Central budgets, heel A/B, canonical mirror, new monoarticular overlay, and coincidence were not silently added to these historical reconstructions. |
| New monoarticular central overlay | Final graph export, individual current operators, domain guards, and fresh repeat exist for separate s3k/Syn6 full-template proposals. | Not yet the same whole-controller composition as the previous row; no phase/walking inheritance. |
| New Syn6 single-IN repair | Isolated full-recipe clone: 794 cells; two heel→InF conductances become zero (3,037→3,035 pairs). The subsequent explicit canonical-right/graded/current/coincidence composition has 3,052 pairs, 368 passing terminal contracts, exact cache/backend parity, and byte-identical fresh graph/proof. | This named profile does not include the separate 22-muscle central overlay or overwrite old FSA exports. No body, training, or complete phase-response result is inherited. |
| Current integrated Hunt1033 | Actual local-reflex correction, canonical-right HUNT_AF timing budget, IBEXC coincidence, saved paired gain checkpoint, and authoritative command-update/cache checks are composed in `integrated_hunt_design.py`. | PT is retained but its body calibration is explicitly separate from unbound electrical proof. Direct contact reset is retained; heel B is not applied. Actual body pilot results belong to the campaign record. |
| Current central-spiking Gait candidate | Canonical inner graph, graded terminal map, and frozen v4 ionic/event-readout core have a composed electrical fixture. | Optional sensory-spiking conversion is a separate component; do not silently stack it. It is not a fully spiking or literal all-GLIF design. |

Current Hunt profile: `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\integration\PROFILE.json`. Its source is the saved `aligned_support_ladder_v1` condition and its explicit `GAIN_PARAMETERS.npz`, rather than an arbitrary mixture of newer winners. The timing mode is `canonical_right_role_budget`, with separate E/F 0.5 nA budgets. These budgets differ from the local signed **0.5 mV** equilibrium contracts; nA and mV are not interchangeable.

The saved Hunt controller cleared MN conductances and restored an incomplete reflex set on every command update, leaving 92 homonymous II→IIX→MN outputs zero. The opt-in local correction writes the authoritative runtime buffers and retains those 92 paths after updates. This is a real implementation repair, not evidence that the historical model had no length-sensitive HUNT_AF inputs. The four RF/HAM cells per leg are two integration outputs plus two inhibitory relays; sartorius, gracilis, and TFL remain separate. PT is an engineering task-balance pathway, not merely a renamed myostatic reflex.

HUNT_AF and the phase-bank cells are implementation feature/basis names, not identified biological neuron classes. The four-part training-method source is Hunt's professional-meeting lecture: A. J. Hunt (2019, February 15), *Real-time control of dynamic walking in a dog-like robot*, 2nd Intel Neuromorphic Research Community Meeting. The implementation's functional basis, PT recruitment, and phase-bank fitting remain engineering realizations; a four-part method citation does not establish literal reproduction of a published connectome or completion of the physical five-stage curriculum.

## 6. Other families: scoped disposition

| Family | Selected interpretation and limit |
| --- | --- |
| W2L108 graded / 704 hybrid-spiking | Existing 216-edge functional graph, 12 homonymous II channels plus four retained phase-II paths, physical contact bindings, and exact source identities remain auditable. A mislabeled right knee Renshaw cell is a logical alias, not missing circuitry. Graded and event/ionic implementations need their own actual clocks and current operators. Supported performance is configuration-specific. |
| Public rat120 | Authored precursor with 16 muscles; preserve source asymmetries and its recorded right-side latch. It is not an exact released final-2023 implementation. |
| New rat122 | Explicit bilateral biarticular integration reconstruction with 216 mirrored edges, preserving and mirroring actual left integration architecture. Its Renshaw current budget is an engineering choice; electrical proof is not a bilateral air/body pass. |
| Curated rat134 | Audit-only, excluded from training: unilateral prototype additions and six mislabeled adapters actually reading Ia prevent blanket biological/pathway approval. Preserve its earlier results as diagnostics. |
| CAD262 | Original ideal-tension controller is hybrid. The later two-pressure-BPA pair retains 52 hybrid routes, not a fully pneumatic SNS walker. Actual IaIN/IIIN/RC populations are absent; dormant builder options are not implemented subsystems. Mechanical repair remains a prerequisite for new stance claims. |
| SNS-Simscape graded / spiking reflexes | Actual RC/graded chemical ports and fixed-reset LIF/additive-capped event synapses are tested under their own functions. Do not relabel the latter adaptive-threshold GLIF or reset-to-peak Eq.27/29. |
| Reduced Deng14 / Rybak14 | Distinct neural mechanism demonstrations; preserve each source's ionic dynamics and unit convention. No body result is required for their limited neural scope, and no walker result is inherited. |

Author-choice/evidence index: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\NON_GAIT_AUTHOR_CHOICE_RECORD.md`. Deng source distinctions and the published 2019 ankle-phase mismatch are documented in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\LITERATURE_MODEL_SOURCE_AUDIT.md`: copying an EXT label is not anatomical moment-arm or stance-phase validation.

## Compact review checklist

- [x] Historical full/pruned recipes separated from later standing snapshots.
- [x] Heel path laterality, two-hop sign, duplicate-pair overwrite, and pulse-versus-level semantics resolved from actual code.
- [x] Central versus local Ia/II/Ib routes distinguished; source omissions and new engineering rules labelled separately.
- [x] Syn6 NMF channels, Di Russo analytic primitives, and inverse versus forward evidence separated.
- [x] Composed electrical profiles distinguished from isolated/static additions; authoritative Hunt runtime II repair identified.
- [x] All-family scope, retained ionic/spiking differences, source preservation, and rat134 exclusion stated without a blanket pass.
- [ ] Actual whole-controller body-domain measurements, hierarchy waveforms, and acceptance decisions: use each newly authorized pilot's own results.
- [ ] New central/heel policies beyond the selected integrated condition: require an explicit combined graph/current/phase test before adoption.
- [ ] Body-calibrated PT, biarticular/non-sagittal torque/phase coverage, and task-specific holdouts remain separate from isolated equilibrium checks.

Historical builders, saved winners, source models, and manuscripts remain unchanged. The separately authorized Syn6 repair changes exactly two conductances in an explicit new compiled candidate; its 128-step electrical proof and fresh repeat do not advance body physics or train parameters.
