# Connectome and optimization audit

## Scope and evidence

Read-only inspection of the two local source PDFs, the embedded v3.3 template library, the sidecar template library, current builders, and W2L project files. No controller, project, optimizer, or dissertation source changed. Inventory results are in `inventory.json`. Paper equations and tables were extracted directly from the PDFs; Di Russo Figure 4 and Shinohara Table A.1 were also rendered and visually inspected. Prior campaign outcomes below are reported evidence, not newly reproduced results.

Ben selected separate Di Russo and Shinohara models. They must remain separate from the existing mixed-literature production models.

## Findings

1. **The corrected W2L is missing from the editor's embedded library.** `connectome_templates.json` contains `bilateralrg2`; the HTML's embedded library does not. Every shared entry matches exactly. The editor loads the embedded library first, so normal use presents the earlier `bilateralrg` topology, not the corrected October 18 topology. `make_editor_templates.py`'s menu mapping also has no `bilateralrg2` entry. Adding that entry and synchronizing both libraries is required before the editor can accurately represent the latest model.
2. **There is no Di Russo template in the embedded library.** `replication/dirusso_connectome_draft.json` is a specification, not a runnable SNS controller or an editor node/edge network. Four phase cells and six NMF synergies are neither the paper's five primitives nor an exact equivalent.
3. **The Shinohara editor entry is Ben's corrected export, not an exact bilateral executable replication.** `make_editor_templates.py` loads `ben_shinohara_20260924.json`, with two symbolic afferent nodes expanded for inspection. Its parent canvas has 20 nodes and 37 edges, seven MNs, and incomplete opposite-side context. The two afferent subviews list constituent sensors with no internal edges. The generator explains the autogenic fan-out in prose; that explanation is not a compiled per-muscle wiring contract.
4. **The W2L project has one neural page.** The main, control20261020, and tuned20261020 `.aproj` files each contain one Subsystem named `Neural Subsystem`, with one valid DiagramXml page. Their whole XML and embedded page XML parse successfully. This confirms the organizational work Ben requested is still needed. It does not establish successful GUI opening or locomotion.
5. **The initial LIT_CIRCUIT_AUDIT verdict table is historical.** It describes an earlier representative build. Later sections, current builders, the October 20 audit, and current runtime configurations must be distinguished. For example, current `build_network.py` contains IIX, IIIN, and IBIN branches behind `full_rules`, contradicting a reading of the older audit as a current absence claim. The production anchor is explicitly a `full_rules=0` configuration. Supported code, enabled topology, and optimized winner are three different things.
6. **The current base builder does not exactly implement Di Russo.** Renshaw cross-inhibition spans other same-side RC pools, not only antagonist RCs. IaINs receive PF phase drive; the Ib-positive branch receives RG-E drive; runner sensory gains are phase/speed gated. The paper deliberately excludes reflex state-machine gating. Antagonists are classified from primary muscle groups, whereas Di Russo Table 1 explicitly includes biarticular roles. The initial IBIN/IaIN wiring checks whether the target was already created, but `_complete_in_mutual` subsequently fills the missing directed edges when full_rules is enabled. This is an existing correction, not an outstanding one-way-wiring finding. These are static findings and a test requirement, not newly measured behavioral effects.
7. **The old conflict audit overstates genuine conflicts.** Its label normalization collapses sensory Ia and IaIN, own and antagonist MNs, ipsilateral and contralateral nodes, and explicit versus intervening inhibitory pathways. Its eight sign-conflict rows are candidate collisions, not eight demonstrated biological contradictions. Gains from different voltage scales, conductance units, and population normalizations cannot be compared by a raw greater-than-two ratio.

## Template interpretation

| Editor entry or family | What it represents | Compliance conclusion |
|---|---|---|
| `rules` | Ben's expanded rule panel | Rule reference; not a complete bilateral walking controller, and not independently certified against every audit row. |
| `rules_motifs` | Older automatically generated motifs | Partial rule examples. |
| `biped2xcpg` | Connectivity audit mined from the Deng-derived AnimatLab model | Strongest organizational precedent for subnetworks, but generator documentation explicitly flags inferred signs. Not an exact Di Russo replication. |
| `w2laproj`, `bilateralrg`, `li` | Different AnimatLab controllers or earlier mines | Source-specific topology. `bilateralrg` predates the corrected W2L. Missing rule families are not evidence that the source paper was implemented incorrectly. |
| `walker_v9`, `walker_v10`, `walker_s1/s2/s3` | Historical SNS snapshots | Not current production proof and not exact paper replications. |
| `walker_s3k`, `_flat` | Compact and flat production schematics | Four PF windows and mixed sources; not Di Russo or exact Shinohara. Representative afferents and RCs are not a full census. |
| `walker_s3k_pruned`, `_flat` | Explicitly pruned variant | Cannot claim compliance with removed contact, interleg, and Ib pathways. Walking performance does not restore those mechanisms. |
| `w2lvar`, `_flat` | W2L-style PF organization on the 92-actuator human model | Extensive reflex interneurons and RCs. It is distinct from the physical W2L model and remains a mixed-literature architecture. |
| `syn6`, `_flat`, `synergy6` | Six-channel synergy organization | Not the five equally spaced Di Russo primitives; channel count alone is insufficient evidence of equivalence. |
| `rybak` | Rybak-derived rule sketch | Assess against the specified Rybak model and documented reductions, not every reflex in Di Russo. |
| `shevtsova` | Ben's laminar RG export | Rhythm/interlimb motif, not a complete sensory-motor reflex controller. Separate SNS template verifications exist with documented reductions and calibrated constants. |
| `shinohara` | Ben's edited Shinohara schematic | Relevant target reference; incomplete bilateral executable specification. |
| `spiking_mirror` | Hybrid spiking reinterpretation | A different neuron model; topology resemblance does not establish dynamics equivalence. |
| sidecar-only `bilateralrg2` | Corrected W2L transcription | Current W2L runtime path, absent from embedded editor. Source-sign and conductance-mapping deviations remain documented. |

No entry can presently be certified as satisfying the entire audit and exactly reproducing either requested paper. That is not a claim that every entry is biologically invalid: the audit combines distinct species, tasks, and levels of abstraction.

## Real source conflicts

### Di Russo Renshaw polarity

The Methods on PDF page 6 state mutual inhibition between antagonist RCs. Figure 4(e)'s caption on PDF page 8 and Table 1's caption on PDF page 9 describe reciprocal excitation. The repository's audit and draft choose mutual inhibition. Retain that recorded interpretation for the separate implementation, label the source inconsistency, and do not characterize the PDF as unambiguous. The same issue affects any assertion of exact source replication.

### Shinohara RG-E to PF target

Visual inspection of PDF page 31 confirms Table A.1 prints RG-E -> PF-F at 0.7, with no value under PF-E. The prose describes corresponding RG/PF excitation, and the repository's draft chooses RG-E -> PF-E. This is a real table/text conflict, not an extraction artifact. Ben ruled: use RG-E -> PF-E and flag the printed-table conflict. Preserve both readings in provenance and test the alternative as a sensitivity condition, not as a hidden optimizer topology switch.

### Rules with different scope

Di Russo's signed primitive weights and ungated reflexes, Shinohara's excitatory central afferents and voltage-dependent sodium inactivation, and Deng-style phase-gated pathways should not be silently merged. A fixed tau_h workaround for the stock SNS backend is not permission to replace Shinohara's published tau_h(V) = 320 + 320/cosh((V+35)/15) ms. A bell-tau SNS backend already exists in `replication/shevtsova_sns.py` and provides a reusable mechanism, with its update-order limitation explicitly documented.

## Optimization assessment

**Parts are appropriate, but the current record does not justify an overall validation claim.** Staged tuning, fresh-process winner replay, per-term diagnostics, ablations, and local basin probes are appropriate. Optimizer choice is secondary to a correct task, implementation, objective, and held-out validation.

- Di Russo section 2.3 first fits neural outputs to reference simulated excitations with recorded state trajectories. It then optimizes dynamically consistent stability and physiology, and finally metabolic performance. Kinematic curve fitting alone is not that procedure. The paper explicitly says imitation alone does not yield stable gait.
- Shinohara section 2.3 optimizes seven PF-to-MN weights and six feedback parameters on step-to-step belt levels {-2, 0, +2} cm, and tests sudden holes after training. Its Appendix C also re-optimizes feedback-ablated models. Both frozen-parameter ablations and re-optimized ablations answer useful, different questions. A good flat-ground score does not validate the ground-loss mechanism.
- The current corrected-W2L campaign fails its own strict basin bar: the report gives at best 9/12 versus required 11/12. Report that as marginal phase robustness, not robust walking.
- The October 20 drive-wean ladder is invalid: `knobs_of` passes free `te` and `tf` while `dscale` is only bookkeeping in this campaign. Its report explicitly withdraws the weaning claim. A real ladder must constrain effective drive and prevent the optimizer from compensating through another drive parameter. This does not retroactively invalidate a different, properly scaled campaign.
- The tuned AnimatLab transfer report records a fall and latched RG, not successful walking. ContactCount-based sensing cannot implement the MuJoCo force-normalization knob `cref`; tonic drives and conductances require simulator-specific calibration. Matching numbers or wiring does not establish matching dynamics.
- Preserve raw height, tilt, support forces, contact events, forward progress, both limb cycles, duty, and phase error. Post-fall bursts, frozen runs, in-place stepping, and harness-supported progression need distinct labels. No harness-free walking is established by this inspection.
- Finalist validation should include multiple initial states and perturbation seeds, held-out surfaces and speeds, longer horizons, timestep convergence, and effective-parameter sensitivity. Tune physiology parameters within justified bounds; do not let optimization silently delete a required paper mechanism.

## Implementation order and subsystem contract

1. Freeze existing production code and its recorded gates; use separate `dirusso` and `shinohara` selectors and complete parameter files.
2. Build one flat, typed, per-muscle edge manifest per model with side, muscle, target muscle, receptor versus IN identity, sign, source, equation/table, units, delay, and interpretation status. Validate required AND forbidden edges.
3. Generate matching diagram hierarchies from that manifest: left/right RG, commissural coordination, PF, hip/knee/ankle sensory-motor modules, and per-muscle Ia/II/Ib/RC motifs. Keep all cross-boundary endpoints explicit. Never infer executable neuron connections from a collapsed box-to-box wire.
4. For W2L, use the corrected main `.aproj` as the structural source and retain tuned/control files as explicitly different parameter variants. Add GUI subsystems in a new reviewable copy; preserve GUIDs, types, gains, stimuli, and the runnable `.asim`. Validate every page and compare the flat circuit before and after reorganization. Then check GUI loading and replay behavior.
5. Di Russo's phase oscillator and five primitives belong outside the reflex leaky-integrator network as in the paper. Standard conductance-based SNS synapses are not numerically identical to its signed weighted-sum leaky integrators. Any SNS implementation must preserve those equations or clearly label the approximation. Preserve receptor delays, biarticular classifications, and the separately summed balance signal.
6. Shinohara needs its own RG/PF/MN sodium currents, published bell tau_h, per-muscle encoding, and central targets. Its 2D cat geometry, forelimb support, and treadmill are not the human or W2L robot model; deploying the controller on either is a documented adaptation, not an exact neuromechanical replication.
7. Run neural unit/dynamic checks, coupled model checks, then the paper-specific staged optimization. Do not launch a large study until hard validity gates and objective-ranking tests pass.

## Per-directive status

- Locate latest W2L project: confirmed folder and main/tuned/control variants.
- Separate Di Russo and Shinohara: recorded Ben's choice; no hybrid introduced.
- Audit editor against literature: inventory and source-level findings complete; no full dynamic compliance certificate claimed.
- Resolve conflicts: real source inconsistencies identified; Shinohara target resolved by Ben as PF-E, with the table conflict retained.
- Add runnable paper models: not implemented in this audit.
- Reorganize AnimatLab and matching MuJoCo subnetworks: not implemented in this audit.
- Assess optimization: concrete invalid-wean, robustness, transfer, and task-validation limitations identified.
- Verification: 24 embedded entries, no top-level duplicate labels or dangling endpoints; all shared sidecar entries equal; corrected sidecar entry absent from embed; three W2L projects and their three DiagramXml pages parse. No simulation rerun or dissertation change claimed.

## Scope clarification

Ben requested a big-picture assessment of obvious mis-wiring and credit restraint before implementation. No new controllers, large optimization runs, or model reorganizations were started. Separate models remain the requested architecture.


