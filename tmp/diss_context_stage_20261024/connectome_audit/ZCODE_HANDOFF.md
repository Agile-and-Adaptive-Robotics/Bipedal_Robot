# ZCode handoff: connectomes, W2L, and syn6 PF

## Scope

Ben requested documentation and stopped further investigation for now. Do not treat this handoff as authorization to launch optimization, restructure models, or edit dissertation sources. No models, controllers, optimizer code, or dissertation files were changed in this session. `AUDIT.md` and `inventory.json` in this folder contain the detailed inspection results. The audit is source-level inspection, not a new simulation verification.

## Ben's decisions

- Implement Di Russo 2023 and Shinohara 2025 as SEPARATE models, not a silent hybrid.
- Shinohara: use RG-E -> PF-E at 0.7 and flag the printed Table A.1 conflict. Direct visual inspection of PDF page 31 confirms the printed 0.7 is under PF-F; the prose and existing draft favor PF-E. This is a real source inconsistency.
- Recently updated W2L lives in `Neuromechanical_Models/Walker_2_Layer_CPG_BilateralRG/`. Distinguish the corrected main `.aproj` from the `control20261020` and `tuned20261020` variants.
- Ben wants more neural subnetworks like the Deng AnimatLab controller, with matching understandable organization in MuJoCo-SNS and reusable subsystem figures. This organizational work has NOT been done.

## Priority unresolved issue: syn6 timing

Ben challenged the explanation that the six synergies can be organized mainly as E/F-driven channels with staggered time constants. He observed that the primitives used to find the synergies appeared to include antagonist muscle groups active at DIFFERENT parts of the gait. Preserve this concern. Do not present the current family assignment as validated physiology or infer that Ben endorsed the architecture.

The explanation given to Ben described the CURRENT IMPLEMENTATION, not proof that it preserves the source primitives. Six synergy channels are not necessarily six independently reproduced temporal primitives. Conversely, common RG input does not alone prove that their final waveforms are identical. This requires measured waveform comparisons.

### What was directly checked

Source: `spinal/build_network_syn6.py`, especially `_load_synergy_tables`, `_build_pf_s6`, and `_wire_muscle`.

- Per leg: six PF output populations `PF_S1` through `PF_S6`, plus `PF_IN_E` and `PF_IN_F` inhibitory relays. Twelve PF output cells bilaterally; sixteen including those four PF relays. The stock four-output-cell PF is a DIFFERENT architecture.
- The loaded W matrices have shape 43 muscles x 6 channels per leg. Positive entries map to excitatory PF-to-MN conductances through the Eq-18 mapping. These are spatial muscle weights, not the complete temporal reconstruction.
- Current family assignment on BOTH sides: S1 F, S2 E, S3 E, S4 F, S5 F, S6 E.
- Family assignment uses the bilateral mean of phase-profile stance fractions when the FSA profiles are available. Current fractions S1..S6: 0.195713, 0.952928, 0.922159, 0.246330, 0.410065, 0.889709. These are the builder's profile-derived ratios, not percentages of contact time or proof of biological stance membership.
- S5 falls in the mixed band (0.35, 0.65): receives both RG-E and RG-F drive, weighted by that fraction. Other channels receive their single-family RG drive.
- S5 is excluded from reciprocal PF family inhibition by default, as a documented latch workaround. `pf_mixed_recip > 0` restores its participation. This exception needs to remain visible in diagrams and interpretation.
- PF response-time multipliers S1..S6: right [1.6, 1.2, 0.9, 1.25, 0.9, 0.6]; left [1.6, 1.2, 0.6, 1.25, 0.9, 0.9]. These multiply `TAU['pf']`; they are not independent phase oscillators or prescribed onset times.
- The inspected stage-5 JSON tunes scalar RG/PF and feedback parameters; it does not contain a free per-muscle six-channel weight matrix.
- Default params for `pf_mixed_recip`, `iain_pf_gate`, and `ib_e_central` were all zero in the inspected environment. This is NOT a claim about every later winner or command-line configuration; use its authoritative loader and record effective settings.

### Investigation to do if Ben resumes

1. Trace exactly which source signals were decomposed: muscle excitations, activations, back-solved neural drive, or another quantity. Read `synergy_model.py`, `fsa_backsolve.py`, and their recorded generation settings. Do not conflate spatial W with temporal coefficients H or measured neural waveforms.
2. Compare the original six temporal coefficients/primitives with the actual six PF voltages and resulting MN activations, on each side, using the SAME gait-phase convention. Determine whether cross-side channel indices were actually matched before bilateral averaging; NMF components have permutation and scaling ambiguities.
3. Show antagonist muscle membership and its source phase, including biarticular roles. Quantify onset, peak, duration, within-family correlation, and reconstruction error. Do not reduce each channel to its stance-average label before examining the source waveform.
4. Test whether RG family drive plus first-order filtering/reciprocal inhibition retains the required early/late stance and swing timing. Distinguish temporal phase offsets from changes in response time. Check S5 separately and measure the effect of its missing default reciprocal inhibition.
5. Only after that evidence, propose whether to retain, reassign, or redesign the PF architecture. Ben has NOT authorized a specific replacement. Di Russo's five primitives are a separate paper model, not an automatic replacement for syn6.

## Wiring discrepancies worth addressing first

- Corrected W2L `bilateralrg2` is in `connectome_templates.json`, but absent from the editor HTML's embedded library and menu mapping. Shared entries match. The editor prefers the embedded library and can show the older W2L.
- AnimatLab `R Hip MN flx RE -> R Hip MN ext RE` carries an excitatory type, while SNS transcribes inhibition under Ben's earlier ruling. Report the discrepancy; do not silently pick a new sign.
- Corrected AnimatLab crossed V3 outputs use the shared 2.749 uS class versus the earlier dedicated 0.1 uS class. The October 20 transfer report flags strong crossed excitation as a latching suspect. This is documented evidence, not a fresh causal test.
- Base `build_network.py` Renshaw cross-inhibition spans other same-side RC pools, not just antagonist pairs, unlike the Di Russo interpretation. This statement concerns the base builder; syn6 explicitly has no RC populations in its documented architecture.
- Do NOT repeat a claim that IaIN/IBIN mutual inhibition is simply one-way: the base builder has `_complete_in_mutual` to fill missing directed edges when full_rules is enabled. The audit was corrected after finding this pass.
- Di Russo's Methods say antagonist RC mutual inhibition, while captions describe reciprocal excitation. Existing repo interpretation follows the Methods. Keep the inconsistency flagged.

## Optimization cautions already supported by the record

- October 20 corrected-W2L drive ladder did not constrain effective drive: te/tf stayed freely optimized while dscale was bookkeeping. Its own report withdraws the weaning claim. Do not generalize this defect to separate campaigns without checking their implementation.
- Corrected W2L's reported best basin was 9/12 against its required 11/12, with phase-lock failures. No new basin trial was run here.
- The tuned MuJoCo-to-AnimatLab transfer report records latching and falling. ContactCount versus force-normalized sensing, neuron dynamics, and effective conductances differ; numeric gains do not transfer unchanged.
- Di Russo's imitation -> dynamic stability -> metabolic stages and Shinohara's randomized-surface training plus held-out ground-loss evaluation are distinct tasks. Kinematic fit or a good supported stepping score alone validates neither paper's full mechanism.

## Files left by this session

- This handoff.
- `AUDIT.md`: detailed findings, template interpretation, source conflicts, and implementation/validation recommendations.
- `inventory.json`: 24 embedded template inventories and three W2L project/page parse results. Only top-level template duplicate labels and dangling endpoints were checked; do not claim recursive or behavioral verification.
- Repo `tmp/connectome_audit/`: scratch text extraction of the two PDFs and rendered source pages. These are intermediate inspection files, not model changes.

No commits or pushes. Ben commits through GitHub Desktop. No production baselines were rerun because production code was unchanged and Ben requested a bounded inspection.
