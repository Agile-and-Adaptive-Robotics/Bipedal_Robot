# Session 1: neural-controller handoff for the defense deck

**Machine: easteregg2. Repository: `D:\GitHub\Bipedal_Robot`.** This handoff is for Ben's defense-deck session on this same machine. Ben directs the content. Current work is controller design and verification; training and curriculum execution await his explicit OK. Compliant tendons, joint damping changes, and muscle damping changes are separate optional experiments, each requiring his OK. No manuscript or Overleaf edits are authorized by this handoff.

**Latest review, November 2:** read `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\GAIT_DISSERTATION_RECIPE_REVIEW.md` for the six actual dissertation-recipe reconstructions and their passing static FSA checks. Circuit repairs remain separate candidates; they have not all been composed and dynamically validated as one runnable controller. The updated DEFENSE_DESIGN_DECISIONS.md explains why Di Russo is an architecture-specific comparison, not a universal wiring specification. Full Hunt/curriculum training is still on hold.

## Read these first, in order

1. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DISSERTATION_TRIAL_PROVENANCE.md` — essential correction: dissertation trials differ from the later standing exports. Start here before describing missing pathways or drawing historical controllers.
2. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DISSERTATION_FIRST_PASS_FIX.md` — passed numerical checks and remaining design/wiring issues. Its historical-trial correction takes precedence over earlier standing-model generalizations.
3. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DISSERTATION_MODEL_MAP.md` — submitted model names versus later research variants. The submitted public rat has 120 neurons; Rat134 is a later diagnostic prototype.
4. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\ROUTE_BUILDER_EXPORT_CHECK.md` — exact EB475WS4 drawing files, hashes, actual routes, and correspondence with easteregg2.
5. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DEFENSE_DESIGN_DECISIONS.md` — committee-facing design answers and corrections; use the provenance correction above when discussing historical trials.
6. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\EXPERIMENT_RECORD.md` — separate, continuing experiment record. `CHATGPT_REPORT.md` is a shorter inter-session handoff.

## Critical distinctions for slides and committee answers

- **Historical full s3k:** Appendix D's −186.4 baseline/ablation result traces to `prune_matrix.py`, v10 base multipliers plus `reports_20260923\s3k_trial34_full_params.json`, and `_curriculum.set_stage(5)`. Central Ia, flexor II, extensor II, and extensor Ib gains are **0.39859, 0.38390, 0.09300, and 0.24355**. `full_rules=1`, IaIN is enabled, and Renshaw is set to 0.5. It is incorrect to say these dissertation trials lacked central feedback merely because the later standing export does.
- **Historical syn6:** the −200.4 result traces to `curriculum_syn6_stage5.json`; heel/toe gains are nonzero, pooled central Ib gain is zero, and the general central Ia/II implementation is absent. The heel-to-both-ipsilateral-RG-interneurons concern therefore matters for this recipe too. Trace signed chains: InE inhibits RG-F; InF inhibits RG-E.
- **Later standing exports:** s3k410, syn6, and w2lvar JSONs under `Gait2392_standing_families` describe named later standing conditions. They are not universal diagrams of every dissertation run. s3k/w2lvar standing settings suppress PF-to-MN expression. Historical source reconstruction with current code must be labeled as reconstruction, rather than a byte-identical saved historical graph.
- **What the new checks establish:** local functional gains, units, clocks, physical input domains, actual signed currents, parameter/cache consumption, symmetry policies, and accelerated/reference agreement. These passes do not establish that every intended physiological pathway exists or that a whole curriculum has passed. The complete design gate remains closed while wiring issues are resolved.
- **Rat models:** original public precursor120 remains preserved. New bilateral-biarticular122 retains the left integration subsystem and adds its actual right counterparts; electrical checks pass. Rat134 contains unilateral diagnostic additions with mismatched sensory labels and disabled physical outputs; it is not adopted as a trained literature controller.

## Zero conductance: describe implementation evidence carefully

For a finite graded conductance synapse, `I_syn = g * activation * (E_rev - V_post)`: an actually zero effective conductance contributes zero synaptic current. Ben has observed response differences and recalls a ZCode rule that removes zero-conductance connections. Its use across all historical implementations has **not** been established.

Current code provides two concrete details:

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\exact_fast_eval.py` constructs sparse indices using `np.nonzero(compiled.g_max_non)` and copies their conductances into `self.gmax`. Edges already zero when that backend is constructed are omitted.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\sparse_pairwise_eval.py` steps using that cached `gmax`. Changing only the dense matrix after construction can leave the cached value unchanged. This is a possible implementation explanation, not an established diagnosis of Ben's earlier observation.

Distinguish zero sensor input, zero effective conductance, a removed edge, a stale cache, and a topology rebuilt with different state/indexing. A gain omitted at construction may require rebuilding to enable the pathway. Spiking ports may also have separate conductance state. Use actual-current and disabled-path controls for the selected runtime; do not claim either universal equivalence or an intrinsic effect of a truly zero graded conductance.

## Design, scaling, and schematic references

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\MODEL_DIMENSION_AND_CONTROL.md` — six floating-base coordinates plus 17 internal coordinates; reduced-controller transfer, control coverage, and apparatus contribution. Its torque-rank diagnostic uses raw unpatched CVT3 mechanics.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\GAIT_TOPOLOGY_SOURCE_REVIEW.md` — actual default/standing/full-rules construction, explicit inhibitory populations, heel laterality, and conditional central routes. Read alongside the historical-trial correction.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\GRADED_GAIT_DESIGN.md` and `GAIT_NEURAL_COMPOSITION_REVIEW.md` — electrical FSA and scaling evidence; stated phase budget applies to encoders, not the complete neuromechanical loop.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\SPIKING_FAMILY_REVIEW.md` — separate legacy spiking comparison, new central-spiking design, and optional six-pool sensory conversion. Mixed versions must be identified accurately.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DESCENDING_POSTURE_AUDIT.md` and `POSTURAL_COMMAND_DESIGN.md` — actual standing POST inputs, CAD's hybrid allocator, and Hunt/PT task feedback versus the proposed descending antagonist posture organization.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\HEEL_STANCE_POLICY_REVIEW.md` — two separate optional heel-routing hypotheses with static signed-current proofs. They are not adopted historical corrections or demonstrated dynamic resets.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\SYN6_K6_K8_DESIGN_REVIEW.md` and `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_results\rank8\fsa_rank8_report.md` — matched K6/K8 inverse reconstruction, not a forward-walking comparison. More NMF factors do not establish more biological synergies or completed RG-to-PF generation.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\FUNCTIONAL_SUBNETWORK_OTHER_MACHINE_PROMPT.md` — editable schematic requirements: actual model/configuration, named synapses, readable functional subnetworks, and explicit source-versus-proposed wiring.

## Existing review media and measured support

These are saved research demonstrations. Keep each GIF attached to its own configuration and result; later FSA profiles have not been body-tested in these recordings.

| Media | Correct short description |
| --- | --- |
| `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\Gait2392_Hunt_supported_stepping.gif` | Supported stepping in place; 5 right/4 left sustained mesh-cleared swings; measured upward support 63.73% BW. Protected knees. Same controller's earlier standing stage remains unresolved. |
| `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\Gait2392_central_spiking_standing.gif` | Central-spiking Gait2392 supported standing; mean upward support about 0.5944% BW after the separate mass-denominator correction. Lateral/rotational guidance remains. |
| `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\cad\CAD_BPA_PAIR_SUPPORT085.gif` | CAD robot with two neural pressure-controlled BPAs plus 52 ideal-tension hybrid routes; measured upward support 87.337% BW. |

Full media provenance and still-image paths:

- `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\GAIT2392_REVIEW.md`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\GAIT2392_REVIEW_MANIFEST.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\cad\README.md`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\cad\REPRODUCTION_MANIFEST.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\GAIT_V4_SUPPORT_MASS_CORRECTION.json` — correction of the fixed-ground contribution to the standing BW denominator; raw forces/trajectories unchanged. Use this over older rounded percentages in media notes.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\GUIDED_WALKER_MEDIA.md` — W2L guided-walker media and context.

## Manuscript snapshot and presentation constraints

Ben's supplied PDF exists at `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf` (18,315,312 bytes; file timestamp October 31, 12:20 p.m.). Source saved during the submission pass is `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\session1_submission_20261031\source`; full downloaded archive is its sibling `Overleaf_source.zip`. Do not extract these over the working manuscript.

Use concise demonstrations: a few cycles aligned across RG, PF, MN, and joint kinematics; readable functional subnetworks; actual support forces; and explicit stimulus timing. Separate recorded results, electrical design checks, and proposed wiring. Emphasize the implemented toolchain, scoped measured acceleration, transfer methodology, and the concrete corrections. Do not invent unmatched engine timings or present a neural-only speed-up as a complete-loop speed-up.

Scope-specific raw records live under `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design` and `fsa_design_20261101`. No new training, curriculum, tendon, or damping experiment should be launched by the deck builder.
