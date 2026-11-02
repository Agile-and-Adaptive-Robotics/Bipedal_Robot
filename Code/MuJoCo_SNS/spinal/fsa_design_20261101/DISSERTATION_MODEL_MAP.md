# Dissertation models and later research variants

Session 1, easteregg2, 2026-11-01. Read-only check of the archived submitted chapter sources; no manuscript changes.

Submitted source reference: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\session1_submission_20261031\source\60-sns.tex` and `95-AppendixD.tex`. These archival copies identify the model names below. They are distinct from the mutable local ProofFinal working tree.

| Submission description | Current audit treatment |
| --- | --- |
| Original Li–Deng/W2L shared-RG AnimatLab walker | Historical source and recorded suspended/fall tests; preserve its original organization. |
| Newer bilateral W2L with contact and flexor/extensor type-II feedback | Prioritize the actual W2L108 source and its recorded supported tests. |
| Gait2392 full controller, W2L organization, and six-synergy organization | Prioritize s3k, w2lvar, and syn6. Each dissertation result has its own settings; later standing conditions do not replace those results. |
| Graded and spiking reflex/SNS-Simscape demonstrations | Prioritize actual neuron and synapse block implementations and their own unit conventions. |
| Reduced Rybak two-level circuit | Neural mechanism example with 13 recorded checks, Appendix D; add this actual source to the design registry. |
| Reduced Deng-derived fourteen-neuron circuit | Neural mechanism/stimulus reference, distinct from the public rat body. |
| Public biarticular rat precursor | **120 neurons, sixteen muscles, fixed pelvis**, stated in chapter 6 lines 305–306. Audit this original source first. |
| CAD-derived robot, including a 262-neuron controller | Preserve the submitted ideal-reference-tension/hybrid example. Later two-BPA support085 results are separate research conditions. |

**Rat134 is a later curated receptor-geometry variant, named for its 134-neuron count.** It is not the 120-neuron public precursor described in the submitted chapter. The earlier AUDIT_ANCHORS.md selected current research reference conditions, not submission versions. Likewise central-spiking dynamic_readout_v4 and its standing checkpoint are later research conditions, not a replacement for the dissertation's spiking comparisons.

Priority: finish the appropriate design, physical-time/signal scaling, and pathway checks for the submitted implementations first; then finish the other active variants. All training and curriculum execution await Ben's explicit OK. Optional tendon compliance, joint damping changes, and muscle damping changes each remain separate, unapproved options.
