# CHANGES MANIFEST — AnimatLab realizer stage (2026-10-16, egg)

Campaign stage: phase 9 robot-track realization. Repo
`D:\GitHub\Bipedal_Robot`, branch KneeTestSetup_BenBo_stw, HEAD b497b67b.
NO commits, NO pushes.

## Working tree

**No tracked file was modified by this stage.** Verified post-run:
`git -C D:\GitHub\Bipedal_Robot status --porcelain Neuromechanical_Models`
→ empty. Ben's `Walker_2_Layer_CPG_BilateralRG` folder (and every other model
folder) is untouched; the V0 control asim in the reports folder is
byte-identical to the source (`fc /b` passed). The only repo-tree additions
are the untracked campaign folder below (plus what earlier stages created).

## Created (all untracked, campaign-scoped)

`Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab\`
- `W2L_ANIMATLAB_REALIZATION_20261016.md` — stage report
- `CHANGES_MANIFEST_20261016_w2l_animatlab.md` — this file
- `build_winner_variant_20261016.py` — variant builder (d51a GUID prefix;
  built-in XML + page-CDATA verification)
- `analyze_run_20261016.py`, `probe_crossings.py` — chart analyzers
- `edits_log.json` — machine-readable log of every edit
- `V0\` control asim + chart byproducts + analysis.json
- `V1\` drive-only probe asim + byproducts + analysis.json
- `V2\` winner-set asim + **aproj (HEADLESS-RUN ONLY, GUI open-check pending
  Ben)** + byproducts + analysis.json
- `V3\` feedback-only probe asim + byproducts + analysis.json

Temporary ssh-side scratch (`%TEMP%\rg_before/after.txt`, `v2_rerun.log`)
deleted after the determinism check.

## Edits carried by the variants (V2 = the full set; V1/V3 subsets)

1. 8 ExternalStimuli (gamma 2.4/3.0 nA on the 4 RG halves; alpha 0 nA on
   L/R Hip+Knee MN ext) — asim `<Stimulus>` + aproj `<Stimulus>` clones.
2. SynapseType SynAmp: "c1 Commissural Inhibit" 2.749→1.1354, "V3
   Commissural Excite" 0.1→0.041304 (asim; aproj MaxSynapticConductance
   mirrors).
3. New type "RG Excite (comm drive, campaign)" SynAmp 1.1354; the 4
   RG→{c1,V3} connexions/links repointed to it.
4. New type "Afferent HC Excite (Contact, campaign)" SynAmp 0.011753; the 20
   contact-neuron→extensor-layer connexions/links repointed to it.
5. 4 ContactCount adapter Gain C: 20 (shipped, SI-amps unit bug) →
   7.1708 nA design-intent × campaign encoder ratio — see report §4.

No page/drawing CDATA was edited (verify: page count check passed; drawings
untouched by construction — stimuli have no drawings, type repointing is
attribute-only).
