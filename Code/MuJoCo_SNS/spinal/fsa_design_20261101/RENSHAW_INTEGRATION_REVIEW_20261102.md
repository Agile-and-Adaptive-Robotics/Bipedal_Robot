# Optional RC integration: released, default off

Session 1, easteregg2. The parent explicitly released integration after the ground-v2 study and fresh gates completed. **Only the reviewed optional caller block was then applied; its default is None/off.** Controller, balance, physics, gains, and training sources remain unchanged. The distinctly identified RC-enabled candidate now has its own fresh stages 1–4 acceptance, and stage 5 passes at contact-PF scale 0.25; full-strength stage 5 fails clearance, and the selected-scalar stage-4 basin passes 8/8. Final S4 and quarter-strength S5 repeats have byte-identical traces and models; the next lower-lift test fails clearance.

The exact [review patch](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/physical_eval_autogenic_rc_v1.patch) inserts only an optional `cfg.renshaw_autogenic_profile` block after existing FSA composition and before SwingDiagnostics/native-step capture or later PT body binding. It requires the reviewed integrated Hunt source, invokes the existing separate RC installer, and records actual helper/profile file paths and SHA256 values in `capture` and `source_hashes`. The existing aligned-capture provenance already receives that dictionary. Missing/None does not import or install the helper.

**Archived pre-RC source:** `c598899ec1e2025180c18062e84791f8c0e5e05a11d6847434304338306d29f8`, preserved byte-for-byte in [pre_rc_physical_eval.py](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_patch_preparation_v2/pre_rc_physical_eval.py>) with `PRE_RC_SOURCE_MANIFEST.json`. The exact reviewed diff was applied only after explicit release. Installed `physical_eval.py` SHA256 is **`494ff3dafa0e8a1b20d01de7167273468f9ac5088a82d8db133b9bb091141602`**. The [preparer](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/prepare_renshaw_integration.py) itself never patches the runtime; its `--baseline` mode confirms that actual on-disk source equals exactly the proposed block plus the archive.

Actual old/new constructor checks cover stages 1–4. None/default has identical membrane/ionic fields, initial states, and 128-step electrical trajectories; existing gain buffers remain exact. The active optional candidate persists through commands 0/0.5/1/0, passes existing local/timing/coincidence checks, records both hashes, and retains the motor collateral under deafferentation. No body integration occurs. [Version-2 constructor/provenance proof](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_patch_preparation_v2/OPTIONAL_PATCH_PROOF.json>).

## Dense peer inhibition changes the operating point

Every one of the92 RC cells receives **45 distinct other RC inputs on the same side**, each 0.5 µS, E=−5 mV, activation window0–5 mV. There is no RC self edge or cross-side peer edge. Each leg has 46 RC cells. The 4,140 peer pairs are therefore a dense same-side circuit, not one weak recurrence per joint. Their endpoints and gains were preserved.

For the complete MN/RC subsystem, drive each MN with 2.5 nA, while other controller inputs are silent. Let `u` and `v` denote the common steady MN and RC voltages. With the new autogenic gain `g=0.0260416667 µS`, the actual source equations give

\[
u(v)={2.5+g(v/5)(-5)\over1+g(v/5)},\qquad
0=-v+(u/5)(8-v)+45(0.5)(v/5)(-5-v).
\]

The solution is **u=2.493701679 mV, v=0.161372543 mV**. Autogenic feedback is **−0.006298321 nA/MN**, not the isolated pair's−0.10 nA; all-peer inhibition into each RC is **−3.748067147 nA**. The full 184-cell,4,324-pair source-derived SNS fixture reproduces those values within4.45×10⁻¹⁵ mV and matches dense versus accelerated sparse output exactly. This is the complete retained MN/RC operator with other circuits silent, not the whole1,033-cell controller or a body claim.

Abort-only voltage guards checked the 3 s electrical fixture without clipping. MN values stayed0.111111…2.493702 mV; RC values stayed−0.174011…0.183879 mV, including permitted inhibitory undershoot relative to rest. Continuous invariant bounds are MN[−5,2.5] and RC[−5,8] for this input condition. The worst saturated per-cell Euler conductance factor is1.6333, so2 ms is **not a global convex-update/invariance guarantee over arbitrary saturated states**. The observed bounds alone are not a stability proof.

## Full184-state derivative and sampled poles

For every actual operator, the complete Jacobian retains both shunting and presynaptic driving-force derivatives:

\[
A_{ii}=-{G_i+\sum_jg_{ij}s_j\over C_i},\qquad
A_{ij}={g_{ij}s'_j(E_{ij}-V_i)\over C_i}\quad(i\ne j).
\]

At the measured coactive equilibrium, **all 184 continuous poles have negative real part**, spanning−848.290694…−23.361778 s⁻¹. The actual2 ms Euler matrix `I+0.002A` has spectral radius **0.953276445**. Its fastest/common-mode pole is **−0.696581387**, inside the unit circle. Finite differences of the actual184-cell SNS one-step map match the analytic sampled matrix within **1.78×10⁻¹¹**. Thus local stability holds for this operating point; it was not inferred from isolated two-cell poles.

All matrices, ordered names, equilibrium, and eigenvalues are preserved in [COACTIVE_LINEARIZATION.npz](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_patch_preparation_v2/COACTIVE_LINEARIZATION.npz>); the proof JSON includes every pole. This is not full-controller, contact, or gait stability.

## Disposition

The inherited per-edge0.5 µS becomes a large aggregate current when copied to45 peers. A declared total-peer current contract would be a **new functional-design choice**, not an automatic biological correction. For example, a hypothetical0.5 nA total inhibitory budget at source/recipient2.5 mV would derive `0.5/(45*0.5*7.5)=0.002962963 µS` per edge. That example was **not installed, optimized, or adopted**. The observed equilibrium is locally stable, so these results do not justify a stability-driven normalization by themselves. The isolated autogenic profile remains correctly labelled as an isolated contract with strongly operating-point-dependent full-circuit expression.

The optional block does nothing unless the named profile is supplied. Parent full-body parity is complete: pre-RC `ground9_s4` versus installed-default-None `ground9_none_s4` has byte-identical aligned NPZ/model MJB and identical assessment/contact/clearance fields. The subsequent RC-enabled candidate was bound and verified independently; it did not inherit the old body pass. [Installed-source/electrical verification](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_installed_v1/OPTIONAL_PATCH_PROOF.json>).

Fresh RC-enabled [S1](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/ground9_rc_s1.result.json>), [S2](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/ground9_rc_s2.result.json>), [S3](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/ground9_rc_s3.result.json>), and [S4](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/ground9_rc_s4.result.json>) pass. The full architecture has **1,033 actual cells; 1,037 names include four logging aliases**. Stage 1 uses its explicit no-PT protocol projection while preserving the motor collateral.

At [figures/ground9_rc_s4](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/figures/ground9_rc_s4>), independent [whole-mesh/once-counted-load checks](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/figures/ground9_rc_s4/aligned_clearance/ALIGNED_CLEARANCE_CHECK.json>) pass with 5/4 clear unloaded swings. [Actual sparse-current closure](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/figures/ground9_rc_s4/SPARSE_CURRENT_REVIEW.json>) passes all 7,700 samples with zero reconstructed voltage-update error. In the [aligned 5–13 s support window](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/figures/ground9_rc_s4/plots/SUPPORT_DECOMPOSITION_5_TO_13.json>), net vertical assistance is **48.8484% BW**, with **15.0322% BW horizontal guidance RMS**. This is supported in-place stepping, not unrestricted balance or walking. No tendon, joint-damping, muscle-damping, contact, or other physics change accompanies the neural repair.

[Quarter-strength stage 5](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/ground9_rc_s5_contact0p25.result.json>) is accepted at `contact_pf_scale=0.25`; its [independent aligned gate](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/figures/ground9_rc_s5_contact0p25/aligned_clearance/ALIGNED_CLEARANCE_CHECK.json>) passes with 5/4 clear unloaded swings, aligned net assistance 48.766126% BW, and horizontal restraint RMS 115.309014 N. [Full-strength stage 5](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/ground9_rc_s5_contact1p0.result.json>) fails with only 2/2 qualifying clear swings despite 5/5 knee cycles. The [stage-4 basin result](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/training/finalist_ground9_rc_basin_v1/BASIN_RESULT.json>) passes all **eight cases perturbing ten selected scalars by ±1%**, with fixed RC gains and 5/4 clear swings in every case. It is not a ten-case test, exhaustive parameter basin, or whole-model robustness proof.

The nominal RC S4 repeat and quarter-strength S5 repeat are complete, with byte-identical aligned NPZ and MJB files. The next 35%-lift quarter-strength case fails with two qualifying swings per foot. Exact hashes are in D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\FINAL_RC_FRESH_HASH_CHECK.json. No completed support-weaning claim is made. The fixed −0.10 mV design target remains an isolated-pair contract, not a claim of a −0.10 mV effect in the full body or of complete physiological Renshaw validation.
