# Named conditions for the all-controller audit

Session 1 designation, **easteregg2**, 2026-11-01. These are audit anchors, not a claim that every anchor has completed Hunt training or locomotion. Do not select a folder by modification time or combine parameters across versions. Preserve the saved source/configuration and its own acceptance criteria. A new FSA design is a separately named opt-in condition.

## Rat134: curated receptor-geometry baseline, with refinement12 as a diagnostic comparator

**Primary model and reference condition:** `Curated134_ReceptorGeometry`, with accepted deafferented A preserved and the declared low-drive receptor-reference condition.

- Model: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\Rat134_ReceptorGeometry_Standalone.asim`
- Editable project: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\Rat134_ReceptorGeometry.aproj`
- Manifest: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\MANIFEST.json`
- Checkpoint proof: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\FINAL_CHECKPOINT_PROOF.json`
- Condition definitions: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\CONDITION_DECISIONS.json`
- Runtime adapter: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\scripts\rat134_receptor_runtime.py`, function `variant_controller()`. The base `rat134_core.controller()` alone is not the curated receptor/geometry condition.

Source standalone SHA256: `e3206341b7f89f3cb97584359ca83f136b3c7881b9fc0cf6d6ecaf9853606888`. It is 134 neurons/274 edges. The reference qualifies 35 informative channels across seven cases, with four fixed low-drive channels and five optional graph-unreachable channels separately recorded. Accepted A is preserved exactly. Keep the strict all-44 failure distinct from this declared reference condition.

**Secondary experimental condition:** `prototype_role_policy_B/refinement12`, nominal trial3, for auditing proposed role policy, consumed gains, and the failure of its holdouts.

- Manifest: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\prototype_role_policy_B\refinement12\MANIFEST.json`
- Record: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\prototype_role_policy_B\refinement12\EXPERIMENT_RECORD.md`
- Exact consumed fields: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\prototype_role_policy_B\refinement12\CONSUMED_FIELD_LEDGER.json`
- Result: `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\prototype_role_policy_B\refinement12\RG_HOLDOUT_AND_A_RESULT.json`

Trial3 passes the nominal four-RG-cell gate but fails both cadence and both amplitude holdouts. It is an intermediate, not accepted Hunt B. PF fitting and actual-body trials were not launched. Do not substitute it for the primary curated model or silently install its gains.

## Central-spiking Gait2392: dynamic_readout_v4 TimingZeroSFA

**Primary saved condition:** the exact `TIMING_CANDIDATE.json` under `Gait_dynamic_readout_v4`, interpreted with the later checkpoint and versioned burst metric.

- Configuration: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\TIMING_CANDIDATE.json`
- Neural checkpoint: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\CENTRAL_CHECKPOINT.json`
- Physical checkpoint: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\PHYSICAL_CHECKPOINT.json`
- Fresh-process proof: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\PHYSICS_FRESH_VERIFICATION.json`
- Versioned neural metric: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\BURST_METRIC_V2.json`
- Source snapshot: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\SOURCE_SNAPSHOT.zip`
- File manifest: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\FILE_MANIFEST.json`
- Record: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\EXPERIMENT_RECORD.md`

Controller class: `TimingZeroSFA`. Parameters: `rg_adapt_jump_nA=0`, `rg_readout_C_rate=5`, `pf_rate_scale=1`, `rg_h_tau_scale=0.9496883348174533`; 16 cells per central population, 0.2 ms internal spiking step, 2 ms neural/physics interval, and 25 ms spike-history constant. Census: 256 central spiking cells and 778 graded cells. This named condition completed the recorded supported-standing protocol (0.6 s pinned preparation plus 7.4 s integration), with 82 arrays exact in a fresh process. Its measured mean upward support was 0.5937% BW, while lateral and rotational guidance remained.

`TIMING_CANDIDATE.json` is a frozen numerical artifact; its original pending metadata is superseded by `CENTRAL_CHECKPOINT.json` and `PHYSICAL_CHECKPOINT.json`. Do not edit the frozen configuration just to synchronize its descriptive status.

`Gait_central_v1` and `Gait_central_inverse_v3` are historical lineage/comparison conditions. Audit their architectural differences when useful; they are not the default physical-result anchor and their settings must not be mixed into v4.

## CAD262: named support085 mixed pressure-pair condition

**Primary reproduction manifest:**

`D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\cad\REPRODUCTION_MANIFEST.json`

It designates the following saved 8 s condition:

- Base model: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\ah001_neural.xml`
- Exact runtime protocol: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped\support085\mixed\mixed\LIVE_PREREGISTRATION.json`
- Result: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped\support085\mixed\mixed\cad.json`
- Assumed pressure-actuator configuration: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped\support085\mixed\ASSUMED_PAIR_CONFIGURATION.json`
- Neural-pair normalization: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped\RECRUITMENT_PREREGISTRATION.json`
- Source manifest: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped\FRESH_SOURCE_MANIFEST.json`
- Measured condition and remaining mechanics: `D:\Users\Ben Bolen\NeuroCampaign_20261031\review_20261101\cad\README.md`

Label it **CAD262, two neural pressure-controlled vas_med BPAs plus 52 ideal-tension hybrid routes, support085**. The pressure pair has zero POST clamp current; the other routes retain their documented allocator and force floor. Measured upward support is 87.3372% BW. Preserve the saved stance/geometry issues in this condition; the later mechanical repair and rigid-loop fixtures are separate diagnostics, not replacement locomotor checkpoints.

This is the designated audit condition via its reproduction manifest, not a newly invented claim that a CAD `CHECKPOINT.json` or all-BPA model already exists.

## Reduced Deng14: fixed-tau neural mechanism reference

**Primary source/condition:** the direct ODE implementation with fixed NaP h time constant, rather than the legacy SNS-NumPy clock or a merged ODE/Simulink parameter set.

- Source: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\deng_cpg_ode.py`
- Saved reference: `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\deng_cpg_ref.mat`
- Separate Simulink implementation: `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos\sns_build_deng_cpg.m`
- Legacy comparison implementation: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\deng_cpg_numpy.py`
- Applicability/unit ledger: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\non_gait_ledgers\reduced_Deng14.json`

The source reference call is `simulate('fixed', t_end=20000.0, dt=0.1, i_pulse=10.0, t_pulse=100.0, pulse_ms=1.0)`. Time arguments here are **milliseconds**: 20 s duration, 0.1 ms step, and a 10 nA RG-extensor pulse starting at 100 ms for 1 ms. The NaP h constant is 350 ms. The circuit has 14 neurons and two hip MN outputs. No designated whole-body/locomotion checkpoint has been identified for this reduced circuit. ODE and Simulink PF-to-MN conductances differ and must remain explicitly versioned.

## Reproduction policy

Use these exact condition files and source manifests. When a live workspace file no longer matches the recorded hash, use the saved source snapshot in an isolated audit copy or record the discrepancy. Do not overwrite the active checkout, rewrite historical acceptance decisions, launch duplicate jobs, or treat the new opt-in FSA adapters as part of an old saved winner.
