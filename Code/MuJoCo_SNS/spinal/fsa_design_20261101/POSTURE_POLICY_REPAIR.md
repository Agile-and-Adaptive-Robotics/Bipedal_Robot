# Optional Hunt standing-input repair

Implemented in the current `physical_eval.py`; both options default to false. Historical configurations and source results remain separate. The physical model, joint ranges, force-law parameters, and acceptance metrics are unchanged.

The selected precise policy is:

```json
{
  "suppress_inherited_posture_mn": true,
  "balance": {
    "replace_inherited": true
  }
}
```

Merge those keys into the complete saved configuration; the nested balance object still needs its existing gains, delay, and bilateral settings. The MN-only policy zeroes **88** slots whose source is exactly `POSTURE` and whose target starts `MN_`. It changes `original_g`, then the live sparse/dense arrays, before final integrated-profile application. Later command-dependent updates therefore retain the zeros. The source neuron remains active, and its two `POSTURE → RG-E` slots retain their exact baseline values. In this actual Hunt1033 composition those two live conductances are already zero; preserving them does not establish a nonzero tonic RG drive. All 92 direct `POST_<act>` applied-current inputs are preserved. The existing `replace_inherited` option independently zeros inherited balance input ports; it does not zero global POSTURE-to-MN connections by itself.

The separate broader diagnostic is `suppress_inherited_posture=true`: it sets only the external `POSTURE` input current to zero before the original neural step. It suppresses every active output of that source and is **not** the selected MN-only repair. In a different builder with active POSTURE-to-RG connections it would also remove their tonic drive; those connections are currently zero in this Hunt composition. Enabling both policies raises an error. If the full-source diagnostic were enabled during an already-running state, its existing source voltage would decay; the policy does not hard-reset neuronal state.

## Why the separation matters

At the original ground calibration pose, the recreated native activation-response matrix exactly reproduces all four saved directional PT calibration responses. An isolated 1 nA BAL-source input, propagated through the actual compiled MN conductances with shunting about a 2.5 mV MN operating point, predicts these native-world CoM acceleration increments:

| Inherited source | Anterior acceleration m/s² | Intended task direction |
| --- | ---: | --- |
| BAL_PF | −0.226286 | Anterior for posterior CoM error |
| BAL_DF | +0.063283 | Posterior for anterior CoM error |

Thus the combined ankle/hip projection has the wrong sagittal sign at that reference. This is a local static operator result, not a closed-loop rollout or a claim about every later contact state. The current reconstruction replays all actual recorded MN synaptic-current groups to at most 1.33e−15 nA.

An additional matched-baseline native check uses MN baseline 0.25 mV, corresponding to the mechanical calibration's activation 0.05, and evaluates the resulting combined activation increments with native nonlinear contact/constraint forward dynamics. It confirms sagittal accelerations of −0.281536 m/s² for BAL_PF and +0.097805 m/s² for BAL_DF. These evaluations integrate no time. Other acceleration components differ substantially from the summed local linear prediction, underscoring the limited scope of the fixed-pose matrix.

There are also two live posture drives: direct `POST_<act> = 6*a_static` and 88 chemical `POSTURE → MN` inputs. The initializer returns an activation target, but the neural mapping is not identity: even direct-only current `6*a_static` at 1 µS leak and a 5 mV activation readout gives `1.2*a_static`, before other inputs and shunting. In the II-only standing verification, inferred soleus-right static target is zero, yet its actual mean activation is 0.2452; global POSTURE supplies 0.7456 nA to that MN. This identifies a substantial additional drive, rather than proving that zero soleus activation would itself be an adequate standing solution.

The repair removes the additional global MN posture drive without adding a runtime inverse-dynamics allocator or cancelling measured synaptic currents. It does not calibrate the remaining direct input to a perfect activation target. Body verification remains necessary.

## Verification and provenance

- Default and explicit-false trajectories match exactly over 128 accelerated electrical steps.
- MN-only masked conductances, dense/sparse slots, and actual MN currents remain exactly zero under commands `0, 0.5, 1, 0`.
- The POSTURE source remains active; both RG-E paths and all direct POST inputs are preserved.
- Learned PF, FF, reflex, and biarticular buffers, all other neural parameters, and archived physical-model arrays remain unchanged.
- The full-source diagnostic gives exactly zero source voltage/current from its zero initial state; conflicting flags reject.
- No body simulation or fitting was launched for this proof.

The existing integrated profile needs no rebind, but the **caller configuration flag is part of the new candidate's identity**. The old profile alone does not encode the new mask; preserve the complete configuration, runtime-source hash, and resulting actual graph. Do not reinterpret an earlier compiled-graph artifact as the masked network.

Evidence: [POSTURE_INPUT_POLICY_PROOF.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/input_policy_v1/POSTURE_INPUT_POLICY_PROOF.json>), [pre-change source manifest](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/input_policy_v1/PRE_MASK_MANIFEST.json>), and [STANDING_PATHWAY_AUDIT.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/STANDING_PATHWAY_AUDIT.json>). Reusable checks: [verify_posture_policies.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/verify_posture_policies.py) and [audit_standing_pathways.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/audit_standing_pathways.py).

## Subsequent supported-standing verifications

The parent ran three conditions with `suppress_inherited_posture_mn=true`, `replace_inherited=true`, and PT gains kp 25 nA/m and kd 10 nA·s/m. All pass the original stage-3 gate. The independent offline audit also finds genuinely improved foot-relative posture, rather than only a favorable fixed-reference metric:

| Condition | Actual contact-hull margin, mean/minimum mm | Mean-centered CoM RMS mm | Horizontal rig RMS N | Actual foot load BW | Downward spring force BW |
| --- | ---: | ---: | ---: | ---: | ---: |
| II-only, coupled initializer | 47.57 / 40.65 | 1.667 | 34.61 | 0.70722 | −0.20858 |
| II-only, contact initializer | 28.42 / 19.18 | 2.156 | 42.25 | 0.72578 | −0.22712 |
| Best23 local targets, coupled initializer | 50.22 / 44.12 | 0.907 | 30.63 | 0.70959 | −0.21092 |

Statistics cover 5–8 s. At every 20 ms sampled late pose, CoM lies inside both the complete projected foot mesh and the native overlapping contact-point hull. All four heel/forefoot regions carry substantial mean load. For best23, right hindfoot/forefoot forces are 146.40/98.24 N, and left hindfoot/forefoot forces are 190.12/111.44 N. Native ankle means are −2.11°/−2.36°, compared with approximately −25° in the earlier failed posture. Its maximum observed knee extensions remain negative, −1.627°/−1.779°. The late CoM change is approximately −2.03 mm anterior-axis displacement across the first/last half-second averages, rather than the large earlier transient drift.

Among these three higher-PT conditions, best23 was selected for repeat checks; the subsequent lower-PT result below supersedes that gain selection. This is evidence for **supported** standing. The prescribed upward force remains +0.5 BW, while the spring pulls downward by approximately 162.35 N, or −0.21092 BW. Their net is +0.28908 BW, but this is opposing apparatus loading, not evidence that only 29% harness assistance is necessary. Lateral and rotational apparatus constraints also remain. No harness-free balance or stepping claim follows from these standing results.

The three conditions differ in local targets and/or initializer, so this table is not a single-factor causal estimate of the mask alone. Their exact configurations, model hashes, original assessments, signed forces, and geometric checks are indexed in [MN_ONLY_STANDING_GEOMETRY_AUDIT.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/input_policy_v1/MN_ONLY_STANDING_GEOMETRY_AUDIT.json>). This agent performed only offline postprocessing of those completed body runs.

## Constant lower-PT standing checkpoint

The completed [pt_low_s3 configuration](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/pt_low_s3.json>) uses constant PT kp 5 nA/m and kd 2 nA·s/m with the selected local-target/posture-policy recipe. It passes the unchanged stage-3 gate and improves the independent geometric assessment further:

- CoM is inside both sampled contact and full-foot mesh hulls throughout the late window. Actual contact-hull margin averages **83.08 mm**, with minimum **81.53 mm**.
- Mean-centered CoM RMS is **0.309 mm**; anterior-axis change between first/last late half-seconds is **+0.245 mm**.
- Horizontal rig-force RMS is **9.126 N**, approximately **1.186% BW**.
- Right hindfoot/forefoot loads average **179.50/65.65 N**, and left loads **232.62/78.83 N**. This is not forefoot-only support.
- Native ankle means are **+1.259°/+1.191°**. Knees remain near the −2° protective bumper; maximum observed extensions are **−1.935°/−1.787°**.

Vertical assistance still has opposing components: +0.5 BW prescribed upward force and −0.22442 BW downward spring force, giving +0.27558 BW net, with feet carrying 0.72311 BW. This remains an apparatus-supported stance. No gain scheduling is justified by a standing failure here: this lower-gain standing condition itself passes. Same-parameter curriculum and fresh/ground bookkeeping belong to the parent experiment record; this offline audit does not independently establish those other stages.

Evidence: [PT_LOW_STANDING_GEOMETRY_AUDIT.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/input_policy_v1/PT_LOW_STANDING_GEOMETRY_AUDIT.json>). No additional body or neuronal integration was performed to produce that audit.
