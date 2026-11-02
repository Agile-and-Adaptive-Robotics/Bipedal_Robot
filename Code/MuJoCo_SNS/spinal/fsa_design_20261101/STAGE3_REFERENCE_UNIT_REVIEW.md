# Stage-3 reference, geometry, and balance-pathway audit

Machine: easteregg2. Read-only assessment of the saved 2 November source, II-only, and weak-all standing verifications. No body steps, neuronal steps, optimization, axis/pose changes, or metric changes were made. The original stage-3 rejection remains in force for all three conditions.

## Main finding

The reported approximately 64 mm “sway” includes a fixed-reference component, but subtracting that component would hide a real assisted-posture failure. The model's feet move anteriorly while its CoM is posterior to them; the rig supplies a substantial anterior force. In the source late window, all measured foot support is in the forefoot regions. This is bounded upright behavior under apparatus assistance, not accepted balanced standing.

| Saved condition | Old maximum ankle-reference distance mm | Late mean reference-offset radius mm | Late mean-centered RMS mm | Mean anterior rig force N | Actual foot load BW | Additional vertical spring support BW |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Source | 64.116 | 46.327 | 9.839 | 112.150 | 0.44716 | 0.05009 |
| II-only +0.02 mV | 63.989 | 47.076 | 9.112 | 110.727 | 0.45256 | 0.04601 |
| Weak-all | 64.513 | 47.280 | 9.558 | 114.298 | 0.45584 | 0.04327 |

“Late” is 5–8 s. The historical maximum uses the full post-warmup 0.6–8 s window, and the historical foot-load gate uses 4–8 s. These different windows must not be conflated. The source's 5–8 s mean offset is **46.308 mm posterior**, plus 1.316 mm native-left; its CoM also changes approximately 28.269 mm anteriorly between the first and last half-seconds of that window. It is not stationary merely because a centered RMS is small. The late horizontal rig-force RMS is 131.987 N, approximately 0.1715 of the historical BW normalization.

These are runner logical times. During the initial 0.6 s held preparation, neural updates run, but the runner calls `mj_forward`, not `mj_step`; actual physics elapsed time does not advance. The first free physical interval ends at logical 0.602 s. Actual activation is taken from aligned capture or the full mechanics diagnostic, never the legacy runner NPZ `act` field, which logs control demand despite its name. Preparation can change neuronal states and control; it must not be described as integrating muscle activation without actual-state evidence.

## 1. Global orientation and actual coordinate signs

The actual archived source model and initial saved pose independently establish anterior as **+world X**:

- Right calcaneus-to-toes vector: `[+0.178790,+0.002931,+0.000157]` m.
- Left calcaneus-to-toes vector: `[+0.178390,+0.012316,+0.000095]` m.
- Pelvis local +X in world coordinates: `[+0.998557,+0.042404,−0.032941]`.
- Torso local +X in world coordinates: `[+0.954125,+0.000977,−0.299408]`.
- Compiled pelvis fixed quaternion is approximately a +90° X rotation, not a 180° facing reversal.

The actual runtime right hip axis is `[0,0,−1]`; the raw converted XML has `[0,0,+1]`. Existing `runner.patch_xml` changes this axis. With anterior established from anatomy, finite differences of the correctly defined downward-thigh angle `atan2(knee_x−hip_x,−(knee_z−hip_z))` give derivatives −0.99920 and −0.99895 for increasing native right/left hip coordinates. Thus increasing the **current runtime** hip coordinate extends the thigh; native −7° is not evidence of anatomical hip extension. Initial world thigh angles are approximately 4.866° and 4.668° anterior to vertical. No coordinate correction is applied by this audit. The existing IK sign adapter and moving-point polynomial conventions must be considered before changing a target or an axis.

The source/runtime hip moving-point equality audit finds ten matching polynomials and no nonzero odd coefficients. If `q_native = −q_OpenSim`, composition `f(−q)` negates only odd coefficients; these particular polynomials are even and unchanged. Therefore an omitted odd-term sign transform is not a demonstrated problem in these ten paths. This does not certify every route or the full force-law alignment.

Initial reconstructed CoM agrees exactly with the logged CoM. Its talus midpoint X equals the recorded reference, −0.0763457535 m. Across sampled saved states, raw and aligned qpos pair exactly. FK versus logged CoM differs by at most 0.397 mm for the source, 0.384 mm for II-only, and 0.402 mm for weak-all. This small discrepancy is consistent with the documented solver-kinematics versus post-integration-qpos timestamp distinction; it cannot explain the large support-geometry separation below.

## 2. Original reference versus actual support geometry

`TaskBalance.bind` sets reference X to the initial talus midpoint and reference Y to initial CoM Y. The runner's historical sway metric uses the same convention:

\[
r_{old}(t)=\|COM_{xy}(t)-[x_{talus,initial},y_{COM,initial}]\|.
\]

It is a distance from a fixed ankle-based reference, not pure oscillatory sway and not distance from the current support polygon. Reconstructing it reproduces every original stage-3 maximum to 1e−12 m.

The source starts with a whole-foot projected mesh extent of approximately X `[−0.13955,+0.13595]` m and CoM inside that extent. Ground alignment has already changed both native ankle angles from zero to **+5.62137°**, with protected knees at −10°; initial foot flattening was performed.

By 5–8 s, the source's actual whole-foot mesh projected X extent averages `[+0.08925,+0.35863]` m, while mean CoM X is approximately **−0.12265 m**. In every 20 ms sampled late pose, CoM is outside even that complete mesh projection; mean outside-edge distance is 0.21557 m. The overlapping native foot-plane contact-point hull is narrower and farther anterior, with mean X span `[+0.29059,+0.33891]` m. It also excludes CoM throughout the late sample. CoM relative to the **moving talus midpoint** averages **−0.27607 m**, rather than the −0.04631 m distance to the old fixed talus X.

The mesh projection is a geometric extent, not a claim that every projected vertex can bear load. The contact-point hull uses actual native overlapping contacts evaluated from the saved pose, not recomputed forces or an inferred CoP. Both are descriptive cross-checks; neither changes acceptance.

Independent saved regional vertical forces average `[0,168.957,0,175.237]` N in right-hindfoot/right-forefoot/left-hindfoot/left-forefoot order. The source is supported entirely through the forefoot regions in that window. Native ankles average −24.928°/−24.840°, and knees approach their pre-neutral protective stop at approximately −1.993°. These are actual postural changes, not merely a static reference-label issue. Both candidate standing verifications show the same general geometric failure.

## 3. PT gains are current gains, with a separate recruitment/transmission chain

The delayed state contains horizontal displacement, filtered velocity, and optionally filtered acceleration. The actual command is

\[
c=\operatorname{clip}[-k_pe_d-k_dv_d-k_aa_d,\,-3,+3]\quad\mathrm{nA}.
\]

Therefore `kp=0.5` is **nA/m**, `kd=0.2` is **nA·s/m**, and `ka` has units nA·s²/m. The configured delay is 0.1 s, velocity/acceleration filtering defaults to 0.04 s, and each added PT source/relay has a 20 ms membrane. These gains are not directly acceleration-PD gains in s⁻² and s⁻¹.

At binding, the code finite-differences native CoM acceleration with a +0.02 activation perturbation about activation 0.05. It solves bounded nonnegative recruitment against each ±1 m/s² directional acceleration target, with vertical/cross-axis effects and activation regularization. Let that dimensionless activation vector be \(b_{ki}\). The actual electrical implementation then assigns

\[
\bar g^+_{ki}=2b_{ki}\;\mu\mathrm S,\qquad
\bar g^-_{ki}=0.25b_{k\oplus1,i}\;\mu\mathrm S.
\]

For source leak \(G_p=1\) µS, its isolated steady voltage is \(U_{p,k}=c_k/G_p\), with window fraction \(\alpha_{p,k}=\operatorname{clip}(U_{p,k}/5\,\mathrm{mV},0,1)\). The inhibitory relay additionally transforms that voltage through its own +8 mV conductance synapse. The resulting MN current is

\[
I_i^{PT}=\sum_k\left[\bar g^+_{ki}\alpha_{p,k}(8-U_i)
+\bar g^-_{ki}\alpha_{relay,k}(-5-U_i)\right].
\]

MN voltage is then divided by 5 mV and clipped to form activation demand; membrane shunting and activation dynamics remain. Consequently the native calibration `response @ b ≈ desired acceleration` does **not** imply that one nA produces one m/s². Even ignoring the inhibitory branch and other shunts, a low-current linearization gives approximately `(0.44–0.64)*b` activation per nA when MN operating voltage ranges from 2.5 to 0 mV. That range is an illustrative operating-point calculation, not a measured whole-network acceleration gain.

The source's actual late anterior command averages **+0.021412 nA**. Recorded `PT_BAL_XP` voltage averages **0.021602 mV**, only 0.00432 of its outgoing activation span; `PT_IN_XP` averages 0.034635 mV. The current chain is therefore operating far below full conductance activation. This motivates inspecting explicit gain scaling with the unchanged source references and limits, but does not justify presuming a corrective outcome or treating `kp` as physical stiffness.

There is also a major stage dependence. Saved source/candidate stage-2 calibrations achieve only approximately 7e−10 to 6e−9 m/s² in the desired directional CoM responses, with maximum basis activation 0.0006593. Stage 3 achieves approximately ±0.984–0.998 m/s² with maximum basis activation 0.3. Both report numerical rank 3, which is not meaningful evidence of air-stage authority at the former absolute scale. Internal muscle forces cannot provide useful free-air whole-CoM acceleration here. The same `kp/kd` recipe thus has very different mechanically calibrated PT→MN matrices in air and on the ground. Passing stage 2 does not test loaded PT balance authority.

## 4. Load-gate coherence

The old assessor requires `load > 0.5`, calculated over 4–8 s. Its source value is 0.49156. The once-counted 5–8 s foot load is lower, 0.44716 BW, because this is a different, later window. The legacy logger's duplicate heel/toe accumulation is already halved in this metric; that known duplication is not by itself a new explanation for the present value.

The source actually receives prescribed upward support of **0.5 BW**, plus mean spring support of **0.05009 BW** over 5–8 s. Their sum with measured foot support is 0.99725 BW. A stationary model receiving approximately 55% vertical assistance should ordinarily load its feet by approximately 45%, not more than 50%. Requiring more than 50% while retaining that assistance favors a transient or a different vertical state, rather than specifying the intended equilibrium.

There is a small normalization distinction: the historical BW is 769.737726 N, including 0.1 kg of fixed-body mass, whereas the moving pelvis-subtree BW is 768.756726 N. This 0.127% effect cannot explain the approximately 5% support mismatch or the posture failure.

A future declared load contract should use matched-window **signed measured** apparatus/explicit forces and the moving-body weight, with vertical acceleration accounted for:

\[
F_{feet,z}+F_{explicit,z}+F_{rig,z}-W_{moving}=M_{moving}\ddot {COM}_z.
\]

It should separately require bilateral load, bounded vertical drift/velocity, and a stated assistance level. This is a proposed specification correction, not retrospective approval of these trials. The current old-gate results remain unchanged.

## 5. Defensible next checks

1. Retain and label the historical ankle-reference metric. Separately report mean setpoint error, variation about that mean, drift, actual moving-foot geometry, and apparatus forces. Do not subtract a fitted offset and call the old trial a pass.
2. Keep explicit current units and the complete PT transmission chain when interpreting corrective gain probes. Larger current gains can affect loaded standing much more strongly than suspended stage 2; they can also change recruitment during the held-pose warmup.
3. If a new stance initialization is considered, solve a declared anatomical/contact task: source joint limits, protected slightly flexed knees, actual flat heel/forefoot geometry, CoM inside the intended support region with an explicit margin, unilateral/friction-feasible contact forces, and the intended measured support. Initial IK alone is not a static muscle-force equilibrium or a controller validation.
4. Verify any selected unchanged controller/stance recipe again through the required air and standing stages and a fresh process before ground stepping. Do not change runtime axes, contact geometry, muscle strengths, or acceptance thresholds as an unreported way to pass.

Evidence: [complete saved-state audit](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/STAGE3_REFERENCE_UNIT_AUDIT.json>), generated by [audit_stage3_reference.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/audit_stage3_reference.py). Actual source equations are in [balance.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/balance.py), historical metrics in [runner.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/runner.py), and stage acceptance in [curriculum_metrics.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/functional_gait2392/curriculum_metrics.py). The JSON indexes exact model/configuration/result hashes, per-window statistics, and all three original outcomes.

## 6. Subsequent bounded PT-gain verifications

The parent ran two explicit corrective conditions with the II-only controller and unchanged reference, physics, and original stage-3 gate. Both pass that historical gate, but the identical offline geometric audit shows that neither establishes adequate foot-relative standing:

| Condition | Original sway mm | CoM behind moving talus mm | Mean distance outside complete foot projection mm | Mean-centered CoM RMS mm | Anterior rig force N | Right/left hindfoot forces N |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| kp 5, kd 2 | 59.45 | 259.60 | 199.13 | 7.25 | 109.95 | 0 / 6.48 |
| kp 25, kd 10 | 50.57 | 145.61 | 83.20 | 2.02 | 84.67 | 0 / 103.78 |

Both have CoM outside the complete projected foot-mesh hull at every sampled late state. The higher-gain condition does improve the physical posture and reduces late mean CoM drift to approximately 0.97 mm, but its right foot remains forefoot-only. Its left forefoot still carries approximately 277.02 N, versus 103.78 N at the left hindfoot; the right forefoot carries approximately 200.53 N. It is a less variable apparatus-supported posture, not demonstrated foot-supported equilibrium.

Critically, the kp 25/kd 10 condition has **downward** vertical spring force throughout 5–8 s: mean −0.25661 BW, range −0.34909 to −0.17972 BW. With prescribed upward +0.5 BW, net apparatus vertical force is +0.24339 BW, and feet carry 0.75524 BW. This does **not** show that the model needs only 24% harness support: the upward support and downward restraint oppose one another, while substantial horizontal assistance remains. The lower-gain spring crosses both signs, approximately −0.19878 to +0.15157 BW, averaging +0.00883 BW.

The old gate outcomes remain valid as results under their original definition, but they are insufficient grounds for promoting either condition to actual ground stepping. Evidence: [PT_GAIN_REFERENCE_AUDIT.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/PT_GAIN_REFERENCE_AUDIT.json>). No new body runs were performed for this postprocessing.
