# Independent CAD evidence review

Reviewed on easteregg2, 16 November 2026. The canonical Chapter 6 and Appendix D were read in full, together with the prior revision records, portable Session 1 index, exact CAD experiment records, retained runner, scoring script, selected configuration, recorded arrays, and local perturbation cases. This audit did not integrate a mechanical model, fit a controller, or change any scientific source.

## Result

The selected `clearance77p5_damped` condition remains supported by the local evidence. The existing manuscript's nominal support, qualifying-swing counts, minimum-knee statement, replay agreement, and five-of-eight local outcome do not require numerical replacement. It is a separate hybrid-control demonstration using ideal reference tensions. It is not a dynamically calibrated 54-BPA result, an unmodified neural-controller result, or harness-free walking.

The primary directory is:

`D:\Users\Ben Bolen\CAD3DOF_20261102\contact_handover_20261103`

The model is:

`D:\Users\Ben Bolen\CAD3DOF_20261102\anatomical_revision\models\toe_base_passive_nominal.xml`

The selected raw and configuration files are `runs\clearance77p5_damped.npz` and `.json`. No later CAD acceptance package was identified in the CAD directory or the newer controller-training records searched. The portable records are accurate for this selected condition; their dated statements were not treated as current evidence by themselves.

## Independent calculations and scoring windows

`CAD_audit.py` reproduces these calculations without invoking any runner. `CAD_AUDIT_RESULTS.json` records exact source hashes, replay comparisons, individual qualifying intervals, and the eight local cases.

| Quantity | Verified value | Source and window |
|---|---:|---|
| Completed integration duration | 27.99999999999536 s | Selected run JSON; requested 28 s |
| Stored sampling | 2,800 samples at 10 ms | Raw NPZ; 0–27.99 s |
| Final-rung scoring samples | 1,200 | `time >= 16`, through 27.99 s |
| Model mass and weight denominator | 75 kg; 735.75 N | Current XML compiled with MuJoCo 2.3.7; gravity `[0,-9.81,0]` |
| Mean upward harness force | 570.20625 N = 77.5% BW | Raw `rig_force[:,1]`; final rung |
| Mean foot-normal load, right/left | 78.5830481/86.9712333 N | Raw `contact_loads`, summed over its two foot regions; final rung |
| Qualifying intervals, right/left | 7/6 | Raw load and mesh-height arrays; final rung |
| Alternations | 12 | Time-ordered qualifying interval onsets; final rung |
| Logged maximum knee angle, right/left | −11.2047449°/−5.1661126° | Raw `knee_coordinate`; whole stored trial |
| Maximum knee angle at every physics step, right/left | −11.2041350°/−5.0995415° | Saved online accumulator in selected JSON; whole 28 s integration |
| Final-rung logged maximum knee angle, right/left | −14.0508627°/−6.0059392° | Raw `knee_coordinate`; final rung |
| Minimum mesh height, right/left | −16.1335400/−12.9055986 mm | Raw `mesh_clearance`; final rung |
| Maximum mesh height, right/left | 35.3954053/53.5396875 mm | Same raw array and window |
| Forward displacement | −74.7140122 mm | First-to-last stored `qpos[:,0]`; whole trial |
| Mean absolute rig moment, x/y/z | 10.7995483/0.9292861/110.0753435 N m | Raw `rig_torque`; final rung |
| Horizontal rig-force RMS, x/z | 1.3505712/2.6951636 N | Raw `rig_force`; final rung |

The interval criteria match `source\measure.py`: whole-foot minimum mesh height above 5 mm, ipsilateral normal load below 7.36 N, contralateral normal load above 15 N, and last-minus-first sampled interval duration at least 80 ms. The saved overall gate additionally requires trial completion, at least six seconds of final-rung observation, at least three intervals on each side, at least five alternations, per-physics-step knee protection, and nonnegative sampled lateral mesh separation. These are engineering acceptance conditions.

The 5.10° minimum-flexion statement is the least flexed knee over the complete integration, specifically the left knee. Its right-side counterpart is 11.20°. The new final-rung plot cannot itself establish the all-trial 2 ms minimum: the raw saved plot arrays have a 10 ms interval. The online every-physics-step accumulator is identified separately instead of conflating its 5.09954° value with the 5.16611° minimum in the saved 10 ms trace. The guard checks predicted angle using positive velocity, so a failed predictive test need not reach the limit in the last recorded angle.

No adopted result was silently renormalized. The selected run's nominal `support=0.7` option is not the measured final support: the spring and cap protocol produce the recorded 77.5% final-rung force. The complete-trial mean support is approximately 81.7761%, including the earlier high-support period, and is not substituted for the final-rung value.

## Model, actuation, and control identity

All eight saved source hashes match the actual model XML, route JSON, frozen controller, dynamic-route implementation, robot kinematics, network builder, network parameters, and exact-fast backend. The frozen runner's SHA256 is `a54c0a7d0f23a497a07fa3b067937b5fd14b0a8b19a0ad3a140aa4db70266ae4`; the model XML's SHA256 is `cb44a27e9d4262786a51938e0bb9a08899556e18851be82e3b8d931c60054789`. The full paths and all hashes are retained in the audit JSON.

Current compiled geometry reproduces the stored final-rung minimum mesh heights at every 10 ms pose to a maximum difference of `3.8163916471489756e-16` m. This verifies the geometry used in the quantitative figure against the current local model without replaying its dynamics. The model has 54 actuators with 30 ms first-order activation, and the runner describes the force law as reference tension capped by the existing actuator capacities. KompareXXVII pressure-actuator sizing and calibrated dynamic BPA laws are not used.

The runner constructs 262 graded cells and 230 input channels. SNS RG activity sets stepping timing, while the engineering allocator supplies posture, stance-load, and mesh-clearance targets. `motor_target_current` explicitly compensates instantaneous background MN synaptic current and imposes the target. It therefore limits ordinary downstream PF/reflex authority, as the manuscript already states. Ia/II/Ib and contact inputs remain in the circuit; their presence does not make the full MN allocation a purely neural controller.

The handover filter is evaluated at the 10 ms allocation interval using a 10 ms time constant. Contact confidence is the smoothstep of force between 1 and 5 N; the normalized phase preference is multiplied by confidence, without renormalizing confidence into a complete load request. The swing-clearance task acts through allocated hip/knee muscle targets. It does not apply an external foot lift. Physical support remains an external pelvis rig with vertical, horizontal, and attitude terms.

## Repeats and local perturbations

All 26 NPZ arrays, not merely the plotted signals, compare exactly between the selected run and each retained fresh-process and reference-backend result. The corresponding JSONs retain the same model XML hash. This audit compares saved results; it does not claim a new dynamic replay.

The eight local raw traces independently reproduce the recorded qualifying-interval counts. Cases 01, 02, 04, 05, and 07 complete and pass. Cases 00, 03, and 06 stop at the predictive knee margin at 18.880, 24.466, and 19.046 s, respectively. The perturbations cover selected shared controller settings, six handover/clearance settings, and three passive-parameter families. They do not perturb every fixed neural matrix entry, the assistance protocol, or the safety bound. This supports limited local verification, not broad robustness or global stability.

The retained 1 ms physics-step test completes with 7/6 clear intervals, but its trajectory is not exact and its least knee-flexion margin is approximately 4.55°. It keeps the 2 ms neural and 10 ms allocation/log intervals. It is an integration-sensitivity test, not another exact replay or a better nominal result. The existing main text appropriately does not substitute it for the selected condition.

## Sampling distinction that must be preserved

The portable 47.23 mm lateral-gap value is a minimum calculated from every fifth 10 ms logged pose (50 ms geometric sampling). Evaluating every 10 ms pose in the same final window gives 46.9349 mm. Both are positive; the sampling difference does not change the accepted gap condition. The existing chapters do not state this gap number, so no prose replacement is needed. If a future revision quotes either value, its geometric sampling must accompany it. This audit does not substitute the denser value for an adopted statistic without disclosure.

## Figure recommendation and provenance

`CAD_supported_traces_20261116.pdf` is a new vector figure generated from all 1,200 final-rung samples. It uses native flexion-negative knee coordinates and one common 16–28 s horizontal axis. Its four panels show bilateral RG flexor-minus-extensor voltage, knee angle, summed normal foot load, and minimum mesh height. It exposes both unloaded clearance and the soft-contact penetrations, making the acceptance and limitations more assessable than the static rendering alone.

The figure is exactly 14 × 15 cm with 10 pt text throughout. Its height was compacted before active placement to reserve caption space in the 19.28 cm manuscript text area. No smoothing, phase shift, amplitude adjustment, or decimation is used. Only stated differences, foot-region summation, and unit conversions are applied. Matplotlib path simplification is disabled. PDF, SVG, PNG, source hashes, transformations, windows, and output hashes are recorded in `CAD_FIGURE_PROVENANCE.json`. The PNG was rendered and visually inspected; panel labels, legend, axes, and threshold lines are readable. Final manuscript scaling remains the coordinating reviewer's responsibility.

Suggested caption:

> Aligned CAD supported-stepping signals during the 16–28 s final-support window: bilateral RG flexor-minus-extensor voltage (A), native flexion-negative knee coordinates (B), summed foot-normal loads (C), and whole-foot minimum mesh heights (D). Blue solid and orange dashed lines identify right and left sides. Dotted thresholds show the −1° knee guard, 7.36 N unloading level, and 5 mm clearance level; the dash-dot line in C marks the 15 N contralateral-load criterion. Qualifying intervals also require the simultaneous opposing load condition and at least 80 ms duration. The 75 kg model receives 77.5% BW mean upward assistance in this window, with horizontal and attitude guidance retained. The complete-trial per-physics-step knee bound is evaluated separately from these 10 ms traces. Actuation uses ideal reference tensions and engineering allocation.

The root agent owns active placement and TeX changes. No existing active asset has been moved or overwritten by this audit.

## Numerical and interpretation ledger for this review

| Existing statement | Disposition | Reason |
|---|---|---|
| 77.5% BW over 16–28 s; 7R/6L; 12 alternations | Retain unchanged | Direct recomputation from raw load, force, and mesh arrays |
| Minimum knee flexion 5.10° | Retain; clarify left side and whole-integration online accumulator if expanded | 2 ms accumulator is distinct from final-rung 10 ms display |
| Reference/fresh exact agreement | Retain unchanged | All 26 saved arrays compare exactly in both pairs |
| Five of eight local cases pass | Retain unchanged | Raw interval counts and retained completion outcomes agree |
| CAD calibrated BPAs or ordinary neural-only control | Continue to exclude | Force law and MN compensation explicitly use ideal tensions and engineering targets |
| 16.13/12.91 mm transient penetrations | Retain unchanged | Direct raw final-window mesh-height minima |
| 47.23 mm lateral gap, if newly quoted | Do not add without its 50 ms sampling definition | All 10 ms poses give 46.93 mm; no existing prose value replaced |
| New trace figure | Add as quantitative evidence if useful | Reveals simultaneous timing, loading, clearance, native knee bounds, and soft-contact limitation |

## Verification limits and errors

No historical compiled MJB or compiled-array digest was identified for this CAD package. Retained source hashes, model XML, exact saved repeat arrays, and direct current-geometry readback establish the documented checks; they do not retroactively create a missing historical compiled-binary hash. Full 2 ms trajectories were not saved, so the every-step extremum is verified as a source-bound online accumulator rather than independently recalculated from omitted samples.

The first audit script attempt failed on importing SciPy because `D:\Anaconda\envs\myo\lib\site-packages\scipy\.libs` raised WinError 5. The final reproducible audit removes that dependency and evaluates linear mesh extrema over all vertices, which gives the same extrema as their convex hull. No software environment was modified. RG period values are explicitly identified as retained source statistics, not freshly recalculated by this script.

No manuscript, model, controller parameter, actuator law, route, physical coefficient, optimizer result, or simulation recording was changed. No simulation, optimization, retraining, commit, push, or external upload was performed.
