# Transferring a reduced locomotor controller to Gait2392

Session1, easteregg2. Design rationale and prerequisites, not a new simulation result.

## What the23coordinates mean

Six coordinates locate and orient the pelvis in space: three translations and three rotations. These are called floating-base coordinates. They are not six joints that connect body segments, and internal muscles do not directly actuate the model's overall translations. Ground reactions and external support participate in changing whole-body motion.

The remaining17 describe joints relative to adjacent segments: seven per leg (three hip rotations, knee flexion, ankle flexion, subtalar inversion/eversion, and toe flexion) and three lumbar rotations. This distinction prevents confusing an anatomically23-coordinate model with a controller that independently stabilizes23motions. Converted moving-path-point follower coordinates are implementation coordinates, not additional anatomical joints.

## Why a small controller can still be useful

A rhythmic task does not require an independent oscillator for every joint. A low-dimensional phase/pattern generator can coordinate many muscles, and feedback, passive mechanics, and ground contact can stabilize deviations. What must be demonstrated is that the assembled system supplies the task's required control, including motions outside the source controller's original sagittal subspace.

The justified transfer is therefore modular: retain an identified rhythm/PF organization, remap its outputs to the target muscle mechanics, add explicitly identified balance/posture/foot-placement pathways where needed, and progressively remove external constraints. Muscle-name matching alone does not establish matching torque signs, moment arms, joint ranges, force capacity, or feedback gains.

## Mathematical statement

For anatomical generalized coordinates q, write

    M(q) qdd + h(q, qdot) = R(q) F(a, length, velocity) + Jc(q)^T lambda_c + tau_harness.

R is the signed muscle moment-arm matrix, contact forces enter through Jc, and the actual external support is tau_harness. The muscle-force law and passive effects must be stated; they can differ between biological muscles, compliant-tendon muscles, and BPAs.

A useful design description is

    u = saturation(W p(phi) + delta_u_feedback(q, qdot, muscle/contact signals)).

This is a design decomposition, not a claim that every current controller literally computes a linear sum in excitation space. Conductance-based SNS circuits implement nonlinear combinations and must be checked against their actual compiled equations. For syn6, the existing fixed nonnegative basis provides the candidate W, while actual RG/PF dynamics produce its time-varying inputs.

At one posture, the local torque authority of that fixed synergy component is approximately

    B_synergy = R(q) (partial F / partial a) (partial a / partial u) W.

Its rank cannot exceed the number of independent synergy inputs. This limits independent instantaneous torque commands, but it is not by itself a test of walking stability. Periodic motion is a coordinated task; contact, passive dynamics, changing posture, and feedback matter. If all corrective commands are forced through the same fixed W, they inherit the same instantaneous span limitation. Additional feedback must actually reach the muscles/directions needed for the added task.

Six primitive channels per leg mean up to12bilateral primitive inputs, not six whole-body inputs. Reconstruction accuracy of muscle activation is different from sufficient torque authority and different from closed-loop stability. Rank/SVD is a useful local diagnostic after scaling torque/coordinate units; it does not replace a constrained stability analysis or perturbation tests.

## Required coverage and staged tests

For each controller/configuration, record each anatomical coordinate, its actual neural/actuator pathways, passive restoring/damping action, harness restraint, and measured task error. State which axes remain free. A spring guide is not a rigid removal of a degree of freedom; its force/torque and deflection must be reported.

1. Establish sagittal coordination under explicitly stated constraints, with actual alternating loading and foot clearance.
2. Add or verify lateral pelvis/trunk balance, hip abduction/adduction, rotation, and appropriate foot control. Use moment arms and signal domains at representative postures, not anatomical names alone.
3. Compare the same saved controller with each added pathway disabled and with restraint reduced separately from vertical body-weight support.
4. Measure out-of-plane angles/velocities, harness forces and torques, foot placement/width, loading, clearance, and perturbation recovery. A small vertical support percentage alone does not establish free3D balance.

For locally linearized stabilization around a posture, check the relevant constrained dynamics. For a periodic gait, use appropriate periodic-orbit/perturbation tests; an equilibrium transfer-function calculation alone does not establish gait stability.

## Literature precedent and defense answer

Song and Geyer explicitly extended a sagittal reflex model to3D by adding lateral hip mechanics, abductor/adductor control, and lateral trunk balance; their later model included functional feedback modules and foot-placement control. This is a close precedent for a modular extension, rather than direct reuse of an unchanged reduced controller.

Their2015human model has eight internal degrees of freedom, so this precedent does not by itself supply all17internal-coordinate pathways needed for our Gait2392 implementation.

- Song and Geyer2013, *Generalization of a Muscle-Reflex Control Model to3D Walking*: https://publications.ri.cmu.edu/storage/publications/2017/11/SongGeyer13cIEEE_EMBC.pdf
- Song and Geyer2015, *A neural circuitry that emphasizes spinal feedback generates diverse behaviours of human locomotion*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4560581/

Suggested spoken answer: "The controller's source model was reduced, but we transferred its modular organization rather than assuming its original gains solve the larger problem. We remap the muscle outputs, scale the signal and actuator domains, identify the motions handled by added feedback, and measure the stabilization supplied by the support apparatus. The present supported trials demonstrate the coordination achieved under those conditions; reducing restraint is the test of the added balance pathways."

## Current action

Control-coverage verification is now required for every family in CAMPAIGN_GATE.json before the next Hunt/curriculum campaign. Existing reported trials remain unchanged. No new3D-walking result is asserted by this design note.

### Static check of the existing syn6 basis

The read-only diagnostic `syn6_torque_span.py` evaluated four postures (0%,25%,50%, and75%) from the existing walking-angle cache. It used the actual native MuJoCo active-force derivative, equality-aware anatomical moment arms, and the unchanged K6 NMF basis. At every posture, the full92-muscle map had numerical rank17; the86muscles represented by the basis also had rank17; projecting through the12bilateral primitive inputs gave rank12.

This identifies a real reduction in independent instantaneous torque directions for that fixed activation basis. It does not mean the complete syn6 controller has only12control inputs: its additional reflex and posture pathways must be evaluated separately, and its conductance-based input/output function is not literally the NMF reconstruction. It also does not mean17independent primitive channels are necessary for a coordinated gait. Contact, passive dynamics, feedback, and the requested task determine what is sufficient.

Full numerical record: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\SYN6_TORQUE_SPAN.json`. This was a static design calculation, without time integration, optimization, training, NMF refitting, or source-model changes.
