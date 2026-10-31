# Canonical 108-neuron W2L contact routing

Session 1, **easteregg2**, 31 October 2026. Read-only audit of the actual canonical runtime, source SHA256 `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`: 108 neurons and 216 chemical connections. No GUI, physical, or controller parameter was changed.

## Main finding

Contact feedback is distributed across **16 direct chemical connections**. Each heel relay drives **seven target populations simultaneously**, including RG, PF, and MN levels. The heel-to-ankle-extensor MN maximum conductance is **50 times** that of the other contact routes. All sixteen paths use contact **level**, not a heel-strike edge detector or a finite reset pulse. This was a candidate for sustained stance reinforcement. The subsequent saved-ground geometry check below finds that the heel inputs are inactive in the inspected baseline window because the model is physically toe-supported. That evidence must guide the ablations; structural fan-out alone does not establish the failure's cause.

### Direct routes, mirrored on both sides

| Physical relay → neural destination | Maximum conductance (uS) | Reversal (mV) | Presynaptic threshold / saturation (mV) | Named type |
|---|---:|---:|---|---|
| Heel → opposite RG flexor HC | 0.01 | −40 | −60 / −40 | `Heel_Contra_RGflx_Excite` |
| Heel → same-side RG extensor HC | 0.01 | −40 | −60 / −40 | `Heel_Ipsi_RGext_Excite` |
| Heel → same-side hip PF extensor | 0.01 | −40 | −60 / −40 | `Heel_Ipsi_HipPFext_Excite` |
| Heel → same-side knee/ankle PF extensor | 0.01 | −40 | −60 / −40 | `Heel_Ipsi_KneePFext_Excite` |
| Heel → same-side hip extensor MN | 0.01 | −40 | −60 / −40 | `Heel_Ipsi_HipMNext_Excite` |
| Heel → same-side knee extensor MN | 0.01 | −40 | −60 / −40 | `Heel_Ipsi_KneeMNext_Excite` |
| Heel → same-side ankle extensor MN | **0.5** | −40 | −60 / −40 | `Heel_Ipsi_AnkMNext_Excite` |
| Toe → same-side ankle extensor MN | 0.01 | −40 | −60 / −40 | `Toe_Ipsi_AnkMNext_Excite` |

These are **direct inputs to the named HC, PF, or MN cells**, not inputs first routed through an inhibitory interneuron. The −40 mV reversal is depolarizing over the usual target-voltage range in these trials; the actual current sign remains `g × (Erev − Vpost)`.

**The ankle MN labels are not their verified mechanical actions.** The route-matched canonical mapping in `runs\canonical_12II\BUILD_LEDGER.json`, combined with the independently verified positive-dorsiflexion coordinate, gives:

| Authored MN label | Actual target muscle name | Verified mechanical action | Moment arm |
|---|---|---|---:|
| `L/R Ank MN ext` | `ankle_L/R_flx` | Dorsiflexion | +0.08 m |
| `L/R Ank MN flx` | `ankle_L/R_ext` | Plantarflexion | −0.04 m |

Consequently, the strong heel-to-`Ank MN ext` path excites the **dorsiflexor**, not the plantarflexor. The unchanged same-label knee-PF-flexor → ankle-MN-flexor connection drives plantarflexion during its active pattern, so its phase relative to knee-flexion swing should be tested explicitly. The existing ankle-PF-crossed option is a diagnostic adaptation of this mechanical binding; it must not be mislabeled as an unchanged reproduction of the published rat model. No phase swap is applied by this audit.

The direct chemical endpoint classes, type parameters, four conversion coefficients, and contact-relay intrinsic parameters are bilaterally matched. This establishes contact-path symmetry, not symmetry of every downstream mechanical or neural response.

## Physical sensor → conversion → graded contact relay

Four distinct source boxes, `foot_L_contact`, `toe_L_contact`, `foot_R_contact`, and `toe_R_contact`, feed the corresponding labelled heel/toe relays. The source heel-labelled boxes are 60 × 80 × 2 mm; toe boxes are 15 × 2 × 80 mm in their serialized local axes. Their parent transforms determine world orientation.

Each PhysicalToNode adapter uses `ContactCount → ExternalCurrent`, polynomial **A=0, B=0, C=1e−9 A/count, D=0**, with gain limits disabled. A binary contact value of one therefore supplies **1 nA**, and zero supplies zero. The relay has rest −60 mV, membrane time constant 10 ms, relative size 1, zero tonic current, and no plateau conductance.

For an isolated continuously active relay, the port's equations approach approximately −59 mV. At that presynaptic voltage, the −60/−40 mV synaptic ramp activates 5% of maximum conductance: approximately 0.025 uS for the heel-to-ankle MN path versus 0.0005 uS for a 0.01 uS path. At a hypothetical postsynaptic −60 mV, those currents are approximately 0.5 and 0.01 nA, respectively. These are equation-based steady-state examples, not measured ground-trial currents.

## Contact level is not a load or strike event

The current MuJoCo runner sets each source contact value to one when the transformed ghost-box bottom is at/below the belt top and its horizontal center lies within the belt extent plus a 100 mm margin. Otherwise it sets zero. It does **not** test a normal-force threshold, differentiate the signal, impose a pulse width, or phase-gate its delivery. The graded relay then filters that level over its membrane time constant.

The actual heel/toe box signal and actual foot vertical forces are separately recorded. A long contact signal can therefore reflect geometric overlap; it should not automatically be read as continuing stance load. Likewise, counting MuJoCo contact manifold points would not be equivalent to this binary source-body contact proxy.

Implementation anchors: `run_virtual_walker.py:89–98` and `native_kernel.py:89–125` in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\`.

## How the signal spreads through interneurons

The direct contact injections act in parallel with the ordinary hierarchy:

- Same-side heel → RG extensor → its RG inhibitory interneuron → opposing RG flexor. RG extensor also drives hip and knee/ankle PF extensor populations.
- Same-side heel → hip or knee/ankle PF extensor → its PF inhibitory interneuron → opposing PF flexor. PF outputs also reach MNs and Ia-related cells.
- Heel/toe → ankle extensor MN → its recurrent Renshaw cell → the recurrent motor/Ia circuit.
- Heel → opposite RG flexor → its RG inhibitory interneuron, commissural c1 pathway, and flexor PF outputs.

Consequently, an extensor MN can receive contact effects directly and indirectly through PF/RG pathways. The audit JSON retains exact first- and second-edge definitions; a path's qualitative net sign should not replace its conductances, reversals, state, or timing. No phase-dependent Ib gate was found in these direct contact connections. Ia/II/Ib receptor pathways are separate and should remain active in contact-only ablations.

## Resolved ablation classes

Exact bilateral connection IDs and parameters are in:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\canonical_12II\CONTACT_ROUTING_AUDIT.json`

The `ablation_edge_groups` dictionary includes each table row separately, plus:

- `heel_ipsilateral_RG_PF_HipKnee_MN_cascade`: ten same-side heel routes, excluding the separate ankle-MN and opposite-RG routes.
- `all_contact_to_MN`: eight direct motor routes.
- `all_contact_to_central_RG_PF`: eight direct central routes.

A useful diagnostic sequence is intact contact, heel ipsilateral-cascade removed, direct contact-to-MN removed, strong heel-to-ankle-MN removed, and all contact routes removed while retaining muscle feedback. Hold the controller, guide geometry, starting state, clock, and non-contact gains fixed. A second experiment can compare continuous level with an explicitly defined heel-strike pulse; that is a changed signal protocol and should be recorded separately from a conductance ablation.

Compare actual stance loading, swing-foot unloading and clearance, bilateral timing, knee limits, and measured support over a common window. Rhythm alone cannot identify which contact path affects loading. No ablation result is asserted by this structural audit.

## Historical connection-name convention supplied by Ben

Ben decoded the authored convention as `{destination}_{Source}{Configuration}{type}{exc/inh}`. His example `RG_HF2i` denotes a connection **to RG**, **from hip flexor**, with **type-II feedback**, and **inhibitory** action. This is a user-supplied key for interpreting the historical labels, not a renaming rule or a substitute for resolving the actual source/destination IDs and reversal potential. Preserve authored names, and check the serialized endpoints and effective type fields. In particular, a central phase-related II route and a newer homonymous II-to-MN excitatory route are different pathways even when they share a muscle/source signal.

## Measured baseline: toe support and restricted dorsiflexion

Forward kinematics of the saved `canonical108_ground_baseline_v1` states over 4–8 s used **orientation-aware box extents**, `abs(R[Y,:])·half_size`. No controller or integration run was added.

- Both heel-box contact duties are **zero**. Toe-box duties are 0.6652 left and 0.0450 right.
- Heel-box bottom gaps are 11.57–100.13 mm left and 10.90–163.04 mm right. The main foot collision boxes remain above the belt, with minimum gaps of 7.91 mm left and 2.13 mm right, while the toes carry contact forces.
- Under recorded foot loading above 20 N, median CoP lies approximately **104.8 mm forward of the main-foot geometry center** on both sides. That is beyond the main foot's +60 mm front edge and in the toe region.
- The heel boxes extend below their corresponding physical sole planes, rather than being mounted too high. The zero heel signals are consistent with the toe-down posture. Existing 2 ms state/logging alignment does not explain the centimeter-scale separation.

Independent geometric calibration confirms **positive raw ankle rotation is dorsiflexion**, with approximately zero anatomical offset: adding +1 raw degree raises geometric dorsiflexion by +0.9999999 degree. Both nominal ankle limits are approximately **−12 degrees plantarflexion to +3 degrees dorsiflexion**. Soft limits permit transient overshoot; the observed 4–8 s maxima are +5.29 degrees left and +6.90 degrees right, not a wider nominal range.

| Time (s) | Side | Ankle DF (deg) | World foot pitch, toe-up positive (deg) | Heel-box gap (mm) | Recorded vertical foot force (N) |
|---:|---|---:|---:|---:|---:|
| 4.0 | L | 1.384 | −20.554 | 23.69 | 20.24 |
| 4.0 | R | 3.119 | −17.434 | 13.67 | 168.00 |
| 4.5 | L | 3.009 | −11.533 | 13.98 | 0 |
| 4.5 | R | −10.896 | −37.518 | 56.00 | 0 |
| 5.0 | L | 1.718 | −18.902 | 22.24 | 0 |
| 5.0 | R | −0.403 | −35.332 | 45.83 | 48.57 |

Holding all other coordinates fixed at 4.0 s would require approximately **21.94 degrees left / 20.55 degrees right DF** to make the main-foot longitudinal axes horizontal. This is a diagnostic of that pose, not a desired standing posture or a recommendation to replace hip/knee/root alignment with ankle motion alone.

The restricted range is **not a MuJoCo converter clamp**. The historical 2021 project at git `23b0e785` already has enabled raw ankle limits −20 to −5 degrees. Its assembled anatomical home is approximately +8.00021 degrees, giving physical limits −12 to +3 degrees. The later geometry builder changes ankle home to zero and shifts the relative limits by +8.00021 degrees, preserving that physical range. It also chooses the initial hip/knee poses that yield foot pitches −20 degrees left and −35 degrees right. Thus, the prior geometry reconstruction changed the starting posture; its comment claiming to keep the feet level is incorrect. The later virtual-walker/guide construction and MuJoCo importer retain these limits.

Detailed evidence in `runs\canonical_12II\`: `HEEL_GEOMETRY_FK_AUDIT.json`, `ANKLE_DORSIFLEXION_FK_AUDIT.json`, `ANKLE_LIMIT_PROVENANCE.json`, and `HISTORICAL_ANKLE_LIMIT_PROVENANCE.json`. These identify the exact source paths, hashes, joint IDs, reference changes, loaded CoP, and sample poses. No physical parameters were changed in this audit.
