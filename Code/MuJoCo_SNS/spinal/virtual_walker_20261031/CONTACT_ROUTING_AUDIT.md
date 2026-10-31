# Canonical 108-neuron W2L contact routing

Session 1, **easteregg2**, 31 October 2026. Read-only audit of the actual canonical runtime, source SHA256 `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`: 108 neurons and 216 chemical connections. No GUI, physical, or controller parameter was changed.

## Main finding

Contact feedback is distributed across **16 direct chemical connections**. Each heel relay drives **seven target populations simultaneously**, including RG, PF, and MN levels. The heel-to-ankle-extensor MN maximum conductance is **50 times** that of the other contact routes. All sixteen paths use contact **level**, not a heel-strike edge detector or a finite reset pulse. This is a concrete candidate for sustained stance reinforcement; ground ablations are needed to establish which paths prevent unloading.

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
