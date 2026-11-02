# Defense controller mathematics and calibration decisions

Gait2392, easteregg2, 2 November 2026. This note describes the current Hunt-PT1033 four-cell RF/hamstring candidate and its version-2 integrated FSA profile. It does not describe every historical S3K, W2Lvar, or Syn6 experiment. The electrical proofs below precede the bounded Hunt-D fit; they do not imply that the fitted controller has passed a body stage. No tendon-compliance or joint/muscle-damping variant is enabled by this neural profile.

## 1. What equations actually run?

Let \(U_i=V_i-V_{\mathrm{rest},i}\) denote a neuron's relative-rest voltage. For the zero-rest cells used here, the conductance equation is

\[
C_i\dot U_i=-G_iU_i+I_{b,i}+I_{\mathrm{app},i}
+\sum_j\bar g_{ij}s_{ij}(U_j)(E_{ij}-U_i)+I_{\mathrm{NaP},i},
\qquad
s_{ij}(U)=\operatorname{clip}\left(\frac{U-E_{lo,ij}}{E_{hi,ij}-E_{lo,ij}},0,1\right).
\]

The first index is postsynaptic. Conductance \(\bar g\) is nonnegative; the reversal and instantaneous driving voltage determine current direction. In this controller, common reversals are +8 mV and −5 mV. A nominally excitatory +8 mV synapse carries outward current if its postsynaptic voltage exceeds +8 mV. Its conductance still shunts other inputs. The 0–5 mV activation window is therefore not a bound on every NaP neuron's voltage.

The implementation uses seconds, mV, nA, µS, and backend capacitance in µF: \(1\,\mu\mathrm S\times1\,\mathrm{mV}=1\,\mathrm{nA}\), and \(C/G\) is seconds. Thus backend \(C=0.03\) µF at \(G=1\) µS means 30 nF and a 30 ms passive membrane time constant. This explicit adapter matters when comparing libraries that use milliseconds and nF.

For the RG's persistent-sodium mechanism,

\[
I_{\mathrm{NaP}}=\bar g_{\mathrm{NaP}}m_\infty(U)h(E_{\mathrm{NaP}}-U),\quad
m_\infty(U)=\frac{1}{1+k_m e^{s_m(e_m-U)}},\quad
\dot h=\frac{h_\infty(U)-h}{\tau_h},\quad
h_\infty(U)=\frac{1}{1+k_h e^{s_h(e_h-U)}}.
\]

The actual saved Hunt-A record supplies \(\bar g_{\mathrm{NaP}}=12\) µS, \(E_{\mathrm{NaP}}=8\) mV, \(k_m=k_h=1\), \(s_m=0.8\,\mathrm{mV}^{-1}\), \(e_m=0\) mV, \(s_h=-2\,\mathrm{mV}^{-1}\), \(e_h=2\) mV, and **fixed** \(\tau_h=0.6313716989889969\) s. These are the fitted source values, not the raw defaults in `params.py`. The fixed-tau backend first advances \(h\) from the old voltage, then uses that updated gate in the explicit voltage update. Replacing fixed \(\tau_h\) with a voltage-dependent gate time constant changes the model.

The compiled sparse implementation retains this neuronal update. Sparse execution and preallocation reduce computation time; they do not shorten biological time constants. For finite admissible states, an edge with \(\bar g=0\) contributes neither current nor shunting. The composed-network zero/delete/rebuild and dense-reference checks agree exactly over their 64-step fixture; this is stronger than checking a diagram's zero label.

Source: [fixed-tau backend](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/build_network.py), [Hunt builder](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/controller.py), and [saved Hunt-A fit](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/runs/step_A/winner.json). These fitted ionic parameters are not a claim of literal reproduction of another paper's membrane parameters.

## 2. Functional targets determine conductances, including shunting

The graded [functional-subnetwork approach of Szczecinski, Hunt, and Quinn (2017)](https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2017.00037/full) motivates specifying a desired transformation before choosing raw weights. The following pool calculation is the implemented engineering endpoint contract, with the actual controller's reversals and activation windows; it is not an assertion that conductance networks are globally linear arithmetic.

Consider one isolated terminal pool. Hold an MN at baseline \(U_0\) with \(I_0=G U_0\), activate the selected sources to measured isolated voltages, and seek \(U_t=U_0+\delta\). Preserve relative compiled weights \(g_j^0\) by setting \(g_j=\lambda g_j^0\). With \(\alpha_j=s_j(U_j)\), equilibrium gives

\[
0=-G\delta+\lambda\sum_jg_j^0\alpha_j(E_j-U_t),\qquad
\boxed{\lambda=\frac{G\delta}{\sum_jg_j^0\alpha_j(E_j-U_t)}}.
\]

The denominator contains the **target** voltage, so the solution accounts for voltage-dependent shunting. It requires a feasible sign and positive finite scale. It preserves meaningful relative allocation within the pool instead of blindly dividing every edge by its fan-in count. If other active inputs establish the baseline, their effective conductance must also be included; the isolated proof deliberately excludes them.

**Worked applied example: `IIX_iliacus_r → MN_iliacus_r`.** The ledger records a 1 nA calibration input to `II_iliacus_r`, producing 1 mV there. The actual upstream relay then produces \(U_{IIX}=1.3333333333333335\) mV: its active conductance is 0.2 µS, with a 1 µS leak and +8 mV reversal. The terminal window is 0–5 mV, hence \(\alpha=0.2666666666666667\). For \(G=1\) µS, \(U_0=2.5\) mV, and \(\delta=+0.5\) mV,

\[
g=\frac{1\times0.5}{0.2666666666666667\,(8-3)}
=0.375\;\mu\mathrm S,
\qquad U_\infty=\frac{2.5+(0.375)(0.2666667)(8)}{1+(0.375)(0.2666667)}=3\;\mathrm{mV}.
\]

The applied ledger value is 0.37499999999999994 µS, obtained by scaling the restored template's 0.4 µS by 0.9374999999999998. That 0.4 µS is the template value before FSA normalization, not proof that the saved body controller previously delivered an active II current of that size. Changing only the II target to +0.02 mV gives \(g\approx0.01368613\) µS by the same equation. Neither value is a measured human synaptic conductance. The 1 nA pulse and the original ±0.5 mV targets are declared calibration choices, not a physical maximum or physiological estimate.

The four pools are homonymous Ia excitation, homonymous II excitation through IIX, autogenic Ib inhibition through IBIN, and pooled reciprocal Ia inhibition through IaIN. Distinct interneurons, reversals, source voltages, and pool membership remain explicit. The simultaneous circuit response is not the sum of the four isolated voltage changes.

An additional optional **Ib/phase coincidence** repair is separate from these four targets. Its illustrative normalized budgets are each \(0.5G\): with +8 mV reversal and zero baseline, either input alone gives \(U=8(0.5)/1.5=2.6667\) mV, below the 2.7 mV output threshold; both give 4 mV. Actual incoming relative load weights and windows are derived from the compiled cells. In the selected profile, the load-high reference is **0.225 nA**, from \(0.5\times1.5\times0.3\times F/F_{max}\) at \(F/F_{max}=1\). It is muscle-load encoding, not a ground-reaction force or universal input bound. Electrical coincidence does not establish appropriate whole-body stance recruitment.

## 3. Time scaling follows a task and its clock, not body mass alone

The design context is the shortest usable **cached** running-stride estimate used in this campaign, \(T=0.559\) s: \(f_1=1/T=1.788909\) Hz. The cached ankle waveform's fourth harmonic supplies the selected 99%-AC-energy context, \(f_d=4/T=7.155635\) Hz. This is neither a universal fastest human stride nor a measured neural bandwidth; it does not establish the original trial's acquisition sampling rate or the bandwidth of every sensor signal.

For a passive encoder \(H(s)=1/(1+s\tau)\), a declared **per-encoder** phase budget \(\beta=10^\circ\) implies

\[
\tau\le\frac{\tan\beta}{2\pi f_d}=3.921848\;\mathrm{ms}.
\]

The actual explicit update runs at unchanged \(\Delta t=2\) ms. If \(a=1-\Delta t/\tau\), input is sampled at the start, and output is timestamped at the end of the update,

\[
H_d(e^{j\theta})=\frac{(1-a)e^{-j\theta}}{1-ae^{-j\theta}},\qquad \theta=2\pi f\Delta t.
\]

Solving the **discrete** phase condition gives the consumed \(\tau=3.896228836665\) ms. The measured encoder lag is 10° at \(f_d\), and 2.508619° at \(f_1\). Additional stale-input age of 0–2 ms adds up to 5.152057° or 1.288014°, respectively; it is not hidden in the 10° allocation. Nor does that allocation cover another relay, the MN, activation dynamics, a delayed balance pathway, or mechanical response. Only eligible passive sensory encoders receive this capacitance; RG NaP gates, PF shape/phase, MN membranes, and muscle dynamics retain their source timing.

Size and speed affect mechanics jointly, as discussed in [Sutton et al.'s 2021 mechanics-space preprint](https://doi.org/10.21203/rs.3.rs-153297/v1). For a local mechanical approximation,

\[
M\ddot q+D\dot q+Kq=\tau,\qquad
\text{relative contributions at }\omega:\;M\omega^2,\ D\omega,\ K.
\]

Even for a geometrically similar gravity-dominated limb, \(J\sim\rho L^5\), \(K_g\sim\rho gL^4\), so \(J\omega^2/K_g\sim L\omega^2/g\); mass alone does not determine the regime or a neuronal time scale. The actual Gait audit instead projects native inertia and forces through the equality map \(E\): \(M_a=E^TM_{native}E\). It derives local damping and restoring-force gradients at declared poses, including apparatus contributions. Negative restoring stiffness and coupled modes cannot be replaced by a positive scalar mass/time rule. That static audit changes no physics.

## 4. Three neural abstractions need different FSA treatments

| Model | Appropriate design/calibration | What is not interchangeable |
| --- | --- | --- |
| Graded network with NaP RG | Preserve ionic states and intrinsic rhythm; derive selected relay bandwidth and conductance-pool functions at declared operating points, then test the coupled system. | Reducing NaP inactivation time is a rhythm-model change, not a free transmission speedup. |
| Pure GLIF carrier | Specify voltage/rate domains, bias, threshold adaptation, membrane and synaptic timing, event law, and a rate-decoding contract. Validate rate and transient approximations. | A threshold-adaptation variable is not the NaP \(h\) gate. Reset-to-peak synapses are not additive-and-capped conductance increments. |
| Analytic phase oscillator and primitives | Specify phase rate/coupling/reset, primitive width/center, and signed recruitment; assess phase and downstream bandwidth. | A raised-cosine primitive is not a passive membrane or a NaP half-center, even if its output resembles a burst. |

The [2020 GLIF extension](https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2020.577804/full) is implemented separately in `table2_glif.py`. Its high-rate construction uses \(\theta^*\approx\theta_0/(1-m/2)\), \(\tau_{mem}=R/(F_{max}\theta^*)\), and, for nonzero threshold adaptation, \(\tau_\theta=\tau_{graded}(1-m/2)\). The public adapter uses \(C_{nF}=1000G_{\mu S}\tau_s\). Those are scoped approximation/design rules, not proof of finite-rate equivalence or permission to overwrite ionic bursters. The current Hunt1033 candidate remains graded; these formulas do not imply that it was converted to spikes.

[Di Russo et al. (2023), Eqs. 6–9](https://infoscience.epfl.ch/server/api/core/bitstreams/4e04db10-e21e-4aea-b772-4a50bf1585aa/content) instead use coupled phase oscillators, contact phase resetting, and five raised-cosine primitives. A primitive has the form \(p_k(\bar\phi)=\tfrac12[1+\cos(\pi(\bar\phi-\mu_k)/\sigma_k)]\) inside its phase window and zero outside; MN input sums signed weighted primitives. Their local reflex implementation is another layer. Production Syn6 is an engineered NMF/neuronal hybrid, not that exact five-primitive controller. The current Hunt candidate has an explicit RF/hamstring four-cell PF circuit per leg; sartorius and gracilis remain separate. Component resemblance is not identity of the complete models.

## 5. What transfer functions and Jacobians can establish

Within a fixed activation-window branch, a passive neuron's small perturbation obeys

\[
C_i\delta\dot U_i=-\left(G_i+\sum_j\bar g_{ij}s^*_{ij}\right)\delta U_i
+\sum_j\bar g_{ij}s'_{ij}(U_j^*)(E_{ij}-U_i^*)\delta U_j+\delta I_i.
\]

Thus \(\tau_{eff}=C_i/(G_i+\sum_j\bar g_{ij}s^*_{ij})\), and the presynaptic gain depends on the postsynaptic operating voltage. NaP adds voltage/gate derivatives and its gate state to the Jacobian. At a window corner, a single smooth derivative is insufficient. A local constant-state model \(\delta\dot x=A\delta x+B\delta u\), \(\delta y=C_y\delta x+D_y\delta u\), gives \(H(s)=C_y(sI-A)^{-1}B+D_y\) only within that operating-point and mode approximation.

A walking state also includes joint position/velocity, muscle activation, relevant filters/delays, neuronal voltage/gates, and contact/protocol state. Along a periodic trajectory it has a **time-varying** linearization \(A(t)\). Return-map/Floquet analysis requires propagation over the stride and the appropriate event/reset sensitivities. Soft contact and the native constraint solver must be differentiated as actually implemented; an ideal hard-impact reset must not be substituted without justification. A stable isolated encoder or a full-rank static muscle moment-arm matrix does not certify that return map.

For the executed code, the corresponding object is the discrete map \(x_{n+1}=\Psi_{\Delta t}(x_n,u_n)\), with \(J_n=\partial\Psi_{\Delta t}/\partial x\) and \(B_n=\partial\Psi_{\Delta t}/\partial u\). At a constant operating point its transfer is \(H_d(z)=C_y(zI-J)^{-1}B+D_y\). One must not substitute \(J=I+\Delta t A\) without checking the actual update order, filters, and solver. The updated NaP gate and the separate neural/physics timestamps are concrete reasons to differentiate the actual map.

Gait2392 has 6 floating-base plus 17 internal anatomical coordinates. Its 105 software coordinates include dependent followers; they are not 105 anatomical DoFs. At the four audited reference poses, the equality-aware native 92-muscle active-force map has rank 17, while the frozen bilateral K6 activation-factor map has rank 12. These are local activation-span results, not the compiled SNS Jacobian or a proof of free 3D balance. Harness restrictions, passive structures, and neural control coverage remain declared separately.

## 6. Why Hunt-D is the current bounded fit

Hunt's four-part training method is attributed here to **A. J. Hunt, “Real-time control of dynamic walking in a dog-like robot,” 2nd Intel Neuromorphic Research Community Meeting, 15 February 2019**, as identified by Ben. In this implementation, A addresses intrinsic rhythm, B expected sensory timing, C feedforward motor recruitment, and D sensory refinement. The later physical curriculum is a separate validation sequence.

With the **same actual saved configuration** (including tonic RG current 1.2 nA), the source controller passes air stages 1 and 2. The new function-equivalent-timing design preserves the stage-1 assessment exactly, but fails stage 2's knee-phase criterion:

| Verification | Knee cycles R/L | Period s | Knee phase error | Period CV | Accepted |
| --- | ---: | ---: | ---: | ---: | --- |
| Source stage 1 | 5/4 | 1.794 | 0.00735 | 0.00365 | Yes |
| Source stage 2 | 5/4 | 1.784 | 0.00562 | 0.00409 | Yes |
| New function-equivalent stage 2 | 5/5 | 1.781 | 0.20895 | 0.24106 | No |

The new stage-2 failure has no recorded crossing or knee hyperextension. Its RG and PF clocks remain regular: RG-E-right onset CV is approximately 0.00056, versus 0.00049 in the source. Knee flexion peaks become irregular downstream. Current reconstruction from the saved full neuronal states matches recorded local-current groups to at most \(3.47\times10^{-18}\) nA. In the 5–14 s right-side RG-E-dominant samples, the new psoas II contribution averages **0.88454 nA**, versus zero in the source, and mean psoas activation rises **0.05567 → 0.19074**. New psoas reciprocal inhibition during RG-F dominance averages **−0.48399 nA**. These are measured/reconstructed currents in the air verification, not claims about normal human stance phases.

Stage 1 deliberately masks timing afferents, sensory reflex outputs, and balance. Its bit-identical result is therefore a useful preservation check, not validation of the new closed-loop feedback. The full network has 1033 neurons; that deafferented protocol omits eight PT neurons and builds 1025. Stage 2 is the first of these two tests to exercise the newly composed feedback.

This evidence selects **Hunt-D local recruitment calibration**, rather than changing an already regular RG clock. It does not yet prove that one isolated pathway alone causes the body failure. The independent coincidence change and simultaneous-input interactions remain possible contributors.

The four positive optimization variables are

\[
d_{Ia},\ d_{II},\ d_{Ib},\ d_{recip};\qquad
\delta=(+d_{Ia},+d_{II},-d_{Ib},-d_{recip}).
\]

Each variable is shared across the bilateral corresponding pool contracts, with each actual pool's conductances derived by Section 2. Shared parameters do not force equal voltages under unequal sensory input or erase antiphase initial conditions. The tested controls are II-only +0.02 mV and weak-all \((+0.1,+0.02,-0.1,-0.1)\) mV. Proposed engineering search bounds are Ia 0.05–0.5 mV, II 0.005–0.15 mV, and inhibitory magnitudes 0.05–0.5 mV; the parent study owns the final consumed bounds and results.

The derived profile sets `parameter_checkpoint=None`, preventing the old reflex-gain snapshot from overwriting these new targets. Preflight verifies all 368 isolated pool contracts, exact preservation of PF/FF/biarticular gain buffers, command persistence, and **frozen timing**. With the original ±0.5 mV targets, omitting that snapshot reproduces the baseline arrays exactly. These proofs validate the fit interface, not the fit outcome.

The subsequent controlled trial lowers only the homonymous II target from +0.5 to +0.02 mV and restores stage-2 acceptance, with knee phase error 0.00748. The independent audit confirms exactly 92 changed conductances, identical physical-model binaries, and identical recorded initial neural/mechanical states before preparation. In right-side RG-E-dominant samples, psoas II current changes 0.88454 → 0.03596 nA, and actual mean activation changes 0.19074 → 0.04161. Those different-trajectory averages include the coupled motion and neural-preparation response; the audit separately verifies the immediate same-state parameter effect. Actual muscle activation remains constant during the held preparation, although neuronal state and control change. This result resolves the immediate air-stage expression failure, not the separate loaded-posture or ground-stepping problem.

The frozen function-equivalent timing budgets are E **0.500991937637077 nA** and F **0.0518116865856527 nA** at the declared 2.5 mV RG operating point; they are not current clamps. A separate balanced 0.5/0.5 nA comparator failed stage 2 more strongly and is not treated as automatically superior. Existing commissural gains remain zero; missing coupling has not been established as this failure's cause. The historical inherited heel-interneuron paths remain masked; the existing direct heel-to-contralateral RG-F event-reset path is retained, without silently stacking a second contact hypothesis.

## 7. Evidence index and scope of acceptance

All source/runtime changes are versioned separately from historical profiles. The v1 archive preserves the exact earlier source; the v2 validator selects matrix elements before subtraction, removing a dense temporary without changing coverage or neuronal mathematics. Fresh electrical comparisons are exact. Original PT binding is checked with zero body steps, including all four nonzero balance directions; this establishes binding/current consumption, not standing or walking stability. The source stage-3 verification itself fails its sway/load criteria, so the air passes do not establish supported walking.

| Evidence | Exact location |
| --- | --- |
| Integrated builder, masks, and verifiers | [integrated_hunt_design.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/integrated_hunt_design.py) |
| Primary v2 profile and source hashes | [PROFILE.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/function_equivalent/PROFILE.json>) |
| Applied local-pool ledger and numeric example | [INTEGRATED_REVIEW.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/function_equivalent/INTEGRATED_REVIEW.json>), `composition.local.design_ledger` |
| Full bound-PT constructor proof | [STATIC_BIND_PREFLIGHT.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/function_equivalent/STATIC_BIND_PREFLIGHT.json>) |
| Four-target API proof | [LOCAL_TARGET_PREFLIGHT.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/LOCAL_TARGET_PREFLIGHT.json>) |
| Actual source/new local-current reconstruction | [LOCAL_CURRENT_DIAGNOSIS.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/LOCAL_CURRENT_DIAGNOSIS.json>) |
| RG/PF versus knee timing | [AIR_PHASE_DIAGNOSIS.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/AIR_PHASE_DIAGNOSIS.json>) |
| Saved body assessments | [source S1](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/source_baseline_s1.result.json>), [source S2](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/source_baseline_s2.result.json>), [new S2](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/configurations/fsa_equivalent_s2.result.json>) |
| Exact validator comparison | [VALIDATOR_PARITY.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/VALIDATOR_PARITY.json>) |
| Frozen v1 source/evidence manifest | [MANIFEST.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v1_frozen/MANIFEST.json>) |
| Encoder adapter, discrete-phase proof, and domains | [graded_profile.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/graded_profile.py), [GRADED_GAIT_DESIGN.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/GRADED_GAIT_DESIGN.md) |
| Separate GLIF mathematical implementation | [table2_glif.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/table2_glif.py), [SPIKING_DESIGN_RECORD.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/SPIKING_DESIGN_RECORD.md) |
| Equality-aware static scaling audit, unchanged physics | [GAIT_MECHANICAL_SCALING_COMPARISON.md](<D:/Users/Ben Bolen/NeuroCampaign_20261031/fsa_design/GAIT_MECHANICAL_SCALING_COMPARISON.md>) |
| Anatomical coverage and static K6 rank | [MODEL_DIMENSION_AND_CONTROL.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/MODEL_DIMENSION_AND_CONTROL.md), [SYN6_TORQUE_SPAN.json](<D:/Users/Ben Bolen/NeuroCampaign_20261031/fsa_design/SYN6_TORQUE_SPAN.json>) |
| Reviewed dissertation family/source choices | [DISSERTATION_WIRING_DECISIONS_20261102.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/DISSERTATION_WIRING_DECISIONS_20261102.md) |
| Controlled II-only refinement and preparation clocks | [CONTROLLED_II_REFINEMENT_CHECK.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/v2/CONTROLLED_II_REFINEMENT_CHECK.json>) |
| Subsequent loaded-posture geometry and PT units | [STAGE3_REFERENCE_UNIT_REVIEW.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/STAGE3_REFERENCE_UNIT_REVIEW.md) |
| Static inherited-path direction and optional duplicate-posture repair | [POSTURE_POLICY_REPAIR.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/POSTURE_POLICY_REPAIR.md) |

For defense questions: the defended contribution is an explicit model-to-code design and calibration process with measurable local contracts, not a universal conductance gain, a mass-derived neuronal clock, or a global walking-stability proof. Shared local feedback targets recover air-stage motor expression in the controlled example; the separate loaded-posture investigation must still establish appropriate foot-relative behavior and apparatus dependence under the unchanged physical acceptance gates.
