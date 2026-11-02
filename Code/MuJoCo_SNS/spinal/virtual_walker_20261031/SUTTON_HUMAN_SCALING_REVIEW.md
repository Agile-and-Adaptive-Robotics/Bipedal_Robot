# Size, cadence, and neural-interface scaling review

Session 1, easteregg2. Read-only controller/reference review; no simulation, optimization, source-reference edit, neural retuning, manuscript edit, or Zotero write. Companion: SZCZECINSKI2020_SEVEN_STEP_CONTROLLER_AUDIT.md. Supplied source is Sutton, Szczecinski, Quinn, and Chiel, *Neural control of rhythmic limb motion is shaped by size and speed*, 2021 Research Square manuscript, DOI10.21203/rs.3.rs-153297/v1.

Source PDF: `D:\Zotero\storage\IYPGFJZ9\Gregory P. Sutton et al. - 2021 - Neural control of rhythmic limb motion is shaped b.pdf`. Figure1D was rendered and visually inspected on PDF p11. Discussion of robotic dynamic similarity is on PDF p10. This is the supplied manuscript version, not an asserted later journal version.

## What size scaling changes

The local, small-angle rotational model is `J*qddot + B*qdot + (K_elastic+K_gravity)*q = torque_active`. Sutton's isometric assumptions give:

| Quantity | Geometric length dependence |
| --- | --- |
| Limb mass | L^3 |
| Rotational inertia J | L^5 |
| Rotational damping B | L^3 |
| Elastic rotational stiffness | L^3 |
| Gravity rotational stiffness | L^4 |
| Moment arm | L |

These are similarity assumptions, not measured constitutive laws for our particular robot. At a chosen posture, gravity can destabilize rather than restore; use the signed local gravity derivative. Active muscle/BPA stiffness, multi-joint coupling, contact, and the harness also affect the effective local dynamics.

For oscillation frequency omega, useful local dimensionless comparisons are `J*omega^2/K_effective`, `B*omega/K_effective`, and `B/(2*sqrt(J*K_effective))` when the effective stiffness is positive. These compare inertial/elastic demand, viscous/elastic demand, and damping ratio. A corresponding available-torque measure is `Fmax*r/(M*g*L)`; delayed sensing and actuation also require normalized delay `delay/T_stride`.

Under these idealized scalings, an elastic-dominated natural timescale grows approximately with L, while a gravity-dominated timescale grows with sqrt(L). They cannot both be preserved by a single mass-based clock multiplier. Real animal allometry and robot materials can differ from isometry. A 75kg total mass alone does not specify leg inertia, stiffness, delay, torque reserve, or an appropriate neural period.

For Gait2392, retain the actual human muscle lengths, moment arms, maximal forces, activation and tendon parameters rather than scaling those a second time from a rat. For the CAD robot, estimate its own quantities from the CAD, routed actuator law, and actual support condition. Include auxiliary follower bodies/armature when evaluating the implemented effective inertia, and distinguish that from original OpenSim physical inertia. The stiff guided-standing condition and a free/contact-changing walker are different mechanical operating conditions.

## Fastest usable cached running stride

The read-only audit loads all43 current caches successfully. Running consists of32 trials/64 cached side cycles. All running waveforms are finite and nonconstant, have periodic=False, and duty fractions0.31–0.46.

- Shortest cached period: **0.559s**, cached `l`, subject11 nominal5m/s trial.
- Stride rate: **1.7889088Hz**; alternating-step equivalent **3.5778175steps/s**, or **214.669steps/min**.
- Paired cached `r` period:0.564s; pair mean0.5615s.
- Pooled running median:0.693250s; right minimum/median0.564/0.700s; left0.559/0.690542s.
- File: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\gait_refs\subject11_Run_50002.npz`.
- SHA256:`726902543a83ed9b66457b14e301708800e1b941be5b1d93ca24a2c2c5890f15`.

This is the shortest **usable cached estimate**, not an independently revalidated raw-recording extremum. Original subject11 IK/GRF files were not located; the generalized library loader does not enforce complete IK-window containment. Cached side assignment is based on earlier loading onset, not verified anatomical force-plate identity. The catalog points to Hamner/Delp's nmbl_running data despite older Arnold labels.100phase samples do not establish raw acquisition sampling rate. The degenerate walking cache ong_selfsel_Init200 left cycle was excluded.

## The fastest stride is not the entire required signal bandwidth

The stored phase grid includes0%and100%; omit the terminal sample and remove DC before the Fourier audit (99phase samples). For the0.559s cached side:

| Joint-angle modulation | 95%AC-energy band | 99%AC-energy band |
| --- | --- | --- |
| Hip | 3.578Hz (harmonic2) | 5.367Hz (harmonic3) |
| Knee | 3.578Hz (harmonic2) | 5.367Hz (harmonic3) |
| Ankle | 5.367Hz (harmonic3) | 7.156Hz (harmonic4) |

Endpoint mismatch exists (about2.99deg hip,2.16deg knee,0.44deg ankle); a linear endpoint-bridge sensitivity calculation leaves all95/99%harmonic counts unchanged. These are joint-angle spectral bands. They do not prove muscle-force, contact-impact, neural, or valve bandwidth. Force/velocity/command targets need their own audit.

For a single first-order transmission filter, phase lag is `atan(2*pi*f*tau)`. At1.78891Hz,20/25/30/80ms yield approximately12.67/15.70/18.63/41.96deg. A pure40ms delay instead gives25.76deg. The multiple stages' combined phase and attenuation matter. At7.156Hz, a first-order filter limited to10deg lag requires tau<=3.92ms. This is an illustrative engineering phase budget, **not a paper-prescribed value** or an instruction to set every PF/CPG state to4ms. Oscillator gates, phase integrators, intended shaping filters, and fast transmission relays perform different functions.

## Required design process across the existing models

1. Define each model's intended behavior and frequency range: standing, human walking/running, rat air stepping, or robot valve control. The0.559s cache is a human-running design input, not a mandatory cadence for every experiment.
2. Characterize actual mechanical J, damping/stiffness, torque capacity, contacts, and assistance. Use local coupled matrices where a scalar joint model is inadequate. Do not change physics merely to make copied neural gains work.
3. Express feedback in explicit physical/normalized ranges (length/l0, velocity/l0 per second, force/Fmax, support/BW, phase). Recompute nonlinear conductances for the intended function, rather than equating range normalization to arbitrary gain division.
4. Specify separate mechanical motion bandwidth, neural coding spike-rate range, spike/readout bandwidth, and actuator/valve bandwidth. Table2 Fmax is spikeHz, not maximal muscle force or strideHz.
5. For graded subnetworks apply the2017 functional method. For actual GLIF pathways derive the2020Table2 bias, membrane/threshold/synapse parameters under its assumptions. For retainedNaP or custom population observers explicitly derive and test the additional mapping instead of claiming literalTable2 completion.
6. Match the analog/digital units. Existing second-based Python Cm=.03 withGm1 is .03microfarad=30nF for a30ms time constant, not.03nF. The nF-based Simulink blocks have a1000conversion. Derive rate-frequency and time units consistently.
7. Validate each pathway's steady gain, transients, saturation, phase lag, and pulse response before full-body fitting. Use shared learned homologous parameters with source-specific known geometry. Then run Hunt's neural training, the physical curriculum, fresh repeats, perturbations, and support measurements separately.

Scharzenberger's long CPG manuscript p12 Eq23 sizes integration current using MAXIMUM active-state duration. That requirement and fastest-signal transmission constraints are complementary, not interchangeable. The ROA paper concerns autonomous-system state-space attraction; a finite +/-1%parameter basin is not that certification. Micklem's reported graded core plus24spiking valve-output neurons is another distinct architecture; its valve pulse timing must not be copied as human stride timing.

## Interpretation of current work

Existing successful graded/spiking demonstrations remain valid in their measured conditions. This audit identifies missing design provenance and signal-domain/bandwidth checks, not a retrospective cancellation of those measured results. No full controller has yet been certified as a complete2020Table2 implementation or a size-invariant75kg-human transfer. Model changes should be proposed and versioned from this ledger, not applied silently across every controller.
