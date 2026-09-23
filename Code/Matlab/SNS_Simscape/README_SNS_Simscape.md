# SNS_Simscape — Synthetic Nervous System + Knee Simscape build (Sept 2026)

Built by ZCode session 2026-09-08; visuals + URDF route + function-subnetwork
figures added 2026-09-09. Goal: wire the knee SNS circuit to a Simscape model of
`09_BA_003.SLDASM`, with an Animatlab/SNS-toolbox-style neuron library drawn in
the journal diagram language Ben specified.

**2026-09-22 LIBRARY REDESIGN (Ben's circuit-language rules, breaking change):**
`SNS_Library.slx` was rebuilt around ONE-INPUT synapses and INTERNAL neuronal
summation (details in "Block appearance conventions" and "The blocks" below):

- `NonSpikingSynapse` = **one input (Vpre) → one output**: it emits the 2-wide
  synaptic signal `[g; g*Esyn]`, and the **postsynaptic neuron** evaluates the
  driving force (`Isyn = sum(g*Esyn) - V*sum(g)` — algebraically identical to
  the old `gmax*Sat(Vpre)*(Esyn-Vpost)`, so all tuning and the numpy units
  reference stay valid). No Vpost sense wire, no Sum block in front of any
  neuron: the neuron has port 1 = `Iapp` [nA] (scalar or vector, summed
  element-wise) and ports `syn1..syn6` for synapses; **unconnected syn ports
  auto-ground to zero** (model diagnostic `UnconnectedInputMsg='none'`, the
  default). Synapse blocks are SMALL and belong right against the neuron they
  synapse onto.
- New `SynSum` junction block sums up to 8 synapse signals into one line, for
  generated big models where a neuron receives >6 synapses (spinal MNs get up
  to 17).
- Icons per Ben: non-spiking neuron = circle + **graded-potential waveform**;
  spiking neuron = circle + **spike-train waveform**; `IaMuscleSpindle` =
  **spindle-shaped capsule**; `IbGolgiTendon` = capsule with braided strands;
  muscles/BPAs = fusiform **with striations**; `MuscleActivation` =
  **pink pentagon**.
- Rebuilt + re-verified on top of it: `sns_units_test_2n` (PASS 6.07e-04 mV,
  same as before), `KneeReflexDemo`, `KneeReflexCircuit` (now a RUNNABLE 1:1
  twin of the demo — the old print-only stub could not run), `BPACPGLegDemo`,
  `BeerCupReflexDemo`, and `SNS_SpinalNetwork.slx` regenerated via the updated
  `sns_build_from_json.m` (network verify PASS 4.2e-06 mV, same as before).
- **Terminology (Ben, 2026-09-22): call things MODELS, not "plants".**
- Fixed the same day: knee animation drew the shank UP alongside the femur
  (leg read as folded 180°) — now flexion swings the shank posteriorly with
  body-fixed muscle insertions; beer animation's biceps origin was floating
  10 cm lateral to the shoulder — now shoulder → forearm; triceps added
  everywhere (model, circuit, animation).

**2026-09-10 session (easteregg2, R2025a):** restyled icons (heavy strokes,
Animatlab colors), added REAL BPA actuator blocks (10/20/40 mm, Ben's Festo
equations) + a Thelen-style biological muscle block, new demos
(BPACPGLegDemo, BeerCupReflexDemo), and a working OpenSim → Simscape route
(`osim_import/mjcf2urdf.py` + `sns_osim_import.m`).

## Folder layout (reorganized 2026-09-10)

```
SNS_Simscape\
  SNS_Library.slx            generated library (11 blocks, R2025a-native)
  sns_build_library.m        builds the library
  sns_build_actuators.m      adds BPA_10/20/40mm + BioMuscle blocks
  sns_test_actuators.m       validates actuator blocks vs Ben's equations
  snsfig.m, sns_draw_circuit.m, sns_function_subnetworks.m   journal figures
  sns_export_diagram.m, sns_osim_import.m, sns_cpg_gait2392.m, sns_urdf_smoke.m
  import_simscape_when_ready.m, check_env.m
  demos\                     demo builders/runners + generated demo .slx
    KneeReflexDemo.slx  BPACPGLegDemo.slx  BeerCupReflexDemo.slx
    sns_build_demo.m / sns_run_demo.m
    sns_build_cpg_demo.m / sns_run_cpg_demo.m
    sns_build_beer_demo.m / sns_run_beer_demo.m
    sns_animate_demo.m (2D mechanism animations)
  results\                   run outputs (.mat)
    pictures\                signal plots + diagnostics
    animations\              GIF mechanism animations
  figures\                   journal figures (PNG/PDF/SVG for LaTeX)
  logs\                      run logs
  dev\                       tuning/debug scratch (tune_cpg*.m, dump_wiring*.m, ...)
  osim_import\               OpenSim -> Simscape (mjcf2urdf.py, gait2392 URDF)
  sw_import\                 SolidWorks COM probes
  urdf_smoke\, simscape_sources\   smoke test + legacy Simscape-lib route
```

The `_R2025a` export copies were deleted 2026-09-10: everything is rebuilt
R2025a-native (verified on easteregg2 R2025a), and R2025b (laptop) opens
R2025a files directly.

**2026-09-20 Rybak restyle:** `SNS_Library.slx` was REBUILT on the laptop
(Rybak pass-through synapse icons, below) and is now R2025b-NATIVE; the
demos' library links still resolve (unchanged block paths), and
`SNS_Library_R2025a.slx` is exported alongside for inspection. R2025a
machines that need to open the demos should simply re-run
`sns_build_library.m` + `sns_build_actuators.m` locally (they regenerate an
R2025a-native library; validation = `sns_test_actuators.m`, PASS 4.9e-10 N).

## Run order (from `SNS_Simscape\`, or run any script by full path)

```matlab
sns_build_library      % 7 neural/mechanical blocks
sns_build_actuators    % + BPA_10mm/20mm/40mm + BioMuscle  -> 11 blocks
sns_build_demo         % KneeReflexDemo.slx
sns_run_demo           % -> results\pictures\sns_demo_results.png
sns_build_cpg_demo     % BPACPGLegDemo.slx
sns_run_cpg_demo       % -> results\pictures\sns_cpg_results.png
sns_build_beer_demo    % BeerCupReflexDemo.slx
sns_run_beer_demo      % -> results\pictures\sns_beer_results.png
sns_draw_circuit       % -> figures\KneeReflex_circuit.{png,pdf,svg}
sns_function_subnetworks
sns_animate_demo       % -> results\animations\*.gif  (in demos\)
```

## Deng 2019 RG/PF layers as separate Simulink files (2026-09-13, VERIFIED)

Per Ben's request: the Deng 2019 (Biomimetics 4(1):21) two-layer CPG — the
architecture of the Animatlab biped port — as separate editable files, with
parameters from Nourse 2023 Tables A4–A7 (`Code\MuJoCo_SNS\spinal\
_nourse2023.txt` is the extracted paper text):

- `SNS_Deng_Library.slx` — **HCNeuron**: persistent-Na half-center neuron
  (Cm 5 nF, Gm 1 µS, Vrest −60 mV, GNa 1.5 µS, ENa 50 mV, m: S 0.2 / E −40 /
  K 1, h: S −0.6 / E −60 / K 0.5, **tau_h FIXED 350 ms**). Rebuild via
  `demos\sns_build_deng_cpg.m` after any change.
- `demos\SNS_Deng_RG.slx` — RG layer: HC_ext/HC_flx + IN-laminated mutual
  inhibition (HC −2.749→ IN −2.749→ HC, Esyn −40/−70 mV, window −60..−25);
  `I_stim` inport; V_ext/V_flx outports.
- `demos\SNS_Deng_PF.slx` — PF layer: hip pair + knee/ankle pair, same
  construction; RG→PF weak exc 0.1 µS; 2 in / 4 out.
- `demos\SNS_Deng_CPGDemo.slx` + `sns_run_deng_demo.m` — ONE 10 nA / 20 ms
  pulse at t = 0.1 s → **continuous oscillation, period 1.938 s (matches the
  numpy ODE reference 1.94 s), RG ext/flx correlation −0.86**, hip MN
  activation alternating via the Deng Fig 6B sigmoid
  (act = 1/(1+e^{0.1532(−70−V)}) − 0.01, MN Vrest −100 mV, PF→MN hip
  2.565/3.632 µS).

**tau_h is the load-bearing choice**: sns_toolbox's tau_h(V) collapses at
depolarized V and QUENCHES this circuit (verified `Code\MuJoCo_SNS\spinal\
deng_cpg_ode.py`: toolbox-tau → flatline, fixed-tau → self-sustained 1.94 s
bursting). Animatlab's port uses tau_h.max as a fixed constant — when
debugging why the Animatlab RG latches instead of oscillating, check its
Na-channel h time-constant handling first. Same lesson applies to the
h-inf/+1 terms: the reciprocal-divider blocks MUST include the +1.

**CORRECTION (2026-09-16, Ben caught it): sns_toolbox 1.5.2 DOES ship
`NonSpikingNeuronWithPersistentSodiumChannel`** (Tutorial 8; constructor:
membrane_capacitance, membrane_conductance, g_ion, e_ion, k_m/slope_m/e_m,
k_h/slope_h/e_h, tau_max_h, name, color). Earlier claims that "the toolbox
can't express persistent-Na dynamics, so the ADAP loop is the only
burst-termination substitute" were WRONG — Tutorial 8 was even executed in
`spinal\sns_tutorials` (2026-09-14). Implications for the ports:

- The numpy/MuJoCo spinal RG half-centers CAN be built as literal Deng-style
  persistent-Na HC neurons with intrinsic burst termination (replacing or
  re-testing the ADAP-loop workaround) — this is the cross-platform
  consistency route: Simulink `SNS_Deng_Library.slx` (fixed tau_h 350 ms)
  and AnimatLab LinearHill (tau_h.max as fixed constant) already work this
  way, so one neuron family could span all three platforms.
- CAVEAT before trusting it: the quenching result above was measured against
  the toolbox tau_h(V) FORMULA as hand-coded in `deng_cpg_ode.py` — not the
  class. Test how the real class's tau_h(V) behaves with tau_max_h at
  depolarized V in the actual Deng circuit first; if it collapses the same
  way, check whether the class allows the fixed-tau_h semantics the working
  ports rely on.
- When porting values across platforms, the standing trap still applies:
  SNS_Library (Simulink/AnimatLab side) vs sns_toolbox (python side) use
  OPPOSITE synapse saturation conventions (ThrPre/Elo).

Solver: fixed-step ode1 @ 0.1 ms (Table A7 dt). The numpy ODE reference is
`deng_cpg_ode.py` (results saved to `results\deng_cpg_ref.mat`); the pointwise
comparison to Simulink is phase-sensitive (onset timing differs), so compare
periods, not samples.

## Tuned spinal network → editable Simulink model (2026-09-12, VERIFIED)

Pipeline that turns the tuned `runner --fitted --best` gait2392 spinal network
(the v4b winner) into an editable Simulink model built from SNS_Library blocks:

1. **Export** (myo env, from `Code\MuJoCo_SNS\spinal\`):
   `python export_network_json.py` → `spinal_net_export.json`
   (410 neurons, 1392 synapses, 376 input ports, 92 MN→actuator outputs,
   composited exactly like `runner --fitted --best`). Companion refs:
   `export_units_ref_2n.py` (2-neuron units test reference) and
   `export_verify_ref.py` (full-network reference, coarse 2 ms Euler + fine
   0.1 ms).
2. **Units test**: `sns_units_test_2n.m` builds the same 2-neuron + synapse
   circuit from `SNS_Library` blocks and matches the SNS-Toolbox numpy
   backend. **PASS: max dev 6.07e-04 mV.** Mapping:
   toolbox `C = tau uF` → block `Cm = 1000*tau nF`; `Gm = 1 uS`; `Vrest = 0`;
   synapse `e_lo/e_hi = 0/5 mV` → `ThrPre = 0 / SlopePre = 5`; `g` uS,
   `Esyn` mV, currents nA.
3. **Generate**: `sns_build_from_json.m` → `results\SNS_SpinalNetwork.slx`
   (~2400 blocks: one NonSpikingNeuron per cell, one one-input
   NonSpikingSynapse per connection landing on the postsynaptic neuron's
   syn1..syn6 port, per-neuron input Muxes onto the Iapp port for the 376
   external input currents, chained SynSum junctions for the neurons with
   more than 6 incoming synapses (146 junctions; MN in-degree reaches 17),
   and a 92-wide Mux of MN `S(V)` outputs → outport `S`, ordered by MuJoCo
   actuator id = ctrl order). All values live in block masks — double-click
   to edit. Regenerate any time the tuning changes.
4. **Verify**: `sns_verify_from_json.m` — **PASS: 4.2e-06 mV** max deviation,
   ALL 410 neurons at t = 0.3 s vs the numpy 2 ms Euler reference
   (DRIVE = 2.5 nA; V_DRIVE Euler-exact to 9 decimals).

**Solver guidance (measured, important):** the network is CHAOTIC — any two
integrators agree only to ~0.4 s (1e-9 @ 0.1 s → 1e-6 @ 0.3 s → O(1) @ 0.5 s;
same effect build_network.py documents as the "reggate_v5_0 lesson"). To
reproduce the tuned production trajectories, run the model **fixed-step
ode1 (Euler) at exactly 0.002 s** — that is bit-compatible with the runner's
numpy stepping. Variable-step ode45 integrates the "true" dynamics but its
trajectory departs from the tuned one after ~0.4 s (and the fine-dt numpy
reference shows the RG period itself shifts ~35% between 0.1 ms and 2 ms
stepping — the v4b tuning is a property of the 2 ms semantics).

Gotchas: `ExternalInput` on `sim` proved unreliable in `-batch` — inject via
a Constant into the demux instead; port-level `DataLogging` (not To
Workspace, which rejects the bus) for logging; Simulink subsystem ports are
referenced numerically (`blk/1`), not by inner port-block names.

## What runs today (MATLAB R2025a and R2025b)

| File | What it is |
|---|---|
| `sns_build_library.m` | builds **`SNS_Library.slx`** — 7 masked blocks with diagram-language icons |
| `sns_build_actuators.m` | adds **BPA_10mm / BPA_20mm / BPA_40mm** (Ben's real Festo equations) + **BioMuscle** (Thelen-style) → 11 blocks |
| `sns_test_actuators.m` | validates the actuator blocks against `festo4.m`/`maxBPAforce.m` + Thelen references |
| `demos\sns_build_demo.m` | builds **`KneeReflexDemo.slx`** — reflex circuit + 1-DOF knee model |
| `demos\sns_run_demo.m` | simulates 5 s, saves results into `results\` |
| `demos\sns_build_cpg_demo.m` + `sns_run_cpg_demo.m` | **`BPACPGLegDemo.slx`** — half-center CPG (mutual inhibition + adaptive inhibition) driving antagonist BPAs on the 1-DOF knee |
| `demos\sns_build_beer_demo.m` + `sns_run_beer_demo.m` | **`BeerCupReflexDemo.slx`** — elbow holds a cup level via Ia/Ib reflex while beer pours in; runs reflex ON vs OFF |
| `demos\sns_animate_demo.m` | replays logged joint motion as 2D mechanism GIFs (`results\animations\`) |
| `sns_draw_circuit.m` | **journal figure**: redraws the demo circuit reader-facing (`figures/KneeReflex_circuit.{png,pdf,svg}`) |
| `sns_function_subnetworks.m` | **journal figure**: 6 arithmetic subnetwork panels (`figures/SNS_function_subnetworks.{png,pdf,svg}`) |
| `snsfig.m` | drawing primitives shared by both figure scripts (neuron/triangle/dot/muscle/box/arrow) |
| `sns_export_diagram.m` | prints Simulink models to `figures/<model>_simulink.{eps,pdf,png}` **as dissertation figures**: per-block 12 pt fonts, annotations normalized 10–12 pt, synapse names hidden (Rybak: a synapse is a connection, not a labelled component), white canvas; styling is applied in-memory only — the .slx files are never modified. EPS = color vector (painters) for `\includegraphics{*.eps}`; Overleaf handles it, but the twin .pdf is the zero-friction option for local pdflatex. |
| `sns_urdf_smoke.m` | minimal URDF → `smimport` smoke test (fails on easteregg2 license, see below) |
| `osim_import/mjcf2urdf.py` + `sns_osim_import.m` | **OpenSim → Simscape**: MyoConverter MJCF (cvt3) → URDF → `smimport`; saves `osim_import/Gait2392_simbody_simscape.slx` |
| `sw_import/sw_probe.py`, `sw_probe2.py`, `sw_mate_probe.m` | SolidWorks COM probes of `09_BA_003.SLDASM` (components + transforms work; mate entities blocked by pywin32 byref bug) |
| `import_simscape_when_ready.m` | one-command `smimport` of the real CAD (URDF **or** Multibody-Link XML) |
| `SNS_Library.slx`, `demos\*.slx` | generated models (open in Simulink, blocks are double-click editable) |
| `figures\`, `results\`, `logs\` | figures, run outputs + animations, run logs |

## Block appearance conventions (Ben, 2026-09-09 + 2026-09-22)

Diagram language follows **Szczecinski et al. 2017 Fig. 2** (the functional
subnetwork paper), Rybak/Shevtsova CPG diagrams, and Animatlab:

| Element | Icon |
|---|---|
| neuron (non-spiking RC) | open circle with a **graded-potential waveform** (smooth depolarizing hump) |
| spiking LIF neuron | open circle with a **spike-train waveform** |
| Ia afferent | **spindle-shaped capsule** (fusiform, tapered ends) labeled "Ia" |
| Ib afferent | capsule with braided collagen strands, labeled "Ib" |
| muscle / BPA / BioMuscle | light-salmon fusiform **with striations** across the belly |
| MuscleActivation | **pink pentagon** (new mask param `A0` = initial activation) |
| **EXCITATORY connection** | synapse block: pass-through axon bar terminating in a **WHITE TRIANGLE, black edges** at the output edge (tip points back toward the presynaptic side) |
| **INHIBITORY connection** | same axon ending in a **SOLID BLACK CIRCLE** at the output edge |

- The `NonSpikingSynapse` icon picks its marker **automatically from the sign of
  `Esyn`** (`Esyn >= 0` → triangle, `< 0` → black dot), so the icon always tells
  the truth about the connection. An un-evaluable `Esyn` expression draws "E?".
- Synapse blocks are SMALL (draw ~40×32) and are placed **immediately left of
  the neuron they synapse onto**, wired into its next free `syn1..syn6` port.
- Shape is the primary code — figures stay readable in grayscale and for
  colorblind readers. Tints are the redundant cue, from the **Okabe-Ito
  CVD-safe palette** (orange = excitatory, blue = inhibitory, green = muscle).
- Keep CIRCLE-icon blocks **square** (width == height): mask icons autoscale to
  the block rectangle, so non-square blocks turn circles into ellipses.
- Simulink gotcha (cost us a debug cycle): mask drawing commands accept
  **numbers only** — no LineSpec strings (`'k-'`), no name-value pairs
  (`'LineWidth'`,…), and `color()` takes a color NAME (`color('black')`), not
  RGB. Fill with `patch(...,[r g b])`, set edge color with `color('black')`
  before `plot`, center labels with `disp()`.

## The blocks (SNS_Library.slx)

Neuron = **non-spiking leaky integrate-and-fire, literally an RC membrane**:
`Cm·dV/dt = Gm·(Vrest − V) + Iapp + Σ g_k·(Esyn_k − V)` — the neuron does the
synaptic summation INTERNALLY. τ_m = Cm/Gm (nF/µS = ms).
Port 1 `Iapp` [nA]: injected current — descending drives and afferent currents
(vector inputs are summed element-wise, so several currents can share it).
Ports `syn1..syn6`: synapse inputs `[g; g·Esyn]` from `NonSpikingSynapse`
blocks; unconnected ports count as zero, so a neuron with fewer synapses needs
no dummy wiring. Outputs membrane `V [mV]` and normalized drive
`S = clip((V−Thr)/Slope, 0, 1)`.

Synapse = **one input, one output**: in = presynaptic `V` [mV]; out =
`[g; g·Esyn]` with `g = gmax·Sat(Vpre)`, `Sat = clip((Vpre−ThrPre)/SlopePre, 0, 1)`.
**Excitatory = Esyn 0 mV, inhibitory = Esyn −72 mV** — an E vs I connection
differs only by the Esyn mask value (and by the icon marker).
`SynSum` = junction for generated models: sums up to eight synapse outputs into
one `[Σg; Σg·Esyn]` line (chainable), for neurons receiving more than 6 synapses.

- `NonSpikingNeuron` — RC membrane with internal synaptic summation (defaults Vrest −52 mV, Gm 0.1 µS, Cm 5 nF ⇒ τ=50 ms; SNS-toolbox conventions)
- `NonSpikingSynapse` — one-in/one-out E/I chemical synapse with auto E/I icon (ThrPre −55 mV, SlopePre 1/mV ⇒ off at rest, graded above)
- `SynSum` — synapse-summing junction (up to 8 inputs, grounded-unused)
- `SpikingLIFNeuron` — spiking LIF with threshold reset (Animatlab-style spiking neuron; scalar/vector Iapp input)
- `IaMuscleSpindle` — stretch + velocity afferent → current (peak 10 nA)
- `IbGolgiTendon` — force afferent → current
- `MuscleActivation` — first-order activation dynamics (τ_act 50 ms; A0 initial activation)
- `BPAForce` — `F = Fmax·A·max(0, epsMax − strain)` placeholder; swap in the Xi-corrected
  MonoPam prediction from Mesh_Optimization when ready

## Function subnetworks (the "how do connections compute" documentation)

Szczecinski, Hunt, Quintin, et al. 2017, *A Functional Subnetwork Approach to
Designing Synthetic Nervous Systems That Control Legged Robot Locomotion*,
Front. Neurorobotics 11:37 (DOI 10.3389/fnbot.2017.00037) — extension to the
generalized LIF neuron: DOI 10.3389/fnbot.2020.577804. Their Fig. 2 caption is
the origin of our marker convention ("triangular synaptic terminations stand
for excitatory inputs, and filled round terminations for inhibitory inputs").

`figures/SNS_function_subnetworks.{png,pdf,svg}` shows all six operations in
our marker language:

| Operation | Wiring | Key params (paper, R = 20 mV) |
|---|---|---|
| addition | two excitatory *signal-transmission* synapses onto one neuron | ΔEs = 194 mV, gs = 115 nS, k_syn = 1 |
| subtraction | excitation + inhibition; weights cancel | gs,2 = −gs,1·ΔEs,1/ΔEs,2 (ΔEs,2 < 0, e.g. −40 mV) |
| division | transmission + *signal modulation* (shunting) | ΔEs,2 = 0 (GABA-like); c_syn = 1/R, gs,2 = 19 µS |
| multiplication | two modulatory synapses **in series** (disinhibitory cascade) | ΔEs = −1 mV, gs = 20 µS; mid-neuron Iapp = R |
| differentiation | same input into two neurons, Cm,1 < Cm,2, subtract (Reichardt) | k_d = Cm,2 − Cm,1, τ_d = Cm,2 |
| integration | mutual **self-disinhibition** line attractor; u into U1, tonic R into U2 | gs,2 = gs,1, ΔEs = −R; ḊU1 = k_i·u |

Mapping to OUR block parameters: the paper's piecewise-linear synapse
`gs·(Vpre−Elo)/(Ehi−Elo)·(ΔEs − Vpost)` equals our
`gmax·Sat(Vpre)·(Esyn − Vpost)` with `gmax ≡ gs`, `ThrPre ≡ Elo`,
`SlopePre ≡ 1/(Ehi−Elo)`, `Esyn ≡ Vrest,post + ΔEs`. With our demo defaults
(Vrest −52 mV, R_pre = 1/0.5 = 2 mV) a k_syn = 1 excitatory synapse is
`gmax = gs·ΔEs/R_pre`; tune `gmax` for gain and `Esyn` for sign, exactly as in
the paper's Table 1 constraints. (Runnable primitive models in Simulink are a
natural follow-up; the figures + this table are the current documentation.)

## KneeReflexDemo.slx

Antagonist BPA pair on a reduced-order 1-DOF knee model (mechanics grouped
into one masked `KneeModel` subsystem; all parameters unchanged from the
2026-09-10 build), with sensory neurons between afferents and synapses:

- Ia(ext) → **Exc** → MN_ext (stretch reflex)
- Ib(ext) → **Inh** → MN_ext (autogenic inhibition)
- Ia(flex) → **Inh** → MN_ext, Ia(ext) → **Inh** → MN_flex (reciprocal inhibition)
- Ia(flex) → **Exc** → MN_flex, Ib(flex) → **Inh** → MN_flex

2026-09-22 architecture rebuild, verified against an independent MATLAB ODE of
the same circuit math (rise 15.8° @ 0.05 s identical; settle band 38–45° with
means within 3°; same antagonist activation alternation): the knee rises from
15° and settles around **43.5° mean (range ~41.6–44.9°)** against the 0.5 N·m
load, with a mild irregular antagonist alternation (~0.5–2 Hz dominant peak
depending on the window; A_ext/A_flex both sweep ~0.26–0.78 around 0.53).
The final SAMPLE of A_ext/A_flex depends on the alternation phase — the
old "~9 Hz limit cycle" note overstated the frequency; compare windowed
means, not final samples. Its physiological significance has not been
evaluated. Tune: `gmax` (loop gain), SN `Cm` (loop delay), `b_knee`, `tauAct`.
Model params (inertia, moment arm, Fmax) are placeholder rig estimates —
replace from CAD mass properties and the Xi-corrected BPA predictions.

## Real BPA actuator blocks (2026-09-10, Ben's equations)

`SNS_Library.slx` now has 11 blocks: the 7 neural/mechanical blocks plus
**BPA_10mm / BPA_20mm / BPA_40mm** and **BioMuscle**.

- **BPA_xx mm** — inputs P [kPa], L [m]; output F [N]. Ben's Festo sfit
  surfaces (from `Functions\festo_lookup_coeffs.json`, same surface
  `festo4.m` evaluates): `Fn = a0*(exp(-a1*rel) - 1) + (P/620)*exp(-a3*rel^2)`,
  rel = contraction/KMAX, F = Fn*Fmax, zero force for rel>=1 and F<0.
  Fmax = `maxBPAforce(Rest,620)` (atan expression) is baked into the internal
  gain — double-click into the block to override. 40 mm default 6000 N (the
  class value; `maxBPAforce.m` says 6398.4 — unverified which Ben wants).
- **BioMuscle** — OpenSim Thelen2003-style Hill muscle, rigid tendon:
  active FL `exp(-(lambda-1)^2/0.45)`, Hill force-velocity
  (A_hill=0.25, v_max in l_opt/s), exponential passive FL (kpe=4, e_pas=0.6).
  Inputs Act, Lmt, Vmt; elementary blocks only (no codegen).
- **Validation** (`sns_test_actuators.m`): BPA blocks match `festo4()*maxBPAforce()`
  to ~1e-10 N over length and pressure sweeps; BioMuscle matches the Thelen
  references to ~1e-14 N (isometric / concentric / passive).
- Gotcha learned the hard way: mask-param DEFAULT strings that reference
  sibling mask params (Fmax referencing Rest) do NOT resolve at sim time —
  block params *inside* the masked subsystem do. And a library must be LOCKED
  for its blocks to be instance-able.

## Demos verified on the LAPTOP R2025b (2026-09-20, and RE-VERIFIED 2026-09-22 on the redesigned library)

2026-09-20 originals (`dev/run_all_demos_laptop_20260920.m`, log
`logs/run_all_demos_laptop_20260920.log`); 2026-09-22 column = rebuilt on the
one-input-synapse architecture:

| Demo | 2026-09-20 | 2026-09-22 rebuild |
|---|---|---|
| KneeReflexDemo | theta final 42.72 deg | settle mean 43.5° (range 41.6–44.9°), same rise/settle/alternation vs independent ODE of the same math |
| BPACPGLegDemo | theta range 9.7..48.0 deg | theta range 9.7..48.0 deg (identical to committed results) |
| BeerCupReflexDemo | max sag ON 5.13 / OFF 9.66 deg — both were STARTUP transients at t≈0.12–0.2 s | starts in equilibrium: ON settles 2.0° vs OFF 7.9°; triceps 0.40→0.32 |
| SNS_Deng_CPGDemo | 10 RG_ext bursts / 20 s, period 1.938 s | unchanged (separate library) |
| sns_units_test_2n | PASS 6.07e-04 mV | PASS 6.07e-04 mV |
| SNS_SpinalNetwork verify | PASS 4.2e-06 mV | PASS 4.2e-06 mV (regenerated) |

**The old beer "max sag" numbers WERE the model settling, not the pour**:
measured from the committed results mat, the 5.13°/9.66° peaks occur at
t = 0.12–0.20 s with activation starting at 0 while the descending drive was
a step sized for the FULL cup — force ramped through the balance point and
overshot. The 2026-09-22 build starts in equilibrium (activation initialized
at the hold value via `MuscleActivation.A0`, drive sized for the empty cup),
so the plotted transient is the pour itself.

Fix shipped (2026-09-20): `sns_run_deng_demo.m` used `corr()` (Statistics
Toolbox — not on the laptop license); replaced with base-MATLAB Pearson. The
Deng sim itself was always fine — only the check crashed.

## Dissertation figure pipeline (2026-09-20, Ben's 10–12 pt / EPS / Rybak ask)

- **`sns_export_diagram.m`** (upgraded): per model writes
  `figures/<model>_simulink.{eps,pdf,png}`. Styling applied in-memory only
  (never saved): every block FontSize → 12 pt (R2025b has NO model-level
  font param — the model `FontSize` set_param errors; blocks carry their
  own), annotations normalized to 10–12 pt, synapse names hidden (Rybak: a
  synapse is a connection, not a labelled component), white canvas.
- **EPS**: `print -depsc` is REFUSED for Simulink systems on R2025b
  ("'epsc' format is not supported with Simulink or Stateflow printing") and
  R2025b ships no ghostscript. Working chain = `print -dpdf -bestfit` →
  **MiKTeX pdfcrop** (user install, `%LOCALAPPDATA%\Programs\MiKTeX\...`) →
  **MiKTeX mgs.exe** (`-sDEVICE=eps2write`) → tight vector EPS. pdfcrop
  cannot write to `C:\` root — the temp crop file lives in the output folder.
  In Overleaf, `\includegraphics{file.eps}` works as-is; locally, prefer the
  twin `.pdf` (identical vector content).
- **`demos/sns_build_circuit_view.m`** builds **`KneeReflexCircuit.slx`** —
  since the 2026-09-22 redesign the demo's top level is already print-clean,
  so the circuit view is simply a runnable 1:1 TWIN of KneeReflexDemo (same
  blocks, mask values, wiring — simulates standalone). The full
  KneeReflexDemo export stays available as the "faithful model" appendix
  figure.
- **`snsfig.m` / `sns_draw_circuit.m`** (vector circuit figure): fonts raised
  to 9.5–11 pt and the canvas set to 16.5 cm = dissertation text width, so
  LaTeX includes it at 100 % scale and the fonts print at face value.

## Demos

| Demo | File | What it shows |
|---|---|---|
| `KneeReflexDemo.slx` | `sns_build_demo.m` / `sns_run_demo.m` | Ia/Ib reflex circuit + antagonist BPAs on a 1-DOF knee model (theta settles ~43° mean with mild antagonist alternation) |
| `KneeReflexCircuit.slx` | `sns_build_circuit_view.m` | RUNNABLE 1:1 circuit view of KneeReflexDemo (identical blocks/values/wiring; simulates standalone) |
| `BPACPGLegDemo.slx` | `sns_build_cpg_demo.m` / `sns_run_cpg_demo.m` | half-center CPG (mutual inhibition + adaptive inhibition via slow Adp neurons) driving antagonist BPAs on the knee; theta sweeps 10-48° (verified identical to the committed 2026-09-20 results) |
| `BeerCupReflexDemo.slx` | `sns_build_beer_demo.m` / `sns_run_beer_demo.m` | elbow holds a cup level while beer pours (0→0.5 kg); BICEPS + TRICEPS BPA_20mm pair; starts in EQUILIBRIUM (no startup transient); Ia stretch reflex + Ib autogenic inhibition + reciprocal Ia inhibition. Reflex ON settles at 2.0° sag vs OFF 7.9°; biceps activation rises 0.39→0.42 while TRICEPS activation FALLS 0.40→0.32 (antagonist inhibition, plotted) |

CPG tuning was done in ODE prototypes first (`tune_cpg.m`, `tune_cpg2.m`):
mutual inhibition 0.8, adaptive current = g_adp(15 nA) x S_adp (LINEAR gain on
the Adp neuron's S output — a synapse there LATCHES the winner because its
driving force collapses at high RG voltage), asymmetric drives 4.0/3.6 nA.
Symmetric drives latch (same symptom as the AnimatLab RG latch).

## OpenSim -> Simscape (COMPLETE + VERIFIED on the laptop, 2026-09-20)

Route: OpenSim gait2392_simbody --MyoConverter--> MJCF (`mjc\gait2392_simbody\
gait2392_simbody_cvt3.xml`) --`osim_import/mjcf2urdf.py`--> URDF -->
`smimport`. The converter handles MuJoCo's default-type hinges, joint anchor
offsets (via `_pre`/`_anchor`/`-jpos` fixed-link chains), the massless muscle
pathpoint slide bodies (kept as links so BPAs can attach at the same points),
and writes proper URDF `<mass>` child elements. Output:
`osim_import/gait2392_simbody.urdf` (149 links / 149 joints).

**STATUS (laptop R2025b, 2026-09-20): `osim_import/Gait2392_simbody_simscape.slx`
IMPORTS, COMPILES, AND SIMULATES.** smimport built 1446 blocks; joint inventory:
85 Prismatic (the massless pathpoint slides) + 20 Revolute Joint (anatomical)
+ 44 Weld (the `_pre`/`_anchor` offset chains) = 149 total. The easteregg2
license blocker no longer applies anywhere that carries SimMechanics.

Two traps fixed on the way (both cost a debugging cycle):

1. **Mesh resolution**: the URDF references `Geometry/*.stl` for the 19
   anatomical Visual blocks. smimport froze those paths at import time, and
   `update` failed with "Geometry/File Name is a file that does not exist"
   for all 19 (a junction created AFTER import does NOT fix an already-saved
   model — the params must be patched). Fix shipped: `osim_import\Geometry`
   is a junction into the cvt3 Geometry folder (gitignored), and the 19 File
   Solid blocks were patched to ABSOLUTE paths via
   `set_param(blk,'ExtGeomFileName',absPath)` — **the mesh param is
   `ExtGeomFileName`** (probed; 'FileName' does not exist).
2. **Crash hazard**: a model whose `set_param('SimulationCommand','update')`
   FAILED leaves Simscape's GUI tree in a state that crashes MATLAB
   (physmod_sm_gui_app_tree.dll access violation) at ANY later touch,
   including process teardown. Always make the compile succeed, or expect the
   crash-at-exit (the log up to that point survives).

Muscles/tendons themselves are NOT imported (URDF has no muscle element) —
they become SNS_Library blocks driving the joints
(`sns_cpg_gait2392.m` is the skeleton for that wiring).

## SolidWorks knee rig -> Simscape (XML route inventoried, 2026-09-20)

`dev/imports2_20260920.m` ran the REAL Multibody-Link XML import
(`Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.xml`) →
**`mdl_knee_rig_xml_imported.slx`**: 96 blocks; joint inventory =
**2 Cylindrical + 3 6-DOF and ZERO Revolute joints — the knee DOF is
MISSING.** The dropped SW Hinge mates (exporter "not supported") became
rigid welds ("unknown constraint" warnings at import), and the KB↔TI
Concentric+Coincident pair did NOT translate into a revolute. Ben's planned
fix (replace the 4 Hinge mates with Concentric+Coincident pairs in SW and
re-export) is still REQUIRED before this route yields a usable model — the
knee angle mate carries the joint state target.
NOTE: the older `mdl_knee_rig_import_tmp_imported.slx` (12 blocks, 0 joints)
is the 1-link sw2urdf STUB from 2026-09-16, not this import.

## SolidWorks -> Simscape status

`sw_import/sw_probe.py` (root-Anaconda python — it has pywin32; the myo env
does NOT) opens `09_BA_003.SLDASM` READ-ONLY via COM and enumerates all 15
components + transforms (ground = `04_02_KB_R_003`, confirmed). SolidWorks
zero-arg methods come back as pywin32 PROPERTIES (use the zcall helper).
Blockers: transforms/mates need typed dispatch (EnsureDispatch refuses
makepy for SW — needs a one-time manual `makepy`), and the sw2urdf GUI export
dialog remains the 2-minute manual path (Ben: File > Export as URDF).

Status of the sw2urdf add-in on this laptop (SW 2025 SP4.1):
- **Not installed.** Installer downloaded: `C:\Users\Ben\Downloads\sw2urdfSetup_1.6.1.exe`.
- ⚠ The official exporter (repo `ros/solidworks_urdf_exporter`, v1.6.1, Nov 2021)
  targets **SolidWorks 2021** (min 2018 SP5). Community reports on SW2024+ are
  mixed (e.g. GitHub issue #147, export dialog vanishing). If it misbehaves on
  SW2025, fallbacks in order:
  1. **Simscape Multibody Link** — it turns out to be **installed and
     registered** on this machine (`HKLM\...\AddIns\{2666BDBF-...}`), just
     disabled at startup. Enable via SolidWorks Tools → Add-Ins, export XML,
     then `import_simscape_when_ready.m` (it accepts the XML too).
  2. Write our own URDF export via the SolidWorks COM API (the knee assembly
     needs only 4 hinge mates + ground, so a purpose-built exporter is small).

Ben's install/export steps (GUI):
1. Run `Downloads\sw2urdfSetup_1.6.1.exe` (close SolidWorks first), restart SolidWorks.
2. Tools → Add-Ins → check **SW2URDF** (if listed as "Start-up", leave enabled).
3. Open `Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.SLDASM`.
4. Tools → **Export as URDF** (v1.6 moved it from File to Tools menu):
   choose **Automatic** export type, pick the reference frame (ground bracket
   origin), accept the auto-detected joints (the 4 Hinge mates), name the main
   revolute joint `knee`, then **Export** → save as
   `...\Knee assembly\09_BA_003.URDF`.
5. Set `inputPath` in `import_simscape_when_ready.m` to that file and run it.

Observed quirk: on the smoke URDF, `smimport`'s returned model name was a
numeric string ("1.540001e+02") rather than the robot name — harmless for
`load_system` via the returned handle, but don't assume the model is named
after the URDF `<robot>`; use the return value.

## CAD work done on 09_BA_003 (via COM API, session 2026-09-08)

- Inspected: 15 components, ground = fixed `04_02_KB_R_003` (knee bracket), 4 hinge
  mates, `Knee Angle` planar-angle mate, PathMates 1/2/3/5 on `Sketch_20mm_route`,
  existing `RotaryMotor1` motion motor, no PDM vault conflicts.
- Exported **`09_BA_003_export.step`** (14.4 MB, AP214) next to the assembly.
- Renders: `09_BA_003_iso.png`, `09_BA_003_front.png`.
- Built **`Tendon_Extensor.SLDPRT`** (2×2 mm section, 541.9 mm) and
  **`Tendon_Flexor.SLDPRT`** (509.4 mm) in the same folder, from the route-sketch
  landmarks: extensor BPA2 endcap (−13.88, −56.10, 0) → tibia insertion (−24.08, −597.94, 0);
  flexor BPA1 endcap (−79.5, −15.6, 103.8) → knee flexor bracket (−43.1, −518.0, 28.2).
  **Not yet inserted into the assembly** (COM interface version mismatch on
  `AddComponent4`'s transform argument). 2-minute GUI insert:
  File → Open the assembly → Insert Components → pick each tendon part → place near
  its endpoints → mate: Concentric (tendon end face ↔ BPA endcap hole) + Concentric
  (other end ↔ bracket hole) + one Coincident/Parallel to lock rotation. `09_BA_003`
  itself was NOT modified or saved.
- Motion study: the assembly already has `RotaryMotor1` + a Results folder — in the
  GUI: Motion Study tab → Calculate → save AVI (API route needs the Motion Study
  Manager add-in, skipped today).

## Gotchas logged for future sessions

- pywin32 + SW 33.4.1: dynamic dispatch works for OpenDoc6 (VARIANT byref ints),
  SaveAs via `Extension`, `GetComponents`, `FeatureManager.GetFeatures`,
  sketch/feature creation. Members that return objects as *properties*
  (FirstFeature, FeatureManager, Extension, GetMathUtility) must be accessed without
  calling. `CreateCornerRectangle` fails with negative corner coords and needs
  `SketchManager.AddToDB = True` for small (few-mm) geometry. Mask params on the
  gen_py-cached interfaces are stale — create MathTransforms via the typed
  `IMathUtility` wrapper (see `sw_lib2.py` pattern in Temp).
- Simulink library building: delete the template In1/Out1 **before** adding ports and
  always set explicit `'Port'` — otherwise port order scrambles; Sum blocks need
  explicit port-count strings (`'++'`, `'+-'`); R2025b requires `Simulink.Mask.create`
  (old `MaskPromptStrings` set_param is gone).
- Mask icon drawing restrictions: see "Block appearance conventions" above.
- `simscape_sources/SNS_lib.slx` is an older generated library kept for reference;
  the canonical build target is `SNS_Library.slx` via `sns_build_library.m`.
