# SNS_Simscape — Synthetic Nervous System + Knee Simscape build (Sept 2026)

Built by ZCode session 2026-09-08; visuals + URDF route + function-subnetwork
figures added 2026-09-09. Goal: wire the knee SNS circuit to a Simscape model of
`09_BA_003.SLDASM`, with an Animatlab/SNS-toolbox-style neuron library drawn in
the journal diagram language Ben specified.

## What runs today (all in this folder, MATLAB R2025b)

| File | What it is |
|---|---|
| `sns_build_library.m` | builds **`SNS_Library.slx`** — 7 masked blocks with diagram-language icons |
| `sns_build_demo.m` | builds **`KneeReflexDemo.slx`** — reflex circuit + 1-DOF knee plant |
| `sns_run_demo.m` | simulates 5 s, saves `sns_demo_results.png` + `.mat` |
| `sns_draw_circuit.m` | **journal figure**: redraws the demo circuit reader-facing (`figures/KneeReflex_circuit.{png,pdf,svg}`) |
| `sns_function_subnetworks.m` | **journal figure**: 6 arithmetic subnetwork panels (`figures/SNS_function_subnetworks.{png,pdf,svg}`) |
| `snsfig.m` | drawing primitives shared by both figure scripts (neuron/triangle/dot/muscle/box/arrow) |
| `sns_export_diagram.m` | prints Simulink models to `figures/<model>_simulink.{png,pdf}` (300 dpi; optional auto-arrange) |
| `sns_urdf_smoke.m` | proves the URDF → `smimport` pipeline works on this license (uses `urdf_smoke/knee_smoke.urdf`) |
| `import_simscape_when_ready.m` | one-command `smimport` of the real CAD (URDF **or** Multibody-Link XML) |
| `SNS_Library.slx`, `KneeReflexDemo.slx` | generated models (open in Simulink, blocks are double-click editable) |
| `sns_demo_results.png`, `sns_demo_results.mat` | latest run output |
| `figures\` | all generated figures (PNG for review, PDF/SVG vector for LaTeX) |

Run order: `sns_build_library` → `sns_build_demo` → `sns_run_demo` →
`sns_draw_circuit` → `sns_function_subnetworks` → `sns_export_diagram`.

## Block appearance conventions (Ben, 2026-09-09)

Diagram language follows **Szczecinski et al. 2017 Fig. 2** (the functional
subnetwork paper), Rybak/Shevtsova CPG diagrams, and Animatlab:

| Element | Icon |
|---|---|
| neuron (non-spiking RC) | open circle, black edge, white fill, "NS" |
| spiking LIF neuron | open circle with spike glyph |
| Ia / Ib afferent | light-gray circle labeled "Ia" / "Ib" |
| muscle (activation, BPA force) | light-green fusiform ellipse |
| **EXCITATORY connection** | **white triangle, black edges — tip points back toward the presynaptic side, flat base at the postsynaptic side (inverted per Ben 2026-09-09)** (+ light-orange backdrop) |
| **INHIBITORY connection** | **solid black circle** (+ light-blue backdrop) |

- The `NonSpikingSynapse` icon picks its marker **automatically from the sign of
  `Esyn`** (`Esyn >= 0` → triangle, `< 0` → black dot), so the icon always tells
  the truth about the connection. An un-evaluable `Esyn` expression draws "E?".
- Shape is the primary code — figures stay readable in grayscale and for
  colorblind readers. Tints are the redundant cue, from the **Okabe-Ito
  CVD-safe palette** (orange = excitatory, blue = inhibitory, green = muscle).
- Keep every masked block **square** (width == height): mask icons autoscale to
  the block rectangle, so non-square blocks turn circles into ellipses.
- Simulink gotcha (cost us a debug cycle): mask drawing commands accept
  **numbers only** — no LineSpec strings (`'k-'`), no name-value pairs
  (`'LineWidth'`,…), and `color()` takes a color NAME (`color('black')`), not
  RGB. Fill with `patch(...,[r g b])`, set edge color with `color('black')`
  before `plot`, center labels with `disp()`.

## The blocks (SNS_Library.slx)

Neuron = **non-spiking leaky integrate-and-fire, literally an RC membrane**:
`Cm·dV/dt = Gm·(Vrest − V) + Σ Isyn`, τ_m = Cm/Gm (nF/µS = ms).
Outputs membrane `V [mV]` and normalized drive `S = clip((V−Thr)/Slope, 0, 1)`.

Synapse (Animatlab + SNS-toolbox convention):
`Isyn = gmax·Sat(Vpre)·(Esyn − Vpost) [nA]`, presyn saturation
`Sat = clip((Vpre−ThrPre)/SlopePre, 0, 1)`.
**Excitatory = Esyn 0 mV, inhibitory = Esyn −72 mV** — an E vs I connection
differs only by the Esyn mask value (and now by the icon marker).

- `NonSpikingNeuron` — RC membrane (defaults Vrest −52 mV, Gm 0.1 µS, Cm 5 nF ⇒ τ=50 ms; SNS-toolbox conventions)
- `NonSpikingSynapse` — E/I chemical synapse with auto E/I icon (ThrPre −45 mV, SlopePre 0.5/mV ⇒ off at rest, graded above)
- `SpikingLIFNeuron` — spiking LIF with threshold reset (Animatlab-style spiking neuron)
- `IaMuscleSpindle` — stretch + velocity afferent → current (peak 10 nA)
- `IbGolgiTendon` — force afferent → current
- `MuscleActivation` — first-order activation dynamics (τ_act 50 ms)
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

Antagonist BPA pair on a reduced-order knee (`I·θ̈ = T_flex − T_ext + T_load − b·θ̇ − K(θ−θ0)`),
with sensory neurons between afferents and synapses:

- Ia(ext) → **Exc** → MN_ext (stretch reflex)
- Ib(ext) → **Inh** → MN_ext (autogenic inhibition)
- Ia(flex) → **Inh** → MN_ext, Ia(ext) → **Inh** → MN_flex (reciprocal inhibition)
- Ia(flex) → **Exc** → MN_flex, Ib(flex) → **Inh** → MN_flex

Demo behavior after correcting the BPA force block to apply activation once
(`F = Fmax*A*max(0, epsMax - strain)`): the knee rises from 15° and reaches
~42.7° against the 0.5 N·m load with graded co-contraction at the final sample
(A_ext≈0.34, A_flex≈0.73). The antagonist loop shows a ~9 Hz alternating
model limit cycle. Its physiological significance has not been evaluated. Tune:
`gmax` (loop gain), SN `Cm` (loop delay), `b_knee`, `tauAct`.
Plant params (inertia, moment arm, Fmax) are placeholder rig estimates — replace
from CAD mass properties and the Xi-corrected BPA predictions.

## CAD → Simscape: the URDF route (Ben's choice, 2026-09-09)

`smimport` consumes URDF natively, so the pipeline is
**SolidWorks → (sw2urdf add-in) → .urdf → `smimport`** — no MathWorks CAD
plugin needed. **The Simscape Multibody "license blocker" was false**: the
license reports `license('test','Simscape_Multibody') = 0` but carries the
product under the legacy feature name `SimMechanics` (= 1), and R2025b
`smimport` runs — proven 2026-09-09 by `sns_urdf_smoke.m` (imported a
femur+tibia revolute-knee URDF, 22 blocks).

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
