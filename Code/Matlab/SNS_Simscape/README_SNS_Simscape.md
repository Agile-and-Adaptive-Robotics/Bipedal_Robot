# SNS_Simscape — Synthetic Nervous System + Knee Simscape build (Sept 2026)

Built by ZCode session 2026-09-08. Goal: export `09_BA_003.SLDASM` into Simscape,
add BPA/tendon CAD representations, and build an Animatlab/SNS-toolbox-style neuron
library (non-spiking RC neurons, E/I synapses, muscle afferents) wired to a knee model.

## What runs today (all in this folder, MATLAB R2025b)

| File | What it is |
|---|---|
| `sns_build_library.m` | builds **`SNS_Library.slx`** — 7 masked blocks |
| `sns_build_demo.m` | builds **`KneeReflexDemo.slx`** — reflex circuit + 1-DOF knee plant |
| `sns_run_demo.m` | simulates 5 s, saves `sns_demo_results.png` + `.mat` |
| `import_simscape_when_ready.m` | one-command `smimport` for when the blockers below clear |
| `SNS_Library.slx`, `KneeReflexDemo.slx` | generated models (open in Simulink, blocks are double-click editable) |
| `sns_demo_results.png`, `sns_demo_results.mat` | latest run output |

Run order: `sns_build_library` → `sns_build_demo` → `sns_run_demo`.

## The blocks (SNS_Library.slx)

Neuron = **non-spiking leaky integrate-and-fire, literally an RC membrane**:
`Cm·dV/dt = Gm·(Vrest − V) + Σ Isyn`, τ_m = Cm/Gm (nF/µS = ms).
Outputs membrane `V [mV]` and normalized drive `S = clip((V−Thr)/Slope, 0, 1)`.

Synapse (Animatlab + SNS-toolbox convention): `Isyn = gmax·Sat(Vpre)·(Esyn − Vpost) [nA]`,
presyn saturation `Sat = clip((Vpre−ThrPre)/SlopePre, 0, 1)`.
**Excitatory = Esyn 0 mV, inhibitory = Esyn −72 mV** — an E vs I connection differs
only by the Esyn mask value.

- `NonSpikingNeuron` — RC membrane (defaults Vrest −52 mV, Gm 0.1 µS, Cm 5 nF ⇒ τ=50 ms; SNS-toolbox conventions)
- `NonSpikingSynapse` — E/I chemical synapse (ThrPre −45 mV, SlopePre 0.5/mV ⇒ off at rest, graded above)
- `SpikingLIFNeuron` — spiking LIF with threshold reset (Animatlab-style spiking neuron)
- `IaMuscleSpindle` — stretch + velocity afferent → current (peak 10 nA)
- `IbGolgiTendon` — force afferent → current
- `MuscleActivation` — first-order activation dynamics (τ_act 50 ms)
- `BPAForce` — `F = Fmax·A·max(0, epsMax − strain)` placeholder; swap in the Xi-corrected
  MonoPam prediction from Mesh_Optimization when ready

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
Plant params (inertia, moment arm, Fmax) are placeholder rig estimates — replace from
CAD mass properties and the Xi-corrected BPA predictions.

## Blockers hit today (both need Ben / MathWorks account)

1. **Simscape Multibody is not in this machine's MATLAB license**
   (`license('test','Simscape_Multibody') = 0`, LME error −5.2: feature absent from the
   named-user license). `smimport` and simulation will not run until the license adds
   Simscape Multibody (or work moves to a machine that has it).
2. **"Simscape Multibody Link" (CAD exporter) is not installed** — only the C-API
   headers exist under `toolbox\physmod\smlink`. It is a separate MathWorks download
   tied to the same entitlement.

### Unblock steps (when licensed)
1. MATLAB → Add-Ons → Get Add-Ons → search **"Simscape Multibody Link"** → Install.
2. In MATLAB: `cd(fullfile(matlabroot,'toolbox','physmod','smlink','mw')); smlink_linksw`
   (registers the add-in with SolidWorks 2025), restart SolidWorks.
3. SolidWorks → Tools → Add-Ins → check **Simscape Multibody Link** (Start-up).
4. Open `09_BA_003.SLDASM` → Simscape Multibody Link tab → **Export** → writes
   `09_BA_003.xml` + STEP parts next to the assembly.
5. Run `import_simscape_when_ready.m` here.

The Hinge mates (Hinge1/2/5/6) become revolute joints; the "Knee Angle" dimensioned
mate carries the joint state target; `Tibia Coordinate System` / `Theta1 Coordinate
System` come across as frames.

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
