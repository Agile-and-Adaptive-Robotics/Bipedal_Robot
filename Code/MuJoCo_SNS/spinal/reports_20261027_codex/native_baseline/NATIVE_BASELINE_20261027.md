# Native Di Russo baseline reproduced on easteregg2 — 2026-10-27

The recovery report's blocker ("native replay blocked by the absence of the
author's custom SegmentalController runtime on this machine") is CLEARED.
The full original stack was built from source on easteregg2 with a portable,
user-scope toolchain (no admin rights required), and the optimized
`spinal.par` replays the author's walking gait on her original
9-DoF, 18-muscle model. This is the first native replay of this baseline in
this program's records.

## Result of record

`sconecmd.exe -e spinal_controller_main.par -r eval_spinal_native.sto`
from the pristine author-repo clone (git `ccf46fad`), 2026-10-27:

| Quantity | Value |
|---|---|
| objective result | 17.2411 |
| step_velocity | 1.16636 m/s (1.174 m/s free-window fit over t ≥ 1 s) |
| step_count | 16 (10 s) |
| forward progression | 11.57 m in 10 s |
| Gait penalty | 0.14253 (threshold 0.05) |
| Effort (Uchida2016 CoT) | 2550.36 / distance 11.57 m |
| Effort (squared activation CoT) | 0.00376 |
| Ankle/Knee/Hip limit penalties | 0.0028 / 0.0091 / 0 |
| HeadStability Y/X | 0.0156 / 0.1479 |
| pelvis height | 0.856–0.918 m (no fall) |
| stance duty (grf_norm_y > 0.1) | 0.548 L / 0.602 R |
| peak GRF | 1.25 L / 1.41 R body weights |
| knee range | −72.5° to +2.2° (R) |
| pelvis tilt range | −11.9° to −5.7° |
| speed | 0.52–0.54 × real time |

Fresh-process determinism: a second evaluation in a new process reproduces
every metric exactly; the two 1338-column trajectories are byte-identical
after the filename header (diff lines = 0).

Trajectory of record: `eval_spinal_native.sto.sto` — 1001 samples ×
1338 columns: joint states, `phase_left/right`, leg states, per-leg
normalized GRF (support forces), per-muscle excitation/activation/force/
length/velocity, every neuron `.input`/`.output`, every synapse
contribution, and the Ong pelvis-tilt PD channels. These are the regression
targets for any adaptation (gait2392 or otherwise). `L5.gv` is the exported
network graph with per-neuron biases.

## How it was built (all user-scope, no admin)

Portable MSYS2 at `D:\GitHub\_zcode_build_tools\msys64`; everything else
under `D:\GitHub\_zcode_build_tools\`:

- gcc 16.2.0 (Rev6), CMake 4.4.4, Ninja 1.13.2, LAPACK 3.12.1,
  Eigen 5.0.1, Boost 1.92.0, TBB 2023.1.0 (MSYS2 packages)
- Simbody 3.5.4 (tag `Simbody-3.5.4`, 73d463d)
- opensim3-scone `ffa32e4` (2021-01-28, tgeijten's patched OpenSim 3.3)
- pagmo2 v2.20.0 (sconelib requires it even though the scenario's
  CmaOptimizer is the spot-based one)
- simgait/SCONE fork `58ea6876` — built headless: sconelib +
  sconeopensim3 + sconecmd (no Qt, no OSG, no GUI)

Exact versions, revisions, and MD5s: `native_build_manifest.txt`.

### MinGW portability patches (documented, minimal, build-system-level)

1. Simbody `Platform/CMakeLists.txt`: skip bundled pthreads-win32 on MinGW
   (no-op PlatformFiles target, empty platform include/link dirs).
2. Simbody root `CMakeLists.txt`: `find_package(BLAS/LAPACK)` and plain
   `pthread` math libs on MinGW.
3. Simbody 3 × `internal/common.h`: `#ifdef _WIN32` dllimport/export gates
   → `&& !defined(__MINGW32__)` (MinGW exports via default visibility).
4. Simbody link: `-Wl,--export-all-symbols`.
5. opensim3-scone: `_USE_MATH_DEFINES`; `IO.cpp` `_mkdir/_chdir` shim;
   `LoadOpenSimLibrary.h` guarded `windows.h`; MinGW math-libs branch in
   root CMakeLists.
6. SCONE fork: xo `mkdir` 1-arg on MinGW; xo `pattern_match` via
   `PathMatchSpecEx` on MinGW (+ shlwapi link); xo `system_tools.h`
   `<typeinfo>` include; xo `get_application_dir` via `GetModuleFileNameA`
   on MinGW (it previously returned empty on Windows-gcc, which broke
   `FindInstallFolder`); sconelib links `tbb12` on MinGW (pagmo headers
   auto-detect TBB); SimTK import libs copied into the OpenSim sdk/lib.
7. pagmo 2.16's FindTBB predates oneTBB naming → used pagmo v2.20.0 with
   `PAGMO_WITH_TBB=OFF` instead (fork uses only stable pagmo-2.x API).

Rebuild everything: cmake+Ninja per component against the install tree;
run with PATH = SCONE build/bin + opensim3/bin + simbody/bin + pagmo2/bin +
mingw64/bin.

## Native controller audit (from source, `SegmentalController.cpp/h`)

Recorded semantics that any adaptation must match:

- Neurons are leaky integrators `tau·dV/dt = −V + I` with **Euler steps at
  a hardcoded 1 ms** (`euler_advance`, dt = 0.001) — the same dt as the
  scenario's `fixed_control_step_size`. Evaluation is **lazy and recursive**:
  pulling one MN's output at time t advances its whole upstream graph in
  1 ms steps; earlier requests are served from `state_history`.
- Sensors: `I = w0 + s·g(t−delay)` with per-muscle delays from
  `delay/neural_delay_testing.txt` (gastroc/soleus/tib_ant 20 ms;
  vasti/bifemsh 10 ms; proximal muscles and pelvis_tilt 5 ms) through
  SCONE `SensorDelayAdapter`s. Ia uses the Procházka spindle
  (`MuscleIaSpindleSensor`) since `use_prochazka = true`; II is normalized
  fiber length; Ib normalized force.
- Activation functions: sensors ReLU (Ia/II/Ib), interneurons and MNs
  **min_max saturation [0,1]**; state clipping at ±50 afterwards.
- Reflex rules (18): Ia monosynaptic excitation + reciprocal IaIN
  inhibition (with IaIN↔IaIN mutual), II excitation/inhibition via IN_II,
  Ib autogenic inhibition (IN_Ib, with IbIN mutual), extensor common-synapse
  Ib facilitation (IN_Ib+), Renshaw: MN→IN_RC, IN_RC→own MN,
  IN_RC→antagonist IN_RC, and **IN_RC→antagonist IaIN** (the antagonist-IaIN
  Renshaw projection).
- Pattern formation: 5 raised cosines (mu = 0.1..0.9, half-width 0.2) per
  side, connected to all same-side MNs with left/right-shared weights
  `P#.ALL.PF-><muscle>.ALL.MN.w`.
- PhaseOscillator (Aoi 2010 style): phi_dot = omega + delta_omega −
  gain·sin(phi_i − phi_j − pi); **heel-strike reset to 0** via Simbody
  discrete event on delayed (50 ms) leg load vs threshold 0.1; swing
  re-sync to `toe_off_percentage·2pi` when `synchronize_swing` (disabled
  here: `synchronize_swing = 0`, `steps_ignore_phase_locking = 100`);
  delta_omega phase-locks omega to the measured step period.
- MN output drives muscles via `AddInput`; the "balance gates" in the
  native controller are the Ong 2019 Geyer-Herr pelvis-tilt PD
  (hamstrings/glut_max/iliopsoas DofReflexes) in the CompositeController.
- Objective (from `import/measures/measure.scone`): GaitMeasure ×100
  (termination height 0.85, velocity band factors 1.09/1.01),
  Uchida2016 cost-of-transport effort with per-muscle slow-twitch ratios ×1,
  squared-activation CoT ×0.5, ankle/knee/hip DofLimit (squared force
  penalty) ×0.1 each, head-acceleration Y ±4.903 / X ±2.452 m/s² ×0.25.
  CMA: lambda 40, sigma 5, seed 0, max 3000 generations, 8 threads.
- 256 parameters in `spinal.par`: 1 omega, 16 state offsets + 2 phase
  offsets, 9 balance, 45 PF weights, 174 reflex weights, 9 length offsets
  (matches the codex parameter inventory).

## Stock SCONE 2.4.4 and Hyfydy attempts (EB475WS4)

`stock_scone_2.4.4_attempt.log`: stock SCONE 2.4.4 initializes its bundled
OpenSim 3.3-2021-01-28 and OpenSim 4.4, but logs the author's
`SegmentalController`/`PhaseOscillator` blocks as "Warning, unused
properties" — the custom controller does not exist outside the simgait
fork, so the model runs muscle-less under the Ong PD alone and falls at
0.243 s. Hyfydy: "Could not initialize Hyfydy: License key has expired
(2026-10-24)". Consequences: the native scenario can only run on the fork
(OpenSim3 backend, as done here). A Hyfydy-backend run would additionally
require the closed-source sconehfd SDK (not present on any lab machine) and
an active Hyfydy license — both are Ben's call.

## MuJoCo speed layer (the requested 20×-class acceleration)

`dirusso_gait2392/dirusso_fast.py` + `verify_fast_parity.py` +
`dirusso_fast_report.json`:

- NUMBA neural step (controller_numba, parity ≤3.4e-16 vs Python),
  vectorized sensor extraction (one actuator read per step, precomputed
  index/range/fmax arrays), one-time model loads (dirusso_v3 re-parsed the
  full XML per joint name per trial — 17+ parses/trial removed).
- **6.2× per 10 s trial** (22.7 s → 3.7 s); with the 8-worker parallel
  pool ≈ **50× search throughput** over the original path. The sparse
  physics backend (SparseForward) is deliberately NOT used: the 2026-10-27
  recovery report showed it fails full-trajectory parity on this program's
  stacks.
- Component-level parity: sensors computed the two ways are IDENTICAL from
  the same state; MuJoCo is bit-deterministic across model instances. The
  closed-loop residual divergence (onset ≈ step 97 of 1000) is ulp-level
  libm differences between the numba and NumPy kernels amplified by the
  falling dynamics of the current default-parameter trial. Protocol
  (unchanged doctrine): fast path screens; every kept result is re-verified
  through the Python-controller path in a fresh process.
- One bug found and fixed during verification: the fast runner initially
  swapped the receptors' velocity/length arguments; the reference-port
  (`_reference_trial`, verified metric-identical to `dirusso_v3.run_trial`)
  keeps the v3 wiring for comparison. NOTE for the go-forward adaptation:
  `dirusso_v3.run_trial` passes `receptors(rii, ria, rib)` with
  ria=length, rii=velocity, i.e. Ia receives the LENGTH-norm and II the
  VELOCITY-norm — swapped relative to controller.py's own docstring and
  the original's Procházka-Ia spindle. This pre-existed this session and
  is flagged for Ben's ruling before the port continues.

## Prerequisites still open

- gait2392 port of the full native controller/objective (now unblocked:
  regression targets + exact semantics above are the reference).
- Hyfydy backend: sconehfd SDK + active license (Ben).
- The rejected somersault search (`optimize_dirusso.py`) remains retired.
