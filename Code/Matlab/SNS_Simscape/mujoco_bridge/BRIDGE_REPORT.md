# MuJoCo ↔ Simulink Bridge — Install & Prove-It Test (EB475WS4)

Task: prove (or disprove) that the official mathworks-robotics/mujoco-simulink-blockset
can drive MUSCLE actuators of our converted gait2392 model and read back muscle
length/velocity/force — before anyone builds on it.

- Date: 2026-09-12
- Machine: EB475WS4 (12 cores, 32 GB, MATLAB R2025a @ D:\Program Files\MATLAB\R2025a, MinGW-w64 gcc 8.1.0 via MATLAB support package, no VS, winget available, no cmake)
- Model: D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml (timestep 0.005 s, 92 muscle actuators, NO sensor section)
- Rule: never modify cvt3.xml; everything new goes in this folder.

## Progress log

- [x] Created working folder + this report.
- [x] Downloaded blockset ZIP (169 KB) -> blockset\mujoco-simulink-blockset-main\ (README, install.m, src/, tools/, blocks/, examples/).
- [x] Read README, install.m, tools\setupBuild.m, tools\build.m, tools\Makefile, src\mj.cpp, src\mj_sfun.cpp.
- [x] Ran install.m with default MJ_VER = **3.3.6** — download + DLL copy + path save all OK
  (log: logs\install_336.log; mujoco.dll 4.4 MB + glfw3.dll into blocks\).
- [x] **3.3.6 CANNOT load our model.** Headless load test with the bundled compile.exe:
  ```
  cd D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody
  ...\blockset\...\lib\win64\mujoco\bin\compile.exe gait2392_simbody_cvt3.xml out.mjb
  -> XML Error: Schema violation: unrecognized attribute: 'collision'
     Element 'option', line 4
  ```
  Our XML has `<option timestep="0.005" collision="predefined"/>` (MyoConverter output for
  mujoco 2.3.7). Decision per plan: switch MJ_VER to 2.3.7 (the version the model was
  converted for; actuator_length/actuator_velocity sensors exist there). The blockset C++
  only uses stable C-API calls (mj_loadXML/mj_makeData/mj_step/mjv_*/mjr_*), nothing
  3.x-only, so compiling against 2.3.7 headers is safe.

### Source-code findings (read BEFORE running anything — answers the mj_step/mj_forward question)

**THE PLANT BLOCK STEPS MUJOCO WITH mj_step — ctrl feeds muscle activation dynamics properly.**

- `src\mj_sfun.cpp`, mdlUpdate (discrete sample hit), lines 432-446:
  ```
  // Step the simulation by one discrete time step. Outputs (sensors and camera) get reflected in the next step
  miTemp->step(uVec);
  ```
  The input port has direct feedthrough = 0 (line 281), i.e. the control signal is read at the
  sample hit and applied during the step (discrete-controller semantics — no algebraic loop).
- `src\mj.cpp`, MujocoModelInstance::step(), lines 189-199:
  ```cpp
  void MujocoModelInstance::step(std::vector<double> u)
  {
      // same memory location will be accessed during gui rendering
      dMutex.lock();
      for (unsigned index = 0; index < u.size(); index++)
      {
          d->ctrl[index] = u[index];
      }
      mj_step(m, d);
      dMutex.unlock();
  }
  ```
  So: input vector -> d->ctrl -> mj_step(m,d). mj_forward is NEVER called in the sim path.
  For dyntype="muscle" actuators this is exactly right: ctrl is the stimulus, activation
  follows first-order dynamics inside mj_step, and force follows activation.

### Interface facts (from mask init + S-function)

- ONE control inport `u`, width = number of actuators (m->nu). Input can be a Simulink Bus
  (one field per actuator, field names = actuator names) or a plain vector [nu x 1] — selected
  by the block's "controlInterfaceType" mask param (mj_maskinit.m lines 51-64). Bus is converted
  to vector inside the block.
- Sensor output is a Simulink Bus: one field per <sensor> in the MJCF, field names = sensor
  names; scalar per dim. (If the model has NO sensors, the sensor port is removed/terminated.)
- The block's sample time is auto-read from the MJCF `<option timestep>` via mj_sampletime
  (returns m->opt.timestep) and written into the block mask. So plant Ts = 0.005 s for our model.
- Rendering (a MuJoCo viewer window) is optional per-block ("Local"/"Global"/None) — headless
  runs are supported (RENDERING_NONE does no GLFW work; the render thread still starts but
  does nothing when there are no windows/cameras... it does call glfwTerminate() at the end,
  which is harmless).

### install.m facts

- Default MJ_VER = '3.3.6' (from https://github.com/deepmind/mujoco/releases), GLFW 3.3.7.
  MJ_VER is a local variable at the top of install.m — editable.
- Downloads into blockset\lib\win64\{mujoco,glfw}\, copies mujoco.dll + glfw3.dll into
  blockset\blocks\, addpath's + savepath's blocks/examples/src/includes, and stores build
  paths in MATLAB prefs (ispref 'mujoco').
- Windows quirk to watch: `savepath` writes to matlabroot\toolbox\local\pathdef.m — may fail
  if that file is not writable without admin.

## PARKED (Ben, 2026-09-12 — stopping point for a commit)

Stopped right at the decision to switch MJ_VER 3.3.6 → 2.3.7. NOT done yet:
the 2.3.7 install, the mex compile (gcc 8.1 vs the README's 12.2+ requirement
is still untested — risk #1), and the three pass/fail tests (a: fire one
muscle, b: two clock rates, c: sensor readback).

### To resume (any session)
1. `blockset\` + `blockset.zip` are gitignored. If missing, re-download from
   https://codeload.github.com/mathworks-robotics/mujoco-simulink-blockset/zip/refs/heads/main
   and extract into `blockset\`.
2. Edit `blockset\mujoco-simulink-blockset-main\install.m` → `MJ_VER = '2.3.7'`,
   run install.m from MATLAB in that folder (mind the savepath permission quirk
   above; if pathdef.m is not writable, addpath manually per session).
3. `tools\setupBuild` (selectedCompilerWin="MINGW") → `mex -setup c++` → `build`.
4. Sensor-patched model copy `gait2392_simbody_cvt3_simbridge.xml` (same folder
   as cvt3.xml so `Geometry\` resolves): jointpos+jointvel for the right knee +
   actuator_length / actuator_velocity / actuatorfrc for one vastus actuator.
   Verify those sensor type names against the mujoco 2.3.7 headers first.
5. Run tests a/b/c; write the verdicts here.

Everything above the PARKED line is verified fact (read from source, not assumed).

## RESUMED & PROVEN (2026-09-12, later session) — ALL THREE TESTS PASS

### What was done

1. **2.3.7 install**: install.m MJ_VER edited to '2.3.7', old
   `lib\win64\mujoco` removed first, install re-run via `matlab\run_install.m`
   (now self-locating). mujoco.dll 2.3.7 verified in both `lib\win64\mujoco\bin\`
   and `blocks\` (log `logs\install_237.log`). savepath worked without admin.
2. **Sensor names CORRECTED**: the MJCF elements in 2.3.7 are
   `jointpos / jointvel / actuatorpos / actuatorvel / actuatorfrc`
   (enum mjSENS_ACTUATORPOS etc. in mjmodel.h). The earlier guess
   "actuator_length / actuator_velocity" in this report was WRONG.
   Verified empirically (`matlab\sensor_name_test.py`, myo env mujoco 2.3.7).
3. **Sensor-patched model**: `matlab\make_simbridge_xml.py` regenerates
   `gait2392_simbody_cvt3_simbridge.xml` next to cvt3.xml (idempotent, cvt3.xml
   untouched): knee_r_pos, knee_r_vel, vas_med_r_len/vel/frc
   (actuator `vas_med_r` = index 28 of 92; joint `knee_angle_r`).
4. **Mex compile**: setupBuild.m switched to MINGW, `mex -setup c++` (MinGW64
   gcc 8.1), `build` → all 4 targets OK (logs `setup_mingw.log`, `build_237.log`).
   gcc 8.1 compiles the blockset fine — risk #1 retired.
5. **LOCAL PATCH to upstream `src\mj.cpp` (required!)**: upstream
   `initData()` calls only `mj_makeData`, which ZERO-fills qpos. For our model
   qpos0 = all zeros and the standing pose exists ONLY as keyframe 0
   (pelvis_ty = 0.95 — the vertical axis in the converted model); bare
   mj_makeData starts the model 0.95 m sunk in the floor and mujoco itself
   warns "unstable" by t=0.025 s. Patch: after mj_makeData, call
   `mj_resetDataKeyframe(m, d, 0)` when `m->nkey > 0` else `mj_resetData`.
   Marked with an AARL LOCAL PATCH comment; rebuilt (log `build_237_patch.log`).
6. **Tests** (`matlab\run_bridge_tests.m`, log `logs\test_abc.log`,
   artifacts `logs\gt_const.{npz,mat}`, `logs\bridge_tests_abc.mat`,
   models `matlab\models\bridge_test_{a,b}.slx`):
   - Ground truth: `matlab\bridge_groundtruth.py` (myo env), keyframe-0 reset,
     ctrl = 1.0 on vas_med_r, 80 steps (0.4 s — the UNCONTROLLED model, one
     muscle fired, goes solver-unstable at ~0.445 s, so stay inside that window).
   - **TEST a (fire one muscle): PASS** — frc 0 → −549.9 N (max |frc| 2005 N
     transient), knee 0 → +87.8° extension. ctrl→activation→force path works.
   - **TEST b (two clock rates): PASS** — 0.001 s From-Workspace sine
     (half-rectified 1 Hz on vas_med_r) → Rate Transition → 0.005 s plant;
     clean 5 ms sensor grid (81 samples), frc follows the stimulus
     (max 143 N), knee range 138.7°.
   - **TEST c (sensor readback vs Python): PASS, BIT-EXACT** — all 5 sensors
     max abs deviation 0.00e+00 vs the Python 2.3.7 ground truth over 80 steps,
     alignment shift 0 (Simulink sample at 0.005k = Python step k; Simulink's
     t=0 row is the pre-step zero state).

### Simulink plumbing gotchas learned (for Part 3 / future wiring)

- `xmlFileRel` mask param fires a callback that resolves it with `which()` —
  pass the ABSOLUTE xml path in BOTH `xmlFileRel` and `xmlFile`.
- `renderingType` 'None' + `rgbOutOption`/`depthOutOption` 'off' gives a
  headless 1-in/1-out block (vector ctrl in, sensor bus out).
- To Workspace refuses the sensor BUS (async-queue quirk in R2025a); use
  port-level DataLogging (`set_param(PortHandles.Outport(1),'DataLogging','on')`)
  and read `out.sigs` (Dataset → element.Values = struct of timeseries).
- Rate Transition block path: `simulink/Signal Attributes/Rate Transition`.
- Plant sample time auto-read from the MJCF (0.005 s); model solver
  FixedStepDiscrete, FixedStep auto (= LCM of rates).

## VERDICT

The bridge WORKS for our pipeline: muscle ctrl in (vector), keyframe start,
sensors out, two-rate SNS-style control compatible, physics identical to the
Python myo-env ground truth (bit-exact, same mujoco 2.3.7). Ready to carry
the SNS network (Part 3) once that exists as a Simulink model.

## EXPERIMENTS (2026-09-13, same machine): first closed-loop SNS x MuJoCo runs

Scripts in `matlab\` (logs in `logs\`, artifacts in `matlab\exp_out\`):

- **E0 `e0_multibody_license.m`**: hand-built Simscape Multibody (World
  Frame + Revolute + Cylindrical Solid, no CAD import) FAILS on THIS
  machine - "a license for Simscape Multibody is not available" fires at
  block-ADD time, so the SimMechanics seat block hits primitives too, not
  just smimport. Ben's two-cylinder-elbow idea is buildable in principle
  (smlib primitives + a force along a frame-to-frame line), but must run
  on the laptop.
- **E1 `e1_sns_reflex_loop.m`** (exp_out\e1_reflex.{mat,png},
  sns_bridge_reflex.slx): first CLOSED-LOOP neuromechanical sim here -
  SNS_Library Ia/II/Ib encoders read the plant's vas_med_r sensors,
  afferent neurons with the network's tuned taus, synapses with the tuned
  conductances (0.6/0.4/0.35 uS) onto MN_vas_med_r, S -> ctrl, 92-wide
  vector via a MATLAB Function block. Result (reflex OFF vs ON, 1.5 s,
  ode45): OFF explodes after the body collapses (|frc| 43 kN transient,
  knee -253 deg, solver unstable at 0.76 s); ON recruits S 0.40 -> 0.97
  from the stretch, bounds |frc| to 469 N (~Fmax scale) and keeps the knee
  ~80 deg straighter - bounded, damped post-collapse behavior.
- **E2 `e2_cpg_open_loop.m`** (exp_out\e2_cpg.{mat,png}, sns_bridge_cpg.slx):
  the FULL 410-neuron SNS_SpinalNetwork as a Model block driving all 92
  muscles, fixed-step ode1 @ 2 ms with the NEW `simbridge2.xml` (timestep
  0.002 + implicitfast - exactly what runner.py rewrites the MJCF option
  to, so plant and network step at the same 2 ms, single-rate). 5 s clean,
  all drives finite, 96% of muscles recruited, knee flexes to -108 deg
  then extensor tone locks it at 0 deg (vas_med |frc| 4.8 kN transient).
  CAVEAT (finding): at bare DRIVE=2.5 with every other input zero, the
  Simulink-Euler realization settles to a TONIC fixed point - no rhythm -
  while the numpy 2 ms reference oscillates: the network sits near a
  bifurcation at this off-operating-point drive, and the production rhythm
  needs the runner's full input schedule (POSTURE biases, BAL terms,
  afferent/MOD gating). Open-loop DRIVE-only is NOT sufficient - matches
  why the runner carries the schedule.
- `gather_exp_const.py` regenerates exp_const.mat (vas_med_r L0 0.1735 m,
  range [0.164,0.247] m, Fmax 1216 N, MN tau 0.03 s, tuned afferent g's,
  u-port indices); `make_simbridge_xml.py` now emits BOTH xml variants.
