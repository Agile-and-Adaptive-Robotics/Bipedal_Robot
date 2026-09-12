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
