# MuJoCo + SNS-Toolbox: gait2392 robot body with BPA muscles

Everything needed to take Ben's `gait2392_robotbody.osim` into MuJoCo and drive
it with SNS neurons and the custom BPA (braided pneumatic actuator) force model.
Set up 2026-09-08.

## 1. Environment (`myoconv` conda env)

One env holds the whole toolchain (MyoConverter, MuJoCo, SNS-Toolbox):

```bash
conda create -n myoconv python=3.10 -y
conda install -n myoconv -c opensim-org opensim=4.5.2 -y
conda run -n myoconv pip install -e C:\Users\Ben\Documents\GitHub\myoconverter
conda run -n myoconv pip install sns-toolbox
```

Version notes (these matter, don't "upgrade"):

- **python 3.10**: myoconverter requires `<3.11`; opensim 4.6 win-64 only ships
  py311+, so we use **opensim 4.5.2** (last win-64 py310 build).
- **numpy 1.25.2 / scipy 1.11.2 / scikit-learn 1.3.0**: opensim's compiled
  extension needs numpy 1.x; newer scipy/sklearn wheels are built against
  numpy 2 and crash on import. This exact trio coexists happily.
- **pyvista 0.41.1 + vtk 9.3.1**: matching pair (a mismatched vtk breaks
  `import pyvista` inside MyoConverter).
- **mujoco 2.3.7**, loguru, lxml, trimesh, networkx, seaborn, fpdf2.
- Always run things with `conda run --live-stream -n myoconv python ...`
  (torch's DLLs need the env's `Library\bin` on PATH; bare python.exe fails).
  `conda run` rejects multi-line `-c` scripts on conda 24.9 — use files.

## 2. Model repair — `gait2392_robot.osim`

`gait2392_robotbody.osim` had two issues (diagnosed 2026-09-08, Ben confirmed):

1. **25 right-side muscle paths are Ben's intentional robot-route edits** (e.g.
   psoas_r via the torso, 9-point ercspn_r). The left side still carried the
   pristine human geometry, so the model was asymmetric.
2. **19 right-side muscles had `appliesForce=false`** (test-time disables).

`repair_robotbody_muscles.py` builds `gait2392_robot.osim`:
mirrors every right-side `GeometryPath` to the left (locations z-negated,
`_r`→`_l` sockets, conditional-point ranges unchanged, moving-pathpoint
x/y splines unchanged and z spline y-values negated — rules verified against
`gait2392_simbody.osim`, whose left side is exactly the mirror of its right),
and re-enables all muscles. Re-run it any time the robotbody file changes:

```bash
python repair_robotbody_muscles.py
```

Original files are never modified. The official human model
(`gait2392_simbody.osim`) is identical to the current
opensim-models/Gait2392_Simbody/gait2392_thelen2003muscle.osim upstream.

## 3. Conversion to MuJoCo — MyoConverter

`convert_to_mujoco.py` converts **both** `gait2392_robot.osim` (robot, Ben's
routes, both legs) and `gait2392_simbody.osim` (pristine human reference) to
MJCF, into `mjc/<model>/`. Bone meshes (.vtp) come from myoconverter's
`Gait2354Simbody/Geometry` example (identical standard OpenSim bone set);
missing decorations (treadmill.vtp) are stripped from a temp copy and a ground
plane is added by the converter.

```bash
# from Solid_Models/OpenSim/Gait2392_Robotbody/  (takes ~30-60 min for both)
conda run --live-stream -n myoconv python convert_to_mujoco.py
```

Final model: `mjc/gait2392_robot/gait2392_robot_cvt3.xml` (validated at default
pose). Load with the bundled keyframe so muscle path constraints are met:

```python
import mujoco
model = mujoco.MjModel.from_xml_path(".../gait2392_robot_cvt3.xml")
data = mujoco.MjData(model)
mujoco.mj_resetDataKeyframe(model, data, 0)
```

Editing muscle routes after conversion = editing `site pos="x y z"` lines in
the MJCF (plain XML, diff-able) — easier than an OpenSim round-trip.

### Local MyoConverter patches (in `C:\Users\Ben\Documents\GitHub\myoconverter`)

Several muscles with Ben's edited routes (semimem_r, glut_med*_r, ercspn_*,
ext_hal_*, intobl_*, extobl_* ...) cannot reach Thelen fiber equilibrium at
extreme joint angles; upstream MyoConverter aborts the whole conversion on the
first failure. Three small patches make it robust (flagged with `[patched]`):

- `optimization/utils/UtilsOpensim.py` — `getMuscleForceMaps` catches
  per-pose equilibrium failures and returns a pose-usable mask.
- `optimization/model_states/OsimMuscleStates.py` — drops failed poses from
  the force maps (numpy import added); muscles whose equilibrium fails at
  every pose (fixed `EQUILIBRIUM_SKIP` set) get placeholder maps without
  calling OpenSim at all (the solver can hang, not just throw).
- `conversion_steps/O2MStep3.py` — placeholder maps skip the force
  optimization (default muscle parameters kept); PSO/diff failures and
  missing `res_opt` no longer crash the run.
- `optimization/utils/UtilsForceOpt.py` — the PSO worker pool gets the model
  **path** instead of the MjModel object: Windows spawn workers cannot pickle
  MjModel, which hung the whole conversion inside `mp.Pool` creation.

Consequence: a handful of human muscles whose routes Ben changed get slightly
degraded Hill-type fits. They are placeholders — the point of this pipeline is
to replace them with BPAs anyway.

Verified results (2026-09-08): both `_cvt3.xml` models load and step cleanly
from their keyframes — `gait2392_robot`: 92 muscle actuators / 92 tendons /
32 bodies; `gait2392_simbody`: 92 / 92 / 46 bodies.

## 4. Custom BPA muscle (`Code/MuJoCo_SNS/`)

| file | purpose |
|---|---|
| `bpa_muscle.py` | The force model, ported from Ben's MATLAB |
| `bpa_mujoco.py` | MuJoCo glue (`BPAMuscleSystem`) |
| `add_bpa_to_mjcf.py` | Append BPA tendon routes + zero-gain actuators to an MJCF |
| `sns_bpa_demo.py` | End-to-end demo: SNS network → BPAs → MuJoCo knee |

`bpa_muscle.py` is an exact port (verified to 9 significant digits against
MATLAB) of:

- `Functions/festo4.m` — normalized force surface (the `sfit` in
  `FestoLookup.mat`): `F = a0*(exp(-a1*rel)-1) + P*exp(-a3*rel^2)`, P
  normalized to 620 kPa, zero above rel=1 and below 0.
- `Functions/maxBPAforce.m` — `Fmax = P*a1*atan(a2*(rest-0.0075)*P)`.
- `Robot_Data/MonoPamDataExplicit_balanceX3.m` — KMAX, Xi0 offset, Xi3 bend
  loss (`delta_L = Xi3*bend*comp^2`), tendon spring rate `Spr()`, and the
  fortz series-stiffness equilibrium (BPA force vs cable/bracket stiffness
  solved with brentq; pass `series_stiffness`, or `xi1/xi2` projected
  compliance via `effective_series_stiffness`, default = tendon rate).

Usage pattern (see `sns_bpa_demo.py`):

```python
muscles = {"knee_flex_bpa": BPAMuscle(name="knee_flex_bpa", diameter=20,
           resting_length=0.155, kmax_length=0.116, fitting_length=0.025,
           bpa_count=2, pressure_max_kpa=620.0)}
system = BPAMuscleSystem(muscles)
system.attach(model)
mujoco.set_mjcb_control(system.cb)     # forces applied inside mj_step
system.set_activation("knee_flex_bpa", 0.8)   # e.g. from an SNS neuron output
mujoco.mj_step(model, data)
```

How it works: each BPA is a MuJoCo `<tendon><spatial>` route with a **zero
gain/bias `<general>` actuator**. MuJoCo applies no force itself but exposes
tendon length/velocity and the exact moment arm (`actuator_moment`);
`BPAMuscleSystem.cb` computes the force from Ben's law each step and writes it
into `qfrc_applied`. The callback owns `qfrc_applied` (zeroes it each step).

MJCF gotchas (mujoco 2.3.7): site-based tendons must be `<spatial>`, not
`<fixed>`; actuators need `biastype="none"` (not `"fixed"`); joints/actuators
with `range`/`ctrlrange` need explicit `limited`/`ctrllimited` or
`<compiler autolimits="true">`; compiler angle default is DEGREES — set
`<compiler angle="radian"/>` or write degree values knowingly.

### SNS-Toolbox side

`sns_bpa_demo.py` builds a two-neuron antagonist network (flexor/extensor
motor neurons with reciprocal inhibition, separate command inputs), compiles
it with the numpy backend (`net.compile(dt, backend="numpy")`), and maps
neuron potentials to BPA activations (V/e_hi clipped to [0,1]). Watch the
units: with the default neuron parameters, `membrane_capacitance=20` gives a
~20 **second** time constant at dt in seconds — the demo uses 0.02 (~20 ms).
Synapses are excitatory by default; pass `reversal_potential=-5` for inhibition.

## 5. Spinal cord network (`spinal/`) — added 2026-09-09

Two-level RG+PF spinal network (SNS-Toolbox, numpy backend) for the
converted `gait2392_simbody` model: 92 MN pools + Ia/II/Ib afferents,
per-leg half-center RG, 4 phase-shifted PF groups per leg, stance-gated Ib
load sharing, descending DRIVE/POSTURE, COM balance inputs, and a solved
standing-posture injection. See `spinal/DESIGN.md` for the architecture,
literature grounding, verified status, and the open-problems list (joint
sign audit, leg-DoF NaNs, pathpoint weld, free balance).

```bash
cd Code/MuJoCo_SNS
D:/Anaconda/envs/myo/python.exe spinal/check_rhythm.py   # rhythm layer alone
D:/Anaconda/envs/myo/python.exe spinal/runner.py         # stand->walk->stand
D:/Anaconda/envs/myo/python.exe spinal/fit_synapses.py   # back-solve scaffold
```

NOTE: on this machine the env is `myo` (`D:/Anaconda/envs/myo`), not
`myoconv` — same stack (py3.10, mujoco 2.3.7, sns-toolbox 1.5.2, opensim
4.4.1). Conda is not on the Git Bash PATH; call the env python directly.

## 6. TODO / next steps

- Map Ben's real BPA configurations (2×20 mm flexor/extensor sets) onto the
  converted gait2392 routes with `add_bpa_to_mjcf.py`, replacing selected
  human muscles per joint.
- Run MyoConverter validation (`validation=True` in `convert_to_mujoco.py`)
  once the BPA replacements are in, if Hill-model torque accuracy matters.
- If a conversion crashes again after muscle edits: delete
  `mjc/<model>/<model>_cvt3.xml` and `Step3_muscleKinetics/` first (stale
  intermediate files are reused otherwise).
