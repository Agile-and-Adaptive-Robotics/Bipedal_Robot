# AH001 robot with lumbar hinge and payload

The MuJoCo model uses the geometry of `10_AH_001.SLDASM`, both spherical hips,
both four-bar knees, a lumbar hinge, and a rigid 0.25 m spine bar. The lower legs
have no ankle joint or separate foot in this CAD assembly. The source CAD and
Simulink files were not saved or modified.

`ah001_with_torso.xml` is the model to inspect. It compiles in MuJoCo 2.3.7 with
15 independent degrees of freedom: six at the floating root, six at the hips,
two at the knees, and one at the lumbar hinge. Its two knee loop constraints
have combined rank four. The companion `ah001_manual_tree.urdf` has three serial
passive hip axes and omits the cut knee edges; use the MJCF for dynamics.

## Adjustable OpenSim routes

`ah001_with_opensim_routes.xml` adds all 92 routes from the symmetric
`gait2392_robot.osim` model. That source preserves Ben's edited right-side
robot paths and mirrors them to the left. The original file is unchanged.
All actuator interfaces generate zero force until a force model and capacities
are assigned.

Edit each point's body-local `position_m` and `enabled` flag in
`opensim_routes.json`, then run `build_opensim_routes.py` and
`audit_opensim_routes.py` with the machine's `myo` Python. The preparer refuses
to overwrite this specification; its `--output` option creates a separate
fresh transfer. Rebuilding the base CAD model does not overwrite attachment edits.

The transfer aligns hip, knee, and distal-shank landmarks. Pelvis alignment uses
the hip centers and approved lumbar hinge. The knee landmark is the mean of
its four pin centers. Torso guides compress onto the square spine bar, and
foot attachments collapse onto adjustable distal-shank cap points. These are
engineering routes rather than measured attachments. Transfer matrices and
the source hash remain in the specification.

The source contains 336 static, 10 conditional, and 8 moving points, with no
wrap objects. All 354 original definitions remain in
`source_route_inventory.json`. Moving-point locations and conditional membership
are **frozen at the source neutral pose** in this initial routing model. Their
functions and conditions need a later dynamic implementation. Source force and
fiber parameters establish provenance, not robot BPA capacity.

The motion audit covers 65 reachable poses: hip changes of ±30°, knee changes
of ±15°, and lumbar changes of ±20° around the CAD reference. Knee derivatives
follow the closed four-bar mechanism. Length-derived moment arms agree with
the tendon Jacobian within 4.665 × 10⁻¹⁰ m. All interfaces produce zero force
under unit control, and total model mass remains 75 kg.

Twenty ankle and toe routes lose joint leverage when foot joints are removed.
Five additional routes disagree with their source neutral sagittal action:
bilateral `add_mag2`, bilateral `grac`, and left `tfl`. Review these attachments
before neural coupling. Seventeen routes change sign over the sampled ranges;
this alone does not establish a defect or a change from the dynamically routed
source. The full record is `opensim_route_motion_audit.json`.

The OpenSim-derived erector routes require 1354.84 N per side to share the
current neutral holding load. This estimate applies to this routed model;
the six-route engineering template below remains separate. Neither estimate
demonstrates actuator capacity.

For a fresh transfer, run `inspect_source_routes.py`,
`prepare_opensim_routes.py`, `build_opensim_routes.py`, and
`audit_opensim_routes.py`. After attachment edits, only the last two are needed.
The routed preview is `ah001_opensim_routes_preview.png`.

## Mass and geometry

| Item | Value | Basis |
|---|---:|---|
| Complete model mass | 75 kg | User specified |
| CAD lower body | 4.307078 kg | Current solid volumes and assumed material densities |
| Spine length | 0.25 m | User specified |
| Spine cross-section | 20 × 20 mm | Adjustable provisional choice |
| Spine mass | 0.12 kg | Solid Onyx, 1200 kg/m³ |
| Payload mass | 70.572922 kg | 75 kg minus lower body and spine |
| Payload center | 100 mm anterior, 75 mm superior to hinge | User confirmed |
| Payload shape | Sphere, radius 100 mm | Adjustable inertia approximation |
| Lumbar hinge position | Pelvis frame (−160, 0, 100) mm | Approximate upper rear center; user approved adjustable position |
| Lumbar axis | Left lateral axis | Sagittal flexion and extension |
| Lumbar range | ±30° | Adjustable provisional choice |

Coordinates are x anterior, y left, and z superior. The CAD-to-MuJoCo rotation
maps CAD (x, y, z) to (x, −z, y). The selected pelvis hinge location is an
engineering approximation rather than a measured mounting feature.

Printed components use the [Markforged Onyx density](https://web-objects.markforged.com/craft/materials/CompositesV5.2.pdf).
Hardware uses provisional densities of 7850 kg/m³ for steel, 8500 kg/m³ for
brass, and 900 kg/m³ for polypropylene. Steel fasteners are simplified as
homogeneous solids, including nylon locknut inserts. Component ownership at
shared pins is recorded in `../spinal/reports_20261027_codex/robot_geometry_audit.json`.
Actuators, tubing, and hardware absent from the CAD are absorbed into the
balancing payload, not individually represented. Update these inputs when
physical mass and mounting data become available.

Each CAD solid's density-one volume and inertia were extracted directly.
An earlier mass query returned the complete assembly mass for every component;
those earlier `robot_geometry_inventory.json` mass fields are invalid and must
not be used. `robot_mass_geometry.json` is the corrected geometric inventory.
Individual body tensors are rotated and combined with the parallel-axis theorem.
SW2URDF exported nine manually assigned rigid mesh groups. Large meshes are
partitioned into assets below MuJoCo's 200,000-face limit without decimation.

## Trunk routes and load

`trunk_routes.json` defines six adjustable straight-line routes for bilateral
erector spinae, internal oblique, and external oblique attachments. These are
engineering candidates, not measured robot attachments or anatomical muscle-path
replicas. They have zero passive stiffness and damping, and no active force
generator is assigned. A single sagittal hinge cannot represent axial rotation
or lateral bending by the obliques.

The neutral payload needs **69.232 N·m** of extension holding torque. The
provisional erector routes each have a 24.821 mm extension moment arm, requiring
**1394.605 N per side** if they share that load equally. This is a force
requirement, not a demonstrated actuator capacity. The actuator force, stroke,
and mounting design must be checked before controller optimization.

The complete model's neutral center of mass is approximately 100.22 mm anterior
and 52.76 mm superior to the hinge. The requested 75 mm height locates the
payload center; lower-body mass shifts the complete model's center downward.

## Verification and regeneration

Run with `D:\Anaconda\envs\myo\python.exe` on easteregg2:

1. Run `audit_robot_geometry.py` in the recovery report folder after extracting
   the current CAD inventory with `InspectAH001.cs`.
2. Compile and run `ExportAH001.cs` against the installed SW2URDF DLL with the
   intended assembly active. Its first argument is `export_groups.json`; its
   second argument is this directory. It exports geometry without saving CAD.
3. Run `build_ah001.py`, then `verify_torso.py`, and `verify_ah001.py` here.

The checks confirm the 75 kg compiled mass, knee mobility and constraint rank,
the analytic lumbar gravity torque, and tendon moment arms against finite
differences. A two-second mechanical exercise with explicit gravity compensation
had a maximum loop error of 8.341 μm and no numerical warnings. This supported
diagnostic does not demonstrate neural control, balance, or walking.

`ah001_preview.png` shows the compiled model. Generated mesh directories remain
local and are ignored by Git; regenerate them before loading on another machine.
