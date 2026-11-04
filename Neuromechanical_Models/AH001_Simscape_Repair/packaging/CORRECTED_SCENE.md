# Corrected packaging scene

Use `output/AH001_Packaging_Corrected.blend`. Earlier files are preserved and are superseded for this visualization.

## Straight-sleeve correction

All BPA bodies are now straight, with their full saved length preserved. Tendons alone follow route guides. Calf sleeves terminate at the first registered foot guide and extend proximally; the proximal leader explicitly shows displacement from the saved origin. This is a display placement, not an approved replacement mount design. Oversize calf actuators can consequently extend above their original origin; their lengths have not been shortened to conceal this.

CalcnA points are converted through their ankle frame into calcnSub before registration. Calf survey coordinates are registered between the knee datum and the shared CAD ankle datum (the source tibia and CAD tibia differ in length). The affine display registration is recorded in `scene_data_corrected.json`. Split sleeve offsets are projected perpendicular to each end section independently. Duplicate endpoint points no longer collapse a sleeve.

`straight_sleeve_audit.csv` records length ranges and first-guide span conflicts. These flags mean the straight body exceeds the first saved guide span, not a certified CAD collision. Flagged route collections carry `CHECK AVAILABLE LENGTH` in Blender. All 84 sleeves pass straightness and length checks at all 25 source poses. Ankle/toe actuators retain their saved neutral-ankle state during the knee-only sweep; no artificial contraction is added.

- All 27 replacement route groups use 20 mm BPAs. The original nBPA counts are retained: 69 parallel channels total. Travel gain is separate and is not interpreted as a BPA multiplier or a physical wheel count. Soleus has two parallel 20 mm BPAs.
- The six knee-crossing routes with saved resting BPA length above 0.50 m use two pinned, end-mounted sections. Each section has half the saved total BPA length, as selected by Ben. This gives 84 sleeve sections across the 69 parallel channels. The central section is Bowden housing/artificial tendon. These are series sections, not four-actuator replacement stacks. The selected split routes are listed in `corrected_configurations.csv`.
- All sleeves use the 20 mm diameter law: uninflated 20 mm, expanding to 38 mm at 25% contraction. Saved contraction is retained as a display input; changing diameter and transmission layout has not been represented as a new force-equilibrium solution. Out-of-calibration saved states are recorded in `corrected_calibration_samples.json`.
- Parallel sleeve centres retain the saved layout direction, with 43 mm spacing (38 mm maximum calibrated OD plus the saved 5 mm bundle margin). This is an explicit preview packaging choice, not a reoptimized mount layout.
- Output tendons share one terminal pin. Bundle offsets move sleeves, not foot insertion coordinates.
- Source calcn and toe frames now share a single registration to the CAD foot. The old route-root folding placed even the same Achilles point at inconsistent world locations. A translation was determined using distal toe extent and dorsal/plantar surfaces, and terminal attachment pins were registered on CAD surfaces. Original source coordinates and each adjustment are recorded in `foot_attachment_registration.csv`. These are visualization registrations, not newly measured hardware attachment points.
- The prior CAD pose contained about 17.20 degrees of right knee flexion and 25.09 degrees of left knee flexion. The new sweep starts at anatomical extension and reaches -120 degrees on both sides. The solver uses the actual femur/tibia long axes, the correct hinge signs for each side, and the four-bar loop closures. `corrected_knee_verification.csv` records all 50 side/pose checks. The ankle orientation is registered to be level at anatomical knee zero, rather than at the old bent reference pose.
- The animation has 289 frames at 24 fps. It moves from 0 to -120 degrees and back. The camera travels 180 degrees about Blender world Z across the entire playback, with linear interpolation.
- The world and viewport background are white. This is a scene setting; no global Blender preferences are changed.

The corrected scene remains a packaging preview. Sheaves, brackets, hardware mounting dimensions, Bowden bend-radius constraints, and new force/pressure solves are not fabricated or inferred from the abstract travel gain. The old erector-spinae route remains available for Ben's separate reverse-pulley redesign.

## Rebuild

Run `prepare_packaging_corrections.py` with the existing myoconv Python, then `build_corrected_blender.py` in Blender. Both scripts derive their paths from the repository. The builder creates a separate scene and writes a new file, preserving the original open model and its user edits. It may also be run from Blender's Python console using the same exec/compile procedure as the original builder.

On this machine, the desktop executable failed in the restricted shell with 0xC0000022. A temporary working copy at `C:\Users\Ben\AppData\Local\Temp\AH001_BlenderRuntime` runs headless successfully; no Windows security settings or file permissions were changed.

No survey MAT/CSV, torque result, optimization output, dissertation source, or original Blender file is changed by this correction pipeline.
