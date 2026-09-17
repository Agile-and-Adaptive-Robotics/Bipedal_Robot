# 09_BA_003 → Simscape Multibody re-export plan (2026-09-16, laptop session)

## Where things stand

- The committed sw2urdf export (`09_BA_003.URDF\urdf\09_BA_003.URDF.urdf`,
  Sep 10) is a **1-link skeleton**: single link `base_link` (anchored at
  component `04_02_KB_R_003-1`, "Tibia Coordinate System") with **zero
  joints**. It proves the pipe (smimport consumed it fine — see below) but is
  not a usable multibody model. A full re-export with links + joints defined
  in the exporter wizard is still needed.
- The MATLAB side is ready and debugged on this machine (R2025b):
  `import_simscape_when_ready.m` now handles the sw2urdf FOLDER structure
  (imports a sanitized temp copy — smimport derives the model name from the
  file name and chokes on the dot in `09_BA_003.URDF`), confirms gravity via
  the `GravityVector` mask param, and saves
  `SNS_Simscape\mdl_knee_rig_import_tmp_imported.slx`.
- License facts (R2025b, this laptop): `license('test','Simscape_Multibody')`
  reads 0 but the legacy `SimMechanics` feature carries it. `smimport` works
  AND hand-built Simscape Multibody models work (`e0_multibody_license.m`
  PASSES here — it fails at block-ADD time on EB475WS4). Ben's two-cylinder
  elbow plant can be built on THIS machine.

## Assembly facts (from README_SNS_Simscape.md, prior session)

- 15 components; ground = fixed `04_02_KB_R_003` (knee bracket).
- 4 hinge mates; `Knee Angle` planar-angle mate; PathMates 1/2/3/5 on
  `Sketch_20mm_route`; an existing `RotaryMotor1` motion motor.
- Knee axis ≈ CAD Z; tibia points −Y (from import_simscape_when_ready.m).
- Coordinate systems of record: "Tibia Coordinate System", "Theta1" (and the
  tibia bracket CS the current export used as base frame).

## sw2urdf wizard recipe (SolidWorks 2025 SP4.1, add-in 1.6.1)

1. Tools > Export > "Export as URDF…" (add-in already installed/enabled).
2. **Links**: group components into rigid links. Suggested partition:
   - `base_link` → femur-side bracket group (everything fixed to the ground
     bracket `04_02_KB_R_003`),
   - `tibia` → tibia-side group (TI_R parts),
   - optionally separate the moving bracket as its own link if you want the
     bracket compliance explicit later.
3. **Joints**: one revolute joint at the knee axis (name it `knee` — the
   importer notes say to keep the name findable). Use the existing hinge-mate
   axes as the joint axes; set limits (use the same range as the MuJoCo knee
   if you want parity, [-120,+10] deg conventionally flexion-negative there —
   note sw2urdf limits are URDF lower/upper on the chosen axis sign).
4. Export to the same folder (it will overwrite the skeleton
   `09_BA_003.URDF\` package — the old export is reproducible, nothing of
   value is lost; the CSV/meshes regenerate).
5. **Known issue**: sw2urdf 1.6.1 on SW2025 has the vanishing-dialog problem
   (upstream issue #147). If the wizard dialog disappears, retry once; the
   fallback is the Simscape Multibody Link add-in (installed + registered,
   enable in Tools > Add-Ins; exports .xml + STEP, which `smimport` also
   consumes — then point `import_simscape_when_ready.m` at the .xml instead).
6. Tendon parts (`Tendon_Extensor/Flexor.SLDPRT`): **keep them OUT of the
   URDF** — they are cable force paths, not rigid bodies. In Simscape they
   should become force elements (line forces between frames) or stay in the
   reduced-order plant. The 2-minute GUI insert into the assembly
   (README_SNS_Simscape.md) is still worthwhile for clearance visualization.

## After export

- Run `import_simscape_when_ready.m` (edit `inputPath` only if the folder
  name changed). It saves `mdl_..._imported.slx` next to the SNS work.
- Then the manual Simulink steps (import script section 4): add Joint
  Actuation on `knee` driven by `KneeReflexDemo` BPAForce outputs; check the
  Tibia/Theta1 frames as sensor frames.
- Pipeline validation on this machine is already green:
  `mujoco_bridge\matlab\laptop_step1_license_import.m` (E0 + MEX check +
  smimport) and `sns_urdf_smoke.m`.

## Audit one-liner (optional, for mate/CS verification)

Open the assembly in SolidWorks first, then:
`"C:/Users/Ben/.anaconda3/python.exe" <scripts>\sw_session.py status`
and use `sw.ActiveDoc` — `OpenDoc6` from COM is still broken by the pywin32
ByRef issue (see solidworks skill). The audit script draft lives in
`%TEMP%\sw_audit_knee_assembly.py`.
