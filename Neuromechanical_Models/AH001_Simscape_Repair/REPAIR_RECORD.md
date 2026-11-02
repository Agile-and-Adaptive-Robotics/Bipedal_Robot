# AH001 Simscape Repair Record — laptop DESKTOP-5Q16KE9, 2026-11-01

Session executing `Code\MuJoCo_SNS\spinal\fsa_design_20261101\LAPTOP_CAD_SIMSCAPE_REPAIR_HANDOFF.md`.
Design/preflight repair only. No training, optimization, curriculum, controller
attachment, commits, or pushes. Original repo sources untouched (one accidental
overwrite by the exporter was reverted — see §9).

**Update (same day, Ben's directions):** lumbar + spine + payload LANDED (§5b);
foot work root-caused (§5c, decision needed); diameter-expanding BPA model
delivered (§5d); payload-split arithmetic corrected (§7 — the earlier
37.25 kg per side was wrong, correct is 35.7539269872/35.7539269873); material
ruling: printed parts solid Onyx (Markforged), fasteners steel (§7).

## 1. Deliverables in this folder

| Path | What it is |
| --- | --- |
| `working/ah001_repair.slx` | **Repaired mechanism** — full-assembly import with the four-bar knee loops physically restored; compiles and integrates (0 s and 0.5 s) |
| `working/ah001_full_import.slx` | Raw smimport of today's CAD-XML export (pre-surgery snapshot) |
| `cad/export/ah001_full.xml(.xml)` + 21 STEP files | Simscape Multibody Link v7.4 export of the CURRENT `10_AH_001.SLDASM` (geometry self-contained; model paths are absolute to this folder) |
| `cad/export/ah001_full_DataFile.m` | smimport data file of the current assembly revision |
| `working/mdl_humanoid_lower_ah001_imported*.slx` + DataFile | Verbatim copies of the pre-existing repo models (audited, unmodified) |
| `lib/SNS_Library.slx` | Copy of `Code/Matlab/SNS_Simscape/SNS_Library.slx` so the OLD models resolve links without putting the originals' folder on the MATLAB path (name shadowing) |
| `scripts/*.m`, `scripts/*.py` | All probes/audits/build steps, re-runnable |
| `logs/*` | Full evidence: environment probe, joint graphs, adjacency, staged build log, mass ledger, free-fall validation, source SHA-256s (`source_hashes_20261101.txt`) |

## 2. Environment (verified at runtime, not via license('test'))

- MATLAB R2025b (25.2.0.2998904); Simscape Multibody 25.2 and **Simscape Multibody Link 25.2 installed and licensed** (scratch Mechanism-Configuration add-in succeeds; `sm_lib` loads).
- SOLIDWORKS Premium 2025 SP4.1 (33.4.1) via COM; Multibody Link add-in GUID `{2666BDBF-...}`, DLL `C:\Program Files\MATLAB\R2025b\bin\win64\cl_sldwks2sm.dll`, ProgID `sldwks2sm.Sldwks2smApp`.

## 3. Audit of the pre-existing SLX (`mdl_humanoid_lower_ah001_imported.slx`)

- `smiData` is already in its model workspace (identical to the DataFile); the easteregg2
  "unresolved symbols" symptom does not occur on this file. Model compiles and runs with
  `SNS_Library` on the path (use `lib/` here).
- It is a **reduced testbed, not the authoritative mechanism**: 2 lumped `knee_hinge`
  revolutes (no four-bar loops), no lumbar joint (spine welded), placeholder brick feet,
  one 61.09 kg `spine_lump` (no 34.76/34.76 split), no ground/contact, no joint targets.
- Payload wired in: 4x `BPA_20mm` (Rest=0.1281, Kmax=0.1057, TendonL=FittingL=0), 4x
  NonSpikingNeuron (Vrest −52, Gm 0.1, Cm 5, Thr −55, Slope 1), 4x NonSpikingSynapse
  (gmax 1.5, Esyn −77, ThrPre −55). This is the knee-pair demo, not CAD262, and is not
  FSA-qualified. Evidence: `logs/probe3_audit_20261101.txt`.

## 4. CAD export from the real assembly

- `10_AH_001.SLDASM` opens **fully resolved on the laptop**: 76/76 components resolved,
  0 suppressed. 71 components under `00_BipedalRobot_Redesign/`; **4 parts reference
  `Solid_Models/Biomimetics_2022-Knee_Test/Knee assembly/`** (pelvis `02_01_PE_001`,
  both `92997A108` setup studs region, knobs — see `logs/sw_census_20261101.json`).
  Not the same-basename assembly the handoff warns about; recorded as a real dependency.
- Headless export recipe (worked): load add-in via `sw.LoadAddIn(DLL)` →
  `GetAddInObject("sldwks2sm.Sldwks2smApp")` → **`SaveAsSimMechXml(path)`** (the
  no-arg `ExportCADXML` blocks on a hidden modal and is unusable headless). A
  "Export Complete" dialog appears; dismiss by posting `WM_COMMAND/IDOK` to the
  `#32770` window. Export XML: 76 instances, 21 parts, 180 mate primitives
  (78 coincident, 76 concentric, 24 parallel, 2 angle), no error-report file.
- smimport: 75 File Solids, 39 Revolute / 2 Spherical / 2 Cartesian / 4 Cylindrical /
  2 Planar joints + 1 Angle Constraint + 1 floating 6-DOF. smimport warnings: the
  KT↔KB constraint set on BOTH sides "could not be mapped to a joint → rigid connection
  added", plus six dowel-pin welds at the pelvis (legitimately rigid).

## 5. Repair surgery on `ah001_repair.slx` (all staged, compile-checked, logged)

The physical four-bar per knee (matching `ah001_neural.xml`: thigh + KT bracket,
BL back link, FL front link, shank + KB bracket):

1. **Removed the two spurious Cartesian joints** (KT to import net) — these were the
   compile blockers ("base and follower ports rigidly connected") and are not part of
   the intended topology; their dangling transform chains deleted.
2. **Deleted the two smimport-added KT↔KB rigid lines** — legal but wrong: they freeze
   the four-bar. (MuJoCo drops the same redundant axle relation.)
3. **Replaced the two KB↔FL Cylindrical joints with Revolute pins** (in-place
   `replace_block`; both block types rotate about base-frame +X, so the pin axis is
   preserved; the cylindrical's stray slide DOF removed).
4. **Removed the degenerate `Parallel` angle constraint** (both sides rigidly connected
   after loop surgery; parallelism is enforced physically by the four pins).
5. Hip knob/stud spin welds (x09/Revolute6..9 → identity transforms): **attempted and
   reverted** — triggers "Spherical Joint 'Spherical1' has base and follower ports
   rigidly connected" (a hidden rigid path pelvis↔knob2/femur-L exists in the import;
   tracing it is the bounded follow-up F-2). The four spins stay as physically-real
   loose thread rotations with negligible inertia.

Final joint census (Stage F): 1x 6-DOF (floating base), 2x Spherical (hips),
41x Revolute (8 four-bar pins + 33 fastener spins), 2x Cylindrical + 2x Planar
(root fastener spins). Structural DOF: 6 (base) + 3+3 (hips) + 1+1 (knees) = 14,
matching the MuJoCo legs; +3 lumbar when the torso spec lands (17 total per handoff).

## 6. Validation performed

- `sim(StopTime=0)` and `sim(StopTime=0.5)`: both succeed on `ah001_repair.slx`.
- Free fall 0.5 s: base point drops **1.22583 m** (ideal ½·9.80665·0.25 = 1.2258) —
  mass model and gravity [0 −9.80665 0] are active and correct in magnitude; all knee
  pin ranges ≤ 1e-12 deg during the fall (closed loops mutually consistent, no
  assembly drift). Evidence: `logs/freefall_check_20261101.txt`.
- **Open item (orientation)**: the base translation appears along joint-Px, not Py —
  the CAD-root↔world orientation still needs the handoff's signed basis-vector check
  (the import's base chain is World → root Transform (smiData.RigidTransform(169) =
  identity) → 6-DOF → LOCKNUT-14 chain). Fix = one explicit adapter rotation at the
  root Transform, then re-run the free-fall probe (expect the drop in Py). Not
  silently rotated — see F-1.

## 7. Mass ledger — CURRENT assembly differs from the old record

Today's `10_AH_001.SLDASM` is a **material/mass revision of the same geometry**
(identical part CoMs, lighter masses; root constraint changed from `RootGround` to a
floating `SixDof`). Per-part totals (kg, from the 75 File Solid instances):

| part | mass each | instances |
| --- | --- | --- |
| 02_01_PE pelvis | 1.687704409 | 1 |
| 03_01/03_03 FB femur-back | 0.197888 / 0.197885 | 2 |
| 03_02/03_04 FH femur-front | 0.073958 / 0.073959 | 2 |
| 04_01/04_03 KT knee-top | 0.219099 / 0.219091 | 2 |
| 04_02/04_04 KB knee-bottom | 0.098580 | 2 |
| 04_05 BL link | 0.00523 | 2 |
| 04_06 FL link | 0.00567 | 2 |
| 05_01/05_02 TI tibia | 0.192589 | 2 |
| knob 6630T310 | 0.028406 | 2 |
| fasteners (screws/studs/nuts/dowels/inserts) | — | 61 |

- **Total lower-body CAD mass (no feet, no spine/payload): 3.3721460255 kg.**
- **Material ruling (Ben, 2026-11-01): printed parts = solid Onyx (Markforged),
  fasteners = steel.** The old→new mass ratio on every printed part is EXACTLY
  5/6 with identical CoMs (pelvis 2.02525→1.68770, KT 0.26292→0.21910, tibia
  0.23111→0.19259), while every steel fastener and polypropylene knob is
  byte-identical — i.e., a uniform density change on the printed parts only,
  consistent with Onyx at 1200 kg/m³ replacing a 1440 kg/m³ assignment (1200/1440
  = 5/6 exactly). Geometry unchanged. (Per-part material-name read via COM is
  pending — the call needs VARIANT-arg handling; the ratio evidence is
  definitive for "material change, not geometry".)
- Old record 5.3566638205 kg (incl. feet) belongs to the older, heavier revision and
  must not be mixed with this one.
- Per `torso_spec.json`'s mass rule with TODAY's assembly:
  payload = 75 − 3.3721460255 − 0.12 (spine) = **71.5078539745 kg**, split
  **35.7539269872 / 35.7539269873** at torso/pelvis — LANDED in the model
  (§5b); recomputed again once real feet are attached. (Sanity: the old record
  satisfies the same rule: 5.3566638205 + 0.12 + 69.5233361795 = 75.0000000000.)
  Full per-instance table: `logs/mass_ledger_20261101.txt`.

## 5b. Torso additions (landed 2026-11-01, Ben's direction)

`working/ah001_repair.slx` now carries: pelvis-datum port (branched off the
pelvis Solid's frame net), lumbar **Spherical Joint** at pelvis [-0.16, 0.10, 0] m
with DampingCoefficient = 2 (spec 2 N·m·s/rad), 20 mm square 0.25 m spine bar
(density 1200 → 0.12 kg), torso payload sphere 35.7539269872 kg at [0.10,
0.075, 0] m from the hinge, pelvis payload sphere 35.7539269873 kg at
[0.10, 0, 0] m from the hinge. **Compiles; total mass exactly 75.0000000000 kg**
(`logs/add_torso_20261101.txt`). DOF now 17 (14 legs + 3 lumbar); target with
articulated feet is 23. Frame assumption documented: pelvis datum axes
anatomical per the 2026-10-27 ruling.

## 5c. Feet (root-caused; needs Ben's decision)

Ben's spec: Tylice's foot (`Foot Ankle.SLDASM`) with a **2-DOF
subtalar-calcaneus joint + 1-DOF MTP joint per foot** — welded feet rejected.
Done so far: CAD-XML export (`cad/export_foot/`, error report notes only the
two cam-lock constraints dropped), smimport, and full topology decode:
- subtalar-calcaneus = native `Revolute` (Ankle Interior ↔ Bottom Joint) plus
  `Planar2` (Interior ↔ Exterior) which must be a revolute (mis-mapped pin);
- MTP = the foot_axle_1 hinge (`Revolute2` main-plate↔axle + `Cylindrical`
  front-foot↔axle → revolute);
- native extras: mid breaks (Planar/Planar1/Planar3, Cylindrical1/2, Revolute1/3),
  five independent toe revolutes.
- **Blocker**: the foot's CAD mate web is self-tangled — parts are held by BOTH
  face-mates and pin-mates (e.g., the toe hardware rigidly bridges front and
  mid foot), so the RAW import already contains rigid cycles and reports a
  real "position violation" at assembly. Welding to Ben's 3-DOF spec closes
  more cycles; downgrading the diagnostics would hide a genuine conflict and
  was NOT done. Evidence: `logs/repair_foot_20261101.txt`.
- **Two fix paths (Ben to choose)**: (A) clean the mate set in a COPY of Ty's
  foot assembly in SOLIDWORKS (remove the redundant face-mates around the
  ankle/toe hardware), re-export — recommended, gives true CAD joint axes;
  (B) hand-build a 4-body foot (housing/spider/rearfoot/front, toes welded to
  front) from the XML instance transforms + native joint frames — no CAD edits,
  but provisional axis provenance.

## 5d. Diameter-expanding BPA (delivered 2026-11-01, Ben's direction)

`scripts/bpa_braid.m` — braid law with inextensible braid, N turns:
`b² = L² + (πN·D)²` → **D(L) = √(b² − L²)/πN** (diameter grows as it
contracts) and virtual-work force `F = P/(4πN²)(3L² − b²)`, tension-only,
active window Lmin = L0(1−Kmax) ≤ L ≤ L0. Self-test PASSED at 300 kPa with
the vas_med-like config (L0=0.2668 m, D0=20 mm, Kmax=0.25, N=1.5 ASSUMED —
replace N with the measured braid turn count): **D 20.00 → 42.45 mm and
F 450.7 → 135.3 N across 0–25% contraction** (constant-P force decreasing
along contraction is the correct McKibben shape). Mounting pattern for the
robot = two frame ports + Internal Force + Transform Sensor + this law
(`scripts/make_bpa_braid_demo.m` holds the wiring; its PS-signal plumbing
needs one interactive port-check pass in the GUI — Transform Sensor outputs
are PS ports, not Simulink outports).

## 8. Honest gate table (handoff acceptance gates)

| Gate | Status | Evidence |
| --- | --- | --- |
| Skills/AGENTS read; hostname; inventory; no source overwrites | PASS | §1, §9 |
| smiData resolution; import verified vs CAD revision | PASS (with revision finding §7) | §3, §4 |
| Anatomical basis / gravity magnitude | PARTIAL — gravity set and validated in magnitude; CAD↔world axis adapter unresolved (F-1) | §6 |
| Hip/lumbar ball axes; all 3 lumbar DOFs; both four-bar closures | four-bars PASS; hips 3-DOF PASS (4 thread spins documented); lumbar NOT yet added (F-3) | §5 |
| Loop residuals / constraint rank / distal-leg ordering / heel-toe grounding / heading / protected flexion | OPEN — needs orientation fix + feet + contact (F-1, F-4, F-5) | §6 |
| Native force/moment/acceleration balance under measured apparatus | OPEN — fixture not built (F-5); free-fall dynamics validated as the base | §6 |
| ideal54 / mixed2-BPA / future all-BPA separation | PASS as documentation (old model payload catalogued; nothing re-labeled) | §3 |
| All-family FSA/scaling/pathway preflight before training | NOT STARTED (out of scope this session; handoff ordering preserved) | — |

## 9. Repo hygiene

- The Multibody Link exporter also wrote STEP files next to the source parts in
  `Solid_Models/Biomimetics_2022-Knee_Test/Knee assembly/` (overwriting 20 tracked
  files from a previous export) and dropped a stray `10_AH_002.xml` there. Both were
  reverted/removed; `git status` is clean except this new untracked folder.
- Nothing to commit from this session unless Ben wants the repair folder in history.
  `working/` copies of the old SLX files are exact copies (hashes in
  `logs/source_hashes_20261101.txt`), so nothing of the old models is lost.

## 10. Bounded follow-ups (in order)

- **F-1 Orientation**: signed basis-vector check of the base chain; set the explicit
  root adapter rotation; re-run the free-fall probe (drop must appear in Py).
- **F-2 Hidden rigid path**: trace pelvis↔knob2↔femur-L rigidity that blocks the
  hip-spin welds (then welds can be retried; expect clean 3-DOF hips).
- **F-3 Torso spec**: lumbar ball at pelvis [-0.16, 0.10, 0] m, 0.25 m spine,
  37.2539/37.2539 kg payload split at the confirmed pelvis-datum frame (§7 numbers).
- **F-4 Feet**: export/attach `Foot Ankle.SLDASM` at a measured distal-tibia mating
  datum (not the mesh-derived center); recompute the payload ledger with feet.
- **F-5 Fixture**: ground contact under both soles, measured apparatus (87.34% BW
  upward support at pelvis COM + guides + rotational guidance, each reported
  separately), native reactions/contacts/accelerations — the equilibrium gate.

## 6. Session 2 (2026-11-02, Ben's directions: do (A); tibia 002; BPA calibration)

**BPA calibration — LANDED.** `scripts/bpa_braid.m` now reproduces Ben's rated
endpoints exactly at 620 kPa (self-test `logs/bpa_braid_caltest_20261102.txt`):
10 mm → 17% contraction, OD 22.5 mm, F 200 N; 20 mm → 25%, OD 38 mm, F 192 N;
40 mm → 32%, OD 75 mm, F 201 N. Naming convention recorded: BPA size =
UNINFLATED diameter; endcaps protrude 22.5 mm per end (mount-to-mount =
L + 45 mm, exposed as `out.mount_length`). The braid-turn count N is solved
per size from the endpoints — replace with measured counts when available.
Model prediction worth checking against Ben's force data: with endpoint
calibration, full-contraction force comes out ~equal (~200 N) for ALL three
sizes at 620 kPa; if the real actuators differ there, add the Festo
coefficients to the force term.

**Foot (A) — toolchain proven, execution staged.** Full SW mate enumeration
of Tylice's foot now works (65 mates; Lock1/Lock2 confirmed type MateLock =
rigid locks of the axle-retention rings — benign in export). Working COM
route recorded in the skills: in-process makepy by ProgID + wrap the raw
dispatch (`M.ModelDoc2(doc._oleobj_)`, features need `M.IFeature(...)`,
assembly ops need `M.AssemblyDoc(...)`); mates are SUB-features of the Mates
folder (GetFirstSubFeature/GetNextSubFeature). Next: map mates to components
(IMateFeatureData wrap), build the suppression list for the redundant
face-mates around the ankle/toe hardware, apply on a copy, re-export.

**Tibia shortening — measurement result CONTRADICTS the premise; ask Ben.**
OpenSim Gait2392 neutral hip→sole (exact, from the osim offset chain + mesh
bbox + converted-model sole): **0.9222 m** (femur 0.4303 + shank 0.4300 +
subtalar 0.0420 + sole 0.0199; ankle height 6.2 cm). Robot current, from the
audit boxes (pose-independent rigid lengths; knee ref = KT/KB box midpoint,
box-corner distal ends, ±20 mm): femur ≈ 0.456 m, tibia ≈ 0.379 m, foot
ankle-joint→sole ≈ 0.062 m → **total ≈ 0.897 m**. The robot appears ~25 mm
SHORT of OpenSim already, before shortening anything. The discrepancy is
within the approximation uncertainty (knee reference choice, foot mount
offset), so no tibia 002 was cut this session — Ben's rule: ask, don't infer.
Needed for a confident Δ: exact four-bar pin axes from smiData (offline) +
the foot housing-top→ankle-axis offset (one SW measurement or per-part boxes
via the Body2 route). The converted Gait2392 mjc skeleton is degenerate at
the knee (MyoConverter collapsed the femur offset) — noted; the osim+mesh
numbers above are the authoritative ones.

## 7. CORRECTION — OpenSim target measured Ben's way (2026-11-02, session 2)

Ben's method (sole contact in the subtalar frame, located relative to the hip)
applied with OpenSim 4.5 Python (`myoconv` env, `gait2392_simbody.osim`,
neutral = default coordinates, OpenSim's own assembly):
- hip center y = 0.88390 (pelvis 0.95 − 0.0661)
- knee (tibia origin) 0.39582 BELOW the hip — **my earlier 0.4303 was wrong**:
  it came from the femur mesh bbox, which includes ~3.4 cm of condyles below
  the knee axis. OpenSim's knee is a joint at 0.3958, not the mesh end.
- shank (knee→ankle) 0.43000, subtalar offset −0.04195 (chain verified:
  all frames axis-aligned at neutral).
- calcaneus origin y = 0.01613; sole point (converted-geometry sole offset
  in the calcn frame, y −0.01990) lands at y = −0.00377 — the model STANDS
  with the sole on the ground plane by construction.

**Corrected target: hip→sole at neutral = 0.8877 m (≈ hip height 0.884 m,
sole on ground).** This supersedes §6's 0.9222. Consequence: robot current
estimate ≈ 0.897 ± 0.02 vs target 0.8877 ± 0.004 → the robot is within noise
of the OpenSim height (Δ = −13 ± 20 mm, i.e., anything between "grow 33 mm"
and "cut 7 mm"). **No tibia 002 should be cut or grown on this basis** — the
robot-side boundaries (four-bar pin axes from smiData; foot housing
top-face→ankle-axis offset) must be measured first; both are queued with the
foot work. Ben's standing directives recorded: tibia MAY grow (002 naming
rule: increment the trailing number); muscle attachment frames must sit at
the proper foot locations (extract the 54-route foot-side sites and carry
them onto the articulated foot's rear/fore segments when attaching).

## 8. FEET ATTACHED — 23 DOF, ISB frame (2026-11-02, overnight assignment, done)

`working/ah001_repair.slx` now carries BOTH articulated feet per Ben's spec
(after the (A)-mate-cleanup path was blocked by the tangled undergrad foot
mates, the feet were built directly from the foot CAD export data — segment
bodies + the CAD's own joint interface frames — the same surgical pattern
that fixed the knees):

- **Per foot: 2-DOF subtalar-calcaneus** (Exterior–rev–Interior–rev–Bottom
  Joint, using the CAD's own mate frames T26/T25 and T12/T11) **+ 1-DOF MTP**
  (main plate–rev–axle1–rev–front foot, T39/T40 + T21/T22). Toes rigid on the
  forefoot, mid foot rigid on the rearfoot. **NO welded feet anywhere.**
- Census: 49 revolute (8 four-bar pins + 8 foot pins + 33 fastener spins),
  3 spherical (2 hips + lumbar), 6-DOF base, hardware spins. **DOF = 23**
  (6 + 6 hips + 2 knees + 3 lumbar + 6 feet).
- **World frame = ISB/OpenSim (x anterior, y superior, z right), VERIFIED**:
  free fall 0.5 s reads Px≈0, Py=−1.22583 m (ideal −1.2258), Pz≈0. Lesson
  learned the hard way: the CAD→world rotation must sit AFTER the 6-DOF
  follower (T_world_cad block) — rotating the root Transform rotates robot
  and joint frame together and changes nothing.
- Anatomical basis of the leg CAD (KT/KB are TOP/BOTTOM brackets per Ben):
  superior = mean tibia axis; lateral = hip-ball line; anterior = sup × lat,
  sign-checked against the four-bar FRONT link (FL) being anterior. Foot CAD
  frame measured from instance map (x width, y up, z anterior-toes — the
  undergrad design does NOT follow the x-forward convention, so everything
  was measured, not assumed).
- Mass: feet 0.3757169027 kg each; payload spheres auto-adjusted to
  35.3782100846 kg per side; **total exactly 75.0000000000 kg**. Compiles;
  integrates 0.5 s cleanly; knee and foot pins static in free fall (loops
  consistent).
- PROVISIONAL (documented, single-block fixes): foot mount position at the
  tibia part datum (no clearance work per Ben's direction); right foot uses
  the LEFT-foot CAD un-mirrored; mount height along tibia to be set when the
  true housing datum is measured (ties to the tibia 002 decision).
