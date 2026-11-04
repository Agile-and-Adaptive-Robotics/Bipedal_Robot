# AH001 BPA Visualization Report and Postmortem — 2026-11-02

Written by ZCode at Ben's direction, for handoff to any future session (Ben:
"maybe chatgpt can salvage something"). Everything below is stated honestly:
what works and is verified, what is broken and exactly why, and what the
Blender option adds. Absolute paths throughout.

## 1. Verdict up front

The **animated** visualization (`robot_bpa27_knee_sweep.html`) is defective in
five specific ways (§4) and should not be shown to anyone. The **static** one
(`robot_bpa27_kompare.html` + `.png`) is coherent but mediocre. The
**underlying data, kinematics, BPA law, and the Simscape robot model are solid
and verified** (§3) — every number in the knee-sweep table is correct even
though the animation that displays it is broken. A salvage attempt should keep
the data layer, throw away the plotly animation mechanism, and fix the
part-placement layer once, properly (§6) — or move the display to Blender (§7).

## 2. Inventory of what exists on disk

All under `C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\`:

| Path | What it is | State |
| --- | --- | --- |
| `working/ah001_repair.slx` | Full Simscape robot: 23 DOF, four-bar knees, lumbar ball, 2-DOF subtalar-calcaneus + 1-DOF MTP per foot, ISB world frame | **VERIFIED** — compiles, free-fall axis test passes (Py=−1.22583 m at 0.5 s, Px≈Pz≈0), 75.0000000000 kg exactly |
| `cad/stl/*.STL` (21 files) | Structural robot meshes exported from the SLDPRTs (pelvis, femur halves, knee brackets, links, tibiae, all key foot parts) | Good (SW SaveAs3, meter-scale confirmed) |
| `cad/cad_to_mujoco_fit.json` | Rigid CAD→MuJoCo-world fit (Procrustes on hip balls + knee refs) | Hips tight; **27 mm max residual** at the knee refs — the root of the "broken leg" |
| `scripts/render_kompare.py` | Static plotly render | Produced the HTML+PNG; visual quality mediocre (§4.1) |
| `scripts/render_kompare_anim.py` | Animated knee-sweep render | **BROKEN** (§4) — do not reuse its animation mechanism |
| `robot_bpa27_kompare.html` + `.png` | Static view, CAD + 27 routes, hover configs | Coherent, ghostly surfaces |
| `robot_bpa27_knee_sweep.html` (35 MB) | The failed animation | Defective (§4) |
| `scripts/bpa_braid.m` | Calibrated BPA law: D(L)=√(b²−L²)/πN, F=P/(4πN²)(3L²−b²); 10/20/40 mm classes reproduce Ben's 620 kPa endpoints (OD 22.5/38/75 mm at 17/25/32% contraction); endcap 22.5 mm per end; `out.mount_length` | **VERIFIED** (self-test log) |
| `scripts/hip_sole_measure.py`, `logs/*` | Measurements, audits, everything else from the session | See `REPAIR_RECORD.md` §1–8 |

Data inputs (modular, swappable when the optimizer re-runs):
`Code\Matlab\Mesh_Optimization\Results\KompareXXVII\route_points.csv` (128
selected-side points, per-point body frame + meters) and
`route_configurations.csv` (per-site nBPA, diameter, G, lengths). The render
scripts take these paths as arguments — **re-optimized routes drop in without
code changes.**

## 3. What is verified correct (keep all of this)

1. **Route kinematics.** Forward kinematics from `ah001_neural.xml` (body
   tree + hinge axes) places all 27 selected routes in the standing ISB world.
   Sanity checks pass: world bbox = standing figure (sole y≈0, spine routes to
   1.23 m); knee sweep produces physiologically sensible length changes:
   `vas_int_r` −220 mm (23% contraction, D→68.2 mm), `rect_fem_r` −156 mm,
   `sar_r` −75 mm, `grac_r` −68 mm, hamstrings near-zero — exactly the
   agonist/antagonist pattern expected. The D(L) expansion values in the same
   table are computed by the calibrated braid law and are correct.
2. **The Simscape model** (§2 table) is the real engineering deliverable of
   this session and is independent of the broken visuals.
3. **The frame mapping** from survey frames to MuJoCo bodies (pelvis, femur→
   thigh_l, tibia→shank_l, calcn*/toes→foot_l, torso*→lumbar_torso) is
   correct at neutral — verified by the bbox and by hip/knee/ankle anchors
   landing on the CAD brackets.
4. **The 27-route set** matches the neural model's installed 54 actuators
   (27 stems × 2 sides) exactly.

## 4. Defects in the delivered visuals, root-caused

### 4.1 Both renders: surfaces look like point-cloud soup
Plotly `Mesh3d` with `opacity≈0.2–0.3`, no `flatshading`, no lighting, over a
bright background renders dense STLs as translucent fog. The pelvis (largest
mesh, hardest decimation: `tris[::stride]`) reads as a cloud of points. It is
a rendering-fidelity failure, not a geometry failure — the STLs themselves are
good. Fix: `flatshading=True`, opacity 0.9+, a dark template, and honest
decimation (quadric decimation, not stride), or stop using plotly for surfaces.

### 4.2 Animated render: leg breaks apart under rotation
I placed each CAD part at its Procructes-fit neutral pose, then under knee
flexion rotated it by "the delta rotation of its assigned MuJoCo body about
that body's world origin." Three compounding errors: (a) the fit carries a
27 mm residual, which a rigid rotation amplifies visibly at the knee; (b) the
part→body assignments were approximations (KT_L→thigh_l, BL_2→shank_l,
FL_2→thigh_l — the four-bar parts straddle the joint and cannot be parented to
one body this way); (c) the correct pivot is the knee hinge anchor in CAD
coordinates, which I never computed. The correct fix is to drive parts from
the **exact CAD mate frames** (the same smiData interface transforms used to
build the Simscape model — they exist, they are exact, and I already used
them for the feet), not from Procrustes-on-bounding-boxes.

### 4.3 Animated render: muscles don't visibly attach
Route endpoints are correct in world coordinates at neutral (§3.1), but the
bodies they should land on are rendered broken/misplaced (§4.2), so endpoints
float off surfaces. Additionally, survey sub-frames (`calcnA`, `calcnSub`,
`femurB`, `torsoB`, `pelvisShift` — 80 of 128 points) were approximated by
their parent frames because their offset transforms live only in the
easteregg2 survey MATs; that adds a few cm of uncertainty on foot/torso
routes. Disclosed when built; still a real approximation.

### 4.4 Animated render: feet missing
`render_kompare_anim.py` includes only the left-leg `LEG` list — the foot
placement block that exists in the static script was never carried over. The
feet are simply absent from that file. (The feet ARE in the static render and
in the Simscape model.)

### 4.5 Animated render: slider does nothing
My hand-rolled plotly animation is wrong twice over: every frame's traces were
added with `visible=(fi==0)` so later frames are invisible, and the
`go.Frame(traces=...)` indices were computed cumulatively as
`len(fig.data) − (per-frame count)` while the figure kept growing, so they
point at the wrong traces; frames built this way carry no data payload, so
`animate` has nothing to apply. I shipped it without once clicking the slider.
This is entirely on me.

## 5. What I actually did this session (summary)

Verified the environment; exported and re-imported the full leg CAD; repaired
the four-bar knee topology in Simscape (compile-verified); added the lumbar
ball, spine, and split payload (75.0000000000 kg exactly); decoded and
attached both articulated feet per Ben's 2+1 DOF spec (23 DOF total, verified
by free-fall axis test); measured the OpenSim neutral hip→sole target (0.8877
m, Ben's frame method, correcting my own earlier mesh-based error); calibrated
the BPA diameter law to Ben's 620 kPa data; produced the static and (broken)
animated route visualizations from the KompareXXVII package. Session log:
`REPAIR_RECORD.md` §1–8.

## 6. Salvage path for the current (plotly) pipeline

1. Keep: CSV data layer, FK, braid law, STLs, frame map (§3).
2. Fix part placement ONCE, exactly: re-export the leg CAD-XML with per-part
   placement in the assembly frame (already have it — `ah001_full.xml.xml`
   InstanceTree), and for articulation rotate parts by the **exact hinge
   frames** (knee axis from the four-bar mate transforms in
   `ah001_full_DataFile.m`, the same source the Simscape repair used). This
   removes the Procrustes residual and the wrong-body parenting entirely.
3. Add the foot block (copy from the static script).
4. Rendering: `flatshading=True`, opacity ≥0.9, dark template, quadric
   decimation ≤8k faces/part.
5. Animation: either (a) proper plotly frames — one trace per entity, frames
   carry the full data arrays, no `visible` juggling; or (b) abandon in-HTML
   animation and ship 7 static pose PNGs from the same loop (simpler, more
   robust, dissertation-usable).
Estimated effort: half a session, mostly the exact hinge-frame fix.

## 7. Blender assessment (Ben's question)

**Verdict: Blender is the right tool for the visuals; Simscape stays the tool
for physics and quantitative interference. Use both on the same data layer.**

What Blender does that plotly cannot:
- Real surfaces: STLs import as meshes with proper materials, smooth shading,
  and viewport/render lighting. No ghost soup.
- Diameter expansion done properly: each BPA as a curve object; tube radius =
  bevel depth driven per frame by D(L(t)) via shape keys or drivers — smooth,
  visible swelling as the joint flexes.
- Articulation: build an armature from the XML body tree (bones at the hinge
  anchors, joint axes set); skin each STL to its bone by the exact CAD mate
  mapping; bake the knee/hip sweep from the same FK. Inverse-kinematics if
  wanted.
- Interference checking: visual overlap is immediate in 3D; Blender's
  boolean/intersection tests give a crude numeric check too — but the
  engineering-grade numbers (contact forces, clearances in mm over a ROM
  sweep) belong in Simscape (Spatial Contact Force blocks between BPA
  capsules and the CAD solids), which the repaired 23-DOF model is ready for.
- Output: stills at any resolution, turntables, MP4 for the defense deck.

Practical notes for this machine: **Blender 5.3.0 Beta is installed and
verified** (2026-11-02): Ben's manual extraction at
`C:/Users/Ben/Desktop/blender-5.3/` was 6622/6627 files complete (the 5
missing were setuptools test-data, now restored — manifest-verified 6627/6627,
zero size mismatches). Exe: `C:/Users/Ben/Desktop/blender-5.3/blender.exe`;
GUI and headless `blender -b --python-expr` both confirmed working (bpy
5.3.0 Beta). Source zip kept at
`C:/Users/Ben/Downloads/blender-5.3.0-beta+v53.ace25df17a10-windows.amd64-release.zip`.
The whole pipeline is automatable from this exe: one script imports the 21 STLs, builds the
armature, creates 27 curve-BPAs from the CSVs, keys the sweep, applies D(L)
radii, and renders stills/frames. The one piece of new engineering is the
same as §6.2 — the exact part→bone mapping from the CAD mate frames; do that
once and both pipelines consume it.

Modularity (Ben's requirement for the optimizer re-run): everything reads the
KompareXXVII-format CSVs; a re-prompted optimizer emits new CSVs; the Blender
scene and the plotly scene regenerate from arguments. No route data is
hard-coded in any scene.

## 9. Skills updated this session (ZCode_Skills repo, uncommitted)

- `solidworks/SKILL.md`: headless Multibody Link export (SaveAsSimMechXml + WM_COMMAND dialog dismissal, exporter file side-effects), early-bound COM dispatch recipe (in-process makepy + M.IFeature wrapping — the thing that unlocked mate enumeration), SW batch-crash mitigation (chunked export with restarts; gmsh STEP import pathologically slow on this machine).
- `matlab/SKILL.md`: the full Simscape-model-surgery section (ReferenceBlock embedded newlines, physical port walking, replace_block port-preservation, weld-cycle hazards, connection-port datum branching, simlog q/w/a/b, script local-function and cell-iteration traps).

## 8. Recommended order of work

1. Fix the exact part→body mapping from CAD mate frames (unblocks everything).
2. Blender pipeline on the fixed mapping (visual deliverable Ben wants).
3. Simscape BPA capsules + contact blocks (quantitative interference, feeds
   the optimizer re-prompt with real packaging data).
4. Re-render plotly statics from the fixed mapping as the quick-look layer.
