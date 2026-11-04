# Reoptimize all 27 robot muscle groups using the BiPulley transmission scheme

Reoptimize all 27 robot muscle groups using the BiPulley transmission scheme, including actual attachment-location and route optimization, physically represented packaging, and both rigid and compliant design cases. Implement and run the campaign; do not stop after producing a plan or another fixed-route length survey.

**Authoritative geometry and existing code**

Read the repository instructions and check the machine identity first. Use easteregg2 for the full campaign if accessible; use this machine for development and bounded verification.

Start from:

- `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim`
- The corresponding stock human `gait2392_simbody.osim`
- `Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json`
- `Opt_run_BiPulley.m`, `Opt_sanity_BiPulley.m`, and `biPulleySpecsFromOpenSim.m`
- `Code/Matlab/Robot_Data/BiPam_pulley.m`
- `Code/Matlab/Mesh_Optimization/BiPam_pulley2_20261022.m`
- Both Routing27 implementations and the `Opt_run_*` family, including cam-assisted optimization and the associated wrapping, exclusion, and route-transition functions.
- The current dissertation's Chapter 30 and geometry appendix for correction-factor definitions, bracket frames, and stiffness assignments.

Resolve the actual current classes and their dependencies; do not assume "BiPulleyData" is the literal class filename.

The old BiPulley driver explicitly optimized only resting length and tendon length while keeping routes fixed. That does not satisfy this request. Show which decision variables now move origins, insertions, intermediate guides, and tackle locations.

The recent Blender packaging files contain display-only attachment registration, affine adjustments, and guessed placement. Do not use those adjustments as anatomical ground truth or optimized results.

**1. Establish correct geometry before optimization**

Extract the original OpenSim path points, their owning bodies, intermediate frames, conditional/moving-point behavior, wrapping objects, and joint transformations. Preserve the complete frame chain.

Compare the extracted routes directly against OpenSim at neutral and representative joint configurations. Pay particular attention to:

- Vastus intermedius/the combined-vasti actuator: the rendered insertion was too far anterior. Check the actuator-map seed preference and patellar routing rather than substituting an arbitrary anterior attachment.
- TFL: the rendered route travelled too far laterally before changing frames.
- CalcnA, calcnSub, toes, and their associated ankle/subtalar/MTP transformations.
- Erector spinae: initialize from the original OpenSim coordinates. Do not preserve the manually wrapped robotbody loop as an imposed transmission topology; that loop predates the formal multi-parallel, multi-sheave scheme. Inspect the original robotbody path and the corresponding stock OpenSim route to distinguish anatomical coordinates from that legacy artificial loop. Document the seed selected.

Do not freeze moving/conditional points across the entire motion range merely because their default-pose coordinates are convenient. Implement their behavior or explicitly replace it with an optimized mechanical guide.

Do not use bounding-box agreement as proof of coordinate correctness.

**2. Use BiPulley mechanics for all 27 groups**

Support monoarticular and biarticular routes through the same verified transmission framework. A route crossing one joint must not acquire a fictitious second joint.

All replacement BPAs are nominally **20 mm**. Optimize/document the actual number of parallel BPAs and the transmission configuration. Keep these distinct:

- Number of parallel BPAs.
- Number of series sections.
- Number and arrangement of physical sheaves.
- Reeving and travel/force ratio.

A configuration calling for three BPAs and two sheaves means three 20 mm BPAs and two physical sheaves, not twelve BPAs or three four-actuator stacks. Never infer physical wheel count from travel gain alone.

Check the transmission using cable-length closure, force balance, and mechanical work. Apply losses and compliance consistently; do not obtain torque or travel benefits by double-counting the transmission ratio.

For each 20 mm sleeve, enforce and verify:

- Diameter = 20 mm at its resting length.
- Diameter = 38 mm at its calibrated fully contracted length.
- A physically consistent diameter–length relation between those endpoints.
- Resting length excludes end fittings; packaging includes both fittings.

Use the validated force model and length-dependent contraction calibration where available. Do not silently substitute a universal 25% contraction limit for all resting lengths. Flag extrapolation outside the experimental domain.

**3. Optimize attachment locations and routing**

Optimize origins, insertions, appropriate via/guide points, resting lengths, tendon lengths, tackle placement, and the discrete transmission configuration.

Attachment locations must move within bounded, anatomically appropriate surface patches on their assigned bones. Define the surface patches and permitted mount standoff in millimeters. Use surface coordinates or projection so a Cartesian search cannot send a mount far beyond the OpenSim bone surface.

Preserve reasonable anatomical attachment regions and mechanical load paths. Do not permit unlimited movement around a bone simply because a candidate improves torque.

Use the established `Opt_run_*` and cam-assisted methods to search routes around obstacles. Verify wrap-point activation, removal, and collinear transitions. A duplicate or eliminated point must never collapse a BPA, create a false direction, or produce a torque discontinuity.

Use actual robot/CAD hardware for structural packaging obstacles, and OpenSim anatomy for route initialization and anatomical bounds. Register the models through documented frames.

**4. Run two genuine optimization branches**

For every muscle, perform:

1. A rigid-body optimization with the relevant compliance terms disabled.
2. A design-factor optimization using Xi₀–Xi₄, with force–deflection–length equilibrium coupled inside every candidate evaluation.

Do not merely reevaluate the rigid optimum using compliant equations. Reoptimize geometry and sizing in the compliant branch as well. Use matched targets and comparable solver budgets, and preserve both results.

Before running, reconcile the requested Xi₀–Xi₄ notation with the dissertation and implementation. Chapter 30 currently describes χ₀–χ₃, including directional bracket stiffness and wrapping loss. Identify whether Xi₄ already exists elsewhere, its physical meaning, units, bounds, and where it enters the equations. Do not invent a fifth factor or silently rename an existing one. If its definition cannot be resolved, state the precise missing definition and continue the independent geometry, rigid-branch, and packaging work.

Assume **two principal attachment brackets**, one at the origin and one at the insertion. Account for both under their actual loads, including the combined loading from parallel BPAs.

Default to Xi₁ stiffness in all three directions only where the dissertation does not prescribe otherwise. Chapter 30 explicitly gives configuration-dependent anisotropic assignments, so check the applicable table and bracket frames. Report the stiffness matrix used for each bracket, the source, and any deliberate generalization to this campaign.

Track additional tackle mounts and Bowden housing stops explicitly. Do not omit their hardware or reactions under the phrase "two brackets," and do not count the same compliance twice.

**5. Optimize packaging from deep layers outward**

Packaging must participate in optimization, not be checked only after selecting a torque optimum.

Establish a dependency order using anatomical layers, initial bone distance, and competing corridors. Search this machine and easteregg2 for the available Atlas of Anatomy and consult it if present. Do not claim to have consulted it unless you actually did.

Use this provisional regional order as a starting point, then verify it against the models and atlas:

- Trunk/deep hip: psoas, iliacus, erector spinae; internal oblique before external oblique.
- Thigh: vastus intermedius/combined vasti and biceps femoris short head before the superficial crossing groups; then adductor magnus, semimembranosus, biceps femoris long head, rectus femoris, gracilis, sartorius, TFL, and the gluteal replacement group as its local corridors permit.
- Lower leg: tibialis posterior, flexor digitorum, and flexor hallucis first; peroneus brevis before peroneus longus; soleus before medial/lateral gastrocnemius. Place tibialis anterior, extensor digitorum, extensor hallucis, and peroneus tertius in their appropriate anterior/lateral corridors.

This list covers the 27 groups but is not a declaration that every regional list is one universal anatomical stack. Resolve regional dependencies and document the final schedule.

Retain the occupied envelopes of optimized inner layers while placing outer layers. After the first pass, perform coordinated refinement so an early poor choice does not permanently block the remaining muscles.

Model sleeves as straight, pinned pressure vessels. Tendons/Bowden cables follow routed paths. Include fittings, brackets, pins, tackle blocks, sheaves, cable exits, housing stops, and service clearance.

**6. Foot and knee tendon requirements**

No BPA body, end fitting, or tackle block may occupy the foot. Tendons may route through and terminate on the foot.

Determine and enforce a minimum tendon length for every foot-actuating route from the required ankle-to-insertion routing and hardware clearance over its motion range. Do not use an arbitrary token tendon length that allows a sleeve to extend onto the foot.

Likewise, define and enforce the minimum artificial-tendon/Bowden span needed across the knee.

For long actuators extending across the knee, consider pinned proximal and distal BPA sections with the central crossing carried by artificial tendon/Bowden cable. Preserve the distinction between total saved BPA length and individual series-section lengths; initialize equal halves where applicable, and explicitly parameterize any subsequently optimized split.

Package the brackets, housing stops, ferrules, cable exits, articulation, and minimum bend radius throughout the knee sweep. Consider common bicycle brake components, but verify their dimensions and load suitability.

Use the robot knee range **0° to −120°**, with no hyperextension. Exercise the ankle, subtalar, MTP, hip, and trunk coordinates relevant to each route; a knee-only animation does not validate all 27 actuators.

**7. Standardize physical tackle designs**

Use current McMaster-Carr catalog data to select components based on expected cable tensions, block reactions, dimensions, and commonality.

Use **at most three standardized tackle design families** for the system. Define their topology, dimensions, sheave count, hardware, load rating, and which muscles use each family. Parameterize genuinely common parts explicitly.

Calculate loads from the transmission equilibrium, not from BPA count or maximum BPA force alone. Include off-design/worst-pose loads and a documented design factor. Verify sheave groove compatibility, cable bend radius, pins, bearings, fasteners, brackets, and housing-anchor reactions.

Record part numbers, manufacturer/catalog links, dimensions, ratings, and retrieval dates. Verify bicycle-component suitability from manufacturer data where available. Do not invent ratings, represent breaking strength as working load, or assume ordinary bicycle hardware can carry every calculated load.

If three families cannot safely cover the demand, report the limiting cases and revise routing/sizing before adding a fourth family.

**8. Objective and constraint structure**

Keep torque capacity and usable motion range as the primary performance objectives, using the established human muscle/group mapping. Human targets must come from the stock human model's routes, not the robot's optimized routes.

Use a small, explicitly weighted torque-shape penalty. Define it separately from capacity/deficit terms so it does not remove the incentive to meet the target torque.

Use a **soft penalty** for slight BPA–BPA interference. BPAs are somewhat compliant, so do not impose a universal hard muscle-to-muscle clearance constraint. Calculate overlap from the pose-dependent inflated envelopes and report its depth and extent. Avoid comparisons between unrelated poses of different muscles; evaluate physically simultaneous configurations.

Retain appropriate hard constraints for rigid hardware/bone interference, prohibited foot occupancy, transmission feasibility, valid actuator lengths, required tendon spans, and rated hardware loads.

There must be **no off-axis torque penalty**, including inherited terms disguised as another regularizer. Record secondary joint torques for inspection without penalizing them. Structural mount reactions still matter for sizing and compliance; they are not off-axis joint-torque penalties.

Print the objective terms separately, including weights and normalization.

**9. Verification, execution, and deliverables**

Before the full campaign, demonstrate that perturbing each attachment/routing variable changes the evaluated geometry and relevant objective/constraint outputs. Confirm that bounds have nonzero extent and that the optimizer does not silently restore seed coordinates.

Verify representative monoarticular, biarticular, trunk, foot, and long knee/Bowden cases. Check rigid-limit recovery, both bracket contributions, parallel load sharing, transmission work consistency, diameter endpoints, and wrapping transitions.

Then launch the full, checkpointed 27-muscle campaign on the appropriate machine. Preserve existing results. Use explicit versioned outputs, reproducible seeds, solver budgets, and resume behavior.

Deliver:

- Both rigid and compliant optimized solutions for all 27 groups, including infeasible cases and exit status.
- Seed-versus-final origins, insertions, and guides in their owning frames, with displacement in millimeters.
- Per-pose sleeve length, diameter, tendon length, force, torque, bracket deflection, transmission state, and hardware reactions.
- The three tackle families and their parts/load assignments.
- Packaging envelopes, overlap statistics, hard-constraint margins, and unresolved fit issues.
- A visualization-ready export containing exact transforms, straight sleeve endpoints, tendon paths, tackle geometry, bracket/housing-stop placement, and component identities. Do not export a generic polyline and leave the renderer to guess which portions are BPA.
- Neutral and representative-motion comparisons against the original OpenSim routes, explicitly including vastus intermedius, TFL, erector spinae, and the foot routes.
- Separate torque, BPA length, relative contraction, and moment-arm plots, plus a concise campaign report documenting what changed and what remains unresolved.

Do not report completion based on a summary CSV alone. Show that attachment locations were actually searched, packaging affected optimization, and both mechanical branches were run.

---

*Superseded/corrected by `HARDWARE_AND_ZCODE_HANDOFF.md` (2026-11-03): 10 mm and 20 mm BPAs allowed with new contraction limits (16.5% / 24%), corrected target groups (vas_int, semimem, add_mag3, glut_max3 aggregate, psoas/iliacus topology), pelvis-frame z=0 centerline constraint, Xi0–Xi3 only (Xi4 waived), and hardware candidate families S/M/L.*
