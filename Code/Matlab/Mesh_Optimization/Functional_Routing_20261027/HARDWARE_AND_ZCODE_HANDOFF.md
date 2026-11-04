# Hardware candidates and corrected optimization contract

Catalog research: 2026-11-03 Pacific. Candidate dimensions are suitable for Blender envelopes. These are not released, load-rated complete tackle designs. The robot Blender scene has not been silently resized from this preliminary screen.

## Hardware for Blender

`hardware_candidates.json` is the machine-readable source. `../output/AH001_Hardware_Candidates.blend` contains three candidate sheave families, loose sheaves, two-sheave block envelopes, and Bowden dimensional references. Append its collections into a robot scene. Blue sheaves use verified catalog envelope dimensions; orange side plates are proposed packaging allowances. Axles, groove profiles, and housing terminals are schematic, not manufacturing drawings. Two sheaves in a displayed block do not establish a transmission ratio without a complete reeving diagram.

| Family | McMaster part | OD × width, mm | Shaft bore, mm | Rope diameter, mm | Published capacity |
|---|---|---|---|---|---|
| S | 9466T44 | 31.75 × 6.35 | 4.7625 | 1.5875 | 90 lbf / 400 N |
| M | 9466T41 | 38.10 × 7.9375 | 9.525 | 2.38125 | 150 lbf / 667 N |
| L | 3434T78 | 63.50 × 12.70 | 12.70 | 3.175 | 675 lbf / 3003 N |

Source: [McMaster catalog](https://www.mcmaster.com/products/cable-sheaves/). The table supplies sheave capacities, not an assembled block WLL. Confirm the capacity's load convention, cyclic duty, and cable compatibility before structural selection. Do not treat catalog pulley diameter as the verified wire pitch diameter. Each family can share parts across repeated blocks; do not introduce a fourth family merely to accommodate an unverified route.

Two-sheave axial stack envelopes, excluding axle head/nut and mounting features: S 21.70 mm, M 24.875 mm, L 36.40 mm. These include the two sheave widths, a proposed 1 mm inter-sheave gap, 1 mm clearance at each outer face, and two proposed plates (3/3/4 mm respectively). Planar block allowance is sheave OD + 4 mm. Plate thicknesses and this guard allowance require structural and cable-retention design.

### Bowden references

- [Shimano Y80098411](https://bike.shimano.com/en-NA/products/service-and-upgrade-parts/pdp.P-Y80098411.html): 1.6 mm inner cable, 2050 mm stock length. Its product listing does not establish a robot assembly WLL.
- [Jagwire manufacturer catalog](https://jagwire.com/files/general/2021_Jagwire_AM_Catalog_LowRes.pdf): KEB-SL compressionless brake housing is the bicycle candidate; use its 5 mm OD envelope. BOT192-BJ alloy end caps fit 5 mm housing and have 6 mm OD. Ferrule length is not verified: the Blender library explicitly uses a 12 mm placeholder. A normal cable-tip crimp is not a load-bearing anchor.
- [Sava C277PL, item 100714](https://www.savacable.com/galvanized-steelhdpe-lining-push-pull-casing-1958): 3.556 mm casing ID, 7.0358 mm metal casing OD, **7.747 mm finished OD**, and specified mating cable diameter 2.3622 mm. Do not use the bare casing OD as its finished packaging size. Although the nominal casing ID exceeds 1/8 inch, the specified mating cable does not: this is not a verified 3.175 mm cable housing.
- [Sava control-cable discussion](https://www.savacable.com/blog/push-pull-cable-assemblies-everything-you-need-to-know) discusses bend radii spanning 2–8 inches depending on construction. It does not supply a product-specific minimum for the selected assemblies. The bicycle asset's 75 mm radius is a preview assumption, not a vendor specification. Evaluate routing sensitivity and obtain the actual rated assembly's radius.

Complete Bowden packaging needs the flexible housing, appropriate liner, two ferrules or rated terminals, two housing stops, a cable anchor at each appropriate mechanical end, adjustability, guide/strain relief, and room for articulation. A two-section BPA route needs explicit attachment pins plus these stops; a black tube joining sleeve tips is insufficient.

Use bicycle components as dimensional candidates for individual low-tension branches. Do not place a combined high-force BPA bundle on bicycle housing simply because it fits a 5 mm envelope. A larger rated pull-only control assembly may be needed for high-load knee routes. No such high-load assembly has been selected or rated here.

## Load sizing that the optimizer must supply

Calculate single-BPA force from `Code/Matlab/Functions/maxBPAforce.m` and the validated force–contraction law. At 620 kPa and 450 mm resting length, the repository Fmax equations give approximately 441 N for 10 mm and 1314 N for 20 mm. These are zero-contraction capacities, not forces throughout the sweep.

For the verified reverse-pulley convention, ideal cable tension is `T = n * F / G`; losses modify this consistently with work balance. Do not assume conventional lifting mechanical advantage if the reverse transmission increases travel. For a sheave turning the cable through angle beta, equal cable tension gives radial reaction `R = 2*T*sin(beta/2)`; use the full vector sum when leg tensions differ. At a 180-degree turn, `R = 2*T`.

Example: three 20 mm BPAs at their 450 mm Fmax, combined into a G=2 reverse transmission, would put approximately 1971 N on an ideal shared cable and 3942 N radial load on a reversing sheave. That exceeds even the L candidate's 3003 N catalog capacity if interpreted as radial capacity, before additional reserve. Therefore these three families cannot be assigned indiscriminately to shared bundles. Explore repeated load-separated branches using common families, feasible reeving, or a stronger replacement for one family. Do not approve a 675 lbf sheave for a 3942 N demand by comparing only cable tension against its label.

Size each pin for shear/bending and each plate for bearing, net section, and bracket bending. For a centered double-shear pin: `tau = 2*R/(pi*d^2)`. Ideal equal two-plate bearing stress is `R/(2*t*d)`. Include lever arms and attachment bolt loads. Neither formula establishes fatigue life or the assembly WLL. Check cable and termination strength with their documented deratings. Reevaluate Bowden tension, housing compression, stop reactions, friction, and compliance; housing stops react cable tension, while routing guides can react a vector sum.

Open items: final reeving and per-branch loads, vendor confirmation of sheave capacity convention, cable fatigue/bend limits, complete cable termination ratings, housing compression rating, exact ferrule length and stop geometry, axle retention, and mount geometry. These are explicit unresolved inputs, not permissions to fabricate ratings.

## BPA sizing and target corrections for ZCode

These instructions supersede the previous all-20-mm prompt.

1. Allow nominal 10 mm and 20 mm BPAs; do not use obsolete 40 mm records. For the initial 10 mm screen, take BPA resting length equal to the replaced muscle's longest path length and compare target MIF against approximately half the 10 mm Fmax at that length. Gracilis is an intended small-actuator example. This is a seed-sizing heuristic, not a replacement for torque/strain/packaging optimization.
2. Use maximum contraction assumptions **16.5% for 10 mm** and **24% for 20 mm**: `Lmin = 0.835*L0` and `Lmin = 0.76*L0`. Preserve 20/38 mm diameter endpoints for 20 mm and verify the 10 mm inflated diameter from its own calibration. Never apply the 20 mm diameter law to a 10 mm sleeve.
3. `vas_int` replaces the combined `vas_med_r`, `vas_lat_r`, and `vas_int_r` torque. Expect several parallel 20 mm BPAs. Explore a cam-assisted transmission with terminal attachment at the patellar-tendon attachment on the tibia, or inferior to it within bounded surface patches. Cam profile, pivot, cable contact, bracket and insertion location are optimization variables; do not use the current renderer's excessively anterior point as ground truth.
4. `semimem` targets the summed, signed torques of `semimem_r` and `semiten_r` at each same configuration.
5. `add_mag3` targets the combined adductor group: add_mag1/2/3, add_long, and add_brev in the present stock mapping. Inspect any additional muscles literally named as adductors in the authoritative model and report them; do not silently pull independently assigned groups into this target. Preserve the alias between the legacy add_mag3 name and the current adductor actuator.
6. The historical `glut_max3_r` replacement targeted all gluteus maximus, medius, and minimus compartments. Retain that aggregate target as the baseline. A three-actuator alternative, each targeting the max/med/min compartments with the same index, is only exploratory. Verify torque-direction compatibility and coverage before adopting it. Such a change would replace one group by three, changing 27 groups to 29; report that consequence. It is not approved merely from matching compartment suffixes or summed MIF.
7. Ben specifies that psoas terminates onto iliacus, and both contribute to hip flexion. Represent this as a shared attachment/load-transfer topology, not a fictitious direct bone insertion for psoas or a serial force multiplier. Retain each actuator's force contribution, trace reactions at their junction, and avoid counting iliacus force twice.
8. **No attachment point may cross the pelvis-frame centerline z=0.** Assign side ownership from the seed/model, transform every candidate attachment into the pelvis frame at every evaluated pose, and enforce `sideSign*z_pelvis >= 0` as a hard constraint. This includes moving attachment points and mechanical mount anchors. Do not flip signs or take absolute coordinates after optimization. Resolve any seed on the centerline explicitly; do not let it switch sides.
9. Use Xi0–Xi3, not Xi4. Run genuine rigid and compliant reoptimizations, with both bracket contributions, the applicable dissertation stiffness assignments, movable bounded surface attachments, inner-to-outer packaging, soft BPA–BPA overlap, no off-axis torque penalty, no BPAs/tackle on the foot, and required foot/knee tendon spans as stated in the preceding optimization prompt.
10. Use the hardware JSON and envelopes as dimensional candidates. Optimize discrete family assignment and repeated branches from actual loads; unresolved ratings cannot be assumed infinite. The three candidate families are not certified final designs.

## Coarse force screen

`bpa_force_screen.csv` and `force_screen_provenance.json` sample each existing human target group on a 5-point-per-relevant-coordinate grid, using knee 0 to -120 degrees. The reported maximum is only a sampled lower bound on the true longest path, and the existing actuator map is used unchanged. All 27 entries are included; the exploratory gluteal split is not applied. For groups, the length screen takes the longest sampled constituent path; final robot routes must be sized from their own feasible geometry.

The gracilis screen gives MIF 162 N, sampled Lmax 0.5134 m, and Fmax10 445.3 N; half is 222.6 N, so 10 mm passes Ben's initial screen. The combined vasti give MIF 4530 N, sampled maximum constituent length 0.2742 m, and Fmax20 1228.2 N. Applying the same half-Fmax convention to 20 mm solely as an initial direct-force estimate yields eight BPAs; this is not an optimized count and does not account for a cam or tackle. Torque targets must be summed at each common pose, not created by summing each muscle's separate peak torque.
