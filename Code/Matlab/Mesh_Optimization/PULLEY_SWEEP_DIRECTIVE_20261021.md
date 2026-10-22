# Pulley transmission sweep directive (Ben Bolen, 2026-10-21)

Standing directive for the next optimizer campaign. Recorded verbatim-intent from
Ben; execute under the same rules as the 2026-10-21 doctrine round (records of
record untouched, sizes re-solved per configuration, regression gates).

## 1. The configuration sweep

Sweep the optimizer over the full transmission grid:

- **nBPA ∈ {1, 2, 3}** — BPAs in parallel (or mirrored) feeding the tackle,
  from three down to one.
- **tackleLineParts G ∈ {1, 2, 3}** — the tackle's moving line parts ("up to 3
  tackles" in Ben's words; the correct term is tackle line parts / rope parts,
  the class field `TackleLineParts`, a.k.a. purchase; travel gain = G,
  force gain = nBPA/G).

This covers, in Ben's words: "3 BPAs in parallel down to one, down to one using
a pulley, and then multiple BPAs in parallel using up to 3 tackles" — i.e., the
direct-drive ladder (3/2/1 BPA at G=1), the single-BPA-plus-tackle column, and
the parallel-pairs-and-triples with tackle columns. The 2026-10-21 doctrine
round already covered (2,G2), (2,G4), (4,G4); this directive completes the
grid at G up to 3 with nBPA up to 3.

Reminder (proven 2026-10-21): actuator sizes MUST be re-solved per
configuration; a direct-drive-sized actuator runs slack at G ≥ 2.

## 2. Moving wrap point (routingMode 'moving_via')

When the wrapping point moves, the bracket deflections must be calculated from
the BPA muscle origin, and the increased internal force accounted for. The
tackle changes the load split (per-BPA force F on the BPA side, nBPA·F/G in the
tendon), and the moving-via geometry rotates/translates the deflection load
path with joint angle, so the compliance solve must reference the BPA origin
and carry the tackle-amplified internal tension through the equilibrium.

## 3. Bowden branch (routingMode 'bowden')

If a bowden cable arrangement is used:
- Give the housing reasonable packaging space (the class carries envelope
  fields `BowdenBossDia` = 7 mm boss and `BowdenRunClearance` = 5 mm as
  constructor defaults — revisit with real hardware dimensions).
- Calculate the deflection of the bowden bracket under the transmitted loads.
- Buttress the bracket stiffness with realistic materials, assuming aluminum
  or CroMoly steel as in bicycles (bicycle brake/shifter cable standards;
  e.g., 6061-T6 aluminum or 4130 chromoly).

## 5. Configuration-selection heuristic (Ben, 2026-10-21, later message)

The solver for multiarticular joints and for single or multiple muscles replaced by
single or multiple BPAs should operate on a fuzzy escalation logic:

- **Lacking in torque?** → more BPAs in parallel.
- **Lacking in range of motion?** → more tackle line parts.

That is, the configuration ladder is chosen by which resource the site is short of,
per site and per pose region, not swept blind.

- **Hip adduction/abduction coupling:** Ben believes limit code already exists that
  limits hip adduction/abduction based on hip flexion/extension. NOT found in the
  bi-pulley builder — `biPulleySpecsFromOpenSim.m` pins secondary coordinates
  (hip adduction/rotation, subtalar, mtp) at zero, and `Add_Mag_Mesh_Opt.m` only
  sweeps pure adduction/abduction. Locate Ben's coupling code (or write it) before
  the hip sites' requirement envelopes are trusted.
- **Extreme-pose corners:** poses like full flexion combined with abduction may be
  achievable by humans but the solver must NOT be over-penalized for failing to
  reach them — soft-penalize or cap the requirement envelope at extreme combined
  poses rather than enforcing full coverage there.

## 6. OpenSim muscle-path artifact audit (before trusting requirement envelopes)

Known issues with the Gait2392-lineage models, to audit per muscle before trusting
required-torque curves at range extremes (Ben asked to check the OpenSim skill and
docs; the skill has nothing — these are from the OpenSim docs and forums):

1. Gait2392 represents muscle wrapping over bone / retinacula with CONDITIONAL
   VIA POINTS (especially at the knee) rather than true wrapping surfaces —
   paths can jump and moment arms can step at the activation angles.
   (OpenSim docs: Gait 2392 and 2354 Models page.)
2. Wrapping-solver artifacts: muscles wrapped on cylinders can flip sides during
   motion (moment-arm discontinuities), and customized paths can penetrate bone
   (simtk forum threads 10885, 20467).
3. GeometryPath results changed across OpenSim versions (post-4.1 changes;
   opensim-core issue #4330 possible body-force/moment-arm bug) — recompute
   requirements on the pinned opensim 4.4.1 and record the version.

Audit step for the sweep campaign: for each site, visually verify the OpenSim
muscle path across the full joint range used for the requirement (no bone
penetration, no via-point jumps at the graded angles), and flag any site whose
requirement curve inherits an artifact — the soft extreme-pose penalty of §5
covers the corners this audit cannot rescue.

## 4. Deliverables

Per the doctrine round: torque surfaces per muscle/muscle-group, margin
surfaces over the configuration grid, quotable RoM/peak tables per site,
driver copies + dated result mats/CSVs in Results/, SURFACES_INDEX.md updated.

Provenance: Ben's directive given 2026-10-21 after the doctrine round
(workflow dwfrun-60eba61a + pulley_doctrine_post_20261021); recorded by the
agent session so future optimizer campaigns execute it without re-derivation.
