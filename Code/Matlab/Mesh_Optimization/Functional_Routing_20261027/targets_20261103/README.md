# targets_20261103 — corrected 27-group target contract (provenance)

`target_contract.json` is the corrected target contract for the Functional Routing
campaign. It implements HARDWARE_AND_ZCODE_HANDOFF.md items 1–9 (2026-11-03),
which supersede the all-20-mm language in MASTER_DIRECTIVE_20261103.md. Built by
the Target-builder workflow on easteregg2, 2026-11-03. Generated artifacts and
verification scripts live in `tmp/functional_routing_20261103/targets/`.

## What was read (evidence for every field)

- `HARDWARE_AND_ZCODE_HANDOFF.md` items 1–10 + the "Coarse force screen" section
  (BPA sizing corrections, contraction limits, centerline rule, Xi4 waiver).
- `MASTER_DIRECTIVE_20261103.md` item 8 context and item 5 provisional layer order.
- `Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json` — group
  membership, primary dofs, MIF grounding, rom_table_deg. Used unchanged.
- `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim` — the stock
  human model; all MIFs and all path-length grids come from it (opensim 4.4.1,
  `D:/Anaconda/envs/myo/python.exe`).
- `Code/Matlab/Functions/maxBPAforce.m` — Fmax law (10/20 mm, 620 kPa).
- `buildKneeFlexorContext20mm.m` (pick 77) and `buildKneeExtContext20mm.m`
  (pick 32) — Xi provenance and bracket/stiffness conventions.
- `Neuromechanical_Models/AH001_Simscape_Repair/packaging/hardware/
  bpa_force_screen.csv` + `force_screen_provenance.json` — LOCATED; used as the
  L0/seed source of record and independently cross-checked (below).

## What was run (all checks executed 2026-11-03 on easteregg2)

1. **Stock MIF census** (`extract_stock_model_20261103.py`): parsed
   gait2392_simbody.osim, 92 muscles. Group sums match the contract's `mifN`:
   semimem 1698 (1288+410), vas_int 4530 (1365+1871+1294), add_mag3 2268
   (381+343+488+627+429), glut_max3 4867 (nine compartments). **Adductor
   census:** exactly five adductor-named right-side muscles exist
   (add_mag1/2/3_r, add_long_r, add_brev_r); NO other literally-named adductor
   → `extraAdductors = []`.
2. **Xi pick verification** (`verify_xi_picks_20261103.m`, MATLAB -batch,
   checkcode clean): `minimizeFlxPin10_results.mat` row 77 (xCols 4/5/6) →
   Xi0 = 0.00393737 m, Xi1 = 39978.7 N/m, Xi2 = 14733.8 N/m;
   `minimizeExt10mmX3_results.mat` row 32 (xCols 5/6/7/8) → Xi0 = −0.00637802 m,
   Xi1 = 39978.7, Xi2 = 14733.8, Xi3 = 0.158289. Lock-check flexor Xi1/Xi2 ==
   extensor Xi1/Xi2: TRUE (one shared pair, advisor requirement).
3. **Screen cross-check** (`crosscheck_screen_20261103.py`, opensim 4.4.1 API,
   5 points per relevant coordinate, knee swept 0..−120 deg, other coordinates
   from the map rom_table_deg clamped to model ranges):
   - all 27 stock group MIFs match the CSV `mif_N` column exactly;
   - all 27 `Fmax10_N`/`Fmax20_N` columns reproduce from the maxBPAforce.m
     formulas at the CSV coarse length (< 0.01 N difference);
   - all 27 ten-mm screen outcomes agree (PASS: sar 156, grac 162, per_tert
     180, ext_hal 162 N — the CSV column matches the ask's
     MIF < 0.5·Fmax10(L0) rule on every row);
   - all 27 CSV `initial_n20_half_Fmax` columns reproduce as
     ceil(MIF/(0.5·Fmax20)) — confirming the CSV's half-Fmax convention;
   - recomputed path lengths match the CSV within 0.05 mm on 20/27 groups;
     deltas on 7 groups up to 20.1 mm (grac 0.4934 vs 0.5134 m) from grid
     convention (this product grid vs the CSV's sampling); no screen outcome
     changes. Rows flagged `mismatch vs CSV conventions: none`.

## Decisions embedded in the contract

- **mifN basis** = stock-model group sums (per the corrected prompt), NOT the
  master actuator MIF. This matters for map ids 15 (tib_post: stock 1588 vs
  master BPA 3549, the documented copy-paste artifact) and 24 (glut_max3: stock
  4867 vs master BPA 5047, the 3.7 %-under all-nine reading). The CSV already
  used the stock basis; the contract follows it.
- **Seed counts** follow the corrected prompt's formula n20 = ceil(MIF/Fmax20(L0))
  (vas_int 4, soleus 3, glut_max3 4, add_mag3 2, semimem 2, ercspn 3, all
  others 1). The CSV's half-Fmax20 alternative (vas_int 8, soleus 6, glut_max3
  8, …) is recorded per group in `sizing_data.n20_csv_half_Fmax20_alternative`
  and flagged as NOT adopted. Both are initial estimates only; cam and tackle
  effects are excluded, and counts are re-solved in the campaign.
- **L0 seeds** = CSV coarse_Lmax (longest sampled constituent path), a sampled
  lower bound; the contract says so in every sizing string and requires
  re-derivation from the robot's feasible geometry.
- **10 mm seeds** (sar, grac, per_tert, ext_hal): 1× 10 mm each, screen pass at
  handoff item 1's rule. The 10 mm inflated-diameter calibration was NOT
  located in this ask (unresolved[3]).
- **layerOrder** = master item 5 provisional order, with one adjustment
  grounded in handoff item 7 (psoas 27 + iliacus 26 solved as a coupled
  junction pair first). The geometry extraction's layering table was NOT
  available on easteregg2 (searched the campaign folder,
  `tmp/functional_routing_20261103/`, and the packaging folder), so no
  bone-distance adjustments could be applied — unresolved[0]. Rule carried:
  deepest layer optimizes against bone; every later group treats all earlier
  layers' occupied envelopes as hard obstacles; coordinated refinement after
  the first pass.
- **Gluteal split** is exploratory only; adopting it would change 27 → 29
  groups and is NOT approved by this contract (`glutealConsequence`).
- **Psoas/iliacus** carry the handoff item 7 shared-attachment/load-transfer
  topology in both groups' targetDef: no fictitious direct bone insertion for
  psoas, no serial force multiplier, iliacus force not counted twice.
- **Xi** = Xi0–Xi3 only (Xi4 waived, handoff item 9), sources as recorded in
  `xiSource` (builders pick rows 77/32; mats of record re-verified in MATLAB).

## Unresolved (mirrored in the contract's `unresolved[]`)

0. Geometry-extraction layering table absent — provisional layer order carried
   unadjusted.
1. n20 count-convention conflict (prompt formula vs CSV half-Fmax convention) —
   documented, both values recorded per group.
2. L0 seeds are sampled lower bounds (grid deltas up to 20.1 mm; no outcome
   changes) — final resting lengths come from feasible robot geometry.
3. 10 mm sleeve inflated-diameter calibration not located.
4. Hardware families S/M/L and Bowden references are dimensional candidates
   only; ratings unresolved per the handoff's "Open items" — no load rating is
   assumed by this contract.
