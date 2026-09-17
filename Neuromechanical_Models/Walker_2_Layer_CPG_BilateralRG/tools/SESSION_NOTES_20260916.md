# Walker_2_Layer_CPG_BilateralRG — session notes (2026-09-16, easteregg2)

Working copy of `..\Walker_2_Layer_CPG` (original untouched). Everything below
is in BOTH the .aproj (GUI-visible, XML-valid, GUI-verified zero error dialogs)
and the `_Standalone_modern.asim` (test vehicle, runs headless via
AnimatSimulator.exe) unless noted. Deterministic GUID prefix: `cafe____-0000-4000-8000-…`.

## Build chain (run in order from backup_v0)

1. `build_rg.pl` — R RG half-center (R RG ext/ext IN/flx/flx IN, ipsilateral to
   R Hip/Knee PF), removed the 4 crossed L-RG→R-PF connexions, Stimulus_2
   (10 nA, 0-10 ms, into R RG flx → legs start neurally antiphase), asymmetric
   start pose (femur_L Z=+12°, femur_R Z=−12°, tibia Y=−28° both), 4 chart
   columns. Air: 11/12 bursts, 0.443 s, clean antiphase.
2. `build_comm.pl` — Shinohara et al. 2025 (bioRxiv 2025.11.11.687930)
   commissural paths: S RG flx→`S c1`→(inh) contra RG flx;
   S RG ext→`S V3`→(exc, weak) contra RG ext IN.
3. `patch_comm_types.pl` — dedicated types: "c1 Commissural Inhibit"
   (Equil −70, SynAmp 2.749) and "V3 Commissural Excite" (Equil −40, SynAmp 0.1).
   **REQUIRED: with RG-Excite (SynAmp 2.749) the V3 path latches both RGs.**
   Air after fix: 11/12 bursts, 0.451 s, c1/V3 burst phasic with their HCs.
4. `build_aff.pl` + inline 0.01 tune — flexor Ia + NEW II (stretch receptor →
   adapter `DataTypeID=II` → `S J flx II` relay) → exc flexor HCs (PF flx + RG
   flx); extensor Ib → exc extensor HCs. Type "Afferent HC Excite" (Equil −40,
   **SynAmp 0.01** — 0.1 latches in air; 0.05 slows to ~0.8 s; 0.01 keeps
   0.477 s antiphase). 24 links.
5. `build_contact.pl` — `L/R heel contact` (foot_*_contact body) + `L/R toe
   contact` (toe_*_contact body) neurons via PhysicalToNode adapters
   (`SourceDataType=ContactCount`, Target ExternalCurrent, Gain C=20), each
   exciting ipsilateral extensor layers: RG ext, Hip/Knee PF ext, Hip/Knee MN
   ext (20 links). Air: inert (no contact), rhythm unchanged. Ground: works.
6. `ground_test.pl` — ground variant generator: Root Freeze→False, Root y
   0.993→−0.15, pelvis shelf box (top y=−0.28) = partial support, Contact.txt
   diagnostic chart. Deliverable: `Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim`.

Backups: `tools/backup_v0` (pristine copy) → `v1_bilateralRG` → `v2_commissural`
→ `v3_preafferent` → `v4_precontact`.

## Hard-won mechanics (do not relearn)

- **AddFlow drawing format (decoded 2026-09-16 evening, verify_handles.pl):**
  a page drawing `<Link Org="N" Dst="M">`'s Org/Dst = **0-based index of the
  endpoint NODE entry in the CDATA's interleaved (nodes+links) file order**.
  Verified: all 208 original synapse drawings match exactly. Consequences:
  (a) appending shapes at the END never shifts existing indexes (safe);
  (b) deleting a mid-file drawing shifts everything after it and silently
  re-docks all later arrows (the original "MyLink cast" disaster and a second
  20-arrow drift both came from this — fix_drawings.pl now normalizes every
  drawing's Org/Dst after any structural change);
  (c) cloning a template drawing without recomputing Org/Dst renders the new
  arrow ON TOP of the template's arrow (Ben's "no links visible" bug — the
  new functional links showed in the tree, but their drawings duplicated old
  arrows). fix_drawings.pl rebuilds all appended drawings with computed
  endpoints: final state 296/296 endpoint-exact (208 synapses + 88 adapters).
- Physical bodies (muscles, SRs, foot/toe contact boxes) ARE drawn as page
  nodes — adapter links dock to them like any node.
- **Effective synaptic strength in the ASIM = SynapseType SynAmp.** The
  per-connexion `<G>` is IGNORED by AnimatSimulator (G=1e-4 vs 0.15 gave
  bit-identical runs). In the APROJ the per-Link `<SynapticConductance>` and
  the type's `<MaxSynapticConductance>` mirror it — set all three consistently.
  Tuning = edit the TYPE's SynAmp/MaxSynapticConductance.
- **Never delete drawing `<Link>` entries from the page CDATA** — other links'
  Org/Dst handles break → GUI error "Unable to cast MyLink to Lassalle.Flow.Node".
  Orphan drawing links (drawing kept, functional link deleted) are tolerated.
  This is why the 4 crossed-link arrows still show stale on the page.
- APROJ uses `<SynapticTypeID>`/`<Text>`/attribute-triplets (`Value/Scale/Actual`);
  ASIM uses `<SynapseTypeID>`/`<Name>`/plain values; rotations are DEGREES in
  aproj, RADIANS in asim.
- Cloned blocks: re-roll EVERY `<ID>` (Ca channels, Gain children) or the sim
  fails "same ID twice". `reroll_ids()` in each builder.
- In `s///` replacements, `"".expr.""` is LITERAL text (no /e) — this wrote
  `Top="".(560+…).""` into the CDATA and the GUI died with "Name cannot begin
  with '.'". Interpolate via a real variable.
- GroundPlane/WalkingPath are TOP-LEVEL bodies (siblings of Root, depth 1) —
  unfreezing/lowering Root does not move the ground. Root was Freeze=True at
  y=0.993; toes ≈1.17 m above the path (path top y=−0.8936).
- GUI verification loop: `tools/gui_err_text.ps1 <aproj>` (launches AnimatLab2,
  dumps any Error window text, kills). WinForms error windows are NOT #32770
  dialogs — enumerate by title.
- Analysis helpers (in %TEMP%): rganalyze.pl (burst/period/antiphase),
  vmean.pl, connmap.pl (full connexion map). Repo copies could be made later.

## Final page layout (after fix_drawings.pl, 2026-09-16 evening)

- R RG half-center mirrors the L RG block exactly, shifted down 150 px
  (L RG at x≈13-94 / y≈578-654; R RG at the same x, y≈728-804).
- c1/V3: 2×2 block in the empty right band (x≈1580-1760, y≈200-320).
- II chains: 4 rows in the right band (y≈400-680): adapter at x≈1580,
  II relay at x≈1670 (L Hip / R Hip / L Knee / R Knee top to bottom).
- Contact groups: bottom row (y≈1354), each a body→adapter→neuron triplet:
  L heel (80/170/260), L toe (340/430/520), R heel (600/690/780),
  R toe (860/950/1040).
- The 4 stale crossed-arrow drawings are DELETED (safe once every drawing's
  Org/Dst is recomputed afterwards, which the fixer does).

## Behavior summary (5.1 s runs)

| config | period | L/R RG ext bursts | notes |
|---|---|---|---|
| single RG (baseline) | 0.443 s | 12/12 | legs antiphase via crossed wiring |
| bilateral RG + poses | 0.443 s | 12/12 | identical L timing; R mirrors |
| + commissural (fixed) | 0.451 s | 11/12 | c1/V3 phasic |
| + afferents @0.01 | 0.477 s | 10/11 | in air |
| + contact layer (air) | 0.466 s | 11/11 | contact neurons silent |
| ground + shelf (no fb) | 0.451 s | 11/12 | toe duty L 22% R 40% |
| ground + all feedback | 0.466 s | 11/11 | toe duty L 16% R 38%, feet 4%/2% |

Ground observations: rhythm survives contact in every config; the walker leans
on the shelf (R-loaded asymmetry); stepping is irregular without balance —
expected, no posture/balance layer exists. Contact gain C=20 and the
"Afferent HC Excite" SynAmp are the tuning knobs for stance reinforcement.

## REMAINING (Ben's task list item 5) — RG/PF/MN subnetwork pages

Structure discovered: `<NervousSystem>` holds ONE `<Node>` = NeuralModule
("Neural Subsystem") containing SynapseTypes + flat neuron `<Node>`s + Links +
one `<DiagramXml>` CDATA page. The Biped (`Biped_2xCPG_wSubs.aproj`) has 15
pages: each page belongs to a child subsystem Node carrying its OWN `<Links>`
+ `<DiagramXml>` (page names LH_RG, LH_HipZ, LH_K&A Pattern Formation, …MH
Motoneuron, etc.). Recipe: create 3 subsystem child Nodes ("RG Layer", "PF
Layer", "MN Layer"), partition the flat neuron Nodes + Links into them, give
each a page CDATA built from the existing drawing entries (keep node
Left/Top; rebuild link entries fresh with template Org/Dst — stale endpoints
are tolerated). GUI-check after. Neurons are NOT nested in subsystems in the
Biped (they stay in one flat `<Nodes>`); subsystems own links + pages.

Also open: II + contact-neuron columns are not in the .aform charts (they are
not charted; add clones of an RG column if wanted); the 4 stale crossed-link
drawings; the aform `Rhythm_Generator.aform` has c1/V3 columns only.
