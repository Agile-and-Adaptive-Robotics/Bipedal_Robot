# ACTIVE WORK — AnimatLab .aproj neural surgery (started 2026-09-09, CONTINUES at 8am PT session)

Ben's asks, in progress: (1) complete the RH-side feedback/connections, (2) add
gastrocnemius, biceps femoris long head, rectus femoris, semimembranosus to each leg and
drive them from the existing pattern-formation layers, (3) connect LH and RH with a rhythm
generator layer. Ben also wants TFL eventually (needs a frontal-plane abduction DOF the Li
legs don't have yet) and the LH_HipZ/RH_HipZ labels renamed LH_Hip/RH_Hip (cosmetic, not
done). **Ben added some missing connections manually in the GUI on 2026-09-09 evening and
reports more are still missing — auditing drawn-arrow completeness against the standard
section is the FIRST task of the next session.**

## Current file state

- `Neuromechanical_Models/Biped_2xCPG_wSubs/Biped_2xCPG_wSubs.aproj` — GUI-verified clean
  open on 2026-09-09 (zero dialogs; earlier builds threw 1–15 dialogs per open).
- 20 muscles (12 original + 8 new biarticular: gas/bflh/semimem/rf per leg, Gait2392 Fmax
  applied: 2241/896/1288/1169 N), 122 neurons, 261 synapses, grids off everywhere.
- Completed: RH half-center + coordination wiring (mirrored from LH), RH knee/ankle
  muscle-drive adapters repointed off Renshaw cells onto MNs, LH↔RH commissural RG
  inhibition (4 synapses via OffPages on the top page), 2 nA tonic on L RG ext,
  full Deng-style chains for the 8 new muscles (MN, afferent Ia, Ia E, Ib, Muscle+SR
  nodes, 3 adapters each), Gait2392-derived Fmax. Gas Ia→Ia correctly retargeted to the
  ankle Ia F interneuron.

## Known remaining gaps (next session, in order)

1. **Missing drawn arrows** — Ben found more after the fix pass. Likely cause: links
   whose two endpoints never shared a page were left functional-but-undrawn. Audit each
   page's drawn arrows against the standard section; draw the missing ones via new
   OffPage instances on a page containing both endpoints.
2. **RG does NOT oscillate** — 2–10 nA tonic on L RG ext latches (L ext depolarized,
   everything else suppressed, even over 20 s). Sweep drive 20–40 nA; try the paper's
   10 nA 1 ms pulse; apply Deng Table A2 PF→MN conductances (hip 2.565/3.632, knee
   4.93/1.516, ankle 4.054/4.522 µS) instead of the 0.5 µS defaults.
3. **16 placeholder attachments** — new-muscle origin/insertion points were copied from
   anatomically similar sites; Ben must position them. Then set RestingLength = TSL + OFL
   (Gait2392: Gas 0.45 m, BFlh 0.435 m, Semimem 0.439 m, RF 0.424 m), LengthTension
   window = TSL+0.5·OFL to TSL+1.5·OFL, Kse ≈ Fmax/(0.033·TSL).
4. No Renshaw cells on new MNs; no II afferents on the new hip-spanning muscles.
5. TFL addition — needs a frontal-plane abduction DOF added to the Li legs first.
6. LH_HipZ/RH_HipZ → LH_Hip/RH_Hip label rename (cosmetic, pending Ben's OK).

## Toolchain (in `Neuromechanical_Models/Biped_2xCPG_wSubs/tools/`)

Full deterministic pipeline, in this order from a pristine baseline (git
`6437441:Neuromechanical_Models/Biped_2xCPG_wSubs/AnimatLab_backups/Biped_2xCPG_wSubs.aproj.bak_20260909_preRH`
is the pre-surgery original): `build_A2.pl` (RH wiring) → `build_B2.pl` (muscles) →
`repair_dup_ids2.pl` (re-GUID child IDs) → `fix_gas_ia3.pl` (Gas Ia retarget) →
`rebuild_pages_v3.pl` (regenerate all 15 page drawings from the standard section; drops
stale, dedups, recomputes all Org/Dst, grids off, exact counts) → `move_links_by_page.pl`
(place generated links into the fragment whose page draws them) → `set_fmax2.pl`
(Gait2392 Fmax). Verify with `verify_placement.pl` (link placement), `verify_final.pl`
(battery), `shape_diff.pl` (element-shape diff vs a reference file).
**If any input changes, re-run the WHOLE chain — never patch a patched file.**

## CRITICAL format gotchas (each cost hours — do not relearn)

- The .aproj has TWO layers: the standard section (functional model: neurons, links,
  muscles — full GUIDs) and per-page `<DiagramXml><![CDATA[...]]></DiagramXml>` AddFlow
  drawings (GUI only, objects referenced by `<Tag>`). The whole-file XML parse does NOT
  validate CDATA contents — validate every page CDATA individually as its own XML doc.
- Subsystem fragments are SCATTERED top-level blocks after `</NeuralModules>` (not
  nested); hierarchy is via SubSystemID references. A drawn link must live in the same
  fragment as the page that draws it, else the GUI shows nothing and clicking throws
  "No item with ID". A page's OffPage instances are that page's local references.
- Depth-counting `<Node>`/`</Node>` MUST (a) skip `<![CDATA[...]]>` spans (they contain
  drawing `</Node>` closers with no matching opens in that counting scheme) and (b) count
  both `<Node>` and `<Node attr...>` opens. Getting this wrong silently redirects
  insertions into page CDATAs → "GetAttribute(Int32 i) out of range" pop-ups (one per
  broken page).
- Page rebuilds must emit, in order: prologue (through `<AddFlow ...>`), **`<Version>`
  header** (inside the body, before the first node — dropping it throws the same
  GetAttribute error), nodes, links, tail (fillers + `</AddFlow></Diagram></Root>`).
  Missing tail = "An error occurred while deserializing the xml data."
- Cloning any block: re-GUID the block's own ID AND every child-object ID inside
  (StimulusTension, LengthTension, Gain, CaActivation, CaDeactivation, PID) or the C++
  sim throws "Attempted to add an object with the same ID twice". Strip inherited
  InLinks/OutLinks. Escape `&` as `&amp;` in Text. Patch `Value` AND `Actual`
  attributes together. Regex replacements that build XML MUST use the /e flag.
- The recurring silent-corruption warning: lvalue `substr($x,$p,length($seg)) = $seg`
  with a GROWN $seg eats trailing bytes (always use the original range length).
- **Verification protocol (MANDATORY before handing the file to Ben):** launch
  `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe` with
  the .aproj path as argument, wait ~15 s, inspect the app state for an "Error" dialog
  window. Each broken page throws its own dialog at load, so zero dialogs ⇒ all 15 pages
  deserialized. Status bar must show "Load project complete". Then `taskkill /IM
  AnimatLab2.exe /F` — NEVER save from the GUI (it overwrites the file with stale memory).
  **Never deploy while Ben may have AnimatLab open** — his open session can't see file
  changes and a save from it clobbers the deployed fixes.
- Run headless sims via `bin\AnimatSimulator.exe <path>\file.asim` (.asim only, never
  .aproj). Chart files (e.g. "Rhythm Generator.txt") land in the sim file's folder; chart
  `<EndTime>` caps collection independently of SimEndTime.

## Background

The .aproj is Deng 2019 (Biomimetics 4(1):21) two-layer CPG (RG→PF→MN, Ia/Ib/II
afferents, Renshaw) ported onto the Li biped biomech; Ben+Connor's walker paper documents
the lineage. The old `_Standalone.asim` exports are stale; Ben exported a fresh
`Biped_2xCPG_wSubs_Standalone.asim` on 2026-09-09 (runs clean headless; contains a
Rhythm Generator chart with RH_RG columns + L Hip charts). Tuning ladder: kinematics →
virtual-walker ground walking → stable walking → transitions (details in AGENTS.md and
the animatlab skill's walker-cpg-architecture.md).

