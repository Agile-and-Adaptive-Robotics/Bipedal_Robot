# ACTIVE WORK CONTINUATION — AnimatLab .aproj wiring repair (paused 2026-09-10 ~13:30 PDT)

**Read this fully before touching anything. The previous session (ZCode, Sep 9-10) burned
Ben's trust with repeated partially-broken fixes; the file history is in git, nothing is
lost, but the discipline rules below are mandatory.**

## Where things stand

File: `Neuromechanical_Models/Biped_2xCPG_wSubs/Biped_2xCPG_wSubs.aproj`
Ben's last GUI save is committed (cd59155-era + later "7ccae6b minimizers, mujoco,
animatlab work"). The working tree has ZCode's subsequent wiring-fix patches — some
correct, some incomplete. **Do not build on the working tree; understand it first.**

### What is DONE and verified
- All 15 page CDATAs are individually well-formed XML (validated via PowerShell
  `[xml]` on each extracted page — the whole-file [xml] parse does NOT validate CDATA).
- Every page has exactly one `<Version>1.5.0.1</Version>` header, exact declared
  Nodes/Links counts, `ShowGrid=False`, closing `</AddFlow></Diagram></Root>`.
- No duplicate child-object IDs, no dangling synapse endpoints, 261 synapses,
  20 muscles, grids off everywhere.
- GUI-verified: opened in AnimatLab2.exe on this machine (2026-09-09/10 sessions);
  the version BEFORE the Sep-10 relocate/merge patches opened with ZERO error dialogs
  (Ben confirmed "ok, no errors"). The CURRENT working tree added 12 new synapse blocks
  that are misplaced (below) — expect ~2 pages to error until fixed.

### What the 12 misplaced blocks are
ZCode's wiring_fix.pl added 12 synapse objects (8 cross-joint Ia: Gas/BFlh/Semimem/RF
IaE-OffPage -> Knee MN per side; 4 commissural RG inhibitory links in the NS fragment)
but inserted each block just BEFORE the target subsystem fragment's closing `</Node>` —
i.e. AFTER the fragment's `</Links>` close. The GUI's AddFlow page loader reads only the
fragment's `<Links>...</Links>` list, so these 12 synapses are loaded in the sim
(connexions are flat) but INVISIBLE on the diagram pages, and `relocate2.pl` then moved
8 of them into NEW duplicate `<Links>` lists placed before `<Nodes>` (wrong too).
Current expected state: LH_Knee fixed (4 blocks merged into its real Links list — verify);
NS and RH_Knee still have misplaced/duplicate-list blocks.

## THE FIX (one careful pass)

For each of the three fragments — `LH_Knee Motoneuron`, `RH_Knee Motoneuron`,
`Neural Subsystem`:
1. Locate the fragment's `<Node>` block (its subsystem node, found via
   `<Text>…fragment name…</Text>` then `rindex '<Node>'`), compute its full span with a
   CDATA-opaque balanced scan (skip `<![CDATA[...]]>`; count BOTH `<Node>` and
   `<Node attr...>` opens; `</Node>` exact).
2. Inside that span, find the fragment's OWN `<Links>…</Links>` list (the first one whose
   position is before the fragment's `<Nodes>` open — layout is Links-then-Nodes).
3. Collect the 7b1a-prefixed synapse blocks that belong to this fragment:
   - LH_Knee: 4 synapses with OriginID = the LH knee-page OffPage instances of
     L_Gas/L_BFlh/L_Semimem/L_RF Ia E (b1a00053/b1a00078/b1a00103/b1a00128 full GUIDs);
     DestIDs = LH_Knee MN E (e7f0f45b…) for the first three, LH_Knee MN F (8b3d2c0c…)
     for RF. (RH twins: b1a00051/76/201/226 OffPages; dests 22970258 / 716e9049.)
   - NS: the 4 commissural synapses (origin/dest = the 4 new OffPage node GUIDs, all
     `7b1a0000/0001/0002/0003`-era ids that start `7b1a` — find them by
     `<LinkedNodeID>` = 7635ff71 / b82ffa13 / 095d4bd1 / e04a7f4b).
4. Move those blocks INSIDE the fragment's `<Links>` list (before `</Links>`), and remove
   any duplicate `<Links>` list containing 7b1a blocks that sits between `</Links>` and
   `<Nodes>` (the relocate2 artifact). Also make sure the 4 OffPage standard NODES
   (LinkedNodeID → RG neuron ids) sit in the NS fragment's `<Nodes>` list — if they are
   in the duplicate list, move them into the real Nodes list.
5. Update the fragment's `<AddFlow Nodes="…" Links="…">` declared counts to the true
   element counts of its Nodes/Links lists.

## Verification (ALL of it, every time)
1. Whole file: PowerShell `[xml](Get-Content -Raw)` parses.
2. Each of the 15 page CDATAs extracted and individually `[xml]`-parsed — ALL must pass.
3. No `<Links>` list containing `<ID>7b1a` blocks left outside the real Links position;
   no duplicate `<Links>` lists within any fragment.
4. Connectivity audit (resolve OffPage LinkedNodeIDs + Node IDs to Text labels):
   synapse count and origin→dest pairs match the pre-fix inventory (261+ additions).
5. GUI: `taskkill /IM AnimatLab2.exe /F` first (Ben may have it open — CONFIRM with him),
   then launch `AnimatLab2.exe <path>` and check app state for Error dialogs. Each broken
   page throws its own dialog at load; zero dialogs = pass. Status bar "Load project
   complete". Then close via taskkill WITHOUT saving.
6. Hand to Ben for the final GUI check and commit.

## Wiring reference (what SHOULD be connected — verified against Deng 2019 + Ben's edits)
- Gas MN ← K&A PF-E OffPage (ankle page: 548c27d1 full-id L / f7eed3a5 R) — stance
  plantarflexion synergy (Ben retargeted this himself from PF F; his edit is in his save).
- BFlh/Semimem MN ← Hip PF-E + K&A PF-F OffPages; RF MN ← Hip PF-F + K&A PF-E OffPages.
- Feedback per new muscle: Ia afferent → Ia E IN → reciprocal inhibition of antagonists
  (Gas Ia → AnkleZ-Dorsi MN + Knee MN F per Ben; BFlh/Semimem Ia → Knee MN F; RF Ia →
  Knee MN F), Ia-Ia mutual, Ib autogenic (own MN). Cross-joint Ia→knee MN links are the
  12 blocks this fix relocates.
- Commissural RG: L RG ext –| RH_RG E, L RG flx –| RH_RG F, both directions
  (4 synapses in NS Links list, OffPages in NS Nodes list).
- Ankle terminology renamed: Ext→Planta, Flx→Dorsi (56 labels done).
- Deng 2019 Biomimetics 4(1):21 is the architecture spec (shared knee-ankle PF: PF-F
  drives knee flexor + ankle dorsiflexor in tandem; PF-E drives knee extensor + ankle
  plantarflexor). DF pairs with knee FLEXION (swing); PF with knee EXTENSION (stance).

## After the file is clean
- Ben re-exports the Standalone asim; RG oscillation is then re-tested (previous
  "latch" result was measured on a half-wired connexion set — invalid; re-sweep drive
  2–40 nA and the 10 nA 1 ms pulse kick from scratch).
- Ben positions the 16 new-muscle attachments (currently anatomical placeholders);
  then RestingLength = TSL+OFL, LengthTension window TSL+0.5·OFL..TSL+1.5·OFL,
  Kse ≈ Fmax/(0.033·TSL).
- Open items: no Renshaw on new MNs; no II afferents on new hip muscles; PF→MN
  conductances all 0.5 µS (Deng A2: hip 2.565/3.632, knee 4.93/1.516, ankle 4.054/4.522);
  TFL needs an abduction DOF; LH_HipZ→LH_Hip rename pending.

## Hard rules (do not skip)
- NEVER save from the AnimatLab GUI — taskkill /F instead; a GUI save re-serializes and
  can drop objects.
- NEVER deploy while Ben has AnimatLab open — his stale in-memory copy overwrites fixes
  on save. Confirm AnimatLab is closed first.
- NEVER validate CDATA pages via the whole-file parse only — extract and validate each.
- NEVER regex-patch XML without /e flags and exact-range substr replaces — the
  short-offset lvalue-substr bug and the un-/e replacement bug each corrupted the file
  once already (both fixed in the tools/ scripts).
- If a fix adds objects: re-GUID block id AND every child-object id; strip inherited
  InLinks/OutLinks; escape & as &amp;; patch Value AND Actual attributes together.
