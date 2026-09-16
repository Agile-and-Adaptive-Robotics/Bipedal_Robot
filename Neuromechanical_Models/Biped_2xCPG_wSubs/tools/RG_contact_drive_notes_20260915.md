# RG oscillation — root-cause findings (2026-09-15 session, EB475WS4)

Why the RG half-center never oscillates, and what the working reference model
actually does. All evidence gathered headless (AnimatSimulator.exe + chart .txt),
scripts live in `C:\Users\Ben Bolen\AppData\Local\Temp\` (build_contact.pl is the
reusable harness builder).

## The working reference: Walker_2_Layer_CPG (CoMorrow branch)

`origin/AddingStepSensor_CoMorrow_stw:Neuromechanical Models/Walker_2_Layer_CPG/Walker_2_Layer_CPG.aproj`
(same folder also has K_A synergy + Walker_Connor .aproj; the three standalone .asim
files on that branch are byte-identical PHYSICAL-ONLY exports — no neurons, useless).

W2L's RG subcircuit (L side, R mirrored) — NO tonic stimulus anywhere:

- L RG ext / L RG flx: NonSpiking, InitialThreshold **+50 mV** (NOT −55), GMaxCa 1.5,
  CaAct mid −40 slope 0.2, CaDeact mid −60 slope −0.6 τ350 — **identical to ours**.
- ext→flx and flx→ext direct: "RG to RG Excite" (Equil −40, type G 0.1), connexion G **0.5**.
- ext→extIN: "RG Excite" (Equil −40, 2.749), G 0.5; extIN→flx: "RG Inhibit" (Equil −70, 2.749), G 0.5 (same on the flx side).
- **Drive = ground contact, not tonic**: spiking neurons (InitialThreshold **−55** — that's
  where the −55 belongs) "L/R_foot ground contact", "L/R_toe ground contact", fed by
  PhysicalToNodeAdapters from contact RigidBodies, with Sigmoid gains.
  - L_foot → L RG ext, SpikingChemical "Depolarizing IPSP" (Equil −50), connexion G **6** (strongest in circuit)
  - L_foot → R RG flx, same type G 6 (**contralateral**: loaded foot drives the OTHER leg's swing flexor)
  - L_foot → L_foot self G 4; foot → own knee-flx/plantarflex "inhibit" interneurons G 6–8
  - "L_CPG inhabit" → L RG ext IN, "R_CPG inhabit" → L RG flx IN, G 1 (interleg modulation)
- Connexion gains: ours are G=5 where W2L uses G=0.5 (10× stronger).

## What was proven on OUR model (Biped_2xCPG_wSubs_Standalone.asim harness)

1. **Tonic-driven oscillation is impossible with these neurons.** Matrix over tonic
   {2..60 nA} × RG G {5, 0.5} × InitialThresh {50, −55}: every case LATCHES (winner-take-all),
   never alternates. The RG's plateau (GMaxCa 1.5) holds the ON cell up; with constant tonic
   there is no escape. W2L doesn't oscillate from tonic either — the drive is REMOVED when
   the foot unloads; the sensory pattern is the pacemaker.
2. **The −55 mV InitialThreshold "fix" on our 8 RG neurons was WRONG** (class mix-up with
   W2L's spiking contact neurons). Reverted 2026-09-15; .aproj RG back to +50 mV.
3. **The biped in the standalone export HANGS**: Root body has `<Freeze>True</Freeze>`
   (mass 30991 g ≈ 31 kg box, Position y=1.02). Pelvis welded mid-air → feet never reach
   the ground → ContactCount = 0 on ALL 41 bodies for the whole sim (verified by charting
   ContactCount directly). This is why "RG doesn't oscillate" runs showed a frozen pose.
4. **Contact sensing works once the pelvis is freed** (Root Freeze→False in a harness copy):
   toes hit the ground, and plain `<Adapter><Type>PhysicalToNode</Type>` elements with
   `SourceDataType=ContactCount`, SourceID = toe body GUID, TargetID = neuron,
   TargetDataType=ExternalCurrent, Gain polynomial C = drive scale DO inject current
   (RG visibly driven at collapse). The dedicated `<Type>Contact</Type>` ContactAdapter
   class exists in the DLLs (binds a ContactSensor DEVICE + RigidBody + FieldPairs) but was
   not needed. ContactSensor DEVICE clones (with Bell FieldGain A=0 — degenerate) are not
   required for the generic-adapter route.
5. First harness runs with freed pelvis: biped collapses (no posture tone) → one big contact
   pulse → static. With ext tonic 20 nA + contact drive 20 nA/contact: RG ext engages
   (duty 95%, 93% antiphase vs the transient) but no sustained stepping yet. Next lever =
   enough PF/MN tone to hold a standing pose, then let contact asymmetry drive swings.

## Reusable harness

`%TEMP%\build_contact.pl <gainC> <ve> <vf>`: loads the standalone .asim, frees the Root,
ensures ContactSensor devices on toe_L/R_contact, adds 4 contact adapters
(toe→own RG ext, toe→contralateral RG flx) at gainC A/count, sets RG tonics (ve ext, vf flx),
runs headless, analyzes Rhythm Generator.txt (bursts/duty/period/antiphase).

AnimatLab gotchas (new): chart DataColumns need VALID GUIDs (invalid → whole chart silently
writes 0 bytes); SimEndTime == chart EndTime → end-of-sim flush never fires (keep sim end
longer); an adapter's Gain child object needs its own fresh GUID.

## Tuning ladder for a W2L-style walker (next sessions)

1. Standing: enough RG/PF drive that the freed biped holds a pose on its feet
   (tone on extensors; possibly PF layer directly). Start from the standing keyframe.
2. Verify symmetric contact → co-contraction (both feet loaded).
3. Perturb (or let physics tip) → unloaded foot drives contralateral swing.
4. Reduce our RG connexion Gs 5 → 0.5 (W2L values) once the pattern exists, then re-tune.
5. Port the winning config to the .aproj via the tools pipeline (standard + page CDATA),
   GUI-verify clean open, then Ben commits.

## UPDATE 2026-09-16: W2L air test PASSES; the starter-kick recipe

Ben exported `Walker_2_Layer_CPG\Walker_2_Layer_CPG_Standalone_modern.asim` from the GUI.
Headless run (AnimatSimulator.exe): VERIFIED OSCILLATION — L Hip MN ext/flx burst in clean
antiphase ~1 Hz (ext −48..−55 mV while flx −83..−87, then swap), PF follows, hip/knee swing
freely (suspended). The drive that starts it: `Stimulus_1` = Tonic Current **10 nA for
0.01 s at t=0** into one RG neuron — a kickoff pulse, nothing sustained. The oscillator
then runs endogenously at W2L's connexion gains (G 0.5). Foot contact (the .aproj wiring)
phase-locks it for ground walking but is NOT needed for air-stepping.

Revised understanding of the earlier "tonic always latches" matrix: latching happens under
SUSTAINED drive. The W2L recipe = brief kickoff pulse + intrinsic plateau/accommodation
dynamics + W2L gain balance. Biped air-test port list (all .aproj-visible, GUI-verified):
(1) per-side RG (Ben: ONE RG per side driving that side's PF layers; weak L↔R coupling
LATER — get each side working first; the export's L side still has old-pattern drive),
(2) RG connexion Gs 5 → 0.5 per W2L table, (3) replace
sustained 2 nA tonic with a 10 ms kickoff stimulus (a Stimuli object — visible in GUI),
(4) optionally asymmetric flexed/extended start. Also this session: repaired
Biped_2xCPG_wSubs.aproj — a prior patch script's stdout (473 bytes) had been written into
the file ahead of <Project>; stripped, XML verified. Ben's rules: never touch
Walker_2_Layer_CPG files; every neural change must land in HIS .aproj (GUI-visible), no
asim-only "backdoor" edits; no physics-engine/binary changes needed (modern asim runs
stock Vortex+Double).
