# W2L tuned corrected-connectome winner -> AnimatLab drop test (2026-10-20, egg)

Ben's ask (2026-10-20): "Did you get the harness to drop in AnimatLab or remove
it and get it to walk after it's tuned?" + his GUI observation that in AnimatLab
"the walker obviously fell down and was trying to walk in circles on its side."

Stage: port the TUNED corrected-connectome winner
(`reports_20261020\w2l_contact2\round_best_r030.json`; study
curr_w2l_contact2_r030, obj20 91.36, m30 fresh survive, 31/28 cycles, fwd
+1.20 m in MuJoCo) back into the AnimatLab Walker_2_Layer_CPG_BilateralRG model
and run the drop protocol headless, 30 s. NO commits; Ben's originals untouched
(`git status` on the model folder shows only the 4 new untracked variant files
below).

## 1. VERDICT (honest)

**No drop-then-walk in AnimatLab, under any condition tested.** Five headless
30.1 s runs (Vortex, same engine as the GUI project; dt 1 ms):

| run | connectome | drive | falls? | RG rhythm (flexor -55 mV bursts / 30 s) |
|---|---|---|---|---|
| **stock30** (old asim wiring, only window extended) | old | 2 kickoffs | **YES** - pelvis 0.99 -> 0.10 m by t=1.5 s | **63/63, period 0.479/0.480 s, 43 bursts in the last 10 s - keeps STEPPING while lying down** (RootZ wanders +-1.5 m, X drifts -0.15 m/s) |
| **control_no2na** (corrected, Ben literal, no intrinsic bias) | corrected | kickoffs | YES - 0.107 m by t=1.21 s | **0 - RG flat** (std 0.05 mV over 20 s) |
| **control** (corrected, Ben literal incl. his 2 nA L RG ext bias) | corrected | kickoffs + 2 nA bias | YES - 0.107 m by t=1.21 s | 0 - RG flat (2 nA bias RULED OUT as the rhythm killer) |
| **tuned_nodrive** (winner feedback ratios, no tonic) | corrected + tuned ratios | kickoffs | YES - 0.110 m by t=1.24 s | 0 - RG flat (body still crawls +1.85 m forward at 0.05 m/s while down) |
| **tuned** (THE WINNER, full set) | corrected + tuned ratios | kickoffs + te 5.26781 / tf 7.5507 nA tonic | YES - 0.074 m by t=1.29 s | **0 - LATCHED**: frozen fixed point from t~0.25 s (L RG ext -55.99 mV, L RG flx -58.18 mV, constant to 30 s); hip flexor muscle tension ramps 210 -> 350 N, MN ext -94 mV inhibited = legs latched in flexion |

Read:
1. **Ben's GUI failure mode is reproduced and explained**: the old-connectome
   model free-dropped falls within ~1.5 s and then keeps stepping on its side
   (63/63 bursts) while the body wanders in circles - exactly "fell down and
   trying to walk in circles on its side". The old 5.1 s records of marching
   never measured posture; with Root telemetry the free drop was always a
   fall. This confirms, in-model, the standing ruling that the AnimatLab
   originals need the walker/harness analog - free standing is structurally
   out for this geometry.
2. **The tuned winner does not transfer as a working gait**: its te/tf tonic
   drive LATCHES the AnimatLab NonSpiking RG half-centers (the documented
   2026-09-15 tonic-latch doctrine; the 2026-10-16 V1 probe latched already at
   te 2.4/tf 3.0). The MuJoCo/SNS half-centers tolerate the same numeric
   drive - the two stacks' RG parameterizations differ, so only FAMILY RATIOS
   were ever transferable, not the drive itself.
3. **The corrected connectome has no free-running rhythm in AnimatLab** (old
   wiring: 63 bursts; corrected wiring: 0, with or without Ben's 2 nA bias,
   with or without the tuned ratios), while its SNS twin free-runs
   (smoke_w2l.py --template bilateralrg2: period 2.997 s, antiphase r -0.782).
   PRIME SUSPECT FLAGGED FOR BEN: in the corrected aproj the crossed
   V3 -> contralateral RG-ext-IN edges were re-typed onto the SHARED
   "RG Excite" type (2.749 uS) - 27x their old dedicated "V3 Commissural
   Excite" class (0.1 uS) and ~34x the MuJoCo twin's v3_weak family
   (0.08 uS). Strong bilateral crossed excitation is a classic
   sync/latch path in the NonSpiking RG. Second suspect: the F1-flagged
   edge "R Hip MN flx RE -> R Hip MN ext RE" now carries the excitatory
   V3-clone type. Both are Ben's literal aproj state, ported faithfully;
   if he retypes the V3-contra edges back to a weak dedicated class
   (~0.1 uS, or comm-scaled), the free rhythm may return - that is the
   first thing to try before any further AnimatLab tuning.

## 2. THE MAPPING (winner knob -> AnimatLab edit)

Both files edited in lockstep under the SAME GUIDs (except the 2 new dedicated
comm types + new stimuli/connexions, deterministic prefix d71a):

| knob (value) | edit |
|---|---|
| te 5.26781 / tf 7.5507 nA | 4 TonicCurrent stimuli "Gamma drive {L,R} RG {ext,flx}", 0 -> 30.1 s |
| ch_hc 0.47826 (x5.97825 vs 0.08 default) | the 5 weak Heel_Ipsi_* cascade types 0.01 -> 0.0597825 uS |
| ch_ank 0.91176 (x0.91176) | Heel_Ipsi_AnkMNext 0.5 -> 0.45588 uS; Toe_Ipsi_AnkMNext 0.01 -> 0.0091176 uS (both families ride contact_ank in the MuJoCo net) |
| ch_contra 1.10224 (x13.778) | Heel_Contra_RGflx 0.01 -> 0.13778 uS |
| comm 0.70048 | NEW dedicated types "RG Inhibit (c1 commissural, tuned 20261020)" and "RG Excite (V3 contra, tuned 20261020)" at 2.749x0.70048 = 1.92562 uS carrying the 2+2 commissural OUTPUT edges (the shared RG Inhibit/RG Excite types stay 2.749 uS; the RG->{c1,V3} DRIVE edges are NOT comm-scaled - the campaign's gains dict leaves comm_c1/comm_v3 at 1.0). Note the stock "c1 Commissural Inhibit" type is UNUSED in the corrected aproj (the working c1 edges carry shared RG Inhibit) - that is why dedicated clones were needed |
| ctn_amp 4.60959 nA | 4 ContactCount adapter Gain C -> 4.60959 nA (asim 4.60959e-009 A; aproj Value/Scale/Actual triplet) |
| cref 129.28324 N | **NOT TRANSFERRED** (same ruling as the 2026-10-16 report): the aproj adapters are ContactCount-based - no force input, and no force-based SourceDataType exists anywhere in the repo's models. The port runs the contact encoder SATURATED (any contact = full ctn_amp); the winner's sub-saturation linear force scaling does not transfer. FLAGGED for Ben |
| ibnA 0.16671 | all 12 per-muscle Tension (Ib) adapter Gain C scaled x0.16671 (8x 5.062 pA/N -> 8.43886e-13, 3x 13.842 -> 2.3076e-12, 1x 2.963 -> 4.93962e-13; Ben's asymmetric hip values preserved) |
| L RG ext intrinsic 2 nA | the corrected APROJ carries it; the stale asim export had 0 -> asim SYNCED to 2e-009 (the MuJoCo net has no such field - known aproj<->MuJoCo delta, documented) |

## 3. THE ASIM GAP (why the asim needed the corrected-connectome delta)

Ben's two standalone asims predate the whole 2026-10-18 corrected connectome
(0 of the 20 dedicated types). The tuned asim was built from the battle-tested
Ground_Standalone.asim + the exact connexion delta vs the corrected aproj
(computed at build time from w2l_aproj_mine_20261018.json + a live re-mine):
**39 retypes, 6 adds (ankle contact + crossed flexor), 10 deletes (toe
cascade)** - the documented 7-point correction - then the tuned values on top.
Neuron sets are GUID-identical (96/96, names match). Post-build audit: 204/204
connexions resolve to the corrected aproj's types (with the comm
substitutions), all tuned SynAmps present, 4+4 gamma stimuli correct, 4 contact
C and 12 Ib C patched. Charts extended to [0,30] s at 1 ms; new "Root Track"
chart (Root WorldPositionX/Y/Z + hip_L/hip_R JointRotation - the drop
telemetry; WorldPositionX is a proven GUI SelectedDataType and worked headless).

## 4. FILES

In Ben's model folder (`Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\`)
- **all new, untracked; Ben commits or deletes; originals untouched**:
- `Walker_2_Layer_CPG_BilateralRG_tuned20261020.aproj` - the tuned winner,
  GUI vehicle. HEADLESS-RUN ONLY, GUI OPEN-CHECK PENDING BEN.
- `Walker_2_Layer_CPG_BilateralRG_tuned20261020.asim` - matching runnable
  export (corrected connectome + tuned values, 30.1 s).
- `Walker_2_Layer_CPG_BilateralRG_control20261020.aproj/.asim` - Ben's
  literal corrected model as a headless asim (first one that exists; the
  evidence for verdict points 1 and 3).

In this campaign folder (`spinal\reports_20261020\w2l_animatlab\`):
- `W2L_ANIMATLAB_TUNED_DROP_20261020.md` (this report)
- `build_tuned_variant_20261020.py` (tuned + `--control` builder),
  `build_stock30_20261020.py` (old-wiring 30 s reference builder),
  `analyze_run_20261020.py` (drop-protocol analyzer),
  `edits_log_20261020.json` (machine-readable edit list of the tuned build)
- `run_tuned\`, `run_tuned_nodrive\`, `run_control\`, `run_control_no2na\`,
  `run_stock30\` - each: asim copy + 17 chart .txt + console.log +
  analysis.json

## 5. REPRODUCE (egg, cmd)

```
D:\Anaconda\envs\myo\python.exe D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261020\w2l_animatlab\build_tuned_variant_20261020.py
D:\Anaconda\envs\myo\python.exe ...\build_tuned_variant_20261020.py --control
D:\Anaconda\envs\myo\python.exe ...\build_stock30_20261020.py
cd /d "D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin"
AnimatSimulator.exe <run_folder>\<asim>            :: one folder at a time
D:\Anaconda\envs\myo\python.exe ...\analyze_run_20261020.py <run_folder> --json <run_folder>\analysis.json
```

## 6. NOT RUN / NOT CLAIMED

- No GUI open-check of either aproj (window-station-less session; labeled
  pending Ben - same protocol as the 2026-10-16 realization).
- No AnimatLab-side re-tuning was attempted (out of scope: the ask was to test
  the tuned winner, not to search for a new one). The two flagged suspects in
  section 1.3 are the recommended starting points.
- The winner's MuJoCo 30 s result (survive, +1.20 m) lives in the MuJoCo
  yield-rig protocol; the AnimatLab model has no rig/harness anywhere (Root is
  a free body; WalkingPath is a frozen rail; body_L_f/Body_11 etc. are debris
  near the origin, not a shelf) - the two protocols are not comparable beyond
  the shared connectome.
- Ben's 2026-10-17 GUI traces predate the 2026-10-18 corrected-connectome
  commits, so his "circling on its side" observation was on a mid-surgery
  state; the stock30 run here shows the same phenomenology on the clean old
  wiring.

- AnimatLab realizer, 2026-10-20
