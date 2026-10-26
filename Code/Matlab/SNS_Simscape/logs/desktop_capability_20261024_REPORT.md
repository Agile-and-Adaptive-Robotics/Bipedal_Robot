# SNS_Simscape desktop capability pass — 2026-10-24 (ZCode, EB475WS4)

Ben's directive: run the SNS-Simscape-Simulink stack to the best of this
machine's ability while the laptop (the only Simscape Multibody license
seat) is down. Evidence levels: everything below is dynamically reproduced
today unless marked otherwise.

## Capability matrix (fresh evidence, 2026-10-24)

| Capability | EB475WS4 (R2025a U1) | easteregg2 (R2025a U1) | laptop (R2025b) |
|---|---|---|---|
| Base Simulink SNS stack | **VERIFIED today** (see below) | documented 2026-09-10 | documented |
| Simscape Multibody checkout (smimport) | blocked (documented) | **blocked — re-verified today**: `sns_urdf_smoke.m` fails at the MechanismConfiguration license callback; both `license('test')` = 0 | works (machine down) |
| mex C compiler | **MinGW64 configured** → MuJoCo bridge rebuild feasible | untested | bridge built there |

## What ran on EB475WS4 today (headless, `matlab -batch`)

Logs: `desktop_capability_20261024_a.log` / `_b.log` in this folder;
runner `dev\run_desktop_capability_20261024.m`.

- `sns_build_library` + `sns_build_actuators`: SNS_Library.slx regenerated
  R2025a-native, 11 blocks (2026-09-22 architecture, side-view icons).
- `sns_test_actuators`: **PASS** — BPA max err 4.929e-10 N, BioMuscle
  1.421e-14 N (documented values reproduced).
- `sns_units_test_2n`: **PASS** — 6.07e-04 mV (documented value exact).
- Demos (fresh processes): KneeReflex theta final 44.87°, A_ext 0.74 /
  A_flex 0.30 (documented settle band); BPACPG theta range 9.7–48.0°
  (identical to committed results); BeerCup ON 2.00° vs OFF 7.85° sag,
  triceps 0.400→0.316 (documented values).
- `sns_build_from_json` + `sns_verify_from_json` on the 410-neuron
  production export: **PASS — all 410 neurons, max dev 4.201e-06 mV,
  V_DRIVE Euler-exact to 9 decimals** (documented value reproduced).
- Winner slice R8 run: **ALL 4 CHECKS PASS** (period 1.7285 s = numpy ref;
  GAMMA duty 0.32→0.66; ALPHA MN 2.055/2.591 mV; heel/toe entrainment).
- `mex -setup C`: MinGW64 configured.

## Known issues found today (not capability gaps)

1. **In-session demo-run artifact**: running a demo immediately after the
   library rebuild in the SAME session fails with Simulink's "Error due to
   multiple causes"; the same demo passes in a fresh process. Fresh-process
   execution is already program doctrine — use it.
2. **Deng CPG — definitive status (probed after the pass).** The builder
   `sns_build_deng_cpg.m` is broken on R2025a (`Invalid Simulink object
   name: 'pair1_s1/2'`, a generated name — the script has evolved past the
   committed .slx set and deletes-then-fails, removing
   `demos\SNS_Deng_RG.slx` when run; restored byte-identical to git HEAD).
   The committed pre-built trio does NOT compile on R2025a either: probe
   (`dev\probe_deng_ports_20261024.m`) shows RG (stim in; V_ext, V_flx
   out) and PF (V_ext, V_flx in; 4 out) are correct 1-wide, but inside
   SNS_Deng_CPGDemo the RG reference presents a 2-wide Output Port 1 into
   PF's 1-wide input, and the reference resolves to the subsystem path
   'SNS_Deng_CPGDemo/RG'. The September verification (period 1.938 s vs
   numpy 1.94 s, corr −0.86) stands as a laptop-era record. Repair =
   bounded task: fix the builder's generated port names on R2025a,
   regenerate the trio, rerun. Today's partial library rebuild is
   preserved at `dev\rebuilt_deng_library_20261024.slx`; the September
   library was restored from git HEAD.
3. **Winner-slice builder path bug**: `sns_build_winner_slice.m` reports
   `SNS_Library_R2025a.slx not found` (the file is in SNS_Simscape root;
   script expects another location). The RUN verification against the
   existing slice is unaffected and passes.

## FIXED + BUILT (2026-10-24/25 overnight, ZCode EB475WS4)

- **Deng trio REBUILT and re-verified** (resolves item 2 above): root
  cause was architectural — the trio predated the 2026-09-22 one-input-
  synapse redesign (Vpost sense line into a port that no longer exists)
  plus one dropped wire (dVsum→memb) found by internal-current probing.
  `sns_build_deng_cpg.m` rewritten for the current architecture; HCNeuron
  now has real syn1..syn4 ports with internal summation and a separate
  `hPowExp` mask param (0.5 = Deng h^Kh semantics, bit-identical; 1 =
  toolbox Tutorial-8 semantics). **Verification: 10 bursts / 20 s, period
  1.938 s (documented value exact), ext-vs-flx correlation −0.857
  (documented −0.86).**
- **Li model built in base-Simulink** (`demos\sns_build_li.m` →
  `demos\SNS_Li_CPG.slx`): all 58 nodes and all 98 connections from the
  GUID-level mine, every synapse carrying ITS OWN per-synapse conductance
  (1–8 µS) per Ben's per-target ruling. Stance CPGs are a LABELED
  adaptation (graded cell + slow Adp, tau = documented tau_Ca 7.5 s);
  swapping in the spiking cell is local to build_stance_cpg(). Smoke run
  (`dev\run_li_smoke_20261024.m`): both stance CPGs oscillate
  (V −52..−16.8 mV, period ≈ 2.0 s per side), inter-leg correlation
  −0.271 under the scripted antiphase heel script.
- **Tutorial-8 Simulink twin: PASS** (`dev\sns_tutorial8_twin_20261024.m`
  + `results\tutorial8_ref.mat`): winner/loser structure and steady
  states EXACT (20.000 / −0.010 mV), loser pointwise 0.595 mV, winner
  excursion within 2.2 mV (7%) and 5 ms of the toolbox reference. Two
  toolbox semantics pinned empirically along the way (probes
  `probe_tutorial8_semantics/_synapse`): the persistent-Na neuron rests
  at Er+R by construction (Gna formula), and the tutorial pair is a
  winner-take-all latch (kick timing decides the winner; a sustained
  1 s drive INVERTS the outcome via h-inactivation — the AnimatLab
  RG-latch lesson, independently reproduced in the twin). Toolbox
  synapse defaults: e_lo/e_hi 0/20; h exponent 1 with Kh as sigmoid
  coefficient (vs Deng's h^Kh).
- **Deng 2019's ACTUAL AnimatLab project ACQUIRED** (ModelDB 264844,
  25.6 MB zip → `dev\deng2019_modeldb\`): the paper's K_A-synergy + RG
  .aproj set — the definitive per-synapse source for the per-target
  Deng rebuild (Ben's order). Mining + headless AnimatLab run = next.

## 2026-10-25 morning additions (ZCode EB475WS4)

- **Deng 2019 paper project MINED** (`dev\mine_deng2019_20261025.py` →
  `deng2019_mine_20261025.json`): 90 non-spiking neurons, 12 muscles,
  12 stretch receptors, 196 synapse edges, 58 named synapse types with
  DISTINCT per-path conductances — the paper itself is per-target:
  PF→MN per joint = hip flx 3.6325 / hip ext 2.6563 / knee flx 1.5165 /
  knee ext 4.9296 / ankle flx 4.5222 / ankle ext 4.0538 µS. Doctrine
  evidence: 114 of 196 edges carry the UNTOUCHED 0.5 µS per-instance
  GUI default while the real values live at type level
  (`dev\check_deng2019_instances_20261025.py`) — for AnimatLab-sourced
  models the TYPE conductance is the effective strength; the mine
  records BOTH per Ben's ruling. My Simulink Deng twin updated to the
  paper's hip values (3.6325 / 2.6563 µS) and re-verified: period
  1.938 s, corr −0.857, unchanged.
- **Tutorial-2 Simulink twin: PASS**
  (`dev\sns_tutorial2_twin_20261025.m` + `results\tutorial2_ref.mat`):
  5-neuron exc/disynaptic-inh/modulatory network — ALL five finals
  match the toolbox numpy backend to ≤1.1e-05 mV. Toolbox default
  neuron confirmed Gm 1 µS / Cm 5 nF (τ 5 ms) by hand-derivation from
  the reference finals.
- **Air-stepping model (curriculum stages 1–2) and v2 leg/cup
  animations**: in progress at report time (`demos\sns_build_air_step.m`
  + `demos\sns_animate_demo_v2.m`).
- **AnimatLab Human-Stretch-Reflex tutorial project: DEAD LINK** —
  animatlab.com's `Arm_Flexion_Model.zip` 404s at every URL variant; the
  video remains at /Portals/0/Media/HumanStretchReflex.mp4. ModelDB's
  Szczecinski stretch-reflex entry is the remaining fetch path.

## Files changed on this machine (for Ben's commit decision)

- `SNS_Library.slx`, `SNS_Library_R2025a.slx`, demo .slx files, Deng
  .slx files, `results\SNS_SpinalNetwork.slx`: regenerated R2025a-native
  (the documented path for R2025a machines). Pre-pass backups:
  `dev\backup_pre_20261024\`.
- New: `dev\run_desktop_capability_20261024.m`, `dev\debug_demo_fail_20261024.m`,
  this report.

## Next unlock (sanctioned by capability, not yet run)

Rebuild the MuJoCo↔Simulink bridge mex on this desktop (MinGW64 present;
win64 MuJoCo libs ship in `mujoco_bridge\blockset\...\lib\win64\mujoco`).
Then SNS Simulink networks can drive MuJoCo physics here without the
laptop.
