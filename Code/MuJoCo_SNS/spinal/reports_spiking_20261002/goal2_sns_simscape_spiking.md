# GOAL 2 — SNS_Simscape spiking campaign (2026-10-02, easteregg2, MATLAB R2025a U1)

Baselines reproduced, spiking capability added ADDITIVELY to `SNS_Library.slx`,
spiking twins of the KneeReflex and BeerCup demos built and compared, one
representative spiking RG→MN subnetwork built and verified (partial alternation,
honestly scoped). The committed demos and units tests reproduce their baseline
numbers EXACTLY with the spiking blocks in the library.

All MATLAB runs: `"D:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch "..."` headless.

---

## 0. Machine + prerequisites

- Host easteregg2, repo `D:\GitHub\Bipedal_Robot`, branch `KneeTestSetup_BenBo_stw`.
- The committed `SNS_Library.slx` + `demos\*.slx` are **R2025b-native** (laptop
  2026-09-22 saves); R2025a REFUSES to load them. Regenerated locally per the
  README's documented route (this modifies those .slx in the git tree — expected,
  listed in §7): `"D:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch "run('D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal2_regen_r2025a.m')"`
  → `sns_build_library` + `sns_build_actuators` (ACTUATOR VALIDATION PASSED,
  max err 4.929e-10 N) + the four demo builders.

## 1. TASK A — baselines (all committed metrics reproduced)

Command: `matlab -batch "run('D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal2_run_baselines.m')"`
(logged values, 2026-10-02; mat: `goal2_baselines.mat`)

| Model | Committed metric | This run | Match |
|---|---|---|---|
| `sns_units_test_2n` | PASS 6.07e-04 mV | PASS, dV_A 6.07e-04 / dV_B 1.22e-04 mV | exact |
| `KneeReflexDemo` | rise 15.8 deg @ 0.05 s; settle 43.5 deg (41.6–44.9) | 15.79 deg; 43.5 deg (41.6–44.9), final 44.87 | exact |
| `KneeReflexCircuit` | runnable twin | identical numbers to the demo | exact |
| `BPACPGLegDemo` | theta range 9.7..48.0 deg | 9.7..48.0 deg (2 hyst switches — the demo is a single sweep, corr(Ae,Af)+0.09, it never cycled) | exact |
| `BeerCupReflexDemo` | ON 2.0 / OFF 7.9 deg sag; A_bi 0.39→0.42; A_tri 0.40→0.32 | ON 2.00 / OFF 7.85 deg; 0.392→0.417; 0.400→0.316 | exact |

## 2. TASK B — spiking blocks (ADDITIVE) + verification

### 2.1 New blocks (git history check: NO pre-redesign spiking synapse existed anywhere in the 8 committed versions of `sns_build_library.m` — only the `SpikingLIFNeuron`, which stays untouched)

Added by **`sns_build_spiking.m`** (opens the existing library, appends, never
touches committed blocks; idempotent; run AFTER `sns_build_library` +
`sns_build_actuators`, same pattern as `sns_build_actuators.m`). Library now
14 blocks. Semantics ported from sns-toolbox 1.5.2 `connections.py` /
`backends.py` (Szczecinski et al. 2017/2020; the doctrine the library header cites):

1. **`SpikingNeuron`** — spiking LIF with INTERNAL synaptic summation (the
   committed `SpikingLIFNeuron` has only an Iapp port; this one is the
   toolbox-faithful drop-in with `syn1..syn6` ports like `NonSpikingNeuron`).
   Between spikes the membrane is the identical RC circuit; `V >= Vth` → spike
   output + reset to `Vreset` (fixed threshold = toolbox m=0, increment=0;
   adaptive threshold documented as future work). Ports: Iapp, syn1..6; out V, spike.
2. **`SpikingSynapse`** — spiking chemical synapse, event-wired (input = the
   presynaptic SPIKE line). `dg/dt = −g/tau_syn` between events; per spike
   (rising edge) `g ← min(gmax, g + ginc)` — the toolbox `conductance_increment`
   semantics, saturating at gmax. Output `[g; g·Esyn]` plugs into ANY neuron's
   syn port (spiking or non-spiking).
3. **`HybridSpikingSynapse`** — the spiking→NON-SPIKING synapse (Tutorials 3/5/9
   pattern: spikes above, analog voltages into the plant). Same conductance
   dynamics, but input = presynaptic **V** line (wired exactly like a
   `NonSpikingSynapse`); the synapse detects the spike itself at
   `Vpre >= ThrPre` (AnimatLab semantics). Its mask documents the ThrPre trap
   (see 2.3).

Implementation notes (probes `goal2_probe*.m` in this folder):
- Integrator port order with external IC + rising reset is **xdot=1, reset=2,
  IC=3** (all 6 permutations tested; correct row's signature 5.000→4.304→4.950→4.304).
- The threshold→reset and IC←state wirings are **algebraic loops**; both are
  broken with the Integrator **STATE PORT** (`ph.State` handle — it is not a
  numbered outport), Simulink's documented LIF pattern.
- Units mapping documented in the script header: toolbox `reversal_potential`
  (relative) = block `Esyn − Vrest_post`; `C[uF] = Cm[nF]/1000`; `tau_syn` in ms.

### 2.2 Verification — `sns_units_test_spiking.m` (NEW; PASS)

Command: `matlab -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape'); sns_units_test_spiking"`

| Test | Check | Result |
|---|---|---|
| 1a | single spike → g = ginc·exp(−(t−ts)/tau_syn) exact | max dev **6.246e-07 µS** |
| 1b | saturating accumulation over a 50 ms pulse train (gmax 0.55, ginc 0.5) | max dev **4.630e-05 µS** |
| 2 | hybrid 2-neuron circuit (SpikingNeuron A → HybridSpikingSynapse → NonSpikingNeuron B) vs a MATLAB Euler reference implementing the toolbox update equations VERBATIM (backends.py:131–181, dt=1e-6 s) | spikes **62/62**, max spike-time dev **11.75 µs**, max V_B dev **1.21e-02 mV**, windowed mean-g dev **4.03e-04 µS** |
| 3 | spiking→spiking chain (A → SpikingSynapse → slower SpikingNeuron C) | A 37 / C 14 spikes in 0.3 s |

**Cross-platform check vs the ACTUAL sns-toolbox numpy backend**
(`goal2_cross_check_spiking.py`, `D:\Anaconda\envs\myo\python.exe`): same circuit,
toolbox dt=1e-6 s:

```
toolbox numpy backend: 62 spikes, ISI 8.0470 ms, V_B end 1.71012 mV
Simulink             : 62 spikes, ISI 8.0472 ms, V_B end 1.71093 mV
diffs: count MATCH, |mean ISI| 0.00019 ms, |V_B end| 0.00081 mV  -> CROSS-CHECK PASS
```

(Simulink's ISI 8.0472 ms equals the analytic value 5·ln(5) exactly; the toolbox
backend at dt=1e-5 carries its own 7.2 µs/spike Euler bias — measured, not a
block error.)

### 2.3 Traps banked (each cost a debug cycle)

- **pchip interpolation across jump discontinuities overshoots** (solver logs a
  point exactly AT each event) — use `linear` and stay ≥10 ms from jumps in
  event-sampled comparisons.
- **A spiking presynaptic V output never renders values ≥ Vth** (it resets at
  threshold): `HybridSpikingSynapse.ThrPre` must sit BETWEEN Vrest_pre and Vth_pre
  and comfortably BELOW Vth_pre. ThrPre = Vth fires NOTHING (verified: mean-g dev
  0.2885 = full cap); ThrPre ≈ 0.94·Vth keys the jump ~1 ms before the true spike.
- The toolbox is **unit-agnostic**: with dt in seconds, `time_constant` must be
  in seconds (my first cross-check passed 100 = 100 s and pinned g at the cap —
  diagnosed from the 1.846 mV asymptote, fixed to 0.1 s).
- Mask default `ThrPre = −55` (below Vrest −52) saturates graded synapses ON at
  rest — the demos override to −45; my first RG build pinned the half-centers
  at −63 mV until I matched the demo convention.

### 2.4 Regression (the ask's gate)

Re-ran §1's full baseline suite + both units tests with the spiking blocks in
the library — **every number identical** (§1 table + `SPIKING UNITS TEST: PASS`
in the same batch). Library loads; `sns_units_test_2n` PASS 6.07e-04 mV; all
four demos reproduce exactly.

## 3. TASK C — spiking demo twins (NEW `*_Spiking.slx`; Ben's demos untouched)

Both twins are COPY-AND-PATCH builds: `save_system` a copy of the committed
demo, then surgically replace ONLY the sensory-neuron + synapse layer (same
knee/elbow model, afferents, MNs, muscles, logging). MNs stay NON-SPIKING
(hybrid doctrine); afferent currents drive SPIKING interneurons whose
`HybridSpikingSynapse` outputs land on the same MN syn ports.

### KneeReflexDemo_Spiking (`sns_build_knee_spiking.m` / `sns_run_knee_spiking.m`)

Gains matched by ḡ = ginc·f·tau_syn ↔ baseline gmax·Sat(Vpre) (3 hand
iterations; final ginc: exc_ext 0.087, exc_flex 0.126, Ib 0.078, recip
0.108/0.072 µS/spike, tau_syn 10 ms, gmax 0.3; documented in the builder):

| Metric (5 s run) | Baseline | Spiking twin |
|---|---|---|
| theta at 0.05 s (rise) | 15.79 deg | **15.79 deg** (identical) |
| settle mean, last 2 s | 43.5 deg | 44.4 deg |
| settle range | 41.6–44.9 deg | 41.9–46.7 deg |
| A_ext / A_flex (window means) | 0.53 / 0.53 | 0.58 / 0.56 |
| antagonist alternation | ~0.5–2 Hz (README) | 7.4 Hz dominant peak |
| IN firing at settle | (graded, subthreshold V) | Ia_ext 95 Hz; Ib 0 Hz |

Honest differences: the spiking twin settles 0.9° higher with a slightly hotter,
faster antagonist alternation — the rate-coding adds synaptic filtering
(tau_syn) that the graded synapses do not have (tau_syn 30 ms produced a 9°
4.7 Hz limit cycle; 10 ms ships). `Ib = 0 Hz` matches the baseline (Ib afferent
current ≈ 0.9 nA is below threshold in BOTH models — that pathway is a
transient brake during the rise only).

### BeerCupReflexDemo_Spiking (`sns_build_beer_spiking.m` / `sns_run_beer_spiking.m`)

Equilibrium start carries over EXACTLY (at t=0 all afferent currents ≈ 0 → INs
below threshold → ḡ = 0, same as the baseline's Sat = 0; same desc drives, same
A0). Gains: ginc = 1.14 × baseline gmax (the active Ia_biceps graded synapses
are saturated during sag; Ia_bi fires ~114 Hz at 2° sag; kReflex-scaled exactly
like the baseline so ON/OFF A/B works).

| Metric (10 s pour) | Baseline | Spiking twin |
|---|---|---|
| max sag after 2 s, reflex ON | 2.00 deg | **2.02 deg** |
| max sag after 2 s, reflex OFF | 7.85 deg | **7.85 deg** (identical) |
| reflex modulation (OFF−ON) | 5.85 deg | 5.83 deg |
| A_bi over the pour (ON) | 0.392→0.417 | 0.392→0.447 |
| A_tri over the pour (ON) | 0.400→0.316 | 0.400→0.254 |
| Ia_bi firing (ON) | (graded) | 166.6 Hz mean |

The spiking reflex holds the cup at essentially the same level; its biceps
response and antagonist inhibition are stronger (first-order rate-coding
nonlinearity vs graded Sat).

## 4. TASK D — SNS_SpinalNetwork assessment + representative subnetwork

**Full 410-neuron conversion: NOT attempted, scoped as future work.** Scale
(from README + `sns_build_from_json.m`): 410 NonSpikingNeurons, 1392 synapses,
376 input ports, 92 MN outputs, 146 chained SynSums (MN in-degree to 17),
~2400 blocks. Every one of the 1392 gains was tuned for GRADED voltages
(MN drive S(V) = clip((V−Thr)/Slope), V in mV above −55); a spiking conversion
is a full re-tune, not a parameter port — the 2026-09-13 session proved this
network sits near a bifurcation where even BLAS summation-order changes the
trajectory. Building it un-tuned would produce a non-walking model and prove
nothing.

**Representative subnetwork BUILT + verified (partial): `SNS_SpikingRG_MN.slx`**
(`sns_build_rgmn_spiking.m` / `sns_run_rgmn_spiking.m`, 12 s run):
spiking RG half-centers (drives 4.5/2.5 nA) with IN-LAMINATED mutual
inhibition (Deng/Biped_2xCPG convention: RG spikes → hybrid synapse → fast
non-spiking IN → graded inhibition onto the other RG) + spike-rate Adp cells
(slow non-spiking, hybrid input) for burst termination, driving a 3-MN
non-spiking pool through hybrid synapses (tau_syn 20 ms).

- **VERIFIED — the hybrid motif works**: the MN pool voltage is GRADED,
  modulating over a ~5 mV depth (−52.0..−46.9 mV) and tracking its RG's rate
  envelope (corr +0.2..+0.4); MN_F1 tracks RG_F. Spikes above, analog
  voltages into the plant — demonstrated end to end.
- **PARTIAL — alternation**: antiphase rate modulation (envelope corr −0.49),
  RG_E bursts 4×/12 s, but RG_F keeps a low-rate floor through E's bursts
  (no discrete F bursts). Diagnosis: direct spike-to-spike inhibition could
  not hold the loser below threshold (4–8 Hz co-firing, measured); the
  laminated route fixed dominance but the Adp engagement threshold vs actual
  firing rates still needs calibration for clean burst handoff.
- 8 parameter iterations are recorded in the builder comments (what failed
  and why) — the honest engineering trail.

## 5. Claims + exact repro commands

| # | Claim | Evidence | Repro |
|---|---|---|---|
| 1 | All committed baselines reproduce on R2025a after the documented local regeneration | `goal2_baselines.mat` + §1 table | `matlab -batch "run('D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal2_run_baselines.m')"` |
| 2 | 3 spiking blocks added additively; committed tests/demos still pass EXACTLY | §2.4 final regression batch | `matlab -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape'); sns_build_spiking; sns_units_test_spiking"` + claim 1's command |
| 3 | Spiking units test PASS: analytic decay 6.2e-07 µS; toolbox-equations match 62/62 spikes, 11.75 µs, 1.2e-02 mV | `sns_units_test_spiking` printout | claim 2's command |
| 4 | Cross-platform match with the actual toolbox numpy backend: ISI diff 0.19 µs, V_B diff 0.81 µV | `goal2_cross_check_spiking.py` printout | `D:\Anaconda\envs\myo\python.exe D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal2_cross_check_spiking.py` |
| 5 | Knee spiking twin: rise identical (15.79°), settle 44.4° vs 43.5° | `goal2_knee_spiking.mat/.png` | `matlab -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/demos'); sns_build_knee_spiking; sns_run_knee_spiking"` |
| 6 | Beer spiking twin: ON sag 2.02° vs 2.00°, OFF 7.85° identical | `goal2_beer_spiking.mat/.png` | `matlab -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/demos'); sns_build_beer_spiking; sns_run_beer_spiking"` |
| 7 | RG→MN hybrid motif verified (graded MN V, ~5 mV depth, tracks RG rate); burst alternation partial | `goal2_rgmn_spiking.mat/.png` | `matlab -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/demos'); sns_build_rgmn_spiking; sns_run_rgmn_spiking"` |

(`matlab` = `"D:\Program Files\MATLAB\R2025a\bin\matlab.exe"`.)

## 6. Not done / out of scope (honest)

- Full 410-neuron SNS_SpinalNetwork spiking conversion — future work (§4).
- Adaptive-threshold spiking neuron variant (toolbox m, theta_increment,
  threshold_floor) — documented in the block description as future work; the
  shipped block is the fixed-threshold (m=0) case.
- Clean burst alternation in the RG→MN subnetwork (§4 partial).
- Transmission delay on spiking synapses (toolbox `transmission_delay`, integer
  steps) — not implemented; a Transport Delay on the spike line covers it if
  ever needed. Noted in the block docs.
- Cosmetic: a `Parameter precision loss` warning fires for `ginc` gains
  (ufix32 quantization, abs err ~6e-12) — numerically negligible, not chased.
- The AnimatLab-side spiking baselines in this same folder (metrics_*.json,
  goal3_animatlab_spiking.md) belong to OTHER goals/agents, not this report.

## 7. File inventory

**Pre-existing files MODIFIED (all by re-running their own committed builders —
no hand edits):**
- `Code/Matlab/SNS_Simscape/SNS_Library.slx` — R2025a regeneration + 3 spiking
  blocks appended (git will show it modified; the laptop's R2025b-native copy is
  what's committed).
- `Code/Matlab/SNS_Simscape/demos/{KneeReflexDemo,KneeReflexCircuit,BPACPGLegDemo,BeerCupReflexDemo}.slx`
  — R2025a regeneration only (required: R2025a cannot load R2025b models).

**NEW (this goal):**
- Library: `sns_build_spiking.m`, `sns_units_test_spiking.m` (SNS_Simscape root).
- Demos: `sns_build_knee_spiking.m`, `sns_run_knee_spiking.m`,
  `sns_build_beer_spiking.m`, `sns_run_beer_spiking.m`,
  `sns_build_rgmn_spiking.m`, `sns_run_rgmn_spiking.m`,
  `KneeReflexDemo_Spiking.slx`, `BeerCupReflexDemo_Spiking.slx`,
  `SNS_SpikingRG_MN.slx` (all R2025a-native).
- Results: `results/sns_units_spk_hybrid.slx`, `results/units_ref_spiking_sim.mat`,
  `results/sns_units_2n.slx` (test artifacts).
- This folder: `goal2_*.{m,py,mat,png}` incl. the four probe scripts + RG-MN
  debug script (kept as the debugging record).

**NOT this goal's changes** (other agents in the same campaign, present in the
tree): `Code/MuJoCo_SNS/spinal/build_network.py` (modified),
`build_network_spiking.py`, `calibrate_spiking.py`, `spiking_calibration.json`,
`spinal_run_spkbase_*.npz`, and the goal-1/3 files in this folder
(`metrics_*`, `baseline_graded_ranges.json`, `goal3_animatlab_spiking.md`,
`models/`, `tools/`, `logs/`).
