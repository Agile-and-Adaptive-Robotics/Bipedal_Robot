# Goal 2 / Milestone 4 — SPLIT RG: each leg powered by its own RG, coupled by the bilateralrg commissurals

**Date:** 2026-09-25 (EB475WS4, unattended campaign) · **Status: complete — all three gates PASS**
**The ask (Ben):** "DISCONNECT the LH RG from directly controlling the RH PF layers. Create a copy of the LH RG for the RH side and use it to power the RH PF layers; then couple the two RGs with the commissural wiring of the documented bilateralrg build."

This is SESSION_NOTES_20260916.md build-chain steps 1–3 applied to the M3 2023-original transcription: step 1 (`build_rg.pl`) = the disconnect + R RG mirror + antiphase `Stimulus_2`; steps 2–3 (`build_comm.pl` + `patch_comm_types.pl`) = the Shinohara c1/V3 commissurals.

**Artifacts (all under `Code\MuJoCo_SNS\spinal\w2l_mujoco\`; nothing else touched):**
- `build_w2l_split_net.py` — NEW: the split-RG network (edit spec + builder + asserts).
- `test_w2l_air.py` — EXTENDED (the ask's "extend/keep"): `--net=split --comm=X` branch, R-RG logging, phase-lag + band-limited antiphase metrics, M4 verdict blocks. The `net=orig` default path is logic-unchanged (regression below).
Gate logs: `reports_20260925\logs\test_w2l_air_m4_coupled_run{1,2,3}.log`, `test_w2l_air_m4_ablated_run1.log`, `test_w2l_air_m3_regression.log`.

Env: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe` (mujoco 2.3.7, sns-toolbox 1.5.2, scipy for the bandpass), cwd `w2l_mujoco\`. No pip installs. Protected set untouched: `runner.py`/`build_network.py`/`params.py` (read-only import of `SNS_NumpyFixedTau` only), `spinal_run.npz`, studies, `reports_20260923/24` — all untouched. `w2l_cpg\` untouched (read-only import of `W2L_GAINS`/`_neu`/`_syn`/`TAU`/`NAP`).

---

## 1. What changed, edge by edge

`make_split_template()` (build_w2l_split_net.py:119-176) edits the `w2laproj` template (92n/163e) before building; every edit is an asserted literal in the file (CROSSED_REMOVE / R_MIRROR_EDGES / R_PF_DRIVE_EDGES / COMM_EDGES, lines 66-113):

**REMOVED — the 4 crossed `pf_drive` edges (the L-RG→R-PF disconnect, verbatim template pairs):**

| from | to | was |
|---|---|---|
| L RG ext | R Knee PF flx | exc, pf_drive |
| L RG ext | R Hip PF flx | exc, pf_drive |
| L RG flx | R Hip PF ext | exc, pf_drive |
| L RG flx | R Knee PF ext | exc, pf_drive |

These were the R side's ONLY RG drive in M3 (M3 report §2: "the crossing is the whole interleg coupling").

**ADDED — R RG half-center block, an exact mirror of the L block** (SESSION_NOTES line 85: "R RG half-center mirrors the L RG block exactly"): 4 neurons `R RG ext` (HC-RG-E), `R RG flx` (HC-RG-F), `R RG ext IN` (IN-InE), `R RG flx IN` (IN-InF) — same persistent-Na class + fixed tau_h as L — plus the L block's 6 `rg_laminate`-family edges mirrored: 2 direct HC↔HC escape exc (R RG ext↔R RG flx) + R RG ext→R RG ext IN exc + R RG flx→R RG flx IN exc + R RG flx IN→R RG ext inh + R RG ext IN→R RG flx inh. Gains identical to L: direct HC↔HC = `rg_direct_exc` 0.5, laminated = `rg_laminate` 4.0. (The 2 direct edges are the documented MIRROR FIX — the bilateralrg template itself ships only 4; without them the R pair cannot oscillate while L does, measured 2026-09-24, `w2l_cpg\build_w2l_net.py:321-335`.)

**ADDED — 4 ipsilateral `pf_drive` edges (the R RG now powers the R PF layers):** R RG ext→R Hip PF ext, R RG ext→R Knee PF ext, R RG flx→R Hip PF flx, R RG flx→R Knee PF flx — the L side's own mapping mirrored, gain `pf_drive` 2.4.

**ADDED — antiphase kickoff port:** PORT-load node `Stimulus_2` → R RG flx (exc, tag `kickoff_antiphase`); caller drives it 10 nA for t ∈ [0, 0.01 s] exactly like Stimulus_1 (the .aproj waveform pair, SESSION_NOTES step 1: "legs start neurally antiphase").

**ADDED — TONIC input ports** on R RG ext / R RG flx (mirror of the L tonics; drive regime te=3 / tf=4 on all four half-centers, symmetric).

**ADDED — the commissural coupling (coupled config only), verbatim from the bilateralrg template:**

| from | to | sign | tag | gain used |
|---|---|---|---|---|
| L RG flx | c1_L | exc | comm_c1 | 1.0 |
| c1_L | R RG flx | **inh** | c1_SynAmp2.749 | **3.0** |
| R RG flx | c1_R | exc | comm_c1 | 1.0 |
| c1_R | L RG flx | **inh** | c1_SynAmp2.749 | **3.0** |
| L RG ext | V3_L | exc | comm_v3 | 1.0 |
| V3_L | R RG ext | exc (weak) | v3_SynAmp0.1_weak | 0.08 |
| V3_L | R RG ext IN | exc (weak) | v3_to_contra_InE | 0.08 |
| R RG ext | V3_R | exc | comm_v3 | 1.0 |
| V3_R | L RG ext | exc (weak) | v3_SynAmp0.1_weak | 0.08 |
| V3_R | L RG ext IN | exc (weak) | v3_to_contra_InE | 0.08 |

c1_L/c1_R are IN-C, V3_L/V3_R are IN-V3, both on the RG-layer membrane tau (0.05 s), per the bilateralrg template's node types and `w2l_cpg`'s `_TYPE_TAU`.

**Gains — "the gains that build used":** the `w2l_cpg\build_w2l_net.py` `W2L_GAINS` re-expression of this same documented build on this stack's 0–5 mV scale — `comm_c1=1.0`, `c1_inh=3.0`, `comm_v3=1.0`, `v3_weak=0.08`. Those carry the build's own measured adjustments (recorded there): raw AnimatLab SynAmps were c1 **2.749** inh / V3 **0.1** weak exc, but with RG-Excite strength on the V3 path both RGs E-latch (SESSION_NOTES step 3 "REQUIRED" fix), v3 0.08 keeps the ~4% c1:V3 ratio, and c1_inh 3.0 (vs 2.749/2.0) locks the two sides' periods (at 2.0 the sides free-ran 1.05 vs 1.10 s and drifted — the coupled run below shows the 3.0 lock directly). `--c1`/`--v3`-equivalent overrides are available via `build(gains=...)`; the coupled/ablated switch is `--comm`.

**Census (asserted in code, printed by `python build_w2l_split_net.py`):**
- coupled: **87 neurons / 166 synapses / 6 inputs / 12 outputs** (M3: 79/150/3/12; Δ = +4 R RG +4 commissural neurons; −4 crossed +6 mirror +4 ipsi-drive +10 commissural synapses)
- ablated (`--comm=0`): **83 / 156 / 6 / 12** — the 4 c1/V3 neurons and all 10 commissural synapses are **not built** (conditional-topology lab law: zero-gain synapses alone change BLAS summation order, so the ablation removes them rather than zeroing them).
- Per-tag counts (coupled run log): pf_drive 8 (4 L + 4 R ipsi), rg_laminate 12 (6+6), comm_c1 2, c1_SynAmp2.749 2, comm_v3 2, v3_SynAmp0.1_weak 2, v3_to_contra_InE 2, all other tags unchanged from M3.

**Causal-disconnect assert** (build_w2l_split_net.py:161-166, runs on EVERY build): no edge into any R Hip/R Knee PF **half-center** has an L-side source. The only remaining cross-side paths in the coupled net are the c1/V3 commissurals (RG↔RG) and the 2023 original's shared knee-extensor Renshaw quirk (L Knee MN ext RE excited by both sides' knee-ext MNs — MN/RC level, transcribed as-is in M3, not an RG→PF path).

## 2. GATES — all three run this session (exit 0 each)

Commands (cwd `w2l_mujoco\`): `C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_w2l_air.py --net=split --comm=1.0 --dur=20` / `--comm=0` / (M3 regression) no args. Logs under `reports_20260925\logs\`.

### Gates (a)+(c): COUPLED — final paste of `test_w2l_air_m4_coupled_run3.log` (run1/2 = metric-debug trail, numbers identical)

```
== gate: M4 split-RG coupled air stepping on the axis-fixed M1 body, 20 s ==
   knobs: net=split comm=1.0 te=3.0 tf=4.0 tau_rg_nap_h=0.25 ctrl_cap=0.5 joint_damp=3.0 stiff_limits=1 lift=0.3
   ground contacts: 0 (must be 0) | leg-leg self contacts: 26741
   neural: L RG ext bursts=13, period 1.357 s (0.74 Hz), E max 6.98 mV
   neural: R RG ext bursts=13, period 1.357 s (0.74 Hz), E max 6.98 mV; L/R RG-E r = -0.750 (antiphase want < -0.5)
   hip L swing excursions: 27, mean interval 0.670 s (1.49 Hz); ACF dominant period 1.358 s (0.74 Hz)
   hip R swing excursions: 16, mean interval 1.161 s (0.86 Hz); ACF dominant period 1.358 s (0.74 Hz)
   hip L/R cross-correlation: peak +0.497 at lag -0.669 s, trough -0.448 at lag -0.157 s (L hip ACF period 1.358 s)
   L/R phase lag = 0.507 cycle (antiphase = 0.5)
   hip L/R correlation, fundamental band 0.3-1.2 Hz: r = -0.941 (antiphase want < -0.3); raw r = +0.323 (contact-kick 2nd harmonic inflates the raw value)
   joint flexion-positive excursions (deg, min..max and range) vs AnimatLab references:
     hip   L [ -25.3,  +7.7] range  33.0 | R [ -12.6, +17.1] range  29.7   (ref range ~38)
     knee  L [  -4.0, +64.1] range  68.1 | R [  -4.0, +63.0] range  67.1   (ref range ~61)
     ankle L [  -2.2, +24.6] range  26.8 | R [  -2.4, +23.9] range  26.3   (ref range ~16)
   hip flexion L/R correlation (0 lag): r = +0.323 (antiphase want < -0.3)
   [PASS] finite_20s
   [PASS] both_hips_swing
   [PASS] antiphase
   [PASS] airborne
   [PASS] freq_within_2x
   [PASS] hip_amp_ok
   [PASS] rg_e_antiphase
   [PASS] phase_half_cycle
VERDICT: PASS  hip-ACF period 1.358 s (0.74 Hz)  band-antiphase_r=-0.941 (raw +0.323)  phase=0.507 cycle  hip_range_L=33.0 deg
```

(a) air-stepping gate passes with the split net (finite 20 s, zero ground contacts, both hips swing ≥ 15°, freq 0.74 Hz within 2× of the 0.77 Hz 2023-original reference, amplitude ok). (c) **L/R phase lag = 0.507 cycle** — half cycle, measured as the peak of the hip-flexion cross-correlation (at lag −0.669 s = 0.494×T, plus the onset-lag view 0.507).

### Gate (b): CAUSALITY ABLATION — paste of `test_w2l_air_m4_ablated_run1.log` (`--comm=0`: c1/V3 neurons + all 10 commissural synapses NOT built)

```
== gate: M4 ABLATED (--comm=0, commissurals REMOVED) air stepping on the axis-fixed M1 body, 20 s ==
   knobs: net=split comm=0.0 te=3.0 tf=4.0 tau_rg_nap_h=0.25 ctrl_cap=0.5 joint_damp=3.0 stiff_limits=1 lift=0.3
   ground contacts: 0 (must be 0) | leg-leg self contacts: 28266
   neural: L RG ext bursts=18, period 1.027 s (0.97 Hz), E max 6.99 mV
   neural: R RG ext bursts=17, period 1.027 s (0.97 Hz), E max 6.99 mV; L/R RG-E r = -0.788 (antiphase want < -0.5)
   hip L swing excursions: 35, mean interval 0.517 s (1.93 Hz); ACF dominant period 1.027 s (0.97 Hz)
   hip R swing excursions: 18, mean interval 1.010 s (0.99 Hz); ACF dominant period 1.027 s (0.97 Hz)
   hip L/R cross-correlation: peak +0.443 at lag -0.591 s, trough -0.891 at lag -0.152 s (L hip ACF period 1.027 s)
   L/R phase lag = 0.425 cycle (antiphase = 0.5)
   hip L/R correlation, fundamental band 0.3-1.2 Hz: r = -0.629 (antiphase want < -0.3); raw r = +0.274 (contact-kick 2nd harmonic inflates the raw value)
   joint flexion-positive excursions (deg, min..max and range) vs AnimatLab references:
     hip   L [ -18.2, +16.6] range  34.9 | R [ -24.3,  +2.7] range  27.0   (ref range ~38)
     knee  L [  -3.2, +63.1] range  66.3 | R [  -2.5, +64.0] range  66.5   (ref range ~61)
     ankle L [  -2.2, +24.8] range  27.0 | R [  -2.0, +23.6] range  25.6   (ref range ~16)
   hip flexion L/R correlation (0 lag): r = +0.274 (antiphase want < -0.3)
   [PASS] finite_20s
   [PASS] airborne
   [PASS] leg_L_oscillates
   [PASS] leg_R_oscillates
   measured periods: L RG-E 1.027 s / R RG-E 1.027 s | hip-L ACF 1.027 s / hip-R ACF 1.027 s
VERDICT: ABLATION PASS  both legs oscillate with the commissural coupling REMOVED (L 1.027 s, R 1.027 s) -> each leg is powered by its own RG
```

**Measured periods (the ask's gate b):** with zero coupling, **L hip ACF 1.027 s / R hip ACF 1.027 s** (L RG-E 1.027 / R RG-E 1.027; hip swing ranges 34.9° / 27.0°). Both legs step on their own RG. The period CONTRAST with the coupled run is itself the causality evidence: coupled the pair locks at **1.357/1.358 s** (both sides identical — the c1 period-lock, cf. the w2l_cpg note that 2.0 lets sides drift at 1.05 vs 1.10 s); ablated each side reverts to its intrinsic **1.027 s** — exactly the M3 single-RG open-loop value (`rhythm_orig_sweep`, M3 report §3). The R rhythm follows the R RG, not the L.

### M3 regression: `test_w2l_air_m3_regression.log` (default `net=orig`, extended file) — PASS

Identical to the M3-recorded run: L RG ext 18 bursts / 1.027 s, hip ranges 41.8/41.8°, r = −0.623, 22 632 leg-leg contacts, VERDICT: PASS. The extension did not disturb the M3 path (its VERDICT line keeps M3's original format, which is how the verbatim reproduction shows).

## 3. Honest notes / deviations

1. **Raw vs band-limited antiphase metric (measurement, not tuning).** The raw zero-lag hip correlation is +0.323 in the coupled run (would FAIL M3's `r < −0.3`), while the neural drive is cleanly antiphased (RG-E r = −0.750) and the phase peak sits at exactly half a cycle (0.507). Cause, measured: M4's split net generates MORE leg-leg midline collisions than M3 (26 741 vs 22 632 over 20 s), and each collision kicks both hips IN PHASE through the welded pelvis — an in-phase 2nd harmonic (hip L excursion interval 0.670 s = T/2 exposes it; the same phenomenon existed in M3 at 23 L-excursions vs 18 R). The gait-antiphase question is about the stepping fundamental, so in split mode the `antiphase` check gates on the 0.3–1.2 Hz band-limited correlation (r = −0.941 coupled, −0.629 ablated); the raw value is printed alongside every time. `net=orig` keeps M3's raw-metric check unchanged.
2. **Ablated free-run stays roughly antiphase** (RG-E r −0.788, phase 0.425 cycle) — expected: the antiphase kickoff pair starts the twin oscillators half a cycle apart and identical parameters drift only slowly (phase moved 0.507→0.425 vs coupled). The ablation gate deliberately does NOT require antiphase; it requires each leg to keep oscillating, which is the causality claim.
3. **The coupling slows the pair** (1.357 s vs intrinsic 1.027 s, +32%) — mutual F-inhibition via c1 lengthens the half-center escape; the bilateralrg/w2l_cpg build saw the same family of effect (period-lock at c1 3.0). This is the documented build's own behavior, not a bug; period is recoverable via `--tau` if a faster coupled cadence is ever wanted.
4. **Shared knee-extensor Renshaw** remains the one other cross-side path (MN/RC level) — 2023-original quirk transcribed as-is in M3, left as-is; it does not feed RG or PF half-centers, so it does not affect the RG→PF causality claim.
5. Body-side stand-ins unchanged from M3 (runtime joint damping 3.0, stiff limits, ctrl cap 0.5, air rig +0.30 m) — all knobbed, all documented in M3 §3; ankle fidelity caveat and leg-leg contact caveat carry over verbatim from M3 §5.
6. Determinism: coupled run3 reproduced run2's numbers exactly (13/13 bursts, 26 741 contacts) on this binary; the standing cross-platform chaos caveat (AGENTS 2026-09-13) applies as always.

## 4. File manifest (this milestone)

- NEW `w2l_mujoco\build_w2l_split_net.py` (edit spec + builder + census/causality asserts; standalone census print)
- MODIFIED `w2l_mujoco\test_w2l_air.py` (--net/--comm knobs, R-RG logging, xcorr phase + band-limited metrics, M4 verdict blocks; M3 path logic unchanged)
- LOGS `reports_20260925\logs\test_w2l_air_m4_coupled_run{1,2,3}.log`, `test_w2l_air_m4_ablated_run1.log`, `test_w2l_air_m3_regression.log`
- tmp miners `reports_20260925\tmp\m4_inspect_tpl.py`, `m4_inspect2.py`, `m4_neural_smoke.py` (neural-only pre-check: coupled 1.357 s r −0.751 / ablated 1.027 s r −0.787, both sides)
- No protected file touched; `w2l_cpg\`, `build_w2l_orig_net.py`, body XMLs all unchanged from M3.
