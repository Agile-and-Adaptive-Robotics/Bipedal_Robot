# Cross-platform bipedal campaign — consolidated report (2026-10-01)

**Tracks:** AnimatLab validation (in AnimatLab, easteregg2) · MuJoCo port of the AnimatLab
walkers (easteregg2) · SNS Simscape rebuild + validation (laptop) · BPA actuation
air-stepping + standing (MuJoCo half easteregg2, Simscape half laptop).
**Reporter:** ZCode subagent "Publisher" (this workflow). Every headline claim below was
checked by the Publisher against the four track reports AND the artifacts/logs they cite,
all present on the laptop; two claims were additionally re-run this session (§4). Verdict
lines are quoted verbatim from the logs of record. Campaign context: `INDEX.md` (same
folder) — this report consolidates only the four cross-platform tracks, not the easteregg2
retune/robustness/validation job families the INDEX tracks separately.

> **EB475WS4 dissertation track: NOT RUN.** Nothing in these four tracks executed on
> EB475WS4 — the AnimatLab/MuJoCo work ran on easteregg2 and the Simscape/BPA-demo work on
> the laptop (Simscape lives on the laptop per Ben's 09-30 call; easteregg2 is R2025a and
> cannot load the R2025b models — INDEX.md item 4). No artifact of this campaign exists on
> EB475WS4. **Sync:** the campaign's changed files are UNCOMMITTED on the laptop and
> easteregg2 working trees (INDEX.md "Code changes"; the shared `w2l_mujoco/` files
> hash-identical on both, per `mujoco_port.md` §2d). Commit + push to GitHub `origin` from
> whichever machine has the change, then pull on the others — GitHub origin is the single
> source of truth; there is no machine-to-machine push path in use.

---

## 1. Track summaries — what ran, evidence, verdict

### 1.1 AnimatLab validation (in AnimatLab itself) — `animatlab/VALIDATION.md` — VERDICT: PASS (11/11 targetable metrics)

What ran (easteregg2, after a `query user` colleague check showing only `ben bolen`, Disc):
the three asims already run headless this campaign + a fresh headless rerun of **Li's own
model** from his unmodified asim in a clean workdir, all traces parsed by a documented
adaptive-burst rule (`parse_animatlab_traces.py` → `animatlab/parse_results_20261001.json`,
on the laptop).

Evidence (all re-read by the Publisher from the parse JSON and the rerun logs):

- Li rerun: `simulator_stdout.log` = `starting sim` … `Simulation stopped. Time: 10.021`
  (6.7 s wall); stderr OSG-noise-only; `DataTool_7.txt` 4,037,905 B with a byte-identical
  t=0 row (`height 1.02, distance -3.454`) to his 2023 reference.
- BilateralRG Ground: L/R RG-ext bursts **12/11**, periods **0.4502 / 0.4668 s** (parse
  JSON: `"bursts": 12, "mean_period_s": 0.4502` / `"bursts": 11, "mean_period_s": 0.4668`);
  L/R lag **0.529 cycle** (xcorr r 0.9762; zero-lag r −0.6473).
- W2L modern: joints FFT **2.2494 Hz** on hip/knee/ankle; flexor HC post-T0 max
  **−56.41 mV, 0 crossings of −55 mV**; ext/flx r **−0.9959** (parse JSON `flexcheck` +
  `corr_post`).
- Li rerun vs his 2023 reference file (both parsed, both in the JSON): foot periods
  **1.2992/1.3327 s** (ref 1.3005/1.3363), duty **0.443/0.555** (ref 0.447/0.552), hip FFT
  0.7778 Hz both, height 0.9522–1.0202 m, distance end 2.9136 m (ref 2.9302) ⇒ **0.637 m/s**.
- `Biped_2xCPG_wSubs_Standalone` measured **latched** post-T0 (L RG ext constant −60.57 mV,
  one 4.55 s "burst"; hip ptp exactly 0) → **excluded as a stale pre-reorg export**, not
  counted as a failure.

**Verdict: 11 of 11 targetable metrics PASS**; the single ± is BilateralRG L-onsets 12 vs
recorded 10 (detector-rule sensitivity; the substantive quantity, period, matches ≤ 3 %).

### 1.2 MuJoCo port of the AnimatLab walkers — `mujoco_port.md` — VERDICT: stand PASS / walking NOT YET (partial)

What ran (easteregg2, WMI-detached, colleague check first): Li drop-protocol gate
(default + 6-knob sweep, `test_li_stepping.py`) and the W2L ground gate converted to Ben's
drop protocol (`test_w2l_ground.py --phase=both`), plus a re-run of `test_w2l_air.py`.
Three pre-existing breakages fixed first (stale 166-synapse censuses in the split/aff
builders after the 09-30 V3 correction — now 87/164, independently rebuilt by the
Publisher this session, §4; duplicate root freejoint injected by `write_ground_xml()`;
laptop tree 2 files behind, synced + hash-verified).

Evidence (verbatim from the fetched logs, all in `tmp/li_fetch/`):

- Li, all 7 protocol variants: `L_foot: 0 stance episodes, period nan s … VERDICT: not yet`.
  Default + sweep (support 0.5/0.85/1.0, clear 0.02/0.06, hold 2.0) all hang at pelvis
  **1.001–1.007 m** (e.g. `li_sw_support05.log`: `pelvis height: min 1.001 m, end 1.004 m`);
  the two earlier release-style runs buckled to **0.099 / 0.138 m**. Diagnosis: the CPG
  runtime pose (knees −64…−67° flexion) keeps the feet clear of the platform under the
  retained PD.
- W2L full gate (`w2l_ground_protocol.log`): `M6 VERDICT: stand=FULL; walk=FALLS/NO-GAIT`.
  Stand: `finite=True fall=no … tilt mean/max 2.0/2.1 deg … [PASS-FULL] unrigged_stand` at
  all 4 rig scales, harness ~50 %. Walk: `heel duty L 0.00 R 0.00 | forward disp +0.001 m`,
  `neural: L RG-E bursts 11 (period 1.539 s), R RG-E bursts 10; L/R RG-E r -0.776`,
  `contact encoder evidence: heel SN max 0.00 mV, toe SN max 2.57 mV` — platform reached,
  but **toe-only loading**.
- `test_w2l_air.py` re-run (`w2l_air_rerun.log`): rhythm reproduces (`L RG ext bursts=18,
  period 1.027 s`, `hip flexion L/R correlation (0 lag): r = -0.414`) but
  `ground contacts: 28890 (must be 0) … [FAIL] airborne … VERDICT: FAIL airborne` — the
  09-30 freejoint rebuild broke the air gate's suspended-body assumption (reported, not
  fixed).

**Verdict: the protocol conversion works mechanically (no falls, rhythm survives the drop),
ground walking is blocked at the foot–platform interface** — Li never touches, W2L touches
toes-only.

### 1.3 SNS Simscape rebuild + validation (laptop) — `sns_simscape.md` — VERDICT: PASS

What ran: `matlab -batch … build_sns_w2l_cpg_20260930` (MATLAB R2025b) built
`results/SNS_W2L_CPG.slx` (301,949 B — size confirmed on disk): a clearly-labeled
**representative core** of the W2L 2-layer CPG (4 persistent-Na RG half-centers with FIXED
tau_h 250 ms in a new `NapHC` block, laminated mutual inhibition, c1/V3 crossed
commissurals, PF E/F, 4 MNs, heel-contact drive; 42 synapses, gains verbatim from
`build_w2l_net.py` `W2L_GAINS`). The existing `sns_build_from_json.m` pair targets the
406-neuron gait2392 net and cannot express NaP half-centers — the fresh-core route was
correct.

Evidence (`Code/Matlab/SNS_Simscape/logs/laptop_20260930/sns_w2l.log`, confirmed on disk):

```
SNS_W2L_CPG core: 4 NaP half-centers + 22 graded cells (incl. 2 heel contact encoders) + 42 synapses
W2L_SIMSCAPE PASS | period L 1.000 vs 1.000 s (0.0%) R 1.000 vs 1.000 s (0.0%) | antiphase r -0.525 (ref -0.525) | bursts 10/10
```

Beyond the gate (`w2l_trace_compare.json`, re-read): all four RG traces cross-correlate
r = 0.99996–0.99998 at a −2 ms lag (one solver step — the known h-update-order
difference); RMSE 0.099–0.111 mV over 2–12 s; L RG-E burst onsets 2.088…11.088 s
(identical grid in both engines, ≤ 2 ms apart per onset). Honest findings recorded in the
report: the net **latches under tonic+kickoff alone** (R RG-E −1.154 mV, zero bursts —
re-proven by the Publisher this session, §4) so both engines used the smoke contact-train
protocol; one leak-sign bug (R side → 8.2e61 mV) was caught and fixed before the passing
run. Honest inventory: Renshaw, all afferent chains, the 2nd PF pair per side, toe contact,
and 8 of the 12 MNs are absent from the core; no plant connected.

**Caveat (Publisher-verified):** `w2l_numpy_ref.json`'s hardcoded `protocol` label reads
`"tonic E 2 / F 3 nA + kickoff, no contact"` — a stale string; the script itself DOES
inject the 1 Hz heel trains (`w2l_numpy_ref_20260930.py:81-85`) and its stdout prints the
true protocol. The science is unaffected; do not cite the JSON label.

### 1.4 BPA actuation: air-stepping + standing — `bpa_actuation.md` — VERDICT: MuJoCo half PASS / Simscape half tuning-blocked

What ran (easteregg2, WMI-detached after colleague check): the split-RG W2L walker with
all 12 Hill muscles replaced by `bpa_muscle.py` BPAs (zero-gain general actuators on the
same tendon routes), `test_bpa_stepping.py --dur=20 --hold=6 --lower=2 --kmax-frac=0.6
--cap=0.6 --count=2`. Simscape half (laptop): `matlab -batch "… sns_run_cpg_demo"`.

Evidence (`campaigns/20260930/bpa_run/bpa_gate.log`, re-read in full this session):

```
AIR knee L: 3 flexion excursions, period 1.542 s, flexion range 63.5 deg [-1.8,+61.8]
AIR knee R: 4 flexion excursions, period 1.543 s, flexion range 63.2 deg [-1.3,+62.0]
RETAINED phase (11.5 s of it): pelvis z min 0.916 m, end 0.916 m
RETAINED vertical load: mean L 19.8 N / R 28.1 N; PEAK L 880 N / R 1372 N; stance episodes L 7 / R 5 (harness carries 70% of 411 N)
BPA pressure: peak 372 kPa (ceiling 620), mean-of-active 318 kPa
VERDICT: BPA ACTUATION PASS (air-stepping + retained-support stand on BPA muscles)
```

Key engineering findings (report-quoted): `festo4()` explodes for negative strain → BPA
rest length sized to the MAX sampled route length over a 7³ pose grid (naive sizing hit
13,096 N knee forces, sim dead at t=0.136 s); force-servo lowering replaced the open-loop
z_target lowering (measured 0 N foot load over 11.5 s — the march pose keeps the feet
clear); stiff joint limits required (ankle otherwise blows through its range to +115°,
ankle_flx force 159,975 N). Documented deviations: `kmax_frac=0.6` exceeds a physical
Festo DMSP-20 KMAX ≈ 0.2 (routes were designed for Hill muscles — the real fix belongs to
the Xi/mesh-optimization program); quiet latched standing fails because in this ported
body the **"knee extensor" route produces flexion-sign torque**
(`actuator_moment[knee_L_ext, knee]` = −0.0553 N·m/N at 0°, −0.0427 at −60°; flexor
+0.0171…+0.0957 — report-quoted probe, not re-run by the Publisher) — the passing mode is
the harness-retained MARCH (7/5 stance episodes, 880/1372 N peaks). The gate run used the
same census-fixed builder (§1.2; Publisher re-verified this session).

Simscape half (`tmp/bpa_track/bpacpg_run.log`, re-read): `CPG: 2 hysteresis switches in
10 s (period ~ 0.155 s, freq ~ 6.46 Hz)` / `CPG run OK: theta range 9.7..48.0 deg` — exit
0, RG half-centers anticorrelated r = −0.701, BPA forces ext 0–66 N / flex 0–30 N, but
**0 sustained knee cycles** at stock params: a tuning problem (drive/adaptation/epsScale),
not a missing-component problem. Library inventory: physical `BPA_10/20/40mm` blocks
(P[kPa], L[m] → F[N]) + normalized `BPAForce` exist; missing for a full Simscape BPA
walker are hip/ankle routing, the full PF→12-MN wiring, and ground afferents.

---

## 2. Cross-engine comparison table

Same quantity, three engines. **Operating-point caveat:** these are different net variants
(2023 original vs modern vs split-RG vs representative core) at different drives, so the
period spread reflects protocol, not fidelity. The two apples-to-apples chains are
(a) AnimatLab rerun vs Li's own 2023 reference file, and (b) the Simulink core vs its
numpy full-net reference.

| Quantity | AnimatLab (validated in situ) | MuJoCo port (easteregg2) | SNS Simscape (laptop) |
|---|---|---|---|
| W2L-family RG burst period | BilateralRG 0.4502 / 0.4668 s endogenous; W2L modern 0.4434 ± 0.0006 s | W2L-2023 air 1.027 s; split-RG under drop protocol 1.539 s; **BPA-driven 1.542 s** | numpy full net 1.000 s = the 1 Hz heel-train drive (net latches without contact — §4); Simulink core 1.000 s (0.0 % vs numpy) |
| L/R (or E/F) antiphase | BilateralRG L/R lag 0.529 cycle (xcorr r 0.976); W2L ext/flx r −0.9959 | ground L/R RG-E r −0.776 (11/10 bursts); air hip L/R r −0.414 | numpy r −0.525; Simulink core r −0.525; traces r ≈ 1.000 at −2 ms |
| Li walker outcome | rerun walks: periods 1.2992/1.3327 s, duty 0.443/0.555, 0.637 m/s, height 0.952–1.020 m = his 2023 file | **0 stance episodes in all 7 protocol variants** (hangs 1.001–1.007 m; buckles to 0.099/0.138 m without retained support) | Li net not rebuilt this session |
| Stance / ground contact | Li: real gait (duty 0.44/0.56); BilateralRG foot-contact duty 3.0–3.7 %, toe_R 28 % (brief touchdowns) | W2L: stand=FULL (tilt ≤ 2.1°, harness ~50 %) but heel duty 0.00, toe-only, +0.001 m; **BPA march: 7/5 stance episodes, peak foot loads 880/1372 N** | no plant (neural core only); BPACPGLegDemo 0 sustained knee cycles (forces 0–66 N) |
| BPA actuation | n/a (AnimatLab muscles) | **PASS**: 63° knee swing in air + retained-support stand, 372 kPa ≤ 620 kPa ceiling | demo runs clean (r −0.701) but 0 sustained cycles at stock params |

**Reading:** the neural rhythm itself ports cleanly across all three engines (every
antiphase measure is negative and roughly half-cycle; the Simulink↔numpy match is
essentially exact). What does NOT yet cross engines is closed-loop ground locomotion:
AnimatLab's Li walks at 0.64 m/s, MuJoCo's ports produce no gait cycles (heel-starved W2L,
never-landing Li), and the BPA march is the closest any engine gets to load-bearing BPA
stepping.

---

## 3. Next actions, ranked by what unblocks BPA-driven walking

1. **Get heel contact in MuJoCo (both walkers).** The single gate between "stands/marches"
   and "walks": W2L reaches the platform but loads toes-only (`heel SN max 0.00 mV`,
   heel duty 0.00), Li never touches (hangs 1.001–1.007 m). Levers: contact-gated lowering
   or a PD z-target below keyframe z by the knee-flexion margin (Li; `mujoco_port.md` §3
   next-levers); larger `--sink`, the `acap` ankle-PF trim, and a check of which rebuilt
   foot geom actually bears (W2L; §5.4); extend the BPA run's force-servo lowering (which
   already found 123 N of foot load) to seek heel-specific load. The BPA march already has
   7/5 stance episodes — heel loading is the missing ingredient for contact-driven cycles.
2. **Body-fidelity fix (both walkers; prerequisite for harness-free stance):** stiff
   solimp on the 8 hinges + Kpe as tendon stiffness (goal-3 §3.2 recipe: ≥ 1 kg virtual
   tendon mass or 0.5 ms timestep). Evidence: Li release-style runs buckle to 0.099/0.138 m;
   probe_hold showed extensors-only hold collapses through knee/ankle limits 0.5 s after
   release.
3. **Knee-extensor route sign (BPA body):** the "extensor" route produces flexion-sign
   torque (−0.0553 N·m/N at 0°) — any extensor-latch quiet stand is architecturally wrong
   on this body; correct the route or tune balanced co-contraction.
4. **Simscape BPA sustained oscillation tuning** (drive/adaptation/epsScale on
   BPACPGLegDemo) — the only blocker for the Simscape half of the stretch goal; all
   required block types exist.
5. **Repair `test_w2l_air.py`'s airborne gate** (pin/lift the root for air runs) so the
   air rhythm (18 bursts, 1.027 s, r −0.414) is citable as a PASS again.
6. **Route redesign toward physical Festo KMAX ≈ 0.2** (Xi/mesh-optimization program) —
   `kmax_frac=0.6` is a valid simulation demonstration, not a buildable muscle design for
   these routes.
7. **Commit + push the campaign** (changed files uncommitted on both machines; EB475WS4
   then pulls from GitHub origin — see header).

---

## 4. Evidence-verification appendix (Publisher, this session)

Everything below was read or executed by the Publisher on the laptop during this ask:

- **Track reports + context read in full:** `animatlab/VALIDATION.md`, `mujoco_port.md`,
  `sns_simscape.md`, `bpa_actuation.md`, `INDEX.md` (all under
  `Code/MuJoCo_SNS/spinal/campaigns/20260930/`).
- **Logs read in full, verdict lines match byte-for-byte:**
  `bpa_run/bpa_gate.log`; `Code/Matlab/SNS_Simscape/logs/laptop_20260930/sns_w2l.log`;
  `tmp/li_fetch/li_drop.log` (3 blocks: 0.099/0.138/1.003 m) and all six sweep logs
  `li_sw_{support05,support085,support100,clear002,clear006,hold20}.log` (pelvis-height
  lines grepped: 1.001–1.007 m, every one); `tmp/li_fetch/w2l_ground_protocol.log`;
  `tmp/li_fetch/w2l_air_rerun.log`; `tmp/bpa_track/bpacpg_run.log`;
  `animatlab/Li_rerun/simulator_stdout.log`.
- **JSON artifacts re-read:** `animatlab/parse_results_20261001.json` (BilateralRG 12/11
  bursts @ 0.4502/0.4668 s, lag 0.529 cycle; W2L −56.41 mV / 0 crossings / r −0.9959 /
  2.2494 Hz; Li rerun 1.2992/1.3327 s, duty 0.4432/0.5551, distance end 2.9136 m vs
  reference 1.3005/1.3363 s, 0.4474/0.5515, 2.9302 m; Biped_2xCPG latched — one 4.55 s
  "burst", hip post ptp 0.0); `w2l_numpy_ref.json` (1.000 s / r −0.525 / 10+10 bursts /
  stale protocol label); `w2l_trace_compare.json` (RMSE 0.099–0.111 mV, xcorr r
  0.99996–0.99998 @ −2 ms, onsets 2.088–11.088 s).
- **Deliverables confirmed on disk:** `SNS_W2L_CPG.slx` (301,949 B),
  `SNS_W2L_CPG_run_20260930.mat`, `dev/build_sns_w2l_cpg_20260930.m`,
  `w2l_mjcf_bpa.xml`, `make_bpa_walker.py`, `test_bpa_stepping.py`,
  `sns_w2l_traces.png`, `Li_rerun/DataTool_7.txt` (4,037,905 B).
- **Publisher re-runs** (`C:\Users\Ben\.anaconda3\envs\myoconv\python.exe
  tmp\reporter_probe_20261001.py`, laptop, this session — output verbatim):
  - `PROBE_A SPLIT_CENSUS neurons=87 synapses=164` — the census fix both MuJoCo-side
    tracks depend on, independently rebuilt.
  - `PROBE_B TONIC_ONLY window 2-12s L RG-E max 4.737 mV bursts 0 | R RG-E max -1.154 mV
    bursts 0` → `LATCH (no rhythm)` — reproduces the Simscape track's latch finding at the
    exact quoted value (−1.154 mV), confirming the contact-train protocol was necessary.
- **Source check:** `w2l_numpy_ref_20260930.py:81-85` injects the alternating 1 Hz heel
  trains while line 115 writes the stale "no contact" label — basis of the §1.3 caveat.

**Not independently re-run by the Publisher:** the easteregg2 simulations themselves
(verified via their fetched logs only), the MATLAB builds/runs (verified via
`sns_w2l.log` / `bpacpg_run.log` and the on-disk `.slx`/`.mat`), the colleague-check
`query user` calls (quoted in the track reports), and the BPA moment-arm / sizing-incident
probes (quoted parsed numbers in `bpa_actuation.md`; no preserved logs).
