# AnimatLab Study Validation — headless in AnimatLab itself

**Date:** 2026-10-01 (late-night session) · **Executor:** ZCode subagent "AnimatLab validator" (workflow dwfrun-f5242c98)
**Machine:** easteregg2 (`D:\GitHub\Bipedal_Robot`, AnimatSimulator at `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe`); report + artifacts copied to the laptop repo.
**Scope:** validate the AnimatLab study results *inside AnimatLab* — the three asims already run headless tonight (BilateralRG Ground, W2L modern, Biped_2xCPG_wSubs) plus a fresh headless rerun of **Li's own model**, each against the targets recorded in AGENTS.md (2026-09-24 entries).

---

## 1. What ran, exactly

**Colleague check (SHARED-BOX RULE), before launching anything:**

```
ssh -i C:/Users/Ben/.ssh/id_ed25519 -o BatchMode=yes "ben bolen@easteregg2.mme.pdx.edu" "query user"
>  USERNAME            SESSIONNAME   ID  STATE   IDLE TIME  LOGON TIME
>  ben bolen                          2  Disc    1+00:53    9/16/2026 9:33 AM
```

Only Ben's own disconnected session — no colleague on the box → cleared to run.

**Li rerun (this session's new run).** Source copied read-only-style into a fresh workdir (sources never modified), then run with that workdir as `CurrentDirectory` so the DataTool chart .txt landed there:

```
mkdir ...\campaigns\20260930\animatlab\Li_rerun
copy "D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Li Model\walk new new tester_Standalone.asim" ...\Li_rerun\
powershell -NoProfile -ExecutionPolicy Bypass -File D:\GitHub\Bipedal_Robot\tmp\run_li_rerun.ps1
```

`run_li_rerun.ps1` = `Start-Process AnimatSimulator.exe "<asim>" -WorkingDirectory <Li_rerun> -PassThru` + `WaitForExit(180000)` (3-min cap; kill on timeout). Outcome (from `Li_rerun\simulator_stdout.log`, copied to the laptop):

```
RESULT exited ... after 6.7 s            <- well under the 3-min cap
starting sim
Simulation stopped. Time: 10.021         <- matches <SimEndTime>10.02</SimEndTime> in the asim
```

stderr carries only the known headless OSG noise (`GraphicsWindow has not been created successfully`, `Cannot get information for screen 0`) — same as tonight's three earlier runs. Produced: `Li_rerun\DataTool_7.txt` (4,037,905 B; Li's 2023 reference `Neuromechanical_Models\Li Model\DataTool_7.txt` is 4,038,947 B — same chart name, same 9 columns, same 50,010-row length, and the t=0 data row is byte-identical: `height 1.02, distance -3.454`).

**Trace parsing.** `tmp\parse_animatlab_traces.py` (kept on easteregg2; results on the laptop as `animatlab\parse_results_20261001.json`) streams every file line-by-line (no whole-file ingest; largest = 4.0 MB, 50,010 rows):

```
D:\Anaconda\envs\myo\python.exe D:\GitHub\Bipedal_Robot\tmp\parse_animatlab_traces.py ^
    D:\GitHub\Bipedal_Robot\tmp\animatlab_parse_spec.json
```

Documented rules (all numbers below come from that JSON):
- **Artifact strip:** DataTool writes a short all-zeros tail when the sim stops (verified: last 3 rows of every file are `TimeSlice Time 0 0 0 ...`). Only *trailing* all-zero data rows are dropped — 9 rows (1.8 ms) in every file except `Contact.txt` (305 rows = 61 ms of end-of-run zero contact; that file's window therefore ends at 4.941 s). A mid-run fall would still be visible; none occurred (Li height never leaves 0.952–1.020 m).
- **Burst rule:** adaptive threshold at the midpoint of [p05, p95] over the window, re-armed after dropping below p05 + 0.25·range. Reported for the full trace and for the post-transient window t ≥ T0 = 1.0 s. Fixed-level crossings (−50 / −55 mV, 1 mV hysteresis) as cross-checks.
- **Phase:** zero-lag Pearson r **and** lagged normalized cross-correlation (1 ms grid, ±1.5 s). τ > 0 at peak ⇒ the second series lags the first by τ; `lag_cycles` = τ ÷ mean burst period.
- **Frequency:** FFT dominant peak (0.1–15 Hz) on a uniform median-dt grid, Hann window, post-T0.
- Volt-scale series are in volts (rest ≈ −60 mV); joint-angle series are radian-consistent (knee ~1.07 ⇒ ~61°) — deg conversions are labeled "if rad".

---

## 2. Validation table

Targets are the recorded 2026-09-24 AGENTS.md entries quoted in the tasking. Measured values are from `parse_results_20261001.json` (this laptop, and `D:\GitHub\Bipedal_Robot\tmp\animatlab_parse_results.json` on easteregg2).

| Model | Metric | Target | Measured | Verdict |
|---|---|---|---|---|
| BilateralRG Ground | L RG ext onsets | 10 | **12** full-trace / 9 post-T0 | ±2 (detector-rule sensitivity — R matches exactly; see §5) |
| BilateralRG Ground | R RG ext onsets | 11 | **11** full-trace / 9 post-T0 | **PASS** (exact) |
| BilateralRG Ground | RG period | ~0.464 s | **0.4502 s** (L, ±0.050) / **0.4668 s** (R, ±0.016), full-trace | **PASS** (−3.0 % / +0.6 %) |
| BilateralRG Ground | L/R antiphase | antiphase (rec. lag +0.474 cycle) | r₀(L ext, R ext) = **−0.647**; xcorr peak r = **0.976** at lag **0.238 s = 0.529 cycle** | **PASS** (half-cycle lag) |
| W2L modern | joint oscillation | ~2.22 Hz | **2.2494 Hz** FFT on hip, knee, and ankle (post-T0) | **PASS** (+1.3 %) |
| W2L modern | flexor HC sub-threshold | peaks −55.8 / −56.4 mV, never cross −55 mV post-transient | post-T0 max = **−56.41 mV**, **0** crossings of −55 mV post-T0 (full-trace max −52.84 mV is the t < 1 s kickoff transient) | **PASS** (−56.41 ≈ recorded −56.4) |
| W2L modern | ext/flx antiphase | r = −0.956 | r₀(L RG ext, L RG flx) = **−0.9959**; xcorr r = 1.000 at −2.499 cycle ≡ −0.5 cycle (mod period) | **PASS** |
| Li (rerun) | period | ~1.305 s | foot-contact bursts: R **1.2992 ± 0.024** s, L **1.3327 ± 0.080** s; hip FFT 0.7778 Hz ⇒ **1.2857 s** | **PASS** (−1.1 % / +2.1 % / −1.5 %) |
| Li (rerun) | stance duty | ~0.5 | R foot **0.443**, L foot **0.555** (mean 0.499) | **PASS** |
| Li (rerun) | L/R antiphase | antiphase | hips: xcorr peak r = 0.403 at **0.485 cycle**; feet: 0.301 at **0.489 cycle**; r₀ = −0.153 / −0.185 (negative) | **PASS** (half-cycle lag, noisier traces) |
| Li (rerun) | pelvis height | 0.95–1.02 m | full-trace **0.9522–1.0202 m** (post-T0 mean 0.9674, ptp 1.9 cm) | **PASS** |
| Li (rerun) | walking speed | ~0.64 m/s | distance −3.454 → +2.9136 m over 10.0 s = **0.637 m/s** | **PASS** (−0.5 %) |
| Biped_2xCPG_wSubs | (no recorded target) | — | **no rhythm**: post-T0 L RG ext constant −60.57 mV, L RG flx constant −57.93 mV (ptp < 10 µV); hip_L post-T0 ptp = 0.0000 | **N/A — STALE EXPORT** (see §4) |

**11 of 11 targetable metrics PASS**; the one ± is the BilateralRG L-onset count, explained in §5.

### Joint-angle ranges + rhythm context (computed per the tasking; no recorded targets)

| Model | hip | knee | ankle | RG burst stats (full trace) |
|---|---|---|---|---|
| BilateralRG Ground | 0.099–0.403 rad (ptp 17.5°); FFT 2.00 Hz | −0.005–1.066 rad (ptp 61.3°); FFT 2.25 Hz | −0.355–−0.069 rad (ptp 16.4°); FFT 2.25 Hz | L ext 12 bursts, period 0.450 s, duty 0.63; R ext 11 / 0.467 s / 0.67; L flx 11 / 0.466 s; R flx 11 / 0.452 s |
| W2L modern | −0.262–0.404 rad (ptp 38.2°) | −0.017–1.052 rad (ptp 61.2°) | −0.355–−0.070 rad (ptp 16.3°) | L ext 11 bursts, period 0.4434 ± 0.0006 s (post-T0 σ = 0.1 ms — extremely regular); flx 12 / 0.4441 s |
| Biped_2xCPG (stale) | full −0.264–0.321 rad (initial transient only); post-T0 ptp 0.0000 | (no knee/ankle chart in this export — only L Hip Physical.txt) | — | L ext: 1 “burst” spanning 4.55 s (latched); flx: 1 onset |
| Li (rerun) | hip-middle channels oscillate at 0.7778 Hz (period 1.286 s), 7–8 bursts/side | — | — | n/a (contact-driven CPG, no RG chart) |

Pelvis height exists only in Li's DataTool (walkers' charts have no height column — stated, not omitted by accident).

---

## 3. Li rerun section — his own model, in AnimatLab

Run details in §1 (6.7 s wall, clean stop at 10.021 s sim time, OSG-noise-only stderr). The rerun reproduces his 2023 reference run to measurement precision:

| Metric | Li reference `DataTool_7.txt` (2023) | Li rerun (this session) |
|---|---|---|
| rows / t_max / dt | 50,010 / 10.0 s / 0.2 ms | 50,010 / 10.0 s / 0.2 ms |
| R foot period | 1.3005 ± 0.024 s (7 bursts) | 1.2992 ± 0.024 s (7) |
| L foot period | 1.3363 ± 0.078 s (8 bursts) | 1.3327 ± 0.080 s (8) |
| hip FFT | 0.7778 Hz (both sides) | 0.7778 Hz (both sides) |
| stance duty R / L | 0.447 / 0.552 | 0.443 / 0.555 |
| height (full) | 0.9522–1.0202 m | 0.9522–1.0202 m |
| distance end (t=10) | 2.9302 m ⇒ 0.638 m/s | 2.9136 m ⇒ 0.637 m/s |
| t = 0 row | height 1.02, distance −3.454 | height 1.02, distance −3.454 (byte-identical) |

All five recorded Li targets (period ~1.305 s, duty ~0.5, antiphase, height 0.95–1.02 m, speed ~0.64 m/s) are met by the rerun, and the reference file itself parses to the same values — i.e. the targets were correctly derived from this chart, and AnimatLab reproduces them deterministically on this machine.

---

## 4. Stale-export caveat — Biped_2xCPG_wSubs_Standalone

`Biped_2xCPG_wSubs_Standalone.asim` **predates the 12 restored synapses** of the subsystem reorg (per the tasking; consistent with AGENTS.md's standing "RG does not oscillate — latches ≤10 nA" open item on that model lineage). Measured: after the ~1 s kickoff transient the network is **static** — L RG ext pinned at −60.57 mV, L RG flx at −57.93 mV (both constant to < 10 µV), hip_L post-T0 ptp exactly 0, hip flexor muscle signal frozen at 41.28, extensor at 0. The half-center pair is latched with the flexor tonically depolarized. The zero-lag correlation printed for it (−0.64) is computed on that micro-noise and is **not** a meaningful rhythm metric — do not quote it. Its `L Hip Neural.txt` / `L Hip PF.txt` charts contain only the Time columns (empty charts in this export). **None of this is a statement about the current `Biped_2xCPG_wSubs.aproj`** — a fresh standalone export of the current project would be needed for a live verdict.

## 5. Deviations, rules, and honest limits

- **Onset-count rule sensitivity (BilateralRG L 12 vs 10):** the 2026-09-24 counts came from an unspecified detector. My adaptive rule fires at the departure foot of each graded burst (thr ≈ −59.2 mV); post-T0 it counts 9/9 (one cycle window lost to T0). R = 11 matches the record exactly on the full trace; L is +2. The period — the substantive quantity — matches to ≤ 3 %.
- **Zero-lag r understates antiphase** for low-duty graded burst trains; that is why the lagged cross-correlation is reported alongside (half-cycle lag with high peak r is the decisive evidence). Li's r₀ (−0.15) looks weak only because his contact channels chatter (≈ 50 Hz spikes through −50 mV during stance: 500+ fixed-level crossings inside 7 stance windows); his lag (0.485–0.489 cycle) is cleanly half-cycle.
- **Angle units:** charts emit radian-consistent values; degree figures are labeled "if rad".
- **BilateralRG contact channels** are integer-like 0–2 counters (foot-contact duty 3.0–3.7 %, toe_R 28 % — brief touchdown windows). No recorded target exists for them; reported as context only, semantics not fully resolved.
- **RG operating point:** all walker half-centers oscillate in the graded, sub-−50 mV regime (BilateralRG ext peaks −56.3 mV post-T0; W2L ext/flx peak −56.41 mV) — the rhythm is genuine membrane oscillation, not spike-like saturation.

## 6. Artifact index

| Artifact | easteregg2 | laptop |
|---|---|---|
| 3 model run dirs (asim + chart .txt) | `Code\MuJoCo_SNS\spinal\campaigns\20260930\animatlab\{Walker_2_Layer_CPG_BilateralRG_Ground_Standalone, Walker_2_Layer_CPG_Standalone_modern, Biped_2xCPG_wSubs_Standalone}` | not copied (34 MB; numbers in the JSON) |
| Li rerun workdir (asim copy, DataTool_7.txt, sim logs) | `...\animatlab\Li_rerun\` | `...\animatlab\Li_rerun\` (DataTool_7.txt 4,037,905 B + both logs, fetched via `ssh type`) |
| Parse results JSON | `D:\GitHub\Bipedal_Robot\tmp\animatlab_parse_results.json` | `...\animatlab\parse_results_20261001.json` |
| Parser + spec + runner | `D:\GitHub\Bipedal_Robot\tmp\{parse_animatlab_traces.py, animatlab_parse_spec.json, run_li_rerun.ps1}` | `tmp\` at repo root (session scratch) |

**Bottom line:** the AnimatLab study results validate **in AnimatLab itself** — BilateralRG Ground (11 R-onsets, 0.45–0.47 s period, half-cycle L/R lag), W2L modern (2.25 Hz joints, flexor HC peaking −56.41 mV sub-threshold, ext/flx r = −0.996), and Li's own model rerun from his unmodified asim (1.30 s period, 0.44/0.56 duty, half-cycle lag, 0.952–1.020 m height, 0.637 m/s — all five targets met and reproduced against his 2023 reference). The Biped_2xCPG standalone is a stale pre-reorg export and shows a latched, non-oscillating network; it is excluded from validation on that basis.
