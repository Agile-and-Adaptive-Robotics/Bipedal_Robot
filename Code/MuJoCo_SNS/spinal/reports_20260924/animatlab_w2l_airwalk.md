# AnimatLab Walker_2_Layer_CPG headless air-walk verification — 2026-09-24 (EB475WS4)

Task: run Ben's Walker_2_Layer_CPG AnimatLab model headless and verify it walks in air
(Ben: "run the Walker_2_Layer_CPG in animatlab. You will see that it works at walking in
the air. We just didn't tune the feedback.").

**Verdict up front: CONFIRMED for three of the four standalone files tried — but the
details matter.** The 2023 `walk new new tester added 2 axis_Standalone.asim` air-walks
with clear spiking bursts (~0.77 Hz, L/R antiphase). The modern
`Walker_2_Layer_CPG_Standalone_modern.asim` also produces sustained, antiphase
air-stepping (2.22 Hz, hip 38° / knee 61° / ankle 16° swings), but its flexor half-centers
never cross the −55 mV burst threshold — the rhythm there is a ~5 mV sub-threshold
alternation. The BilateralRG modern variant's RG does cross −55 mV (2.16 Hz, antiphase).
Full numbers per run below.

## Environment

- AnimatLab found via registry
  `HKLM\SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\AnimatLab_is1`
  → `InstallLocation = D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\`.
  Runner: `<install>\bin\AnimatSimulator.exe` (51,200 bytes, 2014). It takes NO flags —
  `AnimatSimulator.exe --help` just prints `No simulation file '--help' was found.`
  So the documented bare invocation `AnimatSimulator.exe "<path>\file.asim"` is the only one.
- Each run was launched from the bin dir with stdout+stderr tee'd to a log; every run
  ended cleanly (`starting sim … Simulation stopped. Time: 5.1` / `10.021`). The only log
  noise is two benign OSG `grabFocusIfPointerInWindow() … Access is denied` warnings
  (headless session, no interactive desktop cursor) — no sim errors.
- AnimatLab2.exe (Ben's GUI, PID 23036) was up during the runs and was left untouched;
  sconestudio.exe untouched; **no model file was modified** (the runs only rewrote chart
  `.txt` byproducts beside the asims and created `DataTool_7.txt`).

## What ran, and what the charts show

### 1. `Walker_2_Layer_CPG_Standalone_modern.asim` (the named primary target)

Command: `AnimatSimulator.exe "D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG_Standalone_modern.asim"`
Sim end 5.1 s (`<SimEndTime>5.1`), charts record 0.01–5 s at 0.2 ms.

Chart updates (BEFORE → AFTER mtimes): all 14 chart `.txt` files beside the asim were
rewritten — e.g. `Rhythm Generator.txt` 2026-09-16T16:49:58 → 2026-09-24T04:44:00.607,
`L Hip PF.txt` → 04:44:00.791, … `L Angles.txt` → 04:44:03.426 (full list in
`_w2l_mtimes_before.json` / `_w2l_mtimes_after.json`). All 14 landed **byte-size-identical**
to the 09-16 outputs (deterministic rerun; the BilateralRG rerun below proved exact
reproducibility). Trace_2020* files untouched.

**Flexor bursts at −55 mV (the 09-18 criterion): ZERO** in `Rhythm Generator.txt`
(L RG flx) and `L Hip PF.txt` / `L Knee PF.txt` / `R Hip PF.txt` / `R Knee PF.txt`
flexor columns — each shows exactly 1 downward crossing, which is the initialization
transient (half-centers start at −52.8 mV, relax to rest) and no crossings after.

**But the rhythm is there — it is sub-threshold.** Unit note: the asim's neurons carry
`RestingPot=-60` (mV) and the charts store **volts**, so −55 mV = −0.055 in chart units
(first analysis pass used the wrong unit and read flat-zero everywhere). Post-transient
(t ≥ 0.5 s), every charted neuron oscillates coherently at **2.222 Hz** with ~5–6 mV
swing that peaks just *below* threshold (L RG: peaks −56.4 mV; PFs: peaks −55.8 mV;
FFT peak/rms ≈ 6 — a spectral line, not noise):

- ext/flx alternation within a pool: L RG ext vs L RG flx **r = −0.956**; L Hip PF ext vs
  flx **r = −0.977** (the four chart columns are distinct neuron GUIDs — this is real
  half-center alternation, not a plotting artifact).
- L vs R: L Hip PF flx vs R Hip PF flx **r = −0.994 at 0 lag** (mirror image) →
  **antiphase**; mid-swing (−58.5 mV) crossing lag **−0.492 cycle**, both sides
  0.4432 s median interval, 10 crossings each in 4.5 s.
- Joints actually step (`L Angles.txt`): hip −15.0..+23.2° (38.2° swing), knee
  −0.9..+60.3° (61.2°), ankle −20.4..−4.0° (16.4°), all at 2.222 Hz with 10 half-range
  upcrossings at 0.443–0.445 s median intervals.
- Neural→joint coupling: r = 0.90–0.92, neural leads joint by ~57–89 ms (knee PF flx
  leads knee angle +57 ms; RG ext leads hip +63 ms in the negative-lag sense reported
  by the sweep) — the leg swing is phase-locked to the neural rhythm.

Why so weak: this export has **no tonic drive at all** (every neuron `TonicStimulus=0`)
and its only external stimulus is `Stimulus_1` — a 10 ms, `CurrentOn=1e-8` tonic-current
tick at t=0 on `L RG ext` (asim lines 9626–9644). The 2023 model it derives from runs on
14 neurons held at 5–6 nA tonic (see run 2). The modern network free-runs on whatever the
tiny kick leaves behind, hence millivolt-scale amplitudes that sit under the −55 mV
burst line while still driving real 38–61° leg swings.

### 2. `walk new new tester added 2 axis_Standalone.asim` (2023, same folder)

Command: `AnimatSimulator.exe "D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\walk new new tester added 2 axis_Standalone.asim"`
Ran to 10.021 s; its single chart wrote `DataTool_7.txt` (4,038,947 bytes, created
2026-09-24 04:53; columns: R/L foot & toe ground-contact neurons, R/L "hip middle"
neurons, height, distance).

**This is the cleanest "works in air" demonstration:**

- All six logged spiking neurons burst **8 times** over 10 s; period **1.287–1.305 s
  (~0.77 Hz)**, burst durations 456–667 ms (e.g. R_hip middle: period median 1.302 s,
  burst dur median 509 ms; L_foot ground contact: 1.305 s, 667 ms).
- L/R hip middle **onset lag +0.457 cycle (+596 ms) ≈ antiphase**.
- Body stays suspended: height 1.020 → settles 0.958–0.970 m for the whole 10 s (the
  `min=0` in a naive read is 9 trailing zero-fill rows at t ≥ 10.000 s). It never falls.
- It drifts in x from −3.454 to +2.930 m (+6.4 m) — the asim's kickoff is a
  `ForceInput` 40 N push for 0–1 s on the root body (asim lines 7270–7288); with no
  footing in air the body just glides while the legs step.
- Drive source: 14 neurons with `TonicStimulus` 5–6 nA (the endogenous drive the modern
  export lost), plus the 40 N kickoff.

### 3. `Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim`

Command: `AnimatSimulator.exe "…\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim"`
Ran twice (second run to recover charts overwritten mid-analysis by the Ground run below);
**both runs byte-identical** (same RG onsets to the millisecond) — deterministic.

- The RG flexors DO cross −55 mV: **L RG flx 10 onsets** [0.710 … 4.878 s], period
  **0.4640 ± 0.0025 s (2.16 Hz)**; **R RG flx 11 onsets**, period 0.4641 s.
  Duty only ~5% (brief peaks to −52.4 mV out of a −64.8..−52 mV swing).
- L→R onset lag **+0.474 cycle (+220 ms) ≈ antiphase** (RG chart here logs both sides,
  plus L/R c1 and L/R V3 commissural cells).
- PF flexor half-centers: 0 crossings (same sub-threshold family as run 1).
- Joints step repeatedly at the RG frequency (`L Angles.txt`): hip −6.9..+23.3°,
  knee −0.3..+37.2° (brief flexion bursts), ankle −20.3..−4.1°; half-range upcrossing
  intervals 0.461–0.463 s (hip/knee 10 each, ankle 9).

### 4. `Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim`

Command: `AnimatSimulator.exe "…\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim"`

- Same RG rhythm as run 3 (L 0.4643 s / R 0.4580 s, ~2.15–2.18 Hz, duty ~5%);
  R Hip PF flx shows 3 marginal late onsets (2.264–3.673 s) — likely contact-driven.
- Joints step at 0.467–0.470 s median intervals; knee reaches 61.1° (vs 37.2° in air).
- `Contact.txt`: toe_L mean 0.26 / max 2, toe_R mean 0.49 / max 2, foot_L mean 0.08,
  foot_R mean 0.02 — **intermittent real foot-ground contact**, i.e. this variant steps
  with occasional ground touches, as its name implies.

## Honest caveats

- **Feedback is untuned, per Ben's own words** ("We just didn't tune the feedback"), and
  **there is no balance layer** — nothing here is claimed to walk on the ground
  under control; the air variants hang at ~0.97–1.02 m and step.
- The −55 mV flexor-crossing criterion (09-18 notes) is **silent on the modern W2L
  export**: its rhythm exists but peaks ~0.5–1.4 mV below the threshold. If that criterion
  is reused in scripts, store it in chart units (−0.055; charts are volts, rest −60 mV),
  and be aware an exports-without-tonic model can step legitimately at 0 crossings.
- The modern W2L/BilateralRG exports carry only a token kickoff (10 ms, 1e-8 current tick;
  all-neuron tonic 0) vs the 2023 model's 5–6 nA tonic on 14 neurons + 40 N push. The
  modern RG's ~2.2 Hz millivolt alternation is a different (weaker) regime than the 2023
  model's full spiking bursts at 0.77 Hz.
- Chart columns for RG are LEFT-side-only in `Rhythm Generator.txt` of the
  Walker_2_Layer_CPG modern export (no R RG columns); L/R antiphase for that model was
  measured from the R/L Hip PF charts instead. The BilateralRG RG chart logs both sides.
- File mtimes of the chart `.txt` byproducts beside each asim now reflect these runs
  (2026-09-24 04:44–04:56). Copies of the outputs as-produced are preserved at
  `w2l_modern_charts_20260924T0444\` (14 files) and `bilat_modern_charts_20260924\`
  (15 files, incl. Contact.txt) under this report folder. No `.asim`/`.aproj`/`.aform`
  was modified.

## Artifacts (this folder)

- `animatlab_w2l_airwalk.md` — this report.
- `_w2l_headless_run.log`, `_w2l2023_headless_run.log`, `_bilat_headless_run.log` (+`2`),
  `_bilatground_headless_run.log` — raw runner stdout/stderr.
- `_w2l_mtimes_before.json` / `_w2l_mtimes_after.json` — chart mtime/size snapshot diff.
- `_w2l_mtimes.py`, `_w2l_mtimes2.py` — mtime snapshot helpers.
- `_w2l_analyze_charts.py` — flexor −55 mV burst/antiphase analysis (folder-parameterized).
- `_w2l2023_analyze.py`, `_w2l2023_analyze2.py` — DataTool_7 spike/burst, antiphase,
  body-motion analysis.
- `_bilat_extras.py` — L/R RG onset antiphase + L Angles ranges.
- `_w2l_freq_deepdive.py` — post-transient amplitudes, FFT frequencies, neural↔joint
  correlation (modern W2L).
- `_w2l_twin_check.py` — ext/flx alternation correlations (modern W2L).
- `_w2l_lr_phase.py` — L/R hip PF antiphase (modern W2L).
- `_angles_course.py` — joint time-course / repeated-stepping check.
- `_w2l_chart_targets.py` — L Hip PF chart DataColumn target dump.
- `w2l_modern_charts_20260924T0444\`, `bilat_modern_charts_20260924\` — preserved chart
  outputs.
