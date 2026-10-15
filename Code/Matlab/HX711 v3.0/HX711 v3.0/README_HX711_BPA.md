# HX711_BPA — BPA Force + Pressure Test App (AARL)

Rebuild of the customized Matlab–Arduino–HX711 app (the one with the
encoder-angle metadata and the pressure sensor). One app, two calibration
routes, portable across all machines and GitHub clones.

## Quick start (any machine, any clone location)

```matlab
run("path\to\Bipedal_Robot\Code\Matlab\HX711 v3.0\HX711 v3.0\Start_HX711_BPA.m")
% or, after cd'ing anywhere in MATLAB:
Start_HX711_BPA
```

The launcher adds **its own folder** to the MATLAB path, so the Arduino
add-on package `+arduinoioaddons/+basicHX711` resolves no matter where the
repo was cloned. No absolute paths are used anywhere: the default save
folder is found by walking up from the app file until a `Testing_Data`
folder is seen, then defaulting to `Testing_Data\2026_06_Festo\Flx_20mm`
(falls back to `<app folder>\saved_data` and can always be overridden in
the Save Data tab).

**The app is named `HX711_BPA` on purpose.** An old File-Exchange "HX711"
app (or an Add-Ons install of it) can exist on other machines and would
shadow a generic `HX711` name; a unique name makes that impossible. If you
type `HX711` and something other than this app opens, that is an installed
add-on copy — uninstall it in MATLAB (APPS > Manage Apps) or just always
use `Start_HX711_BPA`.

Prerequisites per machine: MATLAB R2025a or newer, the *MATLAB Support
Package for Arduino Hardware*, and the rig's Arduino. The HX711 library is
uploaded from this folder automatically on first **Connect**.

## Calibration — two routes, both always available

**Route 1, original workflow** (needs hardware):
1. **Tare** — mean raw counts at zero load (Settings tab: number of readings).
2. **Scale Factor** — hang the known weight (Settings tab: weight + unit);
   slope = (mean − tare)/weight.
3. **Calibration** — verification overlay: readings vs a Gaussian with the
   known weight marked. **Raw Read** = single converted reading.

**Route 2, known-factor entry** (no hardware step needed):
- *Known LC Cal* tab — type the zero offset/tare (raw counts) and the
  calibration slope (counts per gram-equivalent), click **Apply Known
  Load-Cell Cal**. A field left blank (NaN) keeps its current value, so
  entering just a new zero offset never touches the scale factor.
- *Pressure Cal* tab — type `a` and `b` in
  `Pressure_kPa = a*Voltage_V + b`, click **Apply Known Pressure Cal**; or
  click **Run 7-Point Pressure Cal** for a guided calibration (regulator
  0–620 kPa; enter the actual gauge kPa at each point, the app fits a/b).

Factors persist across sessions in `hx711_bpa_last_cal.mat` next to the app
(machine-local; gitignored, so each machine keeps its own). Every saved
file embeds the factors used, so provenance is never lost.

### Ben's normal load-cell procedure (and how the app protects it)

1. Hang the load cell unloaded → **Tare** (zero offset #1).
2. Take it down, tie the known weight on, hang it → **Scale Factor**.
3. Remove the weight.
4. Mount the load cell horizontally, tie it to the tibia → **Tare** again
   (zero offset #2 for the mounted orientation).

The app is built around this: **Tare updates only the zero offset.** The
scale factor (counts per gram — a property of the cell, not the mounting)
is never modified by a tare, in the live session or in the saved
calibration cache, and the message after each tare states the preserved
scale explicitly. On the Known LC Cal tab the slope field defaults to
blank/NaN; Apply with a blank slope keeps the current scale factor, so a
tare-only entry can never wipe the calibration. The offline self-test
includes a regression check for exactly this re-zero sequence.

Acquisition unlocks as soon as tare + scale are set by *either* route.

## Dynamic pressure calibration (Pressure Ctrl tab)

For a **free BPA** (not mounted on the leg): give a desired pressure and
the built-in PID reaches it through the valves, using time-proportioning.
The PID output `u` in [−100, +100] % becomes, each control period:

| `u` | Valve state | Pins (Increase / Maintain) |
| --- | --- | --- |
| u ≥ +2 % | **FILL** | High / High |
| u ≤ −2 % | **VENT** | Low / Low |
| otherwise | **HOLD** | Low / High |

(the ±2 % minimum duty avoids valve chatter; the "Increase/Maintain"
buttons in the Valves panel drive the same three states manually).

- **Kp [%/kPa], Ki [%/kPa/s], Kd [%·s/kPa]** — PID gains; **Ctrl Period**
  is the valve-update tick (default 0.10 s). Derivative acts on the
  measured pressure (no setpoint kick); the integrator has anti-windup.
  Gains persist in the machine-local calibration cache like the other
  factors.
- **Run Step Test (Dynamic Cal)** — steps the free BPA from its current
  pressure to the setpoint, plots the response on the Calibration Result
  axes with the ± deadband, then reports **overshoot and undershoot in
  kPa and as % of the step**, 10–90 % rise time, settling time (last time
  outside the deadband), and steady-state error — on the message line and
  in the command window. Each test is saved as `DPC_S##_R##.mat` in the
  save folder (`DPC_Data` = t, kPa, V, duty %, valve state;
  `DPC_Metrics`; `DPC_PID`). Click **Pause** to abort a test early; the
  valves always park on HOLD afterwards.
- **PID servo during Get Data** (checkbox) — the same PID holds the
  setpoint while force data is acquired, one valve update per sample.
- Tuning procedure: start with Kd = 0, raise Kp until the response is
  fast with modest overshoot, add Ki to remove steady-state error, then
  add a little Kd if the approach rings. Compare overshoot/settling
  numbers across step tests (each run saved) to confirm improvements.

## Data acquisition

- Get Data / Pause / Save; unit N, g, kg, kN; session-time limit and
  sample-period controls; Samples caps one run (≤ 750).
- Pressure is read **every sample** (the old app froze it at connect time —
  fixed).
- **Pressure servo** (Data Acquisition tab, default OFF): while acquiring,
  the valves hold the setpoint ± deadband
  (Increase = D11 High + D6 High; maintain = D11 Low + D6 High; vent/decrease
  = D11 Low + D6 Low). On stop the valves park on *maintain*. The three
  **Valves** buttons give the same states manually.
- **Metadata tab**: Knee angle and Load-cell angle [deg] are recorded into
  every saved file (type them from the encoder display before Save).

## Saving

Save writes a pair of files, `<Prefix><Series>_<Run>` in the save folder
(default prefix `FlxTest`; Run auto-increments after each save):

- `.mat` — `Data` (750×6: Time_s, RawHX711_counts, Force_N, Pressure_kPa,
  PressureVoltage_V, Force_SelectedUnit), `Stats` table, `Metadata` struct
  (angles, pins, calibration factors, servo settings, valve logic, MATLAB
  version), `ColumnNames`.
- `.txt` — same two-column format the original app wrote
  (`Time [ s ]`, `Force [ unit ]`), so existing text-file tooling keeps
  working. A Save that overwrites asks first.

## Offline self-test (no hardware)

```matlab
test_HX711_BPA_offline
```

Constructs the app and checks: pre-Connect guards, known-factor apply,
conversion math, MAT+txt save content, run auto-increment, save-folder
resolution, and Clean. Ends with `PASS`/`FAIL` in the command window. Run
it on any machine after pulling.

## What happened to the old apps (2026-09-27)

| File | Status |
| --- | --- |
| `HX711_BPA.m`, `Start_HX711_BPA.m`, `test_HX711_BPA_offline.m` | **the app** |
| `HX711.mlapp` → `HX711_customized_original.mlapp.disabled` | old customized app (UI redesign was never wired up: dead calibration tabs, pressure frozen at connect). Kept for reference; rename back to `.mlapp` to open it in App Designer. |
| `HX711_Pressure.mlapp` → `HX711_Pressure.mlapp.disabled` | the other app in this folder, per Ben's request; kept, disabled. |
| `Code\Matlab\HX711-LoadCell\HX711 v3.0\HX711.mlapp` → `.disabled` | clean File-Exchange original; renamed so no `HX711` name is left on disk to shadow anything. |
| `Code\Matlab\HX711-LoadCell\HX711 v3.0\hx711_redesign_matlab2025\HX711.m` → `.superseded` | earlier redesign draft this app supersedes. |
| Zips (`github_repo.zip`, `HX711 v3.0.zip`, `hx711_redesign_matlab2025.zip`) | untouched archives; not on the MATLAB path. |

**2026-10-03: the whole `Code\Matlab\HX711-LoadCell\` folder (clean File-Exchange
original, redesign draft, archives) was DELETED at Ben's request — its inner
"HX711 v3.0" directory name collided with the live app folder and caused
confusion. Everything in it is recoverable from git history
(`git checkout <commit> -- "Code/Matlab/HX711-LoadCell"`). The window now
also shrink-fits + centers itself on the launch display (`fitToScreen` in
`HX711_BPA.m`): the design grid is 1120x900 px and is uniformly scaled
(component positions + fonts, floor 8 pt) when a screen is shorter, so the
app is usable on laptop displays.**

License: original app copyright 2018 Nicholas Giacoboni (BSD); AARL/BPA
modifications as noted in the app's License tab.
