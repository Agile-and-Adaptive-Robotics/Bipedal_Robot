# BPA Walker Program — Workflow, Architecture, and Test Plan

Revision 2026-09-26 (AARL / Ben Bolen). This is the master program document that ties the
sensors, actuators, simulation stacks, and the physical test rig into one workflow. It is the
contract between the three "hosts" we develop against: **Simscape Multibody + SNS**, **MuJoCo +
sns-toolbox**, and the **Orin Nano hardware loop**.

---

## 0. What already exists (build on, don't rebuild)

| Asset | Location | State |
|---|---|---|
| Spinal SNS network (RG+PF, 92 MN pools, Ia/II/Ib, contact resets, vestibular stage) | `Code\MuJoCo_SNS\spinal\` | walks (s3k winners), curriculum stages 3/4/5 |
| Muscle conversion formulas (Ia ∝ L̇/0.6, II ∝ (L−Lmid)/Lhalf, Ib ∝ F/Fmax) | `spinal\DESIGN.md`, Deng-style figure | audited |
| BPA force law (Festo sfit, maxBPAforce) | `Code\MuJoCo_SNS\bpa_muscle.py`, `SNS_Simscape\SNS_Library.slx` (BPA_10/20/40mm blocks) | validated vs festo4.m |
| Knee rig in Simscape, BPA+SNS driven, antagonistic | `Code\Matlab\SNS_Simscape\mdl_leg_rig_ba003_imported.slx` | working (corr −0.17, 999 N) |
| Lower-body humanoid ah001 in Simscape (73.0 kg, feet, spine, 2 SNS knees) | `Code\Matlab\SNS_Simscape\mdl_humanoid_lower_ah001_imported.slx` | compiles + sims |
| OpenSim bone → Simscape parts pipeline | `Solid_Models\Simscape_Part_Library\OpenSim_Bones\`, `sns_bone_leg_demo.slx` | proven |
| gait2392 27-BPA actuator map (names, Fmax, routes) | `Solid_Models\OpenSim\Gait2392_Robotbody\bpa_actuators_27.json` | extracted |
| Gait reference library (43 refs: Falisse, Arnold, Ong) | `Code\MuJoCo_SNS\spinal\gait_refs\` | loaded via gait_lib_loader |
| Synergy basis (6-synergy NMF from back-solved activations) | `spinal\synergy_basis.npz` | VAF 0.945 |
| Ty's ankle-foot CAD (final = "Newest (7-17-20)", reorganized as `01_00_00_Foot`) | `Solid_Models\Prospective_Parts\01_00_00_Foot\` | parts all present; assembly needs link re-path (no Pack-and-Go) |

---

## 1. Goal ladder

| Goal | Deliverable | Exit criteria |
|---|---|---|
| **G1** — sagittal walker, uniarticular | 2-D-ish sagittal walker in MuJoCo-SNS (first) and SNS-Simscape (parity), 5 uniarticular muscle groups/side (hip flex, hip ext, knee ext, ankle PF, ankle DF) | 10+ steps on level "ground" without harness in sim; falls are clean and diagnosable |
| **G2** — full walker, ah001 + 27 BPAs/side | ah001 Simscape model with all 27 muscle routes/side driven by the spinal SNS | stands untethered in sim; ≥3 gait cycles at ≤0.8 m/s; kine score vs gait refs reported |
| **G3** — design freeze | valve manifold map, sensor BoM + mounts, foot choice (Ty's 01_00_00 or OpenSim-metric foot), wiring schematics, Orin carrier plate | Ben signs off; every BPA has a valve channel + pressure sensor; every joint has encoder + Liquid Wire |
| **G4** — simulation validation | Same controller binary/config runs in Simscape AND MuJoCo with matched afferent encodings; sensitivity sweep (gains ±20 %, delays 5–40 ms) | controller survives both stacks within spec'd tolerance; sim-to-sim delta documented |
| **G5** — physical test (sagittal walker, Orin on-board, computer-in-loop) | Walker walks on a boom; live neural deletions + stimulus injection executed from the lab console; logged dataset | ≥20 successful deletion/stimulus trials logged with synchronized neural + sensor traces |

Order of work: G1 → (G3 partial: manifold + sensors on bench) → G4 → G2 in parallel with G5 build.

---

## 2. The one-interface principle

Everything plugs into **one plant/controller contract** regardless of host:

```
        sensors (physical or simulated)
          │  mm, deg, kPa, m/s², contact bools
          ▼
    ┌──────────────┐   afferent encoding (II/Ia/Ib/contact/vestibular)
    │ afferents.py │ ────────────────────────────────┐
    └──────────────┘                                 ▼
   ┌──────────────────────── SNS network (same wiring in all hosts) ───────────────┐
   │ RG+PF, MN pools, Renshaw, IaIN, IBEXC, VEST … gains = the current best json   │
   └───────────────────────────────┬──────────────────────────────────────────────┘
                                   │ per-BPA activation a∈[0,1] (or valve duty)
                                   ▼
   ┌───────────────┐  fill/vent  ┌────────┐  pressure  ┌────────────────────┐
   │ valve driver  │────────────▶│  BPA   │───────────▶│ plant (Simscape /  │
   │ (Orin PWM)    │             │ muscle │            │ MuJoCo / hardware) │
   └───────────────┘             └────────┘            └────────────────────┘
```

Contract files live in `Code\Hardware\walker_io\` (Python, used on the Orin and as the
reference for both sims). The sim hosts implement the same `Plant` interface (see
`skeletons` in the doc §10) so G4 is a config swap, not a rewrite.

---

## 3. Sensor → afferent mapping (the heart of the workflow)

| Physical sensor | Signal | SNS afferent class | Encoding (from spinal conversion map) | Notes |
|---|---|---|---|---|
| **Joint encoder** (per DOF) | θ, θ̇ | **Ia** (primary spindle) | Ia rate ∝ θ̇ / θ̇_max, sign = stretch of the antagonist group | also joint position state for the force law |
| **Liquid Wire** (per BPA) | muscle length L (real-time) | **II** (secondary spindle) | II rate ∝ (L − Lmid)/Lhalf, tonic component also gates MN bias | Leo's sensor is the only per-muscle length source — trust it over geometry |
| **Pressure transducer** (per BPA line) | P (kPa) | **Ib** (Golgi → force) | Ib ∝ F/Fmax where F = festo4(P, L) | autogenic inhibition + stance-gated load sharing (already in the spinal net as IBEXC) |
| **Insole capacitive array** (develop in house) | heel/toe contact + CoP | **contact mechanosensors** | HEEL_c/TOE_c edges → S2W phase reset of the ipsilateral RG (exists: heel/toe gains) | start with 2 zones (heel/toe), expand to CoP later |
| **IMU** (pelvis ± trunk) | orientation, ω, accel | **vestibular down-command** | VEST cells (stage-4 balance work) drive trunk/hip balance | 2nd IMU + 3D camera correct drift |
| **3D camera** (machine vision) | base pose in world | **drift correction / descending drive** | slow correction of pelvis state estimate; terrain → DRIVE modulation | lab-space only at first; not on the walker |

Rule: **no sensor feeds the SNS raw**. Every signal passes `afferents.py`, which owns the
normalization constants (per-muscle Lmid, Lhalf, Fmax from `bpa_actuators_27.json`), sign
conventions, and units. The same file is the single source of truth for all three hosts.

---

## 4. Muscle & valve mapping

### 27 BPAs per side (from `bpa_actuators_27.json`)

3 trunk + 8 thigh/knee + 3 ankle PF + 5 subtalar (incl. 2 tib-post routes) + 4 toes + 4 hip.
Names, rest lengths, Fmax, and route points are already in the JSON — `afferents.py` loads it
and exposes `muscles[i] -> (name, joint, Fmax, Lmid, Lhalf)`.

### Valve drive

- One Festo manifold channel per BPA line: **fill** (supply) and **vent** (exhaust) solenoids.
  Two 3/2 valves per BPA (or 5/2 per antagonistic pair where plumbing allows).
- **Modes** (both proven by you; the driver supports both):
  - `bang_bang`: valve open/closed with hysteresis on pressure error (Orin inner loop, 100 Hz).
  - `proportional_taped`: fixed short PWM packets (duty = activation) — treat valve as
    proportional orifice per your taping; calibrate duty→pressure-per-second on the bench and
    store the table in `valves.py`.
- **Mapping activation → valve command**: per antagonistic pair (e.g. VASTI vs HAM),
  `duty_fill = a_agonist`, `duty_vent = a_antagonist` with co-contraction floor 0.05.
  Pressure closed loop sits **below** the SNS: the SNS commands activation, Orin's fast loop
  realizes pressure. Never let the pressure loop fight the afferents — Ib sees measured force,
  not commanded.

### G1 muscle set (uniarticular sagittal)

| Group | gait2392 analog | BPAs |
|---|---|---|
| hip flex | iliopsoas | 1 pair? — 1 BPA (agonist only) + passive return at first |
| hip ext | gluteus maximus | 1 BPA |
| knee ext | vasti (VasMed+VasLat+Int) | 2 BPAs |
| ankle PF | soleus | 1 BPA |
| ankle DF | tibialis anterior | 1 BPA |

10 BPAs total → 10 valve channels. Foot: Ty's `01_00_00_Foot` (preferred — real hardware) or
OpenSim `foot.stl` from the part library for the pure-sim G1.

---

## 5. Neural deletions & stimulus injection (the experiment interface)

The spinal SNS is built with **conditional topology** — every pathway is a named gain in
`params.G` (default 0 = absent). The experiment interface reuses exactly that:

- **Deletion** = set a named gain/pool output to 0 live. Examples: `ia_in=0` (kill IaIN
  reciprocal inhibition), `ib_rge=0` (kill load sharing), per-pool MN silence (`pool_off =
  ["vasmed_r","vaslat_r"]`), Renshaw off, commissural off (`--no-interleg` analog).
- **Stimulus injection** = add current to named neurons: `stim = {"HC-RG-E_R": 2.0, "t": 12.0,
  "dur": 0.2}` (nA, time, duration) — the same PRESET/VEST input ports already exist in the
  built network.
- Both are expressed as **ephemeral overrides** on top of the current best-parameter json —
  the base config is never edited. Overrides are logged with the telemetry so every trial is
  exactly reproducible offline in MuJoCo/Simscape.

---

## 6. Computer-in-the-loop protocol (Orin ⇄ lab PC)

Link: **UDP over Ethernet** (direct cable, 1 GbE). Orin = `192.168.10.1` (walker), lab PC =
`192.168.10.2`.

| Stream | Direction | Rate | Format | Content |
|---|---|---|---|---|
| telemetry | Orin→PC | 50 Hz (10 ms batches) | JSON lines | t, phase, all pool activations, all afferents, all pressures, joint θ/θ̇, IMU, contact |
| command | PC→Orin | on demand | JSON | `{"cmd":"override","G":{"ia_in":0}}`, `{"cmd":"stim","neuron":"HC-RG-E_R","I":2.0,"dur":0.2}`, `{"cmd":"valve_override","name":"VASTI_R","duty":0.4}`, `{"cmd":"stop"}` |
| heartbeat | Orin→PC | 10 Hz | 8-byte seq | watchdog: **loss of PC heartbeat ≠ stop** (walker must survive console loss); loss of *sensor* heartbeat = vent-all |

Safety rails (non-negotiable, implemented in `orin_daemon.py`):
1. Any sensor stream stale > 20 ms → neutral valve state (all vent), log.
2. Pressure hard caps per BPA (620 kPa supply reg + software cap) enforced below the SNS.
3. Physical e-stop cuts manifold supply pneumatically, independent of Orin.
4. Walker runs on a boom/harness until G5 exit; boom load logged.

---

## 7. G1 sagittal walker spec (start here)

- **Plant**: planar legs (hip pin, knee pin, ankle pin) — bodies from the OpenSim part library
  (`femur_r.stl`, `tibia_r.stl`, `foot.stl`) at Onyx density, so G1 geometry ≈ final robot.
- **Muscles**: the 5 uniarticular groups/side (§4) as point-force BPAs (same Internal-Force
  construction as the rig demo).
- **Controller**: the spinal SNS at reduced pool set (hip-F, hip-E, knee-E, ankle-PF, ankle-DF
  per side) + contact resets + Ib; walk in MuJoCo first (fast iteration, existing stack),
  port to Simscape for G4 parity.
- **Reference**: Ong self-selected + Falisse walking from `gait_refs` for the kine score.
- **Tuning**: reuse the curriculum idea (air → tethered → ground) so sim training transfers.

---

## 8. G5 physical test plan

### Bench bring-up (per subsystem, walker clamped)
1. Pneumatics: regulator → manifold; every channel crack-test at 200 kPa; duty→dP/dt table per
   BPA (fills `valves.py` calibration).
2. Sensors: encoder direction/sign per DOF; Liquid Wire calibration (length vs ADC, per BPA);
   pressure transducer span/offset; IMU alignment.
3. Afferent sanity: `lab_console.py` live view — stretch a BPA by hand, see Ia/II spike; press
   the foot, see contact reset fire; load the knee, see Ib rise.

### Boom walking (the G5 experiment)
- Walker on an overhead boom (constant-height harness, sagittal freedom + small lateral
  play), over a 3 m walkway.
- Trial matrix (each 20 s, 5 repeats): baseline → single deletions (Ia, II, Ib, contact,
  vestibular, commissural, Renshaw) → stimulus injection (RG-E, RG-F, MN pools at 3 amplitudes)
  → combinations suggested by the deletion results.
- Every trial = one `trial_id`; telemetry + override log land in
  `Testing_Data\walker_trials\<date>_<trial_id>\`.

---

## 9. Schematics

### Pneumatic (per antagonistic pair; ×27 pairs/side as routed)

```
  supply bottle ──regulator (620 kPa)──shutoff──┬────────────┬───────────────
                                                │            │
                                          fill valve A   fill valve B
                                                │            │
                                             BPA agonist  BPA antagonist
                                                │            │
                                          vent valve A   vent valve B
                                                │            │
                                             exhaust (muffler box)
  pressure transducer tees into each BPA line (analog 0.5–4.5 V → ADC)
```

### Electrical / compute (Orin Nano carrier)

```
  Jetson Orin Nano dev kit
   ├─ PWM/GPIO (74HC595 or PCA9685) ── MOSFET boards ── 24 Festo manifold solenoids/side
   ├─ ADS1256/ADS8688 ADC (SPI) ── pressure transducers (10/side) + Liquid Wire dividers
   ├─ SPI/I2C ── IMU (pelvis) + IMU (trunk)
   ├─ GPIO int/quadrature ── joint encoders (hip/knee/ankle, both sides)
   ├─ USB ── insole array front-end (MCU streaming contact/CoP)
   └─ Ethernet ── lab PC (computer-in-the-loop console)

  Power: Li-ion pack → 5 V/20 A buck (Orin) + 24 V (valves) + 12 V (sensors), common ground
  E-stop: hardware latching relay in the 24 V supply line to the manifold shutoff
```

### Compute / experiment flow

```
  Lab PC (console, deletions, plots)        Orin Nano (walker)
   ── UDP command (overrides/stim) ────────▶ SNS 1 kHz loop
   ◀── UDP telemetry 50 Hz batches ───────── afferents → SNS → valves → plant
   ◀── CSV/git-annexed trial logs ──────────▶ same contract in Simscape/MuJoCo for replay
```

---

## 10. Repository map (this program)

```
Code\Hardware\walker_io\          protocol, valves, sensors, afferents, orin_daemon, lab_console
Code\MuJoCo_SNS\spinal\           SNS network + curriculum (G1/G2 controller source of truth)
Code\Matlab\SNS_Simscape\         Simscape hosts: rig + ah001 + bone demo (G4 parity)
Solid_Models\Simscape_Part_Library\OpenSim_Bones\   bone STLs for Simscape File Solids
Testing_Data\walker_trials\       physical-trial logs (G5)
Documentation\Program_Workflow\   this document
```

Skeleton code (marked `SKELETON` in each file header) implements the protocol, the valve/
sensor abstractions with hardware TODOs, the afferent encodings (real formulas, constants
from `bpa_actuators_27.json`), and both ends of the UDP link. Sim hosts implement the same
`Plant` duck-type: `state()` → dict, `apply(action)` → None, `step(dt)`.

---

## 11. Open decisions (Ben)

1. Foot for G1/G2 sims: Ty's `01_00_00_Foot` (needs SolidWorks link re-path + Multibody-Link
   export) vs OpenSim `foot.stl` (available today). Default: OpenSim foot for G1, Ty's foot for
   the physical walker.
2. Valve channel count vs BPA count: 2 channels per BPA (54+ channels) or shared
   fill/vent-pair plumbing (fewer)? Decide with the manifold part number in hand.
3. Liquid Wire analog front-end: which ADC + excitation (Leo's paper hardware?).
4. Insole array: start with commercial dev kit or straight to in-house capacitive cells?
5. R-leg muscle mirror fix in the ah001 model (half-sweep asymmetry noted 2026-09-26).
