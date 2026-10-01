# SNS Simscape rebuild + validation (laptop) — W2L 2-layer CPG

**Date:** 2026-10-01 (laptop DESKTOP-5Q16KE9, MATLAB R2025b, python `myoconv`)
**Verdict: PASS** — `W2L_SIMSCAPE PASS | period L 1.000 vs 1.000 s (0.0%) R 1.000 vs
1.000 s (0.0%) | antiphase r -0.525 (ref -0.525) | bursts 10/10` (log line 6 of
`Code/Matlab/SNS_Simscape/logs/laptop_20260930/sns_w2l.log`).

## Deliverables

| Artifact | Path |
|---|---|
| Simulink model | `Code/Matlab/SNS_Simscape/results/SNS_W2L_CPG.slx` (301,949 B) |
| Builder (new dev m-file) | `Code/Matlab/SNS_Simscape/dev/build_sns_w2l_cpg_20260930.m` |
| Run log | `Code/Matlab/SNS_Simscape/logs/laptop_20260930/sns_w2l.log` |
| Simulink traces + verdict | `Code/Matlab/SNS_Simscape/results/SNS_W2L_CPG_run_20260930.mat` |
| Numpy reference driver + data | `spinal/campaigns/20260930/w2l_numpy_ref_20260930.py`, `w2l_numpy_ref.{json,npz}` |
| Trace comparison | `spinal/campaigns/20260930/w2l_trace_compare_20260930.py`, `w2l_trace_compare.json`, `sns_w2l_traces.png` |

## Step 1 — existing builder search (per ask)

Searched `*.m` + `dev/*.m` under `Code/Matlab/SNS_Simscape` for
json/build/connectome. A builder + net-JSON pair DOES exist:
`sns_build_from_json.m` (root) driven by `spinal/spinal_net_export.json`
(grep hits: `sns_build_from_json.m:3`, `:23`). **It could not be used for the
W2L net:** it targets the 406-neuron tuned gait2392 spinal net (its own header,
`sns_build_from_json.m:1-4`) and its entire block vocabulary is
`SNS_Library/NonSpikingNeuron` + `NonSpikingSynapse`, which cannot express the
W2L rhythm generator — the four persistent-Na RG half-centers with FIXED
tau_h (`NonSpikingNeuronWithPersistentSodiumChannel` +
`SNS_NumpyFixedTau`, `spinal/w2l_cpg/build_w2l_net.py:29-37`). Per the ask's
step 2, a new script-built representative core was built instead, reusing the
pair's proven conventions (units mapping, 2026-09-22 syn-port architecture,
fixed-step 2 ms Euler, `sns_units_test_2n.m` block-parameter style).

## What was built — SNS_W2L_CPG (REPRESENTATIVE CORE, not the full net)

The model annotation and `Description` both label it as the representative
core. Contents (all gains verbatim from `build_w2l_net.py` `W2L_GAINS`;
wiring directions verified against `spinal/connectome_templates.json` key
`bilateralrg` edge-by-edge before building):

- **4 NaP RG half-centers** (`L/R_RG_ext/flx`) — new custom masked block
  `NapHC`: `Cm dV/dt = Gm(Vrest-V) + Iapp + sum(g*Esyn) - V*sum(g) +
  GNa*m*h*(ENa-V)`; `m = 1/(1+kM*exp(sM*(eM-V)))` instantaneous;
  `dh/dt = (hInf-h)/tauH` with **tauH FIXED 250 ms** (AnimatLab/Deng
  semantics = `SNS_NumpyFixedTau`). Params = `spinal/params.py` NAP verbatim
  (GNa 12 uS, ENa 8 mV, m k1/s0.8/e2; h k1/s−2/e3.5), membrane tau 50 ms,
  0–5 mV scale.
- **Laminated mutual inhibition** per side (E→InE exc, InE→F inh, F→InF exc,
  InF→E inh, g 4.0) + **weak direct HC↔HC escape excitation** (g 0.5).
- **Crossed commissural inhibition:** c1 (RG_flx→c1 exc 1.0, c1→*contra* RG_flx
  inh 3.0) and V3 (RG_ext→V3 exc 1.0, V3→*contra* InE exc 0.08) — template tags
  `comm_c1`/`c1_SynAmp2.749`/`comm_v3`/`v3_to_contra_InE`.
- **PF E/F half-centers** per side (RG→PF 2.4; PF cross-IN laminated 4.0, IN
  driven by own PF inhibiting the antagonist — template `pf_cross` edges).
- **4 representative MNs** (PF→MN g 2.0).
- **Heel-contact drive channel** per side: graded contact neuron (tau 40 ms) +
  g 1.0 contact synapses onto ipsilateral RG/PF/MN *extensor* cells
  (representative subset of the template's 20 `contact_C20` edges).
- Totals: **4 NaP half-centers + 22 graded cells + 42 synapses** (printed in
  the log, line 1). Solver ode1 fixed 2 ms = the numpy dt exactly.

Drive protocol (= `smoke_w2l.py` defaults, identical in both engines): tonic
2 nA (E) / 3 nA (F) on all four HCs + the verbatim .aproj antiphase kickoff
pair (10 nA for 10 ms: L RG ext, R RG flx) + alternating 1 Hz heel trains
(20 nA, 30 % width, R phase-delayed 0.5 s), 12 s run with the 2–12 s window
analyzed (10 s of rhythm, per the ask).

**Honest finding that forced the contact trains:** under tonic+kickoff only,
the full numpy net does NOT oscillate — measured: R RG-E latches at max
−1.154 mV over 2–12 s, zero bursts (first run of `w2l_numpy_ref_20260930.py`,
exit 1). The W2L architecture is contact-driven (its AnimatLab template has no
tonic wiring at all; the TONIC ports are the python package's runner
convention, `build_w2l_net.py:295-303`). So "DRIVE" for this net = tonic +
contact trains, delivered identically to both engines.

## Validation numbers

**Numpy reference** (full 95-neuron/208-synapse net, `myoconv` python,
`SNS_NumpyFixedTau`, dt 2 ms):

```
W2L_NUMPY_REF PASS period_L=1.000 period_R=1.000 r=-0.525
bursts: L=10 (period 1.000 s)  R=10 (period 1.000 s)
window 2-12 s | L RG E max 6.831 mV, R RG E max 6.831 mV
```

**Simulink core** (same protocol, ode1 @ 2 ms):

```
SIMSCAPE W2L CORE: window 2-12 s | L RG E max 6.840 mV, R RG E max 6.837 mV
bursts: L=10 (period 1.000 s)  R=10 (period 1.000 s); antiphase r=-0.525
W2L_SIMSCAPE PASS | period L 1.000 vs 1.000 s (0.0%) R 1.000 vs 1.000 s (0.0%) |
                 antiphase r -0.525 (ref -0.525) | bursts 10/10
```

Gate (ask): period within 20 % both sides ✓ (0.0 %), antiphase present in both
✓ (r = −0.525 = reference; both < −0.5).

**Beyond the gate** (`w2l_trace_compare_20260930.py`): all four RG traces
cross-correlate at r = 1.000 with a −2 ms lag (one solver step — the known
update-order difference: numpy evaluates I_Na with the *updated* h, Simulink
ode1 integrates both states from step-start values); pointwise RMSE over the
full 2–12 s window 0.099–0.111 mV (~1.5 % of peak); peaks within 0.13 %;
L RG-E burst onsets **identical** in both engines: 2.09, 3.09, …, 11.09 s.

One real bug was caught and fixed during the build (recorded for honesty):
the NapHC leak Sum was first wired `'+-'` (V−Vrest — positive feedback; R RG-E
exploded to 8.2e61 mV); corrected to the library's `'-+'` (Vrest−V), after
which the run above passed.

## Honest inventory — MISSING vs the full 95-neuron W2L net

- Renshaw cells (RC own/mutual, 3 tag families) — absent.
- All afferent chains: Ia relay + IaIN reciprocal inhibition, II relays, Ib
  autogenic excitation, `aff_HC_excite` afferent→HC edges — absent (the core
  has no sensory inputs at all).
- The 2nd PF pair per side: the full net has hip+knee PF half-center pairs per
  side (4 PF HCs/side + 4 INs); the core has 1 pair + 2 INs.
- Toe contact (`SN-toe`, 10 of the 20 contact edges) — absent; heel only.
- MN outputs: 4 representative MNs vs 12 muscle outputs in the full net
  (hip/knee ext/flx per side); no muscle/plant load connected.
- The `ia_*`/`other`/`rc_*` gain families are therefore unexercised in Simulink.
- Not done: R2025a export of the model (Simscape now lives on the laptop by
  Ben's 09-30 call); MuJoCo co-simulation of the Simulink core; the Li net
  (`build_li_net.py`) was not rebuilt this session.

## Commands (exactly as run)

```
# numpy reference (PASS line above)
C:\Users\Ben\.anaconda3\envs\myoconv\python.exe w2l_numpy_ref_20260930.py
# build + validate + log (exit 0)
"C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch ^
  "addpath('...\SNS_Simscape\dev'); build_sns_w2l_cpg_20260930" ^
  -logfile "logs\laptop_20260930\sns_w2l.log"
# trace comparison + figure
C:\Users\Ben\.anaconda3\envs\myoconv\python.exe w2l_trace_compare_20260930.py
```
