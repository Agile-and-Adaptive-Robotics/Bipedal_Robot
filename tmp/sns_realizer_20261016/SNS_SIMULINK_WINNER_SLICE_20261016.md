# SNS_Simscape winner realization — Simulink slice of round8_best3 (2026-10-16 night)

Stage: Simscape realizer (overnight walk campaign). Machine: easteregg2,
MATLAB R2025a, repo D:\GitHub\Bipedal_Robot (branch KneeTestSetup_BenBo_stw,
working tree — NO commits, NO pushes, per campaign contract).

## Deliverable

`Code\Matlab\SNS_Simscape\demos\SNS_WinnerSlice_R8.slx` (+ generator
`sns_build_winner_slice.m`, verifier `sns_run_winner_slice.m`, README
`README_winner_slice_20261016.md`) — a reduced but faithful slice of the
WINNER human circuit (round8_best3.json = curr_s8b_wean04 trial 35,
drive_scale 0.3807379786697019) realized in the REDESIGNED SNS_Library
architecture (2026-09-22: one-input synapses -> `[g; g*Esyn]`, internal
neuronal summation on syn1..syn6, SynSum for >6 incoming).

Scope (per the night's ask): one RG half-center pair per side with the
L/R alpha+gamma hand-drive inputs, one PF layer, heel/toe afferent resets
for one leg (both legs' reset INs included so the left RG pair sees its
production input class), KINH/IBEXC/LBIN stance machinery, and one knee
MN pair (vas_lat_r / bifemsh_r) driving a demo load = scope-signal
knee-torque proxy (NOMINAL gains 190/60 N·m — placeholders, documented).
NO Simscape Multibody port attempted (explicitly out of scope tonight).

## How the values were obtained (no hand-typing)

`spinal\_slice_probe_20261016.py` loads the winner EXACTLY the way
`wean_tune_20261016.py cmd_one` does (`apply_config(knobs, scale, fb)` =
set_stage(4, base+fb) + knobs x scale + alpha_to_mn), builds the real
winner network (606 neurons / 6908 synapses / 386 inputs on the cvt3
model actuator list), and dumps the 34-neuron / 75-synapse slice subgraph
with all EFFECTIVE values (units mV/nA/uS/nF, the export_network_json.py
mapping) to `slice_export.json`. Effective winner values observed in the
probe log: descend_to_rg_e 1.698797 / f 0.996393, rg_to_pf 2.600515,
rg_mutual_inh 4.0, pf_recip_inh 4.0, pf_to_mn 2.0, posture_to_mn 1.0,
posture_to_rg_e 0.8, rg_nap_h 0.842083 s, c1_gain 0.1 (FORCED by
no_cross=1; the json's 0.4979 is overridden by set_stage), heel_rge
1.482757, toe_rge 0.360005, ib_rge 0.312219, ib_e_central 0.139386,
renshaw 0.5, f1_kneext_inh 0.178469, ib_exc_to_mn 0.6; knob G values
gamma_l 0.51157 / gamma_r 0.30362 / alpha_l 0.32862 / alpha_r 0.37062 /
alpha_to_mn 0.87210; walk drive 2.0773198435110425 nA; RG NaP (g_ion 12,
e_ion 8, m: k 1 / slope 0.8 / e 2, h: k 1 / slope -2 / e 3.5, tau_h
FIXED 842.08 ms, h0 0.000911 read from the compiled backend).

The RG somata are a purpose-built masked `NapNeuron` subsystem following
the redesigned-library port conventions (SNS_Library has no persistent-Na
block; the Deng HCNeuron predates the redesign and lacks syn ports).
Promotion of NapNeuron into SNS_Library.slx is left for Ben. The build
asserts the R2025a library copy really is the redesigned one (7-in/2-out
neuron, 1-in/1-out synapse).

The 4 stacked IBEXC->LBIN edges (audit finding N1, reported-not-fixed) are
merged into one 2.0 uS synapse — algebraically exact (parallel
conductances sum in the backend).

## Empirical findings along the way (probes 1-5, all logged)

1. **The winner RG does NOT free-run at constant drive with zero
   afferents**: both half-centers sit on the NaP plateau (RG_E_r
   5.10..5.19 mV tonic, zero bursts at gamma_r x0.5/x1/x2; probe1) — the
   documented s3k-lineage "schedule/afferent-sustained" property. The
   production rhythm (round8_full3.npz: RG_E_r 6 bursts in the walk
   window, period 1.7284 s) is contact-afferent-driven.
2. A stylized bilateral antiphase foot-contact script (production cadence
   1.7284 s, runner-shaped port currents incl. contact_onset 0.8173 edge
   kicks with 0.25 s decay) ENTRAINS the winner network 1:1 in numpy
   (probe2: RG_E_r 5 bursts, period 1.7285 s, duty 0.315; RG_F antiphase
   0.630). The same script drives the Simulink model (Repeating-Sequence
   blocks, gated to zero before the t=5 s walk onset).
3. **Modulation sweep (probe3, full net)**: gamma_r x4 flips the rhythm
   (duty 0.315 -> 0.643, near-full-swing -0.20..7.05 mV); drive x1.25
   breaks entrainment (3:1 locking); alpha_r x0.5/x2 shifts only the MN
   level (2.902 / 3.424 mV full-net); script cycle 1.4 / 2.1 s entrains
   1:1. These picked the verification levers.
4. **Renshaw network effect**: full-net MN_vas sits at 3.089 mV but the
   2-pool slice at 2.247 mV — in the full winner every RC is inhibited by
   ~46 same-side pools. Probe5 (exact slice-topology stepper on the
   exported subgraph, anchored by the matching RG stats) gives the
   slice-equivalent targets: MN_vas 2.055 / 2.247 / 2.591 mV at alpha
   x0.5 / x1 / x2. Simulink reproduces these EXACTLY (to 3 decimals).

## Verification actually run (headless egg R2025a, exit 0)

`matlab -batch sns_run_winner_slice` — model simulates 22 s at ode1/2 ms
(production Euler semantics), 7 sims total. Observed output:

```
CHECK A base: cycles=5 period=1.7285 duty=0.316 v=[5.028,5.745] | RG_F duty=0.630 v=[-2.610,-1.959]
CHECK A PASS (targets: 5 cyc, 1.7285 s, duty 0.315, vmax 5.75)
CHECK B gamma x4: cycles=5 period=1.7155 duty=0.662 v=[-0.234,7.067]
CHECK B PASS (numpy target: duty 0.643, vmin -0.204 - the rhythm is MODULATED by the GAMMA hand-drive knob)
CHECK C MN_vas mean: a0.5=2.055 a1=2.247 a2=2.591 mV
CHECK C PASS (slice-ref targets 2.055 / 2.591; full-net ref 2.902 / 3.424 for the 46-pool Renshaw network)
CHECK D cycle 2.1: cycles=5 period=2.1005 duty=0.321
CHECK D PASS (heel/toe afferent resets entrain the RG 1:1 to the contact cadence)
=== WINNER SLICE VERIFICATION: ALL CHECKS PASS ===
```

Numpy references: period
1.7285 (slice & full net), duty 0.315/0.316, V-ranges agreeing to
~2e-3 mV; gamma x4 slice-ref duty 0.641, vmin -0.231; alpha means
2.055/2.591. Figure + results mat:
`winner_slice_verify.{png,mat}` (this folder).

## What is NOT in the slice (documented, by design)

Afferent encoders (Ia/II/Ib), BAL_*/VEST cells, the pm_* runner phase
machine, contra_* crossed edges (no_cross=1 makes them 0), the MuJoCo
body, and the runner's MuJoCo standing-solve POST_ biases (ports present,
default 0 — they cannot affect the RG/PF rhythm; MNs project only to
their Renshaw cells). The torque-proxy gains are nominal placeholders.
NOT run tonight: opening the model in the Simulink GUI (headless build +
sim only), Simscape Multibody anything, the laptop R2025b copy of the
model.

## Files (all on easteregg2, uncommitted)

- `Code\Matlab\SNS_Simscape\demos\`: SNS_WinnerSlice_R8.slx (new),
  sns_build_winner_slice.m, sns_run_winner_slice.m,
  README_winner_slice_20261016.md (new)
- `Code\Matlab\SNS_Simscape\README_SNS_Simscape.md` (TRACKED file — short
  section appended; see CHANGES_20261016_sns_simulink.md)
- `Code\MuJoCo_SNS\spinal\`: _slice_probe{,2,3,4,5}_20261016.py (evidence
  + reusable references)
- `Code\MuJoCo_SNS\spinal\reports_20261016\sns_simulink\`: slice_export.json,
  slice_ref.npz, slice_ref2.npz, slice_ref_targets.json,
  winner_slice_verify.{png,mat}, this report

## Reproduce

```
ssh easteregg2
cd /d D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos
matlab -batch sns_build_winner_slice   % rebuild + 0.2 s smoke
matlab -batch sns_run_winner_slice     % 4-check verification
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
D:\Anaconda\envs\myo\python.exe _slice_probe_20261016.py    % re-export
```
