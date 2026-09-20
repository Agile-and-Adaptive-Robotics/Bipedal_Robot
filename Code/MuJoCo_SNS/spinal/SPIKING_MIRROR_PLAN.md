# Spiking-mirror plan for the gait2392 spinal network

Ben, 2026-09-18: "How to make a mirror copy of our existing model (once
you wire the PF layers) but using spiking neurons instead?"  The
non-spiking joint-layer PF wiring landed behind `G["joint_pf"]` on
2026-09-18; this is the plan for the spiking twin.

## Why a hybrid, not an all-spiking clone

* The plant is analog: MuJoCo muscles take a continuous activation, and
  MN pools are rate-coded populations.  Keep the MN layer NON-spiking
  and drive it through `SpikingNonSpikingSynapse` — the standard
  SNS-Toolbox hybrid pattern (Tutorials 5-8): spikes above, analog
  voltages into the plant.
* Spikes matter where timing and event coding matter: afferent encoding,
  RG/PF/IN classes, commissurals.  Those become spiking.

## Cell-class mapping (non-spiking -> spiking mirror)

| current (non-spiking)                 | spiking mirror                                                            |
|---|---|
| afferent encoders Ia/II/Ib (analog current in, V out) | `SpikingNeuron` (LIF, threshold + reset); input current -> firing rate via a linear I->rate map chosen from the existing afferent gain ranges |
| RG-E/F (persistent-Na half-centers, fixed tau_h) | LIF + spike-frequency adaptation, or the toolbox's spiking classes from Tutorials 5-8; FIRST verify what spiking burst machinery 1.5.2 actually ships (the NaP class is non-spiking-only) — if none burst natively, use adapting LIF pairs with mutual inhibition (classic half-center) |
| PF layer HCs (joint-layer or phase cells) | adapting LIF, same E/F lamination topology |
| INs (InE/InF, PF_IN, CIN, IaIN, KINH, RC, LBIN, IB-EXC) | LIF, chemical synapses; Renshaw = inhibitory LIF driven by MN-layer V through a rate encoder |
| MN pools + POSTURE bias | stay `NonSpikingNeuron`; spikes -> `SpikingNonSpikingSynapse` (weight = PSP-equivalent conductance, calibrated so the MN V range and `a = clip(V/5mV,0,1)` map are unchanged) |
| Heel/toe mechanosensors | already event-like: become literal spike sources (contact onset -> spike burst) — this is the cleanest win of the mirror |

## Parameter-mapping rules (v1 starting points, to be calibrated)

1. Voltages: spiking cells get real mV ranges (rest -70, threshold -50,
   reset -60); non-spiking 0-5 mV range maps to rate/charge, not volts.
2. Synapse signs: excitation reversal ~0 mV, inhibition ~ -70 mV;
   weights (µS) start at the non-spiking conductances rescaled so a
   single-PSP steady train reproduces the non-spiking steady-state
   MN current (calibration script, not hand-tuning).
3. Time constants: keep the same ms values; the loop dt must drop to
   0.1-0.5 ms (2 ms cannot resolve spikes) -> expect 4-10x slower
   wall-clock; consider sub-stepping the network inside the 2 ms plant
   step (network x4-20, plant x1) to keep plant cost unchanged.
4. Topology: reuse build_network structure verbatim (same names, same
   conditional gains) in `build_network_spiking.py` so `_fix_check`-
   style topology gates can be mirrored.

## Verification gates (mirror the existing ones)

1. Topology: same connection-class contract as `_fix_check.py`.
2. Rhythm: constant-DRIVE self-sustained alternation (the
   `check_selfsustain` recipe), period within ~20% of the non-spiking
   build at matched DRIVE.
3. Basin: the `basin_gate.py` perturbation protocol on the spiking RG.
4. Air smoke: `_smoke_nap`-style 20 s MuJoCo loop; compare RG period,
   E-duty, knee/hip ranges to the non-spiking run (not bit-exact — the
   gate is regime match, not identity).

## Effort + risks

* 1-2 sessions for the build + rhythm gate; +1 for the air smoke.
* Risks: (a) 1.5.2 may lack a spiking bursting class — fallback is
  adapting-LIF half-centers (well supported, changes burst shape);
  (b) dt cost; (c) all tuning is new — the mirror starts from mapped
  parameters, not from a study; expect a fresh (short) curriculum.
* Payoff: event-coded afferents (real mechanosensor spikes), spike-
  timing-dependent phenomena later (Hebbianness), and a check that our
  non-spiking conclusions are not artifacts of analog coding.
