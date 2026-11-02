# Optional IBEXC phase/load coincidence

`ibexc_coincidence.py` implements an opt-in engineering circuit in the existing
ten `IBEXC_*` neurons of each graded Gait family. It does not change default
builders, add neurons, enable Hunt training, or release the all-family gate.
The policy is a binary coincidence design, not multiplication or a validated
physiological model of positive force feedback.

## Actual source and electrical contract

Canonical s3k and w2lvar compiled gates sum RG_E and Ib excitation. In the
ankle-plantarflexor representative, phase alone gives 0.80 outgoing activation,
and the declared F0-equivalent pooled load alone gives 0.26944. Original syn6
has no RG_E input at these gates: the same load produces 0.26944 with either
phase. The candidate preserves source-relative Ib weights but gives each input
pool a total conductance budget, rather than a fixed conductance per muscle.

For actual gate leak G and excitatory reversal E=8 mV, the phase budget is
0.5G and the summed load budget is 0.5G. The loaded equilibrium is

    U = (Gphase*Ephase + sum(Gload_i*Eload_i)) / (G + Gphase + sum(Gload_i)).

Each input alone has U <= 2.666667 mV; both fully activated give 4 mV.
The candidate outgoing activation window is [3.333333, 3.733333] mV.
Each actual gate's reversals, leak, and strict margins are checked before
committing any changes. Unsupported incoming paths, gate bias, ionic gate
currents, ambiguous input ports, and unmatched bilateral source-current
contracts are rejected. Its clocked gate update must remain monotone.

The load current is not assumed to be 1 or 5 nA. Canonical runner code uses

    uIb = 0.5 * g_ib_effective * clip(force / Fmax, 0, 1.5)
    g_ib_effective = AFF.ib_gain * (1 - 0.7 * stance)

for extensor-class muscles. With the current canonical AFF.ib_gain=1.5 and a
declared stance value of 1, F/Fmax=1 corresponds to **0.225 nA**, or 0.225 mV
at an isolated encoder with leak 1 uS and input weight 1. Only the
Ib-to-IBEXC activation window changes to [0, 0.225] mV. Autogenic Ib-to-IBIN
and all other Ib windows are preserved. `F0` in the probe names means this
caller's normalization reference; it does not assert typical stance force or
identity with every saved controller's source-muscle F0.

Callers must supply an explicit current for each source Ib encoder, including
their actual effective gain and force-normalization convention. The supplied
probe policy is a canonical runner engineering example, not a saved-winner or
Hunt-feature binding. Under full phase and equal normalized encoder loads,
this example starts opening above pooled load 3/7 and reaches full activation
at 0.75. It is not a sensitive detector of arbitrarily small stance loads.
Expected source-force distributions must determine a later physical binding.

## Verification and graph changes

All ten gates per family satisfy the static inequality checks. The largest
fan-in representative is run as an SNS-Toolbox encoder/gate/target circuit,
copying actual compiled parameters, with the graded FSA encoders installed.
The RG_E source is clamped to its declared on/off voltage for this electrical
test; the original NaP rhythm is not replaced or tested as a body controller.

The candidate truth table is 0/0/0/1 for neither/phase/load/both. With phase
turned off after 0.5 s while load remains at 1, 1.5, or 5 times the declared
F0 scale, the gate transmission reaches zero in 24 ms and the delivered
target synaptic current reaches zero in 26 ms. The source current is
additionally multiplied by 10/3 during phase-off to stress the runner's
possible gain increase. Five times F0 is an out-of-range saturation stress,
not an assertion that the runner's 1.5 force-ratio clip is removed.
The target neuron retains its normal membrane decay after current release.

Source-weight allocation, graph changes, and all untouched other Ib routes
are enumerated in each evidence JSON. s3k/w2lvar need no new phase edge.
Canonical-right syn6 explicitly gains ten RG_E-to-IBEXC phase edges. The
candidate otherwise preserves neuron capacitances and dynamics, including
the NaP RG/PF circuits. Actual syn6 source NMF K6 is unchanged.

**w2lvar retains 216 direct central Ib routes bypassing IBEXC.** The inventory
finds no such direct central routes among these Ib sources in canonical s3k
or syn6. All 54/270/54 other Ib routes are preserved, respectively. Closing
this gate therefore does not prove whole-network or biomechanical unloading.

## Consumption order

1. Build the selected graded family normally.
2. For syn6, apply the explicit canonical-right candidate policy.
3. Apply the graded FSA profile, then `apply_coincidence(net, policy)`.
4. Construct the sparse accelerator after the topology change.
5. Reapply the graded FSA profile to synchronize final caches.
6. Check `verify_profile(net)` and `verify_coincidence(net)`, including after
   any existing bilateral projection. In the tested candidates that projection
   changes zero fields.

`apply_coincidence(net, None)` is a no-op. Policies reject unused keys.
Final compiled conductances, reversals, windows, encoder parameters, and input
bindings are fingerprinted. Deliberate later overwrite is detected. Dense
and sparse 128-step neural traces are required to match exactly.

There is no injected Hunt dynamic-update hook. The explicit adapter described
below verifies the actual existing update policy instead. A later changed gain
set or runtime method requires a renewed check; the field verifier detects
overwrites. This module installs no body controller, and no Hunt fitting,
curriculum, or body trial was performed.

Evidence is stored outside the repo at
`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\ibexc_coincidence\`:
three `*_POLICY.json` files, three `*_COINCIDENCE.json` files, and a compact
checkpoint/source snapshot. All-family campaign release remains false.

## Explicit Hunt/PT adapter and actual runtime composition

Policy family `Hunt_PT` is accepted only for the actual `controller.HuntNetwork`
class object and source path, with verified JointIntegrationNetwork,
FunctionalNetwork, and W2LVarnet ancestry, 1033/1027 neurons, 92 MNs, and the
expected RG/HUNT_AF cells. No class/module labels are spoofed. The already
constructed sparse wrapper is supported without adding topology. Gates lacking
active outgoing paths remain masked; their original fields are preserved.
`eligible_load_sources(net, 'Hunt_PT')` supplies the exact current-map keys,
including read-only inactive homologous input contracts so encoding remains
shared even when the historical output mask is asymmetric.

Both saved conditions pass commands 0/0.5/1/0 and 32 original `n.step` calls
after local FSA initialization, the explicit `canonical_right_role_budget`
timing profile, and coincidence installation. Source-equivalent timing controls
are preserved separately. The actual weight-update source writes `gmax` only;
its index-set union has zero overlap with incoming gate edges and writes no
gate windows. It rewrites 54/53 allocated gate-to-MN slots from unchanged,
command-independent `reflex_gain`. Direct parameter/cache checks verify their
consumption after every update. No extra per-step clamp or callback is required.
The unbound PT mechanical object is excluded from these neural-only checks.

Current Hunt/PT has ten active gates. Historical four-cell1027 has nine: its
right trunk-extensor gate remains masked while the left is active. No masked
output, missing homologue, or zero gain is restored by the coincidence adapter.
This historical asymmetry is explicitly recorded, not silently repaired.

The authoritative saved buffers were audited before and after the composition.
Only indices actually consumed by `apply_weights` are compared. Pairing uses
named R/L labels and existing source-declared midline mappings, not invented
frontal/transverse channel swaps. Geometry-dependent physical equivalence and
the unbound PT calibration remain outside that comparison.

| Saved source | PF mismatched pairs/max uS | FF mismatched pairs/max uS | Reflex mismatched pairs/max uS | Biarticular mismatched pairs/max uS |
|---|---|---|---|---|
| Current Hunt/PT1033 | 0/0 | 0/0 | 0/0 | 0/0 |
| Historical four-cell1027 | 96/2.53160 | 93/4.05124 | 132/1.79000 | 6/0.546489 |

The comparison tolerance is 1e-12 uS. Source pair counts are respectively
340/466/337/6. After the local-reflex design, reflex pair count is383; the
current source remains fully matched on these existing paired parameters,
while the historical comparator retains16 reflex mismatches (max0.0584955 uS).
Its PF/FF/biarticular mismatch counts are unchanged. Mirroring timing alone
therefore does not repair or qualify the old controller.

Both sources have144 unpaired shared-lumbar bank slots with no declared
bilateral-midline option. In current Hunt/PT, the31 literally nonzero values
are all floating-point residue (maximum3.30943e-14 uS); none exceeds the
comparison tolerance. This is not presented as meaningful asymmetric drive.
The historical source has38 unpaired values above tolerance, maximum3.84577 uS,
from the existing first-leg bank (left). No midline drive is added automatically.
Both networks already contain all four RF/HAM PF-layer cells per leg.

The 0.225 nA source-current contract remains an explicit canonical-runner
engineering reference, not a body-validated saved-condition force encoding.
These checks establish command-policy composition, not training, full control
coverage, balance, stepping, or any inherited curriculum pass.

Evidence and source-equivalent controls are under
`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_ibexc_composition\`.
The final requested timing policy is in its `canonical_right_role_budget`
subfolder, with one `HUNT_COINCIDENCE.json` per condition. Source weight files,
controller/profile code, and condition hashes are recorded. The three canonical
family regression traces remain identical to the earlier gate implementation.
