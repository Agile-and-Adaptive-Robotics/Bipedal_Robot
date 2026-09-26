# goal4 — VARIANT 1 `w2lvar`: Walker_2_Layer_CPG layout on gait2392 (build + gates)

**Status: COMPLETE (all three gates PASS).** 2026-09-25 12:55–14:15 EB475WS4,
env `myo` (`C:\Users\Ben Bolen\.conda\envs\myo\python.exe`), cwd
`Code\MuJoCo_SNS\spinal`.

## What was built

| file | role |
|---|---|
| `spinal\build_network_w2lvar.py` | NEW — the variant builder (`W2LVarnet`, `build(...)`), all variant gains in a module-local `GAINS` dict |
| `spinal\build_network.py` | ONE minimal default-inert edit: a variant selector at the top of `build()` (lines 1018–1031; shared with the sibling `syn6` variant — merged, and gate A re-ran PASS after the merge) |
| `spinal\reports_20260925\tmp\gate_ab_w2lvar.py` | gates (a)+(b) script |
| `spinal\reports_20260925\tmp\run_gate_c.bat` | gate (c) driver (`AARL_NET=w2lvar`, `AARL_NPZ=w2lvar_smoke.npz`, `runner.py --no-ground --drive 2.5`) |
| `spinal\reports_20260925\tmp\analyze_smoke.py` | npz rhythm-metric analysis |
| `spinal\reports_20260925\logs\gate_c_smoke{,2,3}.log` | run logs (1 = co-latch baseline, 2 = per-muscle-divided afference, 3 = final PASS) |
| `spinal\reports_20260925\w2lvar_smoke2.npz`, `w2lvar_smoke3.npz` | smoke traces (neuro columns incl. the watch names) |

Selection: **env `AARL_NET=w2lvar`** (lazy import inside `build()` — with the
var unset the standard path never touches variant code). The variant applies a
documented `VARIANT_G` overlay to `params.G` (heel_rge 0.5, toe_rge 0.5,
joint_pf 1.0, ia_in 0.5, renshaw 0.5, ib_rge 0); this is required because the
RUNNER itself scales the `HEEL_c`/`TOE_c` port currents by `G["heel_rge"]/
["toe_rge"]` (runner.py:1387–1388) and picks the neuro watch columns by
`G["joint_pf"]` (runner.py:1145). A variant run is therefore a dedicated
process; `optuna_walk.db` and every default-path behavior untouched.

## Architecture (words)

Per side, mirroring Walker_2_Layer_CPG on the 92-actuator body:

1. **One RG** per side: persistent-Na half-centers `RG_E/RG_F`
   (`NonSpikingNeuronWithPersistentSodiumChannel`, params.NAP, **fixed
   tau_h = 0.35 s** via the inherited `SNS_NumpyFixedTau` swap) + `InE/InF`
   laminated mutual inhibition (g = `rg_mutual_inh` 4.0, no direct RG↔RG
   inhibitory synapses) + DRIVE→E/F (1.7/1.4) + POSTURE→E (0.8). This is a
   verbatim copy of the parent `_build_rg` mechanics.
2. **Two PF pairs** per side: `PF_HIP-E/F` (hip) and `PF_KNEE-E/F` (ONE
   knee+ankle **synergy** pair). RG-E drives both E cells, RG-F both F cells
   (g 2.4). Each pair carries its own IN lamination
   (`PF_IN_{HIP,KNEE}-{E,F}`, 4 INs/side) at the drawing's **2.749**.
   The synergy cells drive per-muscle MNs (all 92 stay per-muscle):
   KNEE-E → knee_ext + ankle_pf, KNEE-F → knee_flex + ankle_df — the W2L
   biarticular-gas pattern; PF→MN weights come from the parent's fitted
   table `joint_pf_weights.json` (merged cell delivers the KNEE-E row to
   knee muscles and the ANK-E row to ankle muscles), anatomical-1.0
   fallback via `W2L_GROUP2HC`/`W2L_CROSS`.
3. **Per-muscle motif set** (92×): MN + Ia/II/Ib encoders (runner-fed ports
   `POST_/Ia_/II_/Ib_<act>`), IaIN (reciprocal), IIX (II exc), IIIN
   (II inh), IBIN (Ib autogenic inh), RC (Renshaw), group IBEXC (stance Ib
   reversal, RG-E-gated), KINH available at gain 0.
4. **Commissurals** (replaces the parent C1/V3 block): the Shevtsova set
   V2a→V0V→contra-INI→inh contra-RG-F, V0D→inh contra-RG-F,
   RG-E→V3-E→(small exc contra RG-E, exc contra InE), DRIVE→inh V0V/V0D
   (brainstem alpha), plus the master-rules crossed heel edge
   heel→contra-InF exc.

## Rule-file → edge mapping

Master rules `ben_rules_20260924.json` (90n/99e):

| rule-file edge / group | w2lvar realization | gain |
|---|---|---|
| `HC-RG-E/F ↔ InE/InF` lamination (0.5) | RG→InE exc, InE→contra-F inh, both directions | 4.0 (`rg_mutual_inh`, stack calibration) |
| `RG-E → HC-PF-E_*` (0.1) | RG→PF drive ×2 cells | 2.4 — see deviations D1 |
| `HC-PF-E_x → IN-PF_x → HC-PF-F_x` (2.749) | per-pair IN lamination, both directions | 2.749 literal |
| `heel → IN-InE_2806` (ipsi, 0.5) | heel IN → InE exc | 0.5 |
| `heel → IN-InF_2813` (contra, 0.5) | heel IN → contra InF exc (wired post-build) | 0.5 |
| `heel → IN-PF_2810/2812/2814` (0.5) | heel IN → `PF_IN_HIP-E` + `PF_IN_KNEE-E` exc (2 pairs vs 3 stacks) | 0.5 |
| `toe → IN-PF_dorsiflexion_inhibit (5) → HC-PF-Dorsiflexion (2.749)` | toe IN → TOEDF exc → inh `PF_KNEE-F` (the dorsiflexion-driving HC) | 5.0 / 2.749 |
| `Ib grp → RG-E / InE / PF-E` (0.5) | per-muscle autogenic: each stance-extensor Ib → RG-E, InE, PF_HIP-E, PF_KNEE-E | 0.5 ÷ 26 per muscle (= aggregate, D2) |
| `ia_homo` (2) / `ia_recip` (0.5 via Ia_A) | Ia→MN homo 0.6; Ia→IaIN 0.6, IaIN→antagonist MN inh 0.4, IaIN↔IaIN 0.5, PF-F gate 0.5, RC→IaIN inh | stack-proven calibration |
| `ii_exc` / `ii_inh` (0.5) | II→IIX→same MN exc 0.4; II→IIIN→antagonist MN inh 0.4 | stack calibration |
| `ib_auto` (0.5) | Ib→IBIN→same MN inh 0.35; IBIN↔IBIN 0.5 | stack calibration |
| `ib_rev` (0.5) | Ib→IBEXC_grp 0.5 (RG-E gate 1.0) → MN exc 0.6 | stack calibration |
| `rc` motif | MN→RC 1.0, RC→MN 0.5, RC↔RC 0.5 | G["renshaw"]=0.5 |
| `cross: RG-F→c1→contra RG-F` (0.5) | replaced by the Shevtsova set (below), per task instruction | — |

Shinohara `ben_shinohara_20260924.json`: symbolic `flexor ii/ia` /
`extensor Ib` expanded **autogenically per muscle** (never group-broadcast);
the duplicated gains (0.0007/0.53, 0.000586/0.583) are honored as the Ia-vs-II
subpopulation distinction via the stack's separate Ia/II calibrations;
flexor afferents → RG-F + InF + PF_F trio, extensor Ib → RG-E + InE + PF_E
trio.

Shevtsova `ben_shevtsova_20260924.json` (json gain × 4.0 calibration,
w2l_cpg convention): RG-F→V2a 4.0; V2a→V0V 4.0; V0V→contra INI 2.4;
INI→contra RG-F inh 2.0; RG-F→V0D 2.8; V0D→contra RG-F inh 1.2; RG-E→V3-E
1.4; V3-E→contra RG-E 0.04; V3-E→contra InE 1.0; DRIVE→inh V0V/V0D 2.0.
(Ben's drawing heel/c1 rows above take precedence where they overlap.)

## Gate evidence

**(a) Defaults gate — PASS.** Command:
`"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" reports_20260925\tmp\gate_ab_w2lvar.py`
(AARL_NET unset), run TWICE (before and after the sibling syn6 merge):

```
counts neurons/inputs/synapses = (410, 376, 1186)
GATE A: PASS
```

Protected `spinal_run.npz` untouched throughout (mtime 01:10 AM, checked
after every runner invocation; the smoke wrote `w2lvar_smoke.npz` via
`AARL_NPZ` — the runner's "saved spinal_run.npz" line is hardcoded text).

**(b) Variant gate — PASS.** Same script with `AARL_NET=w2lvar`:

```
variant counts neurons/inputs/synapses = (888, 382, 7430)
delta vs default = (478, 6, 6244)
spot-check names present: 23/23
parent CIN block absent: True
watch aliases -> merged cell: True
compiled backend class: build_network.SNS_NumpyFixedTau
GATE B: PASS
```

(+6 inputs = HEEL_c/TOE_c/LOAD_c × 2 sides. Watch aliases: runner logs
`PF_ANK-F_{r,l}` → the merged `PF_KNEE-F` cell; alias bookkeeping lives in
`idx` only, the compile assert is alias-aware.)

**(c) Rhythm smoke — PASS (after 2 tuning rounds).** Command:
`reports_20260925\tmp\run_gate_c.bat > reports_20260925\logs\gate_c_smokeN.log`
= `runner.py --no-ground --drive 2.5`, 22 s air run (10 s analyzed walk
window ≥ the required 10 s), analysis
`tmp\analyze_smoke.py reports_20260925\w2lvar_smoke3.npz`:

```
RG_E_r/RG_F_r  r = -0.910   bursts 6 / 5     (ipsilateral antiphase)
RG_E_l/RG_F_l  r = -0.916   bursts 5 / 6
RG_E_r/RG_E_l  r = -0.370   RG_F_r/RG_F_l r = -0.726   (L/R alternation)
PF_HIP-E_r/PF_KNEE-F_r r = -0.976             (hip-E vs knee-ankle-F)
runner summary: RG_r 6 cycles, period 1.80 s (0.55 Hz), E-duty 0.69;
                RG_l 5 cycles, 1.80 s, E-duty 0.71; stayed up
joint sweep (r): hip -34..+64, knee -39..+13, ankle -82..+51 deg
```

Tuning history (both fixes are in the module, documented inline):
1. smoke1 (`gate_c_smoke.log`): **co-latch** — RG_E AND RG_F both pinned at
   ~3.0 mV, 0 bursts. Cause: the drawing's afferent-node gains (Ib→E trio
   0.5, flexor→F trio 0.53) were broadcast from EVERY per-muscle afferent
   (26 stance / 19 flexor muscles per side), ~15× the intended aggregate.
   Fix D2: per-muscle edges carry the aggregate ÷ family count.
2. smoke2 (`gate_c_smoke2.log`): rhythm on but L/R in-phase
   (r(RG_E_r,RG_E_l)=+0.86) — the V3-E sync leg (4.0) overpowered the
   crossed F inhibition (0.3/0.28). Fix D3: antiphase edges raised
   (INI 2.0, V0D 1.2), sync leg reduced (V3E→InE 1.0); all four commissural
   cell classes kept.
3. smoke3 (`gate_c_smoke3.log`): PASS numbers above.

## Deviations from the rule files (all deliberate, none silent)

- **D1 — RG→PF drive 2.4, not the drawing's 0.1.** 0.1 is the SynAmp-scale
  value of the contact-driven AnimatLab reference where PF-E additionally
  receives 0.5 of live load afference; on the 0..5 mV toolbox scale the
  identical template edge is calibrated to `G["rg_to_pf"]`=2.4
  (`w2l_cpg\README.md` pf_drive row). A literal 0.1 leaves PF cells
  sub-threshold (V_ss ≈ 0.1 mV).
- **D2 — per-muscle afferent central gains = aggregate ÷ family count**
  (Ib 0.5/26, Ia/II 0.53/19). The rule files' SN-Ia/SN-II/extensor-Ib are
  ONE symbolic aggregate node per family; the autogenic expansion
  duplicates the edge per muscle, so without the division the sum over the
  family is ~26× the drawing (measured: co-latch, smoke1).
- **D3 — Shevtsova antiphase/sync rebalance** (INI-inh 0.075→2.0, V0D-inh
  0.07→1.2, V3E→contra-InE 1.0→1.0×… reduced from 4.0 to 1.0; contra-RG-E
  0.02→0.04): at the raw json ratios the V3-E sync leg wins and both legs
  step together (smoke2). The set is structurally complete; the weights are
  the curriculum's to refine.
- **D4 — toe→dorsiflexion inhibition lands on the merged KNEE-F cell**,
  so knee-flex drive dips with it. Inherent to the ONE-knee+ankle-synergy
  layout this variant was asked to build (Ben's drawing has a separate
  ankle stack; the 3-stack layout is where that edge is exact).
- **D5 — heel→RG direct edges (legacy parent branch) not present**; heel
  reaches the RG only through the laminated InE/InF + PF-layer INs, per the
  master drawing. LBIN/LOAD_c port exists for the runner contract but is
  inert (ib_rge 0; the drawing's Ib→RG-E is per-muscle, D2).
- **D6 — trunk groups get no PF drive** (parity with the parent joint-layer
  build; trunk rides POSTURE/BAL_TRK, and the W2L walker has no trunk).

## Known cosmetic issue

`runner.py`'s final matplotlib panel hardcodes the phase-cell names
(`PF_E2_r`, runner.py:1787 region) → "(plot skipped: 'PF_E2_r')" at the end
of a variant run; the npz itself carries all watch columns. Harmless; a
variant-aware plot list is a later fix.

## Not done (next stages)

- No ground walking, no curriculum runs, no optuna studies — the reserved
  fresh names `curr_w2lvar_s1_air_deaff … s5_walk_contact` were NOT touched
  (this task = build + gates only). If/when the variant enters studies, the
  module-local `GAINS` keys must move into the search space per the JSON
  RULE (they are documented knobs today, not `params.G` keys — a deliberate
  choice to keep `params.py` unedited).
- `_curriculum.py` is NOT parameterized for variants yet (its db path is
  hardcoded `sqlite:///optuna_walk.db`; npz redirect via `AARL_NPZ` works).
- The synergy loses knee-vs-ankle phase differentiation by design (one
  pair drives both); if Ben wants it back, the 3-stack layout is variant-2
  territory (`syn6` sibling slot remains per the campaign plan).
