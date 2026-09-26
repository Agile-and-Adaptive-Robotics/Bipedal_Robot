# Goal 4 — VARIANT 2: the 6-synergy walker network (`syn6`)

Campaign 2026-09-25, EB475WS4. Status: **PARTIAL** — gates (a) and (b) PASS,
gate (c) PASS on the network-only air rhythm, PARTIAL on the full-runner air
run (marginal alternation; two E-latch mechanisms found and fixed, one
residual documented). All numbers below were read back from commands run in
this session; every run command is quoted.

## 1. Summary

`build_network_syn6.py` implements the 6-synergy layout: per side ONE rhythm
generator (the stock `_build_rg` pattern: persistent-Na half-centers +
InE/InF lamination) driving SIX synergy PF layers (`PF_S1..PF_S6_{r,l}`),
one per channel of `synergy_basis.npz` (W_r/W_l [43 x 6], name-keyed).
PF->MN conductances are the measured W mapped through the audited
Szczecinski-2017 Eq-18 function imported from `fsa_backsolve.py`
(`analytical_conductance`, fsa_backsolve.py:132-137). The rule-file dress
(heel stance reset at the PF layer, toe DF-inhibition, Ib load, Shevtsova
commissurals, Shinohara autogenic afferents) carries Ben's connectome JSON
gains. Selection is default-inert: env `AARL_NET=syn6` or `G["syn6"] > 0`
(new params.py key, default 0.0).

Gate results:

| gate | result | evidence |
|---|---|---|
| (a) defaults | **PASS** | `gate_syn6.py a` → `default build: neurons=410 inputs=376 synapses=1186` (expected 410/376/1186 exact) |
| (b) variant | **PASS** | `gate_syn6.py b` → `syn6 build: neurons=794 inputs=382 synapses=3033`, port coverage 382/382, W coverage 86/86 non-pruned, Eq-18 spot check `k=0.5: builder g=0.454545 == closed form` |
| (c) rhythm, network-only | **PASS** | `gate_syn6.py c1` → 12 s, DRIVE=2.5 nA: RG-E 5 peaks (period 2.133 s), E-F corr **−0.989**, amplitude 8.97 mV, all 6 PF channels burst (4.67–5.65 mV); S5 distinct (see §4) |
| (c) rhythm, full runner air | **PARTIAL** | 22 s `runner.py --no-ground` under `AARL_NET=syn6`: E/F events strictly interleaved `FEFEFFE` (peaks E: 5.60/6.00/11.77 s, F: 5.29/5.77/6.15/11.72 s), knee flexes to −19.9°, but the rhythm is slow/irregular (~2–3 s events, then quiet); the gate's ≥5-peak criterion FAILed |

## 2. Files

**New** (all in `Code\MuJoCo_SNS\spinal\`):
- `build_network_syn6.py` — the variant builder (`Syn6Network` subclasses the
  stock `SpinalNetwork`; `build()` drop-in with the same signature/contract).
- `reports_20260925\tmp\gate_syn6.py` — the three gate scripts (sections
  a / b / c1 / c2).
- `reports_20260925\tmp\{inspect_basis,probe_syncount,probe_phase,analyze_c1,analyze_c2,c2_peaks}.py`
  — one-off evidence probes (kept under tmp for this campaign's audit trail).
- `reports_20260925\logs\syn6_net_smoke.npy` (C1 trace log),
  `reports_20260925\logs\syn6_air_smoke.npz` + `syn6_air_smoke.log` (C2 run).

**Modified (default-inert; gate (a) re-run after every edit):**
- `build_network.py` `build()` (build_network.py:1018-1034): variant selector
  block — co-edited with the w2lvar agent's `AARL_NET=w2lvar` branch; added
  `AARL_NET=syn6 or G.get("syn6", 0.0) > 0` → lazy import of
  `build_network_syn6`. Unset/0 → stock path untouched.
- `params.py` G dict: two new default-0 keys, `syn6=0.0` (selector) and
  `syn6_brainstem=0.0` (see §4, E-latch finding 1). JSON RULE satisfied:
  both are params.G keys defaulting 0; any winner json that records full G
  carries them; the existing `--bestN` loaders ignore unknown keys, and
  `ib_rge`/`heel_rge`/`toe_rge` (which gate variant dress edges, §3) are
  long-established keys already saved by v11+ studies.
- `runner.py` `_pf_watch` selection (runner.py:1145-1153): syn6 branch first
  (`PF_S1_r..PF_S4_r` when `G["syn6"]>0` or `AARL_NET=syn6`), then the
  stock joint_pf/phase branches unchanged. Default behavior identical
  (gate (a) exercises the stock build; the logged npz `neuro_names` confirm
  the syn6 branch only under the env).

## 3. Architecture

- **RG**: stock `_build_rg` verbatim pattern — NaP bursters (params.NAP,
  fixed tau_h via `SNS_NumpyFixedTau` inherited from the parent `compile`),
  IN-laminated mutual inhibition (`G["rg_mutual_inh"]`), DRIVE/POSTURE
  edges, optional `rg_weak_exc`.
- **Six PF layers**: `PF_S{k}_{side}`, tau = `TAU["pf"]` × a per-channel
  multiplier. Families from the measured phase profiles
  (`fsa_results/fsa_backsolve.npz {side}_pf_phase_mean`, heel=0/toe-off=50):
  stance-dominant → E (driven by RG_E), swing-dominant → F. **Symmetrized**
  across sides (mean stance fraction) — the measured per-side fractions
  disagree only on S5 (0.581 r vs 0.239 l), exactly the S5/S6 bilateral
  instability flagged in synergy_model.py:435; symmetrized families =
  `[F,E,E,F,F,E]` on BOTH sides (3E/3F). Within a family, channels ranked
  by measured peak phase get tau multipliers interpolated over
  (E: 0.6→1.2, F: 0.9→1.6) — e.g. side r E: S6 0.6 / S3 0.9 / S2 1.2.
  Measured inputs (gate b output): stance_frac=[0.196, 0.953, 0.922, 0.246,
  0.41, 0.89], peak_phase r=[95.0, 40.5, 10.0, 80.5, 26.0, 3.5] %.
- **PF lamination**: stock topology (E cells → PF_IN_E → inh F cells and
  back, `G["pf_recip_inh"]`) over the COMMITTED channels only; the genuinely
  MIXED channel (stance fraction inside `SYN6_MIXED_BAND` = 0.35–0.65; here
  S5 at 0.41) is excluded from lamination and receives BOTH RG drives
  (`G["rg_to_pf"]` × stance fraction and × complement) — see §4, degeneracy.
- **PF → MN (the data)**: for synergy k, muscle act:
  `g = Eq18(W[act,k]) = k·R·Gm/(ΔE − k·R)`, R = E_HI = 5 mV, Gm = 1 µS,
  ΔE = E_REV_EXC = 8 mV — imported from `fsa_backsolve.analytical_conductance`
  (not re-implemented). Edges only where W > 0 (conditional topology); all
  86×6 W entries have k < ΔE/R = 1.6 so none is Eq-18-invalid. Measured:
  152 edges (l) / 159 (r), g ∈ [0.0007, 3.3833] µS (l) / [0.0001, 2.6944] (r).
  The 6 runner-pruned actuators (`runner.PRUNE_MUSCLES`, set-equality
  asserted in gate b) keep their MN + reflex arc (the runner indexes
  `mn_names` for all 92) but receive NO synergy drive — their names are
  absent from W, mirroring the runner's force-zeroing of the same list.
- **Rule dress** (gains from Ben's JSONs, verbatim unless noted):
  heel → InE 0.5 + InF 0.5 + PF_IN_E 0.5 (stance reset AT the PF layer);
  toe → TOEDF 5.0 → (inh 2.749) the dorsiflexion channel (S5, = the
  F channel with the largest ankle_df W mass, both sides); per stance
  muscle Ib → IBEXC 0.5 → MN 0.5 (stance reversal); Shevtsova
  V2a/V0V/V0D/V3-E per side (RG-F→V2a 1.0, RG-F→V0D 0.7, RG-E→V3-E 0.35,
  V2a→V0V 1.0, V0V→contra InE 0.6, V0D→contra RG-F 0.07, V3-E→contra RG-E
  0.02, V3-E→contra InE 1.0); Shinohara autogenic motif per muscle
  (Ia→MN 2.0, Ia→IaIN 1.0 → antagonist MN 0.5 with IaIN↔IaIN 0.5;
  II→IIX 1.0 → MN 0.5; II→IIIN 0.5 → antagonist MN 0.5; Ib→IBIN 1.0 →
  MN 0.5 with IBIN↔IBIN 0.5) — proprioception autogenic per muscle,
  never group-broadcast (Ben's reading rule).
- **Documented deviations**: D1 — RG→PF and PF lamination use the
  S3K-tuned knobs (`G["rg_to_pf"]` 2.4 / `G["pf_recip_inh"]` 4.0), not the
  drawing's W2L-spiking-scale 0.1/2.749; Eq-18 presumes the PF source
  reaches E_HI (fsa analytic_pred peak 1.27 confirms that regime).
  D2 — Shevtsova IniE/IniF/Ini folded into the existing laminated InE/InF;
  the file's asymmetric V3-E→contra-RG-E gains (0.02 vs 0.5) taken as the
  "small exc" 0.02 symmetric. D3 (revised after gate c2) — the runner's
  single DRIVE port already carries the brainstem gamma/alpha role via
  `descend_to_rg_e/f`; the folded gamma/alpha cells are built only when
  `G["syn6_brainstem"] > 0` (default 0; see §4).
- **Runner compatibility**: input ports = every port the runner writes
  (382/382 checked): 8 shared + HEEL_c/TOE_c/LOAD_c per side (stance_fb
  True) + POST/Ia/II/Ib per actuator. `aff_loops`/`vest` False; no
  Renshaw/F1-KINH/AFF loops (not in this variant's dress list). The heel/
  toe/LOAD ports are silent until the existing gain keys `heel_rge`/
  `toe_rge`/`ib_rge` are set > 0 by a study (runner-side scaling), and the
  variant's Ib-load reflex edges are gated on `G["ib_rge"] > 0` (§4,
  E-latch finding 2).

## 4. Findings made during the gates (each fixed or documented)

1. **Brainstem fold E-latches the full net (FIXED, gated).** Gate-c2
   attempt 1 (`AARL_NET=syn6 runner --no-ground`, log
   `logs\syn6_air_smoke.log`): PF_S2/S3 saturated tonically (mean 5.03 mV),
   F channels pinned (−1.34 mV), RG swing 0.43 mV. Cause: folding the
   Shevtsova brainstem gamma/alpha onto the single DRIVE port adds a fixed
   +0.5-gain E/F bias on TOP of `descend_to_rg_e/f` at the runner's
   `walk_drive = 4.0` (runner.py:710) plus POSTURE→RG_E — over-drive tips
   the NaP E half-center into never-release. Fix: `syn6_brainstem` key,
   default 0 (cells + all their edges built only when > 0).
2. **Un-gated Ib-load reflex E-latches the full net (FIXED, gated).**
   Attempt 2 still latched after ~2 cycles. Cause: per-muscle
   Ib→LBIN(0.5)→RG_E/InE/PF-E is FORCE-proportional; in air the standing
   muscle tone already produces real Ib encoder drive (runner writes
   `u[Ib_a] = 0.5·g_ib·force_norm` unconditionally, runner.py:1283), so the
   RG gets a tonic E bias the stock stack does not have (its LBIN pathway
   is gated by `G["ib_rge"]=0`). Fix: the variant's Ib→LBIN and
   LBIN→{RG_E,InE,E-channels} edges are built only when `G["ib_rge"] > 0`
   (existing default-0 key; the curriculum can turn the load reflex on).
3. **Synergy-layer degeneracy (MITIGATED, residual documented).** Gate-c1
   first measurement: all E-family channels collapsed to one waveform
   (pairwise corr +0.98…+0.997, common peak phase ~55.8 % cycle) and all
   F-family to another (~2.2 %) — the tau stagger alone is too weak vs a
   multi-second cycle. Mitigation 1 (blind soft RG drive split) DEADLOCKED
   the F family (faint opposite-family tone held the antagonistic PF
   lamination on — the same co-contraction trap as the stock net's old
   DRIVE→PF tonic term, params.py:88-91 note). Mitigation 2 (soft split
   for mixed channels only) flipped the latch: mixed S5 rode RG_E through
   the E phase and held PF_IN_F on, pinning S2/S3/S6 at −0.11 mV. Final
   form: mixed channels get BOTH drives and are EXCLUDED from lamination;
   result (final c1): S5 peak phase 78.4 % vs 2.0 % for its family
   siblings, |corr(S5, any)| ≤ 0.51 — genuinely distinct. RESIDUAL: the
   committed within-family channels remain near-degenerate in constant-
   drive air (S2/S3/S6 +0.98…+0.997; S1/S4 +0.997). Their MN-side
   projections still differ (distinct W columns), but per-channel delay
   interneurons would be the next lever if Ben wants 6 visually distinct
   traces; not built in this pass.
4. **Left-right synchrony (documented, not a gate criterion).** The
   Shevtsova set is near-symmetric between antiphase (V0V/V0D) and
   synchrony (V3-E) systems; measured RG-E peaks are left/right
   synchronous (identical peak times in c1 and c2). Ground gait will need
   the antiphase side to win — a curriculum tuning surface (and the same
   issue the stock stack saw with `v3_gain`, params.py:224-226).
5. **fsa npz key trap (banked).** The phase grid is saved as
   `{side}_phase_grid` (no "pf"; fsa_backsolve.py save block) while the
   profiles are `{side}_pf_phase_mean` — the builder asserts and would
   otherwise silently fall back to W-mass family assignment (that fallback
   is what gate-b's first run actually measured: peak_phase=[0,10,20,30,
   40,50] was the synthetic fallback, caught because it was too clean).

## 5. Gate commands and raw evidence

All runs: cwd `D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal`, interpreter
`C:\Users\Ben Bolen\.conda\envs\myo\python.exe`.

```
python -m py_compile build_network_syn6.py build_network.py params.py runner.py   → COMPILE_OK
python reports_20260925/tmp/gate_syn6.py a    → default 410/376/1186 PASS (re-run 3x, PASS each time)
python reports_20260925/tmp/gate_syn6.py b    → syn6 794/382/3033 PASS (see §1/§3 numbers)
python reports_20260925/tmp/gate_syn6.py c1   → 12 s rhythm PASS (period 2.133 s, corr −0.989)
python reports_20260925/tmp/gate_syn6.py c2   → runner --no-ground exit 0, npz written;
        RG peaks E [5.6, 6.0, 11.77], F [5.29, 5.77, 6.15, 11.72], interleave FEFEFFE,
        knee_angle_r −19.9..+10.0 deg; ≥5-peak criterion FAIL ⇒ PARTIAL
python reports_20260925/tmp/analyze_c1.py     → full 6x6 corr matrices (S5 mitigation evidence)
python reports_20260925/tmp/c2_peaks.py       → peak/interleave extraction quoted in §1
```

Artifacts: `logs\syn6_net_smoke.npy`, `logs\syn6_air_smoke.npz`,
`logs\syn6_air_smoke.log` (22 s runner log). The three c2 attempts used the
same npz path (overwritten); the numbers quoted are from the FINAL attempt
(ib_rge-gated build) except where labeled as attempt 1/2 in §4.

## 6. Not done / next steps (for the campaign owner)

- Curriculum chain for this variant NOT started (out of the ~110-min
  budget): the five fresh study names reserved are `curr_syn6_s1_air_deaff`,
  `curr_syn6_s2_air_aff`, `curr_syn6_s3_balance`, `curr_syn6_s4_walk_nocontact`,
  `curr_syn6_s5_walk_contact`. Hooks are ready: `G["syn6"]=1` (or env) in
  the variant's copy of `_curriculum.py`; that copy must use its OWN sqlite
  db (the stock db path is hardcoded) and must re-verify npz column
  identity (variant neuro columns: RG unchanged, PF watch = PF_S1..S4_r).
- The full-runner rhythm is slow (~2.6–3 s events vs the 1.23 s reference
  cadence) and dies after ~2 cycles + one late event — tuning surface:
  DRIVE level, `rg_weak_exc`, `rg_nap_h`, POSTURE→RG_E, plus turning on
  `ib_rge`/`heel_rge`/`toe_rge` for ground stages. This is stage-tuning
  work, not builder work; the s3k stack went through the same (v4-era
  "E-duty is architectural" notes apply).
- If Ben wants 6 visibly distinct layer traces, add per-channel delay INs
  (§4.3 residual); if he wants the Shevtsova brainstem in the loop, give
  the runner a second descending port instead of the DRIVE fold (§4.1).
- w2lvar cross-check: both variants now share the selector block in
  build_network.py — a combined `AARL_NET` value is rejected by fallthrough
  (only exact `w2lvar`/`syn6` match); no interference observed (gate a ran
  with the w2lvar branch present).
