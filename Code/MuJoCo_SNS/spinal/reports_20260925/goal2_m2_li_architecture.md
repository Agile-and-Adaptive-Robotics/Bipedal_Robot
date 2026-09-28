# Goal 2 / Milestone 2 — Li's neural architecture in MuJoCo-SNS

**Date:** 2026-09-25 (EB475WS4, unattended campaign) · **Status: partial — architecture transcribed + closed loop runs end-to-end; stepping gate NOT yet passed (verdict "not yet")**
**Artifacts (all new, under `Code\MuJoCo_SNS\spinal\w2l_mujoco\`):**
- `build_li_net.py` — transcribes the `li` template (58 nodes / 98 edges) from `spinal\connectome_templates.json` into sns-toolbox, plus Li's 12 MV relay neurons and the MuJoCo muscle interface. Run standalone for a build census.
- `test_li_stepping.py` — the M2 gate: 20 s closed-loop run on the M1 body (`w2l_mjcf.xml`), measured against Li's reference; prints the VERDICT line.
- This report. Gate logs: `reports_20260925\logs\test_li_stepping_run{1..4}.log`.

Env: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe` (mujoco 2.3.7, sns-toolbox 1.5.2, CONDA_PREFIX set before the mujoco import). No pip installs. No protected file touched (new folder only; default runner/build/params untouched; no study touched). Read-only sources: `Neuromechanical_Models\Li Model\walk tester rearranged.aproj` and `...\walk new new tester added 2 axis_Standalone.asim`.

---

## 1. Li's reference numbers (HIS results, re-measured on this box)

**The named primary source contains no data.** `Li Model\Trace_2020*.txt` (all 5 files) are AnimatLab **logger byproducts, 18 lines each** (`StdXml::Load ... LogFile: F:\AnimatLab\bin\Logs\AnimatLab`) — no chart data. Verified: `python -c` count of lines per file = 18,18,18,18,18; content shown in §1 of the trace probe.

So the targets were generated from **Li's own model**: ran his `walk new new tester added 2 axis_Standalone.asim` headless on this box (method proven 2026-09-24, `reports_20260924\animatlab_w2l_airwalk.md`; this run: `AnimatSimulator.exe <asim>`, exit 0, "Simulation stopped. Time: 10.021", log `reports_20260925\logs\li_asim_run1.log`). It wrote **`Li Model\DataTool_7.txt`** (50010 rows, 0.2 ms, 10 s: contact-neuron voltages, hip-middle, height, distance — this file is Li's asim's own output and is the reference of record for this milestone).

| Li target (DataTool_7.txt, analyzed by `tmp\li_ref_analysis.py` + `tmp\li_ref_episodes.py`) | value |
|---|---|
| Gait period (L_foot stance onsets, median interval) | **1.305 s (0.77 Hz)**; R side 1.290 s |
| Stance duration | L 0.669 s / R 0.614 s → duty ≈ **0.5** |
| L/R phase | antiphase (stance-episode overlap > 20 ms: 3 of 15, and those are brief) |
| Toe contact lags heel | ≈ 0.2 s (foot rollover) |
| Pelvis height | 0.95–1.02 m the whole 10 s — **never falls** |
| Forward speed | −3.454 → +2.930 m = **+0.638 m/s** |
| Neuron firing style | contact/hip cells spike ~113 Hz in stance-bound clusters (−60 mV rest, −55 mV threshold crossings) |

Cross-check: the 0.77 Hz matches the independent 09-24 audit ("2023 original bursts ~0.77 Hz") — same walker family.

## 2. What Li's architecture actually is (mined from his .aproj this session)

Fresh read-only XML dumps (`tmp\dump_li_full3.py`, `tmp\dump_li_adapters.py`; full outputs kept in `tmp\li_full_dump.txt`). Key facts the template's `_note` summarizes, now verified against the file:

- **44 Spiking + 12 NonSpiking neurons** in an IntegrateFire module at **0.2 ms** timestep (module: rest −60 mV, threshold −55 mV (MNs −59 mV), τm 5 ms, spike strength 1, refractory 2 ms; `MaxCaConductance = 0` on every neuron — the note's "tau_Ca 7.5 s burst termination" is **not present in the aproj neurons**; burst turnover in this net is closed-loop through contact switching, not an intrinsic Ca current).
- **Tonics (TonicStimulus, nA):** 6 nA on both `L/R_CPG stance`; 5–6 nA on all 12 PF half-centers; 0 on MNs/INs/contact cells.
- **Synapse types used (EquilibriumPotential):** 'Nicotinic ACh' −10 mV (PF→MN exc), **'Depolarizing IPSP' −50 mV = excitatory despite the name** (contact, hip-II, CPG→PF, commissural-excite), 'Hyperpolarizing IPSP for CPG/IN' −70 mV (all inhibition). Per-connexion conductances 1–8 µS — these ARE the template's per-edge gains (verified per-link: heel→knee-flex-inh 8 µS, hip-swing-PF→knee-ext-inh 8 µS, hip stance PF→knee ext PF 2 µS, MN→MV 3/5/8 µS, everything else 1 µS).
- **12 "MV" NonSpiking relay neurons** (rest −100 mV) between MNs and muscles: MN −Nicotinic→ MV −identity-gain adapter→ muscle StimulusTension. The template's `mn_to_muscle` edges abut exactly these MN→MV synapses.
- **Sensory encoders** (PhysicalToNodeAdapter → ExternalCurrent): heel/toe = Sigmoid on ContactCount (A = 0.7–0.8 nA, B = 1e-8, C = 25); hip angle → `*_hip middle` = Sigmoid(A = 15 nA, B = 1e-8, C = 5, D = 0) = a near-step 0→15 nA crossing hip angle **+3.68°** (= ln(1e-8)/5).
- **Interleg coupling (the whole "Renshaw-free" gait loop):** `L_foot contact → R_CPG inhabit` (6 µS exc) and mirror; `L_CPG stance ↔ R_CPG stance` mutual inhibition (1 µS); `L_CPG stance → R_hip swing PF` excitation (1 µS) and mirror.
- Li's joints: friction `Enabled=False` on all hip/knee/ankle (0.02 coefficient, only toe joints enabled) — **his body is damped through the muscles' LinearHill B (400–800 N·s/m), which the M1 body does not model** (M1 report §3, "B — none"). This matters in §5.

## 3. Realization — `build_li_net.py` (build census, from the standalone run)

```
neurons: 56  synapses: 84  inputs: 20  outputs: 12
synapses per tag: {'pf_cross': 22, 'contact_drive': 10, 'mn_to_mv': 12,
 'rg_commissural': 4, 'pf_to_mn': 24, 'rg_laminate': 2, 'pf_drive': 8, 'ii_aff': 2}
```

= 44 template neurons + 12 MV relays; 72 template synapses + 12 MN→MV; 6 encoder ports (heel/toe/hip × L/R) + 14 tonic slots; 12 MV-voltage outputs = the M1 actuator drives. Every template edge is transcribed with its per-edge µS gain and its (tag, sign)→reversal mapping (`_reversal()`), so any edge traces back to the template/aproj. MUSCLE MAP is by geometry, not aproj label (the aproj's R-side hip muscle names are crossed vs their attachments — M1 report §4c): hip stance→back/back actuator (`hip_L_ext`/`hip_R_flx`), hip swing→front/front (`hip_L_flx`/`hip_R_ext`), knee ext→patella route, plantarflexion→heel route, dorsiflexion→toe route.

### DEVIATIONS (all loud, none silent)

1. **Graded surrogate instead of SpikingNeurons** (the ask's documented-fallback clause, w2l_cpg deviation #5 precedent). sns-toolbox 1.5.2 DOES ship `SpikingNeuron`/`SpikingSynapse` and they run on the numpy backend (probed: `tmp\check_spiking.py/.py` successors; threshold adaptation, spike reset, g_spike all execute). **Rejected on physics, not availability:** in `SNS_Numpy.forward` (site-packages `sns_toolbox/backends.py`, read this session) the spike conductance `g_spike` and its increment `g_increment` are **per-NEURON scalars added to every incoming connection row** — convergent mixed-sign inputs (every PF half-center and MN in Li's net) would mis-sum excitation/inhibition and collapse his 1-vs-8 µS per-synapse differentiation. The graded realization preserves topology, per-synapse g, reversals and tonics exactly, on AnimatLab's own −60 mV / nA / µS scale, treating presynaptic depolarization above rest as the firing-rate proxy (synapse window e_lo = rest −60, e_hi = −45). A window anchored at the −55 mV threshold was tried first (gate run 1) and is dead on this scale — a 1 µS synapse can lift a 1 µS-leak cell only to −55 mV, so nothing ever conducted; Li's cells are rate-coding around threshold, which is what the linear-rate window encodes.
2. **+12 MV relay neurons not in the 58-node template** — aproj-verified (§2); gives MuJoCo a graded actuator drive exactly where AnimatLab had one.
3. **MV→ctrl map is a calibrated soft threshold** (`MV_ON=-40 mV`, saturate at the −10 mV Nicotinic ceiling): stands in for the aproj per-muscle StimulusTension A/B/C/D curves, which were not ported (M1-scope item, same family as M1's Kse/Kpe/B exclusion).
4. **Contact encoder** = count-gated current, `CONTACT_GAIN_NA=15` (the aproj Sigmoid's effective amplitude is not reliably recoverable — AnimatLab's Sigmoid formula convention is not documented in the file; A = 0.7–0.8 nA would leave the cell ~5 mV below its own threshold). **Hip encoder** = hard sigmoid at the aproj-derived +3.68° crossing, amplitude 15 nA; `HIP_SIGN` knob flips convention if the closed loop proves mirrored.
5. **Kickoff:** +3 nA on the R tonic port for 0.2 s (Li's asim ships an inactive Stimulus_1; without a break of symmetry the start pose holds both stance CPGs equally suppressed).
6. **Settle-assist + joint damping (gate-only, runtime):** documented inside `test_li_stepping.py`; motivated by §5.

## 4. The gate — runs, numbers, verdict

Command (each): `C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_li_stepping.py [--dur=N --knob=val]`, logs `reports_20260925\logs\test_li_stepping_run{1..4}.log`.

| run | config | result |
|---|---|---|
| 1 | threshold-anchored synapse window | collapses; synapses never conduct (dead circuit) — motivation for deviation #1's window |
| 2 | linear-rate window + MV_ON map, 8 s | pelvis min −1.44 m (tunnels through the floor); knees +68..77° (through the +60° limit); 0 stance episodes |
| 3 | + 0.5 s settle-assist co-contraction, 8 s | same failure mode (min −1.44 m) |
| 4 | + runtime joint damping 1.5 N·m·s/rad as the muscle-B stand-in, 8 s | improved but still collapses: min −0.84 m, end −0.52 m; knee L +42..66°, ankle L −48..−20° (limit −20); **0 stance episodes** |

**VERDICT (runs 1–4): not yet.** The full 20 s run was not reached in a passing configuration; the 8 s iteration runs above are what was executed (each is the same gate at shorter horizon — the ask's ≥20 s run is the deliverable once a stable configuration exists).

## 5. Honest diagnosis (why it fails — measured, not guessed)

The failure is **body-side, not neural-side**. The diagnostic loop (`tmp\diag_li_loop.py`) shows the neural net alive and behaving (CPG/PF cells sit at −56.8/−56.6 mV, MV relays respond, ctrl engages) while the body crumples at landing: knee flexes to ~76° within 200 ms of the drop, ankle slams through its [−20,−5] limit, and the pelvis tunnels below the plane (contacts then never re-register — 0 stance episodes is a consequence of the fall, not the cause). Three compounding root causes, each traceable to a documented earlier deviation:

1. **Start condition:** the M1 rest pose hovers 2.6–3.1 cm (M1 report §2.8) while Li's model starts feet-on-ground at height 1.02 (his DataTool_7 shows first contact episode at t = 0.095). The M1 drop ends in the documented "crumpled kneel" before any neural drive exists. M1 open item #4 (standing-pose solve) is the real fix.
2. **No muscle damping:** AnimatLab LinearHill B = 400–800 N·s/m per muscle is what damps Li's body (his joint frictions are Enabled=False); MuJoCo 2.3.7 `<muscle>` has no damping attribute (M1 §3). Run 4's 1.5 N·m·s/rad joint damping stand-in helped (knee overshoot 76→65°, ankle −82→−48°) but is far from the ~10² N·s/m the muscles carry.
3. **Rigid tendons:** without Kse, muscle-force transients hit the undamped hinges at full rate.

## 6. Next steps (in order, for the follow-up session)

1. Give the M1 body its muscle damping: MuJoCo 3.x `<muscle damping>` if available on any lab machine, or a custom actuator plugin (the AnimatLab skill's BPA-plugin route), or per-DOF `dof_damping` sized from ΣB·r² (moment-arm-weighted muscle damping — principled, not a knob).
2. Standing-pose solve / spawn settled (M1 open item #4): start the run feet-loaded at height ~0.99 like Li's t=0.
3. Re-tune only then: `TONIC_SCALE` (cells ran ~3 mV under the naive linear-rate estimate once mutual inhibition loads them), `MU_GAIN`/`MV_ON`, `CONTACT_GAIN_NA`, `HIP_SIGN` — all exposed in `LI_KNOBS`, all reachable on the `test_li_stepping.py` command line (`--tonic-scale=…` etc.).
4. Re-run the 20 s gate; compare against §1's table (period 1.30 s, duty ≈ 0.5, antiphase, 0.64 m/s, height ≥ 0.95).

## 7. Byproducts / housekeeping

- `Li Model\DataTool_7.txt` was CREATED by the headless asim run (AnimatLab chart byproduct — same accepted category as the 2026-09-24 runs; no model file touched). It is Li's reference data and is cited by path above.
- `reports_20260925\tmp\` holds the read-only miners (`dump_li_*.py`, `li_ref_*.py`, `diag_li_loop.py`, `check_spiking*.py`) and their outputs; all repo writes this session are `w2l_mujoco\{build_li_net.py, test_li_stepping.py}` + this report + logs.

---

## 8. 2026-09-26 (easteregg2) — RETRY ON THE M3-FIXED BODY: loop still does not close; root cause is STRUCTURAL, not neural

Mandate: one bounded retry of the closed loop on `w2l_mjcf_fixed.xml` (the knee/ankle axes are now
sagittal), ≤2 h diagnosis if it fails, no retuning of Li's gains. **Verdict: still does not walk —
and the failure is now root-caused to the BODY, in two layers, with a causal A/B.**

Runs (env `D:\Anaconda\envs\myo\python.exe`, CONDA_PREFIX=D:\Anaconda\envs\myo; net census on this
machine reproduced §3 exactly: 56/84/20/12):
1. 3 s smoke: `test_li_stepping.py --dur=3` → 0 stance episodes.
2. **20 s gate** (log `logs\test_li_stepping_fixed_run1_20s.log`): 0 stance episodes, duty 0.00,
   heel flags ON 0.0% (L) / 0.2% (R). VERDICT "not yet".
3. 20 s forensics probe `tmp\diag_li_fixed.py` (log `logs\diag_li_fixed_run1.log`).
4. **Welded-vs-freejoint A/B** `tmp\diag_li_freejoint_ab.py` (log `logs\diag_li_freejoint_ab.log`;
   probe model `tmp\li_free_probe.xml` = shipped xml + ONLY `<freejoint name="root"/>`).

### Finding 1 — the M2 "ground body" has a WELDED pelvis: the loop is open by construction

`w2l_mjcf_fixed.xml` (and `w2l_mjcf.xml`) ship the Root body with **no joint at all** (mujoco census:
`njnt = 8`, all leg hinges; Root `body_jntnum = 0`, `body_pos = [-3.454, 0, 0.99298]`). A jointless
body is bolted to the world: the pelvis sat at exactly z = 0.993 m for all 20 s in probe A while the
feet hover 2.6–3.1 cm up (M1 §2.8) and **no leg pose can reach the floor** — the same property the
air gate exploits by lifting the welded Root +0.30 m (`test_w2l_air.py:17,90`). Li's CPG is
**contact-driven** — its stance trigger is heel contact — so with heels that never touch, the SNs
never fire (L_foot ON 0.000 of the run) and no stance episode can exist. This contradicts the
generator's own deviation note, `make_w2l_mjcf.py:49-50`: "Root Freeze=True in the aproj is IGNORED
(the Freeze trap: pelvis must be FREE in MuJoCo)" — the freejoint was never emitted, and
`validate_body.py:52` then enshrined `njnt == 8` as the M1 gate.

### Finding 2 — the gate's height/speed metrics were never pelvis height: they are joint angles

`test_li_stepping.py:147` logs `d.qpos[2], d.qpos[0]` as "pelvis z, x". With a welded root, qpos
holds only the 8 hinge angles (`qposadr` census: 0=hip_L, 1=knee_L, **2=ankle_L**, …) — so every
"pelvis height" printed since run 2 was **ankle_L in radians**: the old "tunnels to −1.44 m" was
ankle_L = −82.5° (the −82° limit blow-through already noted in the code comment), today's
"min −0.639, end 0.988" is ankle_L −36.6°..+56.6°. The body never fell and never translated; the
"speed vs Li 0.64 m/s" line is likewise hip_L drift and can never pass on a welded root. (Fix when
resumed: address the pelvis through a freejoint's qposadr, as probe B does.)

### Finding 3 — with the freejoint added (A/B probe B), the body crumples in <0.5 s: the freejoint is necessary but NOT sufficient

Probe B (identical loop, defaults, zero gain changes): pelvis falls to z = 0.100 m and the 31 kg
Root box RESTS ON THE GROUND (`ground × Root` = 6.6e6 N·steps, the largest contact of the run);
knee_L −44° by t = 0.5 s; heels still never register (L_foot ON 0.000). This is M1 report §"gate
behavior" (passive drop → crumpled kneel) and M2 §5's root causes, unchanged by the axis fix:
- **no muscle damping** (AnimatLab LinearHill B 400–800 N·s/m is what damps Li's body; MuJoCo
  2.3.7 `<muscle>` has no damping; runtime stand-in is 1.5 N·m·s/rad per DOF — measured ~2 orders
  short), **rigid tendons**, and
- **the spawn pose sits OUTSIDE its own limits**: `qpos0 = 0` for every joint vs ankle range
  [−20,−5]° (aproj limits verified via `w2l_source_dump.json`) → 5° outside at t = 0, first breach
  at t = 0.00 s; kN-scale muscles (ankle actuators measured to 1.8–1.9 kN) then blow the soft
  limits through **18.7–99.1°** (ankle 90.9° welded / 99.1° free), the femur boxes grind on each
  other (`femur_L × femur_R` 3.8e6 N·steps, shins cross 0.3–0.6 s). Open sub-question (not
  resolvable from the files here): whether the aproj ankle's joint-zero coincides with the
  assembled rest pose — same family as the knee-range sign M3 had to reconcile by hand
  (`fix_joint_axes.py` "knee RANGE sign reconciliation").

Neural side checked and ALIVE in all runs: MV relays rest at −99 mV (below MV_ON, zero drive) and
reach −23…−35 mV once tonics/hip-II load them, with ctrl > 0 on knee-ext/ankle-flx/hip pools — the
net drives; it just never receives the heel trigger and its body cannot use what it drives.

### What was changed / NOT changed
- Changed: `test_li_stepping.py` now points `MJCF` at `w2l_mjcf_fixed.xml` (the retry's vehicle;
  loud comment in place). New byproducts: `logs\test_li_stepping_fixed_run1_20s.log`,
  `logs\diag_li_fixed_run1.log`, `logs\diag_li_freejoint_ab.log`, `tmp\diag_li_fixed.py`,
  `tmp\diag_li_freejoint_ab.py`, `tmp\li_free_probe.xml` (temp model, not the model of record).
- NOT changed: Li's gains/tonics/encoders/MV map (all verbatim defaults), `w2l_mjcf_fixed.xml`,
  `w2l_mjcf.xml`, `make_w2l_mjcf.py`, `validate_body.py`, any protected path. The W2L air gate was
  not re-run here (taken as given from the ask).

### Order of work when resumed (Ben's call; not started)
1. Emit `<freejoint>` for the Root in `make_w2l_mjcf.py` (its own comment already requires it) and
   fix `validate_body.py`'s `njnt == 8` gate + `test_li_stepping.py`'s qpos addressing together.
2. M1 open items #2/#4 become hard prerequisites: spawn settled (feet loaded, ~0.99 m like Li's
   t = 0) and a muscle-damping stand-in sized from ΣB·r² (moment-arm-weighted), NOT a scalar knob.
3. Reconcile the ankle joint-zero/range at spawn (qpos0 vs [-20,-5]°).
4. Only then re-attempt the 20 s gate against §1's table.
