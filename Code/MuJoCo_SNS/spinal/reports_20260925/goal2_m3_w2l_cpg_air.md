# Goal 2 / Milestone 3 — Walker_2_Layer_CPG 2023 ORIGINAL (single LH RG) air-stepping on the M1 MuJoCo body

**Date:** 2026-09-25 (EB475WS4, unattended campaign) · **Status: complete — GATE PASS** (20 s clean air stepping, alternating L/R, frequency within 2× of the 2023-original reference)
**Artifacts (all NEW, under `Code\MuJoCo_SNS\spinal\w2l_mujoco\`; no protected file touched):**
- `fix_joint_axes.py` — **repairs a real M1 transport defect** (below); emits `w2l_mjcf_fixed.xml` from the aproj ground truth; verification PASS 8/8 lateral axes
- `w2l_mjcf_fixed.xml` — the axis-fixed body (M1's `w2l_mjcf.xml` and `make_w2l_mjcf.py` left UNTOUCHED)
- `w2l_air.xml` — air rig: same body with the welded Root lifted +0.30 m (written by the gate)
- `build_w2l_orig_net.py` — the 2023-original neural architecture transcribed from the `w2laproj` template (92 nodes / 163 edges)
- `test_w2l_air.py` — **the M3 gate** (≥20 s air stepping; exit 0 only on PASS)

Gate logs: `reports_20260925\logs\test_w2l_air_run1.log` (metrics bugs), `test_w2l_air_run2.log` (PASS), `test_w2l_air_run3.log` (**final PASS, exit 0 confirmed**). Probe trail: `reports_20260925\tmp\probe_w2l_actuators.py`, `probe2_directions.py`, `probe3_air_directions.py`, `probe4_zerog_directions.py`, `probe5_kinematic_arms.py` / `probe5b_kinematic_fixed.py`, `probe6_knee_debug.py`, `probe7_world_axes.py`, `rhythm_orig_sweep.py`, `cloop_smoke.py`.

Env: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe` (mujoco 2.3.7, sns-toolbox 1.5.2, `CONDA_PREFIX` set before the mujoco import). No pip installs.

---

## 1. THE HEADLINE FINDING — M1 shipped knee/ankle joints as YAW joints (now fixed)

The M1 body as delivered could not step for a reason no tuning could fix: **the knee and ankle hinge axes were vertical (world −z) instead of lateral**, i.e. both joints rotated about the gravity axis.

- Evidence on the SHIPPED xml (`probe7_world_axes.py`): knee_L world axis `(−0.052, 0.000, −0.999)`, ankle_L `(−0.000, 0.000, −1.000)` — vertical. Hip/toe were correct (`(0, ±1, 0)`).
- Consequence (`probe5_kinematic_arms.py`, kinematic finite differences of tendon length vs joint angle): knee muscles had **0.00–0.36 mm of tendon travel per +10°** of "knee" angle — zero sagittal moment arm. The knees/ankles could produce no flexion torque at all.
- This is what the earlier "passive crumpled kneel" actually was: the legs **twirled** about vertical axes (M1 gate-2's knee "+68°" was yaw; probe3 showed knee rotation moving the toe 11 mm SIDEWAYS and 0.000 vertically). It also reframes M2's failure diagnosis: its "knees +68..77 through the limit / ankle slams" body-side collapse was largely this artifact, on top of the (real but secondary) missing muscle damping.
- **Root cause** (`make_w2l_mjcf.py`, emit_body, ~line 402): the joint frame was composed with the body's **parent-relative** rotation `b.R_al` instead of its **world** rotation `b.R_world`, AND the resulting world-frame vector was written into MuJoCo's **body-local** `axis` attribute. Two errors that cancel for near-identity bodies (femur, toe) and don't for the 90°-rotated tibia/foot. The generator's own DUMP code (line ~521) uses `R_world` correctly — so the M1 report §4b world-axis table was always right and is the reference this fix reproduces.
- **The fix** (`fix_joint_axes.py`, new file): recompute `world_AL = R_world(child) @ R_joint(local) @ e_x` per joint from the aproj (read-only), convert to MuJoCo body-local coordinates with MuJoCo's own body rotation (`axis_local = R_worldbody_MJᵀ @ world_MJ`), write `w2l_mjcf_fixed.xml`. Verification in the same run: **all 8 hinges world-lateral (|y| > 0.99) — "AXIS FIX: PASS"** (run output §4).
- After the fix the same kinematic probe shows real arms: knee flexor 0.29 mm/°, knee extensor 0.96 mm/°, ankle PF 0.63 mm/°, ankle DF 1.38 mm/°, hip flexor 1.73 mm/°, hip extensor 1.38 mm/° — every joint's two actuators with opposite signs (table in §4).
- M1's gate could not see this: it validated masses, ranges, anchors, and the rest pose — never axis directions. Recorded as a standing lesson: **axis directions need a gate.**

## 2. The neural architecture — 2023 ORIGINAL, single LH RG

`build_w2l_orig_net.py` transcribes `connectome_templates.json` key **`w2laproj`** (92 nodes / 163 edges, mined from `Walker_2_Layer_CPG.aproj`) mechanically — it designs no topology:

- **Single LEFT RG**: `L RG ext`/`L RG flx` persistent-Na half-centers (spinal proven NAP set, **fixed tau_h** via `SNS_NumpyFixedTau` — AnimatLab semantics; stock tau_h(V) quenches these), laminated InE/InF inhibition (g 4.0) + direct HC↔HC escape excitation (g 0.5).
- **One RG powers BOTH sides' PF layers via 8 `pf_drive` edges**: 4 ipsilateral (L RG ext→L Hip/Knee PF ext; L RG flx→L PF flx) + **4 CROSSED antiphase** (L RG ext→R Hip/Knee PF **flx**; L RG flx→R PF **ext**). That crossing is the whole interleg coupling — there is no R RG, no commissural neurons.
- PF layers for hip+knee only (HC-PF-E/F + IN-PF lamination, g 4.0 cross); **ankle MNs ride the knee PFs** (`pf_to_mn`, the biarticular routing mined from the aproj).
- **Renshaw**: 11 RCs with the original's shared knee-extensor quirk transcribed as-is (`L Knee MN ext RE` is excited by BOTH sides' knee-ext MNs and inhibits both).
- **Ia chain** (SN-Ia → IN-IaIN → antagonist MN, mutual IaIN inhibition) and **Ib** (SN-Ib → homonymous MN exc) built and asserted but receive zero current — the air condition; nothing sensory reaches RG/PF (Deng doctrine).
- Kickoff: `Stimulus_1 → L RG ext`, 10 nA for 10 ms (template gain 10, aproj waveform) — the endogenous start.
- Conventions/gains verbatim from `w2l_cpg\build_w2l_net.py` (0..5 mV scale, `W2L_GAINS` calibration table; unused contact/c1/V3 keys simply don't occur in this template).
- Build census (asserted in code, run log §4): **79 neurons / 150 synapses / 3 inputs / 12 outputs**; per-tag counts sum exactly: rg_laminate 6, pf_drive 8, pf_cross 16, pf_to_mn 12, rc_own 24, rc_mutual 12, ia_recip 24, ia_mutual 12, ib_auto 12, other 24 (=163 − 12 mn→muscle pass-throughs − 1 stimulus-becomes-port).
- **Muscle map by anatomy, not name** (the aproj's R-hip muscle names are crossed vs their attachment sites, M1 §4c; kinematically confirmed this session, §1 table): `R Hip MN flx → hip_R_ext`, `R Hip MN ext → hip_R_flx`; all other pools map to the same-named actuator (knee flex→*_flx, ankle plantarflexion→*_ext, dorsiflexion→*_flx). With this map **flexion = −qpos on every joint, both sides** — a symmetric convention (flexion-positive degrees below).

## 3. Drive regime and body-side stand-ins (documented, all knobbed)

- **Neural regime** (`rhythm_orig_sweep.py`): tonic is REQUIRED — pure kickoff latches E (E max 4.49 mV, duty 1.00, 0 bursts; the NaP + direct-excitation pair self-holds). Working window: te=3, tf=4, tau_rg_nap_h=0.25 → **endogenous period 1.027 s (0.97 Hz), antiphase r −0.966, duty 0.44** (21 bursts / 21 s). te=2/tf=3 → 1.55 s; tau=0.35 scales ~1.3–2.0 s (the documented period knob). Defaults in the gate = te 3, tf 4, tau 0.25.
- **MN→ctrl**: stack convention `a = clip(V/5 mV, 0,1)`, capped at `ctrl_cap=0.5`.
- **Body stand-ins for the unported LinearHill B (400–800 N·s/m; M1 §3 "B — none")**: runtime joint damping `dof_damping=3.0 N·m·s/rad` + stiffer joint limits (`jnt_solimp [0.9,0.99,...]`, `jnt_solref (0.006,1)`) — the default-soft limits let full-strength 1000–1500 N muscles violate them by 8–65° (measured in probes 2/3). All runtime-only (model arrays, not the XML), all exposed as gate knobs.
- **Air rig**: pelvis welded at standing height (as M1 shipped — see §6 note) + Root lifted +0.30 m in `w2l_air.xml` so no leg pose can reach the plane (leg reach ~0.95 m); ground contact measured 0 all run.

## 4. GATE — final run (paste of `logs\test_w2l_air_run2.log`; command: `C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_w2l_air.py`, cwd `w2l_mujoco\`)

```
== M3 gate: W2L 2023 ORIGINAL (single LH RG) air stepping on the axis-fixed M1 body, 20 s ==
   knobs: te=3.0 tf=4.0 tau_rg_nap_h=0.25 ctrl_cap=0.5 joint_damp=3.0 stiff_limits=1 lift=0.3
   ground contacts: 0 (must be 0) | leg-leg self contacts: 22632
   neural: L RG ext bursts=18, period 1.027 s (0.97 Hz), E max 6.99 mV
   hip L swing excursions: 23, mean interval 0.583 s (1.72 Hz); ACF dominant period 1.027 s (0.97 Hz)
   hip R swing excursions: 18, mean interval 0.994 s (1.01 Hz); ACF dominant period 0.269 s (3.72 Hz)
   joint flexion-positive excursions (deg, min..max and range) vs AnimatLab references:
     hip   L [ -24.5, +17.3] range  41.8 | R [ -24.7, +17.1] range  41.8   (ref range ~38)
     knee  L [  -4.5, +64.0] range  68.5 | R [  -4.1, +64.7] range  68.8   (ref range ~61)
     ankle L [  -2.1, +26.6] range  28.8 | R [  -2.6, +25.0] range  27.5   (ref range ~16)
   hip flexion L/R correlation (0 lag): r = -0.623 (antiphase want < -0.3)
   [PASS] finite_20s
   [PASS] both_hips_swing
   [PASS] antiphase
   [PASS] airborne
   [PASS] freq_within_2x
   [PASS] hip_amp_ok
VERDICT: PASS  period~0.974s (hip L ACF)  antiphase_r=-0.623  hip_range_L=41.8 deg
```

**Reference comparison (honest):**
- **Frequency ~0.97 Hz**: within 2× of the **2023 original (0.77 Hz)** band [0.39, 1.54] ✓; vs the modern W2L air reference 2.22 Hz it is 2.3× slow — outside that band (the modern walker's faster cadence would need tau_rg_nap_h ≈ 0.12–0.15 s, below the 0.25–0.35 range the ask gave; the knob is there if wanted).
- **Hip range 41.8° vs ref 38°** — matched (+10%). **Knee 68.5/68.8° vs 61°** (+12%, and includes ~4–8° of residual limit violation past 60° flexion — bounded, stable). **Ankle 27.5–28.8° vs 16°** (+75% — the ankle is the least faithful joint: its transported range [−20°,−5°] excludes the rest pose, so it operates pinned near one limit and violates it under load; honest gap, see §5).
- **Antiphase r = −0.623** over the full 20 s (weaker than the RG's own −0.97 because leg-leg contact impulses perturb the body traces; excursion counts 23 L / 18 R and the R-hip ACF's fast 0.269 s component are the same contact signature).
- **Neural period in the loop (1.027 s) equals the open-loop sweep (1.027 s) exactly** — the RG is autonomous in the closed loop, i.e. the Deng doctrine (nothing sensory reaches RG/PF) is verified behaviorally, not just by wiring.
- **ground contacts = 0** (the gate's "no ground contact" requirement, measured every step). leg-leg self contacts 22632 (~113/s; the antiphase deep-knee swings bump at midline crossing) — reported, not ideal; they are part of this body's measured behavior.

## 5. Deviations and caveats (all loud)

1. **Axis fix + knee-range negation are body changes vs the M1 artifact.** `w2l_mjcf_fixed.xml` = M1 xml with 8 axes corrected (§1) and knee ranges [0,+60°] → [−60°,0] (flexion = −qpos by muscle-geometry evidence; the aproj's LowerLimit=0 is "straight" so the transported sign was inconsistent with the axes). Hip/ankle ranges left as transported. The one-line generator fix for whoever regenerates: in `emit_body`, use `b.R_world` and emit `axis_local = R_world_MJᵀ @ al2mj(world_AL)`.
2. **Runtime-only stand-ins** (§3): joint damping 3.0, stiffer limits, ctrl cap 0.5 — substitutes for AnimatLab's unported muscle damping B and stiff limit dynamics; not tuned, just measured sufficient.
3. **Ankle fidelity is poor** (range +75% over ref, limit violations, transported rest-pose/range inconsistency 0 ∉ [−20,−5] — an M1 finding re-flagged here; resolving it properly needs the aproj's joint-coordinate offset, not more tuning).
4. **Leg-leg self-contact** during deep antiphase flexion (§4) — no ground contact, but the legs are not perfectly independent in air.
5. **The R-hip crossing** means the template's "R Hip flx/exc" MN pools drive actuators named oppositely (`hip_R_ext`/`hip_R_flx`); the template MUSCLE labels were mined from the aproj so the crossing is in the MUSCLE_MAP table, commented, in `build_w2l_orig_net.py`.
6. **Chaos caveat** (standing, AGENTS 2026-09-13): eval-vs-eval bit-exactness holds only for identical binaries/seeds; the gate is deterministic per run but the numbers above are not a cross-platform contract.

## 6. One more M1-record correction found on the way

The M1 report says the Root ships FREE ("the Freeze trap deliberately NOT copied") — it does not: the shipped `w2l_mjcf.xml` Root has **no joint** (`body_jntadr −1`, `nq=8` in M1's own gate log), because `make_w2l_mjcf.py` only emits joints found in source bodies and the aproj Root has none. The M1 "passive drop" was leg collapse under a fixed pelvis, and M2's `qpos[2]` "pelvis height" was actually `ankle_L` radians (its −1.44 "m" = the −82° limit slam M2 itself reports). For M3 this is good news — the welded pelvis IS the air rig (and matches how the AnimatLab standalone-asim air reference was measured). Left unmodified; the gate lifts the welded Root +0.30 m in a separate air XML.

## 7. What is NOT done / natural next steps

- No ground phase: contact encoders, Ia/Ib/II encoder wiring (the SN ports), and the standing-pose solve are all future work; this build has no sensory input at all by design.
- The `w2l_cpg` bilateralrg variant (real contact-driven stepping) is the sibling port — this milestone only claims the 2023 ORIGINAL architecture.
- If the modern 2.22 Hz cadence is wanted: `--tau=0.15` is available but outside the sanctioned 0.25–0.35 range and untested here.
- Knob defaults chosen by hand (§3 sweep, 15 configs, ~4 min); no optimizer was run — out of scope for the time budget.
