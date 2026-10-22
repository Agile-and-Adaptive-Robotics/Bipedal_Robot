# Co-contraction stand descent — DIRECT experiment (2026-10-21)

Ben's directive (2026-10-21): "Massively activated antagonistic muscle
groups should free-stand — the trick is doing it energy efficiently. So
heavily activate antagonists so the torques balance (legs become rigid
struts), THEN reduce from there." The prior s3q study penalized
co-contraction (metabolic term) BEFORE proving the capability —
backwards. This experiment ran the capability test first.

**Method.** Mechanism = the per-side ALPHA hand-drive cells (uniform
MN-pool drive onto EVERY ipsilateral pool, extensors and flexors alike;
the 2026-10-16 Shevtsova brainstem-analog wiring, all three families
inherit it), driven by the new default-off `runner.py --alpha-stand`
hook (additive; anchor re-verified BIT-EXACT after the edit:
`_anchor_coactstand_20261021.log`). Config per family = its QUIET
stage-3 incumbent (v10 multipliers + renshaw + quiet-winner params,
mechanism knobs 0), `--stand-eval 8 --rig-scale 0` (rig fully OFF;
ligament surrogates stay). Driver `_coact_stand_20261021.py`; results
JSONs `coact_stand_{stock,syn6,w2lvar}_20261021.json` beside this file;
logs `_coactstand_*_20261021.log` in `spinal\`.

## 1. EXISTENCE PROOF — the answer is NO, and the failure is measured

Every family, at the maximum drive the mechanism can deliver (port
current 5 nA saturates the ALPHA cell; conductance raised until all 92
MN pools clip at activation 1.0 — verified `act_sum = 92.0`,
`bal_effort = 0.9993`, co-contraction index 0.9999), FALLS within
~1.5 s of rig release, at rig 0:

| family | MAX lock config | all pools at 1.0? | verdict | first joint divergence (t) | COM −8 cm at |
|---|---|---|---|---|---|
| stock (s3k) | g=2.0, I=5 nA | YES (act_sum 92.0) | FALLS | pelvis_rotation +3.9 deg (t=1.006 s; pelvis_tilt/ankles/hips all by 1.02 s) | t=1.06 s |
| syn6 | g=2.0, I=5 nA | near (min antagonist mean 0.86, act_sum 88.9) | FALLS | pelvis_rotation (t=1.006 s) | t=1.06 s |
| w2lvar | g=4.0, I=5 nA | YES (act_sum 92.0) | FALLS | pelvis_rotation (t=1.006 s) | t=1.06 s |

**THE MECHANISM GAP, measured statically** (`_coact_torque_probe_
20261021.py`: every muscle activation set to 1.0 at the standing pose,
`actuator_moment.T @ actuator_force`):

- **The ankle cannot be locked by any uniform drive.** Net muscle
  torque at full co-activation: ankle_r **−182.7 N·m**, ankle_l
  **−187.4 N·m plantarflexor** (gravity needs ~1 N·m). Components R:
  PF soleus −118.5 + med_gas −61.0 + lat_gas −28.7 (sum −236.0) vs DF
  tib_ant +34.5 + ext_dig +18.8 (sum +53.3) — a **4.4:1 PF:DF torque
  mismatch**. Balancing the ankle would need roughly DF 1.0 vs PF 0.23
  — no uniform, and no co-contraction-of-equal-magnitudes, can do it.
- The knee IS torque-balanced at full co-activation (−8.3 / +11.9 N·m;
  it pins into the +10 deg extension stop and strut-holds). The hips
  are extensor-dominant (R −152.9, L −17.8 N·m — the split-stance
  pose), lumbar balanced (−0.1).
- **Dynamics at MAX (stock, `_coactstand_stock_one_20261021.npz`):**
  rig release at 0.6 s; both ankles snap into plantarflexion within
  0.1 s (−35/−27 deg at t=0.7), reach −44/−66 deg by t=0.9 s; the body
  rises onto its toes, the COM passes the toe pivot, pelvis pitches
  (+42 deg by 0.9 s) and the model cartwheels forward to an INVERTED
  heap by ~2 s (pelvis tilt −220 deg, knees at the +10 extension stop,
  hips 40–70 deg flexed, ankles −50..−70). It ends hanging upside-down
  with its feet pressing the ground (the "feet ~100 %BW" numbers on
  fallen runs are heap artifacts, NOT standing loads).
- Probes: symmetric start pose falls identically (ankles −65/−69 deg
  by t=0.9 s); the coact knob at its max (1.0, no ALPHA) also falls —
  by construction it caps flexor-antagonist drive at the extensor
  POSTURE weight, a much weaker regime (act_sum 54.3, ankles still
  −60 deg by 0.9 s).

So the failure is NOT an activation-to-force ceiling (all pools reach
1.0 and the muscles deliver), and NOT joint range per se — it is the
**agonist/antagonist force mismatch at the ankle** (plus unactuated
pelvis orientation): uniform co-contraction of a 4:1-unequal antagonist
pair produces a large NET plantarflexor torque, which destroys the
support geometry before any descent could start.

## 2. DESCENT CURVE — there is no standing level to descend from

Fresh 8-s runs per level (drive fractions of the 5 nA lock; rig 0,
stock; syn6/w2lvar in the JSONs — identical verdicts). Steady-window
values on fallen runs describe the inverted heap, not a stand:

| level (of lock) | ALPHA current nA | mean antagonist activation (hip/knee/ankle/trunk) | stands? | act_sum (energy proxy) |
|---|---|---|---|---|
| 100 % | 5.0 | 1.00 / 1.00 / 1.00 / 1.00 | FALLS | 92.0 |
| 80 % | 4.0 | 0.98 / 1.00 / 1.00 / 1.00 | FALLS | 91.5 |
| 60 % | 3.0 | 0.90 / 0.94 / 0.99 / 1.00 | FALLS | 89.3 |
| 40 % | 2.0 | 0.79 / 0.84 / 0.94 / 1.00 | FALLS | 83.2 (fresh-process re-run identical: 104.8 %BW / −0.66 / 83.2) |
| 20 % | 1.0 | 0.60 / 0.70 / 0.83 / 0.51 | FALLS | 71.6 |
| 10 % | 0.5 | 0.47 / 0.59 / 0.74 / 0.34 | FALLS | 63.4 |
| 0 % (incumbent) | 0.0 | 0.30 / 0.45 / 0.63 / 0.15 | FALLS | 52.7 |

The 56-s staircase run (lock 8 s then 6 steps down) confirms it: the
model topples during the FIRST lock window (COM global min −0.81 m);
every later window reads the settled heap (COM ≈ −0.67 m, knees at the
+10 stop). **Minimum level that still stands: NONE — the descent curve
for free-standing does not exist for this mechanism.** The energy
question ("how little co-contraction suffices?") is answered upstream:
zero levels stand, so the efficient operating point is undefined until
the ankle balance problem is solved.

## 3. Other models

- **syn6 and w2lvar**: mechanism exists (same inherited ALPHA cells);
  run — same verdict at every level (JSONs above). Their unrigged
  incumbents fall identically.
- **W2L Hill robot and BPA walker: NOT run (noted).** Neither stack has
  the ALPHA/coact mechanism (w2l_mujoco bilateralrg2 / tonic BPA drive,
  no such cells), so there is nothing to descend; their existing
  unrigged results already answer capability (W2L unrigged step regime
  falls at t=1.18 s; BPA body cannot statically stand — COM outside the
  support polygon, per VERDICT_TABLE_20261021.md).

## What this buys the dissertation

The honest capability statement is now TWO-SIDED and measured: (a) the
SNS standing solutions cannot free-stand at ANY uniform co-contraction
level — the blocking mechanism is a static 4.4:1 PF:DF ankle torque
imbalance (plus the unactuated pelvis base), not insufficient drive;
(b) the balance objective's metabolic penalty on co-contraction was
therefore not the binding constraint — the s3q search-order critique is
resolved by capability-first measurement. The next mechanism to add, if
free-standing is the goal, is joint-specific torque balancing at the
ankle (foot geometry / DF-side reinforcement / graded PF:DF ratios as
the NNLS standing solve already computes) — not more uniform drive.

Fresh-process verification: MAX lock re-run standalone reproduces the
matrix rung exactly (act_sum 92.0, sway 0.903, COM −0.641, feet 104.3);
40 % rung re-run reproduces exactly (83.2 / 104.8 / −0.66). Anchor
−161.56754173676563 BIT-EXACT after the runner edit.
