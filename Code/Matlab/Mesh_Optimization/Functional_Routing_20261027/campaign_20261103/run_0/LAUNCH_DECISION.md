# LAUNCH_DECISION.md — run_0_weights_placeholder (2026-11-04, easteregg2)

Ben's rule: the escape criterion is stated here BEFORE anything is launched.

## Escape criterion (as handed to the launcher)

Threshold: extrapolate the full campaign including `parpool(10)` overhead and
MATLAB startup; if the extrapolated wall clock is **≤ 36 h**, launch ALL 27
groups; otherwise launch only the demonstrator slice (psoas, iliacus,
glut_max3, add_mag3, vas_int, semimem + one foot route + erector spinae) and
queue the rest with reasons.

## Measured input (gate of record)

`campaign_20261103/output/gate_r4.log` (04-Nov-2026 05:20:48, machine
EASTEREGG2, GATE PASSED all checks; the ask states the three fresh-MATLAB
gate rounds exited 0 each, and gate_r1..r3 logs are present):

```
evaluations measured: 39 | median 0.128 s | mean 0.177 s
extrapolated full campaign (27 groups x {rigid, compliant} x
  [600 solve + 1 x 150 refine evals] = 40500 evals): 1.99 h
```

(The launcher ask's handed gate numbers — median 0.129 s / mean 0.172 s,
~2 h — match this round within measurement scatter.)

## Actual arithmetic

- Evaluation count: 27 groups × 2 branches (rigid, compliant Xi0–Xi3)
  × (600 solve + 150 refine evals) = **40,500 evaluations**.
- Gate mean applied to all of them: 40,500 × 0.177 s = 7,168.5 s =
  **1.99 h serial** (1.93 h at the ask's 0.172 s mean).
- MATLAB startup + `parpool(10)` startup overhead: bounded at **≤ 0.2 h**
  (observed fresh `-batch` startups ~0.5–1 min in this session; pool spin-up
  on 10 local cores typically 1–3 min).
- Unknown not in the gate mean: per-candidate retained-envelope obstacle
  scans grow for outer layers (1→26 envelopes). Conservative **×10 mean**
  bound (1.77 s/eval): 40,500 × 1.77 s = 71,685 s = **19.9 h SERIAL**, and
  `parfor(10)` per layer (layers ≤ 10 groups on this 10-core box) can only
  reduce wall below the serial bound — pessimistic 2× parallel-efficiency
  loss still ≈ 10 h wall, optimistic ≈ 2–4 h wall.

## Decision

Worst-case bound 19.9 h + 0.2 h overhead ≈ **20.1 h ≤ 36 h** →
**FULL branch: launch ALL 27 groups** (branch "full-27").

Launch-safety condition (not a re-scoping): every case is built by the new
contract-driven builder `fr_cases27.m` and must pass a pre-flight
build + seed evaluation on both branches (`fr_preflight27.m`) BEFORE the
campaign starts. A group that fails pre-flight is QUEUED with its exact
reason and the remaining groups launch; that fact will be reported in the
launch result and in `run_0/manifest.json`.

Weights are PLACEHOLDER-PENDING-BEN (`objective_spec.json`, hash recorded in
the run manifest at launch). Ben's weight ruling triggers a refinement
re-solve FROM CHECKPOINTS, not a rebuild.

## Pre-flight addendum (2026-11-04 06:07, before launch)

`fr_preflight27` PASSED: all 27 contract-driven cases built in 8.1 s and
seed-evaluated finitely on both branches (log: `campaign_20261103/output/
preflight27.log`). 54 seed evaluations: median 0.299 s, mean 0.354 s, max
0.980 s (compliant evals carry the Xi3 fixed point and are the slower
half). Direct serial extrapolation of all 40,500 campaign evaluations at
that mean: 40,500 × 0.354 s = **3.98 h**. The solve-budget evaluations run
the same evaluator the gate measured; retained-envelope obstacle scans are
additive tens-of-ms per candidate (measured code path: ≤26 envelopes ×
≤9 poses × ≤3 capsule/segment distance tests), so a realistic worst bound
is ≈ 2× serial ≈ 8 h, and `parfor(10)` per layer reduces wall below the
serial bound. Decision unchanged and re-confirmed: **full-27 branch**.
Seed-state honesty note for the record: several groups carry nonzero seed
constraint violations (e.g. sar 2.45, flex_hal 0.91) and the human target
curves are the stage-1 ΣFmax·(−dL/dθ) curves (peak-vs-contract-mifN gaps
−94% to −99%, the documented informational basis) — run_0's checkpoints
record these as-is; refinement from checkpoints follows Ben's weights.
