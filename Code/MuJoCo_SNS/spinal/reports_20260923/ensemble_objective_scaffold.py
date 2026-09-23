"""ENSEMBLE-OBJECTIVE SCAFFOLD -- DESIGN SKETCH ONLY (Goal 4, 2026-09-23).

*** NOT WIRED IN. NOT IMPORTED BY ANYTHING. NOT TESTED AS A STUDY. ***
This file is a design sketch for the "perturbed variants trained jointly"
objective (Ben's question). It documents the intended structure so the next
study script can be written from it deliberately. Do NOT run this expecting
an optimization: `objective_ensemble` is complete in form but the plant-
jitter variant needs a runner-side Fmax/damping hook that does not exist
yet (see VARIANT NOTES below).

Design (see reports_20260923/goal4_optuna_review.md section 4):

  score(trial) = mean_i(s_i) - LAMBDA * std_i(s_i),   i = 0..N-1 variants

  variant 0  = nominal (today's --eval run, exactly reproducible)
  variant 1  = plant jitter      (Fmax per muscle group, leg damping, rig k)
  variant 2  = reference/pose    (2nd left-leg reference cycle; pelvis_ty
                                  offset via AARL_PELVIS_TY, which the
                                  runner already reads)

  feasibility = >= K of N variants return a real kine dict (not NaN, not
  the frozen sentinel) -- handled natively by optuna if the sentinel is
  moved into TPESampler(constraints_func=...) instead of the score.

Cost: measured median eval time ~1.0-1.1 min (v10/v8b CSV timestamps), so
N=3 costs ~3 min/trial -> 60 trials ~3 h. The 2-STAGE mode below keeps the
screening stage at today's cost and only ensemble-rescores the finalists.

Wiring points that already exist (read from, don't edit, per study rules):
  - _curriculum.set_stage() already drives AARL_KY / AARL_PELVIS_TY envs
    (_curriculum.py lines 99-110).
  - kine_ref.compare(t, q, neuro, walk_start, ref=..., contact=...) accepts
    a per-call `ref` dict (kine_ref.py line 193), so reference variants are
    a precomputed list built once via kine_ref.load_reference()-style code
    with different onset windows.
  - runner --eval returns the metrics dict (runner.py lines 1507-1549).

VARIANT NOTES (what is missing):
  - Plant jitter needs a runner-side hook to scale actuator gainprm[:,2]
    (Fmax) per variant. TODAY no env/flag exists for that; adding one is a
    runner.py change that must keep default-off = bit-identical (same
    contract as every conditional-topology knob).
  - Rollout-level seeding: the runner has no --seed flag; run-to-run
    divergence currently comes from BLAS summation-order nondeterminism
    (documented in AGENTS.md optimizer-relevant insights / basin_gate
    caveat). If deterministic per-variant seeds are wanted, a seed env must
    be added the same way.
"""
from __future__ import annotations

import numpy as np

# ---- knobs (proposed starting values; see review section 4 rationale) ----
N_VARIANTS = 3          # nominal + plant-jitter + reference/pose
LAMBDA = 0.25           # score units per std unit (kine scores ~-60..-320)
K_FEASIBLE = 2          # >=2 of 3 variants must produce a real kine dict
TOPK_REScore = 8        # 2-stage: ensemble-rescore only this many finalists
TOPK_DELTA = 5.0        # ...within this many score points of the best

# Plant-jitter proposal (variant 1). Fmax groups follow the runner's muscle
# grouping convention (params.W_PF_MN keys are joint groups, but the model
# Fmax lives in actuator gainprm[:,2] -- hence the missing hook above).
PLANT_JITTER = dict(
    fmax_group_rel_sigma=0.07,   # ~+/-7% lognormal per muscle GROUP
    groups=("hip_ext", "knee_ext", "ankle_pf", "ankle_df", "biart"),
    leg_damping_rel_sigma=0.10,  # runner --leg-damping exists
    rig_k_rel_sigma=0.10,        # AARL_KY exists (lateral rig)
    pelvis_ty_mm_sigma=5.0,      # AARL_PELVIS_TY exists
)


def variant_specs(trial_number: int) -> list[dict]:
    """Deterministic per-trial variant table (seeded by trial number).

    Same trial -> same variants (reproducibility); different trials ->
    different draws (domain randomization across the study). SCONE's analog
    is NoiseController{random_seed=...} + initial_state_offset sampling
    (Tutorials 3b/4c, local copies).
    """
    rng = np.random.default_rng(10_000 + trial_number)
    specs = [dict(name="nominal")]  # variant 0: exactly today's eval
    specs.append(dict(
        name="plant",
        fmax_scale={g: float(np.exp(PLANT_JITTER["fmax_group_rel_sigma"]
                                    * rng.standard_normal()))
                    for g in PLANT_JITTER["groups"]},
        leg_damping_scale=float(np.exp(
            PLANT_JITTER["leg_damping_rel_sigma"] * rng.standard_normal())),
        rig_k_scale=float(np.exp(
            PLANT_JITTER["rig_k_rel_sigma"] * rng.standard_normal())),
    ))
    specs.append(dict(
        name="refpose",
        # left leg has a 2nd cycle in the reference GRF (measured:
        # 3 left onsets); pose jitter through the existing env
        ref_cycle_left=1,
        pelvis_ty_offset_mm=float(
            PLANT_JITTER["pelvis_ty_mm_sigma"] * rng.standard_normal()),
    ))
    return specs


def score_ensemble(scores: list[float], feasible: list[bool]) -> float | None:
    """mean - LAMBDA*std over feasible variants; None if infeasible."""
    if sum(feasible) < K_FEASIBLE:
        return None  # caller reports via constraints_func, NOT a sentinel
    s = np.array([x for x, f in zip(scores, feasible) if f], dtype=float)
    return float(s.mean() - LAMBDA * s.std(ddof=1 if len(s) > 1 else 0))


def objective_ensemble(trial):  # pragma: no cover - sketch, needs runner hook
    """Form of the optuna objective (to be adapted into the NEXT study
    script -- this scaffold must not be imported by existing studies).

    Pseudo-outline:
      p = suggest(...)                      # narrowed space (see review 3)
      specs = variant_specs(trial.number)
      scores, feas = [], []
      for spec in specs:
          apply_variant(spec)               # envs today; Fmax hook TODO
          m = R.main(["--eval", "--drive", repr(p["drive"])])
          feas.append(not m["nan"] and m.get("kine") is not None)
          scores.append(float(m["kine_score"]) if feas[-1] else float("nan"))
      trial.set_user_attr("variant_scores", scores)   # survives in the db
      val = score_ensemble(scores, feas)
      if val is None:
          raise optuna.TrialPruned()  # or constraints_func path
      return val
    """
    raise NotImplementedError("design scaffold -- see module docstring")


def rescore_finalists(study_results: dict[int, float]) -> set[int]:
    """2-stage budget saver: trials within TOPK_DELTA of the best, capped
    at TOPK_REScore, get the N-variant treatment; everything else is scored
    nominal-only. Mirrors the basin_gate philosophy (robustness gate on
    finalists) without paying N x on every trial."""
    best = max(study_results.values())
    keep = sorted((t for t, v in study_results.items()
                   if best - v <= TOPK_DELTA),
                  key=lambda t: -study_results[t])[:TOPK_REScore]
    return set(keep)


if __name__ == "__main__":  # smoke: the pure-math parts work
    s = score_ensemble([-160.0, -175.0, -168.0], [True, True, True])
    print("ensemble smoke:", round(s, 3))
    print("finalists smoke:", sorted(rescore_finalists(
        {0: -160.0, 1: -161.0, 2: -200.0, 3: -159.0, 4: -300.0})))
