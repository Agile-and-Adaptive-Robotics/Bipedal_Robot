"""Parallel evaluation pool for the W2L per-synapse campaign (2026-10-26).

Uses multiprocessing to evaluate trials in parallel (one worker per core).
Each worker pre-imports the evaluation stack (test_w2l_ground,
_w2l_wean_eval) and evaluates trials independently. The sparse neural
backend (fast_eval.SparseForward) is available but the full run_phase
is used for correctness — the parallel pool alone gives the first
order-of-magnitude speedup.

Usage:
    from parallel_eval import ParallelCampaign
    pc = ParallelCampaign()
    results = pc.evaluate_batch([knobs1, knobs2, ...])
    # or integrate with optuna:
    pc.evaluate_optuna(study, n_trials=500)

Speedup: ~10x (number of cores) over sequential evaluation.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
W2L_MUJOCO = HERE.parent / "w2l_mujoco"
WEAN = HERE.parent / "reports_20261016" / "w2l_wean"

for p in (str(W2L_MUJOCO), str(WEAN), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)


# ---- worker functions (run in child processes) ----

def _init_worker():
    """Called once per worker process: import the evaluation stack."""
    global _EV, _T
    import test_w2l_ground as T
    import _w2l_wean_eval as EV
    _EV = EV
    _T = T
    # pre-build the network so the first trial doesn't pay the build cost
    T.KV["net"] = "bilateralrg2"


def _eval_one(knobs: dict) -> dict:
    """Evaluate one trial. Returns the metric dict."""
    global _EV, _T
    # make all keys KV-resident (forwarding fix)
    for k, v in knobs.items():
        _T.KV.setdefault(k, float(v))
    return _EV.evaluate(knobs, tag=f"par_{os.getpid()}",
                       csv_path=None)  # no csv from workers


# ---- parent-process controller ----

class ParallelCampaign:
    """Evaluate W2L trials in parallel across all available cores."""

    def __init__(self, n_workers: int | None = None):
        if n_workers is None:
            n_workers = min(mp.cpu_count(), 10)  # cap at 10 (10-core box)
        self.n_workers = n_workers
        self.pool = None

    def start(self):
        """Start the worker pool."""
        self.pool = mp.Pool(self.n_workers, initializer=_init_worker)

    def evaluate_batch(self, knobs_list: list[dict]) -> list[dict]:
        """Evaluate a batch of knob dicts in parallel."""
        if self.pool is None:
            self.start()
        return self.pool.map(_eval_one, knobs_list)

    def evaluate_optuna(self, study_name: str, storage: str,
                        n_trials: int = 500, base_knobs: dict | None = None):
        """Run an Optuna study with parallel evaluation.

        Uses a custom ask-and-tell loop: generates a batch of parameter
        suggestions, evaluates them all in parallel, then feeds results
        back to Optuna. This is the pattern for PSO-style parallel
        exploration even with TPE.
        """
        import optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)

        study = optuna.create_study(
            study_name=study_name, storage=storage,
            direction="maximize", load_if_exists=True,
            sampler=optuna.samplers.TPESampler(seed=20261026))

        # get the search space from the first existing trial or a seed
        if not study.trials and base_knobs:
            study.enqueue_trial(base_knobs)

        batch_size = self.n_workers
        n_done = len(study.trials)

        while n_done < n_trials:
            n_this = min(batch_size, n_trials - n_done)

            # ask: generate a batch of trials
            trials = []
            for _ in range(n_this):
                trial = study.ask()
                trials.append(trial)

            # build knob dicts from the asked trials
            knob_list = []
            for trial in trials:
                # the trial already has suggested params
                knobs = dict(trial.params)
                if base_knobs:
                    for k, v in base_knobs.items():
                        knobs.setdefault(k, v)
                knob_list.append(knobs)

            # evaluate all in parallel
            t0 = time.time()
            results = self.evaluate_batch(knob_list)
            elapsed = time.time() - t0

            # tell: feed results back
            for trial, result in zip(trials, results):
                study.tell(trial, float(result.get("obj", -100)))

            n_done = len(study.trials)
            best = study.best_value if study.trials else float("-inf")
            print(f"  batch: {n_this} trials in {elapsed:.1f} s "
                  f"({n_this/elapsed:.1f} trials/s) | total {n_done} | "
                  f"best {best:.1f}", flush=True)

        return study

    def close(self):
        if self.pool:
            self.pool.close()
            self.pool.join()
            self.pool = None


if __name__ == "__main__":
    # quick test: evaluate 4 trials in parallel vs sequential
    import json
    winner = json.load(open(HERE / "round_w25.json", encoding="utf-8"))
    P = winner["params"]

    base = dict(dscale=1.0, te=3.0, tf=4.0, wsup=0.25, sagittal=50.0,
                dur_walk=10.0, legcoll=0.0, sink=0.06, rig=0.25,
                walker_yield=1.5, **P)

    # sequential
    _init_worker()
    t0 = time.time()
    r0 = _eval_one(base)
    t_seq = time.time() - t0
    print(f"sequential: {t_seq:.1f} s | obj {r0['obj']:.1f} "
          f"survive {r0['survive']}")

    # parallel
    pc = ParallelCampaign(n_workers=4)
    batch = [base] * 4
    t0 = time.time()
    results = pc.evaluate_batch(batch)
    t_par = time.time() - t0
    print(f"parallel (4 workers, 4 trials): {t_par:.1f} s "
          f"({4/t_par:.1f} trials/s)")
    print(f"speedup: {4*t_seq/t_par:.1f}x")
    print(f"parallel obj: {results[0]['obj']:.1f} "
          f"(match: {abs(results[0]['obj'] - r0['obj']) < 0.1})")
    pc.close()
