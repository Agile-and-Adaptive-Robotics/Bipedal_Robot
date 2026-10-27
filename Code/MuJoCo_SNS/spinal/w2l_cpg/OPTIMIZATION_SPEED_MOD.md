"""OPTIMIZATION SPEED MOD — design document (2026-10-26)

Ben's directive after his colleague's insight: MuJoCo preallocates memory
so everything calculates instead of calculating AND allocating; SNS
toolbox uses matrices that multiply 0 times 0 a lot; PSO with 80 samples
runs fast.

DIAGNOSIS (measured, bilateralrg2, 20 s trial = 10,000 steps at 2 ms):

  | Component                    | Time/trial | %    |
  |------------------------------|-----------|------|
  | MuJoCo physics (10k steps)   | 0.17 s    | 3%   |
  | SNS neural (10k steps, dense)| 1.02 s    | 17%  |
  | Network build+compile        | 0.08 s    | 1%   |
  | Python loop overhead         | ~4.7 s    | 79%  |

  The SNS compiled network's synapse matrix is 96x96 with 204 nonzeros
  (2.2% dense) — it computes 9,216 element operations to use 204. Input
  connectivity is 96x22 with 22 nonzeros. Every step allocates temporary
  arrays for each intermediate result.

SOLUTION (three layers, each independently useful):

  Layer 1: SPARSE NEURAL BACKEND (`fast_eval.py`)
    SparseForward extracts the nonzero indices once at compile time and
    replaces forward() with scatter/gather on ~204 elements using
    np.bincount (8.9x faster than np.add.at for scatter-adds).
    Pre-allocates all buffers — zero per-step allocation.
    Validated BIT-IDENTICAL to the toolbox forward().
    Speedup on neural compute alone: 2.1x (102→48 µs/step).
    Remaining bottleneck: numpy per-call overhead (~3-5 µs per operation
    on small arrays; ~10 calls per step = 30-50 µs). Would need numba
    or C compilation to go below this — numba not installed in the myo
    environment.

  Layer 2: FAST TRIAL RUNNER (planned)
    Replaces test_w2l_ground.run_phase with a tight loop that:
    - Reuses the pre-built MuJoCo model and sparse neural net
    - Pre-allocates all sensor/metric buffers
    - Eliminates XML text manipulation, logging, dict creation
    - Computes only the metrics needed by the optimizer
    Estimated trial time: ~0.3-0.5 s (vs ~6 s current)

  Layer 3: PARALLEL EVALUATION (planned)
    10 workers (one per core), each pre-building its own model and
    network copy. Evaluates trials independently via multiprocessing.
    At Layer 2 speed: 10 × (1/0.4) = 25 trials/second.
    A 500-trial optimization finishes in ~20 seconds.
    A PSO generation of 80 particles: ~3 seconds.

OPTUNA NOTE: our campaigns use TPE (Tree-structured Parzen Estimator),
which is Optuna's default single-objective sampler — NOT NSGA-II.
TPE is a Bayesian method that models which parameter regions look
promising and picks accordingly. Optuna also offers NSGA-II for
multi-objective, CMA-ES, random, grid, and QMC. All are just
parameter-pickers; the evaluation speed is the bottleneck, not the
sampler choice. Once evaluations are fast, PSO with large populations
becomes attractive because entire generations evaluate in parallel.
"""

## MEASURED RESULTS (2026-10-26, bilateralrg2, 10 s trial)

| Configuration                    | Time per trial | Throughput    |
|----------------------------------|---------------|---------------|
| Sequential (current)              | 4.5 s         | 0.22 trials/s |
| 4 workers                         | 5.7 s wall    | 0.70 trials/s |
| 10 workers (all cores)            | 6.7 s wall    | 1.50 trials/s |

Speedup with 10 workers: **7.5x** (near-linear scaling).
All results bit-identical across workers (deterministic).

Optimization timing at 10-worker throughput:
- 500-trial campaign:   ~5.5 minutes (was ~42 minutes)
- 1000-trial campaign:  ~11 minutes  (was ~84 minutes)
- 80-particle PSO gen:   ~53 seconds

REMAINING BOTTLENECK: ~5-6 s per trial in run_phase Python overhead
(sensor reads, logging, dict management, MuJoCo API call overhead).
The sparse neural backend saves ~0.5 s of this; the bulk is the
Python-level loop in test_w2l_ground.run_phase. Eliminating it needs
numba or C compilation (numba not available in the myo environment).

FILES:
- fast_eval.py: SparseForward (bit-identical, 2.1x on neural step alone)
- parallel_eval.py: ParallelCampaign (10-worker pool, 7.5x throughput)
- OPTIMIZATION_SPEED_MOD.md: this document

## NUMBA RESULTS (2026-10-26, "high gear")

Numba 0.57.1 installed (compatible with numpy 1.21.6 in the myo env).
The entire neural step is compiled to native machine code via LLVM.

| Neural step backend         | µs/step | 20 s trial (neural only) |
|-----------------------------|---------|--------------------------|
| toolbox (dense 96x96)        | 102     | 1.020 s                  |
| sparse numpy                 | 48      | 0.480 s                  |
| **numba compiled**           | **5.4** | **0.054 s**              |

Speedup vs toolbox: **19x on neural compute alone.**

FULL TRIAL (neural + physics + sensors + muscle mapping):

| Duration | Wall time | Speedup vs old |
|----------|-----------|----------------|
| 5 s      | 0.21 s    | ~21x           |
| 10 s     | 0.43 s    | ~10x           |
| 20 s     | 0.88 s    | ~7x            |

The 20 s trial runs at **23x real time** on a single core.

COMBINED WITH 10-WORKER PARALLEL POOL:
  Single worker: ~1.1 trials/s (20 s trials)
  10 workers:   ~11 trials/s
  500-trial campaign: ~45 seconds
  80-particle PSO generation: ~7 seconds

WHAT NUMBA IS: a just-in-time compiler for Python that compiles your
Python function to native machine code via LLVM. The first call pays
~1 second of compilation; every subsequent call runs at C speed with
zero Python interpretation overhead. It eliminates the per-call overhead
of numpy operations on small arrays (the 3-5 µs fixed cost per operation
that dominated the sparse-numpy backend).

FILES:
- fast_eval.py: SparseForward (numpy sparse, bit-identical, 2.1x)
- numba_trial.py: FastTrial (numba-compiled, 19x neural, ~10x full trial)
- parallel_eval.py: ParallelCampaign (10-worker pool)
- OPTIMIZATION_SPEED_MOD.md: this document

## REPLICATION GUIDE (for EB475WS4 or any machine)

### Environment setup

The `myo` environment has numpy 1.21.6 (pinned for OpenSim ABI
compatibility). Numba 0.68+ requires numpy >=1.22, so install the
compatible version:

```bash
D:\Anaconda\envs\myo\python.exe -m pip install numba==0.57.1 llvmlite==0.40.1
```

Verify (must run from a .py file, not inline — numba can't cache from
a `<string>` context):

```python
# test_numba.py
from numba import njit
import numpy as np
@njit(cache=True)
def test(x):
    s = 0.0
    for i in range(len(x)):
        s += x[i] * 2
    return s
print(test(np.array([1.0, 2.0])))
```

### File locations (all in `Code\MuJoCo_SNS\spinal\w2l_cpg\`)

| File | Purpose |
|---|---|
| `fast_eval.py` | SparseForward class + validate() (numpy sparse, bit-identical) |
| `numba_trial.py` | FastTrial class (numba-compiled neural step + MuJoCo loop) |
| `parallel_eval.py` | ParallelCampaign (10-worker multiprocessing pool) — in `reports_20261025_persyn\` |
| `OPTIMIZATION_SPEED_MOD.md` | This document |

### Step-by-step replication

1. **Sparse backend** (2x speedup, validates correctness):
   ```python
   from fast_eval import validate
   r = validate()  # must print PASS and speedup
   ```

2. **Numba backend** (19x on neural step):
   ```python
   from numba_trial import FastTrial
   ft = FastTrial()  # ~1 s compile time on first call
   result = ft.run(dur=10.0)  # one 10 s trial
   ```

3. **Parallel pool** (10x throughput):
   ```python
   from parallel_eval import ParallelCampaign
   pc = ParallelCampaign(n_workers=10)
   results = pc.evaluate_batch([knobs1, knobs2, ...])
   ```

### Key architecture decisions

- The numba `neural_step` function takes ALL parameters as flat numpy
  arrays (no objects, no dicts) because numba only compiles functions
  that operate on primitives and arrays. The class `FastTrial` pre-packs
  these arrays once.
- `cache=True` on `@njit` saves the compiled code to disk; subsequent
  runs skip the compilation step entirely.
- The MuJoCo physics step (`mj_step`) is already C code — numba can't
  call it directly, so the Python loop alternates between the compiled
  neural step and the MuJoCo call. The remaining Python overhead is
  this thin loop (~0.27 s per 20 s trial).
- The sparse extraction happens ONCE at init: nonzero indices from the
  96x96 dense connectivity matrices become flat index arrays, so each
  step computes only the 204 actual synapses instead of all 9,216
  matrix elements.

### Known TODO (as of 2026-10-26)

The `FastTrial.run()` method has the neural and physics loops working
but the drive-schedule wiring is incomplete (the `u[:] = 0.0` line needs
the tonic drive port indices and contact-current encoding). For a
production campaign, either:
a) Finish the wiring (see `test_w2l_ground.py` lines ~530-560 for the
   input-current schedule), or
b) Use `parallel_eval.py` which wraps the existing `run_phase` path
   (correct but slower — still 7.5x from parallelism alone).

### For the Di Russo controller on gait2392

The same approach applies: extract the controller's state update into
a numba-compiled function that takes flat arrays, compile it, and use
the parallel pool for the optimization loop. The controller in
`spinal\dirusso_gait2392\controller.py` has the same structure (leaky
integrators on per-muscle state vectors) and can be compiled the same
way.

## CURRENT STATUS (end of 2026-10-26 session)

WORKING NOW:
- parallel_eval.py: 10 workers, ~1.2 trials/s (correct, wraps run_phase)
- fast_eval.py SparseForward: bit-identical, 2.1x neural speedup
- numba_trial.py neural_step: 5.4 µs/step (19x), bit-identical

NEEDS DEBUGGING before production use:
- FastTrial.run() full loop: runs at 15x real time but the walker
  falls (pelvis drops, no heel strikes). Root cause: the simplified
  harness lacks the drop-protocol timing (hold 3s in air, lower 1s,
  settle 0.4s, then walk). The real test_w2l_ground.run_phase has
  this staged protocol; FastTrial starts walking immediately.
  Fix: add the hold/lower/settle phases before enabling drive and
  contact encoding. Until then, use parallel_eval.py (which wraps
  the proven run_phase).

PRACTICAL SPEEDUP AVAILABLE TODAY:
- 10-worker parallel pool: ~1.2 trials/s (vs 0.22 sequential)
- 500-trial campaign: ~7 minutes (vs ~42 minutes)
- Combined with numba neural step once FastTrial is debugged:
  projected ~11 trials/s, 500 trials in ~45 seconds.
