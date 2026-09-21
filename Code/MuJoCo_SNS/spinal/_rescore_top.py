"""Re-score the top s3c trials under the FIXED kine_ref v2 (frozen-leg
hole closed). Runs each trial's eval, then scores the npz the runner
just wrote (window 2.0 = the runner's eval walk window)."""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
import numpy as np

import _curriculum as C
import kine_ref as KR
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5

s = optuna.load_study(study_name="curr_s3c_ground",
                      storage="sqlite:///optuna_walk.db")
done = [(t.value, t.number, t.params) for t in s.trials
        if t.value is not None]
done.sort(key=lambda x: -x[0])
rows = []
for old_v, num, params in done[:8]:
    try:
        C.set_stage(3, {**mul, **params})
        import os
        os.environ["AARL_NPZ"] = "spinal_run_rescore.npz"
        m = R.main(["--eval", "--drive", repr(params["drive"])])
        z = np.load("spinal_run_rescore.npz", allow_pickle=True)
        k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                       ref=KR.ref_cached(), contact=z["contact"])
        z.close()   # NpzFile holds the zip open -> blocks the next
        # trial's os.replace (WinError 5)
    except Exception as e:
        print(f"t{num:3d} FAILED: {e}", flush=True)
        rows.append((float("-inf"), old_v, num))
        continue
    if k is None:
        new_v = -320.0
        k = {}
    else:
        new_v = float(k["kine_score"])

    def g(key, fmt="{:.2f}"):
        v = k.get(key)
        return fmt.format(v) if isinstance(v, (int, float)) else "--"
    rows.append((new_v, old_v, num))
    print(f"t{num:3d} old {old_v:8.1f} -> new {new_v:8.1f} "
          f"duty r/l {g('duty_r')}/{g('duty_l')} "
          f"cyc r/l {g('n_cycles_r', '{:d}')}/{g('n_cycles_l', '{:d}')} "
          f"frozen {g('frozen_r', '{!s:}')}/{g('frozen_l', '{!s:}')}",
          flush=True)
rows.sort(key=lambda x: -x[0])
print("\nCORRECTED ranking:")
for r in rows:
    print(f"  t{r[2]:3d} new {r[0]:8.1f} (was {r[1]:.1f})")
with open("_rescore_top_out.json", "w", encoding="utf-8") as f:
    json.dump([{"trial": r[2], "new": r[0], "old": r[1]} for r in rows],
              f, indent=1)
print(f"best under fixed objective: trial {rows[0][2]} ({rows[0][0]:.1f})")
