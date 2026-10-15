"""Replay a study's best trial and save the npz run artifact (for
figures/gifs). Reconstructs each curriculum's exact param merge.

Usage:
  python _replay_best_20261003.py <model> <db> <study> <npz_out> [stage]
    model: s3k | syn6 | w2lvar | pruned | spiking
"""
import io
import json
import os
import sys


class _NoCloseRaw(io.RawIOBase):
    """Raw facade whose close() never closes the real stdout - the
    curriculum modules each wrap sys.stdout at import; when an
    intermediate TextIOWrapper is garbage-collected its close() would
    otherwise close the shared buffer (the 2026-10-03 replay crash)."""

    def __init__(self, raw):
        super().__init__()
        self._raw = raw

    def write(self, b):
        return self._raw.write(b)

    def writable(self):
        return True

    def flush(self):
        return self._raw.flush()

    def close(self):
        try:
            self.flush()
        except Exception:
            pass
        # deliberately does NOT close self._raw


sys.stdout = io.TextIOWrapper(_NoCloseRaw(sys.stdout.buffer),
                              encoding="utf-8", errors="replace")

model, db, study, npz_out = sys.argv[1:5]
stage = int(sys.argv[5]) if len(sys.argv) > 5 else 4

import optuna

st = optuna.load_study(study_name=study, storage="sqlite:///" + db)
done = [t for t in st.trials if t.state.name == "COMPLETE" and
        t.value is not None]
best = max(done, key=lambda t: t.value)
print(f"[replay] {study}: best {best.value:.4f} (trial {best.number})")
print(json.dumps(best.params, indent=1))

os.environ["AARL_NPZ"] = npz_out

if model == "s3k":
    os.environ.pop("AARL_NET", None)
    import _curriculum as C
    import runner as R

    base = dict(json.loads(open("best_walk_params_v10.json",
                                encoding="utf-8").read())["multipliers"])
    base["renshaw"] = 0.5
    if stage in (4, 5):
        s3k = json.loads(open(
            "reports_20260923/s3k_trial34_full_params.json",
            encoding="utf-8").read())["params"]
        base.update({k: v for k, v in s3k.items()})
        if s3k.get("pf_gain") is None:
            base.pop("pf_gain", None)
    p = {**base, **best.params}
    C.set_stage(stage, p)
    m = R.main(["--eval", "--drive", repr(p["drive"])])
    print("[replay] kine", m.get("kine"), "score", m.get("kine_score"))

elif model in ("syn6", "w2lvar"):
    mod = __import__(f"_curriculum_{model}")
    base = dict(json.loads(open("best_walk_params_v10.json",
                                encoding="utf-8").read())["multipliers"])
    if model == "syn6":
        base["renshaw"] = 0.0
        base["syn6"] = 1.0
        base["syn6_brainstem"] = 0.0
    else:
        base["renshaw"] = 0.5
    os.environ["AARL_NET"] = model
    import runner as R

    p = {**base, **best.params}
    mod.set_stage(stage, p)
    m = R.main(["--eval", "--drive", repr(p["drive"])])
    print("[replay] kine", m.get("kine"), "score", m.get("kine_score"))

elif model == "spiking":
    os.environ["AARL_NET"] = "spiking"
    import _curriculum_spiking as CS
    import runner as R

    # the module pins AARL_NPZ at import; re-pin OUR output name
    os.environ["AARL_NPZ"] = npz_out
    base = dict(json.loads(open("best_walk_params_v10.json",
                                encoding="utf-8").read())["multipliers"])
    base["renshaw"] = 0.5
    p = {**base, **best.params}
    CS.CS3.set_stage(stage, p)
    if stage <= 2:
        args = ["--no-ground", "--time", "20",
                "--drive", repr(p["drive"])]
        if stage == 1:
            args += ["--no-afferents", "--no-interleg"]
    else:
        args = ["--eval", "--drive", repr(p["drive"])]
    m = R.main(args)
    print("[replay] kine", m.get("kine"), "score", m.get("kine_score"))

elif model == "pruned":
    os.environ.pop("AARL_NET", None)
    import _curriculum_pruned as CP
    import runner as R

    # the module pins AARL_NPZ at import; re-pin OUR output name
    os.environ["AARL_NPZ"] = npz_out
    argv = CP._config(best.params)
    m = R.main(["--eval", "--drive",
                repr(float(best.params["drive"]))] + argv)
    print("[replay] kine", m.get("kine"), "score", m.get("kine_score"))

else:
    raise SystemExit(f"unknown model {model}")
