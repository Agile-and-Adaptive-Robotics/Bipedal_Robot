"""Literature validation, all walkers (2026-09-30): record a WALK run
for each variant (s3k, s3kpruned combo, w2lvar, syn6) through its
authoritative config, then score it against EVERY gait_refs reference
via kine_ref.compare. Writes gait_validation_20260930.{md,csv}.

Usage: python gait_validate_all.py [variant ...]   (default: all four)
"""
import io
import json
import os
import sys
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import numpy as np

import kine_ref as KR
import params as P
import runner as R


def _jload(fn):
    with open(HERE / fn, encoding="utf-8") as f:
        return json.load(f)


BASE = dict(_jload("best_walk_params_v10.json")["multipliers"])
PRUNE_KEYS = {"heel_rge": 0.0, "toe_rge": 0.0, "contact_onset": 0.0,
              "contra_swing": 0.0, "pm_gain": 0.0, "pm_aff": 0.0,
              "ib_rge": 0.0, "ib_e_central": 0.0, "rg_weak_exc": 0.0}


def config_variant(v):
    os.environ.pop("AARL_NET", None)
    if v == "s3k":
        import _curriculum as CS3
        W = {**BASE, **_jload(
            "reports_20260923/s3k_trial34_full_params.json")["params"]}
        CS3.set_stage(5, dict(W))
        return repr(float(W["drive"])), []
    if v == "s3kpruned":
        import _curriculum as CS3
        W = {**BASE, **_jload(
            "reports_20260923/s3k_trial34_full_params.json")["params"]}
        CS3.set_stage(5, dict(W))
        for k, vv in PRUNE_KEYS.items():
            P.G[k] = float(vv)
        return repr(float(W["drive"])), ["--no-interleg"]
    if v == "w2lvar":
        os.environ["AARL_NET"] = "w2lvar"
        import _curriculum_w2lvar as CV
        W = {**BASE, "renshaw": 0.5, **_jload(
            "curriculum_w2lvar_stage5.json")["params"]}
        CV.set_stage(5, dict(W))
        return repr(float(W["drive"])), []
    if v == "syn6":
        os.environ["AARL_NET"] = "syn6"
        import _curriculum_syn6 as CS
        W = {**BASE, "renshaw": 0.0, "syn6": 1.0, "syn6_brainstem": 0.0,
             **_jload("curriculum_syn6_stage5.json")["params"]}
        CS.set_stage(5, dict(W))
        return repr(float(W["drive"])), []


def load_npz_ref(p: Path):
    z = np.load(p, allow_pickle=True)
    ref = {s: {j: z[f"{s}_{j}"] for j in ("hip", "knee", "ankle")}
           for s in ("r", "l")}
    for k in z.files:
        if "_" not in k[:2] or k.startswith(("duty_", "T_", "mean_", "ds",
                                             "hip_", "knee_", "ankle_")):
            v = z[k]
            ref[k] = v.item() if v.shape == () else v
    for s in ("r", "l"):
        for j in ("hip", "knee", "ankle"):
            ref.pop(f"{s}_{j}", None)
    return ref


def _ensure_stdout():
    """The curriculum loaders each wrap sys.stdout; after several R.main
    calls an inner wrapper gets GC'd and closes the stream ("I/O
    operation on closed file" - killed the syn6 pass 2026-09-30).
    Re-open fd 1 directly if the current stream is dead."""
    try:
        sys.stdout.write("")
    except Exception:
        sys.stdout = io.TextIOWrapper(open(1, "wb", closefd=False),
                                      encoding="utf-8",
                                      errors="replace")


def main():
    variants = sys.argv[1:] or ["s3k", "s3kpruned", "w2lvar", "syn6"]
    refs = sorted((HERE / "gait_refs").glob("*.npz"))
    import csv
    csv_path = HERE / "gait_validation_20260930.csv"
    new_file = not csv_path.exists()
    f_csv = open(csv_path, "a", encoding="utf-8", newline="")
    w = csv.writer(f_csv)
    if new_file:
        w.writerow(["variant", "reference", "kine_score"])
    for v in variants:
        _ensure_stdout()
        # crash-safety: skip variants already fully scored in the csv
        done_refs = set()
        if not new_file:
            with open(csv_path, encoding="utf-8") as fr:
                for row in csv.reader(fr):
                    if row and row[0] == v:
                        done_refs.add(row[1])
        if len(done_refs) >= len(refs) + 1:
            print(f"[val:{v}] already scored, skipping", flush=True)
            continue
        npz = HERE / f"spinal_run_val_{v}.npz"
        os.environ["AARL_NPZ"] = npz.name
        drive, argv = config_variant(v)
        _ensure_stdout()
        print(f"[val:{v}] running walk eval...", flush=True)
        m = R.main(["--eval", "--drive", drive] + argv) or {}
        _ensure_stdout()
        z = np.load(npz, allow_pickle=True)
        t, q, neuro, contact = z["t"], z["q"], z["neuro"], z["contact"]
        w.writerow([v, "kine_score",
                    float(m.get("kine_score") or -25.0)])
        f_csv.flush()
        for p in refs:
            try:
                ref = load_npz_ref(p)
                s = KR.compare(t, np.degrees(q), neuro, 2.0, ref=ref,
                               contact=contact)
                w.writerow([v, p.stem,
                            float(s["kine_score"]) if s else None])
            except Exception as ex:
                w.writerow([v, p.stem, None])
                _ensure_stdout()
                print(f"[val:{v}] {p.stem}: {ex!r}", flush=True)
            f_csv.flush()
        _ensure_stdout()
        print(f"[val:{v}] scored vs {len(refs)} refs", flush=True)
    f_csv.close()
    _ensure_stdout()
    print(f"wrote {csv_path.name}", flush=True)


if __name__ == "__main__":
    main()
