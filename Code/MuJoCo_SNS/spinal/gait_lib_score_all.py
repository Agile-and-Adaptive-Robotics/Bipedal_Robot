"""Score the s3k production walker against EVERY integrated gait-library
reference (2026-09-24 thumb-drive campaign).

Loads all gait_refs/*.npz (33 refs: Falisse predicted walking + 8 Arnold
running subjects x 4 speeds), rebuilds each as kine_ref.REF_CACHE, and
scores the recorded s3k run by kine_ref.compare(). Read-only; writes
goal5_allrefs_scoring.{md,csv} to reports_20260923.
"""
import csv
import io
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import kine_ref as KR

NPZ = HERE / "spinal_run_s3k.npz"
WALK_START = 2.0
REFS = sorted((HERE / "gait_refs").glob("*.npz"))


def load_npz_ref(p: Path):
    z = np.load(p, allow_pickle=True)
    ref = {s: {j: z[f"{s}_{j}"] for j in ("hip", "knee", "ankle")}
           for s in ("r", "l")}
    for k in z.files:
        if "_" not in k[:2] or k.startswith(("duty_", "T_", "mean_", "ds",
                                             "hip_", "knee_", "ankle_")):
            if k not in ("r_hip",):
                v = z[k]
                ref[k] = v.item() if v.shape == () else v
    # keep only scalar/array keys that are NOT side-joint arrays
    for s in ("r", "l"):
        for j in ("hip", "knee", "ankle"):
            ref.pop(f"{s}_{j}", None)
    return ref


def main():
    z = np.load(NPZ, allow_pickle=True)
    t, q, neuro, contact = z["t"], z["q"], z["neuro"], z["contact"]

    # reference of record for the baseline row
    rows = []
    s0 = KR.compare(t, q, neuro, WALK_START, ref=None, contact=contact)
    if s0:
        rows.append(("subject01_of_record",
                     float(s0.get("kine_score", np.nan)),
                     s0.get("n_cycles_r", 0), s0.get("n_cycles_l", 0),
                     s0.get("duty_r", np.nan), s0.get("T_r", np.nan),
                     s0.get("knee_min_r", np.nan)))

    for p in REFS:
        try:
            ref = load_npz_ref(p)
        except Exception as e:
            print(f"{p.stem}: ref load FAILED {e}")
            continue
        old = KR.REF_CACHE
        KR.REF_CACHE = ref
        try:
            s = KR.compare(t, q, neuro, WALK_START, ref=ref,
                           contact=contact)
        except Exception as e:
            print(f"{p.stem}: compare FAILED {e}")
            KR.REF_CACHE = old
            continue
        KR.REF_CACHE = old
        if s is None:
            print(f"{p.stem}: no cycles")
            continue
        ks = s.get("kine_score")
        print(f"{p.stem}: kine_score {ks}")
        rows.append((p.stem,
                     float(ks) if ks is not None else np.nan,
                     s.get("n_cycles_r", 0), s.get("n_cycles_l", 0),
                     s.get("duty_r", np.nan), s.get("T_r", np.nan),
                     s.get("knee_min_r", np.nan)))

    csv_path = HERE / "reports_20260923" / "goal5_allrefs_results.csv"
    csv_path.parent.mkdir(exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ref", "kine_score", "cycles_r", "cycles_l",
                    "sim_duty_r", "sim_T_r", "sim_knee_min_r",
                    "ref_T_r", "ref_duty_r", "ref_knee_min_r"])
        zr = {p.stem: p for p in REFS}
        for name, ks, cr, cl, dr, tr, km in rows:
            rT = rD = rK = ""
            if name in zr:
                zz = np.load(zr[name], allow_pickle=True)
                rT = "%.3f" % float(zz["T_r"])
                rD = "%.2f" % float(zz["duty_r"])
                rK = "%.1f" % float(zz["knee_min_r"])
            w.writerow([name, "%.2f" % ks if ks == ks else "", cr, cl,
                        "%.2f" % dr if dr == dr else "",
                        "%.3f" % tr if tr == tr else "",
                        "%.1f" % km if km == km else "", rT, rD, rK])

    lines = ["# s3k production walker vs the full thumb-drive library "
             "(2026-09-24)", "",
             f"Run: `{NPZ.name}`, walk window from t={WALK_START} s. "
             f"{len(rows)} scored references (lower kine_score = better).",
             "",
             "| ref | kine_score | sim duty/T/knee | ref duty/T/knee |",
             "|---|---|---|---|"]
    for name, ks, cr, cl, dr, tr, km in rows:
        lines.append(f"| {name} | {ks:.1f} | "
                     f"{dr:.2f} / {tr:.2f}s / {km:.0f}deg | |")
    lines += ["",
              "All 32 Arnold trials are RUNNING (duty 0.31-0.46, knee_min "
              "-80..-140 deg across 4 speeds x 8 subjects); Falisse "
              "Case_40 is predicted WALKING. The s3k walker is a WALKER "
              "(duty ~0.69, knee ~-53): it should match walking refs far "
              "better than running refs - the spread quantifies gait-"
              "specificity of the current tuning.",
              "",
              "Training-side use: swap kine_ref.REF_CACHE per trial "
              "(gait_lib_pilot2.py pattern) to build a multi-reference / "
              "multi-gait curriculum."]
    out = HERE / "reports_20260923" / "goal5_allrefs_scoring.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("saved", out, "and", csv_path)


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
