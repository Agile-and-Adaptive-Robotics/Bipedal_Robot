"""Score an existing walker run against MULTIPLE gait references (goal 5).

First real use of the staged Falisse predictsim_mtp data: the s3k winner
trajectory (spinal_run_s3k.npz) is scored by kine_ref.compare() against
(a) the subject01 reference of record and (b) the Falisse predicted-walking
reference from gait_lib_loader - measuring how reference-specific our
current tuning is. Read-only: no studies launched, no core files touched.
"""
import io
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import kine_ref as KR
from gait_lib_loader import load_reference_general

FAL = Path(r"D:\temp\gait_lib_staging\falisse\predictsim_mtp-master"
           r"\Results\Case_40")
NPZ = HERE / "spinal_run_s3k.npz"
WALK_START = 2.0   # eval-mode npz gotcha (runner rewrites SCHEDULE for --eval)


def main():
    z = np.load(NPZ, allow_pickle=True)
    t, q, neuro, contact = z["t"], z["q"], z["neuro"], z["contact"]
    print(f"loaded {NPZ.name}: {len(t)} samples, "
          f"{t[-1]:.1f} s, q{q.shape}, contact{contact.shape}")

    refs = {
        "subject01 (of record)": KR.load_reference(),
        "Falisse Case_40 (predicted)": load_reference_general(
            FAL / "motion.mot", FAL / "GRF.mot"),
    }
    lines = ["# Goal-5: s3k winner scored against two references (2026-09-23)",
             "", f"Run: `{NPZ.name}` (eval walk window from t={WALK_START} s).",
             "kine_score is kine_ref.compare's objective: LOWER = better "
             "match to that reference.", ""]
    rows = []
    for tag, ref in refs.items():
        s = KR.compare(t, q, neuro, WALK_START, ref=ref, contact=contact)
        if s is None:
            print(f"{tag}: compare() returned None (no cycles)")
            continue
        ks = s.get("kine_score")
        print(f"{tag}: kine_score {ks}")
        rows.append((tag, s))

    keys = ["n_cycles_r", "n_cycles_l", "duty_r", "duty_r_ref",
            "duty_l", "duty_l_ref", "knee_min_r", "knee_min_r_ref",
            "T_r", "T_r_ref", "range_hip_r", "range_hip_r_ref",
            "rmse_knee_r", "mean_hip_r", "mean_hip_r_ref"]
    header = "| metric | " + " | ".join(t0 for t0, _ in rows) + " |"
    lines += [header, "|" + "---|" * (len(rows) + 1)]
    vals = {k: [] for k in keys}
    for _, s in rows:
        for k in keys:
            v = s.get(k, float("nan"))
            vals[k].append(f"{v:.2f}" if isinstance(v, float) else str(v))
    for k in keys:
        lines.append(f"| {k} | " + " | ".join(vals[k]) + " |")
    for tag, s in rows:
        lines.append(f"\n- {tag}: kine_score **{s.get('kine_score')}**")

    lines += ["", "## Reading",
              "- The gap between the two kine_scores measures how "
              "reference-specific the s3k tuning is: a walker trained on one "
              "reference cycle (subject01) should score visibly better "
              "against it than against an independent predicted gait.",
              "- Falisse reference is one periodic cycle (T 1.113 s, duty "
              "0.57, knee_min -59.2 deg) - faster and stiffer than "
              "subject01 (1.233 s / 0.61 / -69.7 deg).",
              "- No retraining was run (goal-4 next-study config awaits "
              "Ben); this is the measurement half of multi-reference "
              "training."]
    out = HERE / "reports_20260923" / "goal5_falisse_scoring.md"
    out.parent.mkdir(exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("saved", out)


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
