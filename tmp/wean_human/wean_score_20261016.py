"""Score a walker npz against the reference set the -189.5 of-record
number came from (kine_ref.compare on the subject01 reference of record,
walk window 2.0 s for EVAL-mode npz — gait_lib_score_all.py pattern) plus
chosen gait-library refs (Ong walking speeds).

Usage: python wean_score_20261016.py <npz> <walk_start> <ref-stem...>
ref-stem "subject01" = the of-record reference (kine_ref.load_reference);
anything else = gait_refs/<stem>.npz (gait_lib_score_all.load_npz_ref).
"""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import kine_ref as KR


def load_npz_ref(p: Path):
    z = np.load(p, allow_pickle=True)
    ref = {s: {j: z[f"{s}_{j}"] for j in ("hip", "knee", "ankle")}
           for s in ("r", "l")}
    for k in z.files:
        if k not in (f"{s}_{j}" for s in ("r", "l")
                     for j in ("hip", "knee", "ankle")):
            v = z[k]
            ref[k] = v.item() if v.shape == () else v
    return ref


def main():
    npz = sys.argv[1]
    ws = float(sys.argv[2])
    stems = sys.argv[3:] or ["subject01"]
    z = np.load(npz, allow_pickle=True)
    t, q, neuro, contact = z["t"], z["q"], z["neuro"], z["contact"]
    n_done = int(np.sum(t > 0)) + 1
    for stem in stems:
        if stem == "subject01":
            ref = KR.load_reference()
        else:
            ref = load_npz_ref(HERE / "gait_refs" / f"{stem}.npz")
        s = KR.compare(t[:n_done], q[:n_done], neuro[:n_done], ws,
                       ref=ref, contact=contact[:n_done])
        if s is None:
            print(f"REFSCORE {stem} kine_score=None")
            continue
        scal = {k: v for k, v in s.items()
                if isinstance(v, (int, float, bool))}
        print(f"REFSCORE {stem} " + json.dumps(scal))


if __name__ == "__main__":
    main()
