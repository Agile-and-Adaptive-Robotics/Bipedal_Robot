"""Audit whether E1/E2/F1/F2 form four identifiable temporal basis columns."""
from pathlib import Path

import numpy as np

from fit_pf import pf_windows_from_run


HERE = Path(__file__).parent
PHASES = ("E1", "E2", "F1", "F2")


def main():
    run = np.load(HERE / "spinal_run.npz", allow_pickle=True)
    win = pf_windows_from_run(run, "r", n_bins=20)
    basis = np.column_stack([win[p] for p in PHASES])
    corr = np.corrcoef(basis, rowvar=False)
    singular = np.linalg.svd(basis, compute_uv=False)
    energy = singular ** 2 / np.sum(singular ** 2)
    rank = np.linalg.matrix_rank(basis)
    condition = np.linalg.cond(basis)

    print("PF basis phase summaries (20 bins):")
    for i, phase in enumerate(PHASES):
        v = basis[:, i]
        peak = int(np.argmax(v))
        com = float(np.sum((np.arange(20) + 0.5) / 20 * v) /
                    max(np.sum(v), 1e-12))
        duty = float(np.mean(v >= 0.5 * np.max(v)))
        print(f"  {phase}: peak={peak / 20:.2f} cycle, "
              f"center={com:.3f}, half-max duty={duty:.2f}")

    print("\nPairwise correlation:")
    print("       " + " ".join(f"{p:>7s}" for p in PHASES))
    for i, phase in enumerate(PHASES):
        print(f"  {phase:>3s}  " + " ".join(f"{x:7.3f}" for x in corr[i]))

    print("\nLinear-basis diagnostics:")
    print("  singular values:", " ".join(f"{x:.4f}" for x in singular))
    print("  energy fractions:", " ".join(f"{x:.3f}" for x in energy))
    print(f"  numerical rank: {rank}/4")
    print(f"  condition number: {condition:.2f}")


if __name__ == "__main__":
    main()
