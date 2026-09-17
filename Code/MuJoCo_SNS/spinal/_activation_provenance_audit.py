"""Identify what the ResultsBSolve activation STO actually contains.

Compare it by muscle name against both arrays saved in bsolve_out.npz:
``acts`` is the converted-MuJoCo ridge/NNLS solution and ``so_act`` is the
OpenSim StaticOptimization result parsed by bsolve_ik.py.  This avoids
inferring provenance from the ambiguous ``StaticOptimization`` filename.
"""
from pathlib import Path

import numpy as np


HERE = Path(__file__).parent
STO = (HERE.parents[2] / "Solid_Models" / "OpenSim" /
       "Gait2392_Robotbody" / "ResultsBSolve" /
       "zz_bsolve_StaticOptimization_activation.sto")


def read_sto(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    end = next(i for i, line in enumerate(lines)
               if line.strip() == "endheader")
    names = lines[end + 1].split()[1:]
    rows = np.array([[float(x) for x in line.split()]
                     for line in lines[end + 2:] if line.split()])
    return rows[:, 0], names, rows[:, 1:]


def compare(label, ref, ref_names, sto, sto_names):
    def text(name):
        return name.decode("utf-8") if isinstance(name, bytes) else str(name)
    ref_col = {text(name): i for i, name in enumerate(ref_names)}
    pairs = [(j, ref_col[name]) for j, name in enumerate(sto_names)
             if name in ref_col]
    if not pairs:
        print(f"{label}: no muscle-name matches; first saved names="
              f"{list(ref_col)[:5]}")
        return
    n = min(len(sto), len(ref))
    x = np.column_stack([sto[:n, j] for j, _ in pairs])
    y = np.column_stack([ref[:n, i] for _, i in pairs])
    err = x - y
    corr = np.corrcoef(x.ravel(), y.ravel())[0, 1]
    print(f"{label}: frames={n}, matched={len(pairs)}, "
          f"max_abs={np.max(np.abs(err)):.9g}, "
          f"rmse={np.sqrt(np.mean(err ** 2)):.9g}, corr={corr:.9f}")


def main():
    time, sto_names, sto = read_sto(STO)
    data = np.load(HERE / "bsolve_out.npz", allow_pickle=True)
    print(f"STO: {STO.name}, {sto.shape}, t={time[0]:.3f}-{time[-1]:.3f}")
    compare("MuJoCo ridge/NNLS (bsolve_out acts)", data["acts"],
            data["act_names"], sto, sto_names)
    so_names = [x.decode("utf-8") if isinstance(x, bytes) else str(x)
                for x in data["so_names"]]
    # Parsed OpenSim arrays include their time column name/data in some runs.
    so = data["so_act"]
    if so_names and so_names[0] == "time" and so.shape[1] == len(so_names):
        so_names, so = so_names[1:], so[:, 1:]
    compare("OpenSim SO (bsolve_out so_act)", so, so_names,
            sto, sto_names)


if __name__ == "__main__":
    main()
