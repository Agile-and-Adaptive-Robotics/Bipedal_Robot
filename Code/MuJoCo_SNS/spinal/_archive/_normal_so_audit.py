"""Quantify the explicitly residual-supported normal.mot OpenSim solve."""
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).parents[3]
OUT = (ROOT / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody" /
       "ResultsNormalSO")
STO = OUT / "normal_backsolve_StaticOptimization_activation.sto"


def read_sto(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    end = next(i for i, line in enumerate(lines) if line.strip() == "endheader")
    names = lines[end + 1].split()[1:]
    rows = np.array([[float(x) for x in line.split()]
                     for line in lines[end + 2:] if line.split()])
    return rows[:, 0], names, rows[:, 1:]


def main():
    time, names, values = read_sto(STO)
    residuals = [n for n in names if n.startswith("normal_residual_")]
    optimal = dict(normal_residual_FX=4.0, normal_residual_FY=8.0,
                   normal_residual_FZ=4.0, normal_residual_MX=2.0,
                   normal_residual_MY=2.0, normal_residual_MZ=2.0)
    stats = {}
    for name in residuals:
        control = values[:, names.index(name)]
        generalized_force = control * optimal[name]
        stats[name] = {
            "control_rms": float(np.sqrt(np.mean(control ** 2))),
            "force_or_torque_mean": float(np.mean(generalized_force)),
            "force_or_torque_rms": float(np.sqrt(np.mean(generalized_force ** 2))),
            "force_or_torque_peak_abs": float(np.max(np.abs(generalized_force))),
        }
    muscle_idx = [i for i, name in enumerate(names) if name not in residuals]
    muscle = values[:, muscle_idx]
    payload = {
        "frames": len(time), "time_s": [float(time[0]), float(time[-1])],
        "model": "gait2392_simbody_normal_so.osim (OpenSim 4.6 upgrade; lat_gas_r disabled)",
        "motion": "Tutorial1/normal.mot", "external_loads": "none",
        "residuals": stats,
        "muscle_activation_min": float(muscle.min()),
        "muscle_activation_max": float(muscle.max()),
        "saturated_entries_ge_0p99": int(np.count_nonzero(muscle >= 0.99)),
    }
    (OUT / "normal_so_audit.json").write_text(json.dumps(payload, indent=2),
                                                encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
