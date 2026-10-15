"""Audit the patched MuJoCo plant that runner.py actually simulates.

Unlike _fmax_audit.py (which intentionally exposes raw converter
artifacts), this check compiles runner.apply_harness(), verifies the eight
1-newton actuator repairs, and summarizes agreement with stock OpenSim.
"""
from __future__ import annotations

import mujoco
import numpy as np

import runner as R
from _fmax_audit import OSIM, PRUNE, parse_osim_fmax

def main() -> None:
    raw = mujoco.MjModel.from_xml_path(str(R.MODEL))
    live = R.apply_harness(raw, mujoco.MjData(raw))
    names = [mujoco.mj_id2name(live, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
             for i in range(live.nu)]
    stock = parse_osim_fmax(OSIM)
    fmax = dict(zip(names, live.actuator_gainprm[:, 2]))

    expected = {"ercspn": 2500.0, "intobl": 900.0,
                "extobl": 900.0, "ext_hal": 162.0}
    for base, target in expected.items():
        for side in ("r", "l"):
            name = f"{base}_{side}"
            value = float(fmax[name])
            assert np.isclose(value, target), (name, value, target)
            print(f"{name:10s} {value:7.1f} N  repaired OK")

    active = [n for n in names if n not in PRUNE]
    relerr = {n: abs(float(fmax[n]) - stock[n]) / stock[n] for n in active}
    within15 = [n for n in active if relerr[n] <= 0.15]
    outside15 = sorted((n, relerr[n]) for n in active if relerr[n] > 0.15)
    print(f"active actuators within 15% of stock: {len(within15)}/{len(active)}")
    for name, err in outside15:
        print(f"  outside 15%: {name:12s} {100 * err:5.1f}%")
    print("LIVE PLANT REPAIR AUDIT PASS; stock-tolerance exceptions reported")


if __name__ == "__main__":
    main()
