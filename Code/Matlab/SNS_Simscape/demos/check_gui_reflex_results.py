"""Check fresh single-joint results against the saved response and controls."""
from argparse import ArgumentParser
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat


def read(path):
    run = loadmat(path, simplify_cells=True)
    data = dict(zip(run["signals"], run["samples"].T))
    for value in data.values():
        assert np.isfinite(value).all(), path
    return run, data


def knee(folder, candidate, baseline):
    a, p = read(folder / (candidate + ".mat"))
    b, c = read(folder / (baseline + ".mat"))
    assert np.array_equal(a["time"], b["time"])
    t = a["time"] - a["onset"]
    d = {name: p[name] - c[name] for name in p}
    extension = -np.rad2deg(d["th"])
    peak = int(np.argmax(extension))
    return {
        "peak_extension_deg": float(extension[peak]),
        "time_to_peak_s": float(t[peak]),
        "max_absolute_angle_difference_deg": float(np.abs(extension).max()),
        "force_difference_at_30ms_N": [float(np.interp(.030, t, d[k])) for k in ("F_ext", "F_flex")],
        "voltage_difference_at_30ms_mV": [float(np.interp(.030, t, d[k])) for k in ("V_MN_ext", "V_MN_flex")],
        "max_spindle_current_increase_first_50ms_nA": float(np.max(d["Ia_ext"][(t >= 0) & (t <= .05)])),
        "max_fiber_elongation_first_50ms_mm": float(1000*np.max(d["Lfiber_ext"][(t >= 0) & (t <= .05)])),
    }


def arm(folder, enabled):
    src, reference = read(folder / f"arm_source_{enabled}.mat")
    gui, candidate = read(folder / f"arm_gui_{enabled}.mat")
    assert np.array_equal(src["time"], gui["time"])
    assert set(reference) == set(candidate)
    differences = {key: float(np.abs(candidate[key]-reference[key]).max()) for key in candidate}
    return {"max_abs_difference": differences, "final_mass_kg": float(candidate["mCup"][-1]),
            "final_sag_deg": float(np.rad2deg(candidate["th"][-1]))}


def check(folder):
    k = folder/"knee"
    records = {case: knee(k, case, baseline) for case, baseline in [
        ("ia_pulse", "control"), ("half_ia", "control"), ("fine_ia", "fine_control"),
        ("off_ia", "off_control"), ("torque_pulse", "control"), ("off_torque", "off_control")]}
    full = records["ia_pulse"]
    assert abs(full["peak_extension_deg"]-0.93605) < 0.001
    assert abs(full["time_to_peak_s"]-.167) < .001
    assert all(np.array(full["force_difference_at_30ms_N"])*[1.,-1.] > 0)
    assert all(np.array(full["voltage_difference_at_30ms_mV"])*[1.,-1.] > 0)
    assert records["off_ia"]["max_absolute_angle_difference_deg"] < 1e-5
    assert abs(records["half_ia"]["peak_extension_deg"]-.46766) < .001
    assert abs(records["fine_ia"]["peak_extension_deg"]-full["peak_extension_deg"]) < .0001
    assert records["torque_pulse"]["max_spindle_current_increase_first_50ms_nA"] > .12
    assert records["torque_pulse"]["max_fiber_elongation_first_50ms_mm"] > .11
    arms = {str(g): arm(folder/"arm",g) for g in (0,1)}
    for checks in arms.values():
        assert checks["final_mass_kg"] == .625
        assert max(checks["max_abs_difference"].values()) < 1e-6
    assert abs(arms["0"]["final_sag_deg"]-7.84641) < .001
    assert abs(arms["1"]["final_sag_deg"]-1.99654) < .001
    result = {"knee": records, "arm": arms}
    (folder/"gui_reflex_check.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(f"Six biological comparisons and all nine arm signals passed: {folder}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    check(parser.parse_args().folder)
