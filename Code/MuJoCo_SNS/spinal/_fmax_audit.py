"""Fmax audit: converted MJCF actuator gainprm[2] vs STOCK OpenSim
gait2392_thelen max_isometric_force.

Answers (Ben, 2026-09-13): "Does OpenSim have that as the Fmax? No,
it's ridiculous." -> the 1-newton trunk muscles are a MyoConverter
artifact, NOT OpenSim values. This script lists EVERY mismatch and
dumps the gainprm/biasprm layout of a broken actuator vs a good one.

Usage: python _fmax_audit.py
"""
import io
import re
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import mujoco

import runner as R

OSIM = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
            r"\Gait2392_Robotbody\gait2392_simbody.osim")
PRUNE = {"quad_fem_r", "quad_fem_l", "gem_r", "gem_l", "peri_r", "peri_l"}


def parse_osim_fmax(path: Path) -> dict[str, float]:
    txt = path.read_text(encoding="utf-8", errors="replace")
    out = {}
    for m in re.finditer(
            r"<(Thelen2003Muscle|Millard2012EquilibriumMuscle)\s+name=\"([^\"]+)\"(.*?)</\1>",
            txt, flags=re.S):
        name, body = m.group(2), m.group(3)
        f = re.search(r"<max_isometric_force>\s*([0-9.eE+-]+)\s*</max_isometric_force>",
                      body)
        if f:
            out[name] = float(f.group(1))
    return out


def main():
    osim = parse_osim_fmax(OSIM)
    print(f"stock OpenSim muscles parsed: {len(osim)} "
          f"({OSIM.name})")
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    names = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
             for i in range(model.nu)]
    fmax_mj = model.actuator_gainprm[:, 2].copy()
    print(f"MJCF actuators: {model.nu}, total body mass "
          f"{model.body_mass.sum():.2f} kg")

    # muscle-layout reference: a healthy actuator vs a broken one
    for probe in ("soleus_r", "ercspn_r"):
        i = names.index(probe)
        print(f"{probe:10s} gainprm={np.array2string(model.actuator_gainprm[i], precision=3)}")
        print(f"{'':10s} biasprm={np.array2string(model.actuator_biasprm[i], precision=3)}"
              f"  lengthrange={np.array2string(model.actuator_lengthrange[i], precision=3)}")

    bad, prune, ok = [], [], 0
    for i, nm in enumerate(names):
        ref = osim.get(nm)
        if ref is None:
            print(f"  ?? {nm}: no stock OpenSim muscle found")
            continue
        if nm in PRUNE:
            prune.append((nm, fmax_mj[i], ref))
            continue
        if abs(fmax_mj[i] - ref) > 0.01 * max(ref, 1.0):
            bad.append((nm, fmax_mj[i], ref))
        else:
            ok += 1
    print(f"\nMATCHING: {ok}/{model.nu} within 1%")
    print(f"INTENTIONAL PRUNES (Fmax zeroed in patch_xml, Ben's list): "
          f"{len(prune)}")
    for nm, mj, ref in sorted(prune):
        print(f"  {nm:10s} mj {mj:7.1f} N   stock {ref:7.1f} N")
    print(f"MISMATCHES (converter artifact): {len(bad)}")
    for nm, mj, ref in sorted(bad):
        print(f"  {nm:10s} mj {mj:7.2f} N   stock {ref:7.1f} N "
              f"({ref / max(mj, 0.01):,.0f}x too small)")
    tot_bad_ref = sum(ref for _, _, ref in bad)
    print(f"stock total force capacity of the mismatched set: "
          f"{tot_bad_ref:,.0f} N")
    return bad


if __name__ == "__main__":
    main()
