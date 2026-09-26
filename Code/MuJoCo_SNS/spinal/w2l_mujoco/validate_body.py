# -*- coding: utf-8 -*-
"""
validate_body.py - Milestone 1 validation gate for the AnimatLab W2L -> MuJoCo
body port. Run AFTER make_w2l_mjcf.py. Read-only w.r.t. the source .aproj.

Gates:
  1. model loads under mujoco 2.3.7
  2. 2 s passive drop: no NaN in qpos/qvel, bounded state
  3. all 12 muscle actuators have finite resting (operational) lengths
  4. joint ranges match the aproj within degree rounding
Plus: body-count/mass fidelity readback vs w2l_source_dump.json.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
XML = os.path.join(HERE, "w2l_mjcf.xml")
DUMP = os.path.join(HERE, "w2l_source_dump.json")

os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import mujoco  # noqa: E402
import numpy as np  # noqa: E402

PASS = []
FAIL = []


def check(name, ok, detail=""):
    (PASS if ok else FAIL).append((name, detail))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  " + detail) if detail else ""))


def main():
    with open(DUMP) as fh:
        dump = json.load(fh)

    print("== gate 1: model loads under mujoco %s ==" % mujoco.__version__)
    model = mujoco.MjModel.from_xml_path(XML)
    data = mujoco.MjData(model)
    check("mj_loadXML", True, "nbody=%d nq=%d nu=%d neq=%d ngeom=%d"
          % (model.nbody, model.nq, model.nu, model.neq, model.ngeom))
    check("mujoco version 2.3.7", mujoco.__version__.startswith("2.3.7"),
          mujoco.__version__)

    # ---- structure counts
    n_bodies_expected = 1 + 13   # world + Root + 4 segments/leg * 2 + 2 welded contact boxes/leg * 2
    check("body count", model.nbody == n_bodies_expected,
          "nbody=%d (expected %d incl. world)" % (model.nbody, n_bodies_expected))
    check("8 hinge joints", model.njnt == 8, "njnt=%d" % model.njnt)
    check("12 muscle actuators", model.nu == 12, "nu=%d" % model.nu)
    check("14 tendons (12 muscle + 2 toe springs)",
          model.ntendon == 14, "ntendon=%d" % model.ntendon)

    # ---- gate 2: 2 s passive drop
    print("== gate 2: 2 s passive drop (timestep %g s) ==" % model.opt.timestep)
    mujoco.mj_resetData(model, data)
    n_steps = int(round(2.0 / model.opt.timestep))
    min_foot_z = 1e9
    for i in range(n_steps):
        mujoco.mj_step(model, data)
        if not np.isfinite(data.qpos).all() or not np.isfinite(data.qvel).all():
            check("passive drop finite", False, "NaN at step %d (t=%.3f s)" % (i, data.time))
            break
        min_foot_z = min(min_foot_z, float(data.xpos[model.body("foot_L").id][2]),
                         float(data.xpos[model.body("foot_R").id][2]))
    else:
        finite = bool(np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all())
        maxq = float(np.abs(data.qpos).max()) if finite else float("nan")
        check("passive drop finite (2 s, %d steps)" % n_steps, finite,
              "max|qpos|=%.4g, min foot height=%.4f m" % (maxq, min_foot_z))
        check("drop settles on ground (0 <= foot z < 0.1)",
              0.0 <= min_foot_z < 0.1, "min foot z=%.4f m" % min_foot_z)

    # ---- gate 3: actuator resting lengths finite
    print("== gate 3: 12 actuators finite resting lengths ==")
    lr = model.actuator_lengthrange  # (nu, 2)
    ok = np.isfinite(lr).all() and (lr[:, 1] > lr[:, 0]).all()
    names = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(model.nu)]
    check("actuator_lengthrange finite + lo<hi for all 12", bool(ok),
          "; ".join("%s[%.3f,%.3f]" % (n, lr[i, 0], lr[i, 1]) for i, n in enumerate(names)))
    # tendon lengths at rest, finite
    mujoco.mj_forward(model, data)
    ten_len = data.ten_length
    check("tendon lengths finite at rest", bool(np.isfinite(ten_len).all()),
          "rest lengths: " + " ".join("%.3f" % v for v in ten_len))

    # ---- gate 4: joint ranges vs aproj
    print("== gate 4: joint ranges match aproj (deg->rad) ==")
    dump_joints = {}
    for b in dump["bodies"]:
        if "joint" in b:
            dump_joints[b["joint"]["name"]] = b["joint"]
    all_ok = True
    detail = []
    for j in range(model.njnt):
        jname = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        src = dump_joints.get(jname)
        if src is None:
            all_ok = False
            detail.append("%s: no source" % jname)
            continue
        limited = bool(model.jnt_limited[j])
        if src["enable_limits"] and not (src["lower_deg"] == 0.0 and src["upper_deg"] == 0.0):
            lo_rad = math.radians(src["lower_deg"])
            hi_rad = math.radians(src["upper_deg"])
            dev = max(abs(model.jnt_range[j][0] - lo_rad), abs(model.jnt_range[j][1] - hi_rad))
            ok_j = limited and dev <= 1e-6   # far tighter than 0.5-deg rounding
            detail.append("%s mj=[%.6f,%.6f] src=[%.6f,%.6f] dev=%.2e%s"
                          % (jname, model.jnt_range[j][0], model.jnt_range[j][1],
                             lo_rad, hi_rad, dev, "" if ok_j else " MISMATCH"))
        else:
            ok_j = not limited
            detail.append("%s unlimited (mj limited=%s)" % (jname, limited))
        all_ok = all_ok and ok_j
    check("joint ranges match", all_ok, " | ".join(detail))

    # ---- fidelity readback: masses
    print("== fidelity readback: body masses (aproj vs MJCF) ==")
    src_boxes = {b["name"]: b for b in dump["bodies"] if b["type"] == "Box"}
    worst = 0.0
    for bname in ("Root", "femur_L", "tibia_L", "foot_L", "toe_L",
                  "femur_R", "tibia_R", "foot_R", "toe_R",
                  "toe_L_contact", "foot_L_contact", "toe_R_contact", "foot_R_contact"):
        bid = model.body(bname).id
        mj_mass = float(model.body_mass[bid])
        src_mass = src_boxes[bname]["mass_kg"]
        dev = abs(mj_mass - src_mass)
        worst = max(worst, dev)
        print("  %-16s aproj=%8.4f kg  mjcf=%8.4f kg  dev=%.2e" % (bname, src_mass, mj_mass, dev))
    check("body masses match", worst < 1e-9, "worst dev=%.2e kg" % worst)

    total_src = sum(src_boxes[n]["mass_kg"] for n in src_boxes)
    total_mj = float(model.body_mass[1:].sum())
    print("  total            aproj=%8.4f kg  mjcf=%8.4f kg" % (total_src, total_mj))

    print()
    print("RESULT: %d passed, %d failed" % (len(PASS), len(FAIL)))
    for name, det in FAIL:
        print("  FAILED: %s  %s" % (name, det))
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
