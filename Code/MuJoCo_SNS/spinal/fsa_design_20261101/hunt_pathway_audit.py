"""Read actual command-gated Hunt electrical pathways; never run curriculum.

Uses the saved condition's options/weights, constructs the neural circuit,
and applies existing weight policies at specified commands without body steps.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"
os.environ["AARL_NET"] = "w2lvar"
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import mujoco

HERE = Path(__file__).resolve().parent
S = HERE.parent
H = S / "hunt_training"
sys.path[:0] = [str(H), str(S), str(S / "functional_gait2392")]
import controller as C

MODEL = S.parents[2] / "Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml"


def options_from_physical_condition(cfg):
    stage = cfg["stage"]
    return dict(
        use_timing_afferents=False if stage == 1 else cfg.get("timing_afferents", True),
        axis_feedback=0. if stage == 1 else cfg.get("axis_feedback", 0.),
        balance=None if stage == 1 else cfg.get("balance"),
        weights_directory=cfg.get("weights_directory", str(H / "runs/fitted")),
        timing_scale=cfg.get("timing_scale", 1.),
        reciprocal_gain=cfg.get("reciprocal_gain", 0.),
        commissural_gain=0. if stage == 1 else cfg.get("commissural_gain", 0.),
        biarticular_pair_gain=cfg.get("biarticular_pair_gain", 0.),
        cross_contact_gain=cfg.get("cross_contact_gain", 0.),
        commissural_e_gain=0. if stage == 1 else cfg.get("commissural_e_gain", 0.),
        ipsilateral_heel_reset_gain=0. if stage == 1 else cfg.get("ipsilateral_heel_reset_gain", 0.),
        bilateral_midline_drive=cfg.get("bilateral_midline_drive", False),
        learn_reciprocal=cfg.get("learn_reciprocal", False),
        **({"hip_abduction_goal": cfg["hip_abduction_goal"]} if "hip_abduction_goal" in cfg else {}),
    )


def run(path):
    cfg = json.loads(path.read_text())
    model = mujoco.MjModel.from_xml_path(str(MODEL))
    names = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(model.nu)]
    net = C.build(names, dt=.002, options=options_from_physical_condition(cfg))
    # PT recruitment requires a separately calibrated mechanical operating
    # point. Its writes address only PT_BAL/PT_IN source columns, disjoint
    # from the native afferent classes audited here. Do not invent that body
    # calibration or misreport a static constructor as a full PT assessment.
    unbound_balance = net.__dict__.pop("balance", None)
    c = net.compiled
    categories = ("Ia_", "IaIN_", "II_", "IIX_", "IIIN_", "Ib_", "IBIN_", "IBEXC_")
    scopes = {key: [i for i in net.mn_edges if net.names[c.pre[i]].startswith(key)] for key in categories}
    cases = []
    for gate in (0., .5, 1.):
        net.apply_weights(gate)
        net.sync()
        cases.append({"command_gate": gate, "source_classes_to_MN": {
            key: {"allocated_edges": len(indices),
                  "active_edges": int(np.count_nonzero(c.gmax[indices])),
                  "maximum_g_uS": float(np.max(c.gmax[indices])) if indices else 0.}
            for key, indices in scopes.items()}})
    return {"machine": "easteregg2", "scope": "Static actual neural runtime weight-policy audit; no neural/body integration",
        "condition": str(path), "condition_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_controller": str(Path(C.__file__).resolve()),
        "source_controller_sha256": hashlib.sha256(Path(C.__file__).read_bytes()).hexdigest(),
        "options": options_from_physical_condition(cfg), "neurons": len(c.V),
        "input_count": len(net.inputs), "cases": cases,
        "PT_balance_calibration": "Excluded from this native-afferent mask audit; PT-specific gain writes were not evaluated" if unbound_balance is not None else "not present",
        "interpretation": "An allocated builder edge is not an active runtime path. HUNT_AF timing inputs and PF/central sensory projections need their separate pathway ledger."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.condition.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(args.output), "neurons": result["neurons"], "cases": result["cases"]}))
