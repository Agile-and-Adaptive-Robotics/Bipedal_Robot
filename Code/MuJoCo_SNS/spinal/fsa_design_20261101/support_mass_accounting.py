"""Explicit robot-body mass, excluding unrelated ground/environment bodies."""
import argparse
import hashlib
import json
from pathlib import Path


def selected_robot_mass(model, root_name):
    import mujoco
    root = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, root_name)
    if root <= 0:
        raise ValueError("An explicit non-world robot root body is required")
    selected = []
    for body in range(1, model.nbody):
        ancestor = body
        while ancestor not in (0, root):
            ancestor = int(model.body_parentid[ancestor])
        if ancestor == root:
            selected.append(body)
    return float(model.body_mass[selected].sum()), selected


def renormalize_saved_record(result_path, mechanical_path):
    result = json.loads(result_path.read_text())
    mechanical = json.loads(mechanical_path.read_text())
    mass = mechanical["body_mass_ledger"]
    if abs(result["body_weight_N"] / 9.81 - mechanical["compiled_model_total_mass_kg"]) > 1e-7:
        raise ValueError("The saved denominator does not match this mechanical model")
    weight = mass["moving_body_mass_kg"] * 9.81
    fields = ("mean_ground_load_bw", "mean_upward_rig_support_bw", "mean_downward_rig_force_bw", "mean_net_rig_force_bw")
    rows = {key: {"historical_fraction": result[key],
                  "mean_force_N_recovered_from_record": result[key] * result["body_weight_N"],
                  "selected_robot_weight_fraction": result[key] * result["body_weight_N"] / weight}
            for key in fields}
    return {"machine": "easteregg2", "operation": "saved_data_denominator_correction_only",
        "source_result": str(result_path), "mechanical_mass_ledger": str(mechanical_path),
        "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (result_path, mechanical_path)},
        "historical_denominator_N": result["body_weight_N"], "corrected_robot_weight_N": weight,
        "robot_mass_kg": mass["moving_body_mass_kg"], "excluded_fixed_ground_mass_kg": mass["fixed_to_world_mass_kg"],
        "metrics": rows, "raw_files_unchanged": True,
        "scope": "Same saved forces and trajectories; exclude fixed ground from this Gait robot's BW denominator. No new simulation or acceptance change."}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--result", type=Path, required=True)
    p.add_argument("--mechanical", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(f"Preserve the existing analysis record; choose a new output: {a.output}")
    z = renormalize_saved_record(a.result.resolve(), a.mechanical.resolve())
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(z, indent=2), encoding="utf-8")
    print(json.dumps({key: row["selected_robot_weight_fraction"] * 100 for key, row in z["metrics"].items()}))
