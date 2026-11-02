"""Electrical-only Hunt profile consumption, command-policy, and II probes."""
import argparse
import json
from pathlib import Path

import numpy as np
if __package__:
    from .hunt_pathway_audit import C, MODEL, options_from_physical_condition, mujoco
    from .hunt_graded_profile import install, serializable_ledger, verify_runtime
    from .graded_profile import default_profile, field_digest
else:
    from hunt_pathway_audit import C, MODEL, options_from_physical_condition, mujoco
    from hunt_graded_profile import install, serializable_ledger, verify_runtime
    from graded_profile import default_profile, field_digest


def run(condition, output):
    cfg = json.loads(condition.read_text())
    m = mujoco.MjModel.from_xml_path(str(MODEL))
    acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]
    n = C.build(acts, dt=.002, options=options_from_physical_condition(cfg))
    # Same explicit local-path scope as hunt_pathway_audit: PT mechanical
    # calibration is separate, and PT-only columns are not probed here.
    n.__dict__.pop("balance", None)
    n.apply_weights(1.)
    before = field_digest(n.compiled.c)
    noop = install(n, None)
    assert noop["default_noop"] and before == field_digest(n.compiled.c)
    saved_gain = n.compiled.gmax.copy()
    design = install(n, default_profile("w2lvar"))
    policy = []
    for gate in (0., .5, 1., 0.):
        n.apply_weights(gate)
        n.sync()
        policy.append({"command": gate, **verify_runtime(n)})
    applied_ledger = serializable_ledger(n)
    c = n.compiled
    idx = {name: i for i, name in enumerate(n.names)}
    edges = {(int(pre), int(post)): i for i, (pre, post) in enumerate(zip(c.pre, c.post))}
    kept = []
    for side in ("r", "l"):
        for a, b in (("II_soleus_", "IIX_soleus_"), ("IIX_soleus_", "MN_soleus_")):
            kept.append(edges[(idx[a + side], idx[b + side])])
    g = c.gmax[kept].copy()
    c.gmax[:] = 0.
    c.gmax[kept] = g
    n.sync()
    # This is an actual sparse compiled-network forward pass after explicitly
    # isolating the two existing II pathways, not a replacement toy circuit.
    channels = [idx[p + side] for side in ("r", "l") for p in ("II_soleus_", "IIX_soleus_", "MN_soleus_")]
    trace = np.zeros((1500, 6))
    u = np.zeros(len(n.inputs))
    for side in ("r", "l"):
        u[n.ports["POST_soleus_" + side]] = 2.5
    for step in range(len(trace)):
        t = step * .002
        for side in ("r", "l"):
            u[n.ports["II_soleus_" + side]] = 1. if 1. <= t < 2. else 0.
        c.forward(u)
        trace[step] = c.V[channels]
    baseline = trace[449, [2, 5]]
    plateau = trace[949, [2, 5]]
    delta = plateau - baseline
    if not np.allclose(baseline, 2.5, atol=1e-8) or not np.allclose(delta, .5, atol=1e-6):
        raise AssertionError(f"Actual homonymousII endpoint mismatch: baseline={baseline},delta={delta}")
    if not np.array_equal(trace[:, :3], trace[:, 3:]):
        raise AssertionError("Isolated bilateral II trajectories differ")
    output.mkdir(parents=True, exist_ok=True)
    result = {"machine": "easteregg2", "condition": str(condition), "scope": "electrical_only",
        "neurons": len(c.V), "default_profile_None_exact": True, "command_policy": policy,
        "old_active_II_count": sum(row["old_active_g_uS"] != 0 for row in design["restored_homonymous_II"]),
        "new_II_count": len(design["restored_homonymous_II"]),
        "actual_II_MN_baseline_mV": baseline.tolist(), "actual_II_MN_delta_mV": delta.tolist(),
        "bilateral_II_trace_exact": True, "training_allowed": False, "curriculum_allowed": False,
        "scope_limit": "Local native-reflex and command-policy proof. HUNT_AF, NaP/PF, PT mechanical calibration, other interactions, and full-family design remain separate."}
    (output / "RESULTS.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (output / "APPLIED_PROFILE.json").write_text(json.dumps(applied_ledger, indent=2), encoding="utf-8")
    np.savez_compressed(output / "II_TRACE.npz", t=np.arange(len(trace)) * .002, voltage_mV=trace)
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.condition.resolve(), args.output.resolve())
