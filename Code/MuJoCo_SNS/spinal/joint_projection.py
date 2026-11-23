"""Project prescribed MuJoCo scalar joint equalities onto a pose."""
from weakref import WeakKeyDictionary

import mujoco
import numpy as np


_plans = WeakKeyDictionary()


def _dependency_layers(model, active):
    rows = np.flatnonzero((model.eq_type == mujoco.mjtEq.mjEQ_JOINT) & active)
    dependent = model.eq_obj1id[rows]
    if len(np.unique(dependent)) != len(rows):
        raise ValueError("A prescribed pose requires one active equality per dependent joint")
    layers = []
    while rows.size:
        parents = model.eq_obj2id[rows]
        ready = ~np.isin(parents, model.eq_obj1id[rows])
        if not np.any(ready):
            names = [model.joint(int(j)).name for j in model.eq_obj1id[rows]]
            raise ValueError(f"Cyclic scalar joint equalities cannot be projected directly: {names}")
        group = rows[ready]
        parent = parents[ready]
        has_parent = parent >= 0
        layers.append((group, model.jnt_qposadr[model.eq_obj1id[group]], has_parent,
                       model.jnt_qposadr[parent[has_parent]]))
        rows = rows[~ready]
    return layers


def project_joint_followers(model, data):
    """Set dependent qpos from all five polynomial coefficients in dependency order.

    A joint without a second joint is fixed at its reference plus the constant
    coefficient. Individual inactive equalities are skipped. Explicit geometry
    projection also works when the dynamics solver's global equality flag is off.
    """
    active = data.eq_active if hasattr(data, "eq_active") else model.eq_active
    topology = (active.tobytes(), model.eq_type.tobytes(), model.eq_obj1id.tobytes(),
                model.eq_obj2id.tobytes(), model.jnt_qposadr.tobytes())
    cached = _plans.get(model)
    if cached is None or cached[0] != topology:
        cached = (topology, _dependency_layers(model, active.astype(bool)))
        _plans[model] = cached
    for rows, dependent, has_parent, parent in cached[1]:
        x = np.zeros(len(rows))
        x[has_parent] = data.qpos[parent] - model.qpos0[parent]
        c0, c1, c2, c3, c4 = model.eq_data[rows, :5].T
        displacement = c0 + x * (c1 + x * (c2 + x * (c3 + x * c4)))
        data.qpos[dependent] = model.qpos0[dependent] + displacement
