"""Opt-in Hunt local-reflex design, authoritative at every command update.

Historical Hunt source and saved weights are unchanged. The candidate restores
the existing homonymous II template to the local reflex parameter set, derives
four declared signed-pool functions, and updates the actual runtime gain arrays.
It does not install the ambiguous antagonistic-II policy or change HUNT_AF.
"""
import copy
import numpy as np
if __package__:
    from . import graded_profile as GP
else:
    import graded_profile as GP


def install(net, profile=None):
    if profile is None:
        return {"applied": False, "default_noop": True}
    if hasattr(net, "hunt_graded_fsa"):
        raise ValueError("Hunt profile already installed; use a fresh candidate")
    if hasattr(net, "graded_fsa_ledger"):
        raise ValueError("Use a fresh Hunt candidate, after all historical weight transforms and before another FSA profile")
    c = net.compiled
    required = ("gmax", "post", "pre", "c")
    if not all(hasattr(c, name) for name in required):
        raise ValueError("Actual Hunt sparse wrapper required")
    if not all(hasattr(net, name) for name in ("original_g", "reflex_gain", "reflex_edges", "sync")):
        raise ValueError("Install after actual Hunt weights/command policy construction")
    # Stage the transformation: the original compiled arrays and authoritative
    # gain buffers can be restored if any functional/bilateral check rejects it.
    backup = {key: getattr(net, key).copy() for key in ("original_g", "reflex_gain", "reflex_edges")}
    prior_g = c.gmax.copy()
    prior_learn = net.options.get("learn_reciprocal", False)
    prior_cap = c.c.c_m.copy()
    prior_time = c.c.time_factor_membrane.copy()
    selected = []
    prior_II = []
    for i in net.mn_edges:
        pre, post = net.names[c.pre[i]], net.names[c.post[i]]
        kind = GP.edge_kind(pre, post, c.reversal[i])
        if kind is not None and net.original_g[i] > 0:
            selected.append(int(i))
            if kind == "II_homonymous":
                prior_II.append({"source": pre, "target": post,
                                 "old_active_g_uS": float(c.gmax[i]),
                                 "template_g_uS": float(net.original_g[i])})
    selected = np.asarray(selected, dtype=int)
    if not len(selected) or len(prior_II) != 92:
        raise ValueError("Expected all92 allocated homonymous-II template paths")
    try:
        # The new condition starts local reflex functions from the documented
        # canonical template. Historical fitted zeros/gains remain in originals.
        c.gmax[selected] = net.original_g[selected]
        net.sync()
        design = GP.apply_profile(net, profile, family="w2lvar")
        values = c.gmax[selected].copy()
        net.reflex_gain[selected] = values
        net.reflex_edges = np.union1d(net.reflex_edges, selected).astype(int)
        net.options["learn_reciprocal"] = True
        # apply_weights already writes reflex_gain after clearing MN edges.
        # Updating this parameter set avoids a hidden per-step clamp/wrapper.
        net.hunt_graded_fsa = {
            "schema": "aarl_hunt_graded_fsa_candidate_v1",
            "profile": copy.deepcopy(profile),
            "protected_indices": selected,
            "initial_derived_values": values,
            "restored_homonymous_II": prior_II,
            "design_ledger": design,
            "gain_parameters": "net.reflex_gain; includes homonymousII in net.reflex_edges",
            "scope": "New local-reflex initialization/function; retains NaP/PF/HUNT_AF/PT/geometry policies",
            "unqualified": "No full-family design, Hunt training, balance, or curriculum pass inherited",
        }
        verify_runtime(net)
        return net.hunt_graded_fsa
    except Exception:
        for key, value in backup.items():
            setattr(net, key, value)
        c.gmax[:] = prior_g
        c.c.c_m[:] = prior_cap
        c.c.time_factor_membrane[:] = prior_time
        net.options["learn_reciprocal"] = prior_learn
        for field in ("hunt_graded_fsa", "graded_fsa_profile", "graded_fsa_ledger"):
            net.__dict__.pop(field, None)
        net.sync()
        raise


def verify_runtime(net):
    """Check dynamic authoritative buffers, not a stale one-time matrix hash.

    Gait command and geometry may change other registered weights. Local
    reflex parameters must be explicit and bilaterally tied if later fitted.
    Future fitting must regenerate/verify its declared functional gain targets.
    """
    ledger = net.hunt_graded_fsa
    ids = ledger["protected_indices"]
    c = net.compiled
    if not np.array_equal(c.gmax[ids], net.reflex_gain[ids]):
        raise ValueError("Command policy overwrote the new local-reflex parameters")
    if not np.array_equal(c.c.g_max_non[c.post, c.pre], c.gmax):
        raise ValueError("Sparse/dense Hunt weights are not synchronized")
    if not np.all(np.isfinite(c.gmax[ids])) or np.any(c.gmax[ids] < 0):
        raise ValueError("Invalid local-reflex conductance")
    pairs = {}
    for i in ids:
        key = GP.canonical_edge(net.names[c.pre[i]], net.names[c.post[i]])
        a, b, relation, side = key
        pairs.setdefault((a, b, relation), {})[side] = float(c.gmax[i])
    for key, pair in pairs.items():
        if set(pair) != {"r", "l"} or not np.isclose(pair["r"], pair["l"], rtol=1e-12, atol=1e-14):
            raise ValueError("Unpaired/asymmetrical local-reflex parameter: " + str(key))
    expected_cap = {row["index"]: row["capacitance_backend_microF"]
                    for row in ledger["design_ledger"]["encoders"]}
    if any(c.c.c_m[index] != value for index, value in expected_cap.items()):
        raise ValueError("Encoder capacitance was overwritten")
    return {"local_reflex_pairs": len(pairs), "restored_II_paths": len(ledger["restored_homonymous_II"])}


def serializable_ledger(net):
    value = copy.deepcopy(net.hunt_graded_fsa)
    value["protected_indices"] = value["protected_indices"].tolist()
    value["initial_derived_values"] = value["initial_derived_values"].tolist()
    return value
