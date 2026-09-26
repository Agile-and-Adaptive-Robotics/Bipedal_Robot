"""Gates for goal4 (syn6 variant): (a) defaults, (b) variant counts,
(c1) 12 s network-only rhythm smoke, (c2) full runner air run.

Run (cwd = Code/MuJoCo_SNS/spinal, myo env):
  python reports_20260925/tmp/gate_syn6.py a
  python reports_20260925/tmp/gate_syn6.py b
  python reports_20260925/tmp/gate_syn6.py c1
  python reports_20260925/tmp/gate_syn6.py c2
"""
import io
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]   # .../spinal
sys.path.insert(0, str(HERE))
LOGS = HERE / "reports_20260925" / "logs"
LOGS.mkdir(parents=True, exist_ok=True)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
PY = sys.executable


def actuator_names():
    import mujoco
    xml = (HERE.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
           / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml")
    m = mujoco.MjModel.from_xml_path(str(xml))
    return [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(m.nu)]


def counts(net):
    # g_max_non is [n_neurons x n_neurons]; each synapse occupies one
    # entry (verified vs len(net.net.connections) on the default build:
    # both = 1186)
    g = net.compiled.g_max_non
    return (len(net.idx), len(net.inputs),
            int(np.count_nonzero(g)))


def gate_a():
    import numpy as np
    os.environ.pop("AARL_NET", None)
    import importlib
    import build_network as bn
    import params
    importlib.reload(params)
    importlib.reload(bn)
    assert params.G.get("syn6", 0.0) == 0.0, "syn6 selector must default 0"
    net = bn.build(actuator_names(), dt=params.DT, interleg=True)
    n, i, s = counts(net)
    print(f"GATE-A default build: neurons={n} inputs={i} synapses={s}")
    ok = (n, i, s) == (410, 376, 1186)
    print("GATE-A", "PASS" if ok else "FAIL", "(expected 410/376/1186)")
    return ok


def gate_b():
    import numpy as np
    os.environ["AARL_NET"] = "syn6"
    import build_network as bn
    import build_network_syn6 as s6
    import params
    import runner
    from muscle_map import classify
    # prune-set equality with the runner (module constant, not imported)
    assert s6.PRUNE_MUSCLES == runner.PRUNE_MUSCLES, \
        "syn6 prune set != runner.PRUNE_MUSCLES"
    print("prune set == runner.PRUNE_MUSCLES:", sorted(s6.PRUNE_MUSCLES))
    acts = actuator_names()
    net = bn.build(acts, dt=params.DT, interleg=True)
    n, i, s = counts(net)
    print(f"GATE-B syn6 build: neurons={n} inputs={i} synapses={s}")
    # coverage: every non-pruned actuator has a W row; prunes do not
    names_r = set(s6._load_synergy_tables()["r"]["names"])
    names_l = set(s6._load_synergy_tables()["l"]["names"])
    covered = names_r | names_l
    nonpruned = [a for a in acts if a not in runner.PRUNE_MUSCLES]
    print(f"W coverage: {len(covered)} names cover "
          f"{sum(1 for a in nonpruned if a in covered)}/{len(nonpruned)} "
          f"non-pruned actuators; prunes covered: "
          f"{sorted(runner.PRUNE_MUSCLES & covered)}")
    assert set(nonpruned) == covered, "W basis != non-pruned actuator set"
    # every runner-written port must exist
    need = ["DRIVE", "POSTURE", "BAL_PF", "BAL_DF", "BAL_TRK_EXT",
            "BAL_TRK_FLX", "BAL_LAT_R", "BAL_LAT_L"]
    for sd in net.sides:
        need += [f"HEEL_c_{sd}", f"TOE_c_{sd}", f"LOAD_c_{sd}"]
    for a in acts:
        need += [f"POST_{a}", f"Ia_{a}", f"II_{a}", f"Ib_{a}"]
    missing = [p for p in need if p not in net.inputs]
    print(f"port coverage: {len(need) - len(missing)}/{len(need)} present,"
          f" missing={missing[:5]}")
    # muscle/afferent maps complete (runner indexes every actuator)
    assert set(net.mn_names) == set(acts)
    # family/tau/df-channel evidence
    for sd in net.sides:
        tab = s6._load_synergy_tables()[sd]
        print(f"side {sd}: family={tab['family']}")
        print(f"         stance_frac={np.round(tab['stance_frac'], 3).tolist()}")
        print(f"         peak_phase={np.round(tab['peak_phase'], 1).tolist()}")
        print(f"         tau_mult={net.syn_tau_mult[sd]}")
        print(f"         df_channel=S{None if tab['df_channel'] is None else tab['df_channel'] + 1}")
        gs = net.syn_g_stats.get(sd, [])
        if gs:
            import numpy as np
            print(f"         PF->MN g: n={len(gs)} min={min(gs):.4f} "
                  f"max={max(gs):.4f} uS (all Eq-18 valid: "
                  f"{all(g < np.inf for g in gs)})")
    # spot-check one Eq-18 value against the closed form
    import numpy as np
    k = 0.5
    g_expect = k * 5.0 * 1.0 / (8.0 - k * 5.0)
    g_got, valid = s6._eq18_conductances(np.array([k]))
    print(f"Eq-18 spot check k=0.5: builder g={float(g_got[0]):.6f}, "
          f"closed form g={g_expect:.6f}, valid={bool(valid[0])}")
    ok = not missing and (n, i, s) == counts(net)
    print("GATE-B", "PASS" if ok and not missing else "PASS" if ok else "FAIL")
    return ok


def gate_c1():
    """12 s network-only rhythm at constant DRIVE=2.5 (air; no MuJoCo)."""
    import numpy as np
    from scipy.signal import find_peaks
    os.environ["AARL_NET"] = "syn6"
    import build_network as bn
    import params
    net = bn.build(actuator_names(), dt=params.DT, interleg=True)
    watch = (["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l"]
             + [f"PF_S{k}_r" for k in range(1, 7)]
             + [f"PF_S{k}_l" for k in range(1, 7)])
    ii = {nm: net.idx[nm] for nm in watch}
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = 2.5
    T = 12.0
    n = int(round(T / params.DT))
    log = np.zeros((n, len(watch)))
    for k in range(n):
        V = net.step(u)
        log[k] = [V[j] for j in ii.values()]
    tt = np.arange(n) * params.DT
    print(f"GATE-C1 network-only air rhythm, DRIVE=2.5 nA, {T:.0f} s "
          f"({n} steps)")
    ok = True
    for sd in ("r", "l"):
        e = log[:, watch.index(f"RG_E_{sd}")]
        f = log[:, watch.index(f"RG_F_{sd}")]
        msk = tt >= 2.0                       # skip build transient
        pk_e, _ = find_peaks(e[msk], prominence=0.3)
        pk_f, _ = find_peaks(f[msk], prominence=0.3)
        corr = float(np.corrcoef(e[msk], f[msk])[0, 1])
        swing = float(e[msk].max() - e[msk].min())
        print(f"  side {sd}: RG-E peaks(2-12s)={len(pk_e)} "
              f"amplitude={swing:.2f} mV  E-F corr={corr:+.3f}")
        per = (tt[msk][pk_e[-1]] - tt[msk][pk_e[0]]) / max(len(pk_e) - 1, 1) \
            if len(pk_e) > 1 else float("nan")
        print(f"    mean E period={per:.3f} s "
              f"(peak times {np.round(tt[msk][pk_e], 2).tolist()})")
        ok &= len(pk_e) >= 5 and corr < -0.5
    # PF layer profiles: all active + pairwise-correlation distinctness
    for sd in ("r", "l"):
        pfm = np.array([log[:, watch.index(f"PF_S{k}_{sd}")]
                        for k in range(1, 7)])
        msk = tt >= 2.0
        amps = pfm[:, msk].max(axis=1)
        c = np.corrcoef(pfm[:, msk])
        off = c[np.triu_indices(6, 1)]
        print(f"  side {sd}: PF max mV per channel="
              f"{np.round(amps, 2).tolist()}")
        print(f"    pairwise corr: min={off.min():+.3f} "
              f"max={off.max():+.3f} (n={off.size} pairs)")
        ok &= bool((amps > 0.3).all())
    print("GATE-C1", "PASS" if ok else "FAIL")
    np.save(LOGS / "syn6_net_smoke.npy", log)
    return ok


def gate_c2():
    """Full runner air run (--no-ground) under AARL_NET=syn6."""
    env = dict(os.environ)
    env["AARL_NET"] = "syn6"
    env["AARL_NPZ"] = "reports_20260925/logs/syn6_air_smoke.npz"
    env["PYTHONIOENCODING"] = "utf-8"
    cmd = [PY, "runner.py", "--no-ground", "--no-view"]
    print("GATE-C2 cmd:", " ".join(cmd))
    with open(LOGS / "syn6_air_smoke.log", "w", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), env=env, stdout=fh,
                             stderr=subprocess.STDOUT)
    print(f"runner exit code: {rc}")
    npz = LOGS / "syn6_air_smoke.npz"
    if not npz.exists():
        print("GATE-C2 FAIL: no npz produced; log tail:")
        print("".join(open(LOGS / "syn6_air_smoke.log",
                           encoding="utf-8").read()[-2000:]))
        return False
    import numpy as np
    from scipy.signal import find_peaks
    d = np.load(npz, allow_pickle=True)
    print("npz keys:", sorted(d.files))
    neuro = np.asarray(d["neuro"], dtype=float)
    names = [str(x) for x in d["neuro_names"]] \
        if "neuro_names" in d.files else None
    print("neuro shape:", neuro.shape, "names:", names)
    t = np.asarray(d["t"], dtype=float)
    ok = True
    if names:
        for nm in ("RG_E_r", "RG_F_r"):
            if nm in names:
                v = neuro[:, names.index(nm)]
                msk = t >= 5.0
                pk, _ = find_peaks(v[msk], prominence=0.3)
                print(f"  {nm}: peaks(t>=5s)={len(pk)} "
                      f"swing={v[msk].max() - v[msk].min():.2f} mV "
                      f"first/last={np.round(t[msk][pk[[0, -1]]], 2).tolist() if len(pk) else '[]'}")
                ok &= len(pk) >= 5
        for nm in ("PF_S1_r", "PF_S2_r", "PF_S3_r", "PF_S4_r"):
            if nm in names:
                v = neuro[:, names.index(nm)]
                print(f"  {nm}: max={v.max():.2f} mV "
                      f"mean(t>=5s)={v[t >= 5.0].mean():.3f}")
    q = np.asarray(d["q"], dtype=float) if "q" in d.files else None
    if q is not None:
        knee = q[:, 4]   # KEY_JOINTS[4] = knee_angle_r (deg, per AGENTS)
        print(f"  knee_angle_r: min={knee.min():.1f} max={knee.max():.1f} "
              f"deg over run")
    print("GATE-C2", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "a"
    res = {"a": gate_a, "b": gate_b, "c1": gate_c1, "c2": gate_c2}[which]()
    sys.exit(0 if res else 1)
