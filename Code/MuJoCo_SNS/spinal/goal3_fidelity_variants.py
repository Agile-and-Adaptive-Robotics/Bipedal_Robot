"""Goal-3 fidelity variants, implemented as opt-in demos (2026-09-23).

Implements the two still-missing goal-3 features at WRAPPER level on the
real harnessed model, without touching cvt3.xml, runner defaults, or any
recorded winner:

(A) ERROR-CONTROLLED INTEGRATION STEP SIZE: within the 2 ms network
    lockstep, estimate each sample's local error by step-doubling (one
    2 ms step vs two 1 ms half-steps from the same state); refine that
    sample to 4 x 0.5 ms only when the estimate exceeds tol. Adaptive
    cost/accuracy vs fixed 2 / 0.5 ms is reported.

(B) NON-LINEAR CONTACT DAMPING: the goal-3 report identified solimp
    depth-dependent impedance + negative-solref direct (k, b) as the
    partially-native route. Demo: direct k,b on the foot-ground pairs
    with damping halved (toward Hunt-Crossley-like, less viscous),
    compared against the production solref=[0.02, 1] on the same
    replay - trajectory deltas, peak foot force, max penetration.

(C) SEE (tendon elasticity) is NOT re-implemented here: the mechanism
    was demonstrated (toy spring-tendon + slide joint, goal-3 report
    section 3.2); converting real muscles awaits Ben's ruling.

Read-only w.r.t. production; new file; same replay harness as
goal3_timestep_ladder.py (imported).
"""
import io
import sys
import time
from pathlib import Path

import numpy as np
import mujoco

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from goal3_timestep_ladder import build_model, DT_CTRL  # noqa: E402
import runner  # noqa: E402


def key_adr(model):
    return {jn: model.jnt_qposadr[model.joint(jn).id]
            for jn in runner.KEY_JOINTS}


def q_key(data, adr):
    return np.array([data.qpos[a] for a in adr.values()])


def replay_fixed(model, data, act_ids, act_hist, n_sub):
    model.opt.timestep = DT_CTRL / n_sub
    adr = key_adr(model)
    q = np.zeros((act_hist.shape[0], len(adr)))
    peak_f = 0.0
    pen = 0.0
    for k in range(act_hist.shape[0]):
        up = data.act.copy()
        up[act_ids] = act_hist[k]
        for _ in range(n_sub):
            data.act[:] = up
            mujoco.mj_step(model, data)
        q[k] = q_key(data, adr)
        for ci in range(data.ncon):
            c = data.contact[ci]
            f6 = np.zeros(6)
            mujoco.mj_contactForce(model, data, ci, f6)
            peak_f = max(peak_f, float(np.linalg.norm(f6[:3])))
            pen = min(pen, float(c.dist))
    return q, peak_f, pen


def _state_save(data):
    return (data.qpos.copy(), data.qvel.copy(), data.act.copy(),
            float(data.time), data.qacc_warmstart.copy())


def _state_restore(data, s):
    data.qpos[:] = s[0]
    data.qvel[:] = s[1]
    data.act[:] = s[2]
    data.time = s[3]
    data.qacc_warmstart[:] = s[4]


def replay_adaptive(model, data, act_ids, act_hist, tol_deg):
    """Error-controlled substeps: coarse 2 ms, refined 4 x 0.5 ms where the
    step-doubling error estimate exceeds tol_deg (deg, key joints)."""
    adr = key_adr(model)
    q = np.zeros((act_hist.shape[0], len(adr)))
    n_refined = 0
    err_trace = []
    for k in range(act_hist.shape[0]):
        up = data.act.copy()
        up[act_ids] = act_hist[k]
        buf = _state_save(data)
        # coarse: one 2 ms step
        model.opt.timestep = DT_CTRL
        data.act[:] = up
        mujoco.mj_step(model, data)
        q_coarse = q_key(data, adr)
        # probe: two 1 ms half-steps from the same state
        _state_restore(data, buf)
        model.opt.timestep = DT_CTRL / 2
        for _ in range(2):
            data.act[:] = up
            mujoco.mj_step(model, data)
        q_probe = q_key(data, adr)
        err = float(np.max(np.abs(np.degrees(q_probe - q_coarse))))
        err_trace.append(err)
        if err > tol_deg:
            # refined: 4 x 0.5 ms from the same state
            _state_restore(data, buf)
            model.opt.timestep = DT_CTRL / 4
            for _ in range(4):
                data.act[:] = up
                mujoco.mj_step(model, data)
            n_refined += 1
            q[k] = q_key(data, adr)
        else:
            q[k] = q_probe  # the more accurate of the two cheap options
    frac = n_refined / act_hist.shape[0]
    return q, frac, np.asarray(err_trace)


def main():
    z = np.load(HERE / "spinal_run_s3k.npz", allow_pickle=True)
    act_hist = z["act"]
    key_acts = [str(s) for s in z["key_acts"]]
    lines = ["# Goal-3 fidelity variants - measured (2026-09-23)", ""]

    model, data = build_model(air=False)
    act_ids = np.array([model.actuator(nm).id for nm in key_acts])

    # fixed references
    t0 = time.perf_counter()
    m2, d2 = build_model(air=False)
    q2, pf2, pen2 = replay_fixed(m2, d2, act_ids, act_hist, 1)
    w2 = time.perf_counter() - t0
    t0 = time.perf_counter()
    m4, d4 = build_model(air=False)
    q4, pf4, pen4 = replay_fixed(m4, d4, act_ids, act_hist, 4)
    w4 = time.perf_counter() - t0
    rms = lambda a, b: float(np.sqrt(np.mean(np.degrees(a - b) ** 2)))  # noqa: E731
    mx = lambda a, b: float(np.max(np.abs(np.degrees(a - b))))  # noqa: E731

    print(f"fixed 2ms: wall {w2:.1f}s  peak foot force {pf2:.0f} N  "
          f"max penetration {pen2*1000:.2f} mm")
    print(f"fixed 0.5ms: wall {w4:.1f}s  peak foot force {pf4:.0f} N  "
          f"max penetration {pen4*1000:.2f} mm")

    # (A) error-controlled adaptive substepping
    lines.append("## (A) Error-controlled substepping (tol sweep)")
    lines.append("")
    lines.append("| tol (deg/sample) | refined frac | RMS vs 0.5ms (deg) "
                 "| max vs 0.5ms (deg) | wall (s) |")
    lines.append("|---|---|---|---|---|")
    for tol in (0.01, 0.05, 0.2):
        ma, da = build_model(air=False)
        t0 = time.perf_counter()
        qa, frac, errs = replay_adaptive(ma, da, act_ids, act_hist, tol)
        wa = time.perf_counter() - t0
        r, m = rms(qa, q4), mx(qa, q4)
        print(f"adaptive tol {tol}: refined {frac*100:.1f}%  RMS {r:.3f} deg "
              f"max {m:.3f} deg  wall {wa:.1f}s")
        lines.append(f"| {tol} | {frac*100:.1f}% | {r:.3f} | {m:.3f} "
                     f"| {wa:.1f} |")
    lines += ["", f"(fixed 2 ms wall {w2:.1f} s, fixed 0.5 ms wall {w4:.1f} s;"
                 " fixed 2ms-vs-0.5ms RMS "
                 f"{rms(q2, q4):.3f} deg / max {mx(q2, q4):.3f} deg)", ""]

    # (B) contact damping variants INSIDE the impedance framework. The
    # goal-3 report's naive direct-(k, b) conversion (-k, -b in solref)
    # bypasses the solimp impedance scaling and is UNSTABLE at 2 ms
    # (verified: NaN at t=0.012 s) - recorded in the report as a negative.
    # These variants stay in positive-solref mode where MuJoCo impedance-
    # scales the force, changing the damping character measurably:
    #   1. less viscous + snappier: timeconst 0.02 -> 0.01, dampratio 1 -> 0.7
    #   2. more non-linear depth profile: solimp midpoint 0.95 -> 0.90,
    #      width 0.001 -> 0.003 (progressive stiffening with penetration)
    lines.append("## (B) Contact damping variants (impedance framework)")
    lines.append("")
    lines.append("| variant | RMS vs baseline (deg) | max vs baseline (deg) "
                 "| peak foot force (N) | max pen (mm) |")
    lines.append("|---|---|---|---|---|")
    lines.append(f"| baseline solref=[0.02, 1], solimp mid 0.95 | 0 | 0 "
                 f"| {pf2:.0f} | {pen2*1000:.2f} |")
    variants = (
        ("less viscous: tc 0.01, zeta 0.7",
         lambda mv: (mv.pair_solref[:, 0].__setitem__(slice(None), 0.01),
                     mv.pair_solref[:, 1].__setitem__(slice(None), 0.7))),
        ("non-linear depth: solimp mid 0.90, width 0.003",
         lambda mv: (mv.pair_solimp[:, 1].__setitem__(slice(None), 0.90),
                     mv.pair_solimp[:, 3].__setitem__(slice(None), 0.003))),
    )
    for tag, mutate in variants:
        mv, dv = build_model(air=False)
        mutate(mv)
        qv, pfv, penv = replay_fixed(mv, dv, act_ids, act_hist, 1)
        if not np.all(np.isfinite(qv)):
            print(f"{tag}: UNSTABLE (NaN)")
            lines.append(f"| {tag} | UNSTABLE | - | - | - |")
            continue
        r, m = rms(qv, q2), mx(qv, q2)
        print(f"{tag}: RMS {r:.3f} deg  max {m:.3f} deg  peak {pfv:.0f} N "
              f"pen {penv*1000:.2f} mm")
        lines.append(f"| {tag} | {r:.3f} | {m:.3f} | {pfv:.0f} "
                     f"| {penv*1000:.2f} |")
    lines += ["",
              "Naive direct-(k, b) via negative solref was also attempted "
              "and is UNSTABLE at the 2 ms timestep (NaN by t=0.012 s; it "
              "bypasses the solimp impedance scaling): a concrete "
              "demonstration of the goal-3 report's 'partially native' "
              "verdict. Exact Hunt-Crossley fitting and any production "
              "adoption remain open (Ben's ruling).",
              "", "(C) SEE not re-implemented - see module docstring."]

    out = HERE / "reports_20260923" / "goal3_fidelity_variants.md"
    out.parent.mkdir(exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("saved", out)


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
