"""Goal-3 timestep convergence ladder (2026-09-23).

Question (Hyfydy critique of MuJoCo fixed-step integration): how much does
the walker's trajectory depend on the 2 ms plant timestep?

Method: replay the SAME recorded 20-key-muscle activation history
(spinal_run_s3k.npz 'act', 2 ms sampling, zero-order hold) into the SAME
harnessed model, integrating the plant at 2 / 1 / 0.5 ms (1/2/4 substeps
per control sample). All variants share initial state and control signal;
trajectory divergence is pure integration error. Two configs: ground
(contacts active) and air (no ground, pin_rot rig) to separate
contact-solver error from smooth-dynamics error.

Read-only w.r.t. production: no model XML edits, no runner/params changes;
model.opt.timestep is a runtime override. New file only.
"""
import io
import sys
import time
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import numpy as np
import mujoco

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import runner  # safe: main-guard; we only use its build helpers

NPZ = HERE / "spinal_run_s3k.npz"
OUT_MD = (HERE / "reports_20260923" / "goal3_convergence_ladder.md")
DT_CTRL = 0.002  # control sampling = network step = production 2 ms


def build_model(air: bool):
    model = mujoco.MjModel.from_xml_path(str(runner.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    runner.apply_start_pose(model, data)
    key_pose = runner.capture_pose(model, data)
    if air:
        model = runner.apply_harness(model, data, kxy=1500.0,
                                     no_ground=True, pin_rot=True)
    else:
        model = runner.apply_harness(model, data, kxy=1500.0)
    data = mujoco.MjData(model)
    runner.seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)
    return model, data


def replay(model, data, act_ids, act_hist, n_sub):
    """Zero-order-hold control at DT_CTRL, plant stepped n_sub times."""
    model.opt.timestep = DT_CTRL / n_sub
    key_adr = {jn: model.jnt_qposadr[model.joint(jn).id]
               for jn in runner.KEY_JOINTS}
    nsteps = act_hist.shape[0]
    q = np.zeros((nsteps, len(runner.KEY_JOINTS)))
    qf = np.zeros((nsteps, model.nq))
    up = data.act.copy()  # stays 0 elsewhere
    for k in range(nsteps):
        up[act_ids] = act_hist[k]
        for _ in range(n_sub):
            data.act[:] = up
            mujoco.mj_step(model, data)
        for i, jn in enumerate(runner.KEY_JOINTS):
            q[k, i] = data.qpos[key_adr[jn]]
        qf[k] = data.qpos
    return q, qf


def compare(qa, qb, name_a, name_b):
    d = np.degrees(qa - qb)  # nsteps x 13 key joints
    rms_j = np.sqrt(np.mean(d ** 2, axis=0))
    per = {jn: (float(rms_j[i]), float(np.abs(d[:, i]).max()))
           for i, jn in enumerate(runner.KEY_JOINTS)}
    rms_all = float(np.sqrt(np.mean(d ** 2)))
    mx = float(np.abs(d).max())
    print(f"  {name_a} vs {name_b}: overall RMS {rms_all:.3f} deg, "
          f"max {mx:.3f} deg")
    worst = sorted(per.items(), key=lambda kv: -kv[1][0])[:5]
    for jn, (r, m) in worst:
        print(f"    {jn:18s} RMS {r:7.3f} deg  max {m:8.3f} deg")
    return rms_all, mx, per


def main():
    z = np.load(NPZ, allow_pickle=True)
    act_hist = z["act"]              # (8000, 20)
    key_acts = [str(s) for s in z["key_acts"]]
    print(f"replaying {act_hist.shape[0]} control samples "
          f"({act_hist.shape[0] * DT_CTRL:.1f} s), {len(key_acts)} muscles")

    lines = ["# Goal-3 timestep convergence ladder (2026-09-23)", "",
             "Same recorded 20-muscle activation control (spinal_run_s3k.npz,",
             "2 ms zero-order hold), plant integrated at 2 / 1 / 0.5 ms.",
             "Shared initial state; divergence = integration error only.", ""]

    for air in (False, True):
        tag = "air (no ground contact)" if air else "ground (contacts on)"
        print(f"== config: {tag} ==")
        model, data = build_model(air)
        act_ids = np.array([model.actuator(nm).id for nm in key_acts])
        res = {}
        wall = {}
        for n_sub, label in ((1, "2ms"), (2, "1ms"), (4, "0.5ms")):
            m, d = build_model(air)  # fresh identical start per variant
            t0 = time.perf_counter()
            q, qf = replay(m, d, act_ids, act_hist, n_sub)
            wall[n_sub] = time.perf_counter() - t0
            res[n_sub] = (q, qf)
            fell = qf[-1, 1] if m.nq > 1 else float("nan")  # pelvis height col
            print(f"  {label}: wall {wall[n_sub]:.1f}s, "
                  f"final pelvis_y {fell:.3f} m")
        r1, m1, _ = compare(res[1][0], res[4][0], "2ms", "0.5ms")
        r2, m2, _ = compare(res[2][0], res[4][0], "1ms", "0.5ms")
        lines += [f"## {tag}", "",
                  f"- 2 ms vs 0.5 ms: overall RMS **{r1:.3f} deg**, "
                  f"max {m1:.3f} deg",
                  f"- 1 ms vs 0.5 ms: overall RMS **{r2:.3f} deg**, "
                  f"max {m2:.3f} deg",
                  f"- wall time: 2 ms {wall[1]:.0f}s, 1 ms {wall[2]:.0f}s, "
                  f"0.5 ms {wall[4]:.0f}s", ""]

    OUT_MD.parent.mkdir(exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"saved {OUT_MD}")


if __name__ == "__main__":
    main()
