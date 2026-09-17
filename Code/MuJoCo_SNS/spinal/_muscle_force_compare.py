"""Muscle-force comparison: MuJoCo (rigid-tendon muscle actuator) vs
OpenSim Thelen2003 (compliant tendon) along the subject01 walk trajectory.

OpenSim side: zz_bsolve_MuscleAnalysis_TendonForce.sto (real AnalyzeTool
output along the IK trajectory). MuJoCo side: replay bsolve's qpos
frames in the patched model with act=1 and read actuator_force (rigid
tendon => force = Fmax*F-L*FV at L_mt directly).

Per muscle: R^2 (raw), shape R^2 (mean+scale removed), best lag via
cross-correlation (frames), mean force ratio. Outputs CSV + figure.
"""
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import mujoco

import bsolve_ik as B
import kine_ref
import runner as R

HERE = Path(__file__).parent
OSIM_DIR = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
                r"\Gait2392_Robotbody\ResultsBSolve")
MUSCLES = ["soleus_r", "med_gas_r", "lat_gas_r", "tib_ant_r", "vas_lat_r",
           "rect_fem_r", "semimem_r", "bifemsh_r", "psoas_r", "glut_max2_r",
           "tfl_r", "per_brev_r"]


def read_sto(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    t, names, vals = kine_ref._read_mot(path)
    return t, names, vals


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    z = np.load(HERE / "bsolve_out.npz", allow_pickle=True)
    qpos = z["qpos"]                       # [T, nq] mj qpos per frame
    t_ik = z["t"]
    T = len(t_ik)

    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    model = R.apply_harness(model, data)
    data = mujoco.MjData(model)
    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    aid = {n: i for i, n in enumerate(acts)}

    # --- replay frames, act=1, collect mj force + length ---
    F_mj = np.zeros((T, len(MUSCLES)))
    data.act[:] = 0.0
    for k in range(T):
        data.qpos[:] = qpos[k]
        B.apply_eq_followers(model, data)
        mujoco.mj_forward(model, data)
        data.act[:] = 0.0
        for j, m in enumerate(MUSCLES):
            data.act[aid[m]] = 1.0
            mujoco.mj_forward(model, data)
            F_mj[k, j] = abs(float(data.actuator_force[aid[m]]))
            data.act[aid[m]] = 0.0
    data.act[:] = 0.0

    # --- OpenSim muscle forces (SO's own force output; the MuscleAnalysis
    # TendonForce file ran control-less for many muscles - soleus et al.
    # are zero there, verified by _sto_diag.py) ---
    t_os, os_names, os_vals = read_sto(
        OSIM_DIR / "zz_bsolve_StaticOptimization_force.sto")
    t_act, act_names_os, act_vals = read_sto(
        OSIM_DIR / "zz_bsolve_StaticOptimization_activation.sto")

    # --- replay pass 2: matched activation (MuJoCo driven with OpenSim
    # SO activations - differences now isolate the FORMULATION) ---
    F_mj_so = np.zeros((T, len(MUSCLES)))
    for k in range(T):
        data.qpos[:] = qpos[k]
        B.apply_eq_followers(model, data)
        mujoco.mj_forward(model, data)
        data.act[:] = 0.0
        for j, m in enumerate(MUSCLES):
            a = float(np.interp(t_ik[k], t_act,
                                act_vals[:, act_names_os.index(m) - 1]))
            data.act[aid[m]] = min(max(a, 0.0), 1.0)
        mujoco.mj_forward(model, data)
        for j, m in enumerate(MUSCLES):
            F_mj_so[k, j] = abs(float(data.actuator_force[aid[m]]))
    data.act[:] = 0.0
    F_os = np.zeros((T, len(MUSCLES)))
    for j, m in enumerate(MUSCLES):
        # os_names INCLUDES the time column; vals excludes it -> -1
        col = os_names.index(m) - 1
        F_os[:, j] = np.interp(t_ik, t_os, os_vals[:, col])

    # --- stats ---
    print(f"dt between frames: {np.median(np.diff(t_ik)) * 1000:.1f} ms")
    print(f"{'muscle':12s} {'R2so':>7s} {'shape':>7s} {'lag[f]':>7s} "
          f"{'mj/os':>7s} | {'R2max':>7s} (act=1 capacity curve)")
    rows = []
    for j, m in enumerate(MUSCLES):
        b = F_os[:, j]

        def stats(a):
            ss_res = float(np.sum((a - b) ** 2))
            ss_tot = float(np.sum((b - b.mean()) ** 2))
            r2 = 1.0 - ss_res / max(ss_tot, 1e-9)
            sa = (a - a.mean()) / max(a.std(), 1e-9)
            sb = (b - b.mean()) / max(b.std(), 1e-9)
            cc = np.correlate(sa, sb, mode="full")
            mid = len(sb) - 1
            lag = int(np.argmax(cc[mid - T // 3: mid + T // 3 + 1])
                      - T // 3)
            r2s = float(np.corrcoef(sa, sb)[0, 1] ** 2)
            return r2, r2s, lag

        r2, r2s, lag = stats(F_mj_so[:, j])
        r2m, _, _ = stats(F_mj[:, j])
        ratio = float(F_mj_so[:, j].mean() / max(b.mean(), 1e-9))
        rows.append((m, r2, r2s, lag, ratio, float(b.max()),
                     float(F_mj_so[:, j].max()), r2m))
    print(f"{'':26s}(matched SO activation)        "
          f"(act=1 isometric)")
    for m, r2, r2s, lag, ratio, po, pm, r2m in rows:
        print(f"{m:12s} {r2:7.3f} {r2s:7.3f} {lag:7d} {ratio:7.2f} | "
              f"{r2m:7.3f}")

    with open(HERE / "muscle_force_compare.csv", "w", encoding="utf-8") as f:
        f.write("muscle,R2_matched,R2shape_matched,lag_frames,"
                "mj_over_os_mean,R2_act1\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")

    # --- figure: matched-activation comparison, normalized to own peak ---
    fig, axes = plt.subplots(3, 2, figsize=(9.5, 8.0), sharex=True)
    fig.subplots_adjust(hspace=0.16, left=0.08, right=0.98, top=0.93,
                        bottom=0.07)
    for k, m in enumerate(("soleus_r", "med_gas_r", "vas_lat_r",
                           "semimem_r", "tib_ant_r", "psoas_r")):
        j = MUSCLES.index(m)
        ax = axes[k // 2][k % 2]
        pc = (t_ik - t_ik[0]) / (t_ik[-1] - t_ik[0]) * 100.0
        ax.plot(pc, F_mj_so[:, j] / max(F_mj_so[:, j].max(), 1e-9),
                lw=1.8, color="#d55e00", label="MuJoCo (rigid tendon)")
        ax.plot(pc, F_os[:, j] / max(F_os[:, j].max(), 1e-9), "--",
                color="black", lw=1.6, label="OpenSim Thelen (compliant)")
        ax.set_title(f"{m} (R2 {rows[j][1]:.2f}, lag {rows[j][3]}f)", fontsize=8)
        ax.grid(True, alpha=0.25, lw=0.4)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6.5)
    for ax in axes[-1]:
        ax.set_xlabel("% of analyzed walk", fontsize=8)
    fig.suptitle("Muscle force along subject01 walk at matched OpenSim SO\nactivations (normalized to own peak)", fontsize=10)
    import matplotlib.pyplot as plt
    fig.savefig(HERE / "figures" / "muscle_force_compare.png", dpi=170)
    print("saved figures/muscle_force_compare.png")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401
    main()
