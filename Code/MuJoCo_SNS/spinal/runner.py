"""Run the spinal network on the converted Gait2392 MuJoCo model.

Schedule (params.SCHEDULE): quiet standing (balance feedback on) -> DRIVE
ramps up -> walking -> DRIVE ramps down -> back to quiet standing.

Each step:
  1. read MuJoCo: tendon lengths/velocities/forces per muscle, pelvis COM,
  2. encode afferents (Ia velocity, II length, Ib force) with presynaptic
     phase/speed gating derived from the RG and DRIVE potentials,
  3. advance the SNS network one dt,
  4. motoneuron potentials -> muscle activations (V / E_HI clipped [0,1]),
  5. mj_step.

Usage (from Code/MuJoCo_SNS/spinal, myo env):
    python runner.py [--drive N] [--no-afferents] [--harness K] [--time S]
Outputs spinal_run.npz + spinal_run.png next to this file and prints an
honest summary (did it stand, did it step, did it fall).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import build_network as bn
from muscle_map import EXTENSOR_STANCE_GROUPS, classify
from params import AFF, BAL, DT, E_HI, MOD, SCHEDULE

import mujoco

from scipy.optimize import nnls

HERE = Path(__file__).parent
MODEL = Path(r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
             r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml")

KEY_ACTS = ("vas_lat_r", "med_gas_r", "soleus_r", "tib_ant_r", "psoas_r",
            "semimem_r", "glut_max2_r", "rect_fem_r",
            "vas_lat_l", "med_gas_l", "soleus_l", "tib_ant_l", "psoas_l",
            "semimem_l", "glut_max2_l", "rect_fem_l")
KEY_JOINTS = ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
              "hip_flexion_l", "knee_angle_l", "ankle_angle_l")

WARMUP = 0.6   # s: pose pinned while muscle activations build from zero


def find_foot(model) -> tuple[int, int]:
    ids = []
    for frag in ("toes", "foot", "talus", "calcn"):
        ids = [i for i in range(model.nbody)
               if frag in mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i)]
        if len(ids) >= 2:
            break
    r = [i for i in ids if mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i).endswith("_r")]
    l = [i for i in ids if mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i).endswith("_l")]
    return r[0], l[0]


def drive_posture(t: float, walk_drive: float) -> tuple[float, float]:
    """Piecewise stand->walk->stand schedule -> (DRIVE nA, POSTURE nA).

    POSTURE is now only a small tonic bias on top of the solved standing
    activations (POST_<muscle> inputs carry the real posture pattern).
    """
    s = SCHEDULE
    if t < s["stand1"][1]:
        return 0.0, 1.0
    if t < s["ramp_up"][1]:
        f = (t - s["ramp_up"][0]) / (s["ramp_up"][1] - s["ramp_up"][0])
        return walk_drive * f, 1.0 - 0.4 * f
    if t < s["walk"][1]:
        return walk_drive, 0.6
    if t < s["ramp_down"][1]:
        f = (t - s["ramp_down"][0]) / (s["ramp_down"][1] - s["ramp_down"][0])
        return walk_drive * (1.0 - f), 0.6 + 0.4 * f
    return 0.0, 1.0


def solve_standing_activations(model, data, Fmax: np.ndarray) -> np.ndarray:
    """Static-optimization standing solution at the keyframe pose.

    The torque muscles must supply is what gravity demands MINUS what the
    ground contact and passive forces already provide:
        tau = qfrc_bias - qfrc_constraint - qfrc_passive   (all DoFs;
        ground contact covers pelvis translation). NNLS over non-negative
        activations with a small ridge; muscles with big Fmax and good
        moment arms (soleus, quads) win automatically.
    """
    mujoco.mj_forward(model, data)
    A = data.actuator_moment.T.copy()          # [nv, nu]
    tau = (data.qfrc_bias - data.qfrc_constraint - data.qfrc_passive).copy()
    # tonic co-contraction preload: net extension torque at hips/knees,
    # plantarflexion at ankles, extension at lumbar (signs follow the
    # OpenSim conventions of this model) - standing has real muscle tone,
    # it is not a zero-torque equilibrium
    def jadr(jn):
        return model.joint(jn).qposadr[0]
    for side in ("r", "l"):
        tau[jadr(f"hip_flexion_{side}")] += -25.0   # hip extension
        tau[jadr(f"knee_angle_{side}")] += -30.0    # knee extension
        tau[jadr(f"ankle_angle_{side}")] += 20.0    # plantarflexion
    for jn in ("lumbar_extension",):
        tau[jadr(jn)] += -10.0

    # prior-regularized least squares: fit the torque rows but stay close to
    # a physiological standing pattern (plain NNLS degenerates to one-muscle
    # vertex solutions; a sparsity ridge is exactly what we do NOT want)
    from build_network import _group_weight
    from muscle_map import classify
    from params import POSTURE_OVERRIDE, W_POSTURE
    x_prior = np.array([
        min(0.5, 0.6 * _group_weight(classify(a), W_POSTURE, POSTURE_OVERRIDE))
        for a in (mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
                  for i in range(model.nu))])
    Aff = A * Fmax[None, :] * 1e-2             # scale for conditioning
    lam = 0.5
    Aaug = np.vstack([Aff, lam * np.eye(model.nu)])
    baug = np.concatenate([tau * 1e-2, lam * x_prior])
    x, _ = nnls(Aaug, baug)
    return np.clip(x, 0.0, 1.0)


def apply_harness(model, data, kxy=1.0e6, kz=1.0e6, ky=1.0e6,
                  dxy=2000.0, dy=2000.0):
    """Rigid pelvis rig (v1 default): the pelvis is spring-locked at the
    keyframe pose in all three translations - the classic biped walker test
    rig. Legs swing under it with ground contact; balance is NOT solved.
    Weaker springs (a compliant rehab tether) explode the stiff contact
    model at dt=5 ms; see DESIGN.md open problems.
    """
    import re
    xml_text = MODEL.read_text(encoding="utf-8")
    # stiff joint springs + explicit Euler at dt=5 ms explode; implicitfast
    # handles the joint-damping/spring terms stably
    xml_text = xml_text.replace(
        '<option timestep="0.005" collision="predefined"/>',
        '<option timestep="0.005" collision="predefined" integrator="implicitfast"/>')
    # converted bone inertias violate the triangle inequality on some bodies
    # (singular mass matrix -> QACC explosions); let MuJoCo rebalance them
    xml_text = xml_text.replace(
        '<compiler angle="radian" autolimits="true"/>',
        '<compiler angle="radian" autolimits="true" balanceinertia="true"/>')
    # model repair (mechanics): OpenSim conditional pathpoints were converted
    # to slide joints on MASSLESS bodies (mass 0, inertia 5e-14) -> singular
    # mass matrix, NaN blowups. These points are computed projections in
    # OpenSim, not physical masses; weld them as fixed via-points at the
    # keyframe geometry (MyoSuite-style). 82 DoFs removed; upgrade path is a
    # per-pose spline treatment, see DESIGN.md.
    xml_text = re.sub(r'<joint name="[^"]*-P\d+_[xyz]"[^/>]*/>', '', xml_text)
    # ...and the equality couplings that drove them (polycoef splines)
    xml_text = re.sub(r'<joint joint1="[^"]*-P\d+_[xyz]"[^/>]*/>', '', xml_text)
    # the recorded keyframe has the old 105-DoF size; pose is seeded
    # programmatically instead (capture_pose/seed_pose)
    xml_text = re.sub(r'<key [^/>]*/>', '', xml_text)
    specs = {"pelvis_tx": (kxy, dxy), "pelvis_ty": (ky, dy),
             "pelvis_tz": (kz, dxy)}
    for jn, (k, c) in specs.items():
        ref = data.qpos[model.joint(jn).qposadr[0]]
        pat = rf'(<joint name="{jn}"[^/]*?)/>'
        # keep the existing damping attribute; just add stiffness + springref
        xml_text = re.sub(pat, rf'\1 stiffness="{k}" springref="{ref}"/>',
                          xml_text, count=1)
    # model repair: the converted knee translation slides carry ZERO stiffness
    # (OpenSim couples them kinematically to knee_angle; the MJCF leaves them
    # as free slides with ~2.7 cm of play) -> the tibia/femur linkage buckles
    # under any load. Spring them to their keyframe values to restore the
    # moving-pathpoint coupling mechanically.
    for jn in ("knee_r_translation1", "knee_r_translation2",
               "knee_l_translation1", "knee_l_translation2"):
        ref = data.qpos[model.joint(jn).qposadr[0]]
        pat = rf'(<joint name="{jn}"[^/]*?)/>'
        xml_text = re.sub(
            pat, rf'\1 stiffness="5000" springref="{ref}" damping="100"/>',
            xml_text, count=1)
    # mesh files are referenced relative to the model dir -> must load the
    # patched XML from a file in that directory, not from a string
    import tempfile, os
    fd, tmp = tempfile.mkstemp(suffix=".xml", dir=str(MODEL.parent))
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(xml_text)
    try:
        return mujoco.MjModel.from_xml_path(tmp)
    finally:
        os.unlink(tmp)


def capture_pose(model) -> dict[str, np.ndarray]:
    """Joint name -> keyframe qpos values (for re-seeding reduced models)."""
    pose = {}
    for j in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        adr, narange = model.jnt_qposadr[j], 7 if model.jnt_type[j] == 0 else 1
        pose[name] = model.key_qpos[0][adr:adr + narange].copy()
    return pose


def seed_pose(model, data, pose: dict[str, np.ndarray]):
    mujoco.mj_resetData(model, data)
    for j in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        if name in pose:
            adr = model.jnt_qposadr[j]
            data.qpos[adr:adr + len(pose[name])] = pose[name]


def main(argv):
    walk_drive = 4.0
    harness = 1500.0     # default ON (balance assist); 0 disables
    use_aff = True
    args = list(argv)
    while args:
        a = args.pop(0)
        if a == "--drive":
            walk_drive = float(args.pop(0))
        elif a == "--harness":
            harness = float(args.pop(0)) if args and args[0].replace(".", "").isdigit() else 1500.0
        elif a == "--no-harness":
            harness = 0.0
        elif a == "--no-afferents":
            use_aff = False
        elif a == "--time":
            SCHEDULE["stand2"] = (SCHEDULE["ramp_down"][1],
                                  SCHEDULE["ramp_down"][1] + float(args.pop(0)))

    model = mujoco.MjModel.from_xml_path(str(MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    key_pose = capture_pose(model)
    if harness > 0:
        model = apply_harness(model, data, kxy=harness)
        data = mujoco.MjData(model)
        seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)

    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    net = bn.build(acts, dt=DT)
    aid = {name: i for i, name in enumerate(acts)}
    muscles = net.muscles

    # per-muscle normalization: length range and peak force. The converted
    # actuators are <general class="muscle"> with the muscle parameter vector
    # in gainprm/biasprm; Fmax sits at gainprm[2] (e.g. soleus 3395 N).
    Lrange = model.actuator_lengthrange  # [nu, 2]
    Lmid = Lrange.mean(axis=1)
    Lhalf = np.maximum((Lrange[:, 1] - Lrange[:, 0]) / 2, 1e-3)
    Fmax = np.maximum(model.actuator_gainprm[:, 2], 5.0)
    act_stand = solve_standing_activations(model, data, Fmax)
    print(f"standing solve: nnls x in [{act_stand.min():.2f}, "
          f"{act_stand.max():.2f}], >0.1: {(act_stand > 0.1).sum()} muscles")
    top = np.argsort(act_stand)[::-1][:8]
    print("  top: " + ", ".join(f"{acts[i]}={act_stand[i]:.2f}" for i in top))
    print(f"Fmax range: {Fmax.min():.0f}..{Fmax.max():.0f} N | "
          f"L range sample: {Lrange[aid['soleus_r']]}")

    iport = {p: net.input_index(p) for p in net.inputs}
    foot_r, foot_l = find_foot(model)
    # balance target = ankle axis (talus) x: the neutral standing COM
    ankle_x = []
    for i in range(model.nbody):
        nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i)
        if nm and nm.startswith("talus"):
            ankle_x.append(float(data.xpos[i][0]))
    x_ref = float(np.mean(ankle_x)) if ankle_x else float(data.subtree_com[0][0])
    jadr = {jn: model.jnt_qposadr[model.joint(jn).id] for jn in KEY_JOINTS
            if model.joint(jn).id >= 0}

    dur = SCHEDULE["stand2"][1]
    nsteps = int(dur / DT)
    log_t = np.zeros(nsteps)
    log_act = np.zeros((nsteps, len(KEY_ACTS)))
    log_q = np.zeros((nsteps, len(KEY_JOINTS)))
    log_com = np.zeros((nsteps, 3))
    log_neuro = np.zeros((nsteps, 8))  # DRIVE, POST, RG_E_r, RG_F_r, PF_E2_r, PF_F1_r, RG_E_l, RG_F_l
    bal_neurons = ("BAL_PF", "BAL_DF")

    u = net.make_inputs()
    stance_ext = np.array([muscles[a].groups[0] in EXTENSOR_STANCE_GROUPS
                           for a in acts])
    for k in range(nsteps):
        t = k * DT
        drive, posture = drive_posture(t, walk_drive)
        # --- afferent encoding ---
        if use_aff:
            L = data.actuator_length
            Ldot = data.actuator_velocity
            F = np.abs(data.actuator_force)
            len_norm = np.clip((L - Lmid) / Lhalf, -1.0, 1.0)
            vel_norm = np.clip(Ldot / AFF["ia_vel_ref"], -1.0, 1.0)
            force_norm = np.clip(F / Fmax, 0.0, 1.5)
            v = net.compiled.V
            drive_frac = float(np.clip(v[net.idx["DRIVE"]] / E_HI, 0, 1))
            stance = {s: float(np.clip(v[net.idx[f"RG_E_{s}"]] / E_HI, 0, 1))
                      for s in net.sides}
            mod_ia = 1.0 + MOD["ia"] * drive_frac
            mod_ii = 1.0 + MOD["ii"] * drive_frac
            g_ia = np.where(stance_ext,
                            AFF["ia_gain"] * mod_ia * (0.5 + 0.5 * np.array(
                                [stance[muscles[a].side] for a in acts])),
                            AFF["ia_gain"] * mod_ia)
            g_ii = np.where(stance_ext,
                            AFF["ii_gain"] * mod_ii * (0.3 + 0.7 * np.array(
                                [stance[muscles[a].side] for a in acts])),
                            0.5 * AFF["ii_gain"] * mod_ii)
            g_ib = np.where(stance_ext,
                            AFF["ib_gain"] * (1.0 - 0.7 * np.array(
                                [stance[muscles[a].side] for a in acts])),
                            AFF["ib_gain"])
            for i, a in enumerate(acts):
                s = muscles[a].side
                u[iport[f"Ia_{a}"]] = max(g_ia[i] * vel_norm[i], 0.0)
                u[iport[f"II_{a}"]] = max(AFF["i0_ii"] + g_ii[i] * len_norm[i], 0.0)
                u[iport[f"Ib_{a}"]] = max(0.5 * g_ib[i] * force_norm[i], 0.0)
        # --- balance feedback (standing only; fades with drive) ---
        fade = (1.0 - drive) / max(walk_drive, 1e-6) if BAL["fade_with_drive"] else 1.0
        com = data.subtree_com[0]
        u_pd = BAL["kx"] * (x_ref - com[0]) - BAL["kv"] * data.qvel[0]
        u_ff = BAL["ff"] * np.exp(-max(t - WARMUP, 0.0) / BAL["ff_tau"])
        u_pd = np.clip(BAL["kx"] * (x_ref - com[0]) - BAL["kv"] * data.qvel[0],
                       -BAL["max_current"], BAL["max_current"]) * fade
        u_pd = float(np.clip(u_pd + u_ff, -1.5 * BAL["max_current"],
                             1.5 * BAL["max_current"]))
        u[iport["BAL_PF"]] = max(u_pd, 0.0)
        u[iport["BAL_DF"]] = max(-u_pd, 0.0)
        u[iport["DRIVE"]] = drive
        u[iport["POSTURE"]] = posture
        # solved standing pattern -> motoneuron posture bias (fades to a
        # floor during walking so the stepping pattern can take over)
        stand_frac = 0.20 + 0.80 * (1.0 - min(drive / max(walk_drive, 1e-6), 1.0))
        for i, a in enumerate(acts):
            u[iport[f"POST_{a}"]] = 6.0 * act_stand[i] * stand_frac

        # --- neural step + motor mapping ---
        v = net.step(u)
        for i, a in enumerate(acts):
            data.ctrl[aid[a]] = np.clip(v[net.idx[net.mn_names[a]]] / E_HI, 0.0, 1.0)

        if t < WARMUP:
            # hold the keyframe pose while activations build (qpos untouched,
            # velocities zeroed, activation states integrate via mj_forward)
            data.qvel[:] = 0.0
            mujoco.mj_forward(model, data)
        else:
            mujoco.mj_step(model, data)

        # --- log ---
        log_t[k] = t
        log_com[k] = com
        for i, a in enumerate(KEY_ACTS):
            log_act[k, i] = data.ctrl[aid[a]]
        for i, jn in enumerate(KEY_JOINTS):
            log_q[k, i] = data.qpos[jadr[jn]]
        log_neuro[k] = (drive, posture,
                        v[net.idx["RG_E_r"]], v[net.idx["RG_F_r"]],
                        v[net.idx["PF_E2_r"]], v[net.idx["PF_F1_r"]],
                        v[net.idx["RG_E_l"]], v[net.idx["RG_F_l"]])

    np.savez_compressed(HERE / "spinal_run.npz", t=log_t, act=log_act,
                        q=np.degrees(log_q), com=log_com, neuro=log_neuro,
                        key_acts=KEY_ACTS, key_joints=KEY_JOINTS)

    # ---------------- honest summary ----------------
    pelz = log_com[:, 2]
    xcom = log_com[:, 0]
    fall = pelz.min() < 0.55
    ifall = np.argmax(pelz < 0.55) if fall else -1
    walk_mask = log_t > SCHEDULE["walk"][0]
    n_steps = int(np.sum(np.diff((log_neuro[walk_mask, 3] > 0.8).astype(int)) == 1))
    print(f"\n== summary ==")
    print(f"pelvis/COM height: start {pelz[0]:.2f} m, min {pelz.min():.2f} m, "
          f"final {pelz[-1]:.2f} m  ({'FELL at t=%.1f s' % log_t[ifall] if fall else 'stayed up'})")
    print(f"COM x: min {xcom.min():.2f} max {xcom.max():.2f} final {xcom[-1]:.2f} m")
    print(f"RG_F_r bursts during walk window: {n_steps}")
    print(f"peak ctrl: {log_act.max(axis=0).round(2)} ({KEY_ACTS})")
    print(f"knee_r range: {log_q[:, 1].min():.0f}..{log_q[:, 1].max():.0f} deg, "
          f"ankle_r: {log_q[:, 2].min():.0f}..{log_q[:, 2].max():.0f} deg")
    print("saved spinal_run.npz")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
        for i, a in enumerate(KEY_ACTS):
            ax[0].plot(log_t, log_act[:, i], lw=0.8,
                       label=a if i < 8 or a.endswith("_l") else None)
        ax[0].set_ylabel("activation")
        for i, jn in enumerate(KEY_JOINTS):
            ax[1].plot(log_t, log_q[:, i], lw=0.8,
                       label=jn if i < 3 else None)
        ax[1].set_ylabel("joint angle [deg]")
        ax[1].legend(fontsize=7, ncol=3)
        ax[2].plot(log_t, log_com[:, 0], label="com x")
        ax[2].plot(log_t, log_com[:, 2], label="com z")
        ax[2].set_ylabel("COM [m]"); ax[2].legend(fontsize=7)
        for j, nm in enumerate(("DRIVE", "POSTURE", "RG_E_r", "RG_F_r",
                                "PF_E2_r", "PF_F1_r", "RG_E_l", "RG_F_l")):
            ax[3].plot(log_t, log_neuro[:, j], lw=0.8, label=nm)
        ax[3].set_ylabel("neural [mV]"); ax[3].legend(fontsize=7, ncol=4)
        ax[3].set_xlabel("t [s]")
        for a in ax:
            a.grid(True, alpha=0.3)
        fig.suptitle("gait2392 + spinal SNS: stand -> walk -> stand")
        fig.tight_layout()
        fig.savefig(HERE / "spinal_run.png", dpi=130)
        print("saved spinal_run.png")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main(sys.argv[1:])
