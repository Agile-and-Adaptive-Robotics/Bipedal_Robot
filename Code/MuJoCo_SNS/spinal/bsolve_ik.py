"""IK trajectory -> MuJoCo: muscle validation vs OpenSim + NNLS back-solve.

Implements Ben's staged plan (DESIGN.md 2026-09-11 "NEXT"):
OpenSim IK walking kinematics -> map onto the converted MJCF joints ->
(1) MuJoCo-vs-OpenSim muscle validation along the trajectory (muscle
lengths; OpenSim reference computed by opensim-cmd: MuscleAnalysis +
StaticOptimization on subject01_simbody.osim) ->
(2) per-timestep NNLS activation back-solve: for every IK frame solve
    min ||a>=0:  M(q) (unit_force .* a)  ~=  tau_ID - M(q) f_passive
where tau_ID is MuJoCo inverse dynamics (measured GRF applied at the
feet via Jacobians; SO-style quasi-static: force-velocity ignored),
M(q) = actuator_moment, f_passive the act=0 muscle forces. The NNLS
activations are therefore directly comparable to OpenSim SO's.

Frame mapping: the MuJoCo model is z-up but the converter PRESERVED
OpenSim's coordinate values (pelvis slides load identity - verified
against the keyframe), while world VECTORS (GRF force/point/torque)
remap os(x, y, z) -> mj(x, -z, y)  (the runner's BAL_LAT note
"mujoco -y = opensim +z"). Rotational coordinate SIGNS are
auto-detected in two stages: greedy matching of muscle lengths vs the
OpenSim MuscleAnalysis reference, then an ID+GRF frontal-consistency
sweep (hip adduction/subtalar/pelvis-list/hip-rotation flips have weak
length signatures but huge frontal-torque signatures).

Outputs (in this folder):
    bsolve_out.npz          everything for fit_pf.py / plotting
    bsolve_lengths.png      MuJoCo vs OpenSim muscle-length trajectories
    bsolve_activations.png  NNLS vs OpenSim SO
    bsolve_groups.png       functional-group activation vs gait cycle
    bsolve_report.txt       per-muscle verdict table
    bsolve_osim.log         raw opensim-cmd output
Usage: python bsolve_ik.py [--skip-osim] [--subsample N]
(cwd must be Code/MuJoCo_SNS/spinal; myo env)
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

import runner  # provides MODEL + the full patch/repair stack
from muscle_map import classify, GROUPS

import mujoco
from scipy.optimize import nnls
from scipy.signal import savgol_filter

HERE = Path(__file__).parent
OSIM_DIR = Path(r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody")
IK_MOT = OSIM_DIR / "subject01_walk1_ik.mot"
GRF_MOT = OSIM_DIR / "subject01_walk1_grf.mot"
OSIM_MODEL = "subject01_simbody.osim"          # relative to OSIM_DIR
OSIM_CMD = r"D:\OpenSim 4.3\bin\opensim-cmd.exe"

ROT_DRIVERS = (  # OpenSim rotational coordinates (deg in the .mot)
    "pelvis_tilt", "pelvis_list", "pelvis_rotation",
    "hip_flexion_r", "hip_adduction_r", "hip_rotation_r", "knee_angle_r",
    "ankle_angle_r", "subtalar_angle_r", "mtp_angle_r",
    "hip_flexion_l", "hip_adduction_l", "hip_rotation_l", "knee_angle_l",
    "ankle_angle_l", "subtalar_angle_l", "mtp_angle_l",
    "lumbar_extension", "lumbar_bending", "lumbar_rotation")
TRANS_DRIVERS = ("pelvis_tx", "pelvis_ty", "pelvis_tz")
FIT_JOINTS = (  # dofs the NNLS fits: leg + lumbar joints (pelvis rig-carried)
    "hip_flexion_r", "hip_adduction_r", "hip_rotation_r", "knee_angle_r",
    "ankle_angle_r", "subtalar_angle_r", "mtp_angle_r",
    "hip_flexion_l", "hip_adduction_l", "hip_rotation_l", "knee_angle_l",
    "ankle_angle_l", "subtalar_angle_l", "mtp_angle_l",
    "lumbar_extension", "lumbar_bending", "lumbar_rotation")


def read_mot(path: Path):
    """Generic .mot/.sto parser -> (t [T], colnames(no time), data [T, n])."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    end = next(i for i, ln in enumerate(lines)
               if ln.strip() == "endheader")
    names = lines[end + 1].split()
    rows = []
    for ln in lines[end + 2:]:
        s = ln.strip()
        if not s:
            continue
        try:
            rows.append([float(v) for v in s.split()])
        except ValueError:
            break
    data = np.array(rows)
    return data[:, 0], names[1:], data[:, 1:]


def build_air_model(id_mode=False):
    """The converted model with ALL repairs, no ground contact, free pelvis
    (zero-stiffness rig) - the same repair stack the live runs use.
    id_mode: additionally zero the armature on the equality-coupled
    pathpoint follower dofs. The armature=1.0 was added for FORWARD-sim
    stability; for inverse dynamics it is ~36 kg*m^2 of fictitious
    reflected inertia at the knee that inflates tau = M qacc with
    differentiation noise - the human torque demand does not contain it."""
    m0 = mujoco.MjModel.from_xml_path(str(runner.MODEL))
    d0 = mujoco.MjData(m0)
    mujoco.mj_resetDataKeyframe(m0, d0, 0)
    model = runner.apply_harness(m0, d0, kxy=0, kz=0, ky=0, dxy=1.0, dy=1.0,
                                 no_ground=True, pin_rot=False)
    if id_mode:
        # 1. zero the armature the forward-sim stability patch added to the
        #    equality-coupled follower dofs (fictitious inertia for ID).
        # 2. DISABLE the equalities entirely: at deep knee flexion the
        #    polyfit couplers (fitted near the straight keyframe) generate
        #    constraint wrenches of hundreds of N*m at the knee row
        #    (diag_knee.py: -890 N*m static at a swing pose). With the
        #    equalities off, inverse dynamics is the clean TREE system
        #    moving along the measured trajectory: the followers keep
        #    their poly-projected positions (set in set_frame), their own
        #    dof rows are excluded from the fit, and no constraint force
        #    contaminates the driver rows.
        n = 0
        for i in range(model.neq):
            if model.eq_type[i] != mujoco.mjtEq.mjEQ_JOINT:
                continue
            adr = model.jnt_dofadr[int(model.eq_obj1id[i])]
            model.dof_armature[adr] = 0.005
            n += 1
        model.opt.disableflags |= int(mujoco.mjtDisableBit.mjDSBL_EQUALITY)
        print(f"ID model: {n} follower armatures zeroed, equalities "
              f"disabled")
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    return model, data


def apply_eq_followers(model, data):
    """Project equality-coupled pathpoint dofs exactly onto their drivers
    (mjEQ_JOINT: y = y0 + poly(x - x0), model.qpos0 offsets)."""
    for i in range(model.neq):
        if model.eq_type[i] != mujoco.mjtEq.mjEQ_JOINT:
            continue
        j_dep = int(model.eq_obj1id[i])
        j_ind = int(model.eq_obj2id[i])
        adr_d = model.jnt_qposadr[j_ind]
        adr_p = model.jnt_qposadr[j_dep]
        x = data.qpos[adr_d] - float(model.qpos0[adr_d])
        c = model.eq_data[i][:4]
        data.qpos[adr_p] = float(model.qpos0[adr_p]) + c[0] + c[1] * x \
            + c[2] * x * x + c[3] * x * x * x


def set_frame(model, data, signs, values):
    """Write one IK frame into qpos.

    `values` maps OpenSim coordinate name -> value (rotations in DEGREES,
    translations in METERS, OpenSim world frame). The converter PRESERVED
    OpenSim's coordinate values on the pelvis slides (keyframe check:
    qpos[pelvis_ty]=0.95 <-> pelvis world z 0.95; pelvis_tz=+0.1 moves the
    pelvis to world y=-0.1 = os +z) - so translations load IDENTITY, and
    the z-up remap lives in the body frames. GRF vectors still go through
    osim_to_mj_vec. Rotational signs are auto-detected.
    """
    mujoco.mj_resetDataKeyframe(model, data, 0)
    for name in ROT_DRIVERS:
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
        if jid < 0:
            raise RuntimeError(f"joint {name} missing from MJCF")
        adr = model.jnt_qposadr[jid]
        data.qpos[adr] = np.deg2rad(values[name]) * signs.get(name, 1.0)
    for name in TRANS_DRIVERS:
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
        data.qpos[model.jnt_qposadr[jid]] = values[name]
    apply_eq_followers(model, data)


def osim_to_mj_vec(v):
    """OpenSim world vector -> MuJoCo world vector."""
    return np.array([v[0], -v[2], v[1]])


def grf_at(t_frames, tt, grf_names, grf, pref=""):
    """Interpolated plate force [T,3], cop [T,3], free torque [T,3] in the
    MuJoCo frame."""
    col = {n: i for i, n in enumerate(grf_names)}
    def ip(nm):
        return np.column_stack([np.interp(t_frames, tt, grf[:, col[f"{pref}{nm}{a}"]])
                                for a in "xyz"])
    F, P, Tq = ip("ground_force_v"), ip("ground_force_p"), ip("ground_torque_")
    return (np.array([osim_to_mj_vec(r) for r in F]),
            np.array([osim_to_mj_vec(r) for r in P]),
            np.array([osim_to_mj_vec(r) for r in Tq]))


def length_err(model, data, signs, ik, frames, L_os, act_names):
    """Median per-muscle normalized length RMSE vs the OpenSim reference."""
    col_os = {n: i for i, n in enumerate(L_os["names"])}
    ref_all = L_os["data"]
    sel = np.asarray(frames)
    sel = sel[sel < ref_all.shape[0]]
    L = np.zeros((len(sel), model.nu))
    for kk, k in enumerate(sel):
        set_frame(model, data, signs, {n: ik[n][k] for n in ik})
        mujoco.mj_forward(model, data)
        L[kk] = data.actuator_length
    errs = []
    for i, a in enumerate(act_names):
        j = col_os.get(a)
        if j is None:
            continue
        ref = ref_all[sel, j]
        rng = max(np.ptp(ref), 1e-3)
        errs.append(np.sqrt(np.mean((L[:, i] - ref) ** 2)) / rng)
    return float(np.median(errs)) if errs else 9.9


def act_names_of(model):
    return [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]


def _run_osim_tool(t_ik):
    """Write the AnalyzeTool XML (MuscleAnalysis + StaticOptimization as an
    ANALYSIS - 4.x has no standalone SO tool) and run it via opensim-cmd."""
    t0, t1 = float(t_ik[0]), float(t_ik[-1])
    an = OSIM_DIR / "zz_bsolve_analyze.xml"
    an.write_text(f"""<?xml version="1.0" encoding="UTF-8"?>
<OpenSimDocument Version="20302">
  <AnalyzeTool name="zz_bsolve">
    <model_file> {OSIM_MODEL} </model_file>
    <replace_force_set> false </replace_force_set>
    <force_set_files>  </force_set_files>
    <results_directory> ./ResultsBSolve </results_directory>
    <output_precision> 8 </output_precision>
    <initial_time> {t0:.6f} </initial_time>
    <final_time> {t1:.6f} </final_time>
    <solve_for_equilibrium_for_auxiliary_states> false </solve_for_equilibrium_for_auxiliary_states>
    <maximum_number_of_integrator_steps> 20000 </maximum_number_of_integrator_steps>
    <maximum_integrator_step_size> 1.0 </maximum_integrator_step_size>
    <minimum_integrator_step_size> 1e-08 </minimum_integrator_step_size>
    <integrator_error_tolerance> 1e-10 </integrator_error_tolerance>
    <AnalysisSet name="Analyses">
      <objects>
        <MuscleAnalysis name="MuscleAnalysis">
          <on> true </on>
          <start_time> -infinity </start_time>
          <end_time> infinity </end_time>
          <step_interval> 1 </step_interval>
          <in_degrees> true </in_degrees>
          <muscle_list> all </muscle_list>
          <moment_arm_coordinate_list> all </moment_arm_coordinate_list>
          <compute_moments> true </compute_moments>
        </MuscleAnalysis>
        <StaticOptimization name="StaticOptimization">
          <on> true </on>
          <start_time> -infinity </start_time>
          <end_time> infinity </end_time>
          <step_interval> 1 </step_interval>
          <in_degrees> true </in_degrees>
          <use_model_force_dof_velocity> false </use_model_force_dof_velocity>
        </StaticOptimization>
      </objects>
      <groups/>
    </AnalysisSet>
    <ControllerSet name="Controllers"><objects/><groups/></ControllerSet>
    <external_loads_file> subject01_walk1_grf.xml </external_loads_file>
    <states_file>  </states_file>
    <coordinates_file> subject01_walk1_ik.mot </coordinates_file>
    <speeds_file>  </speeds_file>
    <lowpass_cutoff_frequency_for_coordinates> 6 </lowpass_cutoff_frequency_for_coordinates>
  </AnalyzeTool>
</OpenSimDocument>
""", encoding="utf-8")
    print("opensim-cmd run-tool zz_bsolve_analyze.xml ...", flush=True)
    p = subprocess.run([OSIM_CMD, "run-tool", an.name], cwd=str(OSIM_DIR),
                       capture_output=True, text=True, timeout=1800)
    (HERE / "bsolve_osim.log").write_text(
        f"=== zz_bsolve_analyze rc={p.returncode}\n{p.stdout}\n{p.stderr}",
        encoding="utf-8")
    if p.returncode != 0:
        print(f"  FAILED rc={p.returncode} (bsolve_osim.log); tail:")
        print("\n".join((p.stdout + p.stderr).splitlines()[-15:]))


def run_osim_reference(t_ik, reuse=False):
    """MuscleAnalysis + StaticOptimization reference outputs, parsed."""
    out = OSIM_DIR / "ResultsBSolve"
    if reuse and sorted(out.glob("*MuscleAnalysis_Length.sto")) \
            and sorted(out.glob("*StaticOptimization_activation.sto")):
        print("reusing existing ResultsBSolve outputs")
    else:
        _run_osim_tool(t_ik)
    L_os = acts_os = None
    MA_os = {}
    lens = sorted(out.glob("*MuscleAnalysis_Length.sto"))
    sos = sorted(out.glob("*StaticOptimization_activation.sto"))
    import re
    for f in out.glob("*MuscleAnalysis_MomentArm_*.sto"):
        m = re.search(r"MomentArm_(.+)\.sto$", f.name)
        if m:
            _, nms, d_ = read_mot(f)
            MA_os[m.group(1)] = dict(names=nms, data=d_)
    if lens:
        t, nms, d = read_mot(lens[-1])
        L_os = dict(names=nms, data=d, t=t)
    if sos:
        t, nms, d = read_mot(sos[-1])
        acts_os = dict(names=nms, data=d, t=t)
    print(f"OpenSim reference: lengths {L_os is not None}, "
          f"SO {acts_os is not None}, moment arms {len(MA_os)} joints")
    return L_os, acts_os, MA_os


def frames_qpos(model, data, signs, ik, ik_names, t_ik):
    """qpos trajectory [T, nq] for one sign hypothesis."""
    qd = np.zeros((len(t_ik), model.nq))
    for k in range(len(t_ik)):
        set_frame(model, data, signs, {n: ik[n][k] for n in ik_names})
        qd[k] = data.qpos.copy()
    return qd


def id_cost(model, data, signs, ik, ik_names, qd, qvel, qacc, rows,
            Fr, Pr, Tr, Fl, Pl, Tl, t_ik, sample=4):
    """Frontal-plane ID residual under measured GRF: median |tau| at the
    hip-adduction + subtalar dofs over strongly-loaded frames (the frontal
    signs have weak muscle-length signatures but move the feet across the
    sagittal plane, so the CoP-to-foot geometry - and hence these torques -
    detects a flip decisively)."""
    cal_r = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_r")
    cal_l = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_l")
    loaded = np.flatnonzero(Fr[:, 2] > 400.0)[::sample]
    c_r, c_l = [], []
    jacp = np.zeros((3, model.nv))
    jacr = np.zeros((3, model.nv))
    for k in loaded:
        set_frame(model, data, signs, {n: ik[n][k] for n in ik_names})
        data.qpos[:] = qd[k]
        data.qvel[:] = qvel[k]
        data.qacc[:] = qacc[k]
        data.act[:] = 0.0
        mujoco.mj_inverse(model, data)
        tau = data.qfrc_inverse.copy()
        for (F, P, Tq, bid) in ((Fr, Pr, Tr, cal_r), (Fl, Pl, Tl, cal_l)):
            mujoco.mj_jac(model, data, jacp, jacr, P[k], bid)
            tau -= jacp.T @ F[k] + jacr.T @ Tq[k]
        c_r.append(abs(tau[rows[1]]))    # hip_adduction_r
        c_l.append(abs(tau[rows[8]]))    # hip_adduction_l
    return float(np.median(c_r) + np.median(c_l))


def fd_moments(model, data, qd, rows, delta=1.0e-4):
    """Muscle moment arms at the FIT dofs by central difference dl/dtheta,
    re-projecting the equality-coupled pathpoint followers at every
    perturbed pose.

    WHY NOT actuator_moment: MuJoCo's transmission Jacobian is the RAW
    tendon-length sensitivity w.r.t. the dof and does NOT propagate
    through the eq couplers - at the knee it misses the 36 vastii
    pathpoint followers (and at hip_flexion 5 more), which dominate the
    knee moment arm. OpenSim's arm IS dl/dtheta with the whole geometry
    following, so the FD here (Ben's definition note 2026-09-11) is the
    matching quantity. Returns [T, nrows, nu]."""
    T = qd.shape[0]
    A_fd = np.zeros((T, len(rows), model.nu))
    for k in range(T):
        data.qpos[:] = qd[k]
        data.act[:] = 0.0
        for jj, adr in enumerate(rows):
            for s in (1.0, -1.0):
                data.qpos[:] = qd[k]
                data.qpos[adr] += s * delta
                apply_eq_followers(model, data)
                mujoco.mj_forward(model, data)
                A_fd[k, jj] += s * data.actuator_length / (2.0 * delta)
    return A_fd


def smooth_derivs(qd, dt):
    """Zero-phase 6 Hz lowpass (the same filtering StaticOptimization uses,
    its setup's lowpass_cutoff_frequency_for_coordinates=6) + gradients ->
    (qvel, qacc)."""
    from scipy.signal import butter, filtfilt
    b, a = butter(4, 6.0 / (0.5 / dt))
    qf = filtfilt(b, a, qd, axis=0, padtype="odd", padlen=9)
    v = np.gradient(qf, dt, axis=0)
    return v, np.gradient(v, dt, axis=0)


def main(argv):
    skip_osim = "--skip-osim" in argv
    subsample = 5
    if "--subsample" in argv:
        subsample = int(argv[argv.index("--subsample") + 1])

    # ---------------------------------------------------------------- data
    t_ik, ik_names, ik_vals = read_mot(IK_MOT)
    ik = {n: ik_vals[:, i] for i, n in enumerate(ik_names)}
    for n in ROT_DRIVERS + TRANS_DRIVERS:
        assert n in ik, f"{n} missing from IK mot"
    t_g, g_names, g_vals = read_mot(GRF_MOT)
    T = len(t_ik)
    print(f"IK: {T} frames {t_ik[0]:.2f}..{t_ik[-1]:.2f} s "
          f"({t_ik[1] - t_ik[0]:.4f} s step); GRF {len(t_g)} rows")

    model, data = build_air_model(id_mode=True)
    acts = act_names_of(model)
    Fmax = np.maximum(model.actuator_gainprm[:, 2], 5.0)

    # ------------------------------------------- OpenSim reference analyses
    L_os = acts_os = MA_os = None
    if not skip_osim:
        L_os, acts_os, MA_os = run_osim_reference(
            t_ik, reuse="--reuse-osim" in argv)

    # ------------------------------------------------------- sign discovery
    signs = {n: 1.0 for n in ROT_DRIVERS}
    frames = np.arange(0, T, subsample)
    if L_os is not None:
        for pas in range(2):
            for j in ROT_DRIVERS:
                err = {}
                for s in (1.0, -1.0):
                    signs[j] = s
                    err[s] = length_err(model, data, signs, ik, frames,
                                        L_os, acts)
                signs[j] = min(err, key=err.get)
                if err[1.0] != err[-1.0]:
                    print(f"  sign {j:22s} -> {signs[j]:+.0f} "
                          f"(err {err[1.0]:.4f} vs {err[-1.0]:.4f})", flush=True)
        base = length_err(model, data, signs, ik, frames, L_os, acts)
        print(f"sign search done: median normalized length err {base:.4f}")

    # ----------------------------------------------- frontal sign sweep
    # GRF (needed here and in the ID pass below)
    Fr, Pr, Tr = grf_at(t_ik, t_g, g_names, g_vals, pref="")
    Fl, Pl, Tl = grf_at(t_ik, t_g, g_names, g_vals, pref="1_")
    print(f"GRF vertical peaks: right {Fr[:, 2].max():.0f} N, "
          f"left {Fl[:, 2].max():.0f} N")
    rows = np.array([model.jnt_dofadr[
        mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, j)]
        for j in FIT_JOINTS])
    from itertools import product
    FRONTAL = ("hip_adduction_r", "hip_adduction_l", "subtalar_angle_r",
               "subtalar_angle_l", "pelvis_list",
               "hip_rotation_r", "hip_rotation_l")
    best_cost, best_combo = np.inf, (1.0,) * len(FRONTAL)
    for combo in product((1.0, -1.0), repeat=len(FRONTAL)):
        s_t = dict(signs)
        s_t.update(zip(FRONTAL, combo))
        qd_t = frames_qpos(model, data, s_t, ik, ik_names, t_ik)
        dt_ = float(np.mean(np.diff(t_ik)))
        v_t, a_t = smooth_derivs(qd_t, dt_)
        cost = id_cost(model, data, s_t, ik, ik_names, qd_t, v_t, a_t,
                       rows, Fr, Pr, Tr, Fl, Pl, Tl, t_ik)
        if cost < best_cost:
            best_cost, best_combo = cost, combo
    signs.update(zip(FRONTAL, best_combo))
    print(f"frontal sign sweep: best {dict(zip(FRONTAL, best_combo))} "
          f"(frontal residual {best_cost:.1f} N*m)", flush=True)

    # ----------------------------------------------- full MuJoCo kinematics
    L_mj = np.zeros((T, model.nu))
    qd = frames_qpos(model, data, signs, ik, ik_names, t_ik)
    dt = float(np.mean(np.diff(t_ik)))
    qvel, qacc = smooth_derivs(qd, dt)
    for k in range(T):
        set_frame(model, data, signs, {n: ik[n][k] for n in ik_names})
        data.qvel[:] = 0.0
        data.act[:] = 0.0
        mujoco.mj_forward(model, data)
        L_mj[k] = data.actuator_length
    A_fd = fd_moments(model, data, qd, rows)

    # ------------------------------------------- inverse dynamics + NNLS
    cal_r = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_r")
    cal_l = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_l")
    a_sol = np.zeros((T, model.nu))
    tau_fit = np.zeros((T, len(rows)))
    res_norm = np.zeros(T)
    jacp = np.zeros((3, model.nv))
    jacr = np.zeros((3, model.nv))
    from scipy.optimize import lsq_linear
    for k in range(T):
        set_frame(model, data, signs, {n: ik[n][k] for n in ik_names})
        data.qvel[:] = qvel[k]
        data.qacc[:] = qacc[k]
        data.act[:] = 0.0
        mujoco.mj_inverse(model, data)
        tau = data.qfrc_inverse.copy()
        # measured GRF as external generalized forces (force at CoP + free
        # torque), subtracted from the required actuator force
        for (F, P, Tq, bid) in ((Fr, Pr, Tr, cal_r), (Fl, Pl, Tl, cal_l)):
            mujoco.mj_jac(model, data, jacp, jacr, P[k], bid)
            tau -= jacp.T @ F[k] + jacr.T @ Tq[k]
        tau_fit[k] = tau[rows]
        # muscle force model at this pose (quasi-static: v=0). The converted
        # actuators are dyntype=muscle: data.act IS the activation state and
        # drives actuator_force directly; data.ctrl is only the excitation
        # target consumed by the integration step (inert under mj_forward -
        # diag_force.py: ctrl=1 force 0, act=1 force -2655 N on soleus_r)
        data.qvel[:] = 0.0
        data.act[:] = 0.0
        mujoco.mj_forward(model, data)
        f0 = data.actuator_force.copy()          # passive (act=0)
        data.act[:] = 1.0
        mujoco.mj_forward(model, data)
        f1 = data.actuator_force.copy()          # active+passive
        unit = f1 - f0                            # Fmax*FL(l) per muscle
        B = A_fd[k] * unit[None, :]               # [nrows, nu], dl/dtheta
        rhs = tau[rows] - A_fd[k] @ f0
        # SO-style solve: min ||B a - rhs||^2 + lam^2 ||a||^2 with a in [0,1]
        # (plain NNLS dumps absurd activations into muscles whose FL ~ 0 at
        # these lengths - near-zero-norm columns soaking up residual; the
        # ridge is exactly OpenSim SO's sum-of-squares objective, and the
        # bounds are its activation limits)
        cn = np.linalg.norm(B, axis=0)
        good = cn > 1e-9
        lam = 0.05 * float(np.median(cn[good])) if good.any() else 1.0
        Aaug = np.vstack([B, lam * np.eye(model.nu)])
        baug = np.concatenate([rhs, np.zeros(model.nu)])
        sol = lsq_linear(Aaug, baug, bounds=(0.0, 1.0), tol=1e-10,
                         max_iter=200)
        a_sol[k] = np.where(good, sol.x, 0.0)
        res_norm[k] = np.linalg.norm(B @ a_sol[k] - rhs) / \
            max(np.linalg.norm(rhs), 1e-9)

    over = (a_sol > 1.001).sum()
    print(f"back-solve: mean residual {res_norm.mean():.3f}, entries a>1: "
          f"{over}, act range [{a_sol.min():.2f}, {a_sol.max():.2f}]")
    tq = np.abs(tau_fit)
    print("ID torque medians (N*m): " +
          ", ".join(f"{jn}={tq[:, i].mean():.1f}"
                    for i, jn in enumerate(FIT_JOINTS)
                    if jn.endswith(("_r",)) or jn.startswith("lumbar")))

    # ------------------------------------------------------------ gait events
    # MuJoCo frame: vertical = component 2 (mj z = OpenSim y)
    Fv_r = Fr[:, 2]
    Fv_l = Fl[:, 2]
    print(f"GRF vertical peaks: right {Fv_r.max():.0f} N, left {Fv_l.max():.0f} N")
    on_r = (Fv_r > 50.0).astype(int)
    on_l = (Fv_l > 50.0).astype(int)
    hs_r = t_ik[np.flatnonzero(np.diff(on_r) == 1)]
    duty_r = float(on_r.mean()) if on_r.any() else float("nan")
    duty_l = float(on_l.mean()) if on_l.any() else float("nan")
    print(f"gait: right-stance fraction {duty_r:.2f}, left {duty_l:.2f}, "
          f"{len(hs_r)} right HS onsets at {np.round(hs_r, 2)}")
    if len(hs_r) >= 2:
        cyc = np.diff(hs_r)
        print(f"cycle period {cyc.mean():.3f}+-{cyc.std():.3f} s "
              f"({1 / cyc.mean():.2f} Hz)")

    # phase profile per functional group over the right gait cycle
    grp_prof, grp_prof_l = {}, {}
    if len(hs_r) >= 2:
        t0, t1 = hs_r[0], hs_r[1]
        keep = (t_ik - t0 >= 0) & (t_ik - t0 < t1 - t0)
        ph = (t_ik[keep] - t0) / (t1 - t0)
        bins = np.clip((ph * 20).astype(int), 0, 19)
        for g in GROUPS:
            rm = [i for i, a in enumerate(acts)
                  if (mi := classify(a)) and mi.side == "r"
                  and g in mi.groups and Fmax[i] > 5]
            lm = [i for i, a in enumerate(acts)
                  if (mi := classify(a)) and mi.side == "l"
                  and g in mi.groups and Fmax[i] > 5]
            if rm:
                grp_prof[g] = np.array(
                    [a_sol[keep][bins == b][:, rm].mean() for b in range(20)])
            if lm:
                grp_prof_l[g] = np.array(
                    [a_sol[keep][bins == b][:, lm].mean() for b in range(20)])

    np.savez_compressed(
        HERE / "bsolve_out.npz",
        t=t_ik, qpos=qd, signs=json.dumps(signs),
        acts=a_sol, act_names=np.array(acts), Fmax=Fmax,
        tau_fit=tau_fit, fit_joints=np.array(FIT_JOINTS),
        res_norm=res_norm, L_mj=L_mj,
        L_os_names=L_os["names"] if L_os else np.array([]),
        L_os=L_os["data"] if L_os else np.zeros((0, 0)),
        so_names=acts_os["names"] if acts_os else np.array([]),
        so_act=acts_os["data"] if acts_os else np.zeros((0, 0)),
        grp_prof=(np.stack([grp_prof[g] for g in grp_prof])
                  if grp_prof else np.zeros((0, 0))),
        grp_prof_names=np.array(list(grp_prof)),
        grp_prof_l=(np.stack([grp_prof_l[g] for g in grp_prof_l])
                    if grp_prof_l else np.zeros((0, 0))),
        grp_prof_l_names=np.array(list(grp_prof_l)),
        duty_r=duty_r, duty_l=duty_l, hs_r=hs_r,
    )
    print("saved bsolve_out.npz")

    if L_os is not None:
        compare_and_plot(L_mj, L_os, a_sol, acts_os, acts, Fmax,
                         t_ik, grp_prof, grp_prof_l, A_fd, MA_os, signs)


def compare_and_plot(L_mj, L_os, a_sol, acts_os, acts, Fmax, t_ik,
                     grp_prof, grp_prof_l, A_fd=None, MA_os=None, signs=None):
    FIT_ROW_IDX = {j: i for i, j in enumerate(FIT_JOINTS)}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    T = min(L_mj.shape[0], L_os["data"].shape[0])
    col_os = {n: i for i, n in enumerate(L_os["names"])}
    pruned, unmatched = [], []
    report, rows = [], []
    for i, a in enumerate(acts):
        if Fmax[i] <= 5:                         # pruned (quad_fem/gem/peri)
            pruned.append(a)
            continue
        if a not in col_os:
            unmatched.append(a)
            continue
        mj, ref = L_mj[:T, i], L_os["data"][:T, col_os[a]]
        rng = max(np.ptp(ref), 1e-3)
        rmse = float(np.sqrt(np.mean((mj - ref) ** 2)))
        r = float(np.corrcoef(mj, ref)[0, 1]) if np.ptp(mj) > 1e-6 else 0.0
        rows.append((a, r, rmse * 100, rmse / rng))
    rows.sort(key=lambda x: x[3])
    report.append(f"coverage: {len(acts)} MJCF actuators, {len(pruned)} "
                  f"pruned ({', '.join(pruned)}), {len(unmatched)} with no "
                  f"name match in the OpenSim lengths file "
                  f"({len(L_os['names'])} columns): "
                  f"{', '.join(unmatched) if unmatched else 'none'}")
    report.append(f"{'muscle':14s} {'r':>6s} {'RMSE cm':>8s} {'nRMSE':>7s}")
    for a, r, cm, nr in rows:
        report.append(f"{a:14s} {r:6.3f} {cm:8.3f} {nr:7.3f}")
    rs = np.array([x[1] for x in rows])
    nrs = np.array([x[3] for x in rows])
    report.append(f"\nmedian r {np.median(rs):.3f} | median nRMSE "
                  f"{np.median(nrs):.3f} | muscles r>0.9: {(rs > 0.9).sum()}"
                  f"/{len(rs)} | nRMSE>0.2: {(nrs > 0.2).sum()}")

    col_so, S = {}, T
    if acts_os is not None and len(acts_os["names"]) > 1:
        col_so = {n: i for i, n in enumerate(acts_os["names"])}
        S = min(a_sol.shape[0], acts_os["data"].shape[0])
        rs2 = []
        for i, a in enumerate(acts):
            if Fmax[i] <= 5 or a not in col_so:
                continue
            x, y = a_sol[:S, i], acts_os["data"][:S, col_so[a]]
            if np.ptp(x) > 1e-6 and np.ptp(y) > 1e-6:
                rs2.append((a, float(np.corrcoef(x, y)[0, 1]),
                            float(np.sqrt(np.mean((x - y) ** 2)))))
        rs2.sort(key=lambda x: -x[1])
        r_all = np.array([x[1] for x in rs2])
        report.append(f"\nNNLS vs OpenSim SO: median r {np.median(r_all):.3f} "
                      f"| r>0.7: {(r_all > 0.7).sum()}/{len(r_all)}")
        report.append("best:  " + ", ".join(f"{a}({r:.2f})" for a, r, _ in rs2[:8]))
        report.append("worst: " + ", ".join(f"{a}({r:.2f})" for a, r, _ in rs2[-8:]))

    # moment-arm comparison (FD dl/dtheta vs OpenSim's dl/dtheta; the
    # expected transmission minus-sign is folded in)
    if MA_os and A_fd is not None and signs:
        Tm = min(A_fd.shape[0], min(v["data"].shape[0] for v in MA_os.values()))
        report.append("\nmoment arms vs OpenSim (median r; MuJoCo transmission "
                      "sign folded in):")
        for jn in ("knee_angle_r", "ankle_angle_r", "hip_flexion_r",
                   "hip_flexion_l"):
            if jn not in MA_os or jn not in FIT_ROW_IDX:
                continue
            col_m = {n: i for i, n in enumerate(MA_os[jn]["names"])}
            rs_j = []
            for i, a in enumerate(acts):
                if Fmax[i] <= 5 or a not in col_m:
                    continue
                mj_arm = A_fd[:Tm, FIT_ROW_IDX[jn], i]
                os_arm = -signs.get(jn, 1.0) * MA_os[jn]["data"][:Tm, col_m[a]]
                if np.ptp(mj_arm) > 1e-9 and np.ptp(os_arm) > 1e-9:
                    rs_j.append(np.corrcoef(mj_arm, os_arm)[0, 1])
            if rs_j:
                rs_j = np.array(rs_j)
                report.append(f"  {jn:16s} median r {np.median(rs_j):6.3f} "
                              f"(r>0.8: {(rs_j > 0.8).sum()}/{len(rs_j)}, "
                              f"sign-agree {(rs_j > 0).mean() * 100:.0f}%)")
        # vastii-specific: Ben's "no patella" check - the vastii must
        # actuate the tibia/knee through the 36 pathpoint couplers, so
        # their knee moment arms should match OpenSim's closely
        vastii = [a for a in acts if a.startswith(("vas_med_r", "vas_lat_r",
                                                   "vas_int_r"))]
        if "knee_angle_r" in MA_os:
            col_m = {n: i for i, n in enumerate(MA_os["knee_angle_r"]["names"])}
            ri = FIT_ROW_IDX["knee_angle_r"]
            vs = []
            for a in vastii:
                i = acts.index(a)
                if a not in col_m:
                    continue
                mj_arm = A_fd[:Tm, ri, i]
                os_arm = -signs.get("knee_angle_r", 1.0) * \
                    MA_os["knee_angle_r"]["data"][:Tm, col_m[a]]
                if np.ptp(mj_arm) > 1e-9 and np.ptp(os_arm) > 1e-9:
                    vs.append((a, float(np.corrcoef(mj_arm, os_arm)[0, 1]),
                               float(np.mean(mj_arm)),
                               float(np.mean(os_arm))))
            for a, r, mm, mo in vs:
                report.append(f"  vastii check {a:12s} r {r:6.3f}  "
                              f"mean FD arm {mm:+.4f} m vs OS {mo:+.4f} m")

    # key-muscle length trajectories
    key = ["soleus_r", "med_gas_r", "vas_lat_r", "rect_fem_r", "semimem_r",
           "bifemlh_r", "glut_max2_r", "iliacus_r", "psoas_r", "tib_ant_r",
           "per_long_r", "ercspn_r"]
    fig, axes = plt.subplots(3, 4, figsize=(15, 8), sharex=True)
    for ax, kn in zip(axes.ravel(), key):
        if kn in col_os and kn in acts:
            ax.plot(t_ik[:T], L_mj[:T, acts.index(kn)], label="MuJoCo")
            ax.plot(t_ik[:T], L_os["data"][:T, col_os[kn]], "--", label="OpenSim")
        ax.set_title(kn, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("muscle-tendon length along subject01 IK: MuJoCo vs OpenSim")
    fig.tight_layout()
    fig.savefig(HERE / "bsolve_lengths.png", dpi=130)
    plt.close(fig)

    # activations: NNLS vs SO + group profiles
    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    for ax, kn in zip(axes.ravel(),
                      ["soleus_r", "vas_lat_r", "semimem_r", "glut_max2_r",
                       "med_gas_r", "tib_ant_r", "iliacus_r", "rect_fem_r"]):
        if kn in col_so and kn in acts:
            ax.plot(t_ik[:S], a_sol[:S, acts.index(kn)], label="NNLS (MuJoCo)")
            ax.plot(t_ik[:S], acts_os["data"][:S, col_so[kn]], "--",
                    label="OpenSim SO")
        ax.set_title(kn, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("back-solved activations: per-frame NNLS vs OpenSim SO")
    fig.tight_layout()
    fig.savefig(HERE / "bsolve_activations.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for nm, prof in grp_prof.items():
        ax[0].plot(np.arange(20) * 5, prof, label=nm)
    ax[0].set_title("right-leg group activations vs gait cycle [%]")
    ax[0].set_xlabel("gait cycle [%] (0 = right foot strike)")
    ax[0].legend(fontsize=7, ncol=2)
    for nm, prof in grp_prof_l.items():
        ax[1].plot(np.arange(20) * 5, prof, label=nm)
    ax[1].set_title("left-leg group activations vs right-normalized cycle [%]")
    ax[1].legend(fontsize=7, ncol=2)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "bsolve_groups.png", dpi=130)
    plt.close(fig)

    (HERE / "bsolve_report.txt").write_text("\n".join(report),
                                            encoding="utf-8")
    print("\n".join(report))


if __name__ == "__main__":
    main(sys.argv[1:])
