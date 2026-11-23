"""Run a full physical-curriculum trial on the 54-route pneumatic robotbody."""
import argparse
import io
import json
import os
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
os.environ["CONDA_PREFIX"] = str(Path(sys.executable).parent)

import mujoco
import numpy as np
from numba import njit
from scipy.spatial import ConvexHull
from scipy.optimize import lsq_linear

from run_robotbody_bpa import HERE, SPINAL, JOINTS, CoupledArms, model_and_routes, transmission_limit
from robotbody_survey import configuration
from cad_interface import RoutedBPABank
from actuator_interface import force_at_pressure as sleeve_force, inverse_pressure
from torso_support import TorsoSupport, trunk_gate

sys.path[:0] = [str(SPINAL), str(SPINAL / "w2l_cpg"), str(SPINAL / "hunt_training")]
import build_network as BN
import params as P
from compact_pairwise_eval import CompactPairwiseForward
from curriculum_metrics import assess as original_motion_assessment
from ground_contact_check import check_samples
from regional_contacts import RegionalContacts


DURATION = {1: 16.0, 2: 16.0, 3: 8.0, 4: 16.0, 5: 16.0}
DEFAULT = dict(
    drive=3.15,
    tau=0.27,
    posture=0.18,
    ext=0.9,
    flex=1.2,
    hip=1.0,
    ankle=1.0,
    rg_to_pf=2.93,
    desc_e=1.59,
    desc_f=2.20,
    commissural=0.8,
    static_bias=0.8,
    stretch=1.0,
    speed=1.0,
    tension=1.0,
    heel=0.25,
    toe=0.1,
    stance_feedback=0.2,
    vestibular=0.2,
    kick=2.0,
    support4=0.65,
    support5=0.45,
)
NATIVE_PF = {phase: dict(values) for phase, values in P.W_PF_MN.items()}


@njit(cache=True)
def force_at_pressure(length, pressure, p):
    raw = sleeve_force(length, pressure, p)
    return np.minimum(raw, transmission_limit(length, p[0], p[1], p[2], p[4]))


@njit(cache=True)
def advance_pressures(target, pressure, length, shift_in, shift_out, gain, params, fill, vent, dt, force):
    for j in range(len(target)):
        tau = fill[j] if target[j] > pressure[j] else vent[j]
        pressure[j] = target[j] + (pressure[j] - target[j]) * np.exp(-dt / tau)
        path = shift_in[j] + (length[j] - shift_out[j]) / gain[j]
        bundle = force_at_pressure(path, pressure[j], params[j])
        force[j] = bundle / gain[j]


def foot_geometry(m):
    groups = []
    for side in ("r", "l"):
        shapes = []
        for tag in ("talus", "calcn", "toes"):
            gid = m.geom(f"{tag}_{side}_geom_1").id
            mid = m.geom_dataid[gid]
            start = m.mesh_vertadr[mid]
            vertices = m.mesh_vert[start : start + m.mesh_vertnum[mid]]
            shapes.append((gid, vertices[ConvexHull(vertices).vertices].copy()))
        groups.append(shapes)
    return groups


def foot_height(d, shapes):
    return np.array(
        [
            min(float((v @ d.geom_xmat[g].reshape(3, 3)[2, :]).min() + d.geom_xpos[g, 2]) for g, v in side)
            for side in shapes
        ]
    )


def model_setup(stage, dt, pose_offset=0.0, guide_scale=1.0):
    """Mobile pelvis in all stages. Suspension uses a real torso-applied wrench."""
    _, names, routes, xml = model_and_routes(dt, "ground")
    root = ET.fromstring(xml)
    if stage <= 2:
        root.remove(root.find("contact"))
    m = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding="unicode"))
    d = mujoco.MjData(m)
    mujoco.mj_resetDataKeyframe(m, d, 0)
    arms = CoupledArms(m)
    shapes = foot_geometry(m)
    for side in ("r", "l"):
        j = m.joint("knee_angle_" + side)
        m.jnt_range[j.id] = np.radians([-120.0, -1.0])
        m.jnt_margin[j.id] = np.radians(2.0)
        m.jnt_solref[j.id] = [0.004, 1.0]
        m.jnt_solimp[j.id] = [0.999, 0.9999, 0.001, 0.5, 2.0]
        for name, value in (
            ("hip_flexion", 6.0),
            ("hip_adduction", 0.0),
            ("hip_rotation", 0.0),
            ("knee_angle", -15.0),
            ("ankle_angle", 5.0),
            ("subtalar_angle", 0.0),
            ("mtp_angle", 0.0),
        ):
            j = m.joint(name + "_" + side)
            d.qpos[j.qposadr[0]] = np.radians(
                value + pose_offset * (1.0 if side == "r" else -1.0) if name == "knee_angle" else value
            )
            if name in ("ankle_angle", "subtalar_angle", "mtp_angle", "hip_adduction"):
                m.jnt_stiffness[j.id] += 30.0 if name == "hip_adduction" else 10.0
                m.dof_damping[j.dofadr[0]] += 3.0
        for group, k, c in (("hip_adduction", 450.0, 24.0), ("hip_rotation", 400.0, 24.0)):
            j = m.joint(group + "_" + side)
            m.jnt_stiffness[j.id] += guide_scale * k
            m.dof_damping[j.dofadr[0]] += guide_scale * c
    if stage >= 3:
        for name in ("lumbar_extension", "lumbar_bending", "lumbar_rotation"):
            j = m.joint(name)
            m.jnt_stiffness[j.id] += guide_scale * 300.0
            m.dof_damping[j.dofadr[0]] += guide_scale * 25.0
    arms.project(d)
    mujoco.mj_forward(m, d)
    # The floor is world z=0 in this specific converted model.
    plane = m.geom("ground-plane").id
    if abs(d.geom_xpos[plane, 2]) > 1e-7 or d.geom_xmat[plane].reshape(3, 3)[2, 2] < 0.999:
        raise ValueError("The floor is not the expected z-up plane.")
    height = 0.22 if stage <= 2 else -0.001
    d.qpos[m.joint("pelvis_ty").qposadr[0]] += (height - foot_height(d, shapes).min()) / d.xaxis[
        m.joint("pelvis_ty").id, 2
    ]
    mujoco.mj_forward(m, d)
    return m, d, names, routes, arms, shapes


def calibration_file(path):
    """Estimate extension excursion structurally on robotbody paths before neural tuning."""
    from scipy.stats import qmc

    model, d, names, routes, arms, _ = model_setup(3, 0.0005)
    q0 = d.qpos.copy()
    base = d.actuator_length.copy()
    # This is a finite calibration atlas. Trials still enforce the force-law
    # domain online and stop if motion exceeds what the installation admits.
    bounds = {
        "hip_flexion": (-40.0, 75.0),
        "hip_adduction": (-25.0, 25.0),
        "hip_rotation": (-25.0, 25.0),
        "knee_angle": (-120.0, 0.0),
        "ankle_angle": (-48.0, 48.0),
        "subtalar_angle": (-38.0, 38.0),
        "mtp_angle": (-40.0, 80.0),
        "lumbar_extension": (-10.0, 10.0),
        "lumbar_bending": (-10.0, 10.0),
        "lumbar_rotation": (-10.0, 10.0),
    }
    limits = np.array([bounds[name[:-2] if name.endswith(("_r", "_l")) else name] for name in JOINTS])
    points = qmc.Sobol(len(JOINTS), scramble=True, seed=617).random_base2(12)
    # Include all sign corners for each leg and trunk, since a multijoint route
    # can have its greatest excursion at several simultaneous extremes.
    for indices in (list(range(7)), list(range(7, 14)), [14, 15, 16]):
        corners = np.full((2 ** len(indices), len(JOINTS)), 0.5)
        corners[:, indices] = [[int(c) for c in f"{i:0{len(indices)}b}"] for i in range(len(corners))]
        points = np.concatenate((points, corners))
    max_seen = base.copy()
    min_seen = base.copy()
    for p in points:
        d.qpos[:] = q0
        d.qpos[arms.addresses] = np.radians(limits[:, 0] + p * np.diff(limits).ravel())
        arms.project(d)
        mujoco.mj_kinematics(model, d)
        mujoco.mj_tendon(model, d)
        lengths = d.ten_length[model.actuator_trnid[:, 0]]
        max_seen = np.maximum(max_seen, lengths)
        min_seen = np.minimum(min_seen, lengths)
    output = dict(
        samples=len(points),
        bounds_deg=bounds,
        source="converted Gait2392 robotbody reference pose",
        names=names,
        neutral_m=base.tolist(),
        minimum_m=min_seen.tolist(),
        maximum_m=max_seen.tolist(),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2))
    return output


def installation(names, routes, variant, atlas, reserve, trunk_contraction=None):
    # The original transfer fixture is referenced to its minus-30-degree
    # geometry. A saved atlas uses a distinct, documented operating-point datum.
    if atlas is None:
        ref, _, _, _ = model_and_routes(0.0005, "ground")
        start = mujoco.MjData(ref)
        mujoco.mj_resetDataKeyframe(ref, start, 0)
        mujoco.mj_forward(ref, start)
        return configuration(
            names,
            start.actuator_length.copy(),
            routes,
            variant,
            "transfer",
            0.4,
            {"tib_ant": 0.7, "per_tert": 0.7},
        )
    if names != atlas["names"]:
        raise ValueError("Calibration route order differs from this model.")
    pose = np.array(atlas["neutral_m"])
    longest = np.array(atlas["maximum_m"])
    template = configuration(names, pose, routes, variant, "transfer", 0.4)
    fractions = {}
    for i, name in enumerate(names):
        r = template["actuators"][name]
        family = name[:-2]
        if family in fractions:
            continue
        ids = [names.index(family + "_" + s) for s in ("r", "l")]
        out_travel = max(longest[ids] - pose[ids])
        stroke = r["resting_length_m"] - r["free_length_at_620kPa_m"]
        fractions[family] = max(0.12, (out_travel + reserve) / (r["transmission_gain"] * stroke))
    if trunk_contraction is not None:
        if not 0.0 <= trunk_contraction <= 1.0:
            raise ValueError("The specified trunk mounting datum is a fraction of the documented stroke.")
        for family in ("ercspn", "intobl", "extobl"):
            fractions[family] = trunk_contraction
    result = configuration(names, pose, routes, variant, "transfer", 0.4, fractions)
    result["mount"] = "locomotion atlas"
    result["extra_travel_reserve_m"] = reserve
    if trunk_contraction is not None:
        result["mount"] = "locomotion atlas with separately positioned trunk sleeves"
        result["trunk_reference_contraction"] = trunk_contraction
    return result


class NeuralControl:
    def __init__(self, names, stage, p):
        self.stage, self.p = stage, p
        P.TAU["rg_nap_h"] = p["tau"]
        P.G.update(
            rg_to_pf=p["rg_to_pf"],
            descend_to_rg_e=p["desc_e"],
            descend_to_rg_f=p["desc_f"],
            c1_gain=p["commissural"],
            heel_rge=p["heel"] if stage >= 4 else 0.0,
            toe_rge=p["toe"] if stage == 5 else 0.0,
            ib_rge=p["stance_feedback"] if stage >= 4 else 0.0,
            vest_ext=p["vestibular"] if stage >= 3 else 0.0,
            vest_flex_inh=0.5 * p["vestibular"] if stage >= 3 else 0.0,
        )
        layer_contacts = p.get("layer_contacts", False) and stage == 5
        P.G.update(
            joint_pf=float(layer_contacts),
            heel_pf_layer=p.get("heel_pf", 0.25) if layer_contacts else 0.0,
            toe_df_inh=p.get("toe_pf", 0.25) if layer_contacts else 0.0,
            heel_in_f_exc=p.get("heel_f", 0.15) if layer_contacts else 0.0,
        )
        for phase, row in NATIVE_PF.items():
            for group, w in row.items():
                gain = (
                    p["hip"]
                    if group.startswith("hip_")
                    else p["ankle"]
                    if group.startswith("ankle_")
                    else 1.0
                )
                if group == "knee_ext":
                    gain *= p["ext"]
                if group == "knee_flex":
                    gain *= p["flex"]
                if group in ("hip_flex", "hip_ext", "trunk_flex", "trunk_ext", "ankle_pf", "ankle_df"):
                    gain *= p.get(group + "_scale", 1.0)
                P.W_PF_MN[phase][group] = w * gain
        self.n = BN.SpinalNetwork(
            muscles={s: BN.classify(s) for s in names}, sides=("l", "r"), interleg=stage != 1
        )
        stance_gain = p.get("load_to_vastus_uS", 0.0)
        if stage >= 4 and stance_gain > 0.0:
            # Test whether low vastus pressure during loaded stance reflects
            # recruitment. The existing ipsilateral load interneuron excites
            # the vastus motor pool; valve and muscle dynamics still intervene.
            for side in ("r", "l"):
                self.n.net.add_connection(
                    BN._syn(stance_gain, exc=True), "LBIN_" + side, "MN_vas_med_" + side
                )
        for cell in ("RG_E_r", "RG_F_l"):
            self.n.net.add_input(cell)
            self.n.inputs.append("start_" + cell)
        self.n.compile(0.002)
        self.n.compiled = CompactPairwiseForward(self.n.compiled)
        self.u = self.n.make_inputs()
        self.ports = {name: i for i, name in enumerate(self.n.inputs)}
        self.mn = np.array([self.n.idx["MN_" + n] for n in names])
        self.rg = np.array([self.n.idx["RG_F_" + s] for s in ("r", "l")])
        self.sens = {g: np.array([self.ports[g + "_" + n] for n in names]) for g in ("Ia", "II", "Ib")}
        self.bias_ports = np.array([self.ports["POST_" + n] for n in names])
        self.start_ports = [self.ports["start_RG_E_r"], self.ports["start_RG_F_l"]]
        self.range = P.E_HI

    def update(self, t, initial, length, speed, force, scale, grf, bw, tilt, balance=None):
        p, u, ports = self.p, self.u, self.ports
        u[ports["POSTURE"]] = p["posture"]
        u[ports["DRIVE"]] = 0.0 if self.stage == 3 or t < 0.0 else p["drive"]
        u[self.start_ports] = p["kick"] if self.stage != 3 and 0.0 <= t < 0.18 else 0.0
        if self.stage != 1:
            u[self.sens["Ia"]] = p["speed"] * np.maximum(speed, 0.0)
            u[self.sens["II"]] = (
                0.3 * p["stretch"] * np.maximum(length - initial, 0.0) / np.maximum(0.1 * initial, 0.01)
            )
            u[self.sens["Ib"]] = 0.15 * p["tension"] * force / scale
        for i, s in enumerate(("r", "l")):
            if "VEST_c_" + s in ports:
                u[ports["VEST_c_" + s]] = (
                    balance[3] if balance is not None else max(0.0, (1 if i == 0 else -1) * tilt)
                )
            if self.stage >= 4:
                u[ports["HEEL_c_" + s]] = min(grf[i, 0] / (0.35 * bw), 1.5)
                u[ports["TOE_c_" + s]] = min(grf[i, 1] / (0.5 * bw), 1.5) if self.stage == 5 else 0.0
                u[ports["LOAD_c_" + s]] = min(sum(grf[i]) / (0.6 * bw), 1.5)
        if balance is not None:
            fore, lateral, trunk, _ = balance
            for name, current in (
                ("BAL_PF", fore),
                ("BAL_DF", -fore),
                ("BAL_TRK_FLX", trunk),
                ("BAL_TRK_EXT", -trunk),
            ):
                u[ports[name]] = max(0.0, current)
            for i, name in enumerate(("BAL_LAT_R", "BAL_LAT_L")):
                u[ports[name]] = max(0.0, -lateral if i == 0 else lateral) * min(
                    1.0, grf[i].sum() / (0.15 * bw)
                )
        v = self.n.step(u)
        return 620.0 * np.clip(v[self.mn] / self.range, 0.0, 1.0), v[self.rg].copy()


class PosturalSensors:
    """Drive the existing spinal balance neurons from measured COM and torso motion."""

    def __init__(self, m, d, pelvis, gain):
        self.pelvis = pelvis
        self.torso = m.body("torso").id
        self.gain = gain
        self.reference = self.read(d)
        self.previous = self.reference.copy()

    def read(self, d):
        up = d.xmat[self.torso].reshape(3, 3)[:, 1]
        return np.array([*d.subtree_com[self.pelvis, :2], np.arctan2(-up[0], max(up[2], 1e-6))])

    def sample(self, d):
        measured = self.read(d)
        delta = measured - self.reference
        velocity = (measured - self.previous) / 0.002
        self.previous = measured
        b = P.BAL
        fore = np.clip(
            self.gain * (-b["kx"] * delta[0] - b["kv"] * velocity[0]), -b["max_current"], b["max_current"]
        )
        lateral = np.clip(
            self.gain * (-b["ky_lat"] * delta[1] - b["kv_lat"] * velocity[1]),
            -b["max_current"],
            b["max_current"],
        )
        trunk = np.clip(
            self.gain * (b["kp_trk"] * delta[2] + b["kd_trk"] * velocity[2]), -b["max_trk"], b["max_trk"]
        )
        tone = np.clip(
            self.gain * (b["kp_trk"] * abs(delta[2]) + b["kd_trk"] * abs(velocity[2])), 0.0, b["max_trk"]
        )
        return fore, lateral, trunk, tone


def lift_fraction(stage, p):
    if stage <= 2:
        return 1.0
    return p.get(
        "ground_lift_override", 0.35 if stage == 3 else p["support4"] if stage == 4 else p["support5"]
    )


def posture_inputs(m, d, arms, names, bank, gain, shift_in, shift_out, stage, p, bw):
    """Compute one static, bilateral motor-current preset, then freeze it."""
    imu = m.body("torso" if stage <= 2 else "pelvis").id
    external = np.zeros(m.nv)
    actual_fraction = lift_fraction(stage, p)
    mujoco.mj_applyFT(
        m, d, np.array([0.0, 0.0, actual_fraction * bw]), np.zeros(3), d.xipos[imu], imu, external
    )
    if stage >= 3:
        for s in ("r", "l"):
            heel = m.body("calcn_" + s).id
            toe = m.body("toes_" + s).id
            center = 0.4 * d.xpos[heel] + 0.6 * d.xpos[toe]
            mujoco.mj_applyFT(
                m,
                d,
                np.array([0.0, 0.0, (1 - actual_fraction) * 0.5 * bw]),
                np.zeros(3),
                center,
                heel,
                external,
            )
    path = shift_in + (d.actuator_length - shift_out) / gain
    capacity = np.array(
        [force_at_pressure(l, 620.0, row) / g for l, row, g in zip(path, bank.parameters, gain)]
    )
    if not np.isfinite(capacity).all():
        return np.zeros(m.nu), dict(reason="initial pose exceeds static capacity domain")
    A = arms.update(d).copy() * capacity
    demand = arms.E.T @ (d.qfrc_bias - d.qfrc_passive - external)
    families = sorted({s[:-2] for s in names})
    groups = np.array([families.index(s[:-2]) for s in names])
    paired = np.column_stack([A[:, groups == j].sum(axis=1) for j in range(len(families))])
    matrix = np.vstack((paired, np.eye(len(families)) * 0.5))
    fit = lsq_linear(matrix, np.r_[demand, np.zeros(len(families))], bounds=(0.0, 1.0))
    fraction = fit.x[groups]
    currents = np.zeros(m.nu)
    for j in range(len(names)):
        pres, flag = inverse_pressure(fraction[j] * capacity[j] * gain[j], path[j], bank.parameters[j])
        if flag in (0, 1):
            currents[j] = 5.0 * pres / 620.0
    return currents, dict(
        rms_equilibrium_residual_Nm=float(np.linalg.norm(paired @ fit.x - demand) / np.sqrt(len(demand))),
        precomputed_once=True,
        prescribed_foot_reactions_applied_in_forward_dynamics=False,
        current_nA=currents.tolist(),
    )


def standing_statistics(q, com, load, pelvis_height, bw):
    start = min(len(com) - 1, 2000)
    segment = com[start:]
    forces = load[start:].mean(axis=0)
    return dict(
        bal_tilt_max=float(np.max(abs(q[start:]))),
        bal_sway=float(np.linalg.norm(segment[:, :2] - segment[:1, :2], axis=1).max()),
        bal_foot_load=float(sum(forces) / bw),
        bal_contact_sym=float(forces[0] / max(sum(forces), 1e-12)),
        bal_fell=bool(np.min(pelvis_height) < 0.55),
    )


def evaluate(variant, stage, p, dt, atlas=None, reserve=0.01, pose_offset=0.0, save=None, duration=None):
    start = time.perf_counter()
    duration = DURATION[stage] if duration is None else duration
    ratio = round(0.002 / dt)
    if not np.isclose(ratio * dt, 0.002) or ratio < 1:
        raise ValueError("Physics step must divide two milliseconds.")
    guides = p.get("guidance_scale", 1.0)
    if not 0.0 <= guides <= 1.0 or not 0.0 <= lift_fraction(stage, p) <= 1.0:
        raise ValueError("Ground guide scale and lift must be fractions between zero and one.")
    m, d, names, routes, arms, shapes = model_setup(stage, dt, pose_offset, guides)
    if "contact_time_s" in p:
        if not 0.004 <= p["contact_time_s"] <= 0.04:
            raise ValueError(
                "For the contact sensitivity test, use a positive 4–40 ms constraint time constant."
            )
        m.pair_solref[:, 0] = p["contact_time_s"]
        mujoco.mj_forward(m, d)
    initial = d.actuator_length.copy()
    cfg = installation(names, routes, variant, atlas, reserve, p.get("trunk_reference_contraction"))
    bank = RoutedBPABank(names, cfg, dt)
    rows = [cfg["actuators"][s] for s in names]
    gain = np.array([r["transmission_gain"] for r in rows])
    shift_in = np.array([r["bpa_side_reference_m"] for r in rows])
    shift_out = np.array([r["output_reference_m"] for r in rows])
    force = np.zeros(m.nu)
    speed = np.zeros(m.nu)
    pmax = bank.parameters[:, 3] / gain
    target = bank.pressure.copy()
    net = NeuralControl(names, stage, p)
    regional = RegionalContacts(m)
    pelvis = m.body("pelvis").id
    BW = float(m.body_subtreemass[pelvis] * 9.81)
    air = TorsoSupport(m, d, stage=stage) if stage <= 2 else None
    reference = d.xipos[pelvis].copy()
    refrot = d.xmat[pelvis].reshape(3, 3).copy()
    balance = (
        PosturalSensors(m, d, pelvis, p["postural_feedback_gain"])
        if stage >= 3 and p.get("postural_feedback_gain", 0.0) > 0.0
        else None
    )
    jj = {n: m.joint(n).qposadr[0] for n in ["pelvis_list", "pelvis_tilt", "pelvis_rotation"]}
    vel = np.zeros(6)
    ground = np.zeros((2, 2))
    wrench = np.zeros(6)
    stride = int(round(duration / 0.002))
    i = -1
    rg = np.zeros(2)
    reason = "duration"
    bad = []
    trace = {
        key: np.zeros((stride, ncol))
        for key, ncol in [
            ("q", len(JOINTS)),
            ("grf", 2),
            ("support", 6),
            ("width", 2),
            ("kinematics", 5),
            ("clearance", 2),
            ("rg", 2),
        ]
    }
    if save:
        trace.update(
            {
                key: np.zeros((stride, m.nu), dtype=np.float32)
                for key in ("pressure", "command", "tension", "length")
            }
        )
    peak = np.zeros(m.nu)
    maxcmd = np.zeros(m.nu)
    residual = 0.0
    maxpenetration = 0.0
    max_knee = np.full(2, -1e8)
    knees = [m.joint("knee_angle_" + s).qposadr[0] for s in ("r", "l")]
    feet = [m.body("calcn_" + s).id for s in ("r", "l")]
    tib = [m.body("tibia_" + s).id for s in ("r", "l")]
    frozen = p.get("frozen_post_currents_by_route_nA")
    adjustment = p.get("calibrate_initial_motor_currents", False) and frozen is None
    if frozen is None:
        preset, standing_preset = posture_inputs(
            m, d, arms, names, bank, gain, shift_in, shift_out, stage, p, BW
        )
        net.u[net.bias_ports] = p["static_bias"] * preset
    else:
        if set(frozen) != set(names) or not np.isfinite(list(frozen.values())).all():
            raise ValueError("Frozen motor preset requires one finite current for every named route.")
        net.u[net.bias_ports] = np.array([frozen[name] for name in names])
        standing_preset = dict(source="frozen per-route motor currents", precomputed_once=False)
    if adjustment:
        # At a held posture fit a constant current at each normal POST input
        # to the motor voltage needed for its static pressure. Freeze these
        # inputs before release. Neither the pose nor individual motor
        # voltages are imposed once the forward trial begins.
        desired = np.minimum(5.0, p["static_bias"] * preset)
        if "trunk_voltage_mV" in p:
            for j, name in enumerate(names):
                if name[:-2] in p["trunk_voltage_mV"]:
                    desired[j] = p["trunk_voltage_mV"][name[:-2]]
            if not np.isfinite(desired).all() or not ((0.0 <= desired) & (desired <= 5.0)).all():
                raise ValueError("Initial trunk motor targets must be in the 0–5 mV operating range.")
        for k in range(500):
            target, rg = net.update(-1.0 + k * 0.002, initial, initial, speed, force, pmax, ground, BW, 0.0)
            if k >= 50 and k % 5 == 0:
                volt = net.n.compiled.V[net.mn]
                error = desired - np.maximum(volt, 0.0)
                net.u[net.bias_ports] = np.clip(net.u[net.bias_ports] + 0.14 * error, -20.0, 15.0)
            force.fill(0.0)
            advance_pressures(
                target,
                bank.pressure,
                initial,
                shift_in,
                shift_out,
                gain,
                bank.parameters,
                bank.fill,
                bank.vent,
                0.002,
                force,
            )
        standing_preset["held_pose_target_voltage_mV"] = desired.tolist()
        standing_preset["mean_abs_preload_voltage_error_mV"] = float(
            np.mean(abs(desired - np.clip(net.n.compiled.V[net.mn], 0, 5)))
        )
    standing_preset["frozen_motor_current_nA"] = net.u[net.bias_ports].tolist()
    for h in range(round(0.6 / dt)):
        if h % ratio == 0:
            target, rg = net.update(h * dt - 0.6, initial, initial, speed, force, pmax, ground, BW, 0.0)
        advance_pressures(
            target,
            bank.pressure,
            initial,
            shift_in,
            shift_out,
            gain,
            bank.parameters,
            bank.fill,
            bank.vent,
            dt,
            force,
        )
    for h in range(stride * ratio):
        t = h * dt
        mujoco.mj_step1(m, d)
        if h % ratio == 0:
            target, rg = net.update(
                t,
                initial,
                d.actuator_length,
                d.actuator_velocity,
                force,
                pmax,
                ground,
                BW,
                d.qpos[jj["pelvis_tilt"]],
                balance.sample(d) if balance is not None else None,
            )
        advance_pressures(
            target,
            bank.pressure,
            d.actuator_length,
            shift_in,
            shift_out,
            gain,
            bank.parameters,
            bank.fill,
            bank.vent,
            dt,
            force,
        )
        if not np.isfinite(force).all():
            reason = "sleeve extension outside force domain"
            bad = [
                dict(
                    name=names[j],
                    route_m=float(d.actuator_length[j]),
                    pressure_kPa=float(bank.pressure[j]),
                    additional_unloaded_sleeve_length_needed_m=float(
                        shift_in[j]
                        + (d.actuator_length[j] - shift_out[j]) / gain[j]
                        - bank.parameters[j, 0]
                        - bank.parameters[j, 2]
                    ),
                )
                for j in np.flatnonzero(~np.isfinite(force))
            ]
            break
        d.ctrl[:] = -force
        if air:
            air.apply(m, d)
            wrench[:] = air.wrench
        else:
            mujoco.mj_objectVelocity(m, d, mujoco.mjtObj.mjOBJ_BODY, pelvis, vel, 0)
            k_lateral = 2000.0
            k_forward = 1500.0
            off = reference - d.xipos[pelvis]
            wrench[:2] = np.array([k_forward, k_lateral]) * off[:2] - np.array([220.0, 260.0]) * vel[3:5]
            wrench[:2] *= guides
            wrench[2] = BW * lift_fraction(stage, p)
            rot = d.xmat[pelvis].reshape(3, 3)
            wrench[3:] = guides * (400.0 * (0.5 * np.cross(rot.T, refrot.T).sum(axis=0)) - 40.0 * vel[:3])
            d.xfrc_applied[pelvis] = wrench
            if stage == 3 and 4.0 <= t < 4.15:
                d.xfrc_applied[pelvis, 0] += 35.0
        expected = -d.actuator_moment.T @ force
        mujoco.mj_step2(m, d)
        residual = max(residual, float(np.max(abs(d.qfrc_actuator - expected))))
        if not np.isfinite(d.qpos).all() or abs(d.time - (h + 1) * dt) > 1e-8:
            reason = "integration failure"
            break
        max_knee = np.maximum(max_knee, np.degrees(d.qpos[knees]))
        peak = np.maximum(peak, force)
        maxcmd = np.maximum(maxcmd, target)
        if d.ncon:
            maxpenetration = max(maxpenetration, -float(min(c.dist for c in d.contact)))
        if (h + 1) % ratio:
            continue
        ground = regional.sample(d).copy()
        mujoco.mj_kinematics(m, d)
        mujoco.mj_comPos(m, d)
        i += 1
        if air:
            air.record(d, force[air.ids])
        trace["q"][i] = np.degrees(d.qpos[arms.addresses])
        trace["support"][i] = wrench
        trace["grf"][i] = ground.sum(axis=1)
        trace["rg"][i] = rg
        trace["kinematics"][i] = [
            *d.subtree_com[pelvis],
            d.xpos[pelvis, 2],
            np.degrees(max(abs(d.qpos[jj["pelvis_tilt"]]), abs(d.qpos[jj["pelvis_list"]]))),
        ]
        direction = d.xmat[pelvis].reshape(3, 3)[:, 2]
        trace["width"][i] = [
            (d.xpos[feet[0]] - d.xpos[feet[1]]) @ direction,
            (d.xpos[tib[0]] - d.xpos[tib[1]]) @ direction,
        ]
        if i % 10 == 0:
            trace["clearance"][i] = foot_height(d, shapes)
        if save:
            trace["pressure"][i] = bank.pressure
            trace["command"][i] = target
            trace["tension"][i] = force
            trace["length"][i] = d.actuator_length
        if stage >= 3 and d.xpos[pelvis, 2] < 0.43:
            reason = "fell below 0.43 m pelvis height"
            break
    count = i + 1
    trace = {k: v[:count] for k, v in trace.items()}
    times = 0.002 * np.arange(1, count + 1)
    trace["t"] = times
    complete = count == stride and reason == "duration"
    out = dict(
        stage=stage,
        design=variant,
        duration_s=duration,
        dt_s=dt,
        completed_s=count * 0.002,
        complete=complete,
        reason=reason,
        failed_routes=bad,
        weight_N=BW,
        neurons=len(net.n.idx),
        commanded_routes=int(np.count_nonzero(maxcmd > 0)),
        tension_routes=int(np.count_nonzero(peak > 1e-8)),
        max_force_balance_error_Nm=residual,
        max_penetration_m=maxpenetration,
        max_knee_deg=max_knee.tolist(),
        physical_sleeves=int(sum(r["parallel_count"] for r in rows)),
        contact_time_constants_s=m.pair_solref[:, 0].tolist(),
        guide_scale=guides,
        nominal_lift_BW=lift_fraction(stage, p),
        measured_support_BW=None,
        scored_five_stage_gate=False,
        accepted=False,
        score=-1000.0 + count * 0.002,
        mount=cfg["mount"],
        parameters=p,
        warmup_fixed_pose_s=0.6 + (1.0 if adjustment else 0.0),
        static_neural_bias=standing_preset,
    )
    if count:
        window = (times >= 5.0) & (times < 13.0)
        if not window.any():
            window = times >= max(0.0, times[-1] - 0.5)
        out.update(
            measured_support_BW=float(trace["support"][window, 2].mean() / BW),
            mean_foot_load_BW=(trace["grf"][window].mean(axis=0) / BW).tolist(),
            horizontal_assistance_rms_N=float(
                np.sqrt(np.mean(np.sum(trace["support"][window, :2] ** 2, axis=1)))
            ),
            orientation_assistance_rms_Nm=float(
                np.sqrt(np.mean(np.sum(trace["support"][window, 3:] ** 2, axis=1)))
            ),
            knee_excursions_deg=np.ptp(trace["q"][window][:, [3, 10]], axis=0).tolist(),
            minimum_pelvis_m=float(trace["kinematics"][:, 3].min()),
        )
        if stage >= 3:
            stance_angles = []
            for side, joint in enumerate((3, 10)):
                loaded = window & (trace["grf"][:, side] > 0.05 * BW)
                values = trace["q"][loaded, joint]
                stance_angles.append(
                    [float(values.min()), float(np.median(values)), float(values.max())]
                    if len(values)
                    else None
                )
            out["loaded_knee_min_median_max_deg_RL"] = stance_angles
    if complete and duration == DURATION[stage]:
        kine = trace["kinematics"]
        stats = standing_statistics(kine[:, 4], kine[:, :3], trace["grf"], kine[:, 3], BW)
        metrics = dict(nan=False, t_end=duration, **stats)
        stream = io.BytesIO()
        np.savez(stream, t=times, q=trace["q"], key_joints=np.array(JOINTS))
        stream.seek(0)
        measure = original_motion_assessment(stream, metrics, 3 if stage == 3 else 2)
        frontal = bool((trace["width"][window].min(axis=0) > 0.02).all())
        ok = measure["accepted"] and frontal and max(max_knee) <= 0.0 and residual < 1e-8
        score = measure["score"] - max(0.0, max(max_knee)) * 2.0 - 30.0 * (not frontal)
        out.update(scored_five_stage_gate=True, movement=measure, frontal_gate=frontal, standing=stats)
        if air:
            out["torso"] = air.summary()
            ok = ok and trunk_gate(out["torso"])
            score -= min(50.0, max(0.0, out["torso"]["max_lumbar_excursion_deg"] - 15.0) * 0.5)
        if stage in (4, 5):
            samples = np.zeros((count, 9))
            samples[:, 0] = times
            samples[:, 1] = trace["support"][:, 2]
            samples[:, 7:9] = trace["grf"]
            contact = check_samples(samples, BW)
            out["contact"] = contact
            idx = np.flatnonzero(window & (np.arange(count) % 10 == 0))
            height = trace["clearance"][idx]
            forces = trace["grf"][idx] / BW
            swing = []
            for side in (0, 1):
                mask = (height[:, side] > 0.001) & (forces[:, side] < 0.025) & (forces[:, 1 - side] > 0.05)
                bounds = np.r_[0, np.flatnonzero(np.diff(mask.astype(int))) + 1, len(mask)]
                swing.append(
                    [
                        dict(
                            on_s=float(times[idx[a]]),
                            duration_s=float((b - a) * 0.020),
                            peak_m=float(max(height[a:b, side])),
                        )
                        for a, b in zip(bounds[:-1], bounds[1:])
                        if mask[a] and (b - a) >= 4 and max(height[a:b, side]) >= 0.010
                    ]
                )
            ok = ok and contact["accepted"] and min(map(len, swing)) >= 3 and stats["bal_tilt_max"] < 30.0
            out["loaded_clear_swings"] = swing
            score += (
                12.0 * min(3, min(map(len, swing)))
                + 6.0 * min(3, min(contact["cycles"]))
                + 10.0 * min(sum(contact["mean_foot_bw"]), 1.0)
            )
        out.update(accepted=bool(ok), score=float(score))
    if save:
        save.parent.mkdir(exist_ok=True, parents=True)
        np.savez_compressed(str(save) + ".npz", **trace, names=names, joint_names=JOINTS)
        Path(str(save) + ".installation.json").write_text(json.dumps(cfg, indent=2))
    out["wall_s"] = time.perf_counter() - start
    return out


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--calibrate", type=Path)
    cli.add_argument("--atlas", type=Path)
    cli.add_argument("--stage", type=int, choices=DURATION, default=1)
    cli.add_argument("--design", choices=("baseline", "selected"), default="baseline")
    cli.add_argument("--config", type=Path)
    cli.add_argument(
        "--fit",
        choices=("coarse", "refined", "local"),
        default="refined",
        help="Search round to read from a full study record.",
    )
    cli.add_argument(
        "--rank",
        type=int,
        choices=(0, 1),
        default=0,
        help="Finalist index if --config is a full study record.",
    )
    cli.add_argument("--output", type=Path)
    cli.add_argument("--dt", type=float, default=0.001)
    cli.add_argument("--reserve", type=float, default=0.01)
    cli.add_argument("--duration", type=float)
    cli.add_argument("--pose-offset", type=float, default=0.0)
    args = cli.parse_args()
    if args.calibrate:
        atlas = calibration_file(args.calibrate)
        print("Calibration poses:", atlas["samples"])
        return
    if not args.output:
        cli.error("A trial requires --output.")
    p = dict(DEFAULT)
    if args.config:
        record = json.loads(args.config.read_text())
        if "stages" in record:
            if args.fit == "local":
                if args.stage not in (3, 4):
                    cli.error("The local fit applies to standing or supported stepping.")
                key = "posture_parameters" if args.stage == 3 else "ground_parameters"
                record = {"parameters": record["ground_and_standing_comparisons"][args.design][key]}
            else:
                record = record["stages"][f"{args.design}_{args.stage}"][args.fit]
                if args.fit == "refined":
                    record = record["finalists"][args.rank]
        p.update(record["parameters"])
    atlas = json.loads(args.atlas.read_text()) if args.atlas else None
    result = evaluate(
        args.design,
        args.stage,
        p,
        args.dt,
        atlas,
        args.reserve,
        args.pose_offset,
        args.output,
        args.duration,
    )
    Path(str(args.output) + ".json").write_text(json.dumps(result, indent=2))
    print(
        json.dumps(
            {
                k: result.get(k)
                for k in (
                    "design",
                    "stage",
                    "complete",
                    "completed_s",
                    "reason",
                    "failed_routes",
                    "movement",
                    "accepted",
                    "wall_s",
                )
            }
        )
    )


if __name__ == "__main__":
    main()
