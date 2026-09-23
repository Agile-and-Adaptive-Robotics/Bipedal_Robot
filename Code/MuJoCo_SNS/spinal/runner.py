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
                     [--no-ground] [--leg-damping X] [--view] [--scope]
                     [--realtime] [--adaptive-tol X] [--contact-damp V]
--adaptive-tol X (default 0 = off): error-controlled integration - any
        2 ms step whose step-doubling error estimate exceeds X deg is
        refined to 4 x 0.5 ms (goal-3, 2026-09-23).
--contact-damp V (off by default): V = lessviscous | nonlinear - the two
        stable impedance-framework contact variants (goal-3, 2026-09-23).
--view  opens MuJoCo's interactive 3D window (left-drag orbit, scroll
        zoom, double-click a body to track); sim runs at FULL SPEED,
        close the window to end early. --realtime = 1x playback.
--scope opens the live neural strip-chart window (can combine with
        --view).
Outputs spinal_run.npz + spinal_run.png next to this file and prints an
honest summary (did it stand, did it step, did it fall).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

import build_network as bn
from muscle_map import EXTENSOR_STANCE_GROUPS, classify
from params import AFF, BAL, DT, E_HI, G, MOD, PHASE_RESET, SCHEDULE

import mujoco

from scipy.optimize import nnls

HERE = Path(__file__).parent
# AARL_MODEL env override (2026-09-18): point the whole chain at a new
# converted model (e.g. the scaled subject01 MJCF once it lands) without
# editing code. Unset = the proven default.
import os as _os
MODEL = Path(_os.environ.get(
    "AARL_MODEL",
    r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
    r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"))

KEY_ACTS = ("vas_lat_r", "med_gas_r", "soleus_r", "tib_ant_r", "psoas_r",
            "semimem_r", "glut_max2_r", "rect_fem_r", "sar_r", "grac_r",
            "vas_lat_l", "med_gas_l", "soleus_l", "tib_ant_l", "psoas_l",
            "semimem_l", "glut_max2_l", "rect_fem_l", "sar_l", "grac_l")
KEY_JOINTS = ("pelvis_tilt", "pelvis_list", "pelvis_rotation",
              "hip_flexion_r", "knee_angle_r", "ankle_angle_r",
              "hip_adduction_r", "subtalar_angle_r",
              "hip_flexion_l", "knee_angle_l", "ankle_angle_l",
              "hip_adduction_l", "subtalar_angle_l")

WARMUP = 0.6   # s: pose pinned while muscle activations build from zero

# Ben's prune list: tiny short-rotator muscles, problematic even in OpenSim
PRUNE_MUSCLES = {"quad_fem_r", "quad_fem_l", "gem_r", "gem_l",
                 "peri_r", "peri_l"}


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


def _pm_window(phi, lo=0.62, hi=0.95, edge=0.05):
    """Raised-cosine swing window of the contact-reset phase machine
    (1 inside [lo, hi], cosine ramps in the edge bands, 0 outside)."""
    if phi < lo - edge or phi > hi + edge:
        return 0.0
    if phi < lo:
        return float(0.5 * (1.0 - np.cos(np.pi
                                         * (phi - (lo - edge)) / edge)))
    if phi > hi:
        return float(0.5 * (1.0 - np.cos(np.pi
                                         * ((hi + edge) - phi) / edge)))
    return 1.0


# --- goal-3 fidelity flags (2026-09-23): error-controlled substepping and
# contact-damping variants, OFF by default (adaptive_tol = 0 keeps the
# exact single mj_step call; contact_damp = None leaves solref/solimp
# untouched) so every recorded winner reproduces bit-identically. See
# reports_20260923/goal3_fidelity_variants.md for the measured demos.
_ADAPT_ADR: dict[int, dict[str, int]] = {}


def _adaptive_step(model, data, tol_deg: float) -> None:
    """One 2 ms control step, error-controlled (step doubling): probe
    one 2 ms step vs two 1 ms half-steps from the same state; when the
    key-joint discrepancy exceeds tol_deg (deg), take 4 x 0.5 ms instead.
    ctrl is held across substeps (production ZOH semantics; activation
    states integrate through the finer steps)."""
    adr = _ADAPT_ADR.get(id(model))
    if adr is None:
        adr = {jn: model.jnt_qposadr[model.joint(jn).id]
               for jn in KEY_JOINTS}
        _ADAPT_ADR[id(model)] = adr

    def _save():
        return (data.qpos.copy(), data.qvel.copy(), data.act.copy(),
                float(data.time), data.qacc_warmstart.copy())

    def _restore(s):
        data.qpos[:] = s[0]
        data.qvel[:] = s[1]
        data.act[:] = s[2]
        data.time = s[3]
        data.qacc_warmstart[:] = s[4]

    def _qkey():
        return np.array([data.qpos[a] for a in adr.values()])

    buf = _save()
    model.opt.timestep = DT
    mujoco.mj_step(model, data)
    q_coarse = _qkey()
    _restore(buf)
    model.opt.timestep = DT / 2
    for _ in range(2):
        mujoco.mj_step(model, data)
    q_probe = _qkey()
    err = float(np.max(np.abs(np.degrees(q_probe - q_coarse))))
    if err > tol_deg:
        _restore(buf)
        model.opt.timestep = DT / 4
        for _ in range(4):
            mujoco.mj_step(model, data)


def _apply_contact_damp(model, variant: str) -> None:
    """Contact-damping variants INSIDE the impedance framework (the
    demonstrated stable ones; the naive direct-(k, b) negative-solref
    route is UNSTABLE at 2 ms - see goal3_fidelity_variants.md)."""
    if variant == "lessviscous":
        model.pair_solref[:, 0] = 0.01   # timeconst 0.02 -> 0.01
        model.pair_solref[:, 1] = 0.70   # dampratio 1 -> 0.7
    elif variant == "nonlinear":
        model.pair_solimp[:, 1] = 0.90   # midpoint 0.95 -> 0.90
        model.pair_solimp[:, 3] = 0.003  # width 0.001 -> 0.003
    else:
        raise ValueError(f"unknown --contact-damp variant: {variant}")


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
    # fit ONLY the DoF rows muscles can act on. Pelvis translation rows
    # carry the full body weight (the rig holds them, muscles can't touch
    # them) - fitting them injects a huge un-fittable residual that the
    # ridge spreads over every muscle, collapsing the solve to a useless
    # all<=0.2 pattern. Joint rows only = classic joint-space static opt.
    rows = np.abs(A).max(axis=1) > 1e-9
    A, tau = A[rows], tau[rows]
    # tonic co-contraction preload: net extension torque at hips/knees,
    # plantarflexion at ankles, extension at lumbar (signs follow the
    # OpenSim conventions of this model) - standing has real muscle tone,
    # it is not a zero-torque equilibrium
    def jadr(jn):
        return model.joint(jn).qposadr[0]
    for side in ("r", "l"):
        tau[jadr(f"hip_flexion_{side}")] += -40.0   # hip extension
        tau[jadr(f"knee_angle_{side}")] += -60.0    # knee extension
        tau[jadr(f"ankle_angle_{side}")] += 40.0    # plantarflexion
    for jn in ("lumbar_extension",):
        tau[jadr(jn)] += -20.0

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
    lam = 0.25
    Aaug = np.vstack([Aff, lam * np.eye(model.nu)])
    baug = np.concatenate([tau * 1e-2, lam * x_prior])
    x, _ = nnls(Aaug, baug)
    return np.clip(x, 0.0, 1.0)


def patch_xml(xml_text: str, model, data, kxy=2.0e5, kz=2.0e5, ky=5.0e5,
              dxy=2000.0, dy=3000.0, no_ground=False, leg_damping=None,
              boundmass=0.1, pp_contact=False, pin_rot=False, rig_scale=1.0):
    """All XML repairs/patches, as one function over the raw MJCF text.

    Returns the patched text; the caller must compile it from a file inside
    the model directory (mesh paths are relative). `no_ground` kills ground
    contact for suspended-in-air tests.
    """
    import re
    # stiff joint springs + explicit Euler at dt=5 ms explode; implicitfast
    # handles the joint-damping/spring terms stably
    xml_text = xml_text.replace(
        '<option timestep="0.005" collision="predefined"/>',
        '<option timestep="0.002" collision="predefined" integrator="implicitfast"/>')
    # converted bone inertias violate the triangle inequality on some bodies
    # (singular mass matrix -> QACC explosions); rebalance + bound the
    # massless pathpoint bodies. boundmass 0.01 (MyoSuite's own value) is
    # NOT enough at human muscle forces: with 92 muscles pulling on the
    # equality-coupled pathpoint DoFs the mass matrix goes singular
    # ("Inertia matrix too close to singular" -> NaN, diag_stab.py found
    # it at 0.3-1.0 co-activation). 0.1 survives 100% whole-body
    # co-contraction for 6 s. The pathpoint bodies carry only slide joints,
    # so boundinertia is irrelevant to them.
    xml_text = xml_text.replace(
        '<compiler angle="radian" autolimits="true"/>',
        f'<compiler angle="radian" autolimits="true" balanceinertia="true" '
        f'boundmass="{boundmass}" boundinertia="1e-6"/>')

    # ---- repair 1: hip joint axes are flipped vs the OpenSim conventions
    # (audit_signs.py: glut_max pulls flexion, iliacus/psoas pull extension;
    # same for adduction). Negate the hinge axes so + = flexion / adduction
    # as the muscle_map tables assume. Ranges are symmetric (+-2.094).
    for jn, old, new in (
            ("hip_flexion_r",  'axis="0 0 1"',  'axis="0 0 -1"'),
            ("hip_flexion_l",  'axis="0 0 1"',  'axis="0 0 -1"'),
            ("hip_adduction_r", 'axis="1 0 0"', 'axis="-1 0 0"'),
            ("hip_adduction_l", 'axis="-1 0 0"', 'axis="1 0 0"')):
        pat = rf'(<joint name="{jn}"[^/]*?){re.escape(old)}'
        xml_text, n = re.subn(pat, rf'\1{new}', xml_text, count=1)
        assert n == 1, f"axis flip failed for {jn}"

    # ---- repair 2: prune the tiny short-rotator muscles that are
    # problematic even in OpenSim (Ben): quad_fem, gem, peri (both sides).
    # Zero their force scale (gainprm[2]) and passive (biasprm[2]) so the
    # actuator is inert; sites/tendons stay in place.
    def _zero_fmax(match):
        tag = match.group(0)
        name = re.search(r'name="([^"]+)"', tag).group(1)
        if name in PRUNE_MUSCLES:
            tag = re.sub(r'gainprm="[^"]+"',
                         lambda p: 'gainprm="' + " ".join(
                             v if i != 2 else "0" for i, v in
                             enumerate(p.group(0).split('"')[1].split())) + '"',
                         tag)
            tag = re.sub(r'biasprm="[^"]+"',
                         lambda p: 'biasprm="' + " ".join(
                             v if i != 2 else "0" for i, v in
                             enumerate(p.group(0).split('"')[1].split())) + '"',
                         tag)
        return tag

    xml_text = re.sub(r'<general name="[^"]+" class="muscle"[^/]*/>',
                      _zero_fmax, xml_text)

    # ---- repair 2c (Ben's go 2026-09-13): the converter shipped 8 trunk
    # muscles with ACTIVE Fmax = 1 newton (gainprm[2]) while stock
    # gait2392 has ercspn 2500 N, intobl/extobl 900 N, ext_hal 162 N
    # (verified against gait2392_simbody.osim by _fmax_audit.py; the
    # passive slot biasprm[2] already carries the stock value — only the
    # active scale is broken). Before this fix the IMU trunk controller
    # (BAL_TRK -> ercspn) drove 1-newton muscles and the rig springs did
    # all the trunk work. NOTE: lengthrange="0.01 1" on these 8 is also
    # converter garbage (healthy muscles get their true RoM) - flagged,
    # not fixed here; the realized force is measured in _torque_budget.py.
    TRUNK_FMAX_FIX = {"ercspn": 2500.0, "intobl": 900.0, "extobl": 900.0,
                      "ext_hal": 162.0}   # stock gait2392_thelen2003

    def _fix_fmax(match):
        nonlocal n_fmax_fix
        tag = match.group(0)
        name = re.search(r'name="([^"]+)"', tag).group(1)
        base = name.rsplit("_", 1)[0]
        if base in TRUNK_FMAX_FIX:
            def _sub(p):
                vals = p.group(1).split()
                vals[2] = repr(float(TRUNK_FMAX_FIX[base]))
                return 'gainprm="' + " ".join(vals) + '"'
            new = re.sub(r'gainprm="([^"]+)"', _sub, tag, count=1)
            if new != tag:
                n_fmax_fix += 1
            tag = new
        return tag

    n_fmax_fix = 0
    xml_text = re.sub(r'<general name="[^"]+" class="muscle"[^/]*/>',
                      _fix_fmax, xml_text)
    assert n_fmax_fix == 8, f"trunk Fmax fix matched only {n_fmax_fix}"

    # ---- repair 3: patella mechanism for rect_fem (Ben). Dynamic audit:
    # rect_fem drives the knee into FLEXION (+894 rad/s^2 at full act) - its
    # route lacks the patella wrap the vastii kept. Route it over the
    # vastii's patella-tracking via point instead:
    # P1(origin) -> vas_med-P4 (moving via point on the patella path)
    # -> P3 (tibial tuberosity).
    for side in ("r", "l"):
        pat = (rf'(<spatial name="rect_fem_{side}_tendon">\s*'
               rf'<site site="rect_fem_{side}_rect_fem_{side}-P1"/>\s*'
               rf')<site site="rect_fem_{side}_rect_fem_{side}-P2"/>')
        xml_text, n = re.subn(
            pat, rf'\1<site site="vas_med_{side}_vas_med_{side}-P4"/>',
            xml_text, count=1)
        assert n == 1, f"rect_fem reroute failed for {side}"

    # ---- repair 2b (stability): pathpoint marker spheres must never
    # collide - they are virtual tendon via-points on ~0.01 kg boundmass
    # bodies, and sphere-sphere contacts there produce violent contact
    # forces (F/0.01 kg) that blow up the leg DoFs.
    if not pp_contact:
        xml_text = xml_text.replace('contype="2" conaffinity="2"',
                                    'contype="0" conaffinity="0"')

    # ---- repair 2f (CRITICAL): the equality-coupled "conditional
    # pathpoint" slide joints carry range limits sampled only around the
    # keyframe (MyoConverter artifact). The equality drives them exactly
    # along their polycoef trajectories, which leave those boxes within
    # ~10-30 deg of knee flexion (measured: 14 of 36 coupled dofs outside
    # their range at knee=10 deg, 31 at 30 deg - _knee_sweep.py). The
    # range-vs-equality fight stalls the knee AND slams huge constraint
    # forces into the leg DoFs (leg-DoF NaNs, dynamics not kinematics).
    # The equality alone defines the trajectory -> unlimit the followers.
    # NOTE: the default class sets limited="true" explicitly, so dropping
    # the range attribute alone would compile to range [0,0] and error;
    # limited="false" must be written explicitly. Also add armature to the
    # followers: with the CPG co-contraction the 92-muscle mass matrix goes
    # singular at these tiny-mass dofs ("Inertia matrix too close to
    # singular at DOF 23", t=8.9 s) - boundmass alone (0.1) is not enough.
    n_eq = 0
    for i in range(model.neq):
        if model.eq_type[i] != mujoco.mjtEq.mjEQ_JOINT:
            continue
        j1 = int(model.eq_obj1id[i])
        nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j1)
        pat = rf'<joint name="{nm}"([^/>]*?) range="[^"]*"'
        new_text, n = re.subn(
            pat, rf'<joint name="{nm}"\1 limited="false" armature="1.0"',
            xml_text, count=1)
        if n:
            xml_text = new_text
            n_eq += 1
    assert n_eq >= 30, f"pathpoint range strip matched only {n_eq} joints"

    # ---- repair 2f2 (Ben, 2026-09-13): ENFORCE the driver-joint RoM.
    # All 14 driver hinges shipped limited="false" (ranges inert), so the
    # swing knee hyperextended to +26 deg against the stock +10 deg cap
    # (flexion-negative convention). Flip limited="true" with the
    # CONVERTED/stock ranges (knee [-120,+10], hips +-120, ankle/subtalar/
    # mtp +-90 - the stock gait2392 values themselves). The equality-
    # coupled pathpoint followers stay limited="false" (repair 2f): they
    # are driven along their polycoefs and their sampled ranges are only
    # valid near the keyframe.
    rom_joints = ("knee_angle", "hip_flexion", "ankle_angle",
                  "subtalar_angle", "mtp_angle", "hip_adduction",
                  "hip_rotation")
    n_lim = 0
    for side in ("r", "l"):
        for jn in rom_joints:
            pat = rf'<joint name="{jn}_{side}"([^/>]*?) limited="false"'
            xml_text, k = re.subn(
                pat, rf'<joint name="{jn}_{side}"\1 limited="true"',
                xml_text, count=1)
            n_lim += k
    assert n_lim == 14, f"RoM limit patch matched only {n_lim}/14"

    # ---- repair 2g: knee DRIVER range — KEEP AS CONVERTED. REVERTED
    # 2026-09-10 evening after Ben's correction + direct measurement: the
    # OpenSim gait2392 range [-2.094, 0.1745] is flexion-NEGATIVE and the
    # converter PRESERVED that convention (an earlier note here claimed the
    # axis was flipped; that was based on a kinematics sweep that rotated
    # the knee without letting the coupled translation dofs follow - the
    # tibia's rolling center - which mirrored the foot path; invalid).
    # Ground truth (real-actuator test, unjammed followers): gravity
    # buckles the standing knee NEGATIVE (flexion), semimem/bifemsh/
    # med_gas drive NEGATIVE (flexion, to -77 deg), vas_lat/rect_fem drive
    # POSITIVE (extension, limited to +10 deg). Do not flip this range.

    # ---- repair 2c/2e: contact-set surgery. The converter ships
    # collision="predefined" with 19 explicit ground-vs-body <pair> lines.
    # Predefined pairs IGNORE contype/conaffinity, so muting the ground
    # geom does nothing: air mode must REMOVE the pairs, and on the ground
    # only the FEET should touch (OpenSim gait2392 has no self-collision
    # and no shank/pelvis dragging either).
    cm = re.search(r'<contact>.*?</contact>', xml_text, flags=re.S)
    assert cm, "contact section not found"
    if no_ground:
        xml_text = xml_text.replace(cm.group(0), "")
    else:
        keep = [ln.strip() for ln in cm.group(0).splitlines()
                if "ground-plane" in ln
                and any(f in ln for f in ("talus", "calcn", "toes"))]
        assert len(keep) == 6, f"expected 6 foot-ground pairs, got {len(keep)}"
        newsec = ("<contact>\n    " + "\n    ".join(keep) + "\n  </contact>")
        xml_text = xml_text.replace(cm.group(0), newsec)

    # ---- repair 2d (optional, stability sweep): raise leg-hinge damping
    # above the converted default 0.05
    if leg_damping is not None:
        n_leg = 0
        for side in ("r", "l"):
            for jn in ("hip_flexion", "hip_adduction", "hip_rotation",
                       "knee_angle", "ankle_angle", "subtalar_angle",
                       "mtp_angle"):
                pat = rf'<joint name="{jn}_{side}"[^/]*/>'
                m = re.search(pat, xml_text)
                if not m:
                    continue
                tag = m.group(0)
                if 'damping="' in tag:
                    tag = re.sub(r'damping="[^"]*"', f'damping="{leg_damping}"',
                                 tag)
                else:
                    tag = tag[:-2] + f' damping="{leg_damping}"/>'
                xml_text = xml_text.replace(m.group(0), tag)
                n_leg += 1
        assert n_leg >= 10, f"leg damping patch matched only {n_leg} joints"

    # ---- rig: spring-locked pelvis translations at the keyframe pose,
    # now WITH explicit viscous damping (previously accepted but never
    # written -> an undamped spring; implicitfast alone doesn't damp it)
    specs = {"pelvis_tx": (kxy, dxy), "pelvis_ty": (ky, dy),
             "pelvis_tz": (kz, dxy)}
    # pelvis ORIENTATION springs: always on. The yaw/roll/pitch hinges have
    # zero damping in the converted model; with only foot friction to resist
    # (and left-right asymmetries in any solved standing pattern) the pelvis
    # freely spins about vertical - measured 40-68 rad/s at t=1.4 s on the
    # ground, ending in friction-force blowup. This is a support rig; the
    # springs are the rig's torso stabilize. Softer than the air pin.
    # Air mode (pin_rot) effectively RIGID-clamps the trunk (Ben
    # 2026-09-11: "clamp the head so the whole trunk stops falling"):
    # 150 N·m/rad let the trunk pike 30 deg, 800 still gave 13 deg;
    # 5000 + heavy damping = a clamp. Ground mode stays springy (400).
    rot_k = 5000.0 if pin_rot else 400.0
    rot_c = 300.0 if pin_rot else 40.0
    specs.update({"pelvis_tilt": (rot_k, rot_c),
                  "pelvis_list": (rot_k, rot_c),
                  "pelvis_rotation": (rot_k, rot_c)})
    # trunk (lumbar) + femur long-axis-rotation springs: same rig family.
    # Ben spotted the torso flipping upside down over the lumbar hinge in
    # the viewer (2026-09-10): the lumbar dofs have damping 0.05 and our
    # erector-spinae/oblique drive is far too weak to hold the ~30 kg
    # torso (real anatomy has muscles+ligaments doing this; our spinal
    # pattern does not yet). hip_rotation was free-spinning the legs too.
    # These are RIG elements - wean them alongside the pelvis springs
    # when balance is solved.
    lumb_k = 3000.0 if pin_rot else 300.0
    hipr_k = 300.0 if pin_rot else 100.0
    lumc = 150.0 if pin_rot else 25.0
    for jn in ("lumbar_extension", "lumbar_bending", "lumbar_rotation",
               "hip_rotation_r", "hip_rotation_l"):
        if mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn) >= 0:
            specs.update({jn: (lumb_k if "lumbar" in jn else hipr_k, lumc)})
    # subtalar + mtp + ANKLE passive ligament surrogate: the converted model
    # has NO passive stiffness on the foot joints (OpenSim coordinate stops
    # did not survive conversion) - the foot flopped to 128 deg of subtalar
    # and -97 deg of ankle in the air runs. Weak spring + damping =
    # ligament restraint (human ankle passive stiffness ~ 5-15 N·m/rad).
    for side in ("r", "l"):
        for jn in (f"subtalar_angle_{side}", f"mtp_angle_{side}",
                   f"ankle_angle_{side}",
                   f"hip_adduction_{side}"):   # capsule stiffness: the
                    # frontal plane was a free pendulum (legs splayed
                    # +-25 deg adduction in air - Ben 2026-09-11)
            if mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn) >= 0:
                specs.update({jn: (10.0 if "hip" not in jn else 30.0, 3.0)})
    # weaning (--rig-scale S): scale every RIG element (pelvis support +
    # rotation stabilizers + lumbar/hip-rotation rig), stiffness by S and
    # damping by sqrt(S) (constant damping ratio). The ligament surrogates
    # (subtalar/mtp/ankle/hip_adduction) are MODEL properties, not rig -
    # they stay.
    if rig_scale != 1.0:
        rs, rc = rig_scale, rig_scale ** 0.5
        for jn in ("pelvis_tx", "pelvis_ty", "pelvis_tz",
                   "pelvis_tilt", "pelvis_list", "pelvis_rotation",
                   "lumbar_extension", "lumbar_bending", "lumbar_rotation",
                   "hip_rotation_r", "hip_rotation_l"):
            if jn in specs:
                k, c = specs[jn]
                specs[jn] = (k * rs, c * rc)
    for jn, (k, c) in specs.items():
        ref = data.qpos[model.joint(jn).qposadr[0]]
        pat = rf'<joint name="{jn}"[^/]*?/>'
        m = re.search(pat, xml_text)
        assert m, f"pelvis joint {jn} not found"
        old_tag = m.group(0)
        # build the fully-modified tag FIRST, then replace the original
        # (modifying before searching silently no-ops the replace)
        new_tag = re.sub(r'damping="[^"]*"', f'damping="{c}"', old_tag)
        new_tag = new_tag[:-2] + f' stiffness="{k}" springref="{ref}"/>'
        xml_text = xml_text.replace(old_tag, new_tag)
    return xml_text


def apply_harness(model, data, kxy=2.0e5, kz=2.0e5, ky=5.0e5,
                  dxy=2000.0, dy=3000.0, no_ground=False, leg_damping=None,
                  boundmass=0.1, pp_contact=False, pin_rot=False,
                  rig_scale=1.0):
    """Rigid pelvis rig (v1 default): the pelvis is spring-locked at the
    keyframe pose in all three translations - the classic biped walker test
    rig. Legs swing under it with ground contact; balance is NOT solved.
    Weaker springs (a compliant rehab tether) explode the stiff contact
    model at dt=5 ms; see DESIGN.md open problems.
    """
    xml_text = patch_xml(MODEL.read_text(encoding="utf-8"), model, data,
                         kxy=kxy, kz=kz, ky=ky, dxy=dxy, dy=dy,
                         no_ground=no_ground, leg_damping=leg_damping,
                         boundmass=boundmass, pp_contact=pp_contact,
                         pin_rot=pin_rot, rig_scale=rig_scale)
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


def capture_pose(model, data=None) -> dict[str, np.ndarray]:
    """Joint name -> starting qpos values (for re-seeding reduced models).
    Reads data.qpos when given (the applied start pose), else the model's
    stored keyframe."""
    pose = {}
    for j in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        adr, narange = model.jnt_qposadr[j], 7 if model.jnt_type[j] == 0 else 1
        src = data.qpos if data is not None else model.key_qpos[0]
        pose[name] = np.asarray(src[adr:adr + narange]).copy()
    return pose


# ---- Ben (2026-09-14): START POSE = the model's own "normal"
# experimental pose (gait2392_simbody Coordinates panel / normal.mot;
# values in DEGREES, OpenSim conventions). Our repair-1 hip axis flips
# make OpenSim +flexion/+adduction == MuJoCo +qpos, so these map
# directly. Slight forward tilt, knees slightly flexed, right leg
# forward (hip +24.6), left leg back (hip -16.6, ankle +9.8 dorsi).
START_POSE_DEG = {
    "pelvis_tilt": -1.870,
    "pelvis_list": -0.470,
    "pelvis_rotation": 2.430,
    "hip_flexion_r": 24.610,
    "hip_adduction_r": 1.090,
    "hip_rotation_r": -1.350,
    "knee_angle_r": -3.940,
    "ankle_angle_r": -1.700,
    "subtalar_angle_r": 0.000,
    "mtp_angle_r": 0.000,
    "hip_flexion_l": -16.560,
    "hip_adduction_l": 2.680,
    "hip_rotation_l": 1.280,
    "knee_angle_l": -8.200,
    "ankle_angle_l": 9.810,
    "subtalar_angle_l": 0.000,
    "mtp_angle_l": 0.000,
    "lumbar_extension": 1.870,
    "lumbar_bending": 0.470,
    "lumbar_rotation": -2.430,
}
START_PELVIS_HEIGHT = 0.92   # m. NOT the OpenSim 0.96: measured with
# mesh-vertex foot bottoms (_foot_flat_check.py), our converted foot
# geometry leaves the soles 3-6 cm ABOVE ground at 0.96 (rig-locked ->
# feet dangle in passive plantarflexion = the tiptoe posture). Right
# sole touches at ty 0.903, left (back leg) at ~0.925; 0.92 = contact
# compromise the knees/soft contact absorb.

# 2026-09-21: symmetric double-stance variant (AARL_POSE=symmetric).
# Hypothesis under test: the s3c/s3d/s3e "frozen left leg" is a POSE
# LATCH - normal.mot starts asymmetric (right hip +24.6 FORWARD, left
# -16.6 BACK = pre-planted), and the planted leg's load->extensor-
# afferent->more-extensor loop is unbreakable by any crossed trigger
# searched. Symmetric start gives neither leg the planted role.
START_POSE_DEG_SYMMETRIC = {
    "pelvis_tilt": -1.870,
    "pelvis_list": -0.470,
    "pelvis_rotation": 2.430,
    "hip_flexion_r": 0.0,
    "hip_adduction_r": 1.090,
    "hip_rotation_r": -1.350,
    "knee_angle_r": -10.0,
    "ankle_angle_r": -1.700,
    "subtalar_angle_r": 0.0,
    "mtp_angle_r": 0.0,
    "hip_flexion_l": 0.0,
    "hip_adduction_l": 2.680,
    "hip_rotation_l": 1.280,
    "knee_angle_l": -10.0,
    "ankle_angle_l": -1.700,
    "subtalar_angle_l": 0.0,
    "mtp_angle_l": 0.0,
    "lumbar_extension": 1.870,
    "lumbar_bending": 0.470,
    "lumbar_rotation": -2.430,
}


def _project_followers(model, data):
    """Re-project equality-coupled pathpoint dofs onto their driver polycoefs
    (copy of bsolve_ik.apply_eq_followers; kept local to avoid a circular
    import - bsolve_ik imports runner)."""
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


def apply_start_pose(model, data) -> None:
    """Overwrite the keyframe with the model's own 'normal' pose
    (normal.mot values, degrees; drivers + lumbar + pelvis orientation,
    then re-project the equality followers). Applied AS GIVEN - no sign
    auto-correction: these are the OpenSim model's own values, and the
    IMU-lean heuristic's sign convention is not the coordinate's."""
    jadr = {jn: model.joint(jn).qposadr[0] for jn in START_POSE_DEG
            if model.joint(jn).id >= 0}
    import os as _os
    _pose = _os.environ.get("AARL_POSE", "normal")
    pose_src = (START_POSE_DEG_SYMMETRIC if _pose == "symmetric"
                else START_POSE_DEG)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    for jn, deg in pose_src.items():
        if jn in jadr:
            data.qpos[jadr[jn]] = np.radians(deg)
    # 2026-09-20 (Ben: "lower the walker enough to make contact with the
    # ground"): env override for the pelvis height. The 0.92 default is a
    # compromise that left the LEFT sole ~5 mm short -> the left foot
    # never loaded (0% contact in the s3b winner) and its knee pinned at
    # the +10 extension stop. The rig anchors at whatever height is set
    # here (apply_harness reads it from data), so this single knob sets
    # the standing height the whole run is sprung to.
    import os as _os
    _ty = _os.environ.get("AARL_PELVIS_TY")
    ty = float(_ty) if _ty is not None else START_PELVIS_HEIGHT
    data.qpos[model.joint("pelvis_ty").qposadr[0]] = ty
    _project_followers(model, data)
    mujoco.mj_forward(model, data)
    torso = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "torso")
    up = data.xmat[torso].reshape(3, 3)[:, 1]
    lean = float(np.degrees(np.arctan2(-up[0], max(up[2], 1e-6))))
    print(f"start pose ({_pose}): applied {len(jadr)} coordinates as "
          f"given; pelvis ty {ty:.3f} m (env override"
          f"{'=' + _ty if _ty is not None else ' unset'}); "
          f"IMU-metric torso lean {lean:+.1f} deg (+ = back)",
          flush=True)


def seed_pose(model, data, pose: dict[str, np.ndarray]):
    mujoco.mj_resetData(model, data)
    for j in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        if name in pose:
            adr = model.jnt_qposadr[j]
            data.qpos[adr:adr + len(pose[name])] = pose[name]


def main(argv):
    import mujoco  # noqa: F811 -- the later `import mujoco.viewer` makes mujoco function-local
    walk_drive = 4.0
    harness = 1500.0     # default ON (balance assist); 0 disables
    use_aff = True
    no_ground = False
    leg_damping = None
    view = False
    scope_on = False
    realtime = False
    interleg = True
    eval_mode = False
    rig_scale = 1.0
    straight_start = False
    adaptive_tol = 0.0      # --adaptive-tol X: error-controlled substepping
    contact_damp = None     # --contact-damp lessviscous|nonlinear
    args = list(argv)
    fit_keys = None
    if "--fitted" in args:
        # load the IK/NNLS back-solve refit (fit_pf.py output) BEFORE any
        # params read: full W_PF_MN + W_POSTURE tables from human data
        import json as _json
        with open(HERE / "fitted_walk_params.json", encoding="utf-8") as f:
            fit = _json.load(f)
        import params as _p
        for ph, tbl in fit["W_PF_MN"].items():
            for g, w in tbl.items():
                _p.W_PF_MN[ph][g] = float(w)
        for g, w in fit["W_POSTURE"].items():
            _p.W_POSTURE[g] = float(w)
        # which entries came from the fitted file: --best's pf_gain must
        # scale ONLY these (the params-default trunk entries the json
        # lacks - E1/E2 trunk_ext, F2 trunk_flex, W_POSTURE trunk_ext -
        # are NOT part of the back-solved table and scaling them silently
        # diverged --best from the study's set_params config by kine 0.64;
        # found by the v5 regression gate 2026-09-12)
        fit_keys = ({ph: set(tbl) for ph, tbl in fit["W_PF_MN"].items()},
                    set(fit["W_POSTURE"]))
        args.remove("--fitted")
        print("loaded fitted_walk_params.json (back-solve refit)", flush=True)
    if "--best10" in args:
        # v10 winner file (transient PRESET + ankle trim); same loader
        best_path = HERE / "best_walk_params_v10.json"
        args.remove("--best10")
    elif "--best9" in args:
        # v9 winner file (flat-foot height + amplitude objective);
        # same loader as --best
        best_path = HERE / "best_walk_params_v9.json"
        args.remove("--best9")
    elif "--best8b" in args:
        # v8b winner file (normal.mot start pose retune); same loader
        best_path = HERE / "best_walk_params_v8b.json"
        args.remove("--best8b")
    elif "--best8" in args:
        # v8 winner file (post start-pose retune); same loader as --best
        best_path = HERE / "best_walk_params_v8.json"
        args.remove("--best8")
    elif "--best7" in args:
        # v7 winner file (post-RoM-limit retune); same loader as --best
        best_path = HERE / "best_walk_params_v7.json"
        args.remove("--best7")
    elif "--best6" in args:
        # v6 winner file (optuna_walk_v6.py output); same loader as --best
        best_path = HERE / "best_walk_params_v6.json"
        args.remove("--best6")
    elif "--best5" in args:
        # v5 winner file (optuna_walk_v5.py output); same loader as --best
        best_path = HERE / "best_walk_params_v5.json"
        args.remove("--best5")
    elif "--best" in args:
        best_path = HERE / "best_walk_params.json"
        args.remove("--best")
    else:
        best_path = None
    if best_path is not None:
        # load the optimizer's winning parameters (optuna output) before
        # anything reads the params module
        import json as _json
        with open(best_path, encoding="utf-8") as f:
            _doc = _json.load(f)
        best = _doc["params"]
        # pf_gain lives at the DOCUMENT TOP LEVEL, not in params. The old
        # `if "pf_gain" in best` tested the params dict (which never has
        # it) so the gain-scaling branch NEVER ran: every --fitted --best
        # "reproduction" since v3 silently evaluated a pf_gain=1.0 config
        # (caught by the v5 regression gate 2026-09-12 - the study path
        # scored -61.174, the --best path -61.82, and the state dump
        # showed the raw unscaled fitted tables).
        gain = float(_doc.get("pf_gain", 1.0))
        import params as _p
        _p.TAU["rg_adapt"] = best["rg_adapt"]
        _p.G["descend_to_rg_e"] = best["desc_e"]
        _p.G["rg_to_pf"] = best["rg_to_pf"]
        _p.W_PF_MN["E2"]["ankle_pf"] = best["e2_pf"]
        _p.W_PF_MN["F1"]["ankle_df"] = best["f1_df"]
        _p.W_PF_MN["F1"]["knee_flex"] = best["f1_kf"]
        _p.W_POSTURE["knee_ext"] = best["post_kneext"]
        if "post_hipext" in best:
            _p.W_POSTURE["hip_ext"] = best["post_hipext"]
        if gain != 1.0:
            # v3+ winners store a global gain on the fitted table (the 5
            # knob values above are ALREADY effective - gain must not
            # double-apply to them, so it only scales the other entries).
            # Scale ONLY the fitted-file entries (see fit_keys above):
            # the params-default trunk keys are not part of the table.
            fk = fit_keys if fit_keys is not None else (
                {ph: set(tbl) for ph, tbl in _p.W_PF_MN.items()},
                set(_p.W_POSTURE))
            for ph, tbl in _p.W_PF_MN.items():
                for g in fk[0].get(ph, ()):
                    tbl[g] *= gain
            for g in fk[1]:
                _p.W_POSTURE[g] *= gain
            for kn in ("e2_pf", "f1_df", "f1_kf"):
                key = {"e2_pf": ("E2", "ankle_pf"),
                       "f1_df": ("F1", "ankle_df"),
                       "f1_kf": ("F1", "knee_flex")}[kn]
                _p.W_PF_MN[key[0]][key[1]] = best[kn]
            _p.W_POSTURE["knee_ext"] = best["post_kneext"]
            _p.W_POSTURE["hip_ext"] = best["post_hipext"]
        if "desc_f" in best:
            _p.G["descend_to_rg_f"] = float(best["desc_f"])
        if "e2_adapt" in best:
            _p.PF_SHAPE["E2"] = (_p.PF_SHAPE["E2"][0], float(best["e2_adapt"]))
        if "phase_reset_e" in best:
            _p.G["phase_reset_e"] = float(best["phase_reset_e"])
        if "phase_reset_f" in best:
            _p.G["phase_reset_f"] = float(best["phase_reset_f"])
        if "f1_kneext_inh" in best:
            _p.G["f1_kneext_inh"] = float(best["f1_kneext_inh"])
        if "f1_anklepf_inh" in best:
            _p.G["f1_anklepf_inh"] = float(best["f1_anklepf_inh"])
        if "renshaw" in best:
            _p.G["renshaw"] = float(best["renshaw"])
        if "ankle_post_walk_trim" in best:
            _p.G["ankle_post_walk_trim"] = float(best["ankle_post_walk_trim"])
        for _k in ("heel_rge", "toe_rge", "ib_rge", "ia_in",
                   "contact_onset", "contra_swing",
                   "contra_kinh", "pm_gain", "pm_T",
                   "pm_ws", "pm_add", "pm_aff",
                   "full_rules",
                   # goal2 balance stage (2026-09-23) - JSON RULE
                   "vest_ext", "vest_flex_inh", "vest_prop"):  # JSON RULE
            if _k in best:
                _p.G[_k] = float(best[_k])
        if "joint_pf" in best:
            # 2026-09-18 JSON RULE: joint-layer PF knob reproduces from
            # the winner json (searched only by later studies; default
            # studies leave it 0 = phase-cell PF)
            _p.G["joint_pf"] = float(best["joint_pf"])
        _p.BAL["kx"] = best["kx"]
        walk_drive = float(best.get("drive", walk_drive))
        print(f"loaded {best_path.name} (drive={walk_drive:.2f}, "
              f"phase_reset {_p.G['phase_reset_e']:.2f}/"
              f"{_p.G['phase_reset_f']:.2f})", flush=True)
    while args:
        a = args.pop(0)
        if a == "--drive":
            walk_drive = float(args.pop(0))
        elif a == "--adaptive-tol":
            # goal-3 (2026-09-23): error-controlled integration - 0 (off)
            # keeps the exact single-step path; > 0 refines any 2 ms step
            # whose step-doubling error estimate exceeds tol (deg).
            adaptive_tol = float(args.pop(0))
        elif a == "--contact-damp":
            # goal-3 (2026-09-23): impedance-framework contact variant
            contact_damp = args.pop(0) if args and \
                not args[0].startswith("--") else "lessviscous"
        elif a == "--joint-pf":
            # 2026-09-18: build the joint-layer PF (T1 fit) instead of
            # the phase cells; value 0/absent = phase cells
            import params as _p
            _p.G["joint_pf"] = float(args.pop(0)) if args and \
                args[0].replace(".", "").replace("-", "").isdigit() else 1.0
        elif a == "--harness":
            harness = float(args.pop(0)) if args and args[0].replace(".", "").isdigit() else 1500.0
        elif a == "--no-harness":
            harness = 0.0
        elif a == "--no-afferents":
            use_aff = False
        elif a == "--no-ground":
            no_ground = True
        elif a == "--no-interleg":
            interleg = False
        elif a == "--eval":
            # fast objective evaluation for the optimizer (optuna_walk.py):
            # 16 s schedule (1.0 stand / 1.0 ramp / 11.0 walk / 1.5
            # wind-up / 1.5 stand), quiet, no npz/png writes; main()
            # returns a metrics dict. The long walk window matters: slow
            # gaits (0.3-0.4 Hz) need >= 3 RG-E cycles in-window or
            # kine_ref returns None and the optimizer sees the frozen
            # sentinel (the ground_walk_v8b collapse, 2026-09-14).
            eval_mode = True
            SCHEDULE.update(stand1=(0.0, 1.0), ramp_up=(1.0, 2.0),
                            walk=(2.0, 13.0), ramp_down=(13.0, 14.5),
                            stand2=(14.5, 16.0))
        elif a == "--leg-damping":
            leg_damping = float(args.pop(0))
        elif a == "--rig-scale":
            # wean the support rig: stiffness *S, damping *sqrt(S)
            rig_scale = float(args.pop(0))
        elif a == "--phase-reset":
            # v5 sensory phase-reset gains: HIP_EXT_SIG -> RG-E exc / RG-F
            # inh; HIP_FLEX_SIG -> RG-F exc / RG-E inh (0 0 = v4 behavior)
            import params as _p
            _p.G["phase_reset_e"] = float(args.pop(0))
            _p.G["phase_reset_f"] = float(args.pop(0))
        elif a == "--kneext-inh":
            # phase-3 swing-knee quad suppression gain (0 = absent)
            import params as _p
            _p.G["f1_kneext_inh"] = float(args.pop(0))
        elif a == "--renshaw":
            # Renshaw recurrent inhibition gain (0 = absent, Deng A6)
            import params as _p
            _p.G["renshaw"] = float(args.pop(0))
        elif a == "--straight-start":
            # revert to the converted straight keyframe (no walking start
            # pose override)
            straight_start = True
        elif a == "--view":
            view = True
        elif a == "--scope":
            scope_on = True
        elif a == "--realtime":
            realtime = True
        elif a == "--time":
            SCHEDULE["stand2"] = (SCHEDULE["ramp_down"][1],
                                  SCHEDULE["ramp_down"][1] + float(args.pop(0)))
        elif a == "--stand-eval":
            # goal2 balance-stage eval (2026-09-23; SCONE Tutorial-3a
            # analog): STANDING ONLY, no walk window. Value = total
            # duration in s (default 8). eval_mode metrics then carry the
            # bal_* balance fields (sway radius, tilt envelope, contact
            # symmetry, fall flag) for the stage-4 objective.
            eval_mode = True
            _T = float(args.pop(0)) if args and \
                args[0].replace(".", "").replace("-", "").isdigit() else 8.0
            SCHEDULE.update(stand1=(0.0, _T), ramp_up=(_T, _T),
                            walk=(_T, _T), ramp_down=(_T, _T),
                            stand2=(_T, _T))
        elif a == "--vest":
            # goal2 vestibular-analog tone gain (VEST -> extensor MNs;
            # 0 / absent = topology absent, bit-identical)
            import params as _p
            _p.G["vest_ext"] = float(args.pop(0)) if args and \
                args[0].replace(".", "").replace("-", "").isdigit() else 0.5
        elif a == "--vest-flexinh":
            # goal2 reciprocal flexor inhibition from the VEST cells
            import params as _p
            _p.G["vest_flex_inh"] = float(args.pop(0)) if args and \
                args[0].replace(".", "").replace("-", "").isdigit() else 0.2
        elif a == "--vest-prop":
            # goal2 stance-gated II length-loop boost (SCONE T3a KL analog)
            import params as _p
            _p.G["vest_prop"] = float(args.pop(0)) if args and \
                args[0].replace(".", "").replace("-", "").isdigit() else 0.5

    # state-dump hook (debug): RUNNER_DUMP_STATE=<file> writes every
    # mutable param right after arg parsing - diff two paths to find
    # config divergence (found the v4b --best/set_params split this way)
    import os as _os
    # ---- Ben's connectome spec (connectome_editor.html export) ----
    # Applied BEFORE the network build: rule enabled -> its gain key
    # gets |gain|, disabled -> 0.0; any rule with hops >= 1 sets
    # full_rules (the interneuron-layer topology). You own the wiring;
    # the machine tunes.
    _cg = HERE / "connectome_gains.json"
    if _cg.exists():
        import json as _js
        spec = _js.loads(_cg.read_text(encoding="utf-8")).get(
            "rules", {})
        _any_in = False
        for _rid, rr in spec.items():
            if not rr.get("enabled", False):
                continue
            _gk = rr.get("gain_key")
            if _gk and _gk in _p.G:
                _p.G[_gk] = abs(float(rr.get("gain", 0.0)))
            if rr.get("hops", 0) >= 1:
                _any_in = True
        if _any_in:
            _p.G["full_rules"] = 1.0
        print(f"connectome_gains.json applied "
              f"({sum(1 for r in spec.values() if r.get('enabled'))} "
              f"rules, full_rules={_p.G['full_rules']})", flush=True)
    if _os.environ.get("RUNNER_DUMP_STATE"):
        import json as _js
        import params as _pp
        _js.dump(dict(
            W_PF_MN=_pp.W_PF_MN, W_POSTURE=_pp.W_POSTURE, TAU=_pp.TAU,
            G=_pp.G, PF_SHAPE=_pp.PF_SHAPE, BAL=_pp.BAL, AFF=_pp.AFF,
            MOD=_pp.MOD, walk_drive=walk_drive,
            SCHEDULE={k: list(v) for k, v in _pp.SCHEDULE.items()},
        ), open(_os.environ["RUNNER_DUMP_STATE"], "w", encoding="utf-8"),
            indent=1, default=float)

    model = mujoco.MjModel.from_xml_path(str(MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    if not straight_start:
        apply_start_pose(model, data)
    key_pose = capture_pose(model, data)

    # standing-activation solve ALWAYS runs on a ground-contact copy of the
    # model: in suspended-air configs there is no contact force to carry
    # gravity and the NNLS degenerates (all activations <= 0.2, wrong
    # muscles) - the ground-on solve is the physiologically real standing
    # pattern and gives the legs enough tone not to collapse in air either
    model_solve = apply_harness(model, data, kxy=1500.0)
    data_solve = mujoco.MjData(model_solve)
    seed_pose(model_solve, data_solve, key_pose)
    mujoco.mj_forward(model_solve, data_solve)

    if harness > 0:
        # original semantics: soft x tether (walks forward under it), rigid
        # y/z lock (apply_harness defaults kz=2e5, ky=5e5); in air the
        # pelvis ORIENTATION is pinned too (otherwise asymmetric CPG drive
        # tips the assembly head-over)
        # 2026-09-21: AARL_KY scales the LATERAL rig spring only - the
        # s3g diagnosis: ky=5e5 anchors the pelvis laterally and
        # physically blocks the weight transfer the swing gate needs
        # (right foot never reached 25% BW). ky_scale searched per trial.
        _ky = 5.0e5 * float(_os.environ.get("AARL_KY", "1.0"))
        model = apply_harness(model, data, kxy=harness, ky=_ky,
                              no_ground=no_ground, leg_damping=leg_damping,
                              pin_rot=no_ground, rig_scale=rig_scale)
    elif no_ground or leg_damping is not None or rig_scale != 1.0:
        # no rig, but XML patches needed (suspended-air / damping sweep /
        # full wean): zero-stiffness springs with negligible damping =
        # free pelvis; --rig-scale<1 keeps a scaled support instead
        s = rig_scale
        model = apply_harness(model, data, kxy=harness * s, kz=2.0e5 * s,
                              ky=5.0e5 * s, dxy=2000.0 * s ** 0.5,
                              dy=3000.0 * s ** 0.5,
                              no_ground=no_ground, leg_damping=leg_damping,
                              pin_rot=no_ground, rig_scale=rig_scale)
    data = mujoco.MjData(model)
    seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)
    if contact_damp is not None:
        _apply_contact_damp(model, contact_damp)
        print(f"contact damping variant applied: {contact_damp}", flush=True)

    # interactive 3D viewer (--view): open the window FIRST (right after
    # the model is ready) so it appears within ~2 s of launching - the
    # network build + standing solve that follow take ~30-60 s of silent
    # console time, and a window that only appears at the end looks like
    # a hang. The window shows the standing pose during the build, then
    # plays the run in slow motion (the sim is faster than real time).
    dur = SCHEDULE["stand2"][1]
    nsteps = int(dur / DT)
    t_wall0 = time.perf_counter()

    viewer = None
    scope = None
    if view:
        try:
            print("opening 3D viewer window...", flush=True)
            import mujoco.viewer
            viewer = mujoco.viewer.launch_passive(model, data)
            # side-view camera on the walker (CAUTION-wrapped: field names
            # vary across mujoco versions)
            try:
                viewer.cam.lookat[:] = (0.0, 0.0, 0.85)
                viewer.cam.distance = 3.2
                viewer.cam.azimuth = 90.0     # side view
                viewer.cam.elevation = -8.0
            except Exception:
                pass
            viewer.sync()
            print("viewer open - LEFT-DRAG orbits the camera, scroll "
                  "zooms, right-drag pans (double-click a body to track "
                  "it). Sim runs at full speed; add --realtime for 1x "
                  "playback. Close the window to stop the run early.",
                  flush=True)
        except Exception as e:
            print(f"(viewer unavailable: {e}; running headless)", flush=True)
    if scope_on:
        # live neural strip-chart (cheap: decimated redraws); can be used
        # alone or together with --view
        try:
            import neuro_scope as ns
            scope = ns.NeuroScope(
                nsteps, DT,
                ["DRIVE", "POSTURE", "RG_E_r", "RG_F_r",
                 "PF_E2_r", "PF_F1_r", "RG_E_l", "RG_F_l"],
                list(KEY_ACTS[:8]))
            print("neural scope open (separate window).", flush=True)
        except Exception as e:
            print(f"(neural scope unavailable: {e})", flush=True)

    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    net = bn.build(acts, dt=DT, interleg=interleg)
    aid = {name: i for i, name in enumerate(acts)}
    muscles = net.muscles

    # per-muscle normalization: length range and peak force. The converted
    # actuators are <general class="muscle"> with the muscle parameter vector
    # in gainprm/biasprm; Fmax sits at gainprm[2] (e.g. soleus 3395 N).
    Lrange = model.actuator_lengthrange  # [nu, 2]
    Lmid = Lrange.mean(axis=1)
    Lhalf = np.maximum((Lrange[:, 1] - Lrange[:, 0]) / 2, 1e-3)
    Fmax = np.maximum(model.actuator_gainprm[:, 2], 5.0)
    act_stand = solve_standing_activations(model_solve, data_solve, Fmax)
    print(f"standing solve: nnls x in [{act_stand.min():.2f}, "
          f"{act_stand.max():.2f}], >0.1: {(act_stand > 0.1).sum()} muscles")
    top = np.argsort(act_stand)[::-1][:8]
    print("  top: " + ", ".join(f"{acts[i]}={act_stand[i]:.2f}" for i in top))
    print(f"Fmax range: {Fmax.min():.0f}..{Fmax.max():.0f} N | "
          f"L range sample: {Lrange[aid['soleus_r']]}")

    iport = {p: net.input_index(p) for p in net.inputs}
    foot_r, foot_l = find_foot(model)
    # torso body id for the IMU (vestibular surrogate reads its orientation)
    torso_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "torso")
    lean_prev = 0.0
    # balance target = ankle axis (talus) x: the neutral standing COM
    ankle_x = []
    for i in range(model.nbody):
        nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i)
        if nm and nm.startswith("talus"):
            ankle_x.append(float(data.xpos[i][0]))
    x_ref = float(np.mean(ankle_x)) if ankle_x else float(data.subtree_com[0][0])
    y_ref = float(data.subtree_com[0][1])   # lateral balance reference
    jadr = {jn: model.jnt_qposadr[model.joint(jn).id] for jn in KEY_JOINTS
            if model.joint(jn).id >= 0}

    log_t = np.zeros(nsteps)
    log_act = np.zeros((nsteps, len(KEY_ACTS)))
    log_q = np.zeros((nsteps, len(KEY_JOINTS)))
    log_qfull = np.zeros((nsteps, model.nq))   # full state for renderers
    log_com = np.zeros((nsteps, 3))
    # full neural stack for the plot tool (plot_run.py): descending,
    # RG, PF group cells (right; phase cells OR joint-layer HCs),
    # balance cells (the HIP_*_SIG channels left with the PRESET
    # removal 2026-09-16)
    if G.get("joint_pf", 0.0) > 0.0:
        _pf_watch = ("PF_HIP-E_r", "PF_KNEE-E_r",
                     "PF_KNEE-F_r", "PF_ANK-F_r")
    else:
        _pf_watch = ("PF_E1_r", "PF_E2_r", "PF_F1_r", "PF_F2_r")
    NEURO_NAMES = ("DRIVE", "POSTURE",
                   "RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
                   *_pf_watch,
                   "BAL_TRK_FLX", "BAL_TRK_EXT", "BAL_LAT_R", "BAL_LAT_L")
    log_neuro = np.zeros((nsteps, len(NEURO_NAMES)))
    # 2026-09-20: per-side ground-truth contact (heel+toe normal force, N)
    # and foot height (min body-origin z of calcn/toes, m), logged EVERY
    # step in ground mode. Ben's critique: the eval objective and the
    # figures must see real contact - the old duty was neural and the
    # s3b "winner" never loaded its left foot (0% contact frames).
    log_contact = np.zeros((nsteps, 2))
    log_footz = np.zeros((nsteps, 2))
    foot_body = {"r": [], "l": []}
    bal_neurons = ("BAL_PF", "BAL_DF")

    u = net.make_inputs()
    stance_ext = np.array([muscles[a].groups[0] in EXTENSOR_STANCE_GROUPS
                           for a in acts])
    # v5 phase-reset signal groups (primary OR secondary membership: the
    # hamstrings' hip_ext arm and rect_fem's hip_flex arm count)
    hipext_idx, hipflex_idx = {}, {}
    for s in net.sides:
        hipext_idx[s] = np.array([i for i, a in enumerate(acts)
                                  if muscles[a].side == s
                                  and "hip_ext" in muscles[a].groups], dtype=int)
        hipflex_idx[s] = np.array([i for i, a in enumerate(acts)
                                   if muscles[a].side == s
                                   and "hip_flex" in muscles[a].groups], dtype=int)
    # v11 mechanosensory stance feedback: per-foot heel/toe contact normal
    # forces (N), computed from MuJoCo contacts each step
    foot_regions = ("calcn", "toes")
    BW = float(model.body_mass.sum()) * 9.81
    bodyid_region = {}
    for b in range(model.nbody):
        nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, b)
        if nm:
            for reg in foot_regions:
                if nm == f"{reg}_r":
                    bodyid_region[b] = reg + "_r"
                    foot_body["r"].append(b)
                elif nm == f"{reg}_l":
                    bodyid_region[b] = reg + "_l"
                    foot_body["l"].append(b)
    con_force = np.zeros(6)  # contact force buffer (mj_contactForce)
    # contact-onset kick state (G["contact_onset"] > 0 only; 2026-09-20):
    # per-side decaying transient set at the loading/unloading EDGES of
    # each foot's contact signal. All gating is runner-side - the network
    # topology is untouched, gain 0 changes no numbers.
    onset_prev_loaded = {"r": False, "l": False}
    onset_kick = {"r": 0.0, "l": 0.0}
    onset_decay = float(np.exp(-DT / 0.25))  # 0.25 s transient
    # crossed swing-trigger state (G["contra_swing"] > 0 only;
    # 2026-09-21): contra_prev[a] = was side a loaded last step;
    # contra_kick[b] = decaying swing trigger applied to side b's HEEL
    # port when the OPPOSITE foot a loads
    contra_prev = {"r": False, "l": False}
    contra_kick = {"r": 0.0, "l": 0.0}
    # per-side contact-reset phase machine (G["pm_gain"] > 0 only;
    # 2026-09-21): phase per leg, reset at that foot's loading onset,
    # antiphase-coupled; gates extensor/flexor ctrl in the swing window
    pm_on = (G["pm_gain"] > 0.0 or G["pm_ws"] > 0.0) and not no_ground
    pm_phase = {"r": 0.0, "l": 0.5}
    pm_prev_loaded = {"r": False, "l": False}
    # goal2 proprioceptive balance (2026-09-23): stance-gated boost on the
    # II length loop. Gate flag (not a multiplier-at-0) so the OFF path is
    # the IDENTICAL code path = bit-identical.
    vest_prop_on = G["vest_prop"] > 0.0
    # unilateral deafferentation set (AARL_DEAFF=l|r|l,r; 2026-09-21):
    # muscle actuators of these sides get ZERO afferent drive; their
    # heel/toe/load channels are also zeroed in the stance_fb block
    import os as _os2
    _deaff = set(_os2.environ.get("AARL_DEAFF", "").split(",")) - {""}
    deaff_acts = {a for a in acts if muscles[a].side in _deaff}
    for k in range(nsteps):
        t = k * DT
        drive, posture = drive_posture(t, walk_drive)
        # network state from the PREVIOUS step (or zeros initially) for
        # phase gating in deafferented runs
        v = net.compiled.V
        stance = {s: float(np.clip(v[net.idx[f"RG_E_{s}"]] / E_HI, 0, 1))
                  for s in net.sides}
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
                # 2026-09-21 unilateral deafferentation (AARL_DEAFF=l|r,
                # Ben's architecture question: is the stepping side
                # dependent on the frozen side's afferents, or vice
                # versa?): zero ALL afferent drive from the listed side.
                if a in deaff_acts:
                    u[iport[f"Ia_{a}"]] = 0.0
                    u[iport[f"II_{a}"]] = 0.0
                    u[iport[f"Ib_{a}"]] = 0.0
                    continue
                u[iport[f"Ia_{a}"]] = max(g_ia[i] * vel_norm[i], 0.0)
                # II baseline is phase-gated with its gain: an ungated
                # i0_ii=1.0 nA tonic excites every extensor MN ~0.2 ctrl
                # through swing (phase-aligned means: quads 0.25 in swing,
                # knee parked at +10 deg - diag_phase.py)
                u_ii = AFF["i0_ii"] * (0.3 + 0.7 * np.array(
                    [stance[muscles[a].side] for a in acts]))[i] \
                    + g_ii[i] * len_norm[i]
                if vest_prop_on:
                    # goal2 (SCONE T3a KL analog): scale the stance-gated
                    # II length component (the proprioceptive-balance loop)
                    u_ii = AFF["i0_ii"] * (0.3 + 0.7 * stance[s]) \
                        + g_ii[i] * len_norm[i] \
                        * (1.0 + G["vest_prop"] * stance[s])
                u[iport[f"II_{a}"]] = max(u_ii, 0.0)
                u[iport[f"Ib_{a}"]] = max(0.5 * g_ib[i] * force_norm[i], 0.0)
            # (PRESET hip-signal ports removed with the pathway 2026-09-16;
            # sensory phase reset now travels the per-muscle afferent ->
            # central edges inside the network - nothing to inject here)
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
        # --- lateral balance: stance-gated abductor strategy ---
        # Abductors pull the CoG toward their OWN stance foot: right
        # abductors pull CoG right (mujoco -y = opensim +z, Ben's
        # convention), left abductors pull CoG left (+y). So the CoG being
        # on the +y side demands RIGHT abductor activity (and vice versa),
        # gated by that side's stance (RG-E) so they work under load.
        # dy_tot = desired CoG displacement (+ = toward +y), incl. damping
        dy_tot = (y_ref - com[1]) - BAL["kv_lat"] * data.qvel[1]
        fade_lat = max(fade, 0.5)   # frontal control stays active in gait
        for side in net.sides:
            pull = max(0.0, -dy_tot) if side == "r" else max(0.0, dy_tot)
            u[iport[f"BAL_LAT_{side.upper()}"]] = float(np.clip(
                BAL["ky_lat"] * pull, 0.0, BAL["max_current"])) * \
                stance[side] * fade_lat
        # --- IMU trunk controller (vestibular surrogate, Ben 2026-09-10):
        # torso up-vector sagittal lean (+ = backward limbo lean) with
        # rate damping -> abdominals pull forward / erectors pull back.
        # Body frames are OpenSim-aligned (local +Y = up), so up = R[:,1].
        if torso_id >= 0:
            up = data.xmat[torso_id].reshape(3, 3)[:, 1]
            lean = float(np.arctan2(-up[0], max(up[2], 1e-6)))
            lean_dot = (lean - lean_prev) / DT
            lean_prev = lean
            e = lean - BAL["trk_ref"]
            u[iport["BAL_TRK_FLX"]] = float(np.clip(
                BAL["kp_trk"] * e + BAL["kd_trk"] * lean_dot, 0.0,
                BAL["max_trk"]))
            u[iport["BAL_TRK_EXT"]] = float(np.clip(
                -BAL["kp_trk"] * e - BAL["kd_trk"] * lean_dot, 0.0,
                BAL["max_trk"]))
            # --- goal2 vestibular-analog cells (2026-09-23; net.vest =
            # True only when a vest gain > 0). Signal = RECTIFIED tilt
            # deviation + rate (otolith/canal analog; the directional
            # correction stays with BAL_TRK/BAL_PF/BAL_DF above). SCONE
            # Tutorial-3a model: vestibular BodyPointReflex = torso-point
            # PD with 0.1 s delay driving all major muscles; Di Russo
            # 2023 JNE eq (5) uses the same kp*(theta-theta0)+kv*theta_dot
            # trunk-lean PD. Magnitudes reuse the Ben-tuned BAL_TRK
            # gains + clamp; the VEST cells' 0.1 s membrane tau supplies
            # the vestibular lag. Bilateral symmetric (SCONE
            # symmetric = 1). Tone-only: rises with ANY deviation.
            if net.vest:
                u_vest = float(np.clip(
                    BAL["kp_trk"] * abs(e)
                    + BAL["kd_trk"] * abs(lean_dot),
                    0.0, BAL["max_trk"]))
                for s2 in net.sides:
                    u[iport[f"VEST_c_{s2}"]] = u_vest
        u[iport["DRIVE"]] = drive
        u[iport["POSTURE"]] = posture
        # --- v11 heel/toe contact mechanosensors (audit P1a) ---
        # 2026-09-20: forces computed EVERY step in ground mode (the eval
        # objective and figures consume them); sensor ports fed only when
        # stance_fb built them.
        heel_n = toe_n = None
        if (not no_ground) or net.stance_fb:
            heel_n = {"r": 0.0, "l": 0.0}
            toe_n = {"r": 0.0, "l": 0.0}
            for ci in range(data.ncon):
                c = data.contact[ci]
                for g in (c.geom1, c.geom2):
                    b = model.geom_bodyid[g]
                    reg = bodyid_region.get(b)
                    if reg:
                        mujoco.mj_contactForce(model, data, ci, con_force)
                        if reg.endswith("_r"):
                            heel_n["r"] += max(con_force[0], 0.0)
                            toe_n["r"] += max(con_force[0], 0.0)
                        else:
                            heel_n["l"] += max(con_force[0], 0.0)
                            toe_n["l"] += max(con_force[0], 0.0)
                        break
            if not no_ground:
                log_contact[k, 0] = heel_n["r"] + toe_n["r"]
                log_contact[k, 1] = heel_n["l"] + toe_n["l"]
                for si, s in enumerate(("r", "l")):
                    if foot_body[s]:
                        log_footz[k, si] = min(
                            data.xpos[b, 2] for b in foot_body[s])
        if net.stance_fb:
            heel_sig = {"r": min(heel_n["r"] / (0.35 * BW), 1.5),
                        "l": min(heel_n["l"] / (0.35 * BW), 1.5)}
            toe_sig = {"r": min(toe_n["r"] / (0.50 * BW), 1.5),
                       "l": min(toe_n["l"] / (0.50 * BW), 1.5)}
            load_sig = {"r": min((heel_n["r"] + toe_n["r"])
                                 / (0.60 * BW), 1.5),
                        "l": min((heel_n["l"] + toe_n["l"])
                                 / (0.60 * BW), 1.5)}
            for s in net.sides:
                u[iport[f"HEEL_c_{s}"]] = G["heel_rge"] * heel_sig[s]
                u[iport[f"TOE_c_{s}"]] = G["toe_rge"] * toe_sig[s]
                u[iport[f"LOAD_c_{s}"]] = G["ib_rge"] * load_sig[s]
            # contact-EVENT transients (G["contact_onset"] > 0 only):
            # heel strike (load edge) = brief +kick on HEEL_c/TOE_c (E
            # trigger, F suppression); toe-off (unload edge) = brief
            # -kick (E release, F disinhibition). Replaces the tonic
            # contact level as the rhythm carrier under load.
            if G["contact_onset"] > 0.0:
                for s in net.sides:
                    loaded = load_sig[s] > 0.05
                    if loaded and not onset_prev_loaded[s]:
                        onset_kick[s] = 1.0
                    elif not loaded and onset_prev_loaded[s]:
                        onset_kick[s] = -1.0
                    else:
                        onset_kick[s] *= onset_decay
                    onset_prev_loaded[s] = loaded
                    kick = G["contact_onset"] * onset_kick[s]
                    u[iport[f"HEEL_c_{s}"]] += kick
                    u[iport[f"TOE_c_{s}"]] += kick
            # crossed swing trigger (G["contra_swing"] > 0 only): when
            # side A's foot loads (heel strike), the CONTRALATERAL side
            # B gets a decaying NEGATIVE kick on its HEEL port -> less
            # heel_in_B excitation of RG-E_B and less inhibition of
            # RG-F_B = reset-to-swing for the stance leg. Breaks the
            # frozen-stance latch (the frozen foot's own unloading edge
            # never fires; the opposite foot's loading edge does).
            if G["contra_swing"] > 0.0:
                for s in net.sides:
                    o = "l" if s == "r" else "r"
                    loaded_a = load_sig[s] > 0.05
                    if loaded_a and not contra_prev[s]:
                        contra_kick[o] = -1.0
                    else:
                        contra_kick[o] *= onset_decay
                    contra_prev[s] = loaded_a
                    u[iport[f"HEEL_c_{o}"]] += G["contra_swing"] * \
                        contra_kick[o]
        # per-side contact-reset phase machine (G["pm_gain"]/pm_ws > 0;
        # needs the contact forces, computed above in ground mode).
        # v2: the swing window is LOAD-GATED - the phase holds at 0.55
        # until the contralateral foot carries >= 25% BW, so a leg only
        # "opens" its swing once the other leg actually bears weight.
        if pm_on and heel_n is not None:
            for s in ("r", "l"):
                o = "l" if s == "r" else "r"
                loaded = (heel_n[s] + toe_n[s]) > 20.0
                other_load = heel_n[o] + toe_n[o]
                if loaded and not pm_prev_loaded[s]:
                    pm_phase[s] = 0.0      # heel-strike phase reset
                elif (G["pm_ws"] > 0.0
                      and 0.55 <= pm_phase[s] < 0.62
                      and other_load < 0.15 * BW):
                    pm_phase[s] = 0.55     # wait at the gate
                else:
                    pm_phase[s] = (pm_phase[s]
                                   + DT / max(G["pm_T"], 0.3)) % 1.0
                pm_prev_loaded[s] = loaded
            # gentle antiphase pull (Di Russo eq-7 coupling)
            err = ((pm_phase["l"] - pm_phase["r"]) % 1.0) - 0.5
            pm_phase["l"] = (pm_phase["l"] - 0.8 * DT * err) % 1.0
        # --- v11b semi-closed loops: per-side afferent relay currents ---
        if net.aff_loops:
            ga, gb = PHASE_RESET.get("stance_gate", (0.3, 0.7))
            for s in net.sides:
                ie = hipext_idx[s]
                ifl = hipflex_idx[s]
                # extensor afferent drive (force + length components)
                ext_drive = 0.0
                if ie.size:
                    ext_force = float(np.mean(
                        [min(F[aid[a]] / max(Fmax[aid[a]], 1e-9), 1.0)
                         for a in np.array(acts)[ie]]))
                    ext_len = float(np.mean(len_norm[ie]))
                    ext_drive = max(ext_force, 0.0) + max(ext_len, 0.0) * 0.5
                    ext_drive *= (ga + gb * stance[s])
                u[iport[f"AFF_E_{s}"]] = G["aff_e_rg"] * max(ext_drive, 0.0)
                # flexor afferent drive (velocity + length components)
                flex_drive = 0.0
                if ifl.size:
                    flex_vel = float(np.mean(vel_norm[ifl]))
                    flex_len = float(np.mean(len_norm[ifl]))
                    flex_drive = max(-flex_vel, 0.0) + max(-flex_len, 0.0) * 0.5
                u[iport[f"AFF_F_{s}"]] = G["aff_f_rg"] * max(flex_drive, 0.0)
        # pm v4 afferent disfacilitation (G["pm_aff"] > 0): during a
        # side's swing window, scale its load-afferent input channels
        # by (1 - pm_aff*w) - the re-latch loop (load -> heel/toe/Ib ->
        # RG-E -> MN) must be broken for the swing to actually happen.
        if pm_on and G["pm_aff"] > 0.0:
            for s in net.sides:
                w = _pm_window(pm_phase[s])
                if w > 0.0:
                    f = 1.0 - G["pm_aff"] * w
                    for ch in (f"HEEL_c_{s}", f"TOE_c_{s}",
                               f"LOAD_c_{s}", f"AFF_E_{s}",
                               f"AFF_F_{s}"):
                        if ch in iport:
                            u[iport[ch]] *= f
        # unilateral deafferentation: zero the side's mechanosensor and
        # afferent-relay channels AFTER all injections (wins over kicks)
        if _deaff:
            for s in net.sides:
                if s in _deaff:
                    for ch in (f"HEEL_c_{s}", f"TOE_c_{s}",
                               f"LOAD_c_{s}", f"AFF_E_{s}",
                               f"AFF_F_{s}"):
                        if ch in iport:
                            u[iport[ch]] = 0.0
        # solved standing pattern -> motoneuron posture bias (fades to a
        # small floor during walking so the stepping pattern can take over;
        # a 0.2 floor kept knee-extensor tone high enough to pin the knees
        # straight in the air-stepping test)
        stand_frac = 0.05 + 0.95 * (1.0 - min(drive / max(walk_drive, 1e-6), 1.0))
        # ankle standing-tone trim: the soleus/tib_post posture bias held
        # a ~-45 deg PF ankle set-point through gait (v9 overlay tiptoe
        # offset) - real soleus tonic EMG drops with locomotor drive, so
        # scale that group's POST bias toward ankle_post_walk_trim as
        # drive rises (1.0 = v9 behavior, unchanged)
        walk_frac = min(drive / max(walk_drive, 1e-6), 1.0)
        trim = 1.0 - (1.0 - G["ankle_post_walk_trim"]) * walk_frac
        for i, a in enumerate(acts):
            sf = stand_frac
            if muscles[a].groups[0] == "ankle_pf":
                sf *= trim
            u[iport[f"POST_{a}"]] = 6.0 * act_stand[i] * sf

        # --- neural step + motor mapping ---
        v = net.step(u)
        for i, a in enumerate(acts):
            c = float(np.clip(v[net.idx[net.mn_names[a]]] / E_HI,
                              0.0, 1.0))
            # per-side contact-reset phase machine (G["pm_gain"] > 0
            # only): during this side's swing window, scale extensor
            # ctrl DOWN and flexor ctrl UP — guarantees every leg a
            # swing phase each cycle (the frozen-stance fix)
            if pm_on and heel_n is not None:
                w = _pm_window(pm_phase[a[-1]])
                if w > 0.0:
                    g0 = muscles[a].groups[0]
                    if g0 in ("knee_ext", "ankle_pf"):
                        c *= (1.0 - G["pm_gain"] * w)
                    elif g0 in ("knee_flex", "ankle_df", "hip_flex"):
                        # v3 (2026-09-21): ADD the burst - the old
                        # multiplicative boost on a ~0 ctrl did nothing
                        # (MNs below threshold on the unloaded side)
                        if G["pm_add"] > 0.0:
                            c = min(1.0, c + G["pm_add"] * w)
                        else:
                            c = min(1.0,
                                    c * (1.0 + 0.6 * G["pm_gain"] * w))
                # weight-shift prep (pm_ws): during the OTHER side's
                # stance-prep ramp, the upcoming-swing side's abductors
                # scale DOWN and the upcoming-stance side's scale UP
                if G["pm_ws"] > 0.0 and muscles[a].groups[0] == "hip_abd":
                    o = "l" if a[-1] == "r" else "r"
                    wp = min(max((pm_phase[o] - 0.42) / 0.20, 0.0), 1.0)
                    if wp > 0.0 and 0.42 <= pm_phase[o] < 0.62:
                        if a[-1] == o:
                            c *= (1.0 - G["pm_ws"] * wp)
                        else:
                            c = min(1.0, c * (1.0 + G["pm_ws"] * wp))
            data.ctrl[aid[a]] = c

        if t < WARMUP:
            # hold the keyframe pose while activations build (qpos untouched,
            # velocities zeroed, activation states integrate via mj_forward)
            data.qvel[:] = 0.0
            mujoco.mj_forward(model, data)
        else:
            if adaptive_tol > 0.0:
                _adaptive_step(model, data, adaptive_tol)
            else:
                mujoco.mj_step(model, data)
        if k % 10 == 0:
            if viewer is not None:
                if not viewer.is_running():
                    print("viewer closed - ending run early")
                    break
                viewer.sync()
            # --realtime: hold the loop to the wall clock (1x playback);
            # default is full speed ("just run")
            if realtime and (viewer is not None or scope is not None):
                delay = t_wall0 + (k + 1) * DT - time.perf_counter()
                if delay > 0:
                    time.sleep(delay)

        # ---- instability forensics: report precursors and the first blowup
        if not np.all(np.isfinite(data.qacc)) or not np.all(np.isfinite(data.qvel)):
            bad = np.flatnonzero(
                ~np.isfinite(data.qacc) | ~np.isfinite(data.qvel))
            print(f"\n!! qacc/qvel non-finite at t={t:.3f}, dofs: {bad[:12]}")
            for b in bad[:12]:
                for j in range(model.njnt):
                    da = int(model.jnt_dofadr[j])
                    if da == b:
                        print(f"   dof {b}: joint "
                              f"{mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)}")
                        break
            print(f"   ncon={data.ncon}")
            break
        ma = float(np.max(np.abs(data.qacc)))
        mv = float(np.max(np.abs(data.qvel)))
        if not eval_mode and (mv > 30.0 or ma > 1.0e5):
            j = int(np.argmax(np.abs(data.qvel)))
            print(f"t={t:6.2f} precursor: max|qacc|={ma:.3g} "
                  f"max|qvel|={mv:.3g} (dof {j}), ncon={data.ncon}, "
                  f"max ctrl={float(np.max(data.ctrl)):.2f}")
        if not eval_mode and t >= WARMUP and k % int(0.5 / DT) == 0:
            print(f"t={t:6.2f} max|qacc|={ma:10.1f} max|qvel|={mv:7.2f} "
                  f"ncon={data.ncon:3d} com_z={float(com[2]):.3f} "
                  f"com_x={float(com[0]):+.3f}")

        # --- log ---
        log_t[k] = t
        log_com[k] = com
        log_qfull[k] = data.qpos
        for i, a in enumerate(KEY_ACTS):
            log_act[k, i] = data.ctrl[aid[a]]
        for i, jn in enumerate(KEY_JOINTS):
            log_q[k, i] = data.qpos[jadr[jn]]
        log_neuro[k] = (drive, posture,
                        v[net.idx["RG_E_r"]], v[net.idx["RG_F_r"]],
                        v[net.idx["RG_E_l"]], v[net.idx["RG_F_l"]],
                        *(v[net.idx[c]] for c in _pf_watch),
                        u[iport["BAL_TRK_FLX"]], u[iport["BAL_TRK_EXT"]],
                        u[iport["BAL_LAT_R"]], u[iport["BAL_LAT_L"]])
        if scope is not None:
            scope.update(t, log_neuro[k], log_act[k, :8])

    # atomic-replace save (2026-09-18): stage-3 trial 0 died with
    # OSError 22 opening spinal_run.npz - transient Windows lock on the
    # just-replaced file (same family as the preview-lock problem the
    # figure renderer already solves). Write beside + os.replace retry.
    import time as _time
    # 2026-09-20: AARL_NPZ redirects the run output (default
    # spinal_run.npz) - external locks (Spyder indexer/AV) on the default
    # name killed WinError-5 loops; point batch scripts at a side name.
    _npz_name = _os.environ.get("AARL_NPZ", "spinal_run.npz")
    _tmp = HERE / f"{_npz_name}.{_os.getpid()}.tmp.npz"
    with open(_tmp, "wb") as _fh:
        np.savez_compressed(_fh, t=log_t, act=log_act,
                            q=np.degrees(log_q), qfull=log_qfull,
                            com=log_com, neuro=log_neuro,
                            contact=log_contact, footz=log_footz,
                            key_acts=KEY_ACTS, key_joints=KEY_JOINTS,
                            neuro_names=NEURO_NAMES,
                            cfg=("ground" if not no_ground else "air"))
    _err = None
    for _ in range(40):   # 2026-09-20: 20x0.15s was not enough - WinError 5
        # (AV/indexer lock) killed a rescore run; 40x0.5s = 20 s backoff
        try:
            _os.replace(_tmp, HERE / _npz_name)
            _err = None
            break
        except OSError as _e:
            _err = _e
            _time.sleep(0.5)
    if _err is not None:
        raise _err
    if viewer is not None:
        viewer.close()
    if scope is not None:
        scope.close()

    if eval_mode:
        # objective metrics for optuna_walk.py. Rows after an instability
        # break are zero-filled, so slice to the valid prefix.
        n_done = int(np.sum(log_t > 0)) + 1
        n_done = min(n_done, len(log_t))
        tt = log_t[:n_done]
        i_ws = int(np.searchsorted(tt, SCHEDULE["walk"][0]))
        finite_q = bool(np.all(np.isfinite(log_q[:n_done])))
        finite_c = bool(np.all(np.isfinite(log_com[:n_done])))
        def _mn(a):
            return float(np.nanmin(a)) if finite_q else 0.0
        def _mx(a):
            return float(np.nanmax(a)) if finite_q else 0.0
        walk_slice = log_neuro[i_ws:n_done, 2]
        try:
            import kine_ref
            k = kine_ref.compare(log_t[:n_done],
                                 np.degrees(log_q[:n_done]),
                                 log_neuro[:n_done], SCHEDULE["walk"][0],
                                 ref=kine_ref.ref_cached(),
                                 contact=log_contact[:n_done])
        except Exception:
            k = None
        # ---- goal2 balance metrics (2026-09-23): always computed over
        # the post-warmup prefix; consumed by the stage-4 standing
        # objective (COM sway radius, tilt envelope, contact symmetry,
        # no falls). Zero-length guard for degenerate runs.
        i_bal = min(int(np.searchsorted(tt, WARMUP)), n_done - 1)
        seg_c = log_com[i_bal:n_done]
        seg_qd = np.degrees(log_q[i_bal:n_done])
        seg_ct = log_contact[i_bal:n_done]
        if seg_c.shape[0]:
            _rad = np.hypot(seg_c[:, 0] - x_ref, seg_c[:, 1] - y_ref)
            bal_sway = float(np.max(_rad))
            bal_sway_rms = float(np.sqrt(np.mean(_rad ** 2)))
            bal_tilt_max = float(np.max(np.abs(seg_qd[:, 0])))
            bal_com_z_min = float(np.min(seg_c[:, 2]))
        else:
            bal_sway = bal_sway_rms = bal_tilt_max = bal_com_z_min = 0.0
        _cr = float(np.mean(seg_ct[:, 0])) if seg_ct.shape[0] else 0.0
        _cl = float(np.mean(seg_ct[:, 1])) if seg_ct.shape[0] else 0.0
        bal_contact_sym = float(_cr / (_cr + _cl)) if (_cr + _cl) > 1e-9 \
            else 0.25   # no contact at all: score between 0 (one-foot) and 0.5 (even)
        metrics = dict(
            nan=not (finite_q and finite_c),
            t_end=float(tt[-1]),
            dx=float(log_com[n_done - 1, 0] - log_com[i_ws, 0])
            if n_done > i_ws else 0.0,
            kz=float(np.min(log_com[:n_done, 2])) if finite_c else 0.0,
            tilt_max=_mx(np.degrees(log_q[:n_done, 0])),
            knee_min=_mn(log_q[:n_done, 4]),
            hip_amp=_mx(log_q[:n_done, 3]) - _mn(log_q[:n_done, 3]),
            duty=float(np.mean(walk_slice > 0.5 * max(np.max(walk_slice),
                                                      1e-6)))
            if walk_slice.size else 0.0,
            kine=k,
            kine_score=(k["kine_score"] if k else -25.0),
            burst_r=int(np.sum(np.diff(
                (log_neuro[:n_done, 2]
                 > 0.5 * max(np.max(log_neuro[:n_done, 2]), 1e-6)
                 ).astype(int)) == 1)),
            # goal2 balance fields (stage-4 objective consumes these)
            bal_sway=bal_sway,
            bal_sway_rms=bal_sway_rms,
            bal_tilt_max=bal_tilt_max,
            bal_contact_sym=bal_contact_sym,
            bal_com_z_min=bal_com_z_min,
            bal_fell=bool(bal_com_z_min < 0.55),
        )
        return metrics

    # ---------------- honest summary ----------------
    pelz = log_com[:, 2]
    xcom = log_com[:, 0]
    fall = pelz.min() < 0.55
    ifall = np.argmax(pelz < 0.55) if fall else -1
    walk_mask = ((log_t >= SCHEDULE["walk"][0]) &
                 (log_t <= SCHEDULE["walk"][1]))
    neuro_col = {name: i for i, name in enumerate(NEURO_NAMES)}
    n_steps = int(np.sum(np.diff(
        (log_neuro[walk_mask, neuro_col["RG_F_r"]] > 0.8).astype(int)) == 1))
    print(f"\n== summary ==")
    print(f"pelvis/COM height: start {pelz[0]:.2f} m, min {pelz.min():.2f} m, "
          f"final {pelz[-1]:.2f} m  ({'FELL at t=%.1f s' % log_t[ifall] if fall else 'stayed up'})")
    print(f"COM x: min {xcom.min():.2f} max {xcom.max():.2f} final {xcom[-1]:.2f} m")
    print(f"COM y: min {log_com[:, 1].min():+.3f} max {log_com[:, 1].max():+.3f} m")
    print(f"RG_F_r bursts during walk window: {n_steps}")
    # per-leg rhythm metrics during the walk window (cycle period, E duty)
    if walk_mask.sum() > 100:
        for leg in ("r", "l"):
            e_col = neuro_col[f"RG_E_{leg}"]
            f_col = neuro_col[f"RG_F_{leg}"]
            rge = log_neuro[walk_mask, e_col]
            rgf = log_neuro[walk_mask, f_col]
            on = rge > 0.5 * max(rge.max(), 1e-6)
            rises = np.flatnonzero(np.diff(on.astype(int)) == 1)
            falls_ = np.flatnonzero(np.diff(on.astype(int)) == -1)
            if len(rises) >= 2:
                cyc = float(np.mean(np.diff(log_t[walk_mask][rises])))
                duty = float(on.mean())
                print(f"RG_{leg}: cycles in window {len(rises)}, period "
                      f"{cyc:.2f} s ({1.0 / max(cyc, 1e-9):.2f} Hz), "
                      f"E-duty {duty:.2f}")
    print(f"peak ctrl: {log_act.max(axis=0).round(2)} ({KEY_ACTS})")
    for side in ("r", "l"):
        cols = [i for i, j in enumerate(KEY_JOINTS) if j.endswith(f"_{side}")]
        for i in cols:
            print(f"  {KEY_JOINTS[i]:18s} {np.degrees(log_q[:, i]).min():7.1f} .. "
                  f"{np.degrees(log_q[:, i]).max():7.1f} deg")
    for i, jn in enumerate(KEY_JOINTS[:3]):
        print(f"  {jn:18s} {np.degrees(log_q[:, i]).min():7.1f} .. "
              f"{np.degrees(log_q[:, i]).max():7.1f} deg")
    print("saved spinal_run.npz")

    try:
        import matplotlib
        if scope is None:
            matplotlib.use("Agg")  # headless: no display needed for the png
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
        for i, a in enumerate(KEY_ACTS):
            ax[0].plot(log_t, log_act[:, i], lw=0.8,
                       label=a if i < 8 or a.endswith("_l") else None)
        ax[0].set_ylabel("activation")
        for i, jn in enumerate(KEY_JOINTS):
            ax[1].plot(log_t, np.degrees(log_q[:, i]), lw=0.8,
                       label=jn if i < 3 else None)
        ax[1].set_ylabel("joint angle [deg]")
        ax[1].legend(fontsize=7, ncol=3)
        ax[2].plot(log_t, log_com[:, 0], label="com x")
        ax[2].plot(log_t, log_com[:, 2], label="com z")
        ax[2].set_ylabel("COM [m]"); ax[2].legend(fontsize=7)
        plot_neuro = ("DRIVE", "POSTURE", "RG_E_r", "RG_F_r",
                      "PF_E2_r", "PF_F1_r", "RG_E_l", "RG_F_l")
        for nm in plot_neuro:
            ax[3].plot(log_t, log_neuro[:, neuro_col[nm]], lw=0.8,
                       label=nm)
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
