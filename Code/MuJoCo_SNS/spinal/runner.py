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
                     [--realtime]
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
from params import AFF, BAL, DT, E_HI, MOD, SCHEDULE

import mujoco

from scipy.optimize import nnls

HERE = Path(__file__).parent
MODEL = Path(r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
             r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml")

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
              boundmass=0.1, pp_contact=False, pin_rot=False):
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
                  boundmass=0.1, pp_contact=False, pin_rot=False):
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
                         pin_rot=pin_rot)
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
    args = list(argv)
    if "--best" in args:
        # load the optimizer's winning parameters (optuna_walk.py output)
        # before anything reads the params module
        import json as _json
        with open(HERE / "best_walk_params.json", encoding="utf-8") as f:
            best = _json.load(f)["params"]
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
        _p.BAL["kx"] = best["kx"]
        args.remove("--best")
        walk_drive = float(best.get("drive", walk_drive))
        print(f"loaded best_walk_params.json (drive={walk_drive:.2f})",
              flush=True)
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
        elif a == "--no-ground":
            no_ground = True
        elif a == "--no-interleg":
            interleg = False
        elif a == "--eval":
            # fast objective evaluation for the optimizer (optuna_walk.py):
            # 12 s schedule (1.5 stand / 1.5 ramp / 7.5 walk / 1.5 wind-up),
            # quiet, no npz/png writes; main() returns a metrics dict
            eval_mode = True
            SCHEDULE.update(stand1=(0.0, 1.5), ramp_up=(1.5, 3.0),
                            walk=(3.0, 10.5), ramp_down=(10.5, 11.0),
                            stand2=(11.0, 12.0))
        elif a == "--leg-damping":
            leg_damping = float(args.pop(0))
        elif a == "--view":
            view = True
        elif a == "--scope":
            scope_on = True
        elif a == "--realtime":
            realtime = True
        elif a == "--time":
            SCHEDULE["stand2"] = (SCHEDULE["ramp_down"][1],
                                  SCHEDULE["ramp_down"][1] + float(args.pop(0)))

    model = mujoco.MjModel.from_xml_path(str(MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    key_pose = capture_pose(model)

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
        model = apply_harness(model, data, kxy=harness,
                              no_ground=no_ground, leg_damping=leg_damping,
                              pin_rot=no_ground)
    elif no_ground or leg_damping is not None:
        # no rig, but XML patches needed (suspended-air / damping sweep):
        # zero-stiffness springs with negligible damping = free pelvis
        model = apply_harness(model, data, kxy=0, kz=0, ky=0,
                              dxy=1.0, dy=1.0,
                              no_ground=no_ground, leg_damping=leg_damping,
                              pin_rot=no_ground)
    data = mujoco.MjData(model)
    seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)

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
    # RG, all four PF groups (right), balance cells
    NEURO_NAMES = ("DRIVE", "POSTURE",
                   "RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
                   "PF_E1_r", "PF_E2_r", "PF_F1_r", "PF_F2_r",
                   "BAL_TRK_FLX", "BAL_TRK_EXT", "BAL_LAT_R", "BAL_LAT_L")
    log_neuro = np.zeros((nsteps, len(NEURO_NAMES)))
    bal_neurons = ("BAL_PF", "BAL_DF")

    u = net.make_inputs()
    stance_ext = np.array([muscles[a].groups[0] in EXTENSOR_STANCE_GROUPS
                           for a in acts])
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
                u[iport[f"Ia_{a}"]] = max(g_ia[i] * vel_norm[i], 0.0)
                # II baseline is phase-gated with its gain: an ungated
                # i0_ii=1.0 nA tonic excites every extensor MN ~0.2 ctrl
                # through swing (phase-aligned means: quads 0.25 in swing,
                # knee parked at +10 deg - diag_phase.py)
                u[iport[f"II_{a}"]] = max(
                    AFF["i0_ii"] * (0.3 + 0.7 * np.array(
                        [stance[muscles[a].side] for a in acts]))[i]
                    + g_ii[i] * len_norm[i], 0.0)
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
        u[iport["DRIVE"]] = drive
        u[iport["POSTURE"]] = posture
        # solved standing pattern -> motoneuron posture bias (fades to a
        # small floor during walking so the stepping pattern can take over;
        # a 0.2 floor kept knee-extensor tone high enough to pin the knees
        # straight in the air-stepping test)
        stand_frac = 0.05 + 0.95 * (1.0 - min(drive / max(walk_drive, 1e-6), 1.0))
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
                        v[net.idx["PF_E1_r"]], v[net.idx["PF_E2_r"]],
                        v[net.idx["PF_F1_r"]], v[net.idx["PF_F2_r"]],
                        u[iport["BAL_TRK_FLX"]], u[iport["BAL_TRK_EXT"]],
                        u[iport["BAL_LAT_R"]], u[iport["BAL_LAT_L"]])
        if scope is not None:
            scope.update(t, log_neuro[k], log_act[k, :8])

    np.savez_compressed(HERE / "spinal_run.npz", t=log_t, act=log_act,
                        q=np.degrees(log_q), qfull=log_qfull, com=log_com,
                        neuro=log_neuro,
                        key_acts=KEY_ACTS, key_joints=KEY_JOINTS,
                        neuro_names=NEURO_NAMES,
                        cfg=("ground" if not no_ground else "air"))
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
        metrics = dict(
            nan=not (finite_q and finite_c),
            t_end=float(tt[-1]),
            dx=float(log_com[n_done - 1, 0] - log_com[i_ws, 0]),
            kz=float(np.min(log_com[:n_done, 2])) if finite_c else 0.0,
            tilt_max=_mx(np.degrees(log_q[:n_done, 0])),
            knee_min=_mn(log_q[:n_done, 4]),
            hip_amp=_mx(log_q[:n_done, 3]) - _mn(log_q[:n_done, 3]),
            burst_r=int(np.sum(np.diff(
                (log_neuro[:n_done, 2]
                 > 0.5 * max(np.max(log_neuro[:n_done, 2]), 1e-6)
                 ).astype(int)) == 1)),
        )
        return metrics

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
    print(f"COM y: min {log_com[:, 1].min():+.3f} max {log_com[:, 1].max():+.3f} m")
    print(f"RG_F_r bursts during walk window: {n_steps}")
    # per-leg rhythm metrics during the walk window (cycle period, E duty)
    if walk_mask.sum() > 100:
        for leg, (e_col, f_col) in (("r", (2, 3)), ("l", (6, 7))):
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
