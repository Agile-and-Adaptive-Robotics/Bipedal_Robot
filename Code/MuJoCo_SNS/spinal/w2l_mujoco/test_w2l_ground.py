r"""M6 GATE — GROUND STANDING BALANCE + CONTACT-DRIVEN GROUND WALKING
(w2l_mujoco\test_w2l_ground.py).

Body: w2l_ground.xml (written by this script) = w2l_mjcf_fixed.xml (the
M3 axis-fixed M1 body) with a FREEJOINT added to the Root pelvis — the
milestone-6 ask ("release the pelvis"). The aproj's rest pose is kept
verbatim; it hovers 2.6-3.1 cm (M1 report section 2.8), and under the
drop protocol below the harness holds it at that AIR height.

DROP PROTOCOL (2026-10-01; Ben 2026-09-30: "In AnimatLab they can drop
it WITH the harness" — the virtual walker holds the body in AIR while
the CPG air-steps, then LOWERS it onto the platform and STAYS
ATTACHED; same protocol as test_li_stepping.py): the rest pose spawns
VERBATIM (plates hover 2.6-3.1 cm, M1 section 2.8 — the old spawn-time
--drop grounding of the free-joint z is REMOVED); a feedforward-weight
+ vertical-PD harness on the root free-joint z holds the AIR height for
--hold s (3.0), lowers to contact height (z_air - --drop - --sink) over
--lower s (1.0), then keeps feedforward at --wsup (0.7) of body weight
+ PD at contact height — the feet bear the remaining fraction so the
heel/toe contact encoders fire, and the body cannot buckle (sidesteps
the no-Kse/Kpe body-fidelity blocker). Metric windows start at
hold+lower+0.5 s.

Contact encoders (the ask: "heel/toe from MuJoCo contact forces per the
rules"): the aproj's four dedicated contact plates are the sensors —
heel S = foot S_contact geom, toe S = toe S_contact geom. Per physics
step the normal force of every contact involving a plate is accumulated
(mujoco.mj_contactForce, contact-frame z = normal) and the PORT current
is a saturating linear encoder:

    I_port[S] = ctn_amp * clip(F_plate[S] / cref, 0, 1)   [nA]

(ctn_amp default 4.0 nA = the saturating/strike case from the M5 gate-b
calibration; cref default 50 N ~ 12% body weight.) The Ib group ports
keep the M5 wiring: I = ibnA * mean normalized extensor force.

Support rig (knob --rig=S, walk-phase default now 1.0 = retained): the
HORIZONTAL leash on the COM + a TILT assist (rotation vector of
R_root @ R0^T) via xfrc_applied, stiffness K_eff = S * K0 with
K0 = {kxy 2000 N/m, krot 400 N m/rad} + critical-rate damping. The old
soft-COM-z kz spring is REPLACED by the protocol vertical channel
above (which is NOT scaled by S — the harness is the protocol piece).
S=0 = leash/tilt off only. WHY A RIG IS NEEDED AT ALL (measured this
session, tmp\m6_static_margin.py): at the grounded rest pose the COM
sits 3.09 cm OUTSIDE the front edge of the foot support polygon — the
pose cannot statically stand even rigid, and the source aproj itself
ships the Root with Freeze=True (never a free-standing model).
SETTLE DEVICE (runtime, documented): for t < settle (default 0.4 s)
the rig runs at S=1 while neural tone builds; at t=settle the anchor
is re-captured at the settled AIR pose (inside the hold), then S
blends 1.0 -> target over 0.6 s.

STAND phase (--phase=stand|--phase=both, dur_stand=12 s): drive regime
--sregime=latch (default; te=tf=0 + the 10 nA 10 ms kickoff on
Stimulus_1 only -> the M3-measured extension latch = a stance posture)
or --sregime=step (the full te=3/tf=4 stepping drive). Sway = COM x/y
std+range; tilt = pelvis rotation vs its settled orientation; a FALL =
pelvis centre below 0.60 m or tilt above 45 deg. Verdict FULL if S=0
stands, PARTIAL if only S>0 stands.

WALK phase (--phase=walk|--phase=both, dur_walk=20 s): the M5
afferented split-RG net with REAL heel/toe contact replacing the
scripted scheduler, te=3/tf=4, unrigged by default (--rig=0). Metrics
over >= 20 s: heel-strike cycles + period, duty, flexion ranges, tilt,
pelvis height, forward speed, falls, RG-E bursts. FULL HONESTY: if it
does not walk, the failure is characterized (shuffle / fall / period)
and the top blockers named with evidence.

Run (cwd w2l_mujoco\):
    C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_w2l_ground.py
        [--phase=both --dur_stand=12 --dur_walk=20 --ctn_amp=4.0
         --cref=50 --ibnA=1.0 --te=3.0 --tf=4.0 --tau=0.25 --cap=0.5
         --damp=3.0 --rig=1.0 --sregime=latch --settle=0.4
         --hold=3.0 --lower=1.0 --wsup=0.7 --sink=0.01
         --contact-mode=normal|shear|total   (2026-10-20; default normal
          = z-only plate force, bit-identical; shear/total swap what the
          encoders - and the RA/SA channels - sense)
         --ra=1 --sa=1 --ra-syn=.. --sa-syn=..   (bilateralrg2 only;
          RA/SA mechanoreceptor channels, default off = bit-identical)
         --presyn-lo=.. --presyn-hi=.. --contact-presyn-lo=..
         --contact-presyn-hi=..   (AnimatLab PreSynapticThreshold /
          SaturationLevel analogs; defaults 0/5 = bit-identical)]
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np           # noqa: E402
import mujoco                # noqa: E402

import test_w2l_air as TW    # noqa: E402  (DT_PHY, NSUB, metric helpers)
from test_w2l_air import FIXED, DT_PHY, NSUB, bursts_refractory  # noqa: E402

GROUND = os.path.join(HERE, "w2l_ground.xml")

KV = dict(phase="both", dur_stand=12.0, dur_walk=20.0,
          ctn_amp=4.0, cref=50.0, ibnA=1.0,
          te=3.0, tf=4.0, tau=0.25, cap=0.5, acap=-1.0, damp=3.0,
          rig=1.0, sregime="latch", settle=0.4, heelN=5.0, drop=0.026,
          hold=3.0, lower=1.0, wsup=0.7, sink=0.01,
          # per-side drive (2026-10-16): te_l/te_r/tf_l/tf_r = per-side
          # GAMMA (RG tonic) overrides in nA (nan = symmetric te/tf
          # baseline, bit-identical); alpha_l/alpha_r = per-side ALPHA
          # (MN-level) drive in nA (0 = ports not built, bit-identical).
          te_l=float("nan"), te_r=float("nan"),
          tf_l=float("nan"), tf_r=float("nan"),
          alpha_l=0.0, alpha_r=0.0,
          # leg-leg collision (2026-10-16 sweep): 1 (default) = the
          # shipped MJCF behavior (everything collides with everything;
          # the near-parallel tibia capsules interpenetrate at spawn and
          # MuJoCo wedges them apart with ~3.4 kN, measured
          # reports_20261016\w2l_walk_sweep\_w2l_static_pose.py). 0 =
          # exclude L<->R leg pairs ONLY (runtime contype remask; the
          # source aproj is a PLANAR walker with NO leg-leg contact —
          # this is a fidelity fix toward the source, default-off).
          legcoll=1.0,
          # yielding virtual walker (2026-10-16 sweep): tau in s for the
          # anchor-x EMA (see run_phase; 0 = fixed anchor = bit-identical
          # default). The AnimatLab protocol's WalkingPath MOVES with the
          # robot; a fixed anchor makes sustained progression impossible.
          walker_yield=0.0,
          # feedback-gain passthrough (2026-10-16 wean stage): build()
          # afferent gain overrides + commissural scale. DEFAULTS = Ben's
          # rules values (AFF_GAINS + comm 1.0), so the default build is
          # UNCHANGED by construction (bit-exact; gate re-run post-edit).
          # These are the knobs the curr_s6_drive_robot study tunes while
          # the te/tf drive scale is weaned down.
          g_heel=0.5,        # aff_heel  (heel -> InE/InF + PF-E INs)
          g_toe_in=1.0,      # aff_toe_in  (toe SN -> toe IN)
          g_toe_df=5.0,      # aff_toe_df  (toe IN -> TOEDF IN)
          g_toe_df_mn=2.749, # aff_toe_df_mn (TOEDF IN -> DF MN inh)
          g_ib=0.5,          # aff_ib    (Ib grp -> RG-E/InE/PF-E)
          comm=1.0,          # commissural gain scale (0 = ablation)
          # ---- 2026-10-20: --net=bilateralrg2 runs Ben's CORRECTED aproj
          # connectome (2026-10-18) via w2l_cpg.build_w2l_net instead of the
          # M5 aff lineage. The ch_* knobs are that net's contact-family
          # gains (W2L_GAINS defaults; ctn_amp/cref/ibnA/te/tf/tau/alpha/
          # comm all carry over). Default 'aff' = the M5 lineage,
          # bit-identical to the wean-walk of record.
          net="aff",
          ch_hc=0.08,        # contact_hc: heel -> ipsi extensor cascade
          ch_ank=1.0,        # contact_ank: heel/toe -> ipsi Ank MN ext
          ch_contra=0.08,    # contact_contra_flex: heel -> CONTRA RG flx
          # ---- 2026-10-20 (Ben's mechanoreceptor/presynaptic ask):
          # contact sensing mode. 'normal' (default) = the historical
          # mj_contactForce z-only accumulation (BIT-IDENTICAL); 'shear'
          # = tangential magnitude hypot(fx, fy) only; 'total' =
          # hypot(normal, tangential). Applies to the heel/toe encoder
          # AND the RA/SA channel feeds (same plate force).
          contact_mode="normal",
          # RA (fast-adapting) / SA (slow-adapting) mechanoreceptor
          # channels, bilateralrg2 only: >0 builds the "RA|SA heel|toe
          # {L,R}" neurons+ports (w2l_cpg.build_w2l_net; +4n/+20s/+4p per
          # kind, wired to the heel contact_hc extensor targets). RA port
          # current = ctn_amp*clip(rect(dF/dt)*tau_ra/cref, 0, 1); SA =
          # ctn_amp*clip(F/cref, 0, 1) (the level encoder into the slow
          # membrane). 0 = not built = census bit-identical.
          ra=0.0, sa=0.0,
          # RA/SA synaptic-gain overrides (nan = the W2L_GAINS defaults
          # ra_contact/sa_contact 0.08).
          ra_syn=float("nan"), sa_syn=float("nan"),
          # PRESYNAPTIC THRESHOLD/SATURATION (AnimatLab
          # PreSynapticThreshold/SaturationLevel, per Ben's tuning ask):
          # global defaults 0/5 = the historical hardcodes = BIT-IDENTICAL.
          # contact_presyn_lo/hi override the CONTACT edge families only
          # (nan = inherit the globals).
          presyn_lo=0.0, presyn_hi=5.0,
          contact_presyn_lo=float("nan"),
          contact_presyn_hi=float("nan"))
for arg in sys.argv[1:]:
    if arg.startswith("--") and "=" in arg:
        k, v = arg[2:].split("=", 1)
        try:
            v = float(v)
        except ValueError:
            pass
        KV[k.replace("-", "_")] = v

# rig base stiffnesses (world frame, applied at the Root body via
# xfrc_applied): [kx, ky N/m] horizontal leash on the COM, [kroll,
# kpitch, kyaw N m/rad] tilt assist. The VERTICAL channel is the drop
# protocol harness (feedforward weight fraction + stiff PD, computed
# in run_phase) — it replaced the old soft COM-z kz spring so the
# harness can actually hold the body in AIR (a 2000 N/m spring would
# need a 20 cm sag to carry the 411 N weight). Rationale for keeping
# a leash/tilt at all: the static margin probe (tmp\m6_static_margin.py)
# measured the COM 3.09 cm OUTSIDE the front edge of the support
# polygon at the grounded rest pose — no open-loop controller can
# free-stand this pose, matching the aproj shipping the Root with
# Freeze=True (it was never a free-standing model).
K0_XY, K0_ROT = 2000.0, 400.0
# rig damping is critical-rate: c = 2*sqrt(K*m_eff); m_eff 42 kg (trans),
# I_eff 5 kg m^2 (rot)
M_EFF, I_EFF = 42.0, 5.0

# the aproj contact plates = the heel/toe sensors
HEEL_GEOMS = ("foot_L_contact", "foot_R_contact")
TOE_GEOMS = ("toe_L_contact", "toe_R_contact")
JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
EXT_ACT = {"L": ("hip_L_ext", "knee_L_ext", "ankle_L_ext"),
           "R": ("hip_R_flx", "knee_R_ext", "ankle_R_ext")}  # R-hip crossing
FALL_Z = 0.60          # pelvis centre below this = fall (start ~0.993)
FALL_TILT = 45.0       # deg


def write_ground_xml() -> None:
    """Ground variant of w2l_mjcf_fixed.xml. The 2026-09-30 MJCF rebuild
    made the Root pelvis FREEJOINT-ed at the source (make_w2l_mjcf.py now
    emits <freejoint name="root"/> — the old M6 route INJECTED one here,
    which after the rebuild produced a duplicate joint name). So this is
    now a read-only copy + header comment; no structural edit."""
    txt = open(FIXED, encoding="utf-8").read()
    assert '<freejoint name="root"/>' in txt, \
        "w2l_mjcf_fixed.xml lost its root freejoint (2026-09-30 rebuild)"
    hdr = ("<!-- GROUND VARIANT for milestone 6 (test_w2l_ground.py):\n"
           "     identical to w2l_mjcf_fixed.xml (axis-fixed body, see\n"
           "     fix_joint_axes.py; Root pelvis freejoint INCLUDED since\n"
           "     the 2026-09-30 make_w2l_mjcf.py rebuild). Rest pose kept\n"
           "     verbatim (hovers 2.6-3.1 cm; the drop protocol harness\n"
           "     holds it in air, then lowers it). -->\n")
    open(GROUND, "w", encoding="utf-8").write(hdr + txt)


def rotvec_deg(R: np.ndarray, R0: np.ndarray) -> np.ndarray:
    """Rotation vector (world frame, radians) of R @ R0^T."""
    from mujoco import mju_mat2Quat, mju_quat2Vel
    Rrel = R @ R0.T
    q = np.zeros(4)
    mju_mat2Quat(q, Rrel.reshape(9))
    v = np.zeros(3)
    mju_quat2Vel(v, q, 1.0)
    return v


def run_phase(kind: str, dur: float, rig: float, sregime: str | None):
    """One closed-loop ground run. kind = 'stand'|'walk'. Returns logs."""
    if str(KV["net"]) == "bilateralrg2":
        # Ben's 2026-10-18 corrected aproj connectome (bilateralrg2
        # template; census 96/204). Port/muscle names differ from the M5
        # aff lineage - resolved through the _P/_IDX adapters below.
        sys.path.insert(0, os.path.join(HERE, "..", "w2l_cpg"))
        import build_w2l_net as B
        B.NAP["tau_max_h"] = KV["tau"]
        ra_on, sa_on = KV["ra"] > 0.0, KV["sa"] > 0.0
        gains = dict(contact_hc=KV["ch_hc"],
                     contact_ank=KV["ch_ank"],
                     contact_contra_flex=KV["ch_contra"],
                     c1_inh=3.0 * KV["comm"],
                     v3_weak=0.08 * KV["comm"])
        # presynaptic window (2026-10-20): only inject NON-default keys so
        # the default gains dict (and therefore the compiled net) is
        # EXACTLY the pre-knob build (bit-identical).
        if KV["presyn_lo"] != 0.0:
            gains["presyn_lo"] = KV["presyn_lo"]
        if KV["presyn_hi"] != 5.0:
            gains["presyn_hi"] = KV["presyn_hi"]
        if np.isfinite(KV["contact_presyn_lo"]) or \
                np.isfinite(KV["contact_presyn_hi"]):
            _clo = (KV["contact_presyn_lo"] if np.isfinite(
                KV["contact_presyn_lo"]) else KV["presyn_lo"])
            _chi = (KV["contact_presyn_hi"] if np.isfinite(
                KV["contact_presyn_hi"]) else KV["presyn_hi"])
            for _fam in ("contact", "contact_hc", "contact_ank",
                         "contact_contra_flex"):
                gains[f"{_fam}_presyn_lo"] = _clo
                gains[f"{_fam}_presyn_hi"] = _chi
        if np.isfinite(KV["ra_syn"]):
            gains["ra_contact"] = KV["ra_syn"]
        if np.isfinite(KV["sa_syn"]):
            gains["sa_contact"] = KV["sa_syn"]
        net = B.build(template_key="bilateralrg2", ib_ports=True,
                      alpha_l=KV["alpha_l"], alpha_r=KV["alpha_r"],
                      ra=ra_on, sa=sa_on, gains=gains)
        _knobs_on = (ra_on or sa_on or KV["presyn_lo"] != 0.0
                     or KV["presyn_hi"] != 5.0
                     or np.isfinite(KV["contact_presyn_lo"])
                     or np.isfinite(KV["contact_presyn_hi"]))
        if _knobs_on:
            print(f"[bilateralrg2] knobs ON: census {len(net.idx)}n/"
                  f"{net.n_synapses}s/{len(net.inputs)}p vs default runner "
                  f"build 96n/204s/22p (delta "
                  f"+{len(net.idx) - 96}n/+{net.n_synapses - 204}s/"
                  f"+{len(net.inputs) - 22}p)", flush=True)
        TAU_RA = B.TAU["ra"]
        from build_w2l_orig_net import MUSCLE_MAP
        MUSCLE_KEYS = {mlab: MUSCLE_MAP[mlab]
                       for mlab in net.muscle_outputs}   # template -> MJCF

        def _P(*names):
            for n in names:
                if n in net.inputs:
                    return net.input_index(n)
            raise KeyError(f"no port among {names}")

        def _IDX(*names):
            for n in names:
                if n in net.idx:
                    return net.idx[n]
            raise KeyError(f"no neuron among {names}")

        P_STIM2 = ("Stimulus_2", "Stimulus_2 (10nA)")
        P_HEEL = lambda s: (f"PORT heel {s}", f"{s} heel contact")  # noqa: E731
        P_TOE = lambda s: (f"PORT toe {s}", f"{s} toe contact")     # noqa: E731
        IDX_HEEL = lambda s: (f"heel {s}", f"{s} heel contact")     # noqa: E731
        IDX_TOE = lambda s: (f"toe {s}", f"{s} toe contact")        # noqa: E731

        def _muscle_ctrl(V):
            return {MUSCLE_KEYS[m]: c for m, c in net.muscle_ctrl(V).items()}

        ACTUATORS = list(MUSCLE_KEYS.values())
        IB_PORTS = [[p for p in net.inputs if p.startswith(f"PORT Ib {s}")]
                    for s in ("L", "R")]
        IIB = [[net.input_index(p) for p in pl] for pl in IB_PORTS]
        # RA/SA mechanoreceptor feed list: (port index, kind, plate geom
        # name) for every channel the build created (empty at defaults).
        RA_SA = [(net.input_index(f"{kind} {site} {s}"), kind, plate)
                 for kind, on in (("RA", ra_on), ("SA", sa_on)) if on
                 for site, plates in (("heel", HEEL_GEOMS),
                                      ("toe", TOE_GEOMS))
                 for s, plate in zip(("L", "R"), plates)]
    else:
        if KV["ra"] > 0.0 or KV["sa"] > 0.0:
            raise SystemExit("--ra/--sa mechanoreceptor channels are "
                             "implemented for --net=bilateralrg2 only "
                             "(w2l_cpg.build_w2l_net)")
        RA_SA, TAU_RA = [], None
        import build_w2l_aff_net as B
        B.SPLIT.NAP["tau_max_h"] = KV["tau"]
        # g_* passthrough (2026-10-16 wean): defaults == AFF_GAINS/comm 1.0,
        # so the unedited command line builds the IDENTICAL net (bit-exact).
        net = B.build(comm=KV["comm"], alpha_l=KV["alpha_l"],
                      alpha_r=KV["alpha_r"],
                      gains=dict(aff_heel=KV["g_heel"],
                                 aff_toe_in=KV["g_toe_in"],
                                 aff_toe_df=KV["g_toe_df"],
                                 aff_toe_df_mn=KV["g_toe_df_mn"],
                                 aff_ib=KV["g_ib"]))
        _P = lambda *n: net.input_index(n[0])                     # noqa: E731
        _IDX = lambda *n: net.idx[n[0]]                           # noqa: E731
        P_STIM2 = ("Stimulus_2",)
        P_HEEL = lambda s: (f"PORT heel {s}",)                    # noqa: E731
        P_TOE = lambda s: (f"PORT toe {s}",)                      # noqa: E731
        IDX_HEEL = lambda s: (f"heel {s}",)                       # noqa: E731
        IDX_TOE = lambda s: (f"toe {s}",)                         # noqa: E731
        _muscle_ctrl = net.muscle_ctrl
        IIB = [[net.input_index(f"PORT Ib {s}")] for s in ("L", "R")]
        ACTUATORS = list(net.muscle_outputs)
    # (mode-shared interface above: _P/_IDX/P_* lambdas, _muscle_ctrl,
    # IIB per-side port-index lists, ACTUATORS = MJCF actuator names)

    m = mujoco.MjModel.from_xml_path(GROUND)
    if KV["legcoll"] == 0:
        # L<->R leg collision exclusion (default-off; see KV comment).
        _L = ("femur_L", "tibia_L", "foot_L", "toe_L",
              "toe_L_contact", "foot_L_contact")
        _R = ("femur_R", "tibia_R", "foot_R", "toe_R",
              "toe_R_contact", "foot_R_contact")
        for names, ct in ((_L, 1), (_R, 2), (("Root", "ground"), 3)):
            for nm in names:
                gid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, nm)
                if gid >= 0:
                    m.geom_contype[gid] = ct
                    m.geom_conaffinity[gid] = ct
    d = mujoco.MjData(m)
    m.dof_damping[:] = KV["damp"]
    m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
    m.jnt_solref[:] = np.array([0.006, 1.0])
    # DROP PROTOCOL (Ben 2026-09-30): spawn at the VERBATIM rest pose
    # (plates hover 2.6-3.1 cm, M1 report section 2.8) — the old
    # spawn-time --drop grounding of the free-joint z is REMOVED. The
    # virtual-walker harness below holds the root in AIR, lowers it to
    # the platform over --lower s, and STAYS ON at --wsup weight
    # fraction (drop WITH the harness).
    jadr = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT,
                                           "root")]
    z_air = float(d.qpos[jadr + 2])             # rest-pose root z
    z_contact = z_air - KV["drop"] - KV["sink"]  # lowered target
    W = float(np.sum(m.body_mass)) * (-m.opt.gravity[2])   # N
    HOLD, LOW = KV["hold"], KV["lower"]
    # vertical harness = feedforward weight + stiff PD (the
    # test_li_stepping.py runner-rig idiom): K_V sized so the sag at
    # full weight is < drop/4
    K_V = 4.0 * W / max(KV["drop"], 0.01)
    C_V = 2.0 * float(np.sqrt(K_V * M_EFF))     # ~critical
    # RUNTIME STAND-IN (documented): the transported ankle range
    # [-20,-5] deg EXCLUDES the aproj's own rest pose (ankle 0 deg) —
    # the M3 report section 5 transport inconsistency. At t=0 the stiff
    # limit solver slams both ankles to -4.8 deg within 20 ms
    # (measured max|qacc| 5.6e3, tmp\m6_probe_fall.py) and the ground
    # run never gets past the transient. Widen the upper bound to +5 deg
    # so the rest pose is INSIDE range (the transported -20 deg PF bound
    # is kept); runtime-only, like the joint damping stand-in.
    for anm in ("ankle_L", "ankle_R"):
        jid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, anm)
        m.jnt_range[jid] = np.array([np.radians(-20.0), np.radians(5.0)])

    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in ACTUATORS}
    qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)]
            for j in JNTS}
    gid = {g: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g)
           for g in HEEL_GEOMS + TOE_GEOMS + ("ground",)}
    body_root = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, "Root")
    fmax = {a: max(float(m.actuator_gainprm[aid][2]), 1.0)
            for a, aid in act_ids.items()}

    u = net.make_inputs()
    iS1 = _P("Stimulus_1")
    iS2 = _P(*P_STIM2)
    itE = [net.input_index(f"TONIC {s} RG ext") for s in ("L", "R")]
    itF = [net.input_index(f"TONIC {s} RG flx") for s in ("L", "R")]
    # per-side GAMMA (RG tonic) + ALPHA (MN) resolution (defaults = the
    # symmetric te/tf baseline, no alpha ports -> bit-identical)
    te = (KV["te"] if not np.isfinite(KV["te_l"]) else KV["te_l"],
          KV["te"] if not np.isfinite(KV["te_r"]) else KV["te_r"])
    tf = (KV["tf"] if not np.isfinite(KV["tf_l"]) else KV["tf_l"],
          KV["tf"] if not np.isfinite(KV["tf_r"]) else KV["tf_r"])
    i_alpha = ([(net.input_index(p), KV["alpha_l"])
                for p in net.inputs if p.startswith("ALPHA L ")] +
               [(net.input_index(p), KV["alpha_r"])
                for p in net.inputs if p.startswith("ALPHA R ")])
    iheel = [_P(*P_HEEL(s)) for s in ("L", "R")]
    itoe = [_P(*P_TOE(s)) for s in ("L", "R")]
    iib = IIB                 # per side: list of port indices
    ivH = {s: _IDX(*IDX_HEEL(s)) for s in ("L", "R")}
    ivT = {s: _IDX(*IDX_TOE(s)) for s in ("L", "R")}
    iE = {s: net.idx[f"{s} RG ext"] for s in ("L", "R")}

    settle = KV["settle"]
    blend = 0.6
    # per-actuator drive cap: acap < 0 = uniform cap (the air-gate
    # convention). acap >= 0 caps the ANKLE actuators lower — the M3
    # report section 5 "ankle fidelity is poor" caveat: at the uniform
    # 0.5 cap the open-loop extension latch drives both plantarflexors
    # to heel-off within 100 ms (measured, tmp\m6_probe_closed.py:
    # ctrl 0.50, 1483 N, heel plate 40 -> 0 N) and the support polygon
    # collapses. A stance trim on the ankle PF is the minimal stance
    # posture fix; it is a documented stand-in knob, not silent tuning.
    acap = KV["cap"] if KV["acap"] < 0 else KV["acap"]
    nphy = int(dur / DT_PHY)
    log = dict(t=np.zeros(nphy), q=np.zeros((nphy, 6)),
               pel=np.zeros((nphy, 3)), com=np.zeros((nphy, 3)),
               tilt=np.zeros(nphy), heelF=np.zeros((nphy, 2)),
               toeF=np.zeros((nphy, 2)), heelS=np.zeros((nphy, 2)),
               toeS=np.zeros((nphy, 2)), heelV=np.zeros((nphy, 2)),
               toeV=np.zeros((nphy, 2)), rgE=np.zeros((nphy, 2)),
               heelin=np.zeros((nphy, 2)), toein=np.zeros((nphy, 2)),
               rigFz=np.zeros(nphy))
    c6 = np.zeros(6)
    cf = {g: 0.0 for g in HEEL_GEOMS + TOE_GEOMS}
    cfs = {g: 0.0 for g in HEEL_GEOMS + TOE_GEOMS}   # tangential (shear)
    f_prev = {g: 0.0 for g in HEEL_GEOMS + TOE_GEOMS}  # RA dF/dt memory
    DT_NET = NSUB * DT_PHY                            # 2 ms neural step
    CMODE = str(KV["contact_mode"])                   # normal|shear|total
    cap = KV["cap"]
    S_target = float(rig)
    anchor = None
    rv_prev = None
    fall_t = None
    for i in range(nphy):
        t = i * DT_PHY
        # current rig scale: settle catch (S=1) -> blend -> target
        if t < settle:
            S_now = 1.0
        elif t < settle + blend:
            S_now = 1.0 + (S_target - 1.0) * (t - settle) / blend
        else:
            S_now = S_target
        if (i % NSUB == 0):
            u[:] = 0.0
            if kind == "stand" and sregime == "latch":
                # extension-latch stance drive: tonic E on BOTH RGs (the
                # split net has no single shared RG), no F tonic, no
                # antiphase kickoff -> a steady extension posture, not
                # stepping. (The first smoke run drove NEITHER side here
                # — kickoff alone latches only the L RG in the split net;
                # found in smoke 2, fixed.)
                u[itE[0]] = te[0]
                u[itE[1]] = te[1]
            else:
                u[iS1] = 10.0 if t < 0.01 else 0.0       # antiphase kickoff
                u[iS2] = 10.0 if t < 0.01 else 0.0
                u[itE[0]] = te[0]
                u[itE[1]] = te[1]
                u[itF[0]] = tf[0]
                u[itF[1]] = tf[1]
            # i_a (2026-10-16 sweep fix): NEVER shadow the physics loop
            # index i -- the audit-stage `for i, a_cur in i_alpha` rebind
            # corrupted every even-index log entry AND silently disabled
            # the rig anchor capture (leash/tilt never engaged in alpha
            # runs; body drifted 11 m off-platform). alpha=0 = loop body
            # never runs = bit-identical defaults.
            for i_a, a_cur in i_alpha:    # per-side alpha (no-op default)
                u[i_a] = a_cur
            # Ib group current from ipsilateral extensor muscle force
            for k, s in enumerate(("L", "R")):
                fnorm = np.mean([max(d.actuator_force[act_ids[a]] /
                                     fmax[a], 0.0) for a in EXT_ACT[s]])
                for _p in iib[k]:
                    u[_p] = KV["ibnA"] * fnorm
            # REAL contact encoders (accumulated over the previous step)
            for k in range(2):
                u[iheel[k]] = KV["ctn_amp"] * min(cf[HEEL_GEOMS[k]] /
                                                  KV["cref"], 1.0)
                u[itoe[k]] = KV["ctn_amp"] * min(cf[TOE_GEOMS[k]] /
                                                 KV["cref"], 1.0)
            # RA/SA mechanoreceptor channel feeds (2026-10-20): driven by
            # the SAME plate force as the encoders above, so --contact-mode
            # (normal|shear|total) applies here too. RA = rectified dF/dt
            # across one neural step, scaled so a ramp reaching cref
            # within tau_ra saturates the encoder (stand-in; the 20 ms
            # RA membrane then shapes the transient into a short pulse).
            # SA = the sustained LEVEL (same encoder formula as the base
            # channel; the 200 ms membrane + its own synaptic gain do the
            # differentiating).
            if RA_SA:
                for iP, kind, plate in RA_SA:
                    if kind == "RA":
                        dfdt = max((cf[plate] - f_prev[plate]) / DT_NET,
                                   0.0)
                        u[iP] = KV["ctn_amp"] * min(dfdt * TAU_RA /
                                                    KV["cref"], 1.0)
                    else:
                        u[iP] = KV["ctn_amp"] * min(cf[plate] /
                                                    KV["cref"], 1.0)
                for g in f_prev:
                    f_prev[g] = cf[g]
            V = net.step(u)
            for a, c in _muscle_ctrl(V).items():
                d.ctrl[act_ids[a]] = min(c, acap if "ankle" in a else cap)
            d.xfrc_applied[body_root, :] = 0.0     # rig force is per-step
            # VIRTUAL-WALKER VERTICAL HARNESS (the drop-protocol piece;
            # always on, NOT scaled by rig S): feedforward weight +
            # vertical PD at the protocol target height. Air-hold ->
            # linear descent -> retained support (Ben 2026-09-30).
            if t < HOLD:                      # AIR: feet clear, air-step
                z_t, wff = z_air, 1.0
            elif t < HOLD + LOW:              # LOWER to the platform
                z_t = z_air + (z_contact - z_air) * (t - HOLD) / \
                    max(LOW, 1e-6)
                wff = 1.0
            else:                             # ON platform, harness kept
                z_t, wff = z_contact, KV["wsup"]
            Fz = max(wff * W + K_V * (z_t - float(d.qpos[jadr + 2]))
                     - C_V * d.qvel[2], 0.0)
            d.xfrc_applied[body_root, 2] += Fz
            if anchor is not None and S_now > 0.0:
                # leash (COM xy) + tilt assist, scaled by S (the rig)
                if KV["walker_yield"] > 0.0 and i % NSUB == 0 and i > 0:
                    # YIELDING VIRTUAL WALKER (default-off, 2026-10-16
                    # sweep): Ben's AnimatLab protocol uses the MOVING
                    # WalkingPath body as the walker; the fixed anchor
                    # here pins the COM so sustained progression is
                    # structurally impossible (measured: winner combo
                    # shifts +0.10 m during the drop, then flat 22 s).
                    # Anchor x follows a slow EMA of the COM (tau =
                    # walker_yield s): the leash still damps FAST
                    # deviations (the trolley has inertia) but yields
                    # to slow step thrust. Default 0 = fixed anchor =
                    # bit-identical.
                    a = float(np.exp(-(NSUB * DT_PHY) / KV["walker_yield"]))
                    anchor["com"][0] = a * anchor["com"][0] + \
                        (1.0 - a) * float(d.subtree_com[0][0])
                com = anchor["com"]
                dv = d.subtree_com[0] - com
                F = np.zeros(3)
                F[0] = -S_now * K0_XY * dv[0]
                F[1] = -S_now * K0_XY * dv[1]
                rv = rotvec_deg(d.xmat[body_root].reshape(3, 3),
                                anchor["R0"])
                T = -S_now * K0_ROT * rv
                if rv_prev is not None and i > 0:
                    T -= S_now * 2.0 * np.sqrt(K0_ROT * I_EFF) * \
                        (rv - rv_prev) / DT_PHY
                rv_prev = rv
                F[:2] -= S_now * 2.0 * np.sqrt(K0_XY * M_EFF) * \
                    np.array([d.qvel[0], d.qvel[1]])
                d.xfrc_applied[body_root, 0] += F[0]
                d.xfrc_applied[body_root, 1] += F[1]
                d.xfrc_applied[body_root, 3:] = T
            log["rigFz"][i] = d.xfrc_applied[body_root, 2]
        mujoco.mj_step(m, d)
        if i % NSUB == 0 and t >= settle and anchor is None:
            # re-anchor the rig at the SETTLED pose (feet loaded)
            anchor = dict(com=d.subtree_com[0].copy(),
                          pel=d.xpos[body_root].copy(),
                          R0=d.xmat[body_root].reshape(3, 3).copy())
        # contact forces (this step) for the NEXT encoder update.
        # --contact-mode (2026-10-20, Ben's shear ask): 'normal' (default)
        # = the historical z-only accumulation, BIT-IDENTICAL; 'shear' =
        # tangential magnitude only (the plate as a pure slip sensor);
        # 'total' = hypot(normal, shear). The tangential magnitude is
        # ALWAYS accumulated into cfs/heelS/toeS for the logs.
        for g in cf:
            cf[g] = 0.0
            cfs[g] = 0.0
        for k in range(d.ncon):
            con = d.contact[k]
            for nm, g2 in ((HEEL_GEOMS[0], gid["foot_L_contact"]),
                           (HEEL_GEOMS[1], gid["foot_R_contact"]),
                           (TOE_GEOMS[0], gid["toe_L_contact"]),
                           (TOE_GEOMS[1], gid["toe_R_contact"])):
                if con.geom1 == g2 or con.geom2 == g2:
                    mujoco.mj_contactForce(m, d, k, c6)
                    fn = abs(float(c6[0]))
                    ft = float(np.hypot(c6[1], c6[2]))
                    cfs[nm] += ft
                    if CMODE == "shear":
                        cf[nm] += ft
                    elif CMODE == "total":
                        cf[nm] += float(np.hypot(fn, ft))
                    else:
                        cf[nm] += fn
        log["t"][i] = t
        log["q"][i] = [np.degrees(d.qpos[qadr[j]]) for j in JNTS]
        log["pel"][i] = d.xpos[body_root]
        log["com"][i] = d.subtree_com[0]
        log["heelF"][i] = [cf[HEEL_GEOMS[0]], cf[HEEL_GEOMS[1]]]
        log["toeF"][i] = [cf[TOE_GEOMS[0]], cf[TOE_GEOMS[1]]]
        log["heelS"][i] = [cfs[HEEL_GEOMS[0]], cfs[HEEL_GEOMS[1]]]
        log["toeS"][i] = [cfs[TOE_GEOMS[0]], cfs[TOE_GEOMS[1]]]
        if anchor is None:
            log["tilt"][i] = 0.0
        else:
            log["tilt"][i] = np.degrees(
                np.linalg.norm(rotvec_deg(d.xmat[body_root].reshape(3, 3),
                                          anchor["R0"])))
        if i % NSUB == 0:
            log["rgE"][i] = [V[iE["L"]], V[iE["R"]]]
            log["heelV"][i] = [V[ivH[s]] for s in ("L", "R")]
            log["toeV"][i] = [V[ivT[s]] for s in ("L", "R")]
            log["heelin"][i] = [u[iheel[0]], u[iheel[1]]]
            log["toein"][i] = [u[itoe[0]], u[itoe[1]]]
        else:
            log["rgE"][i] = log["rgE"][i - 1]
            log["heelV"][i] = log["heelV"][i - 1]
            log["toeV"][i] = log["toeV"][i - 1]
            log["heelin"][i] = log["heelin"][i - 1]
            log["toein"][i] = log["toein"][i - 1]
        if fall_t is None and (d.xpos[body_root][2] < FALL_Z or
                               log["tilt"][i] > FALL_TILT):
            fall_t = t
        if not np.isfinite(d.qpos).all():
            print(f"[FAIL] qpos non-finite at t={t:.3f} s")
            log["finite"] = False
            log["fall_t"] = fall_t if fall_t is not None else t
            return log
    log["finite"] = True
    log["fall_t"] = fall_t
    return log


def stand_report(rig: float, log: dict, regime: str) -> dict:
    # window = after the drop (hold + lower + settle transient); the
    # settle/blend window sits INSIDE the air-hold under the protocol
    t0 = max(KV["settle"] + 0.6, KV["hold"] + KV["lower"] + 0.5)
    w = log["t"] >= t0
    t = log["t"][w]
    com = log["com"][w]
    pel = log["pel"][w]
    tilt = log["tilt"][w]
    heelF = log["heelF"][w]
    toeF = log["toeF"][w]
    span = t[-1] - t[0]
    res = dict(rig=rig, regime=regime, span=span, finite=log["finite"],
               fall_t=log["fall_t"],
               sway_x=(float(com[:, 0].std()), float(com[:, 0].ptp())),
               sway_y=(float(com[:, 1].std()), float(com[:, 1].ptp())),
               com_z=(float(com[:, 2].mean()), float(com[:, 2].min())),
               pel_z=(float(pel[:, 2].mean()), float(pel[:, 2].min())),
               tilt=(float(tilt.mean()), float(tilt.max())),
               duty_L=float((heelF[:, 0] > KV["heelN"]).mean()),
               duty_R=float((heelF[:, 1] > KV["heelN"]).mean()),
               harness_pct=float(100.0 * log["rigFz"][w].mean() / 411.0))
    return res


def main() -> int:
    write_ground_xml()
    phase = str(KV["phase"])
    ok = True
    stand_full = False
    stand_rows = []
    walk = None

    if phase in ("stand", "both"):
        print(f"== M6 gate (a): GROUND STANDING BALANCE, afferented walker, "
              f"pelvis FREE, {KV['dur_stand']:.0f} s ==")
        print(f"   knobs: sregime={KV['sregime']} ctn_amp={KV['ctn_amp']} "
              f"cref={KV['cref']} ibnA={KV['ibnA']} damp={KV['damp']} "
              f"settle={KV['settle']} (+0.6 s blend) | DROP PROTOCOL: "
              f"hold={KV['hold']} s -> lower={KV['lower']} s -> retained "
              f"wsup={KV['wsup']} (drop={KV['drop']}+sink={KV['sink']} m) "
              f"| rig kxy={K0_XY:.0f} N/m krot={K0_ROT:.0f} N m/rad | "
              f"fall: pelvis z<{FALL_Z} m or tilt>{FALL_TILT} deg")
        for rig in (0.0, 0.25, 0.5, 1.0):
            log = run_phase("stand", KV["dur_stand"], rig, KV["sregime"])
            res = stand_report(rig, log, KV["sregime"])
            stand_rows.append(res)
            fell = res["fall_t"] is not None
            print(f"   rig S={rig:.2f}: finite={res['finite']} "
                  f"fall={'t=%.2f s' % res['fall_t'] if fell else 'no'} | "
                  f"sway_x std/range {res['sway_x'][0]*100:.1f}/"
                  f"{res['sway_x'][1]*100:.1f} cm | "
                  f"sway_y std/range {res['sway_y'][0]*100:.1f}/"
                  f"{res['sway_y'][1]*100:.1f} cm | "
                  f"COM z mean/min {res['com_z'][0]:.3f}/"
                  f"{res['com_z'][1]:.3f} m | "
                  f"tilt mean/max {res['tilt'][0]:.1f}/{res['tilt'][1]:.1f} deg | "
                  f"heel duty L/R {res['duty_L']:.2f}/{res['duty_R']:.2f} | "
                  f"harness carries {res['harness_pct']:+.0f}% weight")
            if rig == 0.0:
                stand_full = (not fell and res["finite"]
                              and res["tilt"][1] < 15.0
                              and res["sway_x"][1] < 0.10
                              and res["sway_y"][1] < 0.10)
        print(f"   balance metrics table: {len(stand_rows)} rig scales")
        print(f"   [{'PASS-FULL' if stand_full else 'PARTIAL'}] unrigged_stand")

    if phase in ("walk", "both"):
        print(f"\n== M6 gate (b): CONTACT-DRIVEN GROUND WALKING attempt, "
              f"real heel/toe contact forces, {KV['dur_walk']:.0f} s, "
              f"rig S={KV['rig']} ==")
        print(f"   knobs: ctn_amp={KV['ctn_amp']} cref={KV['cref']} N "
              f"ibnA={KV['ibnA']} te={KV['te']} tf={KV['tf']} tau={KV['tau']} "
              f"cap={KV['cap']} acap={KV['acap']} damp={KV['damp']} | "
              f"DROP PROTOCOL: hold={KV['hold']} s lower={KV['lower']} s "
              f"wsup={KV['wsup']} drop={KV['drop']}+sink={KV['sink']} m "
              f"rig S={KV['rig']}")
        wlog = run_phase("walk", KV["dur_walk"], KV["rig"], None)
        t0 = KV["hold"] + KV["lower"] + 0.5   # after the drop settles
        w = wlog["t"] >= t0
        t = wlog["t"][w]
        q = -wlog["q"][w]                    # flexion-positive
        pel = wlog["pel"][w]
        com = wlog["com"][w]
        tilt = wlog["tilt"][w]
        heelF = wlog["heelF"][w]
        toeF = wlog["toeF"][w]
        rg = wlog["rgE"][w]
        span = t[-1] - t[0]
        # heel-strike cycles (rising through 20 N, 0.3 s refractory)
        strikes = {}
        for k, s in enumerate(("L", "R")):
            on = heelF[:, k] > 20.0
            st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
            keep = []
            for x in st:
                if not keep or (x - keep[-1]) * DT_PHY > 0.3:
                    keep.append(x)
            strikes[s] = np.array(keep)
        per = {s: (float(np.median(np.diff(strikes[s])) * DT_PHY)
                   if len(strikes[s]) >= 3 else float("nan"))
               for s in strikes}
        duty = {s: float((heelF[:, k] > KV["heelN"]).mean())
                for k, s in enumerate(("L", "R"))}
        # RG-E bursts + antiphase
        stL = bursts_refractory(rg[:, 0], dt=DT_PHY, min_gap_s=0.4)
        stR = bursts_refractory(rg[:, 1], dt=DT_PHY, min_gap_s=0.4)
        rg_per = (float(np.diff(stL).mean() * DT_PHY)
                  if len(stL) >= 3 else float("nan"))
        rg_r = float(np.corrcoef(rg[::5, 0], rg[::5, 1])[0, 1])
        fwd = float(com[-1, 0] - com[0, 0])
        v = fwd / span
        fell = wlog["fall_t"] is not None
        print(f"   finite={wlog['finite']} | "
              f"fall={'t=%.2f s' % wlog['fall_t'] if fell else 'no'} | "
              f"pelvis z min {pel[:, 2].min():.3f} m "
              f"(start 0.993) | tilt max {tilt.max():.1f} deg | "
              f"harness carries "
              f"{100.0 * wlog['rigFz'][w].mean() / 411.0:+.0f}% weight")
        print(f"   cycles (heel strikes >20 N, window {span:.1f} s): "
              f"L {len(strikes['L'])} (interval {per['L']:.3f} s) | "
              f"R {len(strikes['R'])} (interval {per['R']:.3f} s)")
        print(f"   heel duty L {duty['L']:.2f} R {duty['R']:.2f} "
              f"| forward disp {fwd:+.3f} m ({v:+.3f} m/s)")
        print(f"   flexion ranges (deg): hip L {q[:,0].ptp():.1f} "
              f"R {q[:,3].ptp():.1f} | knee L {q[:,1].ptp():.1f} "
              f"R {q[:,4].ptp():.1f} | ankle L {q[:,2].ptp():.1f} "
              f"R {q[:,5].ptp():.1f}")
        print(f"   flexion min..max: hip L [{q[:,0].min():+.1f},"
              f"{q[:,0].max():+.1f}] knee L [{q[:,1].min():+.1f},"
              f"{q[:,1].max():+.1f}] ankle L [{q[:,2].min():+.1f},"
              f"{q[:,2].max():+.1f}]")
        print(f"   neural: L RG-E bursts {len(stL)} "
              f"(period {rg_per if np.isfinite(rg_per) else float('nan'):.3f} s), "
              f"R RG-E bursts {len(stR)}; L/R RG-E r {rg_r:+.3f}")
        hi_, ti_ = wlog["heelin"][w], wlog["toein"][w]
        amp_ = float(KV["ctn_amp"])
        print(f"   contact encoder evidence: heel SN max "
              f"{wlog['heelV'][w].max():.2f} mV, toe SN max "
              f"{wlog['toeV'][w].max():.2f} mV, heel port current max "
              f"{hi_.max():.2f} nA")
        # encoder duty split (2026-10-20 diagnosis ask): is the channel
        # silent, partially active, or pinned at the ctn_amp ceiling?
        for nm_, x_ in (("heel", hi_), ("toe", ti_)):
            print(f"   {nm_} port current nA (L/R): mean "
                  f"{x_[:, 0].mean():.3f}/{x_[:, 1].mean():.3f}, max "
                  f"{x_[:, 0].max():.3f}/{x_[:, 1].max():.3f} | frac "
                  f"steps ACTIVE (I>0) "
                  f"{(x_[:, 0] > 1e-9).mean():.3f}/"
                  f"{(x_[:, 1] > 1e-9).mean():.3f}, SATURATED "
                  f"(I>=0.995*ctn_amp) "
                  f"{(x_[:, 0] >= 0.995 * amp_).mean():.3f}/"
                  f"{(x_[:, 1] >= 0.995 * amp_).mean():.3f} "
                  f"(ceiling {amp_:.2f} nA, mode {KV['contact_mode']})")
        print(f"   plate shear (tangential) N: heel mean/max "
              f"{wlog['heelS'][w].mean():.1f}/{wlog['heelS'][w].max():.1f}, "
              f"toe mean/max {wlog['toeS'][w].mean():.1f}/"
              f"{wlog['toeS'][w].max():.1f} (sensed mode: "
              f"{KV['contact_mode']})")
        walked = (not fell and wlog["finite"]
                  and len(strikes["L"]) >= 3 and len(strikes["R"]) >= 3
                  and 0.3 <= min(duty.values()) and max(duty.values()) <= 0.9
                  and abs(fwd) > 0.2)
        shuffled = (not fell and wlog["finite"] and not walked
                    and len(strikes["L"]) + len(strikes["R"]) >= 4)
        verdict = ("WALKS" if walked else
                   "SHUFFLES-IN-PLACE" if shuffled else "FALLS/NO-GAIT")
        print(f"   walk verdict: {verdict}")
        walk = dict(finite=wlog["finite"], fall_t=wlog["fall_t"],
                    strikes={s: int(len(strikes[s])) for s in ("L", "R")},
                    per=per, duty=duty, fwd=fwd, v=v,
                    ranges={n: (float(q[:, c0].ptp()), float(q[:, c1].ptp()))
                            for n, (c0, c1) in (("hip", (0, 3)),
                                                ("knee", (1, 4)),
                                                ("ankle", (2, 5)))},
                    tilt_max=float(tilt.max()), pel_z_min=float(pel[:, 2].min()),
                    rg_bursts=(int(len(stL)), int(len(stR))), rg_per=rg_per,
                    rg_r=rg_r, verdict=verdict)
        ok = walked

    # ------------------------------------------------------------- verdict
    if phase == "stand":
        verdict = "PASS (FULL: unrigged stand)" if stand_full else \
            "PARTIAL (rigged stand only)" if any(
                r["fall_t"] is None for r in stand_rows) else "FAIL"
    elif phase == "walk":
        verdict = f"WALK={walk['verdict']}"
    else:
        verdict = (f"stand={'FULL' if stand_full else 'PARTIAL/FAIL'}; "
                   f"walk={walk['verdict'] if walk else 'n/a'}")
    print(f"\nM6 VERDICT: {verdict}")
    return 0 if (phase != "walk" or walked) else 1


if __name__ == "__main__":
    sys.exit(main())
