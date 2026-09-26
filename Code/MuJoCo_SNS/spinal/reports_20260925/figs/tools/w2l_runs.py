r"""FIG RUNNER — replays the 2026-09-25 walker-campaign gate runs and captures
offscreen render frames for the figure set. STANDALONE: imports the milestone
machinery read-only (build_w2l_*_net, test_w2l_air/afferented/ground); the
milestone scripts themselves are NOT modified.

Usage (cwd = this folder):
    C:\Users\Ben Bolen\.conda\envs\myo\python.exe w2l_runs.py m3|m4|m5|m6

Outputs (under figs\tools\):
    data_<run>.npz      logs (t, joint angles deg, RG voltages, ...)
    frames_<run>\*.png  raw offscreen renders (air runs only)

Configs are the canonical gate commands:
    m3: test_w2l_air.py                       (net=orig, te3 tf4 tau.25 cap.5 damp3, 20 s)
    m4: test_w2l_air.py --net=split --comm=1.0
    m5: test_w2l_air_afferented.py --dur=20 --amp=1.5 --pamp=4.0 --pulse=0.5
        --ibnA=1  (control; then the gate-b causal heel-L pulse run, 4 nA x 0.5 s)
    m6: test_w2l_ground.py --phase=walk --dur_walk=20 --rig=1.0 --acap=0.15 --ibnA=0
"""
from __future__ import annotations

import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
# NOTE: do NOT re-wrap sys.stdout here — test_w2l_air/ground wrap it again on
# import and the stacked TextIOWrapper gets GC-closed mid-run (ValueError:
# I/O operation on closed file). reconfigure() mutates in place.
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")

import numpy as np

W2L = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
FIGS = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\figs"
TOOLS = os.path.join(FIGS, "tools")
sys.path.insert(0, W2L)

import mujoco  # noqa: E402
import test_w2l_air as TW  # noqa: E402  (DT_PHY, NSUB, write_air_xml, helpers)

# test_w2l_air/ground REASSIGN sys.stdout with a new TextIOWrapper over the
# SAME underlying buffer; when an older wrapper is released it is gc-closed
# and takes the shared buffer with it ("I/O operation on closed file").
# Keep strong references so no wrapper that wraps the real buffer is ever
# collected mid-run.
_STDOUT_KEEP = [sys.stdout]

from test_w2l_air import DT_PHY, NSUB, acf_period, band_limited  # noqa: E402
from test_w2l_air import bursts_refractory, swing_peaks, write_air_xml  # noqa

JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
DUR = 20.0
FRAME_EVERY = 100          # physics steps -> 0.1 s per frame, 200 frames
RW, RH = 360, 470           # render size


def frame_times(nphy: int) -> np.ndarray:
    return np.arange(0, nphy, FRAME_EVERY) * DT_PHY


class FrameGrabber:
    """mujoco.Renderer wrapper (mujoco 2.3.7 has no Renderer.close)."""

    def __init__(self, model, tag: str):
        self.r = mujoco.Renderer(model, height=RH, width=RW)
        self.cam = mujoco.MjvCamera()
        self.cam.lookat[:] = [-3.40, 0.0, 1.05]
        self.cam.azimuth = 90.0
        self.cam.elevation = -3.0
        self.cam.distance = 2.35
        self.outdir = os.path.join(TOOLS, "frames_" + tag)
        os.makedirs(self.outdir, exist_ok=True)
        self.n = 0

    def grab(self, d, i: int) -> None:
        self.r.update_scene(d, camera=self.cam)
        rgb = self.r.render()
        from PIL import Image
        Image.fromarray(rgb).save(
            os.path.join(self.outdir, "frame_%04d.png" % self.n))
        self.n += 1


def common_body_setup(m, d, cap: float, damp: float = 3.0, stiff: float = 1.0):
    m.dof_damping[:] = damp
    if stiff > 0.0:
        m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
        m.jnt_solref[:] = np.array([0.006, 1.0])
    act_ids = {}
    return act_ids


def run_air(kind: str) -> None:
    """M3 (orig) / M4 (split coupled) closed loop, verbatim drive of
    test_w2l_air.main + frame capture."""
    TE, TF, TAU, CAP, LIFT = 3.0, 4.0, 0.25, 0.5, 0.30
    if kind == "m3":
        import build_w2l_orig_net as B
        B.NAP["tau_max_h"] = TAU
        net = B.build()
    else:
        import build_w2l_split_net as B
        B.NAP["tau_max_h"] = TAU
        net = B.build(comm=1.0)
    SPLIT = kind == "m4"

    TW.KV["lift"] = LIFT
    write_air_xml()
    m = mujoco.MjModel.from_xml_path(TW.AIR)
    d = mujoco.MjData(m)
    common_body_setup(m, d, CAP)

    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)]
            for j in JNTS}
    g_ground = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "ground")

    u = net.make_inputs()
    iE = net.input_index("TONIC L RG ext")
    iF = net.input_index("TONIC L RG flx")
    iS = net.input_index("Stimulus_1")
    if SPLIT:
        iEr = net.input_index("TONIC R RG ext")
        iFr = net.input_index("TONIC R RG flx")
        iS2 = net.input_index("Stimulus_2")

    nphy = int(DUR / DT_PHY)
    tlog = np.zeros(nphy)
    qlog = np.zeros((nphy, 6))
    rgL = np.zeros(nphy)
    rgR = np.zeros(nphy) if SPLIT else np.zeros(1)
    g_contacts = other_contacts = 0
    grab = FrameGrabber(m, kind)
    fi = 0
    ftimes = frame_times(nphy)
    V = np.zeros(1)
    for i in range(nphy):
        t = i * DT_PHY
        if i % NSUB == 0:
            u[iE] = TE
            u[iF] = TF
            u[iS] = 10.0 if t < 0.01 else 0.0
            if SPLIT:
                u[iEr] = TE
                u[iFr] = TF
                u[iS2] = 10.0 if t < 0.01 else 0.0
            V = net.step(u)
            for a, c in net.muscle_ctrl(V).items():
                d.ctrl[act_ids[a]] = min(c, CAP)
        mujoco.mj_step(m, d)
        for k in range(d.ncon):
            con = d.contact[k]
            if con.geom1 == g_ground or con.geom2 == g_ground:
                g_contacts += 1
            else:
                other_contacts += 1
        qlog[i] = [np.degrees(d.qpos[qadr[j]]) for j in JNTS]
        if i % NSUB == 0:
            rgL[i] = V[net.idx["L RG ext"]]
            if SPLIT:
                rgR[i] = V[net.idx["R RG ext"]]
            else:
                rgR[0] = 0.0
        else:
            rgL[i] = rgL[i - 1]
            if SPLIT:
                rgR[i] = rgR[i - 1]
        tlog[i] = t
        if fi < len(ftimes) and i == int(round(ftimes[fi] / DT_PHY)):
            grab.grab(d, i)
            fi += 1
        if not np.isfinite(d.qpos).all():
            print(f"[FAIL] qpos non-finite at t={t:.3f} s")
            return

    # sanity metrics vs the gate reports
    an = slice(int(2.0 / DT_PHY), None)
    flex = -qlog[an]
    st = bursts_refractory(rgL[an], dt=DT_PHY, min_gap_s=0.4)
    per = float(np.diff(st).mean() * DT_PHY) if len(st) >= 3 else float("nan")
    fb, rb = band_limited(flex[:, 0]), band_limited(flex[:, 3])
    r_band = float(np.corrcoef(fb, rb)[0, 1])
    hipL = float(flex[:, 0].ptp())
    print(f"[{kind}] ground contacts {g_contacts} | self {other_contacts} | "
          f"RG-E period {per:.3f} s | hip L range {hipL:.1f} deg | "
          f"band r {r_band:+.3f}")
    if SPLIT:
        rgr_r = float(np.corrcoef(rgL[an][::5], rgR[an][::5])[0, 1])
        print(f"[{kind}] R RG-E r {rgr_r:+.3f} (report: -0.750, period 1.357 s)")
    else:
        print(f"[{kind}] (M3 report: period 1.027 s, hip L 41.8 deg, r -0.623)")

    out = dict(t=tlog, q=qlog, rgL=rgL, ftimes=ftimes,
               g_contacts=g_contacts, other_contacts=other_contacts)
    if SPLIT:
        out["rgR"] = rgR
    np.savez_compressed(os.path.join(TOOLS, f"data_{kind}.npz"), **out)
    print(f"[{kind}] saved data + {grab.n} frames")


def run_m5() -> None:
    """M5: control walk (canonical flags) then the gate-b causal heel-L
    pulse run; frames captured on the PERTURBED run (the pulse moment)."""
    import build_w2l_aff_net as BA
    BA.SPLIT.NAP["tau_max_h"] = 0.25

    import test_w2l_air_afferented as TA
    TA.KV.update(dur=DUR, causal="none", amp=1.5, pulse=0.5, ibnA=1.0)
    A = TA.run_sim()
    an = slice(int(2.0 / DT_PHY), None)
    onA = [t for t in A["onsets"]["L"] if t >= 2.0]
    T_hat = float(np.median(np.diff(onA))) if len(onA) >= 3 else 1.03
    tstar = onA[2] + TA.KV["tfac"] * T_hat
    flex = -A["q"][an]
    st = bursts_refractory(A["rgL"][an], dt=DT_PHY, min_gap_s=0.4)
    per = float(np.diff(st).mean() * DT_PHY) if len(st) >= 3 else float("nan")
    hipL = float(flex[:, 0].ptp())
    heel_pk = float(A["heel"][an].max())
    print(f"[m5 control] period {per:.3f} s | hip L {hipL:.1f} deg | "
          f"heel SN max {heel_pk:.2f} mV | onsets {len(onA)} | "
          f"tstar {tstar:.3f} s   (M5 report: 1.389 s / 33.1 / 1.50 mV / "
          f"t* 6.244 s)")
    np.savez_compressed(
        os.path.join(TOOLS, "data_m5c.npz"),
        t=A["t"], q=A["q"], rgL=A["rgL"], rgR=A["rgR"],
        heelSN=A["heel"], heelIN=A["ibin"] * 0.0,  # placeholder, see below
        onL=np.array(A["onsets"]["L"]), onR=np.array(A["onsets"]["R"]),
        tstar=tstar, ftimes=np.zeros(1))

    # ---- perturbed run with the extra causal heel-L pulse + frames
    import mujoco as mj
    net = BA.build(comm=1.0)
    TA.KV["lift"] = 0.30
    write_air_xml()
    m = mj.MjModel.from_xml_path(TW.AIR)
    d = mj.MjData(m)
    m.dof_damping[:] = 3.0
    m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
    m.jnt_solref[:] = np.array([0.006, 1.0])
    act_ids = {a: mj.mj_name2id(m, mj.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    qadr = {j: m.jnt_qposadr[mj.mj_name2id(m, mj.mjtObj.mjOBJ_JOINT, j)]
            for j in JNTS}
    g_ground = mj.mj_name2id(m, mj.mjtObj.mjOBJ_GEOM, "ground")
    u = net.make_inputs()
    iS1 = net.input_index("Stimulus_1")
    iS2 = net.input_index("Stimulus_2")
    itE = [net.input_index(f"TONIC {s} RG ext") for s in ("L", "R")]
    itF = [net.input_index(f"TONIC {s} RG flx") for s in ("L", "R")]
    iheel = [net.input_index(f"PORT heel {s}") for s in ("L", "R")]
    itoe = [net.input_index(f"PORT toe {s}") for s in ("L", "R")]
    iib = [net.input_index(f"PORT Ib {s}") for s in ("L", "R")]
    iE_L, iE_R = net.idx["L RG ext"], net.idx["R RG ext"]
    iH = [net.idx[f"heel {s}"] for s in ("L", "R")]
    fmax = {a: max(float(m.actuator_gainprm[aid][2]), 1.0)
            for a, aid in act_ids.items()}
    TE, TF, CAP = 3.0, 4.0, 0.5
    PAMP, PDUR = 4.0, 0.5

    sch = {s: dict(last=-9.0, in_burst=False, T=1.03,
                   heel_to=-1.0, toe_from=-1.0, toe_to=-1.0,
                   rmax=1e-9, onsets=[]) for s in ("L", "R")}
    nphy = int(DUR / DT_PHY)
    tlog = np.zeros(nphy)
    qlog = np.zeros((nphy, 6))
    rgL = np.zeros(nphy)
    rgR = np.zeros(nphy)
    heelSN = np.zeros((nphy, 2))
    heelIN = np.zeros((nphy, 2))
    g_contacts = 0
    grab = FrameGrabber(m, "m5")
    fi = 0
    ftimes = frame_times(nphy)
    V = np.zeros(1)
    for i in range(nphy):
        t = i * DT_PHY
        if i % NSUB == 0:
            u[:] = 0.0
            u[iS1] = 10.0 if t < 0.01 else 0.0
            u[iS2] = 10.0 if t < 0.01 else 0.0
            u[itE[0]] = u[itE[1]] = TE
            u[itF[0]] = u[itF[1]] = TF
            for k, s in enumerate(("L", "R")):
                fnorm = np.mean([max(d.actuator_force[act_ids[a]] /
                                     fmax[a], 0.0) for a in TA.EXT_ACT[s]])
                u[iib[k]] = 1.0 * fnorm
            for k, s in enumerate(("L", "R")):
                if t >= 1.5:
                    if sch[s]["heel_to"] > t:
                        u[iheel[k]] = 1.5
                    if sch[s]["toe_from"] <= t <= sch[s]["toe_to"]:
                        u[itoe[k]] = 1.5
            if tstar <= t < tstar + PDUR:
                u[iheel[0]] = PAMP            # the EXTRA causal heel-L pulse
            V = net.step(u)
            for k, (s, ie) in enumerate((("L", iE_L), ("R", iE_R))):
                stt = sch[s]
                stt["rmax"] = max(stt["rmax"], V[ie])
                active = V[ie] > 0.5 * stt["rmax"]
                if active and not stt["in_burst"] and t - stt["last"] > 0.4:
                    if stt["last"] > 0:
                        stt["T"] = 0.5 * stt["T"] + 0.5 * (t - stt["last"])
                    stt["last"] = t
                    stt["onsets"].append(t)
                    stt["heel_to"] = t + 0.02 + 0.38 * stt["T"]
                    stt["toe_from"] = t + 0.62 * stt["T"]
                    stt["toe_to"] = t + 0.92 * stt["T"]
                stt["in_burst"] = active
            for a, c in net.muscle_ctrl(V).items():
                d.ctrl[act_ids[a]] = min(c, CAP)
        mj.mj_step(m, d)
        for k in range(d.ncon):
            con = d.contact[k]
            if con.geom1 == g_ground or con.geom2 == g_ground:
                g_contacts += 1
        qlog[i] = [np.degrees(d.qpos[qadr[j]]) for j in JNTS]
        if i % NSUB == 0:
            rgL[i] = V[iE_L]
            rgR[i] = V[iE_R]
            heelSN[i] = [V[iH[0]], V[iH[1]]]
            heelIN[i] = [u[iheel[0]], u[iheel[1]]]
        else:
            rgL[i] = rgL[i - 1]
            rgR[i] = rgR[i - 1]
            heelSN[i] = heelSN[i - 1]
            heelIN[i] = heelIN[i - 1]
        tlog[i] = t
        if fi < len(ftimes) and i == int(round(ftimes[fi] / DT_PHY)):
            grab.grab(d, i)
            fi += 1
        if not np.isfinite(d.qpos).all():
            print(f"[FAIL] m5 perturbed qpos non-finite at t={t:.3f} s")
            return
    onB = [t for t in sch["L"]["onsets"] if t >= 2.0]
    nxtA = [t for t in onA if t > tstar]
    nxtB = [t for t in onB if t > tstar]
    shifts = [b - a for a, b in zip(nxtA[:4], nxtB[:4])]
    print(f"[m5 perturbed] ground contacts {g_contacts} | frames {grab.n} | "
          f"heel-L max in {heelIN[:, 0].max():.2f} nA | first-onset shifts "
          f"{['%+.0f' % (s * 1e3) for s in shifts]} ms (report: +54/+12/+10/+12)")
    np.savez_compressed(
        os.path.join(TOOLS, "data_m5p.npz"),
        t=tlog, q=qlog, rgL=rgL, rgR=rgR, heelSN=heelSN, heelIN=heelIN,
        onL=np.array(sch["L"]["onsets"]), onR=np.array(sch["R"]["onsets"]),
        tstar=tstar, ftimes=ftimes, g_contacts=g_contacts)


def run_m6() -> None:
    """M6 harness-supported walk (rig S=1, acap 0.15, ibnA 0) via the gate's
    own run_phase; no frames (static figure)."""
    import test_w2l_ground as TG
    _STDOUT_KEEP.append(sys.stdout)   # keep every stdout wrapper alive
    TG.KV.update(dur_walk=DUR, rig=1.0, acap=0.15, ibnA=0.0)
    TG.write_ground_xml()
    log = TG.run_phase("walk", DUR, 1.0, None)
    w = log["t"] >= 1.0
    q = -log["q"][w]
    heelF = log["heelF"][w]
    duty = {s: float((heelF[:, k] > 5.0).mean()) for k, s in enumerate("LR")}
    stL = bursts_refractory(log["rgE"][w, 0], dt=DT_PHY, min_gap_s=0.4)
    per = float(np.diff(stL).mean() * DT_PHY) if len(stL) >= 3 else float("nan")
    keep = dict(t=log["t"], q=log["q"], com=log["com"], pel=log["pel"],
                tilt=log["tilt"], heelF=log["heelF"], toeF=log["toeF"],
                heelV=log["heelV"], toeV=log["toeV"], rgE=log["rgE"],
                rigFz=log["rigFz"], heelVmax=float(log["heelV"][w].max()))
    # save BEFORE printing (a closed stdout must never lose the data)
    np.savez_compressed(os.path.join(TOOLS, "data_m6.npz"), **keep)
    print(f"[m6] finite {log['finite']} fall_t {log['fall_t']} | "
          f"pel z min {log['pel'][w, 2].min():.3f} | "
          f"tilt max {log['tilt'][w].max():.1f} | "
          f"harness {100.0 * log['rigFz'][w].mean() / 411.0:+.0f}% | "
          f"heel duty L/R {duty['L']:.2f}/{duty['R']:.2f} | "
          f"heelF max {heelF.max():.1f} N | RG-E {len(stL)} @ {per:.3f} s")
    print("[m6] report: fall=no, pel z min 0.660, tilt 11.7, harness +39%, "
          "duty 0.00, heel SN 0.00 mV, 14 bursts @ 1.357 s")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else ""
    if which == "m3":
        run_air("m3")
    elif which == "m4":
        run_air("m4")
    elif which == "m5":
        run_m5()
    elif which == "m6":
        run_m6()
    else:
        print("usage: w2l_runs.py m3|m4|m5|m6")
        sys.exit(2)
