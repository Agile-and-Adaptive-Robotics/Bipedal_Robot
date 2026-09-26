r"""W2L AIR GATE — milestones 3 AND 4.

M3 (net=orig, default): the 2023 ORIGINAL single-LH-RG architecture
(build_w2l_orig_net). M4 (net=split): Ben's ask — the L-RG->R-PF crossed
drive edges REMOVED, a mirrored R RG powering the R PF layers ipsilaterally,
the two RGs coupled by the documented bilateralrg commissurals (c1 mutual
inhibition + weak V3 excitation; build_w2l_split_net). --comm=0 ablates the
commissurals (conditional topology: the 4 c1/V3 neurons + 10 synapses are NOT
built), proving each leg is powered by its OWN RG.

Run:
    C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_w2l_air.py [knobs]
    knobs: --net=orig|split  --comm=1.0|0.0  --dur=N --te=N --tf=N --tau=N
           --cap=N --damp=N --stiff=N --lift=N

Body: w2l_air.xml (written by this script) = w2l_mjcf_fixed.xml with the
welded Root lifted +0.30 m -> the pelvis is FIXED at standing height and no
leg pose can reach the ground plane (the ask's "pelvis held + no ground
contact" rig; also how the AnimatLab standalone .asim air-stepping reference
was measured). w2l_mjcf_fixed.xml = w2l_mjcf.xml with the 8 hinge axes
corrected (fix_joint_axes.py; the shipped knee/ankle axes were vertical) and
the knee ranges negated to [-60 deg, 0] (flexion = -qpos; see that module).

Neural drive regime: verbatim kickoff (Stimulus_1 10 nA 10 ms -> L RG ext;
split mode adds Stimulus_2 10 nA 10 ms -> R RG flx, the antiphase kickoff of
SESSION_NOTES build_rg.pl) + tonic te/tf on the RG half-centers (tonic-F is
the escape engine; split mode drives BOTH RGs symmetrically).

VERDICTS
  M3 (net=orig): PASS requires a finite >=20 s run, ZERO ground contacts,
      both hips swinging >= 15 deg in antiphase, period within 2x of a
      reference (2.22 Hz modern W2L air / 0.77 Hz 2023 original).
  M4 coupled (net=split --comm>0): the M3 checks PLUS rg_e_antiphase (L/R
      RG-E Pearson r < -0.5, smoke_w2l convention) PLUS phase_half_cycle:
      the peak positive correlation of hip L vs R flexion sits at
      0.35..0.65 of a cycle (antiphase = half-cycle lag).
  M4 ablated (net=split --comm=0): finite run, zero ground contacts, and
      BOTH legs oscillate autonomously: per side >= 4 RG-E bursts with a
      measurable period, >= 8 hip swing excursions, hip range >= 15 deg,
      finite ACF dominant period. Antiphase is NOT required (free-running
      twin oscillators may drift); measured periods are reported.
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import mujoco  # noqa: E402

FIXED = os.path.join(HERE, "w2l_mjcf_fixed.xml")
AIR = os.path.join(HERE, "w2l_air.xml")

# ------------------------------------------------------------------ knobs
KV = dict(net="orig", comm=1.0,
          te=3.0, tf=4.0, tau=0.25, cap=0.5, damp=3.0, stiff=1.0,
          dur=20.0, lift=0.30)
for arg in sys.argv[1:]:
    if arg.startswith("--") and "=" in arg:
        k, v = arg[2:].split("=", 1)
        try:
            v = float(v)
        except ValueError:
            pass
        KV[k.replace("-", "_")] = v

SPLIT = KV["net"] == "split"
ABLATED = SPLIT and KV["comm"] <= 0.0

DT_PHY = 0.001
NSUB = 2                      # physics steps per 2 ms net step (net dt = DT)
REFS = dict(modern_hz=2.22, orig2023_hz=0.77,
            hip_deg=38.0, knee_deg=61.0, ankle_deg=16.0)


def write_air_xml() -> None:
    txt = open(FIXED, encoding="utf-8").read()
    old = 'pos="-3.454 0 0.99298"'
    assert old in txt, "root pos pattern not found in w2l_mjcf_fixed.xml"
    hdr = ("<!-- AIR VARIANT for milestones 3/4 (test_w2l_air.py): identical\n"
           "     to w2l_mjcf_fixed.xml (axis-fixed body, see fix_joint_axes.py)\n"
           "     except the welded Root is lifted +0.30 m so NO leg pose can\n"
           "     reach the ground plane (leg vertical reach ~0.95 m). Keeps\n"
           "     'no ground contact' provable (ground contacts == 0). -->\n")
    new = 'pos="-3.454 0 %.5f"' % (0.99298 + KV["lift"])
    open(AIR, "w", encoding="utf-8").write(hdr + txt.replace(old, new))


def swing_peaks(sig: np.ndarray, min_gap_ms: int = 400):
    """Onsets of sustained flexion excursions (rising through the mean + 0.3
    std, with a refractory gap so threshold ripples don't double-count)."""
    thr = sig.mean() + 0.3 * sig.std()
    on = sig > thr
    starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    if len(starts) == 0:
        return starts
    keep = [starts[0]]
    for s in starts[1:]:
        if s - keep[-1] > min_gap_ms:
            keep.append(s)
    return np.array(keep)


def bursts_refractory(sig: np.ndarray, dt: float, min_gap_s: float = 0.4):
    """RG-E burst onsets (rising through 0.5 max) with a refractory gap."""
    on = sig > 0.5 * sig.max()
    starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    if len(starts) == 0:
        return starts
    keep = [starts[0]]
    for s in starts[1:]:
        if (s - keep[-1]) * dt > min_gap_s:
            keep.append(s)
    return np.array(keep)


def acf_period(x: np.ndarray) -> float:
    """Dominant oscillation period = lag of the first autocorrelation
    peak beyond lag 0.25 s (robust to waveform distortion)."""
    x = x - x.mean()
    n = len(x)
    ac = np.correlate(x, x, "full")[n - 1:]
    ac /= ac[0]
    lo = int(0.25 / DT_PHY)
    pk = ac[lo:]
    if pk.max() < 0.1:
        return float("nan")
    return float((int(np.argmax(pk)) + lo) * DT_PHY)


def xcorr_phase(a: np.ndarray, b: np.ndarray, tmax: float = 2.0):
    """Cross-correlation r(lag) of a vs b for lags in [-tmax, +tmax] (5 ms
    grid). Returns (lag_of_max, lag_of_min, r_at_max, r_at_min); positive lag
    means b leads a by |lag| (b(t+lag) aligned to a(t))."""
    a = a - a.mean()
    b = b - b.mean()
    n = len(a)
    ac = np.correlate(a, b, "full")          # index n-1+k -> b(t+k) vs a(t)
    lags = (np.arange(len(ac)) - (n - 1)) * DT_PHY
    m = np.abs(lags) <= tmax
    ac, lags = ac[m], lags[m]
    std = (np.sqrt((a * a).mean()) * np.sqrt((b * b).mean()))
    r = ac / (len(a) * std)
    i_max, i_min = int(np.argmax(r)), int(np.argmin(r))
    return float(lags[i_max]), float(lags[i_min]), float(r[i_max]), float(r[i_min])


def band_limited(x: np.ndarray, lo: float = 0.3, hi: float = 1.2) -> np.ndarray:
    """Fundamental-band version of a trace (4th-order Butterworth, zero
    phase). M4's split-RG runs carry MORE leg-leg midline-bump impulses than
    M3 (26741 vs 22632 in the first coupled run) and each bump kicks both
    hips IN PHASE through the welded pelvis (an in-phase 2nd harmonic —
    hip L excursion interval T/2 exposes it). The gait-level antiphase
    question is about the stepping fundamental, so split mode gates on the
    band-limited correlation; the raw one is printed alongside."""
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [lo, hi], btype="band", fs=1.0 / DT_PHY, output="sos")
    return sosfiltfilt(sos, x)


def main() -> int:
    if SPLIT:
        import build_w2l_split_net as B
        B.NAP["tau_max_h"] = KV["tau"]
        net = B.build(comm=KV["comm"])
    else:
        import build_w2l_orig_net as B
        B.NAP["tau_max_h"] = KV["tau"]
        net = B.build()

    write_air_xml()
    m = mujoco.MjModel.from_xml_path(AIR)
    d = mujoco.MjData(m)

    # runtime stand-ins, all documented deviations (M1 report section 3):
    # AnimatLab LinearHill muscles carry B = 400-800 N s/m; MuJoCo 2.3.7
    # <muscle> has no damping. Joint damping + capped ctrl + stiffer joint
    # limits keep the undamped hinges inside their (soft-default) limits.
    m.dof_damping[:] = KV["damp"]
    if KV["stiff"] > 0.0:
        m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
        m.jnt_solref[:] = np.array([0.006, 1.0])
    CAP = KV["cap"]

    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
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

    nphy = int(KV["dur"] / DT_PHY)
    tlog = np.zeros(nphy)
    qlog = np.zeros((nphy, 6))            # hipL kneeL anklL hipR kneeR anklR
    rgE = np.zeros(nphy)
    rgR = np.zeros(nphy) if SPLIT else None
    g_contacts = 0                        # contacts involving the ground plane
    other_contacts = 0                    # leg-leg self contacts
    for i in range(nphy):
        t = i * DT_PHY
        if i % NSUB == 0:
            u[iE] = KV["te"]
            u[iF] = KV["tf"]
            u[iS] = 10.0 if t < 0.01 else 0.0     # verbatim kickoff window
            if SPLIT:
                u[iEr] = KV["te"]
                u[iFr] = KV["tf"]
                u[iS2] = 10.0 if t < 0.01 else 0.0  # antiphase kickoff
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
            rgE[i] = V[net.idx["L RG ext"]]
            if SPLIT:
                rgR[i] = V[net.idx["R RG ext"]]
        else:
            if SPLIT:
                rgR[i] = rgR[i - 1]
        rgE[i] = rgE[i] if i % NSUB == 0 else rgE[i - 1]
        tlog[i] = t
        if not np.isfinite(d.qpos).all():
            print(f"[FAIL] qpos non-finite at t={t:.3f} s")
            return 1

    # ------------------------------------------------------------- metrics
    an = slice(int(2.0 / DT_PHY), None)    # drop the kickoff transient
    flex = -qlog[an]                       # flexion-positive convention
    rg = rgE[an]
    tt = tlog[an]

    tag_mode = ("M4 ABLATED (--comm=0, commissurals REMOVED)" if ABLATED else
                "M4 split-RG coupled" if SPLIT else
                "M3 W2L 2023 ORIGINAL (single LH RG)")
    print(f"== gate: {tag_mode} air stepping on the axis-fixed M1 body, "
          f"{KV['dur']:.0f} s ==")
    print(f"   knobs: net={KV['net']} comm={KV['comm']} te={KV['te']} "
          f"tf={KV['tf']} tau_rg_nap_h={KV['tau']} ctrl_cap={KV['cap']} "
          f"joint_damp={KV['damp']} stiff_limits={KV['stiff']:.0f} "
          f"lift={KV['lift']}")
    print(f"   ground contacts: {g_contacts} (must be 0) | "
          f"leg-leg self contacts: {other_contacts}")

    # neural periods (RG ext bursts, check_rhythm convention + refractory).
    # rg is sampled per PHYSICS step (held between net updates) -> dt = DT_PHY.
    starts = bursts_refractory(rg, dt=DT_PHY, min_gap_s=0.4)
    rg_per = float(np.diff(starts).mean() * DT_PHY) if len(starts) >= 3 else float("nan")
    print(f"   neural: L RG ext bursts={len(starts)}, period {rg_per:.3f} s "
          f"({1.0 / rg_per:.2f} Hz), E max {rg.max():.2f} mV")
    rg_per_R = float("nan")
    rg_r = float("nan")
    if SPLIT:
        rgr = rgR[an]
        startsR = bursts_refractory(rgr, dt=DT_PHY, min_gap_s=0.4)
        rg_per_R = (float(np.diff(startsR).mean() * DT_PHY)
                    if len(startsR) >= 3 else float("nan"))
        rg_r = float(np.corrcoef(rg[::5], rgr[::5])[0, 1])
        print(f"   neural: R RG ext bursts={len(startsR)}, period {rg_per_R:.3f} s "
              f"({1.0 / rg_per_R:.2f} Hz), E max {rgr.max():.2f} mV; "
              f"L/R RG-E r = {rg_r:+.3f} (antiphase want < -0.5)")

    # body: swing periods from hip-flexion excursions + ACF dominant period
    per = {}
    for side, col in (("L", 0), ("R", 3)):
        pk = swing_peaks(flex[:, col])
        per[side] = float(np.diff(pk).mean() * DT_PHY) if len(pk) >= 3 else float("nan")
        ap = acf_period(flex[:, col])
        if np.isfinite(per[side]):
            print(f"   hip {side} swing excursions: {len(pk)}, "
                  f"mean interval {per[side]:.3f} s "
                  f"({1.0 / per[side]:.2f} Hz); ACF dominant period {ap:.3f} s "
                  f"({1.0 / ap:.2f} Hz)")
        else:
            print(f"   hip {side} swing excursions: {len(pk)} (too few); "
                  f"ACF dominant period {ap:.3f} s")

    ds = slice(None, None, 5)   # 5 ms samples
    a, b = flex[:, 0][ds] - flex[:, 0][ds].mean(), flex[:, 3][ds] - flex[:, 3][ds].mean()
    r = float(np.corrcoef(a, b)[0, 1])

    # phase lag: peak positive correlation of L vs R hip flexion, in cycles
    # of the L hip's dominant period (antiphase pairs peak at half a cycle).
    phase = float("nan")
    r_band = float("nan")
    if SPLIT:
        T = acf_period(flex[:, 0])
        lag_pos, lag_neg, r_pos, r_neg = xcorr_phase(flex[:, 0], flex[:, 3],
                                                     tmax=2.0)
        if np.isfinite(T) and T > 0:
            phase = (lag_pos % T) / T
        fb = band_limited(flex[:, 0])
        rb = band_limited(flex[:, 3])
        r_band = float(np.corrcoef(fb, rb)[0, 1])
        print(f"   hip L/R cross-correlation: peak {r_pos:+.3f} at lag "
              f"{lag_pos:+.3f} s, trough {r_neg:+.3f} at lag {lag_neg:+.3f} s "
              f"(L hip ACF period {T:.3f} s)")
        print(f"   L/R phase lag = {phase:.3f} cycle (antiphase = 0.5)")
        print(f"   hip L/R correlation, fundamental band 0.3-1.2 Hz: "
              f"r = {r_band:+.3f} (antiphase want < -0.3); raw r = {r:+.3f} "
              "(contact-kick 2nd harmonic inflates the raw value)")

    print("   joint flexion-positive excursions (deg, min..max and range) "
          "vs AnimatLab references:")
    ranges = {}
    for name, col, ref in (("hip", (0, 3), REFS["hip_deg"]),
                           ("knee", (1, 4), REFS["knee_deg"]),
                           ("ankle", (2, 5), REFS["ankle_deg"])):
        rl = flex[:, col[0]].max() - flex[:, col[0]].min()
        rr = flex[:, col[1]].max() - flex[:, col[1]].min()
        ranges[name] = (rl, rr)
        print(f"     {name:<5} L [{flex[:, col[0]].min():+6.1f},{flex[:, col[0]].max():+6.1f}]"
              f" range {rl:5.1f} | R [{flex[:, col[1]].min():+6.1f},{flex[:, col[1]].max():+6.1f}]"
              f" range {rr:5.1f}   (ref range ~{ref:.0f})")

    print(f"   hip flexion L/R correlation (0 lag): r = {r:+.3f} "
          "(antiphase want < -0.3)")

    # ------------------------------------------------------------- verdict
    finite = True
    both_step = all(np.isfinite(per[s]) and len(swing_peaks(flex[:, c])) >= 8
                    for s, c in (("L", 0), ("R", 3)))
    amp_ok = all((flex[:, c].max() - flex[:, c].min()) >= 15.0
                 for c in (0, 3))
    anti = (r < -0.3)
    airborne = (g_contacts == 0)
    acf_hz = [1.0 / acf_period(flex[:, c]) for c in (0, 3)]

    if SPLIT:
        # ---------------- M4 verdicts ----------------
        if ABLATED:
            # CAUSALITY ABLATION: commissurals removed -> each leg must still
            # oscillate on its OWN RG. Per side: >= 4 RG-E bursts with a
            # measurable period, >= 8 hip excursions, hip range >= 15 deg,
            # finite ACF period.
            def side_ok(col, rgv):
                st = bursts_refractory(rgv, dt=DT_PHY, min_gap_s=0.4)
                return (len(st) >= 4
                        and np.isfinite(np.diff(st).mean() * DT_PHY)
                        and len(swing_peaks(flex[:, col])) >= 8
                        and (flex[:, col].max() - flex[:, col].min()) >= 15.0
                        and np.isfinite(acf_period(flex[:, col])))
            okL = side_ok(0, rg)
            okR = side_ok(3, rgR[an])
            checks = dict(finite_20s=finite, airborne=airborne,
                          leg_L_oscillates=okL, leg_R_oscillates=okR)
            for k, v in checks.items():
                print(f"   [{'PASS' if v else 'FAIL'}] {k}")
            print(f"   measured periods: L RG-E {rg_per:.3f} s / "
                  f"R RG-E {rg_per_R:.3f} s | hip-L ACF "
                  f"{acf_period(flex[:, 0]):.3f} s / hip-R ACF "
                  f"{acf_period(flex[:, 3]):.3f} s")
            if all(checks.values()):
                print(f"VERDICT: ABLATION PASS  both legs oscillate with the "
                      f"commissural coupling REMOVED "
                      f"(L {rg_per:.3f} s, R {rg_per_R:.3f} s) -> each leg is "
                      f"powered by its own RG")
                return 0
            print("VERDICT: ABLATION FAIL "
                  + ", ".join(k for k, v in checks.items() if not v))
            return 1

        # coupled: M3 checks + RG-E antiphase + half-cycle phase lag.
        # antiphase gates on the FUNDAMENTAL-BAND body correlation (raw
        # zero-lag r is inflated by the in-phase leg-leg contact kicks, see
        # band_limited); the raw value is printed alongside.
        rg_anti = (rg_r < -0.5)
        anti_band = (r_band < -0.3)
        phase_ok = np.isfinite(phase) and 0.35 <= phase <= 0.65
        cand_hz = [h for h in acf_hz + ([1.0 / rg_per] if np.isfinite(rg_per) else [])
                   if np.isfinite(h) and h > 0]
        freq_ok = any(0.5 * ref <= hz <= 2.0 * ref
                      for hz in cand_hz
                      for ref in (REFS["modern_hz"], REFS["orig2023_hz"]))
        checks = dict(finite_20s=finite, both_hips_swing=both_step,
                      antiphase=anti_band, airborne=airborne,
                      freq_within_2x=freq_ok, hip_amp_ok=amp_ok,
                      rg_e_antiphase=rg_anti, phase_half_cycle=phase_ok)
        for k, v in checks.items():
            print(f"   [{'PASS' if v else 'FAIL'}] {k}")
        if all(checks.values()):
            print(f"VERDICT: PASS  hip-ACF period {1.0 / acf_hz[0]:.3f} s "
                  f"({acf_hz[0]:.2f} Hz)  band-antiphase_r={r_band:.3f} "
                  f"(raw {r:+.3f})  phase={phase:.3f} cycle  "
                  f"hip_range_L={ranges['hip'][0]:.1f} deg")
            return 0
        print("VERDICT: FAIL " + ", ".join(k for k, v in checks.items() if not v))
        return 1

    # ---------------- M3 verdict (net=orig; unchanged logic) ----------------
    cand_hz = [h for h in acf_hz + ([1.0 / rg_per] if np.isfinite(rg_per) else [])
               if np.isfinite(h) and h > 0]
    freq_ok = any(0.5 * ref <= hz <= 2.0 * ref
                  for hz in cand_hz
                  for ref in (REFS["modern_hz"], REFS["orig2023_hz"]))
    checks = dict(finite_20s=finite, both_hips_swing=both_step, antiphase=anti,
                  airborne=airborne, freq_within_2x=freq_ok, hip_amp_ok=amp_ok)
    for k, v in checks.items():
        print(f"   [{'PASS' if v else 'FAIL'}] {k}")
    if all(checks.values()):
        print(f"VERDICT: PASS  period~{acf_hz[0]:.3f}s (hip L ACF)  "
              f"antiphase_r={r:.3f}  "
              f"hip_range_L={ranges['hip'][0]:.1f} deg")
        return 0
    print("VERDICT: FAIL " + ", ".join(k for k, v in checks.items() if not v))
    return 1


if __name__ == "__main__":
    sys.exit(main())
