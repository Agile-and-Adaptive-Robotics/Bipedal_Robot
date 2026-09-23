"""Six-synergy neuromechanical model (Goal 7, 2026-09-23).

First-pass neuromechanical model driven by the six per-leg NMF synergies
previously identified in this project: ``fsa_backsolve.py`` selected a
shared rank of 6 NMF components per leg (smallest rank with centered
VAF >= 0.90; left needed 6, right 5) on ``bsolve_out.npz['acts']`` --
the ``bsolve_ik.py`` back-solved activations from subject01_walk1 IK +
measured GRF.  See fsa_results/fsa_backsolve_report.md and DESIGN.md
("the six NMF synergies ... S1--S6").

This module does NOT touch runner.py / params.py / build_network.py /
_curriculum.py (owned by another agent); it imports runner read-only for
MODEL path + patch/harness helpers.

Stages
------
1. ``basis``  - write the canonical W (muscle x synergy) / H (synergy x
   time) layout to ``synergy_basis.npz``, extracted VERBATIM from the
   existing k=6 basis of record (fsa_results/fsa_backsolve.npz
   ``{side}_gain`` / ``{side}_source``).  Nothing is re-factorized.
2. ``replay`` - replay the reference gait through ONLY the six synergy
   channels: W fixed, per-timestep NNLS solves the synergy coefficients
   onto the reference muscle activations ("H drives").  Reports
   uncentered VAF (synergy-literature convention) and centered R^2,
   overall and per functional muscle group (muscle_map primary group),
   plus the correlation between the NNLS coefficients and the stored
   NMF H traces.
3. ``demo``   - minimal MuJoCo OPEN-LOOP replay: muscle ctrl driven by
   the NNLS synergy coefficients (interpolated onto the 2 ms sim grid,
   2 s clip looped), suspended rig (no ground, pelvis pinned incl.
   rotation - runner.apply_harness(no_ground=True, pin_rot=True)).
   Standing is not attempted; the suspended configuration is chosen
   deliberately (open-loop activation replay has no balance closure).
   Unique timestamped npz output; joint tracking vs the reference IK
   is reported for the last clip loop.
4. ``note``   - print the walker-integration design note (synergy ->
   MN-pool mapping via W, six channels as phase-windowed drives).

Run (cwd = Code/MuJoCo_SNS/spinal, myo env)::

    set CONDA_PREFIX=C:\\Users\\Ben Bolen\\.conda\\envs\\myo
    C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe synergy_model.py
    # optional: --stage basis|replay|demo|note  --seconds 6
"""
from __future__ import annotations

import argparse
import datetime as _dt
import io
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

HERE = Path(__file__).parent
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

FSA_NPZ = HERE / "fsa_results" / "fsa_backsolve.npz"
FSA_REPORT = HERE / "fsa_results" / "fsa_backsolve_report.md"
BSOLVE_NPZ = HERE / "bsolve_out.npz"
BASIS_NPZ = HERE / "synergy_basis.npz"
SIDES = ("r", "l")
K = 6                      # synergies per leg (basis of record)
WARMUP_S = 0.5             # pose-hold while activation ramps in (s)
NAN_GUARD = True


def _vaf_unc(target: np.ndarray, pred: np.ndarray) -> float:
    """Uncentered VAF (Tresch/Cheung synergy convention)."""
    sse = float(np.sum((target - pred) ** 2))
    return 1.0 - sse / max(float(np.sum(target ** 2)), 1e-12)


def _r2(target: np.ndarray, pred: np.ndarray) -> float:
    """Centered R^2 (== 'centered VAF' in fsa_backsolve terms)."""
    sse = float(np.sum((target - pred) ** 2))
    sst = float(np.sum((target - target.mean()) ** 2))
    return 1.0 - sse / max(sst, 1e-12)


# ------------------------------------------------------------------ stage 1
def stage_basis() -> dict:
    """Extract the k=6 basis of record into the canonical W/H layout."""
    if not FSA_NPZ.exists():
        raise FileNotFoundError(
            f"{FSA_NPZ} not found - the six-synergy basis of record is "
            "missing; regenerate with fsa_backsolve.py (bsolve_out.npz "
            "must exist first, from bsolve_ik.py)")
    fsa = np.load(FSA_NPZ, allow_pickle=True)
    payload = {
        "provenance": (
            "k=6 per-leg NMF synergies extracted from "
            "fsa_results/fsa_backsolve.npz ({side}_gain/{side}_source), "
            "produced by fsa_backsolve.py on bsolve_out.npz['acts'] "
            "(bsolve_ik.py back-solve, subject01_walk1 IK + measured "
            "GRF, 121 frames 0.500-2.500 s). Rank 6 = smallest shared "
            "rank with centered VAF>=0.90 (left 6, right 5)."
        ),
        "time": fsa["time"],
    }
    report = {}
    for side in SIDES:
        gain = fsa[f"{side}_gain"]            # [K, M] spatial weights
        source = fsa[f"{side}_source"]        # [T, K] normalized coeffs
        names = np.asarray([str(x) for x in fsa[f"{side}_names"]])
        target = fsa[f"{side}_target"]        # [T, M] analysis target
        W = np.ascontiguousarray(gain.T)      # [M, K] muscle x synergy
        H = np.ascontiguousarray(source.T)    # [K, T] synergy x time
        rec = source @ gain                   # == H.T @ W.T
        payload[f"W_{side}"] = W
        payload[f"H_{side}"] = H
        payload[f"muscle_names_{side}"] = names
        payload[f"target_{side}"] = target
        report[side] = {
            "n_muscles": int(W.shape[0]),
            "r2_stored_H": _r2(target, rec),
            "vaf_stored_H": _vaf_unc(target, rec),
        }
    with open(BASIS_NPZ, "wb") as fh:
        np.savez_compressed(fh, **payload)
        fh.close()
    for side in SIDES:
        print(f"[basis] side {side}: W {payload[f'W_{side}'].shape}, "
              f"H {payload[f'H_{side}'].shape}, "
              f"{report[side]['n_muscles']} muscles, "
              f"stored-H replay R2={report[side]['r2_stored_H']:.4f} "
              f"VAF={report[side]['vaf_stored_H']:.4f}")
    print(f"[basis] wrote {BASIS_NPZ.name}")
    report["path"] = str(BASIS_NPZ)
    return report


# ------------------------------------------------------------------ stage 2
def _group_of(name: str) -> str:
    import muscle_map
    info = muscle_map.classify(name)
    return info.groups[0] if info else "unmapped"


def stage_replay() -> dict:
    """Per-timestep NNLS synergy replay + reconstruction quality."""
    basis = np.load(BASIS_NPZ, allow_pickle=True)
    raw = np.load(BSOLVE_NPZ, allow_pickle=True)
    raw_names = [str(x) for x in raw["act_names"]]
    raw_acts = np.asarray(raw["acts"], dtype=float)
    report = {"per_group": {}, "per_side": {}}
    all_t, all_p, all_tr = [], [], []
    for side in SIDES:
        W = basis[f"W_{side}"]                     # [M, K]
        H = basis[f"H_{side}"]                     # [K, T]
        names = [str(x) for x in basis[f"muscle_names_{side}"]]
        target = basis[f"target_{side}"]           # [T, M]
        T, M = target.shape

        # (a) stored-NMF-H replay (the basis as shipped)
        rec_h = H.T @ W.T
        # (b) per-timestep NNLS: min_{c>=0} ||W c - a(t)||^2 per t
        C = np.zeros((T, K))
        for i in range(T):
            C[i], _ = nnls(W, target[i])
        rec_n = C @ W.T

        raw_idx = [raw_names.index(n) for n in names]
        raw_side = raw_acts[:, raw_idx]

        side_rep = {
            "vaf_stored_H": _vaf_unc(target, rec_h),
            "r2_stored_H": _r2(target, rec_h),
            "vaf_nnls": _vaf_unc(target, rec_n),
            "r2_nnls": _r2(target, rec_n),
            "vaf_nnls_vs_raw": _vaf_unc(raw_side, rec_n),
            "r2_nnls_vs_raw": _r2(raw_side, rec_n),
        }
        # coefficient agreement: NNLS c(t) vs stored NMF H(t)
        corr = []
        for k in range(K):
            a, b = C[:, k], H[k, :]
            if a.std() > 1e-12 and b.std() > 1e-12:
                corr.append(float(np.corrcoef(a, b)[0, 1]))
            else:
                corr.append(float("nan"))
        side_rep["nnls_vs_H_corr"] = corr
        report["per_side"][side] = side_rep
        all_t.append(target)
        all_p.append(rec_n)
        all_tr.append(raw_side)

        # per-muscle then per-group (primary group)
        per_m_r2 = 1.0 - np.sum((target - rec_n) ** 2, axis=0) / \
            np.maximum(np.sum((target - target.mean(axis=0)) ** 2,
                              axis=0), 1e-12)
        per_m_vaf = 1.0 - np.sum((target - rec_n) ** 2, axis=0) / \
            np.maximum(np.sum(target ** 2, axis=0), 1e-12)
        groups = {}
        for j, name in enumerate(names):
            groups.setdefault(_group_of(name), []).append(j)
        for gname, idx in sorted(groups.items()):
            report["per_group"].setdefault(gname, {})[side] = {
                "n": len(idx),
                "r2_mean": float(np.mean(per_m_r2[idx])),
                "vaf_mean": float(np.mean(per_m_vaf[idx])),
            }
        print(f"\n[replay] side {side} ({M} muscles, {T} frames)")
        print(f"  stored-H replay : VAF {side_rep['vaf_stored_H']:.4f} "
              f"R2 {side_rep['r2_stored_H']:.4f}")
        print(f"  NNLS replay     : VAF {side_rep['vaf_nnls']:.4f} "
              f"R2 {side_rep['r2_nnls']:.4f}")
        print(f"  NNLS vs raw acts: VAF {side_rep['vaf_nnls_vs_raw']:.4f} "
              f"R2 {side_rep['r2_nnls_vs_raw']:.4f}")
        print("  NNLS-vs-H coefficient r: "
              + " ".join(f"S{k+1}={c:.3f}" for k, c in enumerate(corr)))

    pooled_t = np.hstack(all_t)
    pooled_p = np.hstack(all_p)
    pooled_raw = np.hstack(all_tr)
    report["overall"] = {
        "vaf_nnls": _vaf_unc(pooled_t, pooled_p),
        "r2_nnls": _r2(pooled_t, pooled_p),
        "vaf_nnls_vs_raw": _vaf_unc(pooled_raw, pooled_p),
        "r2_nnls_vs_raw": _r2(pooled_raw, pooled_p),
    }
    print(f"\n[replay] pooled both sides ({pooled_t.shape[1]} muscles): "
          f"VAF {report['overall']['vaf_nnls']:.4f} "
          f"R2 {report['overall']['r2_nnls']:.4f}; vs raw acts "
          f"VAF {report['overall']['vaf_nnls_vs_raw']:.4f} "
          f"R2 {report['overall']['r2_nnls_vs_raw']:.4f}")
    print("\n[replay] per functional group (primary), NNLS replay:")
    print(f"  {'group':10s} {'n_r':>3s} {'VAF_r':>6s} {'R2_r':>6s} "
          f"{'n_l':>3s} {'VAF_l':>6s} {'R2_l':>6s}")
    for gname, per in report["per_group"].items():
        r = per.get("r", {"n": 0, "vaf_mean": float("nan"),
                          "r2_mean": float("nan")})
        l = per.get("l", {"n": 0, "vaf_mean": float("nan"),
                          "r2_mean": float("nan")})
        print(f"  {gname:10s} {r['n']:3d} {r['vaf_mean']:6.3f} "
              f"{r['r2_mean']:6.3f} {l['n']:3d} {l['vaf_mean']:6.3f} "
              f"{l['r2_mean']:6.3f}")
    return report


# ------------------------------------------------------------------ stage 3
def stage_demo(seconds: float = 6.0, leg_damping: float = 0.5) -> dict:
    """Suspended (no-ground) open-loop MuJoCo replay of the NNLS synergy
    drive.  Standing/balance is NOT attempted - open-loop activation
    replay has no feedback path; the pelvis rig (no_ground=True,
    pin_rot=True) is the suspended configuration.

    ``leg_damping`` passes through to runner.apply_harness (repair 2d,
    the project's own air-stability knob).  0 = converted default 0.05.
    Without it the replay goes non-finite at ~1.8 s (rect_fem patella
    follower singular under unloaded deep-knee-flexion drive)."""
    import mujoco
    import runner                      # read-only: MODEL + helpers

    basis = np.load(BASIS_NPZ, allow_pickle=True)
    raw = np.load(BSOLVE_NPZ, allow_pickle=True)
    raw_names = [str(x) for x in raw["act_names"]]

    # NNLS coefficients + reconstruction per side (same solve as replay)
    ctrl_side, included = {}, {}
    for side in SIDES:
        W = basis[f"W_{side}"]
        target = basis[f"target_{side}"]
        T = target.shape[0]
        C = np.zeros((T, K))
        for i in range(T):
            C[i], _ = nnls(W, target[i])
        ctrl_side[side] = {"t": basis["time"], "A": C @ W.T,
                           "names": [str(x) for x in
                                     basis[f"muscle_names_{side}"]]}
        included[side] = set(ctrl_side[side]["names"])

    t_a = ctrl_side["r"]["t"]
    period = float(t_a[-1] - t_a[0])

    # ---- model: suspended rig (runner's own --no-ground configuration)
    m0 = mujoco.MjModel.from_xml_path(str(runner.MODEL))
    d0 = mujoco.MjData(m0)
    mujoco.mj_resetDataKeyframe(m0, d0, 0)
    model = runner.apply_harness(m0, d0, no_ground=True,
                                 pin_rot=True, leg_damping=leg_damping
                                 if leg_damping else None)
    data = mujoco.MjData(model)
    runner.apply_start_pose(model, data)

    aid = {mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i): i
           for i in range(model.nu)}
    missing = [n for n in raw_names if n not in aid]
    if missing:
        raise RuntimeError(f"actuators missing from model: {missing}")
    zeroed = [n for n in raw_names
              if n not in included["r"] and n not in included["l"]]

    # per-actuator interp sources (None -> ctrl 0)
    src = {}
    for side in SIDES:
        A = ctrl_side[side]["A"]
        for j, name in enumerate(ctrl_side[side]["names"]):
            src[name] = A[:, j]

    jadr = {jn: model.joint(jn).qposadr[0] for jn in runner.KEY_JOINTS}
    pelvis_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY,
                                   "pelvis")
    ref_q = np.degrees(raw["qpos"][:, [jadr[jn]
                                       for jn in runner.KEY_JOINTS]])

    dt = float(model.opt.timestep)
    n_steps = int(round(seconds / dt))
    log_t = np.zeros(n_steps)
    log_q = np.zeros((n_steps, len(runner.KEY_JOINTS)))
    log_com = np.zeros((n_steps, 3))
    log_ctrl = np.zeros((n_steps, model.nu))
    broke_at = None

    for k in range(n_steps):
        t = k * dt
        # looped clip time on the analysis grid
        tau = (t - WARMUP_S) % period
        ta = t_a[0] + tau
        ramp = min(max(t / WARMUP_S, 0.0), 1.0)
        data.ctrl[:] = 0.0
        for name, col in src.items():
            data.ctrl[aid[name]] = ramp * float(
                np.interp(ta, t_a, col))
        if t < WARMUP_S:
            data.qvel[:] = 0.0            # runner's warmup pattern
            mujoco.mj_forward(model, data)
        else:
            mujoco.mj_step(model, data)
        if NAN_GUARD and not (np.all(np.isfinite(data.qacc))
                              and np.all(np.isfinite(data.qvel))):
            broke_at = t
            print(f"[demo] !! non-finite qacc/qvel at t={t:.3f} - "
                  "stopping log here")
            n_steps = k + 1
            log_t, log_q = log_t[:n_steps], log_q[:n_steps]
            log_com, log_ctrl = log_com[:n_steps], log_ctrl[:n_steps]
            break
        log_t[k] = t
        log_com[k] = data.subtree_com[pelvis_id]
        for i, jn in enumerate(runner.KEY_JOINTS):
            log_q[k, i] = data.qpos[jadr[jn]]
        if k % int(0.5 / dt) == 0:
            print(f"[demo] t={t:6.2f} com_z={log_com[k, 2]:.3f} "
                  f"knee_r={np.degrees(log_q[k, 4]):+7.1f} "
                  f"knee_l={np.degrees(log_q[k, 9]):+7.1f}")

    # ---- tracking vs reference IK on the LAST full clip loop
    track = {}
    if n_steps * dt > WARMUP_S + period and broke_at is None:
        sl = slice(int((log_t[-1] - period) / dt), n_steps)
        tt = log_t[sl]
        ta = t_a[0] + ((tt - WARMUP_S) % period)
        ref_i = np.clip(np.searchsorted(t_a, ta), 0, len(t_a) - 1)
        ref = ref_q[ref_i]
        sim_deg = np.degrees(log_q[sl])
        for i, jn in enumerate(runner.KEY_JOINTS):
            if jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
                      "hip_flexion_l", "knee_angle_l", "ankle_angle_l"):
                rmse = float(np.sqrt(np.mean(
                    (sim_deg[:, i] - ref[:, i]) ** 2)))
                a, b = sim_deg[:, i], ref[:, i]
                r = float(np.corrcoef(a, b)[0, 1]) \
                    if a.std() > 1e-9 and b.std() > 1e-9 else float("nan")
                track[jn] = {"rmse_deg": rmse, "r": r}
        print("[demo] open-loop tracking vs reference IK "
              f"(last {period:.2f} s loop):")
        for jn, m in track.items():
            print(f"  {jn:16s} RMSE {m['rmse_deg']:6.1f} deg  "
                  f"r={m['r']:+.3f}")

    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = HERE / f"synergy_replay_air_{stamp}.npz"
    with open(out, "wb") as fh:
        np.savez_compressed(fh, t=log_t, q=np.degrees(log_q),
                            com=log_com, ctrl=log_ctrl,
                            key_joints=runner.KEY_JOINTS,
                            ref_q_deg=ref_q,
                            zeroed_actuators=np.array(zeroed, dtype=object),
                            cfg="air-openloop-synergy6",
                            leg_damping=leg_damping)
        fh.close()
    status = f"broke at t={broke_at:.2f}" if broke_at is not None \
        else "clean (no non-finite states)"
    print(f"[demo] wrote {out.name} ({n_steps} steps, "
          f"{n_steps * dt:.2f} s, {status})")
    print(f"[demo] zero-ctrl actuators (pruned, outside basis): "
          f"{', '.join(zeroed)}")
    rom = {jn: (float(np.degrees(log_q[:, i].min())),
                float(np.degrees(log_q[:, i].max())))
           for i, jn in enumerate(runner.KEY_JOINTS)}
    return {"npz": out.name, "steps": n_steps, "broke_at": broke_at,
            "seconds": n_steps * dt, "tracking": track,
            "rom_deg": rom, "zeroed": zeroed}


# ------------------------------------------------------------------ stage 4
NOTE = """DESIGN NOTE - wiring the six synergies into the walker
---------------------------------------------------------
Route (from the s3k decision-tree verdict, DESIGN.md 2026-09-21/22:
output/input gating levers are exhausted; the PATTERN LAYER must
change - 'motor primitives with antagonist groups', with the pm_*
contact-reset phase machine as the ready-made phase source):

1. Basis W (synergy_basis.npz, per leg [43 muscles x 6]): each column
   is the MN-pool weighting of one synergy channel. In SNS terms a
   synergy channel becomes a PF-population -> MN-pool projection:
   PF source neuron(s) -> one MN pool per W row with W > 0. The
   conductance mapping V_MN = E_HI * a with g = k*R*Gm/(dE - k*R)
   (Szczecinski et al. 2017 Eq. 18) is already implemented and
   audited in fsa_backsolve.py ({side}_g_fit / {side}_g_analytic).

2. Six channels as PHASE-WINDOWED drives: H (synergy x time) is the
   measured drive waveform per channel. Cycle-fold H (heel strike=0,
   toe-off=50%) and replace each channel with a raised-cosine window
   fitted to its phase profile (Di Russo et al. 2023 eq. 8 shape);
   the pm_* phase machine already provides phase 0-1 per side. This
   is exactly the 'phase-driven signed primitives' build spec the
   decision tree recommends (LIT_CIRCUIT_AUDIT.md section 10),
   except the weight matrix is the MEASURED W rather than searched.

3. What this demo establishes: the six channels alone reconstruct
   the reference activations at the reported VAF, and an open-loop
   suspended replay produces alternating-sagittal motion without
   any neural rhythm generator. What it does NOT establish: ground
   gait, balance, or closed-loop stability (no feedback path); the
   activation-level ctrl bypasses MN membrane dynamics and the
   afferent circuit entirely.

Open items for Ben: (a) whether synergy channels map 1:1 to PF
populations or pairs of them share half-centers (DESIGN.md's
'phase/PF semantics' warning: a synergy is not a PF neuron);
(b) S5/S6 are not bilaterally phase-stable - keep or drop per side;
(c) whether to phase-window H per channel or replay measured H.
"""


def stage_note() -> str:
    print(NOTE)
    return NOTE


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", default="all",
                    choices=["all", "basis", "replay", "demo", "note"])
    ap.add_argument("--seconds", type=float, default=6.0,
                    help="open-loop demo duration (s)")
    ap.add_argument("--leg-damping", type=float, default=0.5,
                    help="leg-hinge damping for the suspended demo "
                         "(runner repair 2d; 0.5 default = the clean-run "
                         "configuration; 0 = converted default 0.05, which "
                         "goes non-finite at ~1.8 s)")
    args = ap.parse_args(argv)

    if args.stage in ("all", "basis"):
        stage_basis()
    if args.stage in ("all", "replay"):
        stage_replay()
    if args.stage in ("all", "demo"):
        stage_demo(args.seconds, args.leg_damping)
    if args.stage in ("all", "note"):
        stage_note()
    print("\nDone.")


if __name__ == "__main__":
    main()
