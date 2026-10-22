"""Co-contraction stand descent experiment (Ben's 2026-10-21 directive).

Ben: "Massively activated antagonistic muscle groups should free-stand -
the trick is doing it energy efficiently. So heavily activate antagonists
so the torques balance (legs become rigid struts), THEN reduce from
there." The prior s3q study penalized co-contraction BEFORE proving the
capability - backwards. This script runs the DIRECT experiment:

  1. EXISTENCE PROOF: drive ALL antagonist muscle pairs bilaterally at
     high co-activation through the per-side ALPHA hand-drive cells
     (uniform MN-pool drive, extensors AND flexors alike; runner.py
     --alpha-stand, new 2026-10-21, default-off) at rig_scale 0 (rig
     OFF). Does the model free-stand with the legs as struts? If it
     buckles even at MAX co-contraction, report the failing joint.
  2. DESCENT CURVE: from the standing lock, reduce the co-contraction
     drive in steps and record stands/falls, feet-load %BW, COM height,
     sway, and mean summed muscle activation (the energy proxy). The
     minimum level that still stands = the efficient operating point.

Baseline config per family = its QUIET stage-3 winner (the s3p trial-0
seed recipe: v10 multipliers + renshaw + quiet winner params, mechanism
knobs 0) - the same incumbent the verdict table reports.

Usage (myo env, cwd = spinal):
  python _coact_stand_20261021.py matrix <stock|syn6|w2lvar>   # phases 1+2
  python _coact_stand_20261021.py stair  <stock|syn6|w2lvar>   # staircase run
  python _coact_stand_20261021.py one <family> <alpha_to_mn> <current_nA> <rig> [dur]
      # single fresh-process run (verification rungs); prints one JSON line
Results JSON: reports_20261020/balance/coact_stand_<family>_20261021.json
"""
import io
import importlib
import json
import os
import sys
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

MODE = sys.argv[1] if len(sys.argv) > 1 else "matrix"
FAMILY = sys.argv[2] if len(sys.argv) > 2 else "stock"

CFG = {
    "stock": dict(mod="_curriculum", quiet_json="curriculum_stage3.json",
                  net=None),
    "syn6": dict(mod="_curriculum_syn6",
                 quiet_json="curriculum_syn6_stage3.json", net="syn6"),
    "w2lvar": dict(mod="_curriculum_w2lvar",
                   quiet_json="curriculum_w2lvar_stage3.json", net="w2lvar"),
}[FAMILY]

OUTDIR = os.path.join("reports_20261020", "balance")
NPZ = f"_coactstand_{FAMILY}_20261021.npz"

os.environ["AARL_NPZ"] = NPZ
if CFG["net"]:
    os.environ["AARL_NET"] = CFG["net"]
else:
    os.environ.pop("AARL_NET", None)
# incumbent pose: mechanism offsets OFF (quiet winners did not use them)
for _e in ("AARL_STANCE_W", "AARL_KNEE_BIAS", "AARL_PELVIS_TY"):
    os.environ.pop(_e, None)

import numpy as np  # noqa: E402
import runner as R  # noqa: E402  (also pins AARL_MODEL default)
import params as P  # noqa: E402
from muscle_map import classify  # noqa: E402

# body weight for %BW numbers (unpatched model masses == patched)
import mujoco  # noqa: E402
_m = mujoco.MjModel.from_xml_path(str(R.MODEL))
BW = float(_m.body_mass.sum()) * 9.81
del _m

ANTAG_PAIRS = (("hip_flex", "hip_ext"), ("knee_flex", "knee_ext"),
               ("ankle_df", "ankle_pf"), ("trunk_flex", "trunk_ext"))
KEY_JOINTS = R.KEY_JOINTS


def load_incumbent():
    C = importlib.import_module(CFG["mod"])
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    base = dict(prev["multipliers"])
    base["renshaw"] = 0.5
    qw = json.loads(open(CFG["quiet_json"], encoding="utf-8").read())["params"]
    p = {**base, **qw}
    C.set_stage(3, p)
    P.G["coact"] = 0.0
    return p


def run_one(alpha_to_mn, current, rig=0.0, dur=8.0, tag=""):
    """One standing eval; returns (metrics, analysis-dict). alpha_to_mn
    None = no --alpha-stand flag at all (pure incumbent reference)."""
    if alpha_to_mn is not None:
        P.G["alpha_to_mn"] = float(alpha_to_mn)
    npz = f"_coactstand_{FAMILY}_{tag}_20261021.npz"
    os.environ["AARL_NPZ"] = npz
    t0 = time.time()
    args = ["--stand-eval", repr(float(dur)), "--rig-scale", repr(float(rig))]
    if current is not None:
        args += ["--alpha-stand", repr(float(current))]
    m = R.main(args)
    wall = time.time() - t0
    a = analyze_npz(npz, m, tag or f"g{alpha_to_mn}_I{current}_rig{rig}")
    a["wall_s"] = round(wall, 1)
    a["alpha_to_mn"] = alpha_to_mn
    a["alpha_nA"] = current
    a["rig_scale"] = rig
    return m, a


def analyze_npz(npz_path, m, tag):
    """Steady-window (t>=4 s) analysis of the saved run: per-antagonist-
    group mean activation, energy proxy, joint envelopes, and buckling
    forensics if the run fell."""
    d = np.load(npz_path)
    t, com, q = d["t"], d["com"], d["q"]          # q in deg, KEY_JOINTS
    ctrl, contact = d["ctrl"], d["contact"]
    n = int(np.sum(t > 0)) + 1
    n = min(n, len(t))
    t, com, q, ctrl = t[:n], com[:n], q[:n], ctrl[:n]
    contact = contact[:n]
    i4 = min(int(np.searchsorted(t, 4.0)), n - 1)
    seg = slice(i4, n)
    # actuator names in model order + group means (primary membership)
    import mujoco as _mj
    _mm = _mj.MjModel.from_xml_path(str(R.MODEL))
    names = [_mj.mj_id2name(_mm, _mj.mjtObj.mjOBJ_ACTUATOR, i)
             for i in range(_mm.nu)]
    del _mm
    groups = {}
    for j, nm in enumerate(names):
        mi = classify(nm)
        if mi is None:
            continue
        groups.setdefault(mi.groups[0], []).append(j)
    gmean = {g: float(np.mean(ctrl[seg, idx].mean(axis=0)))
             for g, idx in groups.items() if idx}
    act_sum = float(np.mean(np.sum(ctrl[seg], axis=1)))
    # antagonist-pair summary
    pairs = {}
    for fl, ex in ANTAG_PAIRS:
        if fl in gmean and ex in gmean:
            pairs[f"{fl}|{ex}"] = (round(gmean[fl], 3), round(gmean[ex], 3))
    # joint envelopes over the steady window (deg)
    jenv = {KEY_JOINTS[i]: (round(float(q[seg, i].min()), 1),
                            round(float(q[seg, i].max()), 1))
            for i in range(q.shape[1])}
    # buckling forensics: (a) first time COM drops 8 cm below its 1-s
    # value; (b) EARLIEST per-joint divergence after warmup release -
    # each KEY_JOINT's first crossing of +/-3 deg from its t=1.0 s value
    # and first crossing of +/-10 deg (which joint ran away first).
    com_ref = float(com[min(int(np.searchsorted(t, 1.0)), n - 1), 2])
    drop = np.flatnonzero(com[:, 2] < com_ref - 0.08)
    forensics = None
    i1 = min(int(np.searchsorted(t, 1.0)), n - 1)
    div = []
    for i in range(q.shape[1]):
        ref = q[i1, i]
        d3 = np.flatnonzero(np.abs(q[i1:, i] - ref) > 3.0)
        d10 = np.flatnonzero(np.abs(q[i1:, i] - ref) > 10.0)
        if d3.size:
            div.append((float(t[i1 + d3[0]]), KEY_JOINTS[i],
                        round(float(q[i1 + d3[0], i] - ref), 1),
                        float(t[i1 + d10[0]]) if d10.size else None))
    div.sort()
    if drop.size or div:
        k0 = int(drop[0]) if drop.size else i1
        deltas = q[k0] - q[i1]
        order = np.argsort(-np.abs(deltas))
        forensics = dict(
            t_buckle_s=round(float(t[k0]), 2) if drop.size else None,
            com_z_ref=round(com_ref, 3),
            com_z_at=round(float(com[k0, 2]), 3) if drop.size else None,
            first_divergences=[(nm, dt, ddeg,
                                (round(t10, 2) if t10 is not None else None))
                               for dt, nm, ddeg, t10 in div[:6]],
            top_movers_at_buckle=[(KEY_JOINTS[i], round(float(deltas[i]), 1))
                                  for i in order[:5]],
        )
    feet_n = 0.5 * float(np.mean(contact[seg, 0] + contact[seg, 1]))
    return dict(
        tag=tag, t_end=float(m.get("t_end", t[-1])),
        fell=bool(m.get("bal_fell", False)),
        feet_pctBW=round(100.0 * feet_n / BW, 1),
        com_z_mean=round(float(np.mean(com[seg, 2])), 3),
        com_z_min=round(float(np.min(com[seg, 2])), 3),
        sway=round(float(m.get("bal_sway", np.nan)), 3),
        tilt_max=round(float(m.get("bal_tilt_max", np.nan)), 1),
        act_sum=round(act_sum, 1),
        cocon_act=round(float(m.get("bal_cocon_act", np.nan)), 2),
        group_mean={g: round(v, 3) for g, v in gmean.items()},
        antag_pairs=pairs, joints=jenv, forensics=forensics,
    )


def main_matrix():
    p = load_incumbent()
    print(f"[coact:{FAMILY}] incumbent quiet winner rig="
          f"{p.get('rig_scale')} BW={BW:.1f} N; matrix starting", flush=True)
    results = {"family": FAMILY, "bw_N": round(BW, 1),
               "incumbent": {k: p[k] for k in p if k in
                             ("vest_ext", "vest_flex_inh", "vest_prop",
                              "rig_scale", "drive")}}

    # -- phase 0: pure incumbent reference at rig 0 (no alpha at all)
    m, a = run_one(None, None, rig=0.0, tag="ref_incumbent_rig0")
    results["ref_incumbent_rig0"] = a
    _report("REF rig0 no-alpha", a)

    # -- phase 1: existence matrix - saturating 5 nA port current (the
    # ALPHA cell's own V saturates at ~5 mV; beyond this the current
    # cannot add drive), rising uniform MN conductance.
    matrix = []
    for g in (0.5, 2.0, 4.0, 8.0):
        m, a = run_one(g, 5.0, rig=0.0, tag=f"lock_g{g}")
        matrix.append(a)
        _report(f"LOCK g={g} I=5nA", a)
        # stop when every antagonist pair is near-ceiling
        mins = [min(v) for v in a["antag_pairs"].values()]
        if mins and min(mins) >= 0.85:
            print(f"[coact:{FAMILY}] lock reached at g={g} "
                  f"(min antagonist mean {min(mins):.2f})", flush=True)
            break
    results["existence_matrix"] = matrix
    # the working lock = the LAST tested (largest) conductance
    g_lock = matrix[-1]["alpha_to_mn"]
    results["g_lock"] = g_lock

    # -- phase 2: descent in drive at fixed lock conductance
    descent = []
    for frac in (1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.0):
        cur = 5.0 * frac
        m, a = run_one(g_lock, cur, rig=0.0, tag=f"desc_{frac:.1f}")
        descent.append(a)
        _report(f"DESC {int(frac*100)}% (I={cur:.2f} nA)", a)
    results["descent"] = descent

    out = os.path.join(OUTDIR, f"coact_stand_{FAMILY}_20261021.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"[coact:{FAMILY}] wrote {out}", flush=True)


def main_stair():
    p = load_incumbent()
    results = json.loads(open(os.path.join(
        OUTDIR, f"coact_stand_{FAMILY}_20261021.json"),
        encoding="utf-8").read()) \
        if os.path.exists(os.path.join(
            OUTDIR, f"coact_stand_{FAMILY}_20261021.json")) else {}
    g_lock = results.get("g_lock", 4.0)
    P.G["alpha_to_mn"] = float(g_lock)
    levels = [5.0 * f for f in (1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.0)]
    dur = 8.0 * len(levels)      # ~4.6 s dwell per level + lock lead-in
    spec = ",".join(repr(v) for v in levels)
    os.environ["AARL_NPZ"] = f"_coactstand_{FAMILY}_stair_20261021.npz"
    t0 = time.time()
    m = R.main(["--stand-eval", repr(dur), "--rig-scale", "0.0",
                "--alpha-stand", spec])
    wall = time.time() - t0
    # per-dwell-window analysis
    d = np.load(os.environ["AARL_NPZ"])
    t, com, contact, ctrl, q = (d["t"], d["com"], d["contact"], d["ctrl"],
                                d["q"])
    n = min(int(np.sum(t > 0)) + 1, len(t))
    dwell = dur / len(levels)
    wins = []
    for i, lv in enumerate(levels):
        # runner switches level i in at t = i*dwell; skip 0.8 s of
        # transition after each switch
        lo = min(int(np.searchsorted(t[:n], i * dwell + 0.8)), n - 1)
        hi = min(int(np.searchsorted(t[:n], (i + 1) * dwell)), n - 1)
        seg = slice(lo, max(hi, lo + 1))
        feet_n = 0.5 * float(np.mean(contact[seg, 0] + contact[seg, 1]))
        wins.append(dict(
            level_nA=lv,
            feet_pctBW=round(100.0 * feet_n / BW, 1),
            com_z_mean=round(float(np.mean(com[seg, 2])), 3),
            com_z_min=round(float(np.min(com[seg, 2])), 3),
            act_sum=round(float(np.mean(np.sum(ctrl[seg], axis=1))), 1),
            knee_deg=[round(float(np.min(q[seg, 4])), 1),
                      round(float(np.max(q[seg, 4])), 1)],
        ))
    entry = dict(g_lock=g_lock, dur_s=dur, dwell_s=round(dwell, 2),
                 fell=bool(m.get("bal_fell", False)),
                 com_z_global_min=round(float(m.get("bal_com_z_min", 0)), 3),
                 windows=wins, wall_s=round(wall, 1))
    results["staircase"] = entry
    out = os.path.join(OUTDIR, f"coact_stand_{FAMILY}_20261021.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"[coact:{FAMILY}] staircase: {json.dumps(entry)}", flush=True)
    print(f"[coact:{FAMILY}] wrote {out}", flush=True)


def main_one():
    g = float(sys.argv[3])
    cur = float(sys.argv[4])
    rig = float(sys.argv[5])
    dur = float(sys.argv[6]) if len(sys.argv) > 6 else 8.0
    load_incumbent()
    m, a = run_one(g, cur, rig=rig, dur=dur, tag="one")
    a["raw_metrics"] = {k: v for k, v in m.items()
                        if not isinstance(v, dict)}
    print("ONE " + json.dumps(a, default=str))


def _report(hdr, a):
    print(f"[coact:{FAMILY}] {hdr:24s} fell={int(a['fell'])} "
          f"feet={a['feet_pctBW']}%BW comz={a['com_z_mean']} "
          f"sway={a['sway']} actsum={a['act_sum']} "
          f"pairs={json.dumps(a['antag_pairs'])}", flush=True)
    if a["forensics"]:
        print("    diverge: "
              + json.dumps(a["forensics"]["first_divergences"]), flush=True)
        bk = {k: a["forensics"][k] for k in ("t_buckle_s", "com_z_ref",
                                             "com_z_at",
                                             "top_movers_at_buckle")}
        print("    buckle: " + json.dumps(bk), flush=True)


def main_probe():
    """Symmetric-pose MAX-lock probe (mechanism-gap diagnosis): the
    incumbent start pose is a split stance (R hip +24.6 deg forward, L
    -16.6 back); the symmetric pose removes that asymmetry."""
    os.environ["AARL_POSE"] = "symmetric"
    p = load_incumbent()
    results = json.loads(open(os.path.join(
        OUTDIR, f"coact_stand_{FAMILY}_20261021.json"),
        encoding="utf-8").read())
    g = float(sys.argv[3]) if len(sys.argv) > 3 else float(results["g_lock"])
    cur = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    m, a = run_one(g, cur, rig=0.0, tag=f"sym_g{g}_I{cur}")
    _report(f"PROBE sym g={g} I={cur}", a)
    results.setdefault("probes", []).append(a)
    out = os.path.join(OUTDIR, f"coact_stand_{FAMILY}_20261021.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"[coact:{FAMILY}] wrote {out}", flush=True)


def main_coactmax():
    """The OTHER mechanism named in the directive, at its max: the coact
    knob (POSTURE->flexor-antagonist drive, capped by definition at the
    extensor posture weight) at 1.0, no ALPHA, rig 0."""
    p = load_incumbent()
    P.G["coact"] = 1.0
    m, a = run_one(None, None, rig=0.0, tag="coactmax_rig0")
    _report("PROBE coact=1.0 no-alpha", a)
    results = json.loads(open(os.path.join(
        OUTDIR, f"coact_stand_{FAMILY}_20261021.json"),
        encoding="utf-8").read())
    results["probe_coact_max"] = a
    out = os.path.join(OUTDIR, f"coact_stand_{FAMILY}_20261021.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"[coact:{FAMILY}] wrote {out}", flush=True)


if __name__ == "__main__":
    if MODE == "matrix":
        main_matrix()
    elif MODE == "stair":
        main_stair()
    elif MODE == "one":
        main_one()
    elif MODE == "probe":
        main_probe()
    elif MODE == "coactmax":
        main_coactmax()
