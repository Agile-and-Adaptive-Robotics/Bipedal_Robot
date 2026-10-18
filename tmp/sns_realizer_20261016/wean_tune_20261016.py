"""Human-model feedback-tuning + DRIVE-WEAN campaign driver (2026-10-16,
overnight walk campaign phase: Ben's wean ladder).

Base = campaign anchor A (bit-exact -161.56754173676563): v10 multipliers
+ renshaw 0.5 + s3k trial-34 winner params + full_rules=0, set_stage(4,..)
+ the four L/R hand-drive knobs (gamma_l/r, alpha_l/r, alpha_to_mn) set
AFTER set_stage, exactly the DRIVE_VERIFY section-6 recipe the drive sweep
used (drive_sweep_one_20261016.py lines 26-39).

WEAN SEMANTICS (Ben's strategy, campaign handoff "ladder ~1.0 -> <=0.7 ->
<=0.4"): the weaned quantity is the HAND-DRIVE crutch only — the knob
currents u[GAMMA_c_s]=drive*G[gamma_s], u[ALPHA_c_s]=drive*G[alpha_s] are
scaled by drive_scale (the s3k descending DRIVE itself stays). Effective
knob = knob * drive_scale.

ROUND IDS (2026-10-16 evening restart):
  1  curr_s6b_drive_human  — eval window, 20 keys. SUPERSEDED: a bug made
     the 13 FB keys INERT (set_stage re-loaded s3k values; only the 5
     drive knobs varied). Its winner (trial 87, -146.03) is a VALID
     walker but only drive-tuned.
  2  curr_s7_wean07 — run under the same bug: a scale-only sweep
     (feedback NOT re-tuned). Invalid as the wean re-tune.
  4  curr_s6c_prodwindow — same bug, aborted at 4 trials.
  5  curr_s6d_feedback — REDO of round 1 with FB injection FIXED
     (fb merged into base BEFORE set_stage). Eval window, 20 live keys.
  6  curr_s7b_wean07 — REDO of round 2 (drive_scale<=0.7, FB re-tuned,
     seeded from round 5).
  7  curr_s6e_prodwindow — REDO of round 4 (production 22 s schedule
     walk window scored, seeded from round 5 winner + bootstrap).

Anti-lottery: every trial's FULL per-term kine diagnostics (rmse/mean/
range per joint per leg, T, period_cv, lag_rl, ds, duty, knee_min,
cycles, frozen flags) are logged to round{r}_w{worker}.jsonl.

Usage (cwd = Code\\MuJoCo_SNS\\spinal, python = D:\\Anaconda\\envs\\myo):
  wean_tune_20261016.py one <spec.json>          # single eval, prints RESULT
  wean_tune_20261016.py gate                     # anchor + bootstrap gates
  wean_tune_20261016.py round <r> <n> <worker>   # optuna worker
  wean_tune_20261016.py best <r> [rank]          # extract winner json
spec.json (one): {"tag":..., "knobs":{...}, "scale":1.0, "fb":{...},
                  "overrides":{G-key:value,...}, "winner":<best.json path>,
                  "npz": abs path, "jsonl": abs path}
"""
import json
import os
import sys
import time

# NOTE: no sys.stdout TextIOWrapper here — importing _curriculum would wrap
# it a second time and the first wrapper's GC closes the real stdout
# ("I/O operation on closed file", caught in the first gate run). The
# launchers set PYTHONIOENCODING=utf-8 instead (drive-sweep pattern).

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "wean_human")
os.makedirs(OUT, exist_ok=True)

DB = "sqlite:///optuna_walk.db"
BASE_DRIVE = 2.0773198435110425          # s3k winner drive (anchor A)
ANCHOR_BASELINE = -161.56754173676563    # knobs-0 gate value
BOOTSTRAP_EVAL = -173.0970129592704      # sweep row g1.5_a1_d0
DRIVE_KEYS = ("gamma_l", "gamma_r", "alpha_l", "alpha_r", "alpha_to_mn")
BOOT_KNOBS = dict(gamma_l=1.5, gamma_r=1.5, alpha_l=1.0, alpha_r=1.0,
                  alpha_to_mn=0.5)
STUDIES = {1: "curr_s6b_drive_human", 2: "curr_s7_wean07",
           3: "curr_s8_wean04", 4: "curr_s6c_prodwindow",
           5: "curr_s6d_feedback", 6: "curr_s7b_wean07",
           7: "curr_s6e_prodwindow", 8: "curr_s8b_wean04"}
SCALE_RANGE = {1: None, 2: (0.30, 0.70), 3: (0.05, 0.40), 4: None,
               5: None, 6: (0.30, 0.70), 7: None, 8: (0.05, 0.40)}
MODE = {1: "eval", 2: "wean", 3: "wean", 4: "prod",
        5: "eval", 6: "wean", 7: "prod", 8: "wean"}
PREV = {2: 1, 3: 2, 6: 5, 7: 5, 8: 6}    # seeding round for wean/prod

# feedback keys searched every round (all have set_stage(4)/set_params
# loaders; the four knobs load via the post-set_stage recipe)
FB_KEYS = ("heel_rge", "toe_rge", "ib_rge", "ia_in", "ib_e_central",
           "ii_e_central", "ii_f_central", "ia_f_central",
           "phase_reset_e", "phase_reset_f", "contact_onset",
           "f1_kneext_inh", "f1_anklepf_inh")

# round-1/5 ranges: every one contains the anchor/s3k incumbent value
R1 = {
    "gamma_l": (0.75, 2.5), "gamma_r": (0.75, 2.5),
    "alpha_l": (0.5, 1.5), "alpha_r": (0.5, 1.5),
    "alpha_to_mn": (0.25, 1.0),
    "heel_rge": (0.0, 1.0), "toe_rge": (0.0, 0.75),
    "ib_rge": (0.0, 1.25), "ia_in": (0.0, 1.25),
    "ib_e_central": (0.0, 0.75), "ii_e_central": (0.0, 0.75),
    "ii_f_central": (0.0, 0.75), "ia_f_central": (0.0, 0.75),
    "phase_reset_e": (0.0, 1.25), "phase_reset_f": (0.0, 2.0),
    "contact_onset": (0.0, 1.25),
    "f1_kneext_inh": (0.25, 2.0), "f1_anklepf_inh": (0.0, 1.5),
}

MODS = {}   # C, P, R module refs, filled by apply_config


def load_base():
    base = dict(json.loads(open("best_walk_params_v10.json",
                                encoding="utf-8").read())["multipliers"])
    base["renshaw"] = 0.5
    s3k = json.loads(open("reports_20260923/s3k_trial34_full_params.json",
                          encoding="utf-8").read())["params"]
    base.update(s3k)
    base["full_rules"] = 0.0
    return base


def apply_config(knobs, scale=1.0, fb=None, overrides=None):
    """set_stage(4, base+fb) then the four knobs (x scale) + alpha_to_mn.

    fb MUST carry the searched feedback keys (merged into base BEFORE
    set_stage — the 17:10 bug: set_stage re-loaded s3k values and the
    sampled FB keys were never applied, making rounds 1/2/4 FB-inert)."""
    import _curriculum as C
    import params as P
    import runner as R
    MODS.update(C=C, P=P, R=R)
    base = load_base()
    if fb:
        base.update({k: float(v) for k, v in fb.items()})
    C.set_stage(4, base)
    eff = {}
    for k in ("gamma_l", "gamma_r", "alpha_l", "alpha_r"):
        v = float(knobs[k]) * float(scale)
        P.G[k] = v
        eff[k] = v
    P.G["alpha_to_mn"] = float(knobs.get("alpha_to_mn",
                                         BOOT_KNOBS["alpha_to_mn"]))
    eff["alpha_to_mn"] = P.G["alpha_to_mn"]
    for k, v in (overrides or {}).items():
        P.G[k] = float(v)
    return eff


def flat_metrics(m):
    """runner eval dict + full per-term kine diagnostics, scalars only."""
    out = {}
    for k, v in m.items():
        if isinstance(v, (int, float, bool)) or v is None:
            out[k] = v
    for kk, v in (m.get("kine") or {}).items():
        if isinstance(v, (int, float, bool)) or v is None:
            out["k_" + kk] = v
    return out


def run_eval(tag, knobs, scale=1.0, fb=None, overrides=None, npz=None):
    apply_config(knobs, scale=scale, fb=fb, overrides=overrides)
    if npz:
        os.environ["AARL_NPZ"] = npz
    else:
        os.environ.pop("AARL_NPZ", None)
    m = MODS["R"].main(["--eval", "--drive", repr(BASE_DRIVE)])
    res = dict(tag=tag, scale=scale, knobs=dict(knobs), fb=fb or {},
               overrides=overrides or {}, ts=time.strftime("%H:%M:%S"),
               **flat_metrics(m))
    return res


# ---------------------------------------------------------------- one
def cmd_one():
    spec = json.loads(open(sys.argv[2], encoding="utf-8").read())
    if "winner" in spec:
        w = json.loads(open(spec["winner"], encoding="utf-8").read())
        knobs = {k: float(w["knobs"][k]) for k in DRIVE_KEYS}
        fb = {k: float(w["params"].get(k, 0.0)) for k in FB_KEYS}
        scale = float(spec.get("scale", w.get("scale", 1.0)))
    else:
        knobs = dict(BOOT_KNOBS)
        knobs.update(spec.get("knobs", {}))
        fb = spec.get("fb")
        scale = float(spec.get("scale", 1.0))
    t0 = time.time()
    res = run_eval(spec.get("tag", "one"), knobs, scale=scale, fb=fb,
                   overrides=spec.get("overrides"), npz=spec.get("npz"))
    res["wall_s"] = round(time.time() - t0, 1)
    line = json.dumps(res)
    if spec.get("jsonl"):
        with open(spec["jsonl"], "a", encoding="utf-8") as f:
            f.write(line + "\n")
    print("RESULT " + line, flush=True)


def cmd_spec_eval(tag, knobs, scale):
    spec_path = os.path.join(OUT, "_spec_tmp.json")
    json.dump(dict(tag=tag, knobs=knobs, scale=scale),
              open(spec_path, "w", encoding="utf-8"))
    import subprocess
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.abspath(__file__),
                        "one", spec_path], capture_output=True, text=True,
                       cwd=HERE, env={**os.environ,
                                      "PYTHONIOENCODING": "utf-8"})
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("RESULT ")]
    if not line:
        print(p.stdout[-3000:] if p.stdout else "")
        print(p.stderr[-3000:] if p.stderr else "")
        raise RuntimeError("subprocess eval failed")
    res = json.loads(line[-1][7:])
    res["wall_s"] = round(time.time() - t0, 1)
    return res


# ---------------------------------------------------------------- gate
def cmd_gate():
    r0 = cmd_spec_eval("gate_anchor", {k: 0.0 for k in DRIVE_KEYS}, 1.0)
    ok0 = (r0["kine_score"] == ANCHOR_BASELINE)
    r1 = cmd_spec_eval("gate_bootstrap", dict(BOOT_KNOBS), 1.0)
    ok1 = abs(r1["kine_score"] - BOOTSTRAP_EVAL) < 1e-9
    print(f"GATE anchor {r0['kine_score']!r} exact={ok0} "
          f"bootstrap {r1['kine_score']!r} match={ok1} "
          f"wall {r0['wall_s']}/{r1['wall_s']} s", flush=True)
    if not (ok0 and ok1):
        sys.exit(3)


# ---------------------------------------------------------------- round
def fb_range(key, v):
    """Adaptive range around incumbent v: contains v."""
    if v <= 1e-6:
        return (0.0, 0.3)
    return (max(0.0, 0.25 * v), 1.75 * v)


def cmd_round():
    r = int(sys.argv[2])
    n = int(sys.argv[3])
    worker = sys.argv[4]
    import optuna
    mode = MODE[r]

    if r in (1, 5):
        base0 = load_base()
        inc_fb = {k: float(base0.get(k, 0.0)) for k in FB_KEYS}
        inc_knobs = dict(BOOT_KNOBS)
    else:
        prev = json.loads(open(os.path.join(
            OUT, f"round{PREV[r]}_best.json"), encoding="utf-8").read())
        inc_fb = {k: float(prev["params"].get(k, 0.0)) for k in FB_KEYS}
        inc_knobs = {k: float(prev["knobs"][k]) for k in DRIVE_KEYS}
    seed_scale = 1.0 if SCALE_RANGE[r] is None else SCALE_RANGE[r][1]

    npz = os.path.join(OUT, f"trial_w{worker}.npz")
    jsonl = os.path.join(OUT, f"round{r}_w{worker}.jsonl")

    if mode == "prod":
        # production-window round: the 22 s default schedule leaves the
        # tracked spinal_run.png churned by every trial (not
        # AARL_NPZ-redirectable) — no-op the figure save so concurrent
        # workers cannot race on it (WinError 32) or churn the repo file.
        import matplotlib.figure as _mf
        _mf.Figure.savefig = lambda self, *a, **k: None

    def prod_run_metrics():
        """Default-schedule run + production-window kine (walk 5-15 s)."""
        import numpy as np
        import kine_ref
        os.environ["AARL_NPZ"] = npz
        MODS["R"].main(["--drive", repr(BASE_DRIVE)])
        d = np.load(npz, allow_pickle=True)
        t, q = d["t"], d["q"]
        neuro, com, ct = d["neuro"], d["com"], d["contact"]
        n_done = int(np.sum(t > 0)) + 1
        finite = bool(np.all(np.isfinite(q[:n_done]))
                      and np.all(np.isfinite(com[:n_done])))
        zmin = float(np.nanmin(com[:n_done, 2])) if finite else 0.0
        i_ws = int(np.searchsorted(t[:n_done], 5.0))
        try:
            kk = kine_ref.compare(t[:n_done], q[:n_done], neuro[:n_done],
                                  5.0, ref=kine_ref.ref_cached(),
                                  contact=ct[:n_done])
        except Exception:
            kk = None
        return dict(
            nan=not finite, t_end=float(t[n_done - 1]),
            dx=float(com[n_done - 1, 0] - com[i_ws, 0]),
            kz=zmin,
            tilt_max=float(np.nanmax(np.abs(q[:n_done, 0]))) if finite
            else 0.0,
            knee_min=float(np.nanmin(q[:n_done, 4])) if finite else 0.0,
            kine=kk, kine_score=(kk["kine_score"] if kk else -25.0),
        )

    def objective(trial):
        sug = {}
        if r in (1, 5):               # full 20-key space around anchor
            for k in DRIVE_KEYS:
                sug[k] = trial.suggest_float(k, *R1[k])
        elif mode == "prod":          # drive knobs around round-5 winner
            for k in DRIVE_KEYS:
                sug[k] = trial.suggest_float(k, *fb_range(k, inc_knobs[k]))
        scale = 1.0
        if SCALE_RANGE[r] is not None:
            scale = trial.suggest_float("drive_scale", *SCALE_RANGE[r])
        rng = R1 if r in (1, 5) else {k: fb_range(k, inc_fb[k])
                                      for k in FB_KEYS}
        for k in FB_KEYS:
            sug[k] = trial.suggest_float(k, *rng[k])
        knobs = ({k: sug[k] for k in DRIVE_KEYS}
                 if r in (1, 5) or mode == "prod" else inc_knobs)
        fb = {k: sug[k] for k in FB_KEYS}
        apply_config(knobs, scale=scale, fb=fb)
        if mode == "prod":
            m = prod_run_metrics()
        else:
            os.environ["AARL_NPZ"] = npz
            m = MODS["R"].main(["--eval", "--drive", repr(BASE_DRIVE)])
        if m["nan"]:
            score = -400.0
        elif m.get("kine") is None:
            score = -320.0
        else:
            score = max(float(m["kine_score"]), -315.0)
            if m["kz"] < 0.62:
                score -= 20.0
            if m["tilt_max"] > 40.0:
                score -= 10.0
            # BILATERAL GATE (2026-10-16, the aborted curr_s6 lesson +
            # campaign contract "frozen runs never count"): a trial whose
            # left OR right leg yields <3 contact cycles sits 60 below
            # every bilateral trial. Ordering among bilateral trials is
            # unchanged (raw kine_score + the established adjustments).
            kk = m["kine"]
            if not (int(kk.get("n_cycles_r", 0) or 0) >= 3 and
                    int(kk.get("n_cycles_l", 0) or 0) >= 3):
                score -= 60.0
        P = MODS["P"]
        row = dict(round=r, worker=worker, number=trial.number,
                   value=score, ts=time.strftime("%H:%M:%S"),
                   params=dict(trial.params),
                   eff_knobs={k: float(P.G[k]) for k in
                              ("gamma_l", "gamma_r", "alpha_l", "alpha_r")},
                   metrics=flat_metrics(m))
        with open(jsonl, "a", encoding="utf-8") as f:
            f.write(json.dumps(row) + "\n")
        print(f"[r{r} w{worker} t{trial.number}] {score:.3f}", flush=True)
        return score

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", storage=DB,
                                study_name=STUDIES[r],
                                load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(
                                    # per-worker seed: identical seeds made
                                    # the 4 workers sample identical random
                                    # startup trials (seen in round 1's
                                    # first 13 rows - duplicated evals)
                                    seed=61 + 10 * r + int(worker),
                                    n_startup_trials=8))
    if len(study.trials) == 0 and worker == "0":
        if mode == "prod":
            # seed 1 = round-5 winner (anchors its own production score);
            # seed 2 = the bootstrap combo (production-verified bilateral
            # in the drive sweep, kine -169.21).
            seed = {**inc_fb, **{k: inc_knobs[k] for k in DRIVE_KEYS}}
            base0 = load_base()
            seed2 = {**{k: float(base0.get(k, 0.0)) for k in FB_KEYS},
                     **BOOT_KNOBS}
        else:
            seed = dict(inc_fb)
            if r in (1, 5):
                seed.update({k: BOOT_KNOBS[k] for k in DRIVE_KEYS})
            else:
                seed["drive_scale"] = seed_scale
        study.enqueue_trial(seed)
        print("SEEDED " + json.dumps(seed), flush=True)
        if mode == "prod":
            study.enqueue_trial(seed2)
            print("SEEDED2 " + json.dumps(seed2), flush=True)
    study.optimize(objective, n_trials=n, gc_after_trial=True)
    print(f"WORKER_DONE r{r} w{worker}", flush=True)


# ---------------------------------------------------------------- best
def cmd_best():
    r = int(sys.argv[2])
    rank = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    import glob
    # select from the jsonl rows (they carry per-term diagnostics): the
    # winner must be BILATERAL (>=3 contact cycles each leg, no NaN);
    # among those, max objective value. rank>1 = n-th best (for full-run
    # production-schedule testing of several candidates).
    rows = []
    for p in glob.glob(os.path.join(OUT, f"round{r}_w*.jsonl")):
        rows.extend(json.loads(l) for l in open(p, encoding="utf-8"))
    elig = [x for x in rows
            if int((x["metrics"].get("k_n_cycles_r") or 0)) >= 3
            and int((x["metrics"].get("k_n_cycles_l") or 0)) >= 3
            and not x["metrics"].get("nan", False)]
    if len(elig) < rank:
        raise SystemExit(f"NOT ENOUGH BILATERAL TRIALS in round {r} "
                         f"({len(elig)} < {rank}, {len(rows)} rows total)")
    bt = sorted(elig, key=lambda x: -x["value"])[rank - 1]
    params = load_base()
    params.update({k: float(v) for k, v in bt["params"].items()})
    scale = float(bt["params"].get("drive_scale", 1.0))
    if r in (1, 5) or MODE[r] == "prod":
        knobs = {k: float(bt["params"][k]) for k in DRIVE_KEYS}
    else:
        prev = json.loads(open(os.path.join(
            OUT, f"round{PREV[r]}_best.json"),
            encoding="utf-8").read())
        knobs = {k: float(prev["knobs"][k]) for k in DRIVE_KEYS}
    out = dict(round=r, study=STUDIES[r], trial=bt["number"],
               score=float(bt["value"]), raw_kine_score=bt["metrics"].get(
                   "kine_score"),
               scale=scale, knobs=knobs, params=params,
               eff_knob_currents={k: knobs[k] * scale for k in
                                  ("gamma_l", "gamma_r", "alpha_l",
                                   "alpha_r")},
               note=("wean round; loader recipe: base=dict(this params); "
                     "set_stage(4, base); for k in (gamma_l,gamma_r,"
                     "alpha_l,alpha_r): P.G[k]=knob*scale; "
                     "P.G[alpha_to_mn]=knobs[alpha_to_mn]; or simply "
                     "wean_tune_20261016.py one with spec "
                     "{\"winner\": this file}"))
    path = os.path.join(OUT, f"round{r}_best{rank if rank > 1 else ''}.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print("WROTE " + path + " score " + repr(bt["value"])
          + " trial " + str(bt["number"]) + " rows " + str(len(rows)),
          flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "one":
        cmd_one()
    elif cmd == "gate":
        cmd_gate()
    elif cmd == "round":
        cmd_round()
    elif cmd == "best":
        cmd_best()
    else:
        raise SystemExit("unknown cmd " + cmd)
