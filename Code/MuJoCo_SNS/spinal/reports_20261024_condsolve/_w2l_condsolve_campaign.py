r"""W2L CONDUCTANCE-CLASS SOLVE (Batch BW Tier 2, 2026-10-24, egg; Ben:
"shouldn't the Optuna optimizer solve for all conductances?" — yes, at the
SYNAPSE-CLASS level).

Net: --net=bilateralrg2 (Ben's 2026-10-18 corrected aproj; census 96/204)
via w2l_cpg.build_w2l_net. Signs and topology are NOT search dims (they
carry Ben's rulings; see deviation_table.md in this folder).

DRIVE IS PINNED: te=3.0, tf=4.0 nA — Ben's literal corrected-aproj
defaults, fixed in every trial (the free-drive compensation flaw is what
invalidated the 2026-10-20 ladder; dscale stays 1.0).

Search space = the distinct CONDUCTANCE CLASSES:
  existing absolute knobs (contact2 campaign, same bounds for
  comparability): ch_hc (heel->ipsi extensor cascade, weak class),
  ch_ank (heel/toe->ipsi ankle MN), ch_contra (heel->contra RG flx),
  ctn_amp/cref (contact encoder amplitude/saturation), ibnA (Ib group
  current), comm (c1 3.0*comm + V3 0.08*comm commissural classes);
  NEW class scales (test_w2l_ground KV s_*, default 1.0 = the documented
  W2L_GAINS anchor = bit-identical default build, proven 2026-10-24):
  s_rg_lam (rg_laminate 4.0), s_rg_exc (rg_direct_exc 0.5),
  s_pf_drive (2.4), s_pf_cross (4.0), s_pf_mn (2.0), s_rc (Renshaw pair),
  s_ia (Ia chain trio), s_aff_hc (aff_HC_excite 0.08).
NO per-edge search (204 free synapses = lottery-bait; the program's
anti-lottery doctrine). JSON RULE: the winner json carries every knob;
the loader branch is test_w2l_ground's KV/gains plumbing (CLI --k=v sets
KV; the s_* branch injects anchored products).

Metric + survival gate VERBATIM from the campaigns of record
(_w2l_wean_eval imported untouched, exactly like _w2l_contact2_campaign).
Per-term diagnostics as optuna user_attrs every trial + csv rows
(anti-lottery). Winner verified in a FRESH PROCESS before any claim;
acceptance = the strict basin gate (_basin_condsolve.py, 12 fresh
trials, +/-1% on every knob, bar 11/12) + inter-leg phase metrics.

Study: curr_w2l_condsolve in this folder's optuna_w2l_condsolve.db.
NO wean-ladder claims in this batch (Ben's rule).

Usage:
    D:\\Anaconda\\envs\\myo\\python.exe _w2l_condsolve_campaign.py [--n=80]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WEAN = os.path.normpath(os.path.join(
    HERE, "..", "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import _w2l_wean_eval as EV       # noqa: E402  (metric + evaluate, verbatim)
import optuna                     # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)

EV.T.KV["net"] = "bilateralrg2"   # the corrected net is a CONSTANT here

DB = "sqlite:///" + os.path.join(HERE, "optuna_w2l_condsolve.db").replace("\\", "/")
SEED = 20261024
DUR_TUNE, DUR_VERIFY = 20.0, 30.0
PY = sys.executable
STUDY = "curr_w2l_condsolve"

# PINNED DRIVE (not in SPACE): Ben's literal corrected-aproj defaults.
TE_PIN, TF_PIN = 3.0, 4.0

SPACE = dict(
    # contact2's absolute class knobs (same bounds = comparable record)
    ch_hc=(0.01, 1.0),      # heel -> ipsi extensor cascade (weak class)
    ch_ank=(0.1, 2.0),      # heel/toe -> ipsi ankle MN ext
    ch_contra=(0.01, 2.0),  # heel -> contra RG flx (the antiphase edge)
    ctn_amp=(0.5, 6.0),     # contact-encoder amplitude (nA)
    cref=(20.0, 150.0),     # contact-encoder saturation force (N)
    ibnA=(0.0, 3.0),        # Ib group sensor current (nA)
    comm=(0.2, 2.0),        # commissural scale (c1 3.0*comm, V3 0.08*comm)
    # NEW: synapse-class scales over W2L_GAINS anchors (1.0 = anchor)
    s_rg_lam=(0.25, 4.0),   # rg_laminate 4.0
    s_rg_exc=(0.25, 4.0),   # rg_direct_exc 0.5
    s_pf_drive=(0.25, 4.0),  # pf_drive 2.4
    s_pf_cross=(0.25, 4.0),  # pf_cross 4.0
    s_pf_mn=(0.25, 4.0),    # pf_to_mn 2.0
    s_rc=(0.25, 4.0),       # rc_own_exc 1.0 + rc_inh 0.5
    s_ia=(0.25, 4.0),       # ia_sn_exc 0.6 + ia_inh 0.4 + ia_mutual 0.5
    s_aff_hc=(0.25, 4.0),   # aff_HC_excite 0.08
)
SCALES = tuple(k for k in SPACE if k.startswith("s_"))
_ONE = {k: 1.0 for k in SCALES}

# Seed 1 = Ben's literal defaults (documented anchors, scales at 1.0).
# Seeds 2/3 = the contact2 r100/r030 winners translated to the PINNED drive
# (their class/encoder gains, scales 1.0): does the solved contact encoding
# transfer at Ben's own drive?
SEEDS = [
    dict(ch_hc=0.08, ch_ank=1.0, ch_contra=0.08,
         ctn_amp=4.0, cref=50.0, ibnA=1.0, comm=1.0, **_ONE),
    dict(ch_hc=0.1423, ch_ank=0.7184, ch_contra=0.113,
         ctn_amp=4.4267, cref=145.7223, ibnA=0.1808, comm=1.2462, **_ONE),
    dict(ch_hc=0.4783, ch_ank=0.9118, ch_contra=1.1022,
         ctn_amp=4.6096, cref=129.2832, ibnA=0.1667, comm=0.7005, **_ONE),
]


def knobs_of(p: dict, dur: float) -> dict:
    return dict(EV.BASE,
                dscale=1.0,
                te=TE_PIN, tf=TF_PIN,       # PINNED - never a search dim
                **{k: float(p[k]) for k in SPACE},
                dur_walk=dur)


def fresh_verify(params: dict) -> dict:
    kn = knobs_of(params, DUR_VERIFY)
    tmp = os.path.join(HERE, "_finalist_knobs.json")
    json.dump(kn, open(tmp, "w"), indent=1)
    r = subprocess.run(
        [PY, os.path.join(HERE, "_verify_finalist.py"), tmp],
        capture_output=True, text=True, timeout=1800)
    for line in r.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    raise RuntimeError(f"fresh verify failed:\n{r.stdout[-800:]}\n{r.stderr[-800:]}")


def main() -> int:
    n = 80
    for a in sys.argv[1:]:
        if a.startswith("--n="):
            n = int(a[4:])
    study = optuna.create_study(
        study_name=STUDY, storage=DB, direction="maximize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=SEED))
    have = len(study.trials)
    if have == 0:
        for sp in SEEDS:
            study.enqueue_trial(sp)
    need = max(n - have, 0)
    if need == 0 and have == 0:
        need = 1

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi)
             for k, (lo, hi) in SPACE.items()}
        m = EV.evaluate(knobs_of(p, DUR_TUNE),
                        tag=f"{STUDY}#{trial.number}",
                        csv_path=os.path.join(HERE, "condsolve_trials.csv"))
        for k, v in m.items():
            if not k.startswith("_") and v is not None:
                try:
                    trial.set_user_attr(k, v)
                except (TypeError, ValueError):
                    pass
        return float(m["obj"])

    if need > 0:
        study.optimize(obj, n_trials=need)
    best = study.best_trial
    out = dict(study=STUDY, n_trials=len(study.trials), pinned_drive=dict(te=TE_PIN, tf=TF_PIN),
               params={k: round(float(best.params[k]), 5) for k in SPACE},
               obj20=float(best.value))
    m30 = fresh_verify(best.params)
    out["m30_fresh"] = {k: m30.get(k) for k in EV.METRIC_KEYS}
    out["survive30"] = bool(m30.get("survive"))
    path = os.path.join(HERE, "condsolve_best.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print(f"== CONDSOLVE best obj20={out['obj20']:.1f} "
          f"survive30(fresh)={out['survive30']} -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
