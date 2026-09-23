"""One-stop walker status: optuna studies, top trials, log milestones,
orphan processes, newest npz files. Read-only.

Usage:
  status.py                     # everything
  status.py --study NAME        # + top trials of one study (dynamic param cols)
  status.py --spinal PATH       # override spinal dir (or env AARL_SPINAL)
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402

# Column priority for the top-trials table: curriculum stage-3 key order first,
# then anything else the trials actually carry (never stale hardcoded keys).
KEY_PRIORITY = [
    "drive", "pf_gain", "ky_scale", "pm_gain", "pm_T", "pm_ws", "pm_add",
    "pm_aff", "contra_swing", "contra_kinh", "ia_f_contra_f", "c1_gain",
    "v3_gain", "full_rules", "no_cross", "ib_rge", "ia_in",
    "v3_to_ibexc", "ankle_post_walk_trim", "contact_onset", "pelvis_ty",
    "f1_anklepf_inh", "rg_adapt", "desc_e", "desc_f", "rg_to_pf", "rg_nap_h",
]


def studies_section(spinal: Path) -> dict[str, object]:
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    url = wcommon.db_url(spinal)
    sums = optuna.get_all_study_summaries(url)
    print("== optuna studies ==")
    if not sums:
        print("  (none)")
    for s in sums:
        best = "-"
        if s.best_trial is not None and s.directions[0].name == "MAXIMIZE":
            best = f"{s.best_trial.value:.3f} @ t{s.best_trial.number}"
        elif s.best_trial is not None:
            best = f"{s.best_trial.value:.3f} @ t{s.best_trial.number} (MIN)"
        print(f"  {s.study_name:34s} trials={s.n_trials:4d} best={best}")
    return {s.study_name: s for s in sums}


def top_section(spinal: Path, study_name: str, top: int = 8) -> None:
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.load_study(study_name=study_name, storage=wcommon.db_url(spinal))
    done = [t for t in study.trials if t.value is not None and t.state.name == "COMPLETE"]
    if not done:
        print(f"== {study_name}: no complete trials ==")
        return
    done.sort(key=lambda t: t.value, reverse=True)  # curriculum maximizes
    print(f"== {study_name}: {len(done)} complete, best {done[0].value:.3f} (trial {done[0].number}) ==")
    rows = done[:top]
    keys = [k for k in KEY_PRIORITY if any(k in t.params for t in rows)]
    extra = []
    for t in rows:
        for k in t.params:
            if k not in keys and k not in extra:
                extra.append(k)
    keys += extra[: max(0, 10 - len(keys))]
    hdr = "  trial    value    " + " ".join(f"{k[:9]:>9s}" for k in keys)
    print(hdr)
    for t in rows:
        cells = []
        for k in keys:
            v = t.params.get(k)
            cells.append(f"{v:9.3f}" if isinstance(v, float) else f"{str(v)[:9]:>9s}")
        print(f"  t{t.number:<6d} {t.value:8.3f}  " + " ".join(cells))


def logs_section(spinal: Path, n_logs: int = 3) -> None:
    print("== curriculum logs (newest) ==")
    import glob
    logs = []
    for pat in ("curriculum_*.log", "curr_s*.log", "run_*.log"):
        logs += [Path(g) for g in glob.glob(str(spinal / pat))]
    logs = sorted({p.resolve() for p in logs}, key=lambda p: p.stat().st_mtime, reverse=True)
    if not logs:
        print("  (none)")
        return
    now = time.time()
    for lg in logs[:n_logs]:
        age = wcommon.fmt_age(now - lg.stat().st_mtime)
        print(f"  {lg.name}  ({lg.stat().st_size // 1024} KB, {age})")
        for line in wcommon.tail_milestones(lg, 2):
            print(f"    | {line[:150]}")


def procs_section() -> None:
    print("== walker processes ==")
    procs = wcommon.running_walker_processes()
    if not procs:
        print("  none running")
    for pid, name, cmd in procs:
        print(f"  [{pid}] {name}: {cmd}")


def npz_section(spinal: Path, n: int = 5) -> None:
    print("== newest spinal_run*.npz ==")
    files = sorted(spinal.glob("spinal_run*.npz"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not files:
        print("  (none)")
        return
    now = time.time()
    for f in files[:n]:
        st = f.stat()
        print(f"  {f.name:34s} {st.st_size / 1e6:7.1f} MB  {wcommon.fmt_age(now - st.st_mtime)}")


def main() -> int:
    wcommon.utf8_stdio()
    ap = argparse.ArgumentParser()
    ap.add_argument("--study")
    ap.add_argument("--spinal")
    ap.add_argument("--top", type=int, default=8)
    args = ap.parse_args()

    spinal = Path(args.spinal) if args.spinal else wcommon.resolve_spinal()
    wcommon.banner(f"spinal={spinal}")
    try:
        studies_section(spinal)
        if args.study:
            top_section(spinal, args.study, args.top)
    except Exception as e:  # db missing / locked / optuna import failure
        print(f"  (optuna read failed: {type(e).__name__}: {e})")
    logs_section(spinal)
    procs_section()
    npz_section(spinal)
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    raise SystemExit(main())
