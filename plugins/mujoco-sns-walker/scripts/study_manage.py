"""Optuna study lifecycle for the walker curriculum: show / seed / archive / purge.

  study_manage.py show NAME [n]          deep dive: top n trials with params
  study_manage.py seed NAME TRIAL        rewrite curriculum_stage3.json from a trial
  study_manage.py archive NAME...        dump all trials to dated json (NO delete)
  study_manage.py purge NAME... --yes    archive FIRST, then delete study from db

purge is destructive: it refuses without --yes and always archives before
deleting; the archive path is printed. All other actions are read-only
except `seed` (one json rewrite).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402


def load_study(spinal: Path, name: str):
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    return optuna.load_study(study_name=name, storage=wcommon.db_url(spinal))


def show(spinal: Path, name: str, n: int) -> None:
    s = load_study(spinal, name)
    trials = [t for t in s.trials if t.value is not None and t.state.name == "COMPLETE"]
    trials.sort(key=lambda t: t.value, reverse=True)
    print(f"== {name}: {len(trials)} complete ==")
    for t in trials[:n]:
        ps = " ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                      for k, v in t.params.items())
        print(f"  t{t.number:<4d} v={t.value:9.3f}  {ps}")


def archive(spinal: Path, name: str) -> Path | None:
    try:
        s = load_study(spinal, name)
    except Exception as e:
        print(f"{name}: not present ({type(e).__name__}) - nothing to archive")
        return None
    dump = [{"number": t.number, "value": t.value, "params": t.params,
             "state": t.state.name} for t in s.trials]
    out = spinal / f"curriculum_{name}_archive_{time.strftime('%Y%m%d')}.json"
    out.write_text(json.dumps(dump, indent=1), encoding="utf-8")
    best = max((t for t in s.trials if t.value is not None), key=lambda t: t.value, default=None)
    bb = f"best {best.value:.3f} @ t{best.number}" if best else "no values"
    print(f"{name}: archived {len(dump)} trials ({bb}) -> {out.name}")
    return out


def purge(spinal: Path, name: str, confirmed: bool) -> None:
    import optuna
    if not confirmed:
        print(f"{name}: purge needs --yes (it DELETES the study from the db after archiving)")
        return
    path = archive(spinal, name)
    if path is None:
        return
    optuna.delete_study(study_name=name, storage=wcommon.db_url(spinal))
    print(f"{name}: DELETED from db (archive kept: {path.name})")


def seed(spinal: Path, name: str, trial_no: int) -> None:
    s = load_study(spinal, name)
    t = next((x for x in s.trials if x.number == trial_no), None)
    if t is None or t.value is None:
        raise SystemExit(f"trial {trial_no} not found / no value in {name}")
    payload = {"stage": 3, "score": t.value, "trial": trial_no,
               "study": name, "params": t.params}
    out = spinal / "curriculum_stage3.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    pf = t.params.get("pf_gain")
    cs = t.params.get("contra_swing")
    extra = ""
    if pf is not None:
        extra += f" pf {pf:.2f}"
    if cs is not None:
        extra += f" cs {cs:.2f}"
    print(f"curriculum_stage3.json -> {name} t{trial_no} value {t.value:.3f};{extra}")


def main() -> int:
    wcommon.utf8_stdio()
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["show", "seed", "archive", "purge"])
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--trial", type=int, help="seed: trial number (required)")
    ap.add_argument("-n", "--top", type=int, default=10)
    ap.add_argument("--yes", action="store_true", help="purge confirmation")
    args = ap.parse_args()

    spinal = wcommon.resolve_spinal()
    wcommon.banner(f"{args.action} spinal={spinal}")

    if args.action == "seed":
        if args.trial is None:
            ap.error("seed needs --trial N")
        if len(args.studies) != 1:
            ap.error("seed takes exactly one study")
        seed(spinal, args.studies[0], args.trial)
        return 0

    for name in args.studies:
        if args.action == "show":
            show(spinal, name, args.top)
        elif args.action == "archive":
            archive(spinal, name)
        elif args.action == "purge":
            purge(spinal, name, args.yes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
