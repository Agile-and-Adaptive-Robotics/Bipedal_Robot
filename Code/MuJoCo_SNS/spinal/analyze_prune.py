"""Prune-matrix analysis (Ben 2026-09-28): prune_results_<variant>.jsonl
-> prune_analysis.md + stdout tables.

Verdict rule per cut (component C prunable iff removing it costs
nothing anywhere):
  AIR/AIRDEAFF: air_score delta >= -5 (rhythm quality preserved)
  WALK:         kine_score delta >= -3
  STAND:        no new falls AND bal_sway delta <= +0.01 m
  PUSH:         no new falls on any axis AND mean bal_push_sway
                delta <= +0.01 m
  PRUNE if all pass; KEEP otherwise (the failing modes name themselves).
Usage: python analyze_prune.py
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace") if hasattr(sys.stdout, "buffer") else sys.stdout

VARIANTS = ["s3k", "w2lvar", "syn6"]
PUSHES = ["+x", "-x", "+y", "-y"]
CUTS = ["noaff", "interleg", "contact", "ia", "ii", "ib", "renshaw",
        "rgweak"]

AIR_TOL = -5.0
WALK_TOL = -3.0
SWAY_TOL = 0.01


def _load(v):
    cells = {}
    fn = f"prune_results_{v}.jsonl"
    if not os.path.exists(fn):
        return cells
    with open(fn, encoding="utf-8") as f:
        for line in f:
            try:
                r = json.loads(line)
            except Exception:
                continue
            if "error" in r:
                cells[(r["config"], r["mode"], r.get("push", ""))] = r
                continue
            cells[(r["config"], r["mode"], r.get("push", ""))] = r
    return cells


def _stand_score(m):
    if not m or m.get("bal_fell") or m.get("nan"):
        return -150.0
    return (100.0 - 400.0 * m.get("bal_sway", 1.0)
            - m.get("bal_tilt_max", 90.0)
            - 40.0 * abs(0.5 - m.get("bal_contact_sym", 0.25)))


def _walk_score(m):
    if not m or m.get("nan"):
        return -400.0
    if not (m.get("kine_components") or {}):
        # no cycles detected: the runner's -25 sentinel must sit BELOW
        # every genuine walker (curriculum -320 discipline)
        return -320.0
    ks = m.get("kine_score", -25.0)
    s = max(float(ks), -315.0)
    if m.get("kz", 1.0) < 0.62:
        s -= 20.0
    if m.get("tilt_max", 0.0) > 40.0:
        s -= 10.0
    return s


def main():
    lines = ["# Prune-matrix analysis (easteregg2, 2026-09-28)", ""]
    for v in VARIANTS:
        cells = _load(v)
        if not cells:
            lines += [f"## {v}: NO RESULTS", ""]
            continue
        n_err = sum(1 for r in cells.values() if "error" in r)
        lines += [f"## {v}  ({len(cells)} cells, {n_err} errors)", ""]
        ref = {}
        for mode in ["AIR", "AIRDEAFF", "WALK", "STAND"]:
            ref[mode] = cells.get(("full", mode, ""))
        ref_push = {ax: cells.get(("full", "PUSH", ax)) for ax in PUSHES}
        # header table
        lines += ["| cut | AIR d | AIRDEAFF d | WALK d | STAND d | "
                  "sway d (m) | push fell | push sway d (m) | verdict |",
                  "|---|---|---|---|---|---|---|---|---|"]
        def _pair(mode, cut_):
            """(cut_metrics, full_metrics) or (None, None) if either
            cell is absent - absent reads as n/a, NOT as a failure."""
            c = cells.get((cut_, mode, ""))
            f = ref.get(mode)
            if c is None or f is None or "error" in c:
                return None, None
            return c, f

        for cut in CUTS + (["combo"] if v == "s3k" else []):
            if (cut, "AIR", "") not in cells and \
               (cut, "STAND", "") not in cells:
                continue
            row, fail = [], []

            def _num(val, tol, name, fmt="{:+.1f}", worse="low"):
                if val is None:
                    row.append("n/a")
                    return
                row.append(fmt.format(val))
                if (val < tol) if worse == "low" else (val > tol):
                    fail.append(name)

            c_air, f_air = _pair("AIR", cut)
            da = None if c_air is None else \
                (c_air.get("air", {}).get("air_score") or -999) - \
                (f_air.get("air", {}).get("air_score") or -999)
            _num(da, AIR_TOL, "air")

            c_ad, f_ad = _pair("AIRDEAFF", cut)
            if c_ad is None and cut == "combo":
                row.append("n/a")   # combo deliberately skips AIRDEAFF
            else:
                dd = None if c_ad is None else \
                    (c_ad.get("air", {}).get("air_score") or -999) - \
                    (f_ad.get("air", {}).get("air_score") or -999)
                _num(dd, AIR_TOL, "air-deaff")

            c_w, f_w = _pair("WALK", cut)
            dw = None if c_w is None else \
                _walk_score(c_w.get("metrics", {})) - \
                _walk_score(f_w.get("metrics", {}))
            _num(dw, WALK_TOL, "walk")

            c_s, f_s = _pair("STAND", cut)
            if c_s is None:
                ds = dsway = newfall = None
            else:
                st, stf = c_s.get("metrics", {}), f_s.get("metrics", {})
                ds = _stand_score(st) - _stand_score(stf)
                dsway = (st.get("bal_sway", 9.0) if st else 9.0) - \
                        (stf.get("bal_sway", 9.0) if stf else 9.0)
                newfall = bool(st.get("bal_fell")) and \
                    not bool(stf.get("bal_fell"))
            _num(ds, -10.0, "stand")
            _num(dsway, SWAY_TOL, "stand", fmt="{:+.3f}", worse="high")
            if newfall:
                fail.append("stand")

            c_pushes = [cells.get((cut, "PUSH", ax)) for ax in PUSHES]
            f_pushes = [ref_push.get(ax) for ax in PUSHES]
            if all(c is None for c in c_pushes):
                row.extend(["n/a", "n/a"])
            else:
                pf = psw = 0.0
                newpf = False
                for pm, pr in zip(c_pushes, f_pushes):
                    if pm is None:
                        continue
                    prm = pr.get("metrics", {}) if pr else {}
                    pmm = pm.get("metrics", {})
                    if pmm.get("bal_fell") and not prm.get("bal_fell"):
                        newpf = True
                    pf += 1 if pmm.get("bal_fell") else 0
                    psw += (pmm.get("bal_push_sway", 9.0) or 9.0) - \
                           (prm.get("bal_push_sway", 9.0) or 9.0)
                psw /= max(1, sum(1 for c in c_pushes if c is not None))
                row.append(f"{int(pf)}/4")
                _num(psw, SWAY_TOL, "push", fmt="{:+.3f}", worse="high")
                if newpf:
                    fail.append("push")

            verdict = "PRUNE" if not fail else ("KEEP (" +
                                                ", ".join(fail) + ")")
            lines.append("| " + cut + " | " + " | ".join(row) +
                         " | **" + verdict + "** |")
        lines.append("")
        # FULL reference summary
        f_air = (ref["AIR"] or {}).get("air", {})
        f_walk = (ref["WALK"] or {}).get("metrics", {})
        f_stand = (ref["STAND"] or {}).get("metrics", {})
        lines += [f"FULL reference: air {f_air.get('air_score')} "
                  f"(rises {f_air.get('rises')}, "
                  f"period {f_air.get('period') and round(f_air['period'], 2)}s), "
                  f"walk kine {f_walk.get('kine_score')} "
                  f"(duty {f_walk.get('duty') and round(f_walk['duty'], 2)}, "
                  f"kz {f_walk.get('kz') and round(f_walk['kz'], 2)}), "
                  f"stand score {round(_stand_score(f_stand), 1)} "
                  f"(sway {f_stand.get('bal_sway') and round(f_stand['bal_sway'], 3)}, "
                  f"fell {f_stand.get('bal_fell')})",
                  ""]
    out = "\n".join(lines)
    with open("prune_analysis.md", "w", encoding="utf-8") as f:
        f.write(out + "\n")
    print(out)


if __name__ == "__main__":
    main()
