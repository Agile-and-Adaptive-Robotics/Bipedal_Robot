import json, glob, os

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"

for path in sorted(glob.glob(SP + r"\prune_results_*.jsonl")):
    rows = [json.loads(l) for l in open(path) if l.strip()]
    print("=" * 60)
    print(os.path.basename(path), "rows:", len(rows))
    from collections import defaultdict
    grp = defaultdict(list)
    for r in rows:
        grp[(r.get("config"), r.get("mode"))].append(r)
    for (cfg, mode), rs in sorted(grp.items(), key=lambda kv: str(kv[0])):
        ks = [r["metrics"].get("kine_score") for r in rs if r["metrics"].get("kine_score") is not None]
        fell = sum(1 for r in rs if r["metrics"].get("bal_fell"))
        bilat = sum(1 for r in rs if r["metrics"].get("kine_components", {}).get("bilateral"))
        duty = [r["metrics"].get("kine_components", {}).get("duty") for r in rs
                if r["metrics"].get("kine_components", {}).get("duty") is not None]
        kmin = [r["metrics"].get("kine_components", {}).get("knee_min") for r in rs
                if r["metrics"].get("kine_components", {}).get("knee_min") is not None]
        line = "  cfg=%s mode=%s n=%d fell=%d bilateral=%d" % (cfg, mode, len(rs), fell, bilat)
        if ks:
            line += " kine[min=%.3f max=%.3f]" % (min(ks), max(ks))
        else:
            line += " kine=none"
        if duty:
            line += " duty[max=%.3f]" % max(duty)
        if kmin:
            line += " knee_min[min=%.1f]" % min(kmin)
        print(line)
        scored = [r for r in rs if r["metrics"].get("kine_score") is not None]
        if scored:
            best = min(scored, key=lambda r: r["metrics"].get("kine_score", 0))
            kc = best["metrics"].get("kine_components", {})
            print("    best: kine=%.3f knee_min=%.1f duty=%.3f ncyc_r=%s ncyc_l=%s T_r=%s bilateral=%s argextra=%s" % (
                best["metrics"]["kine_score"], kc.get("knee_min", float("nan")), kc.get("duty", float("nan")),
                kc.get("n_cycles_r"), kc.get("n_cycles_l"), kc.get("T_r"), kc.get("bilateral"),
                best.get("argv_extra")))

for path in sorted(glob.glob(SP + r"\robust_results_*.jsonl")):
    rows = [json.loads(l) for l in open(path) if l.strip()]
    print("=" * 60)
    print(os.path.basename(path), "rows:", len(rows))
    for r in rows:
        mt = r["metrics"]
        ps = mt.get("bal_push_sway")
        print("  variant=%s dty=%s mode=%s push=%-3s fell=%s tilt_max=%.2f sway_rms=%.4f push_sway=%s kine=%.2f" % (
            r.get("variant"), r.get("dty"), r.get("mode"), r.get("push"),
            mt.get("bal_fell"), mt.get("bal_tilt_max", float("nan")),
            mt.get("bal_sway_rms", float("nan")),
            ("%.4f" % ps) if ps is not None else "n/a",
            mt.get("kine_score", float("nan"))))
