import json, glob, csv, collections

D = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\tmp\gate_20261016"

def near(v, t, tol=1e-6):
    return v is not None and abs(v - t) < tol

for pat in ("round8_w*.jsonl", "round6_w*.jsonl"):
    rows = []
    for f in sorted(glob.glob(D + "\\" + pat)):
        for ln in open(f, encoding="utf-8", errors="replace"):
            ln = ln.strip()
            if ln.startswith("{"):
                rows.append(json.loads(ln))
    print(f"\n===== {pat}: {len(rows)} trials =====")
    sent = collections.Counter()
    for r in rows:
        s = r.get("kine_score")
        if r.get("nan"): sent["nan(-400?)"] += 1; continue
        if near(s, -25.0): sent["frozen -25"] += 1; continue
        if near(s, -320.0): sent["no-cycles -320"] += 1; continue
        if near(s, -60.0): sent["bilateral-gate -60"] += 1; continue
        sent["genuine"] += 1
    print("score buckets:", dict(sent))
    bil = [r for r in rows if r.get("k_bilateral") and not r.get("nan")
           and r.get("kine_score") is not None and abs(r["kine_score"] + 25.0) > 1e-6
           and abs(r["kine_score"] + 320.0) > 1e-6 and abs(r["kine_score"] + 60.0) > 1e-6]
    bil.sort(key=lambda r: r["kine_score"])
    print(f"bilateral non-sentinel trials: {len(bil)}; top 5:")
    for r in bil[:5]:
        print(f"  tag={r.get('tag')} score={r['kine_score']:.5f} cyc={r.get('k_n_cycles_r')}/{r.get('k_n_cycles_l')} "
              f"duty={None if r.get('k_duty_r') is None else round(r['k_duty_r'],3)}/{None if r.get('k_duty_l') is None else round(r['k_duty_l'],3)} "
              f"scale={r.get('scale')} nan={r.get('nan')}")
    # minimum-cycle winners check
    mins = [r for r in bil if (r.get("k_n_cycles_r") or 0) <= 3 or (r.get("k_n_cycles_l") or 0) <= 3]
    print(f"bilateral trials sitting AT the 3-cycle gate minimum: {len(mins)}")
    # scale distribution of top 10
    print("top-10 scales:", [round(r.get("scale", -1), 3) for r in bil[:10]])

print("\n===== wean_trials.csv (robot) =====")
with open(D + r"\wean_trials.csv", newline="", encoding="utf-8", errors="replace") as fh:
    rd = csv.DictReader(fh)
    rows = list(rd)
print("columns:", rd.fieldnames)
print("total rows:", len(rows))
studies = collections.Counter(r.get("study", "?") for r in rows)
print("by study:", dict(studies))
r060 = [r for r in rows if r.get("study", "").endswith("r060")]
print(f"\nr060 rows: {len(r060)}")
def fnum(r, k):
    try: return float(r[k])
    except Exception: return None
r060v = [(fnum(r, "obj"), r) for r in r060]
r060v.sort(key=lambda x: -(x[0] if x[0] is not None else -1e9))
print("r060 top 5 by obj:")
for obj, r in r060v[:5]:
    print(f"  trial={r.get('trial')} obj={obj} survive={r.get('survive')} cyc={r.get('cycL')}/{r.get('cycR')} phase={r.get('phase_med')} slope={r.get('slope')}")
dead = sum(1 for r in r060 if fnum(r, "obj") is not None and fnum(r, "obj") <= 5.001)
print(f"r060 dead-sentinel rows (obj<=5): {dead}")
